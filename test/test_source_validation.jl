using Test
using vSmartMOM
using vSmartMOM.CoreRT

@testset "Source requests preserve their requested physics" begin
    params = read_parameters(joinpath(@__DIR__, "..", "config", "quickstart.yaml"))
    model = model_from_parameters(params)
    external = model_from_parameters(params; external_solar=true)
    lin_forward, lin_model = model_from_parameters(LinMode(), deepcopy(params))
    ngas = size(lin_model.τ̇_abs[1], 1)
    rpv = CoreRT.rpvSurfaceScalar(0.1, 1.0, 0.0, 0.1)
    sif = SurfaceSIF(SIF₀=fill(0.01, 1, 1))

    @testset "Every solve entry rejects ignored beam geometry and extra beams" begin
        solves = (
            s -> rt_run(model; sources=s),
            s -> rt_run_ss(model; sources=s),
            s -> rt_run_ss_exact(model; sources=s),
            s -> rt_run_streams(model; sources=s),
            s -> rt_run_toa(external; sources=s),
            s -> rt_run_atmosphere(model; sources=s),
            s -> rt_run(lin_forward, lin_model, 0, ngas, 1; sources=s),
        )
        for solve in solves
            @test_throws r"does not match the model solar zenith angle" solve(SolarBeam(sza=75))
            @test_throws r"Multiple SolarBeam sources" solve(SolarBeam() + SolarBeam())
        end
        blackbody = BlackbodySource(1500, [12987.0]; pol_n=1)
        @test_throws r"Multiple SolarBeam sources" rt_run(model; sources=SolarBeam() + blackbody)

        # Stored sources and explicit overrides use the same validation.
        bad_model = model_from_parameters(params; sources=SolarBeam(sza=75))
        @test_throws r"does not match the model solar zenith angle" rt_run(bad_model)
        @test rt_run(bad_model; sources=SolarBeam()).toa == rt_run(model).toa

        # Interior observers take a separate driver branch.
        interior_params = deepcopy(params)
        interior_params.obs_alt = 5.0
        interior = model_from_parameters(interior_params)
        @test_throws r"does not match the model solar zenith angle" rt_run(interior; sources=SolarBeam(sza=75))
        @test_throws r"Multiple SolarBeam sources" rt_run(interior; sources=SolarBeam() + SolarBeam())
    end

    @testset "Matching geometry and summed irradiance remain valid" begin
        reference = rt_run(model).toa
        @test rt_run(model; sources=SolarBeam(sza=params.sza)).toa == reference
        @test rt_run(model; sources=SolarBeam(F₀=fill(2.0, 1, 1))).toa ≈ 2reference
        @test iszero(rt_run(model; sources=NoSource()).toa)
        for FT in (Float32, Float64)
            surface = LambertianSurfaceScalar(FT(0.15))
            @test isnothing(CoreRT.validate_source_requests(SolarBeam(sza=60.1), FT(60.1), surface))
            @test_throws ArgumentError CoreRT.validate_source_requests(SolarBeam(sza=NaN), FT(60.1), surface)
            for SourceFT in (Float32, Float64)
                @test isnothing(CoreRT.validate_source_requests(
                    SolarBeam(sza=SourceFT(35.1)), FT(35.1), surface))
                @test_throws ArgumentError CoreRT.validate_source_requests(
                    SolarBeam(sza=SourceFT(35.101)), FT(35.1), surface)
                @test_throws ArgumentError CoreRT.validate_source_requests(
                    SolarBeam(sza=SourceFT(Inf)), FT(35.1), surface)
            end
        end
    end

    @testset "Exact atmospheric single scattering honors source irradiance" begin
        reference = rt_run_ss_exact(model)
        doubled = model_from_parameters(deepcopy(params); sources=SolarBeam(F₀=fill(2.0, 1, 1)))
        dark = model_from_parameters(deepcopy(params); sources=NoSource())
        @test rt_run_ss_exact(doubled) ≈ 2reference
        @test iszero(rt_run_ss_exact(dark))
        @test rt_run_ss_exact(dark; sources=SolarBeam()) == reference
        @test_throws r"atmosphere-only" rt_run_ss_exact(model; sources=SolarBeam() + sif)
        @test_throws r"atmosphere-only" rt_run_ss_exact(model; sources=ThermalEmission())
        @test_throws r"atmosphere-only" rt_run_ss_exact(model;
            sources=SolarBeam() + SurfaceSIF(SIF760=0.0, wavenumber_cm1=[12987.0]))
        @test rt_run_ss_exact(model; sources=SolarBeam() + SurfaceSIF()) == reference
        polarized_params = deepcopy(params)
        polarized_params.polarization_type = vSmartMOM.Scattering.Stokes_IQU()
        polarized = model_from_parameters(polarized_params)
        @test_throws r"requires an unpolarized solar source" rt_run_ss_exact(
            polarized; sources=SolarBeam(F₀=reshape([1.0, 0.1, 0.0], 3, 1)))
    end

    @testset "Unsupported SIF surfaces fail in forward, linearized and cache replay" begin
        rpv_params = deepcopy(params)
        rpv_params.brdf[1] = rpv
        rpv_model = model_from_parameters(rpv_params)
        for solve in (s -> rt_run(rpv_model; sources=s),
                      s -> rt_run_ss(rpv_model; sources=s),
                      s -> rt_run_atmosphere(rpv_model; sources=s))
            @test_throws r"SurfaceSIF.*unsupported" solve(SolarBeam() + sif)
        end
        # Surface changes at solve time must also validate the linearized path.
        lin_forward.surfaces[1] = rpv
        nsurf = CoreRT.surface_parameter_count(rpv)
        @test_throws r"SurfaceSIF.*unsupported" rt_run(
            lin_forward, lin_model, 0, ngas, nsurf; sources=SolarBeam() + sif)

        cache = rt_run_atmosphere(model; sources=SolarBeam() + sif)
        @test_throws r"SurfaceSIF.*unsupported" rt_run_surface(cache, rpv)
        @test rt_run_surface(cache, model.surfaces[1])[1] == rt_run(model; sources=SolarBeam() + sif).toa
        @test all(rt_run(model; sources=SolarBeam() + sif).toa .> rt_run(model).toa)

        # A zero prescribed source is allowed; a zero retrieval amplitude still
        # has a nonzero derivative and needs an implemented injection method.
        for zero_sif in (SurfaceSIF(), SurfaceSIF(SIF₀=zeros(1, 1)))
            @test isnothing(CoreRT.validate_source_requests(SolarBeam() + zero_sif, params.sza, rpv))
        end
        retrievable = SurfaceSIF(SIF760=0.0, wavenumber_cm1=[12987.0])
        @test_throws r"SurfaceSIF.*unsupported" rt_run(rpv_model; sources=SolarBeam() + retrievable)
        prepared = prepare_sources(SolarBeam() + retrievable, Float64, 1, 1, Array)
        @test_throws r"SurfaceSIF.*unsupported" CoreRT.validate_source_surface(prepared, rpv)
        @test !supports_surface_sif(rpv)
        @test supports_surface_sif(model.surfaces[1])
    end
end
