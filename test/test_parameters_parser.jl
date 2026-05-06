using Test
using vSmartMOM
using vSmartMOM.IO
using vSmartMOM.CoreRT
using CanopyOptics

function _minimal_parameter_dict(; surface = ["LambertianSurfaceScalar(0.1)"],
                                   float_type = "Float64",
                                   architecture = "CPU()",
                                   quadrature_type = "GaussQuadHemisphere()")
    return Dict(
        "radiative_transfer" => Dict(
            "spec_bands" => ["[13000.0, 13000.1]"],
            "surface" => surface,
            "quadrature_type" => quadrature_type,
            "polarization_type" => "Stokes_I()",
            "max_m" => 1,
            "Δ_angle" => 0.0,
            "l_trunc" => 20,
            "depol" => 0.0,
            "float_type" => float_type,
            "architecture" => architecture,
        ),
        "geometry" => Dict(
            "sza" => 30.0,
            "vza" => [0.0],
            "vaz" => [0.0],
            "obs_alt" => 1000.0,
        ),
        "atmospheric_profile" => Dict(
            "T" => [290.0],
            "p" => [1000.0, 0.0],
            "profile_reduction" => 1,
        ),
    )
end

@testset "valid minimal dict" begin
    params = parameters_from_dict(_minimal_parameter_dict())

    @test params.float_type === Float64
    @test params.architecture isa vSmartMOM.Architectures.CPU
    @test params.brdf[1] isa LambertianSurfaceScalar{Float64}
    @test params.spec_bands == [[13000.0, 13000.1]]

    metal_params = parameters_from_dict(_minimal_parameter_dict(
        float_type = "Float32",
        architecture = "MetalGPU()",
    ))
    @test metal_params.float_type === Float32
    @test metal_params.architecture isa vSmartMOM.Architectures.MetalGPU
end

@testset "TOML file input" begin
    mktempdir() do dir
        path = joinpath(dir, "params.toml")
        open(path, "w") do io
            write(io, """
                [radiative_transfer]
                spec_bands = ["[13000.0, 13000.1]"]
                surface = ["LambertianSurfaceScalar(0.1)"]
                quadrature_type = "GaussQuadHemisphere()"
                polarization_type = "Stokes_I()"
                max_m = 1
                "Δ_angle" = 0.0
                l_trunc = 20
                depol = 0.0
                float_type = "Float64"
                architecture = "CPU()"

                [geometry]
                sza = 30.0
                vza = [0.0]
                vaz = [0.0]
                obs_alt = 1000.0

                [atmospheric_profile]
                T = [290.0]
                p = [1000.0, 0.0]
                profile_reduction = 1
            """)
        end

        params = parameters_from_file(path)
        @test params.float_type === Float64
        @test params.brdf[1] isa LambertianSurfaceScalar{Float64}
        @test read_parameters(path).spec_bands == [[13000.0, 13000.1]]
        @test_throws ArgumentError parameters_from_yaml(path)
    end
end

@testset "surface dispatch parser" begin
    surf = vSmartMOM.IO.parse_surface_str("CoxMunkSurface(wind_speed=5.0)", Float32)

    @test surf isa CoxMunkSurface{Float32}
    @test surf.wind_speed === Float32(5)
    @test_throws ArgumentError vSmartMOM.IO.parse_surface_str("UnknownSurface(1.0)", Float64)
    @test_throws ArgumentError vSmartMOM.IO.parse_surface_str("LambertianSurfaceScalar()", Float64)
    @test_throws ArgumentError vSmartMOM.IO.parse_surface_str("CoxMunkSurface()", Float64)
end

@testset "analytic phase-function aerosol parser" begin
    cfg = _minimal_parameter_dict()
    cfg["radiative_transfer"]["polarization_type"] = "Stokes_IQU()"
    cfg["radiative_transfer"]["max_m"] = 3
    cfg["radiative_transfer"]["l_trunc"] = 8
    cfg["atmospheric_profile"]["profile_reduction"] = -1
    cfg["atmospheric_profile"]["p"] = [100.0, 1000.0]
    cfg["scattering"] = Dict(
        "aerosols" => [Dict(
            "τ_ref" => 0.02,
            "phase_function" => "SyntheticPolarizedHenyeyGreensteinPhaseFunction(g=0.45, polarization_fraction=0.6)",
            "ssa" => 0.92,
            "p₀" => 550.0,
            "σp" => 200.0,
        )],
        "r_max" => 1.0,
        "nquad_radius" => 16,
        "λ_ref" => 0.55,
        "decomp_type" => "NAI2()",
    )

    params = parameters_from_dict(cfg)
    aer = params.scattering_params.rt_aerosols[1]
    @test aer.phase_function isa vSmartMOM.Scattering.SyntheticPolarizedHenyeyGreensteinPhaseFunction
    @test aer.ϖ == 0.92

    cfg_dict_phase = deepcopy(cfg)
    cfg_dict_phase["scattering"]["aerosols"][1]["phase_function"] = Dict(
        "type" => "HenyeyGreensteinPhaseFunction",
        "g" => 0.35,
    )
    delete!(cfg_dict_phase["scattering"]["aerosols"][1], "ssa")
    cfg_dict_phase["scattering"]["aerosols"][1]["ϖ"] = 0.85
    params_dict_phase = parameters_from_dict(cfg_dict_phase)
    aer_dict_phase = params_dict_phase.scattering_params.rt_aerosols[1]
    @test aer_dict_phase.phase_function isa vSmartMOM.Scattering.HenyeyGreensteinPhaseFunction
    @test aer_dict_phase.phase_function.g === 0.35
    @test aer_dict_phase.ϖ == 0.85

    model = model_from_parameters(params)
    @test model.aerosol_optics[1][1].ω̃ ≈ 0.92
    @test length(model.aerosol_optics[1][1].greek_coefs.β) == params.l_trunc
    @test all(model.τ_aer[1][1, :] .>= 0)
    @test sum(model.τ_aer[1][1, :]) ≈ aer.τ_ref
    @test_throws ArgumentError model_from_parameters(LinMode(), params)

    cfg_no_trunc = deepcopy(cfg)
    cfg_no_trunc["radiative_transfer"]["truncation"] = "NoTruncation()"
    params_no_trunc = parameters_from_dict(cfg_no_trunc)
    model_no_trunc = model_from_parameters(params_no_trunc)
    @test length(model_no_trunc.aerosol_optics[1][1].greek_coefs.β) == params_no_trunc.l_trunc
    @test model_no_trunc.aerosol_optics[1][1].fᵗ == 0

    cfg_no_trunc_lmax = deepcopy(cfg)
    cfg_no_trunc_lmax["radiative_transfer"]["truncation"] = "NoTruncation(l_max=5)"
    model_no_trunc_lmax = model_from_parameters(parameters_from_dict(cfg_no_trunc_lmax))
    @test length(model_no_trunc_lmax.aerosol_optics[1][1].greek_coefs.β) == 5

    bad_cfg = deepcopy(cfg)
    bad_cfg["scattering"]["aerosols"][1]["phase_function"] = "UnknownPhaseFunction(g=0.2)"
    @test_throws ArgumentError parameters_from_dict(bad_cfg)

    bad_microphysics_cfg = deepcopy(cfg)
    bad_microphysics_cfg["scattering"]["aerosols"][1]["μ"] = 0.1
    @test_throws ArgumentError parameters_from_dict(bad_microphysics_cfg)
end

@testset "canopy clumping parser" begin
    cfg = _minimal_parameter_dict()
    cfg["canopy"] = Dict(
        "LAI" => 2.0,
        "clumping" => 0.65,
    )
    params = parameters_from_dict(cfg)
    canopy = params.brdf[1]
    @test canopy isa CanopySurface{Float64}
    @test canopy.canopy_clumping isa CanopyOptics.ConstantClumping{Float64}
    @test CanopyOptics.clumping_index(canopy.canopy_clumping, 0.5) == 0.65
    @test canopy.canopy_quadrature.n_azimuth == CanopyOptics.CanopyQuadrature().n_azimuth

    cfg = _minimal_parameter_dict(float_type = "Float32")
    cfg["canopy"] = Dict(
        "LAI" => 2.0,
        "clumping" => Dict("type" => "chen_leblanc", "Omega0" => 0.6, "c" => 1.5, "e" => 2.0),
        "n_leaf_quadrature" => 12,
    )
    params = parameters_from_dict(cfg)
    canopy = params.brdf[1]
    @test canopy.canopy_clumping isa CanopyOptics.ChenLeblancClumping{Float32}
    @test CanopyOptics.clumping_index(canopy.canopy_clumping, Float32(1)) ≈ Float32(0.6)
    @test canopy.canopy_quadrature.n_leaf == 12
    @test canopy.canopy_quadrature.n_azimuth == CanopyOptics.CanopyQuadrature().n_azimuth
end

@testset "environment path expansion" begin
    env_key = "VSMARTMOM_TEST_LUT_DIR"
    old_value = get(ENV, env_key, nothing)

    try
        ENV[env_key] = "/tmp/vsmartmom-luts"
        @test vSmartMOM.IO._expand_env_path("\${ENV:$(env_key)}/O2.jld2") == normpath("/tmp/vsmartmom-luts/O2.jld2")
        ENV[env_key] = "/tmp/vsmartmom-luts/O2.jld2"
        @test vSmartMOM.IO._expand_env_path("\${ENV:$(env_key)}") == normpath("/tmp/vsmartmom-luts/O2.jld2")
        @test vSmartMOM.IO._expand_env_path("relative/O2.jld2") == "relative/O2.jld2"

        delete!(ENV, env_key)
        @test_throws ArgumentError vSmartMOM.IO._expand_env_path("\${ENV:$(env_key)}/O2.jld2")
    finally
        if old_value === nothing
            delete!(ENV, env_key)
        else
            ENV[env_key] = old_value
        end
    end
end

@testset "invalid configuration raises ArgumentError" begin
    cfg = _minimal_parameter_dict(float_type = "Float16")
    @test_throws ArgumentError parameters_from_dict(cfg)

    cfg = _minimal_parameter_dict(architecture = "TPU()")
    @test_throws ArgumentError parameters_from_dict(cfg)

    cfg = _minimal_parameter_dict(quadrature_type = "BadQuad()")
    @test_throws ArgumentError parameters_from_dict(cfg)

    cfg = _minimal_parameter_dict()
    delete!(cfg["geometry"], "sza")
    @test_throws ArgumentError parameters_from_dict(cfg)
end
