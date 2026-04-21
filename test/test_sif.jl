# Phase 3a — SIF injection + data-loader smoke tests
# =========================================================================
#
# Verifies:
#   - inject_surface_SIF! adds the expected (1/π) * SIF₀ * exp(-τ/μ)
#     contribution to R_SFI for a Lambertian surface, m=0.
#   - Non-Lambertian surfaces and m>0 leave the source untouched.
#   - SIF loaders (load_sif_spectrum, load_ficus_reflectance,
#     build_sif_source) run without error on the bundled fixtures.

using vSmartMOM, vSmartMOM.CoreRT, vSmartMOM.InelasticScattering
using Test

@testset "Phase 3a SIF injection + loaders" begin
    params = parameters_from_yaml("test_parameters/Phase1b_RRS_761-764nm.yaml")
    params.architecture = vSmartMOM.Architectures.CPU()
    model = model_from_parameters(params)
    FT = CoreRT.float_type(model)
    pol_type = CoreRT.polarization_type(model)
    nSpec = length(model.atmosphere.spec_bands[1])

    @testset "(1/π)·SIF₀ injection magnitude" begin
        SIF_I = FT(0.01)

        rs0 = InelasticScattering.noRS{FT}()
        rs0.F₀ = zeros(FT, pol_type.n, nSpec); rs0.F₀[1, :] .= 1
        rs0.SIF₀ = zeros(FT, pol_type.n, nSpec)
        R0 = CoreRT.rt_run_test_ss(rs0, model, 1)[1]

        rs1 = InelasticScattering.noRS{FT}()
        rs1.F₀ = zeros(FT, pol_type.n, nSpec); rs1.F₀[1, :] .= 1
        rs1.SIF₀ = zeros(FT, pol_type.n, nSpec); rs1.SIF₀[1, :] .= SIF_I
        R1 = CoreRT.rt_run_test_ss(rs1, model, 1)[1]

        δ = R1[1, 1, 50] - R0[1, 1, 50]
        # Expected: (1/π)·SIF_I·exp(-τ_atm/μ_view). τ_atm ≈ 0.025 Rayleigh-only
        # at 763 nm, μ_view=1 at vza=0 → attenuation ≈ 0.975.
        expected = SIF_I / π
        # The attenuation is bounded by [0.95, 1] for this scenario.
        @test 0.9  * expected <  δ <  1.0 * expected
        # Q/U untouched (isotropic unpolarized SIF).
        for ipol in 2:pol_type.n
            @test R1[1, ipol, 50] ≈ R0[1, ipol, 50] atol=1e-8
        end
    end

    @testset "m>0 SIF is zero (isotropic)" begin
        # Inject into m=0 only — verify j₀⁻ unchanged for m=1 via direct call.
        arch = vSmartMOM.Architectures.CPU()
        FT = Float32
        SIF₀ = zeros(FT, 3, 10); SIF₀[1, :] .= FT(0.01)
        brdf = CoreRT.LambertianSurfaceScalar(FT(0.0))

        quad = model.quad_points
        Nquad = quad.Nquad
        NquadN = Nquad * pol_type.n
        added = CoreRT.make_added_layer(
            InelasticScattering.noRS{FT}(), FT, vSmartMOM.Architectures.array_type(arch),
            (NquadN, NquadN), 10)

        # m=0: injection should hit
        fill!(added.j₀⁻, zero(FT))
        CoreRT.inject_surface_SIF!(brdf, added, 0, pol_type, SIF₀, arch)
        @test maximum(abs, added.j₀⁻) > 0
        @test added.j₀⁻[1, 1, 1] ≈ FT(2) * SIF₀[1, 1]

        # m=1: no injection
        fill!(added.j₀⁻, zero(FT))
        CoreRT.inject_surface_SIF!(brdf, added, 1, pol_type, SIF₀, arch)
        @test all(iszero, added.j₀⁻)
    end

    @testset "non-Lambertian surfaces: inject no-op" begin
        arch = vSmartMOM.Architectures.CPU()
        FT = Float32
        SIF₀ = zeros(FT, 3, 10); SIF₀[1, :] .= FT(0.01)

        quad = model.quad_points
        Nquad = quad.Nquad
        NquadN = Nquad * pol_type.n
        added = CoreRT.make_added_layer(
            InelasticScattering.noRS{FT}(), FT, vSmartMOM.Architectures.array_type(arch),
            (NquadN, NquadN), 10)
        fill!(added.j₀⁻, zero(FT))

        # rpv surface — not Lambertian
        rpv = CoreRT.rpvSurfaceScalar{FT}(FT(0.1), FT(0.1), FT(0.7), FT(-0.1))
        CoreRT.inject_surface_SIF!(rpv, added, 0, pol_type, SIF₀, arch)
        @test all(iszero, added.j₀⁻)
    end

    @testset "SIF₀ = nothing dispatch no-op" begin
        arch = vSmartMOM.Architectures.CPU()
        FT = Float32
        quad = model.quad_points
        NquadN = quad.Nquad * pol_type.n
        added = CoreRT.make_added_layer(
            InelasticScattering.noRS{FT}(), FT, vSmartMOM.Architectures.array_type(arch),
            (NquadN, NquadN), 10)
        added.j₀⁻ .= FT(0.7)

        brdf = CoreRT.LambertianSurfaceScalar(FT(0.0))
        CoreRT.inject_surface_SIF!(brdf, added, 0, pol_type, nothing, arch)
        @test all(added.j₀⁻ .≈ FT(0.7))
    end

    @testset "load_sif_spectrum smoke" begin
        ν, j = vSmartMOM.load_sif_spectrum()
        @test issorted(ν)
        @test length(ν) == length(j)
        @test all(j .≥ 0)
        @test 11_000 < ν[1] < 13_000   # 640 nm ≈ 15625, 850 nm ≈ 11765
        @test 15_000 < ν[end] < 16_000

        # Rescale-off path still returns something positive.
        ν2, j2 = vSmartMOM.load_sif_spectrum(rescale_to_peak=false)
        @test length(ν2) == length(ν)
        @test maximum(j2) > 0
    end

    @testset "load_ficus_reflectance smoke" begin
        λ, R = vSmartMOM.load_ficus_reflectance()
        @test length(λ) == length(R)
        @test 0.59 < λ[1] < 0.61
        @test 0.79 < λ[end] < 0.81
        @test all(0 .≤ R .≤ 1)
    end

    @testset "build_sif_source populates RS_type.SIF₀" begin
        rs = InelasticScattering.noRS{FT}()
        rs.F₀ = zeros(FT, pol_type.n, nSpec); rs.F₀[1, :] .= 1
        rs.SIF₀ = zeros(FT, pol_type.n, nSpec)
        ν_sif, jSIF = vSmartMOM.load_sif_spectrum()
        ν_model = collect(model.atmosphere.spec_bands[1])
        vSmartMOM.build_sif_source(rs, ν_model, ν_sif, jSIF; pol_component=1)
        @test all(rs.SIF₀[1, :] .> 0)
        @test all(rs.SIF₀[2, :] .== 0)  # untouched
    end
end
