# =========================================================================
# Lambertian surface-flavour consolidation tests
# =========================================================================
#
# Exercises the unified `create_surface_layer!` scaffold shared by the four
# Lambertian flavours (Scalar / Spectrum / Legendre / Spline) plus the
# `surface_albedo` kernels. Verifies:
#
#   1. Cross-flavour equivalence of the surface AddedLayer at m=0 and m=1:
#        • Spectrum with a constant ρ vector  ==  Scalar (bit-identical)
#        • Legendre with only the constant coeff == Scalar (bit-identical)
#      plus the agreed t / r / j conventions (t⁺⁺=t⁻⁻=I, r⁺⁻=0,
#      j₀⁺ = attenuated direct beam at m=0, everything zero at m>0).
#   2. LambertianSurfaceSpectrum end-to-end through `rt_run` equals the
#      equivalent Scalar config (bit-identical R).
#   3. Regression for the old copy-paste bug where the m>0 branch zeroed
#      r⁻⁺ twice and never r⁺⁻ — a pre-seeded r⁺⁻ must come back zero.
#   4. Spectrum albedo-length validation.
#
# All Lambertian flavours are exact at m=0, but Rayleigh forces the m loop
# up to 2, so the m>0 branch is live and must be well-defined.

using vSmartMOM
using vSmartMOM.CoreRT
using vSmartMOM.InelasticScattering
using Test
using LinearAlgebra
using Interpolations

const _LS = vSmartMOM.CoreRT

@testset "Lambertian surface flavours" begin

    arch     = vSmartMOM.Architectures.CPU()
    arr_type = vSmartMOM.Architectures.array_type(arch)

    # ---------------------------------------------------------------------
    # A small real model just to borrow a consistent QuadPoints / pol_type.
    # ---------------------------------------------------------------------
    params = parameters_from_yaml("test_parameters/PureRayleighParameters.yaml")
    params.architecture = arch
    model    = model_from_parameters(params)
    FT       = CoreRT.float_type(model)
    pol_type = CoreRT.polarization_type(model)
    quad     = model.quad_points
    Nquad    = quad.Nquad
    NquadN   = Nquad * pol_type.n

    nSpec  = 8
    albedo = FT(0.23)
    τ_sum  = arr_type(collect(range(FT(0.05), FT(0.4), length=nSpec)))

    # Helper: fresh AddedLayer, run create_surface_layer!, return it.
    function build_surface(brdf, m; SFI=true)
        added = CoreRT.make_added_layer(
            InelasticScattering.noRS{FT}(), FT, arr_type, (NquadN, NquadN), nSpec)
        CoreRT.create_surface_layer!(brdf, added, SFI, m, pol_type, quad, τ_sum, arch)
        return added
    end

    scalar   = CoreRT.LambertianSurfaceScalar(albedo)
    spectrum = CoreRT.LambertianSurfaceSpectrum(fill(albedo, nSpec))
    legendre = CoreRT.LambertianSurfaceLegendre([albedo])   # only constant P₀ term

    # ---------------------------------------------------------------------
    # surface_albedo kernels
    # ---------------------------------------------------------------------
    @testset "surface_albedo kernels" begin
        @test _LS.surface_albedo(scalar, τ_sum) == albedo
        @test _LS.surface_albedo(spectrum, τ_sum) == fill(albedo, nSpec)
        # Legendre P₀ ≡ 1, so a single constant coeff is flat == coeff.
        @test all(_LS.surface_albedo(legendre, τ_sum) .≈ albedo)

        # Spectrum length mismatch must throw.
        bad = CoreRT.LambertianSurfaceSpectrum(fill(albedo, nSpec + 1))
        @test_throws DimensionMismatch _LS.surface_albedo(bad, τ_sum)
    end

    # ---------------------------------------------------------------------
    # m = 0 conventions + cross-flavour bit-equality
    # ---------------------------------------------------------------------
    @testset "m=0 conventions" begin
        a0 = build_surface(scalar, 0)

        # t⁺⁺ = I (pass-through), t⁻⁻ = 0 (opaque lowest boundary), r⁺⁻ = 0.
        Id = arr_type(Matrix{FT}(I, NquadN, NquadN))
        @test a0.t⁺⁺[:, :, 1] == Id
        @test all(iszero, a0.t⁻⁻)
        @test all(iszero, a0.r⁺⁻)
        # r⁻⁺ is the (positive) reflectance — nonzero.
        @test maximum(abs, a0.r⁻⁺) > 0
        # j₀⁺ = attenuated direct beam (NOT zero — consumed by interaction_hdrf!).
        @test maximum(abs, a0.j₀⁺) > 0
        # j₀⁻ = reflected beam — nonzero.
        @test maximum(abs, a0.j₀⁻) > 0
    end

    @testset "Spectrum(const) == Scalar  (m=0, bit-identical)" begin
        a_sc = build_surface(scalar, 0)
        a_sp = build_surface(spectrum, 0)
        @test a_sp.r⁻⁺ == a_sc.r⁻⁺
        @test a_sp.r⁺⁻ == a_sc.r⁺⁻
        @test a_sp.t⁺⁺ == a_sc.t⁺⁺
        @test a_sp.t⁻⁻ == a_sc.t⁻⁻
        @test a_sp.j₀⁺ == a_sc.j₀⁺
        @test a_sp.j₀⁻ == a_sc.j₀⁻
    end

    @testset "Legendre(const) == Scalar  (m=0, bit-identical)" begin
        a_sc = build_surface(scalar, 0)
        a_le = build_surface(legendre, 0)
        @test a_le.r⁻⁺ == a_sc.r⁻⁺
        @test a_le.r⁺⁻ == a_sc.r⁺⁻
        @test a_le.t⁺⁺ == a_sc.t⁺⁺
        @test a_le.t⁻⁻ == a_sc.t⁻⁻
        @test a_le.j₀⁺ == a_sc.j₀⁺
        @test a_le.j₀⁻ == a_sc.j₀⁻
    end

    # Spline with a flat interpolator should also equal Scalar at m=0.
    @testset "Spline(const) == Scalar  (m=0, bit-identical)" begin
        wl   = collect(range(FT(0), FT(1), length=nSpec))
        itp  = extrapolate(interpolate((wl,), fill(albedo, nSpec),
                           Gridded(Linear())), Flat())
        spline = CoreRT.LambertianSurfaceSpline(itp, wl)
        a_sc = build_surface(scalar, 0)
        a_pl = build_surface(spline, 0)
        @test a_pl.r⁻⁺ == a_sc.r⁻⁺
        @test a_pl.j₀⁻ == a_sc.j₀⁻
        @test a_pl.j₀⁺ == a_sc.j₀⁺
    end

    # ---------------------------------------------------------------------
    # m > 0 conventions: everything zero; t⁺⁺ = I (pass-through), t⁻⁻ = 0 (opaque).
    # ---------------------------------------------------------------------
    @testset "m>0 conventions (all flavours)" begin
        Id = arr_type(Matrix{FT}(I, NquadN, NquadN))
        for brdf in (scalar, spectrum, legendre)
            a = build_surface(brdf, 1)
            @test all(iszero, a.r⁻⁺)
            @test all(iszero, a.r⁺⁻)
            @test all(iszero, a.j₀⁺)
            @test all(iszero, a.j₀⁻)
            @test a.t⁺⁺[:, :, 1] == Id
            @test all(iszero, a.t⁻⁻)
        end
    end

    # ---------------------------------------------------------------------
    # Regression: m>0 must zero r⁺⁻ (old code zeroed r⁻⁺ twice, never r⁺⁻).
    # ---------------------------------------------------------------------
    @testset "r⁺⁻ zeroing regression (m>0)" begin
        added = CoreRT.make_added_layer(
            InelasticScattering.noRS{FT}(), FT, arr_type, (NquadN, NquadN), nSpec)
        # Pre-seed BOTH reflectances with garbage; the m>0 branch must clear both.
        added.r⁻⁺ .= FT(7)
        added.r⁺⁻ .= FT(9)
        CoreRT.create_surface_layer!(scalar, added, true, 1, pol_type, quad, τ_sum, arch)
        @test all(iszero, added.r⁻⁺)
        @test all(iszero, added.r⁺⁻)   # would FAIL with the old copy-paste typo
    end

    # ---------------------------------------------------------------------
    # End-to-end: Spectrum(const) rt_run == Scalar rt_run (bit-identical R).
    # ---------------------------------------------------------------------
    @testset "Spectrum end-to-end == Scalar (rt_run)" begin
        nS = length(model.atmosphere.spec_bands[1])
        a  = FT(0.3)

        p_sc = parameters_from_yaml("test_parameters/PureRayleighParameters.yaml")
        p_sc.architecture = arch
        p_sc.brdf = CoreRT.AbstractSurfaceType[CoreRT.LambertianSurfaceScalar(a)]
        R_sc, _ = rt_run(model_from_parameters(p_sc))

        p_sp = parameters_from_yaml("test_parameters/PureRayleighParameters.yaml")
        p_sp.architecture = arch
        p_sp.brdf = CoreRT.AbstractSurfaceType[CoreRT.LambertianSurfaceSpectrum(fill(a, nS))]
        R_sp, _ = rt_run(model_from_parameters(p_sp))

        @test R_sp == R_sc          # bit-identical
        @test maximum(abs, R_sc) > 0
    end
end
