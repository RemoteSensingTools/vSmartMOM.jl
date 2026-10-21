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
#      plus the agreed t / r / j conventions (t⁺⁺=I, t⁻⁻=0, r⁺⁻=0,
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
using YAML

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

    @testset "Legendre coefficient Jacobians vs finite differences" begin
        coeff = FT[0.22, 0.035, -0.018]
        p = parameters_from_yaml("test_parameters/PureRayleighParameters.yaml")
        p.architecture = arch
        p.brdf = CoreRT.AbstractSurfaceType[CoreRT.LambertianSurfaceLegendre(coeff)]
        lin_model, optical_lin = model_from_parameters(LinMode(), p)
        result = rt_run(lin_model, optical_lin, 0,
                        size(optical_lin.τ̇_abs[1], 1), length(coeff))

        @test result.layout.n_surface == length(coeff)
        @test CoreRT.surface_parameter_count(p.brdf[1]) == length(coeff)
        δ = FT(1e-5)
        for k in eachindex(coeff)
            cp = copy(coeff); cp[k] += δ
            cm = copy(coeff); cm[k] -= δ
            pp = deepcopy(p); pp.brdf = CoreRT.AbstractSurfaceType[CoreRT.LambertianSurfaceLegendre(cp)]
            pm = deepcopy(p); pm.brdf = CoreRT.AbstractSurfaceType[CoreRT.LambertianSurfaceLegendre(cm)]
            Rp = rt_run(model_from_parameters(pp)).toa
            Rm = rt_run(model_from_parameters(pm)).toa
            fd = (Rp .- Rm) ./ (FT(2) * δ)
            col = CoreRT.surface_range(result.layout)[k]
            @test result.toa_jacobian[:, :, :, col] ≈ fd rtol=FT(2e-4) atol=FT(2e-7)
        end

        @test_throws ArgumentError rt_run(lin_model, optical_lin, 0,
                                           size(optical_lin.τ̇_abs[1], 1), 1)
    end
end

# B8 (docs/dev_notes/RELEASE_READINESS_AND_CLAUDE_HANDOFF_2026-10-16.md): a
# combined-band solve applies one lower boundary to every band. It used to take
# the first band's surface silently; differing or band-resolved surfaces must
# now be rejected.
@testset "Combined-band solves require a common surface" begin
    d = YAML.load_file(joinpath(pkgdir(vSmartMOM), "config", "quickstart.yaml"))
    d["radiative_transfer"]["spec_bands"] = ["[12987.0, 12988.0]", "[13100.0, 13101.0]"]
    d["radiative_transfer"]["surface"] = ["LambertianSurfaceScalar(0.15)",
                                          "LambertianSurfaceScalar(0.3)"]
    differing = model_from_parameters(read_parameters(d))
    rs = noRS{Float64}()
    @test_throws r"solved together with a single lower boundary" CoreRT.rt_run(rs, differing, [1, 2])
    @test_throws r"solved together with a single lower boundary" CoreRT.rt_run_ss(rs, differing, [1, 2])
    @test_throws r"solved together with a single lower boundary" CoreRT.rt_run_test_ms(rs, [1], differing, [1, 2])
    # Each band alone is unaffected.
    @test CoreRT.rt_run(rs, differing, [2]).toa ≈ rt_run(differing; i_band = 2).toa

    # Equal plain-value surfaces from separate entries are the same boundary,
    # and the elastic combined solve reproduces the separate band solves.
    d["radiative_transfer"]["surface"] = ["LambertianSurfaceScalar(0.15)",
                                          "LambertianSurfaceScalar(0.15)"]
    common = model_from_parameters(read_parameters(d))
    separate = cat(rt_run(common; i_band = 1).toa, rt_run(common; i_band = 2).toa; dims = 3)
    @test CoreRT.rt_run(rs, common, [1, 2]).toa ≈ separate rtol = 1e-12

    # A shared spectral albedo spanning both bands gives each band its own
    # albedo; it must match separate solves with those albedos.
    spanning = LambertianSurfaceSpectrum([0.15, 0.15, 0.3, 0.3])
    common.surfaces[1] = spanning
    common.surfaces[2] = spanning
    combined = CoreRT.rt_run(rs, common, [1, 2]).toa
    @test combined[:, :, 1:2] ≈ separate[:, :, 1:2] rtol = 1e-12
    @test combined[:, :, 3:4] ≈ rt_run(differing; i_band = 2).toa rtol = 1e-12
    one_band = LambertianSurfaceSpectrum([0.1, 0.2])
    common.surfaces[1] = one_band
    common.surfaces[2] = one_band
    @test_throws r"albedo values" CoreRT.rt_run(rs, common, [1, 2])

    # Legendre albedos and canopies cannot span bands.
    legendre = CoreRT.LambertianSurfaceLegendre([0.2, 0.1])
    common.surfaces[1] = legendre
    common.surfaces[2] = legendre
    @test_throws r"LambertianSurfaceLegendre" CoreRT.rt_run(rs, common, [1, 2])
    canopy = CanopySurface(; soil = LambertianSurfaceScalar(0.15), LAI = 2.0)
    common.surfaces[1] = canopy
    common.surfaces[2] = canopy
    @test_throws r"CanopySurface" CoreRT.rt_run(rs, common, [1, 2])
end
