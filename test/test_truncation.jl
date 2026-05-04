# Phase-function truncation invariants
#
# Validates the AbstractTruncationType hierarchy added for the
# Sanghavi & Stephens 2015 (JQSRT 159, 53–68) δ-m / δ-fit family.
#
# Tests cover:
#   1. NoTruncation is the identity map on AerosolOptics.
#   2. Canopy invariant — δBGE on a smooth canopy phase function should
#      give the same RT result as NoTruncation, because canopy phase
#      functions have no sharp forward peak (paper Eq. 8 with f_tr → 0
#      collapses to the identity).
#   3. δ-m absorption budget: τ·(1 − ω) is invariant under δ-scaling
#      (algebraic consequence of Eqs. 8 in the paper). This test is
#      written against the proper τ/ω rescaling and is currently
#      `@test_skip`'d because the rescaling lines are still commented
#      out in `truncate_phase` — flip back to `@test` once they're
#      restored.

using Test
using vSmartMOM
using vSmartMOM.CoreRT
using vSmartMOM.Scattering
using CanopyOptics

@testset "Truncation invariants" begin

    @testset "NoTruncation is identity on AerosolOptics" begin
        # Synthetic AerosolOptics with arbitrary GreekCoefs.
        l = 8
        g = Scattering.GreekCoefs(
            randn(l), randn(l), randn(l),
            randn(l), randn(l), randn(l))
        aero = AerosolOptics(greek_coefs=g, ω̃=0.95, k=1.7, fᵗ=0.0)
        same = Scattering.truncate_phase(NoTruncation(), aero)
        @test same === aero          # actually returns the same object
        @test same.greek_coefs === g

        same2 = Scattering.truncate_phase_lowconf(NoTruncation(), aero)
        @test same2 === aero
    end

    @testset "l_max accessor" begin
        @test Scattering.l_max(NoTruncation()) == typemax(Int)
        @test Scattering.l_max(NoTruncation(; l_max = 32)) == 32
        @test Scattering.l_max(δBGE{Float64}(20, 2.0)) == 20
    end

    @testset "Canopy invariant: δBGE ≡ NoTruncation for canopy" begin
        # PROSPECT-PRO leaf optics at 685 / 800 nm; planophile LD; black
        # soil; LAI = 4 — same setup as
        # sandbox/canopy_foursail_*_685_800nm.jl.
        leaf = CanopyOptics.LeafProspectProProperties(
            N=1.5, Ccab=40.0, Ccar=8.0, Canth=0.0, Cbrown=0.0,
            Cw=0.012, Cm=0.009, Cprot=0.0, Ccbc=0.0)
        opti = CanopyOptics.createLeafOpticalStruct(400.0:1.0:2500.0)
        T_grid, R_grid = CanopyOptics.prospect(leaf, opti)
        λgrid = [Float64(v.val) for v in opti.λ]
        gix = [findmin(abs.(λgrid .- λ))[2] for λ in (685.0, 800.0)]
        R_leaf = R_grid[gix]; T_leaf = T_grid[gix]

        function build(trunc_str)
            rt = Dict{String,Any}(
                "spec_bands" => [string("[", join(1e7 ./ (685.0, 800.0), " "), "]")],
                "surface" => ["LambertianSurfaceScalar(0.0)"],
                "quadrature_type" => "RadauQuad()",
                "polarization_type" => "Stokes_I()",
                "max_m" => 8, "Δ_angle" => 2.0, "l_trunc" => 32,
                "depol" => -1.0, "float_type" => "Float64", "architecture" => "CPU()")
            if trunc_str !== nothing
                rt["truncation"] = trunc_str
            end
            params = parameters_from_dict(Dict{String,Any}(
                "radiative_transfer" => rt,
                "geometry" => Dict{String,Any}("sza" => 30.0,
                    "vza" => [30.0, 0.0, 30.0],
                    "vaz" => [180.0, 0.0, 0.0],
                    "obs_alt" => 1000.0),
                "atmospheric_profile" => Dict{String,Any}(
                    "T" => [285.0], "p" => [1012.99, 1013.0],
                    "profile_reduction" => -1)))
            params.brdf[1] = CanopySurface(;
                soil = LambertianSurfaceScalar(0.0), LAI = 4.0, n_layers = 1,
                LAD = CanopyOptics.planophile_leaves2(Float64),
                leaf_reflectance = R_leaf, leaf_transmittance = T_leaf,
                leaf_optics_grid = [685.0, 800.0], grid_unit = :nm)
            model = model_from_parameters(params)
            for τ in model.τ_rayl; fill!(τ, 0); end
            for τ in model.τ_abs;  fill!(τ, 0); end
            for τ in model.τ_aer;  fill!(τ, 0); end
            R_v, _ = rt_run(model)
            return (π / cosd(30.0)) .* R_v[:, 1, :]
        end

        v_default = build(nothing)            # legacy default = δBGE(32, 2.0)
        v_notrunc = build("NoTruncation()")    # explicit NoTruncation
        v_dbge    = build("δBGE{Float64}(32, 2.0)")  # explicit δBGE

        # All three must agree exactly: canopy has no forward peak, so
        # the δBGE retained-fraction c₀ is 1 (f_tr = 0) and Eq. 8
        # collapses to the identity.
        @test v_default ≈ v_notrunc atol = 1e-12
        @test v_default ≈ v_dbge    atol = 1e-12
    end

    # Pending invariant — δ-m / δ-BGE absorption budget (Sanghavi &
    # Stephens 2015 Eq. 8):
    #
    #     τ*·(1 − ω*) = τ·(1 − f_tr·ω)·(1 − ω·(1−f_tr)/(1−f_tr·ω))
    #                 = τ·(1 − f_tr·ω − ω + f_tr·ω)
    #                 = τ·(1 − ω).
    #
    # The current `truncate_phase(::δBGE, ...)` returns the truncated
    # Greek coefficients but does NOT rescale ω̃ or k (lines 231–235 of
    # `truncate_phase.jl` are commented out — see Eq. 8 in the paper).
    # Restore those, then implement this test against a real Mie
    # AerosolOptics from `compute_aerosol_optical_properties`. Marking
    # as @test_skip with a stub assertion until the rescaling lands.
    @testset "δBGE absorption budget τ(1-ω) — pending τ/ω rescaling" begin
        @test_skip false  # placeholder; see comment block above
    end
end
