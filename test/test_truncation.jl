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

    @testset "NoTruncation aerosol passthrough resets fᵗ" begin
        # Raw Mie outputs initialise `fᵗ = 1` as a "untruncated yet"
        # sentinel; downstream `delta_m_forward` interprets a literal 1
        # as "everything is in the forward peak" and zeros out
        # scattering. NoTruncation must therefore return `fᵗ = 0`,
        # not the raw `fᵗ = 1`.
        l = 8
        g = Scattering.GreekCoefs(zeros(l), ones(l), zeros(l),
                                  ones(l), zeros(l), ones(l))
        raw = AerosolOptics(greek_coefs=g, ω̃=0.95, k=1.7, fᵗ=1.0)
        out = Scattering.truncate_phase(NoTruncation(), raw)
        @test out.fᵗ == 0
        @test out.ω̃ == raw.ω̃
        @test out.k == raw.k
        @test out.greek_coefs === raw.greek_coefs
    end

    @testset "NoTruncation lin passthrough resets fᵗ and ḟᵗ" begin
        # Same sentinel issue applies to the linearised
        # (Jacobian) path: `delta_m_truncation_lin` reads `fᵗ` and
        # `ḟᵗ` and would zero the SSA Jacobian if `fᵗ = 1` leaked
        # through. NoTruncation must reset both.
        l = 8
        g  = Scattering.GreekCoefs(zeros(l), ones(l), zeros(l),
                                   ones(l), zeros(l), ones(l))
        lg = Scattering.linGreekCoefs(zeros(2, l), zeros(2, l), zeros(2, l),
                                      zeros(2, l), zeros(2, l), zeros(2, l))
        raw     = AerosolOptics(greek_coefs=g, ω̃=0.95, k=1.7, fᵗ=1.0)
        raw_lin = Scattering.linAerosolOptics(lin_greek_coefs=lg,
                      ω̃̇=zeros(2), k̇=zeros(2), ḟᵗ=ones(2))   # raw lin sentinel
        out, lout = Scattering.truncate_phase(NoTruncation(), raw, raw_lin)
        @test out.fᵗ == 0
        @test all(lout.ḟᵗ .== 0)
        @test out.greek_coefs === raw.greek_coefs
        @test lout.lin_greek_coefs === raw_lin.lin_greek_coefs
    end

    @testset "Explicit truncation is preserved (no silent rebuild)" begin
        # Regression for the P2 finding from codex review: setting
        # `params.truncation = δBGE(40, 5.0)` with legacy `Δ_angle=2.0`
        # `l_trunc=20` must keep δBGE(40, 5.0); the model-builder must
        # NOT silently rebuild from the legacy fields.
        rt = Dict{String,Any}(
            "spec_bands" => ["[14492 14493]"],
            "surface" => ["LambertianSurfaceScalar(0.1)"],
            "quadrature_type" => "RadauQuad()",
            "polarization_type" => "Stokes_I()",
            "max_m" => 8, "Δ_angle" => 2.0, "l_trunc" => 20,
            "depol" => -1.0, "float_type" => "Float64", "architecture" => "CPU()")
        p = parameters_from_dict(Dict{String,Any}(
            "radiative_transfer" => rt,
            "geometry" => Dict{String,Any}("sza" => 30.0, "vza" => [0.0],
                                            "vaz" => [0.0], "obs_alt" => 1000.0),
            "atmospheric_profile" => Dict{String,Any}(
                "T" => [285.0], "p" => [1012.99, 1013.0], "profile_reduction" => -1)))
        @test p.truncation isa δBGE                  # legacy default
        p.truncation = δBGE{Float64}(40, 5.0)        # user overrides
        resolved = vSmartMOM.CoreRT._resolved_truncation(p, Float64)
        @test resolved isa δBGE
        @test resolved.l_max == 40
        @test resolved.Δ_angle == 5.0
    end

    @testset "String parser whitelist (no eval)" begin
        # parameters_from_dict should match a small allow-list and
        # reject anything else with ArgumentError — no Meta.parse + eval.
        base_rt = Dict{String,Any}(
            "spec_bands" => ["[14492 14493]"],
            "surface" => ["LambertianSurfaceScalar(0.1)"],
            "quadrature_type" => "RadauQuad()",
            "polarization_type" => "Stokes_I()",
            "max_m" => 8, "Δ_angle" => 2.0, "l_trunc" => 20,
            "depol" => -1.0, "float_type" => "Float64", "architecture" => "CPU()")
        function build(trunc_spec)
            rt = copy(base_rt); rt["truncation"] = trunc_spec
            parameters_from_dict(Dict{String,Any}(
                "radiative_transfer" => rt,
                "geometry" => Dict{String,Any}("sza" => 30.0, "vza" => [0.0],
                                                "vaz" => [0.0], "obs_alt" => 1000.0),
                "atmospheric_profile" => Dict{String,Any}(
                    "T" => [285.0], "p" => [1012.99, 1013.0], "profile_reduction" => -1)))
        end
        @test build("NoTruncation()").truncation isa NoTruncation
        @test build("NoTruncation(l_max=32)").truncation.l_max == 32
        @test build("NoTruncation(32)").truncation.l_max == 32
        @test build("δBGE(20, 2.0)").truncation isa δBGE
        @test build("δBGE{Float64}(20, 2.0)").truncation.l_max == 20
        # Whitelist enforcement — refuses anything that isn't a
        # supported constructor shape, even valid Julia. Defends
        # against arbitrary-code execution from untrusted YAML.
        @test_throws ArgumentError build("println(\"pwned\")")
        @test_throws ArgumentError build("run(`whoami`)")
        @test_throws ArgumentError build("Foo()")
    end

    @testset "δ-M Eq. 8 invariants (delta_m_forward)" begin
        # Sanghavi & Stephens 2015 Eq. 8 is implemented downstream in
        # `delta_m_forward`, not inside `truncate_phase` (the
        # commented-out lines at truncate_phase.jl:115-116 / 260-261
        # would *double-apply* the rescaling and are correctly
        # disabled). Test the live function instead.
        using vSmartMOM.CoreRT: delta_m_forward
        Z⁺⁺ = randn(4, 4); Z⁻⁺ = randn(4, 4)

        # Invariant 1 — absorption budget τ(1−ω) is exactly invariant
        # under δ-M scaling, for any (τ, ω̃, fᵗ):
        #   τ_mod·(1 − ϖ_mod) = (1−fᵗω̃)τ · (1 − (1−fᵗ)ω̃/(1−fᵗω̃))
        #                     = τ·(1 − ω̃).
        for τ_aer in (0.1, 0.5, 1.0, 5.0),
            ω̃    in (0.1, 0.5, 0.9, 0.99),
            fᵗ   in (0.0, 0.1, 0.5, 0.9)
            out = delta_m_forward(τ_aer, ω̃, fᵗ, Z⁺⁺, Z⁻⁺)
            @test out.τ * (1 - out.ϖ) ≈ τ_aer * (1 - ω̃) rtol = 1e-12
        end

        # Invariant 2 — scattering optical depth shrinks by exactly
        # `(1 - fᵗ)`: τ_mod·ϖ_mod = τ·ω̃·(1−fᵗ).
        for τ_aer in (0.5, 2.0), ω̃ in (0.3, 0.95), fᵗ in (0.0, 0.3, 0.7)
            out = delta_m_forward(τ_aer, ω̃, fᵗ, Z⁺⁺, Z⁻⁺)
            @test out.τ * out.ϖ ≈ τ_aer * ω̃ * (1 - fᵗ) rtol = 1e-12
        end

        # Invariant 3 — fᵗ = 0 is the identity (NoTruncation limit).
        for τ_aer in (0.1, 1.0, 5.0), ω̃ in (0.05, 0.5, 0.99)
            out = delta_m_forward(τ_aer, ω̃, 0.0, Z⁺⁺, Z⁻⁺)
            @test out.τ ≈ τ_aer
            @test out.ϖ ≈ ω̃
        end

        # Invariant 4 — phase matrices pass through unchanged.
        out = delta_m_forward(0.5, 0.9, 0.3, Z⁺⁺, Z⁻⁺)
        @test out.Z⁺⁺ === Z⁺⁺
        @test out.Z⁻⁺ === Z⁻⁺
    end
end
