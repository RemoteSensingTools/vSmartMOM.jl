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
using Distributions
using Statistics

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

    @testset "δBGE renormalises all six Greek families (DoLP invariance)" begin
        # SS2015 Eq. 8: the δBGE phase matrix is
        #     Zᵗ = (Z − fᵗ·δ(cosΘ−1)·E) / (1 − fᵗ),      c₀ ≡ 1 − fᵗ,
        # so *every* Greek family is renormalised by 1/c₀. Only the diagonal
        # families α/β/δ/ζ additionally lose the forward-peak amount
        # (β_l − c_l), because the Dirac spike is diagonal:
        # γ_l^δ = ϵ_l^δ = 0 (S2014 Eqs. A.5–A.10; cf. SS2015 Eqs. 27c,d for
        # the δ-m analogue γⁿ_l = γ_l/(1−f_tr)).
        #
        # Observable consequence: since γᵗ builds f₁₂ and βᵗ builds f₁₁, and
        # both fits reproduce the untruncated functions outside the exclusion
        # cone, the degree of linear polarisation −f₁₂/f₁₁ must be *invariant*
        # under truncation. Dropping the 1/c₀ on γᵗ dilutes DoLP by exactly
        # c₀ — a ~29 % error for the coarse aerosol below, and ~30 % on the
        # TOA Stokes Q of a τ=0.5 coarse-dust layer.
        function coarse_optics(; μ_r = 1.5, σ_r = 1.3, r_max = 4.0, l_tr = 32,
                                 Δ_angle = 2.0, nquad_radius = 800)
            aero = Scattering.Aerosol(LogNormal(log(μ_r), log(σ_r)), 1.55, 0.005)
            mie = Scattering.MieModel(computation_type = Scattering.NAI2(),
                                      aerosol = aero, λ = 0.55,
                                      polarization_type = Stokes_IQUV(),
                                      truncation_type = δBGE{Float64}(l_tr, Δ_angle),
                                      r_max = r_max, nquad_radius = nquad_radius)
            raw = Scattering.compute_aerosol_optical_properties(mie, Float64)
            return raw, Scattering.truncate_phase(δBGE{Float64}(l_tr, Δ_angle), raw)
        end

        "DoLP = −f₁₂/f₁₁ reconstructed from Greek coefs on a θ grid [deg]."
        function dolp_of(greek, θ)
            sm = Scattering.reconstruct_phase(greek, cosd.(θ))
            return -sm.f₁₂ ./ sm.f₁₁
        end

        # ── coarse aerosol: genuine forward peak, fᵗ ≈ 0.29 ──────────────
        raw, tr = coarse_optics()
        c₀ = 1 - tr.fᵗ
        @test tr.fᵗ > 0.2                       # a real peak to remove

        # Sample well outside the 2° exclusion cone, away from the DoLP
        # sign changes where the ratio is numerically meaningless.
        θ = collect(range(10.0, 180.0, length = 341))
        d_raw = dolp_of(raw.greek_coefs, θ)
        d_tr  = dolp_of(tr.greek_coefs,  θ)
        keep  = abs.(d_raw) .> 0.02
        @test count(keep) > 100                 # the comparison has support
        ratio = d_tr[keep] ./ d_raw[keep]

        # The invariant: DoLP is preserved (residual is δBGE *fit* error at
        # l_tr = 32, not a systematic rescaling).
        @test median(ratio) ≈ 1 atol = 0.05
        @test median(abs.(d_tr[keep] .- d_raw[keep])) < 0.02
        # …and it is emphatically NOT diluted by c₀, which is what the
        # un-normalised γᵗ produced (median ratio ≈ 0.65 pre-fix).
        @test abs(median(ratio) - c₀) > 0.15

        # ── fine aerosol: fᵗ = 0 ⇒ c₀ = 1 ⇒ strict no-op ────────────────
        raw0, tr0 = coarse_optics(μ_r = 0.1, σ_r = 1.5, r_max = 1.0)
        @test tr0.fᵗ < 1e-6
        d0_raw = dolp_of(raw0.greek_coefs, θ)
        d0_tr  = dolp_of(tr0.greek_coefs,  θ)
        keep0  = abs.(d0_raw) .> 0.02
        # 1/c₀ with c₀ = 1 changes nothing: DoLP agreement at roundoff level.
        @test maximum(abs.(d0_tr[keep0] .- d0_raw[keep0])) < 1e-7
    end

    @testset "δBGE linearised path: γ̇ᵗ/ϵ̇ᵗ carry the 1/c₀ chain rule" begin
        # `truncate_phase(δBGE, aero, lin_aero)` must supply derivatives of the
        # *renormalised* coefficients, i.e. including the −xᵗ·ċ₀ term with
        # ċ₀ = ẋβ[:,1]. Verified by central finite differences of the very same
        # function, so the check is independent of the (pre-existing) fact that
        # the linearised fits sum over the full μ range while the two-argument
        # `truncate_phase` restricts them to the `Δ_angle` exclusion cone.
        # Derivative order is (nᵣ, nᵢ, rₚ, σₚ); we finite-difference nᵣ.
        l_tr, Δ_angle = 24, 2.0
        mk(nᵣ) = Scattering.MieModel(
            computation_type = Scattering.NAI2(),
            aerosol = Scattering.Aerosol(LogNormal(log(0.8), log(1.3)), nᵣ, 0.005),
            λ = 0.55, polarization_type = Stokes_IQUV(),
            truncation_type = δBGE{Float64}(l_tr, Δ_angle),
            r_max = 2.5, nquad_radius = 600)
        # truncated Greek coefs via the *linearised* method at refractive index nᵣ
        function trunc_lin_at(nᵣ)
            a, la = Scattering.compute_aerosol_optical_properties(
                LinMode(), mk(nᵣ), Float64)
            return Scattering.truncate_phase(δBGE{Float64}(l_tr, Δ_angle), a, la)
        end

        nᵣ, h = 1.55, 1e-5
        tr0, lin0 = trunc_lin_at(nᵣ)
        trp, _    = trunc_lin_at(nᵣ + h)
        trm, _    = trunc_lin_at(nᵣ - h)

        @test tr0.fᵗ > 0.01                       # non-trivial c₀ = 1 − fᵗ

        # β̇ᵗ is the family that already carried the full chain rule — it is the
        # control. γ̇ᵗ / ϵ̇ᵗ are the ones this fix repaired; pre-fix they were
        # off by both the 1/c₀ factor and the −xᵗ·ċ₀ term.
        for (f, fd) in ((:β, :β̇), (:γ, :γ̇), (:ϵ, :ϵ̇), (:α, :α̇), (:δ, :δ̇), (:ζ, :ζ̇))
            fd_ref   = (getfield(trp.greek_coefs, f) .- getfield(trm.greek_coefs, f)) ./ (2h)
            analytic = getfield(lin0.lin_greek_coefs, fd)[1, :]
            scale    = maximum(abs, fd_ref)
            @test maximum(abs, analytic .- fd_ref) < 1e-3 * scale
        end

        # ḟᵗ = −ċ₀ must be consistent with the same finite difference.
        @test lin0.ḟᵗ[1] ≈ (trp.fᵗ - trm.fᵗ) / (2h) rtol = 1e-3
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
