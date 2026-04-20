# Phase 1b regression gate — RRS forward model against a frozen reference
# =========================================================================
#
# Runs the toy 1-band RRS configuration from Phase1b_RRS_761-764nm.yaml
# (762–765 nm, Stokes_IQU, Float32, CPU) through the full forward path,
# then compares I/Q/U/V and ieR/ieT against the reference JLD2 committed
# alongside this file.
#
# The committed reference was generated on sanghavi-unified itself (i.e.,
# a self-reference, not a sanghavi-worktree cross-check). That means this
# test catches any numerical regression on sanghavi-unified after the
# Phase 1b merge is cemented — but it does NOT (yet) verify "matches
# sanghavi physics to the 6th decimal place." A follow-up commit will
# regenerate this reference on the sanghavi worktree so the comparison
# becomes the sanghavi-authoritative gate the merge plan actually calls for.
#
# Tolerances (from user's 2026-04-19 Phase 1b kickoff Q3 answers):
#  - I, ieI   :  atol = 1e-6, rtol = 1e-6  (6 decimal places, Float32)
#  - Q, U, V  :  atol = 1e-6, rtol = 0     (zero-valued pixels → atol only)
#  - ieQ/U/V  :  atol = 1e-6, rtol = 0
#  - wall-clock: merged ≤ 1.02 × reference wall (2% ceiling)
# Single-pixel fail-hard comparison (no percentile aggregation).
#
# Skipped automatically when CUDA is unavailable AND we don't have time
# for the ~3 min CPU run; pass `PHASE1B_CPU=1` env var to run on CPU anyway.

using vSmartMOM, vSmartMOM.CoreRT
using vSmartMOM.InelasticScattering
using Statistics
using JLD2
using Test

const PHASE1B_RUN_CPU = get(ENV, "PHASE1B_CPU", "0") == "1" || try
    using CUDA; CUDA.functional()
catch; false; end

@testset "Phase 1b RRS regression" begin

if !PHASE1B_RUN_CPU
    @test_skip "PHASE1B_CPU=1 not set and no CUDA — skipping 197s CPU run"
else

    params = parameters_from_yaml("test_parameters/Phase1b_RRS_761-764nm.yaml")
    model = model_from_parameters(params)

    iBand = 1
    FT = Float32
    ν = model.atmosphere.spec_bands[iBand]
    nSpec = length(ν)
    ν̃ = mean(ν)
    effT = (model.profile.vcd_dry' * model.profile.T) / sum(model.profile.vcd_dry)
    n2, o2 = InelasticScattering.getRamanAtmoConstants(ν̃, effT)
    nPol = CoreRT.polarization_type(model).n

    F₀ = zeros(FT, nPol, nSpec); F₀[1, :] .= 1
    SIF₀ = zeros(FT, nPol, nSpec)

    RS_type = InelasticScattering.RRS(
        n2 = n2, o2 = o2,
        greek_raman = InelasticScattering.GreekCoefs(
            [FT(1)], [FT(1)], [FT(1)], [FT(1)], [FT(1)], [FT(1)]),
        fscattRayl  = [FT(1)], ϖ_Cabannes = [FT(1)],
        ϖ_λ₁λ₀ = zeros(FT, 1), i_λ₁λ₀ = zeros(Int, 1),
        Z⁻⁺_λ₁λ₀ = zeros(FT, 1, 1), Z⁺⁺_λ₁λ₀ = zeros(FT, 1, 1),
        i_ref = argmin(abs.(ν .- ν̃)),
        n_Raman = 0,
        F₀ = F₀, SIF₀ = SIF₀)
    CoreRT.getRamanSSProp!(RS_type, 1e7/ν̃, ν)

    t0 = time()
    R_rrs, T_rrs, ieR, ieT = CoreRT.rt_run_test(RS_type, model, iBand)
    wall_merged = time() - t0

    ref = load(joinpath(@__DIR__, "reference", "phase1b_RRS_unified_selfref.jld2"))
    R_ref   = ref["R_rrs"]
    T_ref   = ref["T_rrs"]
    ieR_ref = ref["ieR"]
    ieT_ref = ref["ieT"]
    wall_ref = ref["wall"]

    atol_all = 1e-6
    rtol_I   = 1e-6

    # Stokes I on elastic and inelastic paths — relative tol enabled.
    @test all(isapprox.(R_rrs[:, 1, :],   R_ref[:, 1, :];   atol=atol_all, rtol=rtol_I))
    @test all(isapprox.(ieR[:, 1, :],     ieR_ref[:, 1, :]; atol=atol_all, rtol=rtol_I))
    @test all(isapprox.(T_rrs[:, 1, :],   T_ref[:, 1, :];   atol=atol_all, rtol=rtol_I))
    @test all(isapprox.(ieT[:, 1, :],     ieT_ref[:, 1, :]; atol=atol_all, rtol=rtol_I))

    # Stokes Q / U / (V if present) — absolute tol only; zeros allowed.
    for ipol in 2:nPol
        @test all(isapprox.(R_rrs[:, ipol, :], R_ref[:, ipol, :]; atol=atol_all, rtol=0))
        @test all(isapprox.(ieR[:, ipol, :],   ieR_ref[:, ipol, :]; atol=atol_all, rtol=0))
        @test all(isapprox.(T_rrs[:, ipol, :], T_ref[:, ipol, :]; atol=atol_all, rtol=0))
        @test all(isapprox.(ieT[:, ipol, :],   ieT_ref[:, ipol, :]; atol=atol_all, rtol=0))
    end

    # Wall-clock — 2 % ceiling relative to reference.
    @test wall_merged ≤ 1.02 * wall_ref
    @info "Phase1b wall-clock" merged=wall_merged ref=wall_ref ratio=wall_merged/wall_ref
end

end # @testset
