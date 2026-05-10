# Phase 1b regression gate — RRS forward model vs sanghavi reference
# =========================================================================
#
# Runs the toy 1-band RRS configuration from Phase1b_RRS_761-764nm.yaml
# (762–765 nm, Stokes_IQU, Float32, CPU) through the full forward path,
# then compares I/Q/U and ieR/ieT against the sanghavi-worktree JLD2 reference
# at test/reference/phase1b_RRS_sanghavi_q0.jld2. This is the "matches
# sanghavi physics" gate the merge plan calls for.
#
# Tolerances (2026-04-22 residual closed; greek_cabannes fix in
# compEffectiveLayerProperties.jl):
#  - I, Q, ieI, ieQ  :  atol = 1e-6, rtol = 0.02
#      Mean ratios match sanghavi to <1e-4 after fix; per-pixel max rel-error is
#      ~0.04% on R, ~1.6% on T (single-pixel Float32 accumulation noise).
#      rtol=0.02 gives headroom while still catching any real regression.
#  - T, ieT          :  atol = 1e-6, rtol = 0.02
#  - U, V, ieU, ieV  :  atol = 1e-6, rtol = 0 (zero for this config)
#  - wall-clock gate removed — comparison against sanghavi reference, hardware
#    asymmetry makes a wall ceiling meaningless here.
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

    R_rrs, T_rrs, ieR, ieT = CoreRT.rt_run_test(RS_type, model, iBand)

    ref = load(joinpath(@__DIR__, "reference", "phase1b_RRS_sanghavi_q0.jld2"))
    R_ref   = ref["R_rrs"]
    T_ref   = ref["T_rrs"]
    ieR_ref = ref["ieR"]
    ieT_ref = ref["ieT"]

    atol_all = 1e-6
    rtol_Stokes = 0.02   # see docstring — headroom for FP32 per-pixel noise after greek_cabannes fix
    rtol_U_V    = 0      # U and V expected zero for this config; atol enforced only

    # Stokes I/Q on elastic + inelastic — rtol enabled.
    for ipol in 1:min(2, nPol)
        @test all(isapprox.(R_rrs[:, ipol, :],   R_ref[:, ipol, :];   atol=atol_all, rtol=rtol_Stokes))
        @test all(isapprox.(ieR[:, ipol, :],     ieR_ref[:, ipol, :]; atol=atol_all, rtol=rtol_Stokes))
        @test all(isapprox.(T_rrs[:, ipol, :],   T_ref[:, ipol, :];   atol=atol_all, rtol=rtol_Stokes))
        @test all(isapprox.(ieT[:, ipol, :],     ieT_ref[:, ipol, :]; atol=atol_all, rtol=rtol_Stokes))
    end

    # Stokes U / (V if present) — zero-valued, atol-only check.
    for ipol in 3:nPol
        @test all(isapprox.(R_rrs[:, ipol, :], R_ref[:, ipol, :]; atol=atol_all, rtol=rtol_U_V))
        @test all(isapprox.(ieR[:, ipol, :],   ieR_ref[:, ipol, :]; atol=atol_all, rtol=rtol_U_V))
        @test all(isapprox.(T_rrs[:, ipol, :], T_ref[:, ipol, :]; atol=atol_all, rtol=rtol_U_V))
        @test all(isapprox.(ieT[:, ipol, :],   ieT_ref[:, ipol, :]; atol=atol_all, rtol=rtol_U_V))
    end
end

end # @testset
