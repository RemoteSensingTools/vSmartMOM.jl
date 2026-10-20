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

# Rotational lines of `rs` with their strengths S = vmr σ (ν₀ + Δν)⁴.
function _phase1b_lines(rs, ν₀)
    lines = Tuple{Float64, Float64}[]
    for mol in (rs.n2, rs.o2)
        e = mol.effCoeff
        for (Δ, σ) in zip(vcat(parent(e.Δν̃_RoRaman_coeff_JtoJp2), parent(e.Δν̃_RoRaman_coeff_JtoJm2)),
                          vcat(parent(e.σ_RoRaman_coeff_JtoJp2), parent(e.σ_RoRaman_coeff_JtoJm2)))
            σ > 0 && push!(lines, (Δ, mol.vmr * σ * (ν₀ + Δ)^4))
        end
    end
    return lines
end

"""
Line weights by source offset (`i => w`, kernel convention n₀ = n₁ + i) on the
band grid `ν`, with shifts from the grid point nearest the band mean, for the
legacy placement (S/2 on each neighbouring grid point, as when the reference
was made; `:legacy`) or for linear interpolation (`:linear`). Built here
independently of the package code.
"""
function _phase1b_line_weights(rs, ν, placement)
    n = length(ν)
    k₀ = argmin(abs.(ν .- sum(ν) / n))
    ν₀ = Float64(ν[k₀])
    Δν_grid = Float64.(ν) .- ν₀
    w = Dict{Int, Float64}()
    add!(k, x) = x > 0 && (w[k₀ - k] = get(w, k₀ - k, 0.0) + x)   # i = -(output offset)
    for (Δ, S) in _phase1b_lines(rs, ν₀)
        first(Δν_grid) < Δ < last(Δν_grid) || continue
        if placement === :legacy
            k = argmin(abs.(Δ .- Δν_grid))
            pair = Δν_grid[k] < Δ ? (k, k + 1) : (k - 1, k)
            foreach(kk -> add!(kk, S / 2), pair)
        else
            k = searchsortedlast(Δν_grid, Δ)
            f = (Δ - Δν_grid[k]) / (Δν_grid[k + 1] - Δν_grid[k])
            add!(k, (1 - f) * S)
            add!(k + 1, f * S)
        end
    end
    return w
end

"""
Per-output-point factor that takes the frozen inelastic reference from the
legacy line placement to linear interpolation: the ratio of the fractions of
the line weight whose source point lies inside the band.
"""
function _phase1b_line_placement_correction(rs, ν)
    n = length(ν)
    legacy = _phase1b_line_weights(rs, ν, :legacy)
    linear = _phase1b_line_weights(rs, ν, :linear)
    in_band(w, k) = sum((wi for (i, wi) in w if 1 <= k + i <= n); init = 0.0) / sum(values(w))
    return [in_band(linear, k) / in_band(legacy, k) for k in 1:n]
end

@testset "Phase 1b RRS regression" begin

if !PHASE1B_RUN_CPU
    @test_skip "PHASE1B_CPU=1 not set and no CUDA — skipping 197s CPU run"
else

    params = parameters_from_yaml("test_parameters/Phase1b_RRS_761-764nm.yaml")
    # The reference predates the dispersive anisotropy (decision 13c), which
    # moves the depol: -1 τ_rayl at 763 nm by -0.3 %. Select the constant
    # anisotropy of the reference.
    params.rayleigh_anisotropy = ConstantAnisotropy()
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
    # Line-placement correction of the frozen inelastic reference (data left
    # untouched). The reference put half of each rotational line on each of
    # its two neighbouring grid points; lines are now shared by linear
    # interpolation (v2.2 decision log, decision 7, R7). The total weight is
    # unchanged, and with a flat F₀ the inelastic radiance at output point k is,
    # to first order (τ_Rayleigh varies by ≈ 1.6 % across the band),
    # proportional to the weight whose source k + i_λ₁λ₀ lies inside the band,
    # so the reference is rescaled by the ratio of those in-band fractions.
    # Away from the points where a line's source leaves the band the two
    # placements agree to ~1e-5; with the rescaling they agree to < 1e-3
    # everywhere (uncorrected: up to 12 % at six such points). The weights the
    # code uses are first checked against the independent linear placement.
    linear = _phase1b_line_weights(RS_type, ν, :linear)
    @test sort(collect(keys(linear))) == sort(RS_type.i_λ₁λ₀)
    w_code = Dict(zip(RS_type.i_λ₁λ₀, RS_type.ϖ_λ₁λ₀ ./ sum(RS_type.ϖ_λ₁λ₀)))
    @test all(isapprox(w_code[i], w / sum(values(linear)); rtol = 1e-5) for (i, w) in linear)
    corr = _phase1b_line_placement_correction(RS_type, ν)
    ieR_ref = ref["ieR"] .* reshape(corr, 1, 1, :)
    ieT_ref = ref["ieT"] .* reshape(corr, 1, 1, :)

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
