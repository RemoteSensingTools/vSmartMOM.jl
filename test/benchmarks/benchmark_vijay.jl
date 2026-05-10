# =============================================================================
# benchmark_vijay.jl — Standalone Natraj (Vijay, 1980) polarized-Rayleigh
# benchmark driver.
#
# Replays the `@testset "compare against natraj paper"` block from
# test/test_CoreRT.jl (lines 41-85) but as a plain script so it can be run
# outside the Pkg.test() harness and produce human-readable diagnostics.
#
# Geometry (from natraj.yaml):
#   - τ_Rayleigh = 0.5, depol = 0, Lambertian α = 0
#   - λ = 360 nm, Stokes_IQUV, RadauQuad, max_m = 3, l_trunc = 20
#   - cos(SZA) = 0.2  (SZA ≈ 78.463°)
#   - 16 view cosines μ ∈ {0.02, 0.06, …, 1.00}  (VZA = acosd(μ))
#   - 7 relative azimuths φ ∈ {0°, 30°, …, 180°}
#
# Truth tables I_trues, Q_trues, U_trues (16 μ × 7 φ) come from
# Natraj et al., J. Quant. Spectrosc. Radiat. Transf. (1980). They live in
# test/benchmarks/natraj_trues.jl.
#
# Run:
#   julia --project=test test/benchmarks/benchmark_vijay.jl
# =============================================================================

using vSmartMOM
using vSmartMOM.CoreRT
using Printf
using Statistics

# Use pkgdir so the path is correct whether the script is run directly
# (`julia ... benchmark_vijay.jl`) or pasted/included from a REPL started in
# the repo root.
const BENCH_DIR = joinpath(pkgdir(vSmartMOM), "test", "benchmarks")
const YAML_PATH   = joinpath(BENCH_DIR, "natraj.yaml")
const TRUES_PATH  = joinpath(BENCH_DIR, "natraj_trues.jl")

# Load reference tables into module scope (natraj_trues.jl just assigns
# I_trues / Q_trues / U_trues at top level).
include(TRUES_PATH)

# Test tolerances from test_CoreRT.jl.
const ϵ_I = 0.002
const ϵ_QU = 0.008
const MIN_REL_Q_U = 0.01   # skip |Q|, |U| below this in rel-error

function run_natraj()
    μ   = [0.02, 0.06, 0.10, 0.16, 0.20, 0.28, 0.32, 0.40,
           0.52, 0.64, 0.72, 0.84, 0.92, 0.96, 0.98, 1.00]
    ϕs  = collect(0.0:30.0:180.0)
    τ   = 0.5
    λnm = 360.0

    I_mod = zeros(length(ϕs), length(μ))
    Q_mod = zeros(length(ϕs), length(μ))
    U_mod = zeros(length(ϕs), length(μ))

    # Natraj reference tables (loaded at module scope above).
    I_true = I_trues
    Q_true = Q_trues
    U_true = U_trues

    @info "Loading YAML" YAML_PATH
    params = parameters_from_yaml(YAML_PATH)
    params.spec_bands = [[1e7 / λnm, 1e7 / λnm + 1]]
    params.vza        = acosd.(μ)
    params.sza        = acosd(0.2)

    # Natraj 1980 uses the unphysical idealization depol = 0. natraj.yaml sets
    # `depol: 0.0`; with the 2026-04-24 fix to model_from_parameters.jl this is
    # honoured directly (any params.depol ≥ 0 is taken as the explicit depol;
    # -1 is the auto-from-molecular-constants sentinel).

    t0 = time()
    for (iϕ, ϕ) in enumerate(ϕs)
        params.vaz = repeat([ϕ], length(μ))
        model = model_from_parameters(params)
        model.τ_rayl[1] .= τ
        # rt_run returns radiance factor L = I/F₀; Natraj table is R = π·L.
        R = π .* CoreRT.rt_run(model, i_band = 1)[1]
        I_mod[iϕ, :] = R[:, 1, 1]
        Q_mod[iϕ, :] = R[:, 2, 1]
        U_mod[iϕ, :] = R[:, 3, 1]
        @info "Fourier-harmonic run" iϕ length(ϕs) ϕ elapsed_s = round(time() - t0; digits = 2)
    end

    ΔI = abs.(I_true .- I_mod') ./ I_true                                  # note transpose: trues are (μ × ϕ)
    ΔQ = abs.(Q_true .- Q_mod') ./ Q_true
    ΔU = abs.(U_true .- U_mod') ./ U_true

    # Restrict Q, U rel-err to cells where |model| ≥ 0.01, the convention in
    # the original testset.
    q_mask = abs.(Q_mod') .>= MIN_REL_Q_U
    u_mask = abs.(U_mod') .>= MIN_REL_Q_U

    ΔI_max = maximum(ΔI)
    ΔQ_max = any(q_mask) ? maximum(ΔQ[q_mask]) : NaN
    ΔU_valid = filter(!isnan, ΔU[u_mask])
    ΔU_max = isempty(ΔU_valid) ? NaN : maximum(ΔU_valid)

    println("\n=== Vijay (Natraj 1980) benchmark summary ===")
    @printf("τ_Rayleigh = %.2f,  cos(SZA) = 0.2,  λ = %.1f nm\n", τ, λnm)
    @printf("%-8s %-14s %-9s %-9s\n", "stokes", "max rel err", "tol", "pass?")
    @printf("%-8s %-14.4e %-9.3e %-9s\n", "I",  ΔI_max, ϵ_I,
            ΔI_max < ϵ_I  ? "PASS" : "FAIL")
    @printf("%-8s %-14.4e %-9.3e %-9s\n", "Q",  ΔQ_max, ϵ_QU,
            ΔQ_max < ϵ_QU ? "PASS" : "FAIL")
    @printf("%-8s %-14.4e %-9.3e %-9s\n", "U",  ΔU_max, ϵ_QU,
            ΔU_max < ϵ_QU ? "PASS" : "FAIL")

    println("\n--- per-(μ, ϕ) rel error (Stokes I) ---")
    print("         ")
    for ϕ in ϕs
        @printf("ϕ=%5.1f   ", ϕ)
    end
    println()
    for (iμ, μi) in enumerate(μ)
        @printf("μ=%5.2f  ", μi)
        for iϕ in eachindex(ϕs)
            @printf("%-10.3e", ΔI[iμ, iϕ])
        end
        println()
    end

    return (; I_mod, Q_mod, U_mod, I_true, Q_true, U_true,
              ΔI_max, ΔQ_max, ΔU_max,
              passed_I = ΔI_max < ϵ_I,
              passed_Q = ΔQ_max < ϵ_QU,
              passed_U = ΔU_max < ϵ_QU)
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_natraj()
end
