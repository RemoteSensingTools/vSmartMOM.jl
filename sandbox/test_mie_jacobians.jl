#!/usr/bin/env julia
# ==========================================================================
# Mie Microphysical Parameter Jacobian Test
# ==========================================================================
# Tests the Jacobians for parameters that go through the Mie derivative 
# interpolation chain (nᵣ, nᵢ, rₘ, σᵣ).
#
# These are the parameters that UNIQUELY exercise the Mie derivative code
# path in lin_model_from_parameters.jl (lines 498-511), which had:
#   Bug 20a: k̇sca used k̇*ω̃̇ instead of product rule k̇*ω̃ + k*ω̃̇
#   Bug 20b: ω̃̇ divided by k̇ instead of quotient rule (k̇sca - ω̃*k̇)/k
#   Bug 21:  τ̇_aer missing k_band factor in k_ref normalization term
#
# NOTE: Each perturbation requires a FULL Mie recalculation (~slow), 
# so we use few perturbation sizes but verify convergence.
# ==========================================================================

using vSmartMOM, vSmartMOM.CoreRT
using vSmartMOM.CoreRT: RT_Aerosol
using Distributions, Statistics

println("=" ^ 70)
println("  Mie Microphysical Parameter Jacobian Tests")
println("=" ^ 70)

# ---------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------
function run_lin_rt(params)
    model, lin_model = model_from_parameters(LinMode(), params)
    NAer = length(params.scattering_params.rt_aerosols)
    NGas = size(lin_model.τ̇_abs[1], 1)
    NSurf = 1
    R, T, dR, dT = rt_run(model, lin_model, NAer, NGas, NSurf)
    return R, dR, NAer, NGas, NSurf, model, lin_model
end

function run_fwd_only(params)
    # Use the same linearized path to ensure consistent code paths
    R, _, _, _, _, _, _ = run_lin_rt(params)
    return R
end

function rel_errors(analytic, fd; threshold=1e-12)
    mask = abs.(fd) .> threshold
    if any(mask)
        errs = abs.(analytic[mask] .- fd[mask]) ./ abs.(fd[mask])
        return maximum(errs), mean(errs), median(errs)
    else
        ae = abs.(analytic .- fd)
        return maximum(ae), mean(ae), median(ae)
    end
end

# =====================================================================
# Run base computation (once)
# =====================================================================
yaml_file = "test/test_parameters/JacobianTestFast.yaml"
println("\n[Setup] Running base linearized RT (includes Mie calculation)...")
println("  Config: $yaml_file")
flush(stdout)
params_base = parameters_from_yaml(yaml_file)
R_base, dR_base, NAer, NGas, NSurf, model_base, lin_model_base = run_lin_rt(params_base)
Nparams = NAer * 7 + NGas + NSurf

nᵣ_base = params_base.scattering_params.rt_aerosols[1].aerosol.nᵣ
nᵢ_base = params_base.scattering_params.rt_aerosols[1].aerosol.nᵢ
sd_base = params_base.scattering_params.rt_aerosols[1].aerosol.size_distribution
μ_base = sd_base.μ   # LogNormal log-mean
σ_base = sd_base.σ   # LogNormal log-std

println("  R shape: $(size(R_base)), dR shape: $(size(dR_base))")
println("  NAer=$NAer, NGas=$NGas, NSurf=$NSurf, Nparams=$Nparams")
println("  nᵣ=$nᵣ_base, nᵢ=$nᵢ_base")
println("  LogNormal: μ=$(round(μ_base, sigdigits=6)), σ=$(round(σ_base, sigdigits=6))")
println("  (median radius=$(round(exp(μ_base), sigdigits=4)) μm)")

# =====================================================================
# Level 0: Check Mie derivative interpolation output  
# =====================================================================
println("\n" * "=" ^ 70)
println("  Level 0: Mie Interpolation Output Check")
println("=" ^ 70)

aer_optics = model_base.aerosol_optics[1][1]
lin_aer_optics = lin_model_base.lin_aerosol_optics[1][1]

param_names_mie = ["nᵣ", "nᵢ", "rₘ(μ)", "σᵣ(σ)"]
println("\n  Forward: k ∈ [$(round(minimum(aer_optics.k),sigdigits=4)), $(round(maximum(aer_optics.k),sigdigits=4))], ω̃ ∈ [$(round(minimum(aer_optics.ω̃),sigdigits=4)), $(round(maximum(aer_optics.ω̃),sigdigits=4))]")

println("\n  Derivative ranges (from interpolation, post-fix):")
for ctr in 1:4
    k_min, k_max = extrema(lin_aer_optics.k̇[ctr,:])
    s_min, s_max = extrema(lin_aer_optics.ω̃̇[ctr,:])
    println("    $(param_names_mie[ctr]): k̇=($(round(k_min,sigdigits=4)), $(round(k_max,sigdigits=4))), ω̃̇=($(round(s_min,sigdigits=4)), $(round(s_max,sigdigits=4)))")
end

# =====================================================================
# FD Convergence: nᵣ (param 2)
# =====================================================================
println("\n" * "=" ^ 70)
println("  FD Convergence: nᵣ (param 2, real refractive index)")
println("=" ^ 70)
flush(stdout)

analytic_nr = dR_base[2, :, :, :]
println("  Analytic dR/dnᵣ range: ($(round(minimum(analytic_nr), sigdigits=4)), $(round(maximum(analytic_nr), sigdigits=4)))")

for δ in [1e-2, 3e-3, 1e-3, 3e-4]
    println("  Running FD with δ=$(δ)...")
    flush(stdout)
    params_pert = parameters_from_yaml(yaml_file)
    params_pert.scattering_params.rt_aerosols[1].aerosol.nᵣ = nᵣ_base + δ
    R_pert = run_fwd_only(params_pert)
    fd = (R_pert .- R_base) ./ δ
    max_e, mean_e, med_e = rel_errors(analytic_nr, fd)
    println("    δ=$δ: max=$(round(max_e, sigdigits=4)), mean=$(round(mean_e, sigdigits=4)), median=$(round(med_e, sigdigits=4))")
    flush(stdout)
end

# =====================================================================
# FD Convergence: nᵢ (param 3)  
# =====================================================================
println("\n" * "=" ^ 70)
println("  FD Convergence: nᵢ (param 3, imaginary refractive index)")
println("=" ^ 70)
flush(stdout)

analytic_ni = dR_base[3, :, :, :]
println("  Analytic dR/dnᵢ range: ($(round(minimum(analytic_ni), sigdigits=4)), $(round(maximum(analytic_ni), sigdigits=4)))")

# nᵢ is very small (1e-8), use absolute perturbations
for δ in [1e-4, 1e-5, 1e-6]
    println("  Running FD with δ=$(δ)...")
    flush(stdout)
    params_pert = parameters_from_yaml(yaml_file)
    params_pert.scattering_params.rt_aerosols[1].aerosol.nᵢ = nᵢ_base + δ
    R_pert = run_fwd_only(params_pert)
    fd = (R_pert .- R_base) ./ δ
    max_e, mean_e, med_e = rel_errors(analytic_ni, fd)
    println("    δ=$δ: max=$(round(max_e, sigdigits=4)), mean=$(round(mean_e, sigdigits=4)), median=$(round(med_e, sigdigits=4))")
    flush(stdout)
end

# =====================================================================
# FD Convergence: μ (param 4, log-mean of LogNormal)
# =====================================================================
println("\n" * "=" ^ 70)
println("  FD Convergence: μ (param 4, log-mean of size distribution)")
println("=" ^ 70)
flush(stdout)

analytic_mu = dR_base[4, :, :, :]
println("  Analytic dR/dμ range: ($(round(minimum(analytic_mu), sigdigits=4)), $(round(maximum(analytic_mu), sigdigits=4)))")

for δ_frac in [1e-1, 3e-2, 1e-2, 3e-3]
    δ = abs(μ_base) * δ_frac
    println("  Running FD with δ/|μ|=$(δ_frac) (δ=$(round(δ, sigdigits=4)))...")
    flush(stdout)
    params_pert = parameters_from_yaml(yaml_file)
    params_pert.scattering_params.rt_aerosols[1].aerosol.size_distribution = LogNormal(μ_base + δ, σ_base)
    R_pert = run_fwd_only(params_pert)
    fd = (R_pert .- R_base) ./ δ
    max_e, mean_e, med_e = rel_errors(analytic_mu, fd)
    println("    δ/|μ|=$δ_frac: max=$(round(max_e, sigdigits=4)), mean=$(round(mean_e, sigdigits=4)), median=$(round(med_e, sigdigits=4))")
    flush(stdout)
end

# =====================================================================
# FD Convergence: σ (param 5, log-std of LogNormal)
# =====================================================================
println("\n" * "=" ^ 70)
println("  FD Convergence: σ (param 5, log-std of size distribution)")
println("=" ^ 70)
flush(stdout)

analytic_sig = dR_base[5, :, :, :]
println("  Analytic dR/dσ range: ($(round(minimum(analytic_sig), sigdigits=4)), $(round(maximum(analytic_sig), sigdigits=4)))")

for δ_frac in [1e-1, 3e-2, 1e-2, 3e-3]
    δ = σ_base * δ_frac
    println("  Running FD with δ/σ=$(δ_frac) (δ=$(round(δ, sigdigits=4)))...")
    flush(stdout)
    params_pert = parameters_from_yaml(yaml_file)
    params_pert.scattering_params.rt_aerosols[1].aerosol.size_distribution = LogNormal(μ_base, σ_base + δ)
    R_pert = run_fwd_only(params_pert)
    fd = (R_pert .- R_base) ./ δ
    max_e, mean_e, med_e = rel_errors(analytic_sig, fd)
    println("    δ/σ=$δ_frac: max=$(round(max_e, sigdigits=4)), mean=$(round(mean_e, sigdigits=4)), median=$(round(med_e, sigdigits=4))")
    flush(stdout)
end

# =====================================================================
# Summary
# =====================================================================
println("\n" * "=" ^ 70)
println("  Summary: Mie Parameter Jacobian Test")
println("=" ^ 70)
println("""
  Parameter pathway analysis:
  - Surface albedo: bypasses Mie → tests adding/interaction only
  - τ_ref:          bypasses Mie → tests δ-M + RT kernel chain
  - p₀, σp:         bypasses Mie → tests profile + δ-M + RT kernel
  - nᵣ, nᵢ, μ, σ:   exercises Mie → tests the FULL chain including:
      1. Mie derivative computation (compute_NAI2_lin.jl)
      2. Mie derivative interpolation (lin_model_from_parameters.jl:498-511)
      3. τ̇_aer with k_ref normalization (lin_model_from_parameters.jl:562-568)
      4. δ-M scaling derivatives (createAero)
      5. Core optical property combination (+)
      6. RT kernel linearization

  Fixes applied:
  - Bug 20a: k̇sca product rule (k̇*ω̃ + k*ω̃̇ instead of k̇*ω̃̇)
  - Bug 20b: ω̃̇ quotient rule ((k̇sca - ω̃*k̇)/k instead of k̇sca/k̇)
  - Bug 21:  τ̇_aer k_band factor (.* k_band in k_ref term)
""")
println("=" ^ 70)
println("  All Mie Jacobian tests complete")
println("=" ^ 70)
