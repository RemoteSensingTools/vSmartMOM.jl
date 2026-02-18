#!/usr/bin/env julia
# ==========================================================================
# Focused FD Convergence Test for τ_ref and p₀
# ==========================================================================
# These only perturb τ_ref / p₀ which don't need Mie recalculation,
# so each FD step is fast (only needs rt_run, not the full Mie calc).
# Also includes Level 4: Mie ω̃̇ interpolation sanity check.
# ==========================================================================

using vSmartMOM, vSmartMOM.CoreRT
using vSmartMOM.CoreRT: RT_Aerosol
using Distributions, Statistics

println("=" ^ 70)
println("  FD Convergence Test + Mie Interpolation Diagnostic")
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
    R, dR, _, _, _, _, _ = run_lin_rt(params)
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
println("\n[Setup] Running base linearized RT...")
params_base = parameters_from_yaml("test/test_parameters/JacobianTest.yaml")
R_base, dR_base, NAer, NGas, NSurf, model_base, lin_model_base = run_lin_rt(params_base)
Nparams = NAer * 7 + NGas + NSurf

τ_ref_base = params_base.scattering_params.rt_aerosols[1].τ_ref
p0_base = mean(params_base.scattering_params.rt_aerosols[1].profile)
σp_base = std(params_base.scattering_params.rt_aerosols[1].profile)

println("  R shape: $(size(R_base)), dR shape: $(size(dR_base))")
println("  NAer=$NAer, NGas=$NGas, NSurf=$NSurf, Nparams=$Nparams")
println("  τ_ref=$τ_ref_base, p₀=$p0_base, σp=$σp_base")

# =====================================================================
# Level 4: Mie ω̃̇ Interpolation Diagnostic (no FD needed)
# =====================================================================
println("\n" * "=" ^ 70)
println("  Mie ω̃̇ Interpolation Diagnostic")
println("=" ^ 70)

aer_optics = model_base.aerosol_optics[1][1]
lin_aer_optics = lin_model_base.lin_aerosol_optics[1][1]

param_names = ["nᵣ", "nᵢ", "rₘ", "σᵣ"]

println("\n  Forward Mie properties:")
println("    k range: $(round(minimum(aer_optics.k),sigdigits=4)) to $(round(maximum(aer_optics.k),sigdigits=4))")
println("    ω̃ range: $(round(minimum(aer_optics.ω̃),sigdigits=4)) to $(round(maximum(aer_optics.ω̃),sigdigits=4))")
println("    fᵗ: $(round(aer_optics.fᵗ, sigdigits=6))")

println("\n  Mie derivative ranges:")
for ctr in 1:4
    k_min, k_max = extrema(lin_aer_optics.k̇[ctr,:])
    s_min, s_max = extrema(lin_aer_optics.ω̃̇[ctr,:])
    println("    $(param_names[ctr]): k̇=($(round(k_min,sigdigits=4)), $(round(k_max,sigdigits=4))), ω̃̇=($(round(s_min,sigdigits=4)), $(round(s_max,sigdigits=4)))")
end

# Check the k̇sca interpolation formula
# Code uses: k̇sca_grid = k̇_ext * ω̃̇ (should be k̇_ext * ω̃ + k_ext * ω̃̇)
println("\n  k̇sca formula check:")
println("  (If ratio ≈ 1.0 with small std, the interpolation is self-consistent)")
println("  (If ratio ≠ 1.0, there's a formula error in lin_model_from_parameters.jl:503)")
for ctr in 1:4
    k_dot = lin_aer_optics.k̇[ctr,:]
    ssa_dot = lin_aer_optics.ω̃̇[ctr,:]
    # Correct k̇_sca = k̇*ω̃ + k*ω̃̇  (product rule for k_sca = k*ω̃)
    ksca_dot_correct = k_dot .* aer_optics.ω̃ .+ aer_optics.k .* ssa_dot
    # What the code uses for interpolation: k̇*ω̃̇
    ksca_dot_code = k_dot .* ssa_dot
    mask = abs.(ksca_dot_code) .> 1e-20
    if any(mask)
        ratio = ksca_dot_correct[mask] ./ ksca_dot_code[mask]
        println("    $(param_names[ctr]): correct/code ratio: mean=$(round(mean(ratio), sigdigits=6)), std=$(round(std(ratio), sigdigits=4)), range=($(round(minimum(ratio),sigdigits=4)), $(round(maximum(ratio),sigdigits=4)))")
    else
        println("    $(param_names[ctr]): k̇sca_code ≈ 0")
    end
end

# =====================================================================
# FD Convergence: Surface Albedo
# =====================================================================
println("\n" * "=" ^ 70)
println("  FD Convergence: Surface Albedo (param $Nparams)")
println("=" ^ 70)

albedo_base = 0.05
analytic_albedo = dR_base[Nparams, :, :, :]

for δ in [1e-2, 1e-3, 1e-4, 1e-5, 1e-6]
    params_pert = parameters_from_yaml("test/test_parameters/JacobianTest.yaml")
    params_pert.brdf = [LambertianSurfaceScalar(albedo_base + δ)]
    R_pert = run_fwd_only(params_pert)
    fd = (R_pert .- R_base) ./ δ
    max_e, mean_e, med_e = rel_errors(analytic_albedo, fd)
    println("  δ=$δ: max=$(round(max_e, sigdigits=4)), mean=$(round(mean_e, sigdigits=4)), median=$(round(med_e, sigdigits=4))")
end

# =====================================================================
# FD Convergence: τ_ref  
# =====================================================================
println("\n" * "=" ^ 70)
println("  FD Convergence: τ_ref (param 1)")
println("=" ^ 70)

analytic_tau = dR_base[1, :, :, :]
println("  Analytic range: ($(minimum(analytic_tau)), $(maximum(analytic_tau)))")

for δ_frac in [1e-1, 3e-2, 1e-2, 3e-3, 1e-3, 3e-4, 1e-4, 3e-5, 1e-5]
    δ = τ_ref_base * δ_frac
    params_pert = parameters_from_yaml("test/test_parameters/JacobianTest.yaml")
    params_pert.scattering_params.rt_aerosols[1].τ_ref = τ_ref_base + δ
    R_pert = run_fwd_only(params_pert)
    fd = (R_pert .- R_base) ./ δ
    max_e, mean_e, med_e = rel_errors(analytic_tau, fd)
    println("  δ/τ=$(δ_frac): max=$(round(max_e, sigdigits=4)), mean=$(round(mean_e, sigdigits=4)), median=$(round(med_e, sigdigits=4))")
end

# =====================================================================
# FD Convergence: p₀
# =====================================================================
println("\n" * "=" ^ 70)
println("  FD Convergence: p₀ (param 6)")
println("=" ^ 70)

analytic_p0 = dR_base[6, :, :, :]
println("  Analytic range: ($(minimum(analytic_p0)), $(maximum(analytic_p0)))")

for δ_frac in [1e-1, 3e-2, 1e-2, 3e-3, 1e-3, 3e-4, 1e-4, 3e-5, 1e-5]
    δ = p0_base * δ_frac
    params_pert = parameters_from_yaml("test/test_parameters/JacobianTest.yaml")
    params_pert.scattering_params.rt_aerosols[1].profile = Normal(p0_base + δ, σp_base)
    R_pert = run_fwd_only(params_pert)
    fd = (R_pert .- R_base) ./ δ
    max_e, mean_e, med_e = rel_errors(analytic_p0, fd)
    println("  δ/p₀=$(δ_frac): max=$(round(max_e, sigdigits=4)), mean=$(round(mean_e, sigdigits=4)), median=$(round(med_e, sigdigits=4))")
end

# =====================================================================
# FD Convergence: σp
# =====================================================================
println("\n" * "=" ^ 70)
println("  FD Convergence: σp (param 7)")
println("=" ^ 70)

analytic_sp = dR_base[7, :, :, :]
println("  Analytic range: ($(minimum(analytic_sp)), $(maximum(analytic_sp)))")

for δ_frac in [1e-1, 3e-2, 1e-2, 3e-3, 1e-3, 3e-4, 1e-4]
    δ = σp_base * δ_frac
    params_pert = parameters_from_yaml("test/test_parameters/JacobianTest.yaml")
    params_pert.scattering_params.rt_aerosols[1].profile = Normal(p0_base, σp_base + δ)
    R_pert = run_fwd_only(params_pert)
    fd = (R_pert .- R_base) ./ δ
    max_e, mean_e, med_e = rel_errors(analytic_sp, fd)
    println("  δ/σp=$(δ_frac): max=$(round(max_e, sigdigits=4)), mean=$(round(mean_e, sigdigits=4)), median=$(round(med_e, sigdigits=4))")
end

println("\n" * "=" ^ 70)
println("  All convergence tests complete")
println("=" ^ 70)
