#!/usr/bin/env julia
# ==========================================================================
# Comprehensive Jacobian Diagnostic Test
# ==========================================================================
# Tests the derivative chain at EACH level separately:
#   Level 1: τ_aer derivatives (model construction, before Mie)
#   Level 2: createAero δ-M scaling derivatives (Mie chain)
#   Level 3: Full RT Jacobian with FD convergence study
#   Level 4: Mie ω̃̇ validation
#
# This isolates whether errors come from:
#   (a) Mie derivatives (ForwardDiff + interpolation)
#   (b) δ-M scaling assembly (createAero)
#   (c) RT propagation (elemental/doubling/interaction)
# ==========================================================================

using vSmartMOM, vSmartMOM.CoreRT
using vSmartMOM.CoreRT: RT_Aerosol
using Distributions, Statistics

println("=" ^ 70)
println("  Jacobian Diagnostic Test Suite")
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
    return R, dR, NAer, NGas, NSurf
end

function run_fwd_only(params)
    R, dR, _, _, _ = run_lin_rt(params)
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
# LEVEL 1: τ_aer Derivative Validation (before RT, before δ-M)
# =====================================================================
println("\n" * "=" ^ 70)
println("  LEVEL 1: τ_aer Derivative (Model Construction)")
println("=" ^ 70)

params_base = parameters_from_yaml("test/test_parameters/JacobianTest.yaml")
model_base, lin_model_base = model_from_parameters(LinMode(), params_base)

τ_aer_base = model_base.τ_aer[1]   # [nAer, nSpec, nZ]
τ̇_aer_base = lin_model_base.τ̇_aer[1]   # [nAer, 7, nSpec, nZ]

τ_ref_base = params_base.scattering_params.rt_aerosols[1].τ_ref
println("  τ_ref = $τ_ref_base")
println("  τ_aer shape: $(size(τ_aer_base))")
println("  τ̇_aer shape: $(size(τ̇_aer_base))")

# --- Test 1a: ∂τ_aer/∂τ_ref ---
println("\n--- Test 1a: ∂τ_aer/∂τ_ref ---")

# Analytic: τ̇_aer[1, 1, :, :] = ∂τ_aer/∂τ_ref
analytic_dtau_dtauref = τ̇_aer_base[1, 1, :, :]

# Theory check: ∂τ_aer/∂τ_ref should = τ_aer / τ_ref
expected_dtau_dtauref = τ_aer_base[1, :, :] ./ τ_ref_base
max_e, mean_e, _ = rel_errors(analytic_dtau_dtauref, expected_dtau_dtauref)
println("  τ̇_aer[1] vs τ_aer/τ_ref: max_rel=$(round(max_e, sigdigits=4)), mean_rel=$(round(mean_e, sigdigits=4))")

# FD check with multiple δ values
for δ_frac in [1e-2, 1e-3, 1e-4, 1e-5, 1e-6]
    δ = τ_ref_base * δ_frac
    params_pert = parameters_from_yaml("test/test_parameters/JacobianTest.yaml")
    params_pert.scattering_params.rt_aerosols[1].τ_ref = τ_ref_base + δ
    model_pert, _ = model_from_parameters(LinMode(), params_pert)
    τ_aer_pert = model_pert.τ_aer[1]
    fd = (τ_aer_pert[1, :, :] .- τ_aer_base[1, :, :]) ./ δ
    max_e, mean_e, _ = rel_errors(analytic_dtau_dtauref, fd)
    println("  FD δ=$(δ_frac): max_rel=$(round(max_e, sigdigits=4)), mean_rel=$(round(mean_e, sigdigits=4))")
end

# --- Test 1b: ∂τ_aer/∂p₀ ---
println("\n--- Test 1b: ∂τ_aer/∂p₀ ---")
analytic_dtau_dp0 = τ̇_aer_base[1, 6, :, :]
p0_base = mean(params_base.scattering_params.rt_aerosols[1].profile)
σp_base = std(params_base.scattering_params.rt_aerosols[1].profile)
println("  p₀ = $p0_base, σp = $σp_base")

for δ_frac in [1e-2, 1e-3, 1e-4, 1e-5, 1e-6]
    δ = p0_base * δ_frac
    params_pert = parameters_from_yaml("test/test_parameters/JacobianTest.yaml")
    params_pert.scattering_params.rt_aerosols[1].profile = Normal(p0_base + δ, σp_base)
    model_pert, _ = model_from_parameters(LinMode(), params_pert)
    τ_aer_pert = model_pert.τ_aer[1]
    fd = (τ_aer_pert[1, :, :] .- τ_aer_base[1, :, :]) ./ δ
    max_e, mean_e, _ = rel_errors(analytic_dtau_dp0, fd)
    println("  FD δ=$(δ_frac): max_rel=$(round(max_e, sigdigits=4)), mean_rel=$(round(mean_e, sigdigits=4))")
end

# =====================================================================
# LEVEL 2: Mie Derivative Validation (ω̃̇, ḟᵗ, k̇)  
# =====================================================================
println("\n" * "=" ^ 70)
println("  LEVEL 2: Mie Derivative Validation (ω̃, fᵗ, k)")
println("=" ^ 70)

# Extract Mie properties from base model
aer_optics = model_base.aerosol_optics[1][1]
lin_aer_optics = lin_model_base.lin_aerosol_optics[1][1]

println("  k (extinction) shape: $(size(aer_optics.k))")
println("  ω̃ (SSA) shape: $(size(aer_optics.ω̃))")
println("  k̇ shape: $(size(lin_aer_optics.k̇))")
println("  ω̃̇ shape: $(size(lin_aer_optics.ω̃̇))")
println("  ḟᵗ shape: $(size(lin_aer_optics.ḟᵗ))")

# Print Mie derivative ranges
println("\n  Mie derivative ranges:")
param_names = ["nᵣ", "nᵢ", "rₘ", "σᵣ"]
for ctr in 1:4
    k_min, k_max = extrema(lin_aer_optics.k̇[ctr,:])
    s_min, s_max = extrema(lin_aer_optics.ω̃̇[ctr,:])
    println("    $(param_names[ctr]): k̇=($(round(k_min,sigdigits=4)), $(round(k_max,sigdigits=4))), ω̃̇=($(round(s_min,sigdigits=4)), $(round(s_max,sigdigits=4)))")
end

# FD test: perturb nᵣ and check k, ω̃ derivatives
println("\n--- Test 2a: FD check of k, ω̃ vs nᵣ perturbation ---")
aer = params_base.scattering_params.rt_aerosols[1]
nr_base = aer.aerosol.nᵣ
println("  Base nᵣ = $nr_base")

for δ_frac in [1e-3, 1e-4, 1e-5, 1e-6]
    δ = nr_base * δ_frac
    params_pert = parameters_from_yaml("test/test_parameters/JacobianTest.yaml")
    params_pert.scattering_params.rt_aerosols[1].aerosol = 
        vSmartMOM.Scattering.Aerosol(
            params_pert.scattering_params.rt_aerosols[1].aerosol.size_distribution,
            nr_base + δ,
            params_pert.scattering_params.rt_aerosols[1].aerosol.nᵢ)
    model_pert, lin_model_pert = model_from_parameters(LinMode(), params_pert)
    
    aer_optics_pert = model_pert.aerosol_optics[1][1]
    
    fd_k = (aer_optics_pert.k .- aer_optics.k) ./ δ
    fd_ssa = (aer_optics_pert.ω̃ .- aer_optics.ω̃) ./ δ
    
    # nᵣ is Mie param 1 → index 1 in k̇, ω̃̇
    max_ek, mean_ek, _ = rel_errors(lin_aer_optics.k̇[1,:], fd_k)
    max_essa, mean_essa, _ = rel_errors(lin_aer_optics.ω̃̇[1,:], fd_ssa)
    println("  δ=$(δ_frac): k̇ max_rel=$(round(max_ek, sigdigits=4)), ω̃̇ max_rel=$(round(max_essa, sigdigits=4))")
end

# =====================================================================
# LEVEL 3: Full RT Jacobian FD Convergence Study
# =====================================================================
println("\n" * "=" ^ 70)
println("  LEVEL 3: Full RT Jacobian - FD Convergence Study")
println("=" ^ 70)

# Run base linearized RT (once)
println("\n  Running base linearized RT...")
R_base, dR_base, NAer, NGas, NSurf = run_lin_rt(params_base)
Nparams = NAer * 7 + NGas + NSurf
println("  R shape: $(size(R_base)), dR shape: $(size(dR_base)), Nparams=$Nparams")

# --- 3a: Surface albedo convergence ---
println("\n--- Test 3a: Surface Albedo Convergence ---")
albedo_base = 0.05
isurf = Nparams
analytic_albedo = dR_base[isurf, :, :, :]

for δ in [1e-2, 1e-3, 1e-4, 1e-5, 1e-6]
    params_pert = parameters_from_yaml("test/test_parameters/JacobianTest.yaml")
    params_pert.brdf = [LambertianSurfaceScalar(albedo_base + δ)]
    R_pert = run_fwd_only(params_pert)
    fd = (R_pert .- R_base) ./ δ
    max_e, mean_e, _ = rel_errors(analytic_albedo, fd)
    println("  δ=$δ: max_rel=$(round(max_e, sigdigits=4)), mean_rel=$(round(mean_e, sigdigits=4))")
end

# --- 3b: τ_ref convergence ---
println("\n--- Test 3b: τ_ref Convergence ---")
itau = 1
analytic_tau = dR_base[itau, :, :, :]

for δ_frac in [1e-1, 1e-2, 1e-3, 1e-4, 1e-5]
    δ = τ_ref_base * δ_frac
    params_pert = parameters_from_yaml("test/test_parameters/JacobianTest.yaml")
    params_pert.scattering_params.rt_aerosols[1].τ_ref = τ_ref_base + δ
    R_pert = run_fwd_only(params_pert)
    fd = (R_pert .- R_base) ./ δ
    max_e, mean_e, med_e = rel_errors(analytic_tau, fd)
    println("  δ=$(δ_frac): max_rel=$(round(max_e, sigdigits=4)), mean_rel=$(round(mean_e, sigdigits=4)), median_rel=$(round(med_e, sigdigits=4))")
end

# --- 3c: p₀ convergence ---
println("\n--- Test 3c: p₀ Convergence ---")
ip0 = 6
analytic_p0 = dR_base[ip0, :, :, :]

for δ_frac in [1e-1, 1e-2, 1e-3, 1e-4, 1e-5]
    δ = p0_base * δ_frac
    params_pert = parameters_from_yaml("test/test_parameters/JacobianTest.yaml")
    params_pert.scattering_params.rt_aerosols[1].profile = Normal(p0_base + δ, σp_base)
    R_pert = run_fwd_only(params_pert)
    fd = (R_pert .- R_base) ./ δ
    max_e, mean_e, med_e = rel_errors(analytic_p0, fd)
    println("  δ=$(δ_frac): max_rel=$(round(max_e, sigdigits=4)), mean_rel=$(round(mean_e, sigdigits=4)), median_rel=$(round(med_e, sigdigits=4))")
end

# --- 3d: nᵣ convergence (Mie parameter, full RT chain) ---
println("\n--- Test 3d: nᵣ Convergence (Full RT Chain) ---")
inr = 2  # nᵣ is parameter 2 (after τ_ref)
analytic_nr = dR_base[inr, :, :, :]

for δ_frac in [1e-2, 1e-3, 1e-4, 1e-5]
    δ = nr_base * δ_frac
    params_pert = parameters_from_yaml("test/test_parameters/JacobianTest.yaml")
    params_pert.scattering_params.rt_aerosols[1].aerosol = 
        vSmartMOM.Scattering.Aerosol(
            params_pert.scattering_params.rt_aerosols[1].aerosol.size_distribution,
            nr_base + δ,
            params_pert.scattering_params.rt_aerosols[1].aerosol.nᵢ)
    R_pert = run_fwd_only(params_pert)
    fd = (R_pert .- R_base) ./ δ
    max_e, mean_e, med_e = rel_errors(analytic_nr, fd)
    println("  δ=$(δ_frac): max_rel=$(round(max_e, sigdigits=4)), mean_rel=$(round(mean_e, sigdigits=4)), median_rel=$(round(med_e, sigdigits=4))")
end

# =====================================================================
# LEVEL 4: Mie ω̃̇ Interpolation Sanity Check
# =====================================================================
println("\n" * "=" ^ 70)
println("  LEVEL 4: Mie ω̃̇ Interpolation Diagnostic")
println("=" ^ 70)

# Check if ω̃̇ interpolation formula is correct
# The code uses: k̇sca = k̇_ext * ω̃̇ (at grid points), then ω̃̇(λ) = k̇sca(λ)/k̇(λ)
# Correct formula: k̇sca = k̇_ext * ω̃ + k_ext * ω̃̇, then ω̃̇(λ) = (k̇sca(λ) - ω̃(λ)*k̇(λ))/k(λ)
# Let's check if the current values match what we'd expect

println("  ω̃ range: $(extrema(aer_optics.ω̃))")
println("  k range: $(extrema(aer_optics.k))")
println("  k̇[1] (∂k/∂nᵣ) range: $(extrema(lin_aer_optics.k̇[1,:]))")
println("  ω̃̇[1] (∂ω̃/∂nᵣ) range: $(extrema(lin_aer_optics.ω̃̇[1,:]))")

# If ω̃̇ is correct, then k_sca_dot = k̇*ω̃ + k*ω̃̇ should be smooth
# If ω̃̇ was computed as k̇*ω̃̇_true/k̇ = ω̃̇_true (approximately), 
# it should still be roughly correct
for ctr in 1:4
    k_dot = lin_aer_optics.k̇[ctr,:]
    ssa_dot = lin_aer_optics.ω̃̇[ctr,:]
    # Reconstruct k_sca_dot two ways:
    ksca_dot_correct = k_dot .* aer_optics.ω̃ .+ aer_optics.k .* ssa_dot
    # What the code effectively stores as "k̇sca":
    ksca_dot_code = k_dot .* ssa_dot
    mask = abs.(ksca_dot_code) .> 1e-20
    if any(mask)
        ratio = ksca_dot_correct[mask] ./ ksca_dot_code[mask]
        println("  Param $(param_names[ctr]): k̇sca_correct/k̇sca_code ratio: mean=$(round(mean(ratio), sigdigits=4)), std=$(round(std(ratio), sigdigits=4))")
    else
        println("  Param $(param_names[ctr]): k̇sca_code ≈ 0 everywhere")
    end
end

println("\n" * "=" ^ 70)
println("  All diagnostic tests complete")
println("=" ^ 70)
