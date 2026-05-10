#!/usr/bin/env julia
# Perturbation test for linearized RT Jacobians
# Compares analytic Jacobians (Ṙ) against finite-difference (ΔR/Δp)
# using the same rt_run function for both

using vSmartMOM, vSmartMOM.CoreRT
using vSmartMOM.CoreRT: RT_Aerosol
using Distributions

println("=" ^ 60)
println("Perturbation Test: Jacobians vs Finite Difference")
println("=" ^ 60)

# Helper: run linearized RT and return (R, dR)
function run_lin_rt(params)
    model, lin_model = model_from_parameters(LinMode(), params)
    NAer = length(params.scattering_params.rt_aerosols)
    NGas = size(lin_model.τ̇_abs[1], 1)
    NSurf = 1
    R, T, dR, dT = rt_run(model, lin_model, NAer, NGas, NSurf)
    return R, dR, NAer, NGas, NSurf
end

# Helper: run forward-only RT using linearized path (to get same format)
function run_fwd_only(params)
    model, lin_model = model_from_parameters(LinMode(), params)
    NAer = length(params.scattering_params.rt_aerosols)
    NGas = size(lin_model.τ̇_abs[1], 1)
    NSurf = 1
    R, T, dR, dT = rt_run(model, lin_model, NAer, NGas, NSurf)
    return R
end

# Helper: compute relative error
function rel_error(analytic, fd)
    mask = abs.(fd) .> 1e-15  # avoid division by tiny FD values
    if any(mask)
        return maximum(abs.(analytic[mask] .- fd[mask]) ./ abs.(fd[mask]))
    else
        return maximum(abs.(analytic .- fd))
    end
end

function mean_rel_error(analytic, fd)
    mask = abs.(fd) .> 1e-15
    if any(mask)
        return mean(abs.(analytic[mask] .- fd[mask]) ./ abs.(fd[mask]))
    else
        return mean(abs.(analytic .- fd))
    end
end

using Statistics

# ====================================================================
# Load base parameters
# ====================================================================
println("\n[Setup] Loading parameters...")
params_base = parameters_from_yaml("test/test_parameters/JacobianTest.yaml")

# ====================================================================
# Phase 1: Run baseline linearized RT 
# ====================================================================
println("\n[Phase 1] Running baseline linearized RT...")
R_base, dR_base, NAer, NGas, NSurf = run_lin_rt(params_base)
Nparams = NAer * 7 + NGas + NSurf

println("  R shape:  $(size(R_base))")
println("  dR shape: $(size(dR_base))")
println("  NAer=$NAer, NGas=$NGas, NSurf=$NSurf, Nparams=$Nparams")
println("  Any NaN in R:  $(any(isnan, R_base))")
println("  Any NaN in dR: $(any(isnan, dR_base))")

# Print dR ranges per parameter
println("\n  dR ranges per parameter:")
for ip in 1:Nparams
    vals = dR_base[ip, :, :, :]
    println("    param $ip: range=($(minimum(vals)), $(maximum(vals)))")
end

# ====================================================================
# Phase 2: Surface Albedo Perturbation Test
# ====================================================================
println("\n" * "=" ^ 60)
println("[Phase 2] Surface Albedo Jacobian Test")
println("=" ^ 60)

albedo_base = 0.05  # from YAML
δ_albedo = 1e-4

# Create perturbed params
params_pert = parameters_from_yaml("test/test_parameters/JacobianTest.yaml")
# Modify surface albedo
params_pert.brdf = [LambertianSurfaceScalar(albedo_base + δ_albedo)]

println("  Base albedo: $albedo_base")
println("  Perturbed:   $(albedo_base + δ_albedo)")
println("  δ = $δ_albedo")

R_pert = run_fwd_only(params_pert)

# Finite difference
FD_albedo = (R_pert .- R_base) ./ δ_albedo

# Analytic: surface albedo is the LAST parameter (index Nparams)
isurf = Nparams
analytic_albedo = dR_base[isurf, :, :, :]

println("\n  FD albedo Jacobian range:       ($(minimum(FD_albedo)), $(maximum(FD_albedo)))")
println("  Analytic albedo Jacobian range: ($(minimum(analytic_albedo)), $(maximum(analytic_albedo)))")
println("  Max relative error: $(rel_error(analytic_albedo, FD_albedo))")
println("  Mean relative error: $(mean_rel_error(analytic_albedo, FD_albedo))")

# ====================================================================
# Phase 3: Aerosol τ_ref Perturbation Test  
# ====================================================================
println("\n" * "=" ^ 60)
println("[Phase 3] Aerosol τ_ref Jacobian Test")
println("=" ^ 60)

τ_ref_base = params_base.scattering_params.rt_aerosols[1].τ_ref
δ_τ_ref = τ_ref_base * 1e-3  # 0.1% perturbation

params_pert2 = parameters_from_yaml("test/test_parameters/JacobianTest.yaml")
# Modify τ_ref directly (RT_Aerosol is mutable)
params_pert2.scattering_params.rt_aerosols[1].τ_ref = τ_ref_base + δ_τ_ref

println("  Base τ_ref: $τ_ref_base")
println("  δ = $δ_τ_ref")

R_pert2 = run_fwd_only(params_pert2)

FD_tau = (R_pert2 .- R_base) ./ δ_τ_ref

# Analytic: τ_ref is parameter 1 (first aerosol, first sub-parameter)
itau = 1
analytic_tau = dR_base[itau, :, :, :]

println("\n  FD τ_ref Jacobian range:       ($(minimum(FD_tau)), $(maximum(FD_tau)))")
println("  Analytic τ_ref Jacobian range: ($(minimum(analytic_tau)), $(maximum(analytic_tau)))")
println("  Max relative error: $(rel_error(analytic_tau, FD_tau))")
println("  Mean relative error: $(mean_rel_error(analytic_tau, FD_tau))")

# ====================================================================
# Phase 4: Aerosol Profile p₀ Perturbation Test
# ====================================================================
println("\n" * "=" ^ 60)
println("[Phase 4] Aerosol Profile p₀ Jacobian Test")
println("=" ^ 60)

aer_base = params_base.scattering_params.rt_aerosols[1]
p0_base = mean(aer_base.profile)
σp_base = std(aer_base.profile)

δ_p0 = p0_base * 1e-3

params_pert3 = parameters_from_yaml("test/test_parameters/JacobianTest.yaml")
# Perturb p₀ by shifting the distribution mean (RT_Aerosol is mutable)
params_pert3.scattering_params.rt_aerosols[1].profile = Normal(p0_base + δ_p0, σp_base)

println("  Base p₀: $p0_base")
println("  δ = $δ_p0")

R_pert3 = run_fwd_only(params_pert3)

FD_p0 = (R_pert3 .- R_base) ./ δ_p0

# Analytic: p₀ is parameter 6 (first aerosol, 6th sub-parameter)
ip0 = 6
analytic_p0 = dR_base[ip0, :, :, :]

println("\n  FD p₀ Jacobian range:       ($(minimum(FD_p0)), $(maximum(FD_p0)))")
println("  Analytic p₀ Jacobian range: ($(minimum(analytic_p0)), $(maximum(analytic_p0)))")
println("  Max relative error: $(rel_error(analytic_p0, FD_p0))")
println("  Mean relative error: $(mean_rel_error(analytic_p0, FD_p0))")

println("\n" * "=" ^ 60)
println("All perturbation tests complete")
println("=" ^ 60)
