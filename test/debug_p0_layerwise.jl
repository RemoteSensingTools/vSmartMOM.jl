#!/usr/bin/env julia
# Layer-wise p₀ diagnostic: check which layer contributes the error
using vSmartMOM, vSmartMOM.CoreRT
using vSmartMOM.CoreRT: RT_Aerosol
using Distributions, Statistics
using Printf

const YAML_FAST = "test/test_parameters/JacobianTestFast.yaml"

function run_lin_rt(params)
    model, lin_model = model_from_parameters(LinMode(), params)
    NAer = length(params.scattering_params.rt_aerosols)
    NGas = size(lin_model.τ̇_abs[1], 1)
    NSurf = 1
    R, T, dR, dT = rt_run(model, lin_model, NAer, NGas, NSurf)
    return R, dR, model, lin_model
end

println("Setting up base model...")
params_base = parameters_from_yaml(YAML_FAST)
R_base, dR_base, model_base, lin_model_base = run_lin_rt(params_base)
p0_base = mean(params_base.scattering_params.rt_aerosols[1].profile)
σp_base = std(params_base.scattering_params.rt_aerosols[1].profile)

# Check τ_ref convergence too
println("\n=== τ_ref FD Convergence Test ===")
τ_ref_base = params_base.scattering_params.rt_aerosols[1].τ_ref
analytic_τ = dR_base[1, :, :, :]
for log_δ in [-2, -3, -4, -5, -6]
    δ = τ_ref_base * 10.0^log_δ
    params_pert = parameters_from_yaml(YAML_FAST)
    params_pert.scattering_params.rt_aerosols[1].τ_ref = τ_ref_base + δ
    R_pert, _ = run_lin_rt(params_pert)
    fd = (R_pert .- R_base) ./ δ
    mask = abs.(fd) .> 1e-12
    if any(mask)
        errs = abs.(analytic_τ[mask] .- fd[mask]) ./ abs.(fd[mask])
        @printf("  δ/τ=1e%d: mean_rel=%.4f, max_rel=%.4f\n", log_δ, mean(errs), maximum(errs))
    end
end

# Check σp convergence
println("\n=== σp FD Convergence Test ===")
analytic_σp = dR_base[7, :, :, :]
for log_δ in [-2, -3, -4, -5]
    δ = σp_base * 10.0^log_δ
    params_pert = parameters_from_yaml(YAML_FAST)
    params_pert.scattering_params.rt_aerosols[1].profile = Normal(p0_base, σp_base + δ)
    R_pert, _ = run_lin_rt(params_pert)
    fd = (R_pert .- R_base) ./ δ
    mask = abs.(fd) .> 1e-12
    if any(mask)
        errs = abs.(analytic_σp[mask] .- fd[mask]) ./ abs.(fd[mask])
        @printf("  δ/σ=1e%d: mean_rel=%.4f, max_rel=%.4f\n", log_δ, mean(errs), maximum(errs))
    end
end

# Print layer properties for p₀ 
println("\n=== Layer-level τ̇_sum for p₀ (param 6) ===")
# Access lin_model properties
for (iB, iBand) in enumerate(model_base.params.spec_bands)
    println("Band $iB:")
    for iz = 1:length(lin_model_base.layer_opt_props[iB])
        τ = lin_model_base.layer_opt_props[iB][iz].τ
        τ̇_p0 = lin_model_base.layer_opt_props_lin[iB][iz].τ̇[6, :]
        @printf("  Layer %d: τ=%.4e, τ̇[p₀]=%.4e\n", iz, τ[1], τ̇_p0[1])
    end
end

# Print p₀ derivatives for all params per VZA/λ
println("\n=== Raw values per (VZA, λ) ===")
for ip in [1, 6, 7, 8, 9]
    val = dR_base[ip, 1, 1, 1]
    @printf("  param %d: dR[1,1,1] = %+13.6e\n", ip, val)
end
