#!/usr/bin/env julia
# Focused dR/dp₀ convergence test: central FD vs analytic
using vSmartMOM, vSmartMOM.CoreRT
using Distributions, Printf

const YAML_FAST = "test/test_parameters/JacobianTestFast.yaml"

# Build base model
params = parameters_from_yaml(YAML_FAST)
model, lin_model = model_from_parameters(LinMode(), params)

NAer = length(params.scattering_params.rt_aerosols)
NGas = size(lin_model.τ̇_abs[1], 1)
NSurf = 1

# Compute base R and analytic dR
println("Computing base RT...")
R_base, T_base, dR_an, dT_an = rt_run(model, lin_model, NAer, NGas, NSurf)

nSpec = size(R_base, 3)
nVZA = size(R_base, 1)

println("\n=== Analytic dR[p₀] (param 6) ===")
for s = 1:min(4, nSpec)
    @printf("  s=%d: R=%.8e, dR_an[6]=%.6e, dR_an[1]=%.6e\n", s, R_base[1,1,s], dR_an[6,1,1,s], dR_an[1,1,1,s])
end

# Now do central FD for p₀ at various step sizes
c_aero = params.scattering_params.rt_aerosols[1]
p0_base = mean(c_aero.profile)
σp_base = std(c_aero.profile)
println("\n  p₀=$p0_base, σp=$σp_base")

println("\n=== Central FD convergence for dR/dp₀ ===")
for log_δ in [-2, -3, -4, -5, -6]
    δ = p0_base * 10.0^log_δ
    
    # Forward perturbation
    c_aero.profile = Normal(p0_base + δ, σp_base)
    model_fwd, lin_model_fwd = model_from_parameters(LinMode(), params)
    R_fwd, _, _, _ = rt_run(model_fwd, lin_model_fwd, NAer, NGas, NSurf)
    
    # Backward perturbation
    c_aero.profile = Normal(p0_base - δ, σp_base)
    model_bwd, lin_model_bwd = model_from_parameters(LinMode(), params)
    R_bwd, _, _, _ = rt_run(model_bwd, lin_model_bwd, NAer, NGas, NSurf)
    
    # Restore
    c_aero.profile = Normal(p0_base, σp_base)
    
    for s = 1:min(2, nSpec)
        dR_fd = (R_fwd[1,1,s] - R_bwd[1,1,s]) / (2δ)
        dR_a  = dR_an[6,1,1,s]
        rel = abs(dR_fd) > 1e-30 ? (dR_a - dR_fd) / dR_fd * 100 : NaN
        @printf("  δ=%.1e s=%d: dR_FD=%.6e, dR_AN=%.6e, rel=%.1f%%\n", δ, s, dR_fd, dR_a, rel)
    end
end

# Also do FD for τ_ref as sanity check
println("\n=== Central FD convergence for dR/dτ_ref ===")
τ_ref_base = params.scattering_params.rt_aerosols[1].τ_ref
for log_δ in [-3, -4, -5]
    δ = τ_ref_base * 10.0^log_δ
    
    params.scattering_params.rt_aerosols[1].τ_ref = τ_ref_base + δ
    model_fwd, lin_model_fwd = model_from_parameters(LinMode(), params)
    R_fwd, _, _, _ = rt_run(model_fwd, lin_model_fwd, NAer, NGas, NSurf)
    
    params.scattering_params.rt_aerosols[1].τ_ref = τ_ref_base - δ
    model_bwd, lin_model_bwd = model_from_parameters(LinMode(), params)
    R_bwd, _, _, _ = rt_run(model_bwd, lin_model_bwd, NAer, NGas, NSurf)
    
    params.scattering_params.rt_aerosols[1].τ_ref = τ_ref_base
    
    for s = 1:min(2, nSpec)
        dR_fd = (R_fwd[1,1,s] - R_bwd[1,1,s]) / (2δ)
        dR_a  = dR_an[1,1,1,s]
        rel = abs(dR_fd) > 1e-30 ? (dR_a - dR_fd) / dR_fd * 100 : NaN
        @printf("  δ=%.1e s=%d: dR_FD=%.6e, dR_AN=%.6e, rel=%.1f%%\n", δ, s, dR_fd, dR_a, rel)
    end
end

println("\nDone!")
