#!/usr/bin/env julia
# Quick convergence test for p₀ Jacobian: does FD → analytic as δ → 0?
using vSmartMOM, vSmartMOM.CoreRT
using vSmartMOM.CoreRT: RT_Aerosol
using Distributions, Statistics

const YAML_FAST = "test/test_parameters/JacobianTestFast.yaml"

function run_lin_rt(params)
    model, lin_model = model_from_parameters(LinMode(), params)
    NAer = length(params.scattering_params.rt_aerosols)
    NGas = size(lin_model.τ̇_abs[1], 1)
    NSurf = 1
    R, T, dR, dT = rt_run(model, lin_model, NAer, NGas, NSurf)
    return R, dR
end

function rel_errors(analytic, fd; threshold=1e-12)
    mask = abs.(fd) .> threshold
    if any(mask)
        errs = abs.(analytic[mask] .- fd[mask]) ./ abs.(fd[mask])
        return (max=maximum(errs), mean=mean(errs))
    else
        ae = abs.(analytic .- fd)
        return (max=maximum(ae), mean=mean(ae))
    end
end

# Run base model
println("Setting up base model...")
params_base = parameters_from_yaml(YAML_FAST)
R_base, dR_base = run_lin_rt(params_base)
p0_base = mean(params_base.scattering_params.rt_aerosols[1].profile)
σp_base = std(params_base.scattering_params.rt_aerosols[1].profile)
analytic = dR_base[6, :, :, :]  # param 6 = p₀

println("\n=== p₀ FD Convergence Test (one-sided) ===")
for log_δ in [-1, -2, -3, -4, -5]
    δ = p0_base * 10.0^log_δ
    params_pert = parameters_from_yaml(YAML_FAST)
    params_pert.scattering_params.rt_aerosols[1].profile = Normal(p0_base + δ, σp_base)
    R_pert, _ = run_lin_rt(params_pert)
    fd = (R_pert .- R_base) ./ δ
    err = rel_errors(analytic, fd)
    println("  δ/p₀=1e$(log_δ) (δ=$(round(δ, sigdigits=3))): mean_rel=$(round(err.mean, sigdigits=4)), max_rel=$(round(err.max, sigdigits=4))")
end

println("\n=== p₀ FD Convergence Test (central difference) ===")
for log_δ in [-2, -3, -4, -5]
    δ = p0_base * 10.0^log_δ
    params_p = parameters_from_yaml(YAML_FAST)
    params_p.scattering_params.rt_aerosols[1].profile = Normal(p0_base + δ, σp_base)
    R_p, _ = run_lin_rt(params_p)
    
    params_m = parameters_from_yaml(YAML_FAST)
    params_m.scattering_params.rt_aerosols[1].profile = Normal(p0_base - δ, σp_base)
    R_m, _ = run_lin_rt(params_m)
    
    fd_central = (R_p .- R_m) ./ (2δ)
    err = rel_errors(analytic, fd_central)
    println("  δ/p₀=1e$(log_δ) (δ=$(round(δ, sigdigits=3))): mean_rel=$(round(err.mean, sigdigits=4)), max_rel=$(round(err.max, sigdigits=4))")
end

println("\n=== Raw values comparison (first VZA, first Stokes, all λ) ===")
δ = p0_base * 1e-4
params_p = parameters_from_yaml(YAML_FAST)
params_p.scattering_params.rt_aerosols[1].profile = Normal(p0_base + δ, σp_base)
R_p, _ = run_lin_rt(params_p)
params_m = parameters_from_yaml(YAML_FAST)
params_m.scattering_params.rt_aerosols[1].profile = Normal(p0_base - δ, σp_base)
R_m, _ = run_lin_rt(params_m)
fd_central = (R_p .- R_m) ./ (2δ)

for iλ in 1:size(R_base, 3)
    a = analytic[1, 1, iλ]
    f = fd_central[1, 1, iλ]
    re = abs(a) > 1e-20 ? abs(a - f) / abs(f) : abs(a - f)
    println("  λ=$iλ: analytic=$(round(a, sigdigits=6)), fd_central=$(round(f, sigdigits=6)), rel_err=$(round(re, sigdigits=4))")
end
