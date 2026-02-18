#!/usr/bin/env julia
# Detailed p₀ Jacobian diagnostic: show ALL values per (VZA, λ)
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
    return R, dR
end

# Run base model
println("Setting up base model...")
params_base = parameters_from_yaml(YAML_FAST)
R_base, dR_base = run_lin_rt(params_base)
p0_base = mean(params_base.scattering_params.rt_aerosols[1].profile)
σp_base = std(params_base.scattering_params.rt_aerosols[1].profile)
analytic = dR_base[6, :, :, :]  # param 6 = p₀

nVZA = size(R_base, 1)
nStokes = size(R_base, 2)
nλ = size(R_base, 3)

# Compute FD at multiple step sizes
step_sizes = [1e-1, 1e-2, 1e-3, 1e-4]
fd_results = Dict{Float64, Array}()

for rel_δ in step_sizes
    δ = p0_base * rel_δ
    params_p = parameters_from_yaml(YAML_FAST)
    params_p.scattering_params.rt_aerosols[1].profile = Normal(p0_base + δ, σp_base)
    R_p, _ = run_lin_rt(params_p)
    fd_results[rel_δ] = (R_p .- R_base) ./ δ
    @printf("  Computed FD for δ/p₀ = %.0e (δ = %.3f hPa)\n", rel_δ, δ)
end

# Also compute central difference at 1e-3
δ_cd = p0_base * 1e-3
params_p = parameters_from_yaml(YAML_FAST)
params_p.scattering_params.rt_aerosols[1].profile = Normal(p0_base + δ_cd, σp_base)
R_p, _ = run_lin_rt(params_p)
params_m = parameters_from_yaml(YAML_FAST)
params_m.scattering_params.rt_aerosols[1].profile = Normal(p0_base - δ_cd, σp_base)
R_m, _ = run_lin_rt(params_m)
fd_central_1e3 = (R_p .- R_m) ./ (2δ_cd)
println("  Computed central FD for δ/p₀ = 1e-3")

println("\n" * "="^120)
println("p₀ Jacobian Detailed Comparison (p₀ = $(p0_base) hPa, σp = $(σp_base) hPa)")
println("="^120)

# Print header
@printf("%-6s %-8s | %-14s | %-14s | %-14s | %-14s | %-14s | %-14s\n",
    "VZA", "λ", "Analytic", "FD(δ/p₀=0.1)", "FD(δ/p₀=0.01)", "FD(δ/p₀=1e-3)", "FD(δ/p₀=1e-4)", "CD(δ/p₀=1e-3)")
println("-"^120)

for iλ = 1:nλ
    for iVZA = 1:nVZA
        for iS = 1:nStokes
            a = analytic[iVZA, iS, iλ]
            f1 = fd_results[1e-1][iVZA, iS, iλ]
            f2 = fd_results[1e-2][iVZA, iS, iλ]
            f3 = fd_results[1e-3][iVZA, iS, iλ]
            f4 = fd_results[1e-4][iVZA, iS, iλ]
            fc = fd_central_1e3[iVZA, iS, iλ]
            @printf("VZA=%d  λ=%d    | %+13.6e | %+13.6e | %+13.6e | %+13.6e | %+13.6e | %+13.6e\n",
                iVZA, iλ, a, f1, f2, f3, f4, fc)
        end
    end
end

println("\n" * "="^120)
println("Relative errors: |analytic - FD| / |FD|")
println("="^120)
@printf("%-6s %-8s | %-14s | %-14s | %-14s | %-14s | %-14s\n",
    "VZA", "λ", "FD(δ/p₀=0.1)", "FD(δ/p₀=0.01)", "FD(δ/p₀=1e-3)", "FD(δ/p₀=1e-4)", "CD(δ/p₀=1e-3)")
println("-"^120)

for iλ = 1:nλ
    for iVZA = 1:nVZA
        for iS = 1:nStokes
            a = analytic[iVZA, iS, iλ]
            f2 = fd_results[1e-2][iVZA, iS, iλ]
            f3 = fd_results[1e-3][iVZA, iS, iλ]
            f4 = fd_results[1e-4][iVZA, iS, iλ]
            fc = fd_central_1e3[iVZA, iS, iλ]
            f1 = fd_results[1e-1][iVZA, iS, iλ]
            re1 = abs(f1) > 1e-20 ? abs(a - f1) / abs(f1) : NaN
            re2 = abs(f2) > 1e-20 ? abs(a - f2) / abs(f2) : NaN
            re3 = abs(f3) > 1e-20 ? abs(a - f3) / abs(f3) : NaN
            re4 = abs(f4) > 1e-20 ? abs(a - f4) / abs(f4) : NaN
            rec = abs(fc) > 1e-20 ? abs(a - fc) / abs(fc) : NaN
            @printf("VZA=%d  λ=%d    | %13.4f | %13.4f | %13.4f | %13.4f | %13.4f\n",
                iVZA, iλ, re1, re2, re3, re4, rec)
        end
    end
end

# Also check: which parameters show nonzero analytic derivatives?
println("\n" * "="^80)
println("All parameter derivatives at VZA=1, λ=1")
println("="^80)
for ip = 1:size(dR_base, 1)
    val = dR_base[ip, 1, 1, 1]
    @printf("  param %d: dR = %+13.6e\n", ip, val)
end

# Check the base radiance values
println("\n" * "="^80)
println("Base radiance R values")
println("="^80)
for iλ = 1:nλ
    for iVZA = 1:nVZA
        @printf("  VZA=%d  λ=%d: R = %+13.6e\n", iVZA, iλ, R_base[iVZA, 1, iλ])
    end
end
