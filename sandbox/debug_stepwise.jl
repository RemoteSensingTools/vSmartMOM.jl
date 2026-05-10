#!/usr/bin/env julia
# Step-by-step debug script for linearized RT
# Tests elemental → doubling → interaction chain to find where NaN originates

using vSmartMOM, vSmartMOM.CoreRT

println("=" ^ 60)
println("Step-by-step Linearized RT Debug")
println("=" ^ 60)

# 1. Build models
println("\n[Step 1] Building models...")
params = parameters_from_yaml("test/test_parameters/JacobianTest.yaml")
model_fwd = model_from_parameters(params)
model_lin, lin_model = model_from_parameters(LinMode(), params)

# 2. Compare model inputs
println("\n[Step 2] Comparing model inputs...")
println("  τ_rayl fwd sum: ", sum(model_fwd.τ_rayl[1]))
println("  τ_rayl lin sum: ", sum(model_lin.τ_rayl[1]))
println("  τ_rayl match:   ", isapprox(sum(model_fwd.τ_rayl[1]), sum(model_lin.τ_rayl[1]), rtol=1e-10))
println("  τ_abs  fwd sum: ", sum(model_fwd.τ_abs[1]))
println("  τ_abs  lin sum: ", sum(model_lin.τ_abs[1]))
println("  τ_aer  fwd sum: ", sum(model_fwd.τ_aer[1]))
println("  τ_aer  lin sum: ", sum(model_lin.τ_aer[1]))

# Check lin_model derivatives
println("\n  Linearized model derivative NaN check:")
println("  τ̇_abs  NaN: ", any(isnan, lin_model.τ̇_abs[1]))
println("  τ̇_aer  NaN: ", any(isnan, lin_model.τ̇_aer[1]))
println("  τ̇_abs  Inf: ", any(isinf, lin_model.τ̇_abs[1]))
println("  τ̇_aer  Inf: ", any(isinf, lin_model.τ̇_aer[1]))

# Show τ̇_aer ranges per parameter
println("\n  τ̇_aer per-parameter ranges:")
td = lin_model.τ̇_aer[1]
for ip in 1:size(td, 1)
    for ipar in 1:size(td, 2)
        vals = td[ip, ipar, :, :]
        if any(!iszero, vals)
            println("    aer=$ip, par=$ipar: range=($(minimum(vals)), $(maximum(vals)))")
        end
    end
end

# 3. Run forward RT only
println("\n[Step 3] Running forward-only RT...")
R_fwd, T_fwd, _, _ = rt_run(model_fwd)
println("  R_fwd shape: ", size(R_fwd[1]))
println("  R_fwd[1,1]:  ", R_fwd[1][1,1])

# 4. Run linearized RT (debug prints will come from rt_kernel_lin.jl)
println("\n[Step 4] Running linearized RT (debug traces from rt_kernel)...")
NAer = length(params.scattering_params.rt_aerosols)
NGas = size(lin_model.τ̇_abs[1], 1)
NSurf = 1
Nparams = NAer * 7 + NGas + NSurf
println("  NAer=$NAer, NGas=$NGas, NSurf=$NSurf, Nparams=$Nparams")

result = rt_run(model_lin, lin_model, NAer, NGas, NSurf)

# 5. Analyze results
println("\n[Step 5] Analyzing results...")
println("  Return type: ", typeof(result))
println("  Return length: ", length(result))
for (i, r) in enumerate(result)
    if r isa AbstractArray
        println("  result[$i]: $(typeof(r)), size=$(size(r)), NaN=$(any(isnan, r)), Inf=$(any(isinf, r))")
    else
        println("  result[$i]: $(typeof(r)) = $r")
    end
end

# Check forward R consistency
if result[1] isa AbstractArray && R_fwd[1] isa AbstractArray
    if size(result[1]) == size(R_fwd[1])
        max_diff = maximum(abs.(R_fwd[1] .- result[1]))
        println("\n  Forward R consistency:")
        println("  Max |R_fwd - R_lin|: ", max_diff)
        println("  Max relative diff:   ", max_diff / max(maximum(abs.(R_fwd[1])), 1e-30))
    else
        println("\n  Forward R sizes differ: fwd=$(size(R_fwd[1])) vs lin=$(size(result[1]))")
    end
end

println("\n" * "=" ^ 60)
println("Debug complete")
println("=" ^ 60)
