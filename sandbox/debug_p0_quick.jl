#!/usr/bin/env julia
using vSmartMOM, vSmartMOM.CoreRT
using Distributions, Statistics, Printf

const YAML_FAST = "test/test_parameters/JacobianTestFast.yaml"
params = parameters_from_yaml(YAML_FAST)
model, lin_model = model_from_parameters(LinMode(), params)
NAer = length(params.scattering_params.rt_aerosols)
NGas = size(lin_model.τ̇_abs[1], 1)
NSurf = 1
println("\nRunning linearized RT (debug output from rt_kernel)...")
R, T, dR, dT = rt_run(model, lin_model, NAer, NGas, NSurf)

# Print τ̇ for p₀ at each layer
println("\n=== Layer-level τ̇[p₀] ===")
for iz = 1:length(lin_model.layer_opt_props[1])
    τ = lin_model.layer_opt_props[1][iz].τ[1]
    τ̇_p0 = lin_model.layer_opt_props_lin[1][iz].τ̇[6, 1]
    ϖ̇_p0 = lin_model.layer_opt_props_lin[1][iz].ϖ̇[6, 1]
    @printf("  Layer %d: τ=%.4e, τ̇[6]=%.4e, ϖ̇[6]=%.4e\n", iz, τ, τ̇_p0, ϖ̇_p0)
end

println("\n=== Final dR[p₀] ===")
for iλ = 1:size(R, 3)
    for iVZA = 1:size(R, 1)
        @printf("  VZA=%d λ=%d: R=%.6e, dR[p₀]=%.6e\n", iVZA, iλ, R[iVZA,1,iλ], dR[6,iVZA,1,iλ])
    end
end
