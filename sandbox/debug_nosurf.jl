#!/usr/bin/env julia
# Test with albedo=0 to isolate atmospheric RT from surface interaction
using vSmartMOM, vSmartMOM.CoreRT
using Distributions, Printf

const YAML_FAST = "test/test_parameters/JacobianTestFast.yaml"

# Build base model
params = parameters_from_yaml(YAML_FAST)
# Set albedo to zero to disable surface contribution
params.brdf = [CoreRT.LambertianSurfaceScalar(0.0)]
model, lin_model = model_from_parameters(LinMode(), params)

NAer = length(params.scattering_params.rt_aerosols)
NGas = size(lin_model.τ̇_abs[1], 1)
NSurf = 1

println("=== Albedo=0 test ===")
R_base, T_base, dR_an, dT_an = rt_run(model, lin_model, NAer, NGas, NSurf)

nSpec = size(R_base, 3)
for s = 1:min(4, nSpec)
    @printf("  s=%d: R=%.8e, dR_an[6]=%.6e, dR_an[1]=%.6e\n", s, R_base[1,1,s], dR_an[6,1,1,s], dR_an[1,1,1,s])
end

# FD for p₀
c_aero = params.scattering_params.rt_aerosols[1]
p0_base = mean(c_aero.profile)
σp_base = std(c_aero.profile)
δ = p0_base * 1e-4

c_aero.profile = Normal(p0_base + δ, σp_base)
params.brdf = [CoreRT.LambertianSurfaceScalar(0.0)]
model_fwd, lin_model_fwd = model_from_parameters(LinMode(), params)
R_fwd, _, _, _ = rt_run(model_fwd, lin_model_fwd, NAer, NGas, NSurf)

c_aero.profile = Normal(p0_base - δ, σp_base)
params.brdf = [CoreRT.LambertianSurfaceScalar(0.0)]
model_bwd, lin_model_bwd = model_from_parameters(LinMode(), params)
R_bwd, _, _, _ = rt_run(model_bwd, lin_model_bwd, NAer, NGas, NSurf)

c_aero.profile = Normal(p0_base, σp_base)

println("\n=== Central FD comparison (albedo=0, p₀) ===")
for s = 1:min(4, nSpec)
    dR_fd = (R_fwd[1,1,s] - R_bwd[1,1,s]) / (2δ)
    dR_a  = dR_an[6,1,1,s]
    rel = abs(dR_fd) > 1e-30 ? (dR_a - dR_fd) / dR_fd * 100 : NaN
    @printf("  s=%d: dR_FD=%.6e, dR_AN=%.6e, rel=%.1f%%\n", s, dR_fd, dR_a, rel)
end

# Also τ_ref for reference
println("\n=== τ_ref (albedo=0) ===")
params.brdf = [CoreRT.LambertianSurfaceScalar(0.0)]
τ_ref_base = params.scattering_params.rt_aerosols[1].τ_ref
δ_τ = τ_ref_base * 1e-4
params.scattering_params.rt_aerosols[1].τ_ref = τ_ref_base + δ_τ
model_fwd2, lin_model_fwd2 = model_from_parameters(LinMode(), params)
R_fwd2, _, _, _ = rt_run(model_fwd2, lin_model_fwd2, NAer, NGas, NSurf)
params.scattering_params.rt_aerosols[1].τ_ref = τ_ref_base - δ_τ
model_bwd2, lin_model_bwd2 = model_from_parameters(LinMode(), params)
R_bwd2, _, _, _ = rt_run(model_bwd2, lin_model_bwd2, NAer, NGas, NSurf)
params.scattering_params.rt_aerosols[1].τ_ref = τ_ref_base

for s = 1:min(2, nSpec)
    dR_fd = (R_fwd2[1,1,s] - R_bwd2[1,1,s]) / (2δ_τ)
    dR_a  = dR_an[1,1,1,s]
    rel = abs(dR_fd) > 1e-30 ? (dR_a - dR_fd) / dR_fd * 100 : NaN
    @printf("  s=%d: dR_FD=%.6e, dR_AN=%.6e, rel=%.1f%%\n", s, dR_fd, dR_a, rel)
end

println("\nDone!")
