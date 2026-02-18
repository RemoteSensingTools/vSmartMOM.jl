#!/usr/bin/env julia
# Verify layer-level τ̇ for p₀ and σp against FD — direct function test
using vSmartMOM, vSmartMOM.CoreRT
using Distributions, Printf

# Test 1: Direct test of getAerosolLayerOptProp derivatives
println("=== Test 1: getAerosolLayerOptProp derivative check ===")
p_half = [0.0, 200.0, 400.0, 600.0, 800.0, 1013.0]
p₀ = 800.0
σp = 50.0
total_τ = 1.0

τ_base, dτ_dp₀, dτ_dσp = CoreRT.getAerosolLayerOptProp(LinMode(), total_τ, p₀, σp, p_half)
println("  p₀=$p₀, σp=$σp")
println("  τ_base = $τ_base  (sum=$(sum(τ_base)))")
println("  dτ/dp₀ = $dτ_dp₀ (sum=$(sum(dτ_dp₀)))")
println("  dτ/dσp = $dτ_dσp (sum=$(sum(dτ_dσp)))")

# FD check for p₀
for δ in [1.0, 0.1, 0.01, 0.001]
    τ_pert, _, _ = CoreRT.getAerosolLayerOptProp(LinMode(), total_τ, p₀ + δ, σp, p_half)
    fd = (τ_pert .- τ_base) ./ δ
    rel_errs = [abs(dτ_dp₀[i]) > 1e-15 ? (fd[i] - dτ_dp₀[i]) / dτ_dp₀[i] * 100 : NaN for i in eachindex(fd)]
    @printf("  δ=%.0e: FD=%.6e..., AN=%.6e..., rel_err=%.2f%%...\n", δ, fd[3], dτ_dp₀[3], rel_errs[3])
end

# FD check for σp
println("\n  σp FD check:")
for δ in [1.0, 0.1, 0.01, 0.001]
    τ_pert, _, _ = CoreRT.getAerosolLayerOptProp(LinMode(), total_τ, p₀, σp + δ, p_half)
    fd = (τ_pert .- τ_base) ./ δ
    rel_errs = [abs(dτ_dσp[i]) > 1e-15 ? (fd[i] - dτ_dσp[i]) / dτ_dσp[i] * 100 : NaN for i in eachindex(fd)]
    @printf("  δ=%.0e: FD=%.6e..., AN=%.6e..., rel_err=%.2f%%...\n", δ, fd[3], dτ_dσp[3], rel_errs[3])
end

# Test 2: Full model — compare τ̇_aer for p₀ against FD 
println("\n=== Test 2: Full model τ̇_aer[p₀] verification ===")
const YAML_FAST = "test/test_parameters/JacobianTestFast.yaml"
params = parameters_from_yaml(YAML_FAST)
model, lin_model = model_from_parameters(LinMode(), params)

# Get base τ_aer for first aerosol, first band
τ_aer_base = lin_model.τ̇_aer[1]  # [iaer, 7params, nSpec, nLayers]
i_aer = 1
Nz = size(model.τ_aer[1], 3)

println("  τ_aer shape: $(size(model.τ_aer[1]))")
println("  τ̇_aer shape: $(size(lin_model.τ̇_aer[1]))")

# Print analytic τ̇_aer for p₀ and σp at each layer  
println("\n  Layer-level τ̇_aer (analytic) at spec=1:")
for iz = 1:Nz
    τ_val = model.τ_aer[1][i_aer, 1, iz]
    τ̇_tref = lin_model.τ̇_aer[1][i_aer, 1, 1, iz]
    τ̇_p0 = lin_model.τ̇_aer[1][i_aer, 6, 1, iz]
    τ̇_sp = lin_model.τ̇_aer[1][i_aer, 7, 1, iz]
    @printf("    Layer %d: τ=%.4e, τ̇[τ_ref]=%.4e, τ̇[p₀]=%.4e, τ̇[σp]=%.4e\n", iz, τ_val, τ̇_tref, τ̇_p0, τ̇_sp)
end

# Now perturb p₀ and recompute — need to perturb the profile distribution
c_aero = params.scattering_params.rt_aerosols[1]
p0_base = mean(c_aero.profile)
σp_base = std(c_aero.profile)
println("\n  Aerosol profile: p₀=$p0_base, σp=$σp_base")

δ = p0_base * 1e-4
println("  FD step δ = $δ")

# Perturb p₀
c_aero.profile = Normal(p0_base + δ, σp_base)
model_pert, _ = model_from_parameters(LinMode(), params)

println("\n  Layer-level FD vs Analytic for p₀:")
for iz = 1:Nz
    τ_b = model.τ_aer[1][i_aer, 1, iz]
    τ_p = model_pert.τ_aer[1][i_aer, 1, iz]
    fd = (τ_p - τ_b) / δ
    an = lin_model.τ̇_aer[1][i_aer, 6, 1, iz]
    rel = abs(an) > 1e-30 ? (fd - an) / an * 100 : NaN
    @printf("    Layer %d: FD=%.6e, AN=%.6e, rel_err=%.2f%%\n", iz, fd, an, rel)
end

# Restore profile
c_aero.profile = Normal(p0_base, σp_base)

# Test 3: Full end-to-end FD for R
println("\n=== Test 3: Full R FD check for p₀ ===")
NAer = length(params.scattering_params.rt_aerosols)
NGas = size(lin_model.τ̇_abs[1], 1)
NSurf = 1

R_base, T_base, dR, dT = rt_run(model, lin_model, NAer, NGas, NSurf)
# Recompute with perturbed p₀
c_aero.profile = Normal(p0_base + δ, σp_base)
model_pert2, lin_model_pert2 = model_from_parameters(LinMode(), params)
R_pert, T_pert, _, _ = rt_run(model_pert2, lin_model_pert2, NAer, NGas, NSurf)
c_aero.profile = Normal(p0_base, σp_base)

println("  R shape: $(size(R_base)), dR shape: $(size(dR))")
for iλ = 1:min(4, size(R_base, 3))
    for iVZA = 1:min(3, size(R_base, 1))
        r_fd = (R_pert[iVZA, 1, iλ] - R_base[iVZA, 1, iλ]) / δ
        r_an = dR[6, iVZA, 1, iλ]
        rel = abs(r_an) > 1e-30 ? (r_fd - r_an) / r_fd * 100 : NaN
        @printf("  VZA=%d λ=%d: R_FD=%.6e, R_AN=%.6e, rel_err=%.2f%%\n", iVZA, iλ, r_fd, r_an, rel)
    end
end

println("\nDone!")
