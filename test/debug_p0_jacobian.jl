#!/usr/bin/env julia
# ===========================================================================
# Diagnostic: Isolate p₀ derivative error 
# ===========================================================================

using vSmartMOM, vSmartMOM.CoreRT
using vSmartMOM.CoreRT: RT_Aerosol
using Distributions, Statistics, LinearAlgebra
using Printf

const YAML_FAST = "test/test_parameters/JacobianTestFast.yaml"

# === 1. Load base and perturbed models ===
println("Loading base model...")
params_base = parameters_from_yaml(YAML_FAST)
model_base, lin_model_base = model_from_parameters(LinMode(), params_base)

p0_base = mean(params_base.scattering_params.rt_aerosols[1].profile)
σp_base = std(params_base.scattering_params.rt_aerosols[1].profile)
println("  p₀ = $p0_base, σp = $σp_base")

δ = p0_base * 1e-4  # use smaller step for better FD
println("  δ = $δ (relative: $(δ/p0_base))")

println("\nLoading perturbed model...")
params_pert = parameters_from_yaml(YAML_FAST)
params_pert.scattering_params.rt_aerosols[1].profile = Normal(p0_base + δ, σp_base)
model_pert, _ = model_from_parameters(LinMode(), params_pert)

# === 2. Profile derivative (Stage 1) ===
println("\n=== Stage 1: Profile derivative ===")
p_half = model_base.profile.p_half

τₚ_base, dτₚdp₀, _ = vSmartMOM.CoreRT.getAerosolLayerOptProp(
    LinMode(), 1.0, p0_base, σp_base, p_half)
τₚ_pert, _, _ = vSmartMOM.CoreRT.getAerosolLayerOptProp(
    LinMode(), 1.0, p0_base + δ, σp_base, p_half)

fd_profile = (τₚ_pert .- τₚ_base) ./ δ
for iz in 1:length(τₚ_base)
    if abs(fd_profile[iz]) > 1e-20
        rel = abs(dτₚdp₀[iz] - fd_profile[iz]) / abs(fd_profile[iz])
        @printf("  layer %2d: τₚ=%+.4e, dτ/dp₀ analytic=%+.6e, FD=%+.6e, rel=%.3e\n", 
                iz, τₚ_base[iz], dτₚdp₀[iz], fd_profile[iz], rel)
    end
end

# === 3. τ_aer derivative (Stage 2) ===
println("\n=== Stage 2: τ_aer per-layer ===")
iBand = 1
iaer = 1
nSpec = size(model_base.τ_aer[iBand], 2)
nZ = size(model_base.τ_aer[iBand], 3)
println("  τ_aer shape: iaer=1, nSpec=$nSpec, nZ=$nZ")
println("  τ̇_aer shape: $(size(lin_model_base.τ̇_aer[iBand]))")

for iz in 1:nZ
    for iλ in 1:min(1, nSpec)  # just first wavelength
        base_val = model_base.τ_aer[iBand][iaer, iλ, iz]
        pert_val = model_pert.τ_aer[iBand][iaer, iλ, iz]
        fd_val = (pert_val - base_val) / δ
        ana_val = lin_model_base.τ̇_aer[iBand][iaer, 6, iλ, iz]  # param 6 = p₀
        if abs(fd_val) > 1e-30
            rel = abs(ana_val - fd_val) / abs(fd_val)
            @printf("  [λ=1, z=%d]: τ_aer=%.4e, ∂τ/∂p₀ analytic=%+.6e, FD=%+.6e, rel=%.4e\n",
                    iz, base_val, ana_val, fd_val, rel)
        else
            @printf("  [λ=1, z=%d]: τ_aer=%.4e, FD≈0, analytic=%+.6e\n",
                    iz, base_val, ana_val)
        end
    end
end

# === 4. Combined layer properties (Stage 3) ===
println("\n=== Stage 3: Combined optical properties ===")
using vSmartMOM.InelasticScattering
RS_type = noRS{Float64}()
m_fourier = 0

layer_opt_base, layer_opt_lin_base, _ = vSmartMOM.CoreRT.constructCoreOpticalProperties(
    RS_type, [iBand], m_fourier, model_base, lin_model_base)
layer_opt_pert, _, _ = vSmartMOM.CoreRT.constructCoreOpticalProperties(
    RS_type, [iBand], m_fourier, model_pert, lin_model_base)

param_idx_p0 = 6
nLayers = length(layer_opt_base)
println("  nLayers = $nLayers, nParams = $(size(layer_opt_lin_base[1].τ̇, 1))")

for iz in 1:nLayers
    τ_b = layer_opt_base[iz].τ
    τ_p = layer_opt_pert[iz].τ
    ϖ_b = layer_opt_base[iz].ϖ
    ϖ_p = layer_opt_pert[iz].ϖ
    
    fd_τ = (τ_p .- τ_b) ./ δ
    fd_ϖ = (ϖ_p .- ϖ_b) ./ δ
    
    ana_τ = layer_opt_lin_base[iz].τ̇[param_idx_p0, :]
    ana_ϖ = layer_opt_lin_base[iz].ϖ̇[param_idx_p0, :]
    
    iλ = 1
    fd_τ1 = fd_τ[iλ]
    fd_ϖ1 = fd_ϖ[iλ]
    ana_τ1 = ana_τ[iλ]
    ana_ϖ1 = ana_ϖ[iλ]
    
    rel_τ = abs(fd_τ1) > 1e-30 ? abs(ana_τ1 - fd_τ1)/abs(fd_τ1) : NaN
    rel_ϖ = abs(fd_ϖ1) > 1e-30 ? abs(ana_ϖ1 - fd_ϖ1)/abs(fd_ϖ1) : NaN
    
    @printf("  layer %2d: τ=%.4e, dτ/dp₀  ana=%+.4e FD=%+.4e rel=%.3e\n",
            iz, τ_b[iλ], ana_τ1, fd_τ1, rel_τ)
    @printf("           ϖ=%.4e, dϖ/dp₀  ana=%+.4e FD=%+.4e rel=%.3e\n",
            ϖ_b isa Number ? ϖ_b : ϖ_b[iλ], ana_ϖ1, fd_ϖ1, rel_ϖ)
    
    # Z derivative norm
    Z_b = layer_opt_base[iz].Z⁺⁺
    Z_p = layer_opt_pert[iz].Z⁺⁺
    fd_Z = (Z_p .- Z_b) ./ δ
    Ż = layer_opt_lin_base[iz].Ż⁺⁺
    
    if ndims(Ż) == 4
        ana_Z_slice = Ż[param_idx_p0, :, :, :]
    else
        ana_Z_slice = Ż[param_idx_p0, :, :]
    end
    
    fd_Z_n = norm(fd_Z)
    ana_Z_n = norm(ana_Z_slice)
    rel_Z = fd_Z_n > 1e-30 ? abs(ana_Z_n - fd_Z_n)/fd_Z_n : NaN
    @printf("           |Ż⁺⁺| ana=%.4e FD=%.4e rel=%.3e\n\n", ana_Z_n, fd_Z_n, rel_Z)
end

println("\nDiagnostic complete.")
