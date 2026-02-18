#!/usr/bin/env julia
# Verify the combined layer optical property derivatives (τ̇, ϖ̇) against FD.
# This tests the createAero + operator chain but NOT the RT kernel.
using vSmartMOM, vSmartMOM.CoreRT
using Distributions, Statistics, Printf

const YAML_FAST = "test/test_parameters/JacobianTestFast.yaml"

# Build base model
params = parameters_from_yaml(YAML_FAST)
model, lin_model = model_from_parameters(LinMode(), params)

# Set up RS_type like rt_run does
using vSmartMOM.InelasticScattering: noRS
RS_type = noRS()
FT = Float64
iBand = 1
pol_type = model.params.polarization_type

# Initialize F₀ and bandSpecLim as rt_run does
RS_type.bandSpecLim = []
nSpec = size(model.τ_abs[iBand], 1)
push!(RS_type.bandSpecLim, 1:nSpec)
F₀ = zeros(FT, pol_type.n, nSpec)
F₀[1,:] .= 1.0
RS_type.F₀ = F₀

m = 0
layer_opt, layer_opt_lin, fScatt = CoreRT.constructCoreOpticalProperties(RS_type, iBand, m, model, lin_model)

Nz = length(layer_opt)
nSpec = length(layer_opt[1].τ)
nParams = size(layer_opt_lin[1].τ̇, 1)
println("=== Combined Layer Optical Properties ===")
println("  Nz=$Nz, nSpec=$nSpec, nParams=$nParams")

# Print base values at layer 4 (peak aerosol)
iz = 4
println("\n--- Layer $iz (peak aerosol) ---")
println("  τ = $(layer_opt[iz].τ)")
println("  ϖ = $(layer_opt[iz].ϖ)")
println("  τ̇[6,:] (p₀) = $(layer_opt_lin[iz].τ̇[6,:])")
println("  ϖ̇[6,:] (p₀) = $(layer_opt_lin[iz].ϖ̇[6,:])")
println("  τ̇[1,:] (τ_ref) = $(layer_opt_lin[iz].τ̇[1,:])")
println("  ϖ̇[1,:] (τ_ref) = $(layer_opt_lin[iz].ϖ̇[1,:])")

# Now perturb p₀ and recompute
c_aero = params.scattering_params.rt_aerosols[1]
p0_base = mean(c_aero.profile)
σp_base = std(c_aero.profile)
δ = p0_base * 1e-4

c_aero.profile = Normal(p0_base + δ, σp_base)
model_pert, lin_model_pert = model_from_parameters(LinMode(), params)
RS_pert = noRS()
RS_pert.bandSpecLim = [1:nSpec]
RS_pert.F₀ = copy(F₀)
layer_opt_pert, _, _ = CoreRT.constructCoreOpticalProperties(RS_pert, iBand, m, model_pert, lin_model_pert)
c_aero.profile = Normal(p0_base, σp_base)  # restore

println("\n=== FD vs Analytic for p₀ (param 6), δ=$δ ===")
for iz = 1:Nz
    τ_fd = (layer_opt_pert[iz].τ .- layer_opt[iz].τ) ./ δ
    ϖ_fd = (layer_opt_pert[iz].ϖ .- layer_opt[iz].ϖ) ./ δ
    τ_an = layer_opt_lin[iz].τ̇[6,:]
    ϖ_an = layer_opt_lin[iz].ϖ̇[6,:]
    
    for iλ = 1:nSpec
        τ_rel = abs(τ_an[iλ]) > 1e-30 ? (τ_fd[iλ] - τ_an[iλ]) / τ_an[iλ] * 100 : NaN
        ϖ_rel = abs(ϖ_an[iλ]) > 1e-30 ? (ϖ_fd[iλ] - ϖ_an[iλ]) / ϖ_an[iλ] * 100 : NaN
        @printf("  Layer %d λ=%d: τ̇_FD=%.4e τ̇_AN=%.4e (%.1f%%) | ϖ̇_FD=%.4e ϖ̇_AN=%.4e (%.1f%%)\n",
            iz, iλ, τ_fd[iλ], τ_an[iλ], τ_rel, ϖ_fd[iλ], ϖ_an[iλ], ϖ_rel)
    end
end

# Also check τ_ref (param 1) for comparison
println("\n=== FD vs Analytic for τ_ref (param 1) ===")
τ_ref_base = params.scattering_params.rt_aerosols[1].τ_ref
δ_τ = τ_ref_base * 1e-4
params.scattering_params.rt_aerosols[1].τ_ref = τ_ref_base + δ_τ
model_pert2, lin_model_pert2 = model_from_parameters(LinMode(), params)
RS_pert2 = noRS()
RS_pert2.bandSpecLim = [1:nSpec]
RS_pert2.F₀ = copy(F₀)
layer_opt_pert2, _, _ = CoreRT.constructCoreOpticalProperties(RS_pert2, iBand, m, model_pert2, lin_model_pert2)
params.scattering_params.rt_aerosols[1].τ_ref = τ_ref_base  # restore

for iz = [3, 4, 5]  # only check layers with significant aerosol
    τ_fd = (layer_opt_pert2[iz].τ .- layer_opt[iz].τ) ./ δ_τ
    ϖ_fd = (layer_opt_pert2[iz].ϖ .- layer_opt[iz].ϖ) ./ δ_τ
    τ_an = layer_opt_lin[iz].τ̇[1,:]
    ϖ_an = layer_opt_lin[iz].ϖ̇[1,:]
    
    for iλ = 1:min(2, nSpec)  # just 2 wavelengths
        τ_rel = abs(τ_an[iλ]) > 1e-30 ? (τ_fd[iλ] - τ_an[iλ]) / τ_an[iλ] * 100 : NaN
        ϖ_rel = abs(ϖ_an[iλ]) > 1e-30 ? (ϖ_fd[iλ] - ϖ_an[iλ]) / ϖ_an[iλ] * 100 : NaN
        @printf("  Layer %d λ=%d: τ̇_FD=%.4e τ̇_AN=%.4e (%.1f%%) | ϖ̇_FD=%.4e ϖ̇_AN=%.4e (%.1f%%)\n",
            iz, iλ, τ_fd[iλ], τ_an[iλ], τ_rel, ϖ_fd[iλ], ϖ_an[iλ], ϖ_rel)
    end
end

# Check Z derivatives for p₀ at peak layer
println("\n=== Z derivative check for p₀ at layer 4 ===")
Z_base = Array(layer_opt[4].Z⁺⁺)
Z_pert = Array(layer_opt_pert[4].Z⁺⁺)
Ż_an = Array(layer_opt_lin[4].Ż⁺⁺)

Ż_fd = (Z_pert .- Z_base) ./ δ
Ż_p0 = Ż_an[6,:,:,:]
println("  Ż_FD[1,1,:] = $(Ż_fd[1,1,:])")
println("  Ż_AN[6,1,1,:] = $(Ż_p0[1,1,:])")
for iλ = 1:nSpec
    rel = abs(Ż_p0[1,1,iλ]) > 1e-30 ? (Ż_fd[1,1,iλ] - Ż_p0[1,1,iλ]) / Ż_p0[1,1,iλ] * 100 : NaN
    @printf("  λ=%d: Ż_FD=%.4e, Ż_AN=%.4e, rel=%.1f%%\n", iλ, Ż_fd[1,1,iλ], Ż_p0[1,1,iλ], rel)
end

# Also check the extractEffectiveProps output (τ_sum, τ̇_sum)
println("\n=== τ̇_sum_all check ===")
_, τ_sum_base, τ̇_sum_base = CoreRT.extractEffectiveProps(layer_opt, layer_opt_lin)
_, τ_sum_pert, _ = CoreRT.extractEffectiveProps(layer_opt_pert, [layer_opt_lin[iz] for iz=1:Nz])

for iz = 1:Nz+1
    τsum_fd = (τ_sum_pert[:,iz] .- τ_sum_base[:,iz]) ./ δ
    τsum_an = τ̇_sum_base[6,:,iz]
    for iλ = 1:min(2, nSpec)
        rel = abs(τsum_an[iλ]) > 1e-30 ? (τsum_fd[iλ] - τsum_an[iλ]) / τsum_an[iλ] * 100 : NaN
        @printf("  τ̇_sum at iz=%d λ=%d: FD=%.4e, AN=%.4e, rel=%.1f%%\n", iz, iλ, τsum_fd[iλ], τsum_an[iλ], rel)
    end
end

println("\nDone!")
