## Import packages

using RadiativeTransfer, RadiativeTransfer.vSmartMOM
using Plots

include("6SV1_R_trues.jl")
R_modeled = zeros(6, 3, 3, 16);
parameters = parameters_from_yaml("/home/rjeyaram/RadiativeTransfer/test/6SV1_1.yaml");

## Case 1: τ = 0.10, (λ = 530 nm), ρ = 0.0

# sza, az, vza
azs = [180, 90, 0]
szas = [23.0739, 53.1301, 78.4630]
R_deltas = zeros(3, 3, 16)
R_true = R_trues[1]

for sza_i in 1:3
    for az_i in 1:3
        parameters.spec_bands = [1e7/530 (1e7/530 +1)]
        parameters.vaz = repeat([azs[az_i]], 16)
        parameters.sza = szas[sza_i]
        parameters.brdf = [RadiativeTransfer.vSmartMOM.LambertianSurfaceScalar(0.0)]
        model = model_from_parameters(parameters);
        model.τ_rayl[1] .= 0.1
        R_modeled[1, sza_i, az_i, :] = vSmartMOM.rt_run(model, i_band=1)[:,1,1] / model.quad_points.μ₀
        R_deltas[sza_i, az_i, :] = abs.(R_true[sza_i][az_i] - R_modeled[1, sza_i, az_i, :]) ./ R_true[sza_i][az_i]
    end
end

# Deltas are maximum( abs(R_true-R_model)/R_true )
@show maximum(R_deltas)

plot(R_modeled[1,1,1,:])
plot!(R_true[1][1])


## Case 2: τ = 0.10, (λ = 530 nm), ρ = 0.25

azs = [180, 90, 0]
szas = [0.0001, 36.8699, 66.4218]
R_deltas = zeros(3, 3, 16)
R_true = R_trues[2]

for sza_i in 1:3
    for az_i in 1:3
        parameters.spec_bands = [1e7/530 (1e7/530 +1)]
        parameters.vaz = repeat([azs[az_i]], 16)
        parameters.sza = szas[sza_i]
        parameters.brdf = [RadiativeTransfer.vSmartMOM.LambertianSurfaceScalar(0.25)]
        model = model_from_parameters(parameters);
        model.τ_rayl[1] .= 0.1
        R_modeled[2, sza_i, az_i, :] = vSmartMOM.rt_run(model, i_band=1)[:,1,1] / model.quad_points.μ₀
        R_deltas[sza_i, az_i, :] = abs.(R_true[sza_i][az_i] - R_modeled[2, sza_i, az_i, :]) ./ R_true[sza_i][az_i]
    end
end

# Deltas are maximum( abs(R_true-R_model)/R_true )
@show maximum(R_deltas)

plot(R_modeled[2,1,1,:])
plot!(R_true[1][1])

## Case 3: τ = 0.25, (λ = 440 nm), ρ = 0.0

azs = [180, 90, 0]
szas = [0.0001, 36.8699, 66.4218]
R_deltas = zeros(3, 3, 16)
R_true = R_trues[3]

for sza_i in 1:3
    for az_i in 1:3
        parameters.spec_bands = [1e7/440 (1e7/440 +1)]
        parameters.vaz = repeat([azs[az_i]], 16)
        parameters.sza = szas[sza_i]
        parameters.brdf = [RadiativeTransfer.vSmartMOM.LambertianSurfaceScalar(0.0)]
        model = model_from_parameters(parameters);
        model.τ_rayl[1] .= 0.25
        R_modeled[3, sza_i, az_i, :] = vSmartMOM.rt_run(model, i_band=1)[:,1,1] / model.quad_points.μ₀
        R_deltas[sza_i, az_i, :] = abs.(R_true[sza_i][az_i] - R_modeled[3, sza_i, az_i, :]) ./ R_true[sza_i][az_i]
    end
end

# Deltas are maximum( abs(R_true-R_model)/R_true )
@show maximum(R_deltas)

plot(R_modeled[3,1,1,:])
plot!(R_true[1][1])

## Case 4: τ = 0.25, (λ = 440 nm), ρ = 0.25

## Case 5: τ = 0.50, (λ = 360 nm), ρ = 0.0

## Case 6: τ = 0.50, (λ = 360 nm), ρ = 0.25

