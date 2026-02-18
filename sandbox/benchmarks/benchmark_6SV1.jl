## Import packages

using vSmartMOM, vSmartMOM.CoreRT
using Statistics, Plots

include("6SV1_R_trues.jl")
R_modeled_all = zeros(6, 3, 3, 16);
R_deltas_all = zeros(6, 3, 3, 16);
parameters = parameters_from_yaml("/home/rjeyaram/vSmartMOM/test/6SV1_1.yaml");

## Function to perform test, given conditions

function test_against_6SV1(case_i, azs, szas, λ, τ, ρ)
    for sza_i in 1:length(szas)
        for az_i in 1:length(azs)
            R_true = R_trues[case_i]
            R_modeled = view(R_modeled_all, case_i, :, :, :)
            R_deltas = view(R_deltas_all, case_i, :, :, :)
            parameters.spec_bands = [1e7/λ (1e7/λ + 1)]
            parameters.vaz = repeat([azs[az_i]], 16)
            parameters.sza = szas[sza_i]
            parameters.brdf = [vSmartMOM.CoreRT.LambertianSurfaceScalar(ρ * π)]
            model = model_from_parameters(parameters);
            model.τ_rayl[1] .= τ
            R_modeled[sza_i, az_i, :] = CoreRT.rt_run(model, i_band=1)[:,1,1] / model.quad_points.μ₀
            R_deltas[sza_i, az_i, :] = abs.(R_true[sza_i][az_i] - R_modeled[sza_i, az_i, :]) ./ R_true[sza_i][az_i]
        end
    end

    println("Maximum discrepancy: ", maximum(R_deltas_all[case_i,:,:,:])*100, "%")
    println("Minimum discrepancy: ", minimum(R_deltas_all[case_i,:,:,:])*100, "%")
    println("Average discrepancy: ", mean(R_deltas_all[case_i,:,:,:])*100, "%")
    println("Std dev discrepancy: ", std(R_deltas_all[case_i,:,:,:])*100, "%")
end

## Case 1: τ = 0.10, (λ = 530 nm), ρ = 0.0

test_against_6SV1(1,                               # test case no. 
                  [180, 90, 0],                    # AZ
                  [23.0739, 53.1301, 78.4630],     # SZA
                  530,                             # λ
                  0.1,                             # τ
                  0.0)                             # ρ

## Case 2: τ = 0.10, (λ = 530 nm), ρ = 0.25

test_against_6SV1(2,                               # test case no. 
                  [180, 90, 0],                    # AZ
                  [0.0001, 36.8699, 66.4218],      # SZA
                  530,                             # λ
                  0.1,                             # τ
                  0.25)                            # ρ

## Case 3: τ = 0.25, (λ = 440 nm), ρ = 0.0

test_against_6SV1(3,                               # test case no. 
                  [180, 90, 0],                    # AZ
                  [0.0001, 36.8699, 66.4218],      # SZA
                  440,                             # λ
                  0.25,                            # τ
                  0.0)                             # ρ

## Case 4: τ = 0.25, (λ = 440 nm), ρ = 0.25

test_against_6SV1(4,                               # test case no. 
                  [180, 90, 0],                    # AZ
                  [23.0739, 53.1301, 78.4630],     # SZA
                  440,                             # λ
                  0.25,                            # τ
                  0.25)                            # ρ

## Case 5: τ = 0.50, (λ = 360 nm), ρ = 0.0

test_against_6SV1(5,                               # test case no. 
                  [180, 90, 0],                    # AZ
                  [23.0739, 53.1301, 78.4630],     # SZA
                  360,                             # λ
                  0.50,                            # τ
                  0.0)                             # ρ

## Case 6: τ = 0.50, (λ = 360 nm), ρ = 0.25

test_against_6SV1(6,                               # test case no. 
                  [180, 90, 0],                    # AZ
                  [0.0001, 36.8699, 66.4218],      # SZA
                  360,                             # λ
                  0.50,                            # τ
                  0.25)                            # ρ