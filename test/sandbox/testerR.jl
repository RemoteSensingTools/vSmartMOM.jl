parameters = parameters_from_yaml("test/benchmarks/6SV1_1_simple.yaml");
τ = 1.0
ρ = 0.00
λ = 400.0
parameters.brdf = [vSmartMOM.CoreRT.LambertianSurfaceScalar(ρ)]
parameters.spec_bands = [1e7/λ (1e7/λ + 1)]
model = model_from_parameters(parameters);
model.τ_rayl[1] .= τ

R_SFI_0, T_SFI, ieR_SFI, ieT_SFI, hdr, bhr_uw, bhr_dw = vSmartMOM.CoreRT.rt_run(model)


ρ = 0.20
parameters.brdf = [vSmartMOM.CoreRT.LambertianSurfaceScalar(ρ)]
model = model_from_parameters(parameters);
model.τ_rayl[1] .= τ
R_SFI_025, T_SFI, ieR_SFI, ieT_SFI, hdr, bhr_uw, bhr_dw = vSmartMOM.CoreRT.rt_run(model)

ρ = 0.80
parameters.brdf = [vSmartMOM.CoreRT.LambertianSurfaceScalar(ρ)]
model = model_from_parameters(parameters);
model.τ_rayl[1] .= τ
R_SFI_08, T_SFI, ieR_SFI, ieT_SFI, hdr, bhr_uw, bhr_dw = vSmartMOM.CoreRT.rt_run(model)



# Vijays for μ₀ = 0.46 and μ = 1 and μ = 0.7
# Albedo = 0.0
alb_0   = [0.15888669  0.0683381 ; 0.19466293 0.08885206]
alb_025 = [0.20013929 0.06833818 ; 0.23075985 0.08895853]
alb_08  = [0.34139850 0.06833818 ; 0.35436472 0.08932314]

(π *  R_SFI_0[1:2,1:2]   - alb_0) 
(π *  R_SFI_025[1:2,1:2] - alb_025) 
(π *  R_SFI_08[1:2,1:2]  - alb_08) 

#=
julia> (π *  R_SFI_0[1:2,1:2]   - alb_0)
2×2 Matrix{Float64}:
  0.000410157  -0.000569149
 -0.000265577   0.000167046

julia> (π *  R_SFI_025[1:2,1:2] - alb_025)
2×2 Matrix{Float64}:
 -0.00863807  -0.000569229
 -0.00823166   0.000210581

julia> (π *  R_SFI_08[1:2,1:2]  - alb_08) 
2×2 Matrix{Float64}:
  0.000482343  -0.000569229
 -0.000478472   0.000546427
 =#