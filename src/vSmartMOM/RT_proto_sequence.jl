using Revise
using Plots
using Pkg
# Pkg.activate(".")
using BenchmarkTools
using RadiativeTransfer
using RadiativeTransfer.Absorption
using RadiativeTransfer.Scattering
using RadiativeTransfer.vSmartMOM
using RadiativeTransfer.SolarModel
# using DelimitedFiles
# using Interact
# using Interpolations
# using InstrumentOperator
# using ForwardDiff

R_0 = [0.03916, 0.04235, 0.04378, 0.04593, 0.04940, 0.05438, 0.05816, 0.06541, 0.07639, 0.08773, 0.09549, 0.11842, 0.13619, 0.17836, 0.22393, 0.27069]

R_90 = [0.03916, 0.03922, 0.03930, 0.03952, 0.04021, 0.04208, 0.04409, 0.04894, 0.05748, 0.06698, 0.07367, 0.09388, 0.10977, 0.14790, 0.18964, 0.23441]

R_180 = [0.03916, 0.03634, 0.03533, 0.03414, 0.03316, 0.03381, 0.03557, 0.04073, 0.05062, 0.06186, 0.06981, 0.09387, 0.11280, 0.15826, 0.20833, 0.26423]

R_trues = [R_0, R_90, R_180]

R_modeled = zeros(16, 3)

for i in 1:3 #[0, 90, 180]

    # Load parameters from file
    # parameters = vSmartMOM.default_parameters()
    parameters = vSmartMOM.parameters_from_yaml("/home/rjeyaram/RadiativeTransfer/test/helper/6SV1_1.yaml");

    parameters.vaz = repeat([[0, 90, 180][i]], 16)

    model = model_from_parameters(parameters);

    model.τ_rayl[1] *= 0.1/sum(model.τ_rayl[1])
    # model.τ_rayl[1] .= 0.1/length(model.τ_rayl[1])

    # @assert sum(model.τ_rayl[1]) == 0.1

    R_modeled[:,i] = vSmartMOM.rt_run(model, i_band=1)[:,:,1]

end

plot(R_modeled[:,1])
plot!(R_trues[1])

plot(R_modeled[:,2])
plot!(R_trues[2])

plot(R_modeled[:,3])
plot!(R_trues[3])

# [0.0000, 11.4783, 16.2602, 23.0739, 32.8599, 43.9455, 50.2082, 58.6677, 66.4218, 71.3371, 73.7398, 78.4630, 80.7931, 84.2608, 86.5602, 88.854]

# "/home/rjeyaram/RadiativeTransfer/test/helper/6SV1_1.yaml"

# params_dict = YAML.load_file("RadiativeTransfer/src/vSmartMOM/ModelParameters/DefaultParameters.yaml")

# vSmartMOM.validate_yaml_parameters(params_dict)

# file = "RadiativeTransfer/src/vSmartMOM/ModelParameters/SampleProfile.yaml"
# profile_nc4 = vSmartMOM.read_atmos_profile(parameters.file, parameters.lat, parameters.lon, parameters.timeIndex);
# profile_yaml = vSmartMOM.read_atmos_profile(file);

# Sets all the "specific" parameters
# parameters = vSmartMOM.default_parameters();

# Generates all the derived attributes from above parameters

#model.params.architecture = RadiativeTransfer.Architectures.GPU();

# function run_auto(x)
#     (x)
#     (x[1] isa ForwardDiff.Dual)
#     (x[1])
#     (typeof(parameters.μ))
#     parameters.μ  = [x[1]]
#     parameters.σ  = [x[2]]
#     parameters.nᵣ = [x[3]]
#     parameters.nᵢ = [x[4]]
#     model = model_from_parameters(parameters);
#     # model.architecture = RadiativeTransfer.Architectures.CPU()
#     model.params.architecture = RadiativeTransfer.Architectures.CPU()
#     R = vSmartMOM.rt_run(model, i_band=1)
#     return R[1,1,:]
# end

# x = [1.3, 2.0, 1.3, 0.000001];

# # @time run_auto(x)
# @time dfdx = ForwardDiff.jacobian(run_auto, x);

# a = 10

# g = x -> ForwardDiff.jacobian(run_auto, x);

# ν_grid = collect((1e7/778):0.015:(1e7/755))

# T = 5777
# black_body = planck_spectrum(T, ν_grid)

# solar = readdlm("RadiativeTransfer/src/solar_merged_20160127_600_26316_100.out")

# solar_idx_start = argmin(abs.(solar[:, 1] .- minimum(ν_grid)))
# solar_idx_end   = argmin(abs.(solar[:, 1] .- maximum(ν_grid)))

# solar_new = solar[(solar_idx_start-10):(solar_idx_end+10), :]

# itp = LinearInterpolation(solar_new[:, 1], 
#                           solar_new[:, 2])

# solar_interpolated = itp.(black_body[:,1])


# plot(black_body[:,1], (solar_interpolated .* black_body[:,2]))
# plot!(black_body[:,1], black_body[:,2])

# solar_out = solar_interpolated .* black_body[:,2]

# earth_out = solar_out .* R[1,1,:]

# file = "/home/cfranken/oco2_L1bScND_18688a_180105_B8100r_180206190633.h5"

# ils_Δ, ils_in, dispersion = InstrumentOperator.read_ils_table(file, "/home/rjeyaram/RadiativeTransfer/src/vSmartMOM/ils_oco2.json");

# ils_pixel = 

# oco2_kernel = VariableKernelInstrument(ils_pixel, ν, ind_out)

# plot!(1e6 ./ νs[argmin(abs.((νs ./ 100) .- solar[end,1])):(end-1)], solar_interpolated * 1e13)

# 1 ./ νs[argmin(abs.((νs ./ 100) .- solar[end,1])):(end-1)]

# a = ones(1500)
# a[380:(end-1)] .= solar_interpolated

# plot(1:1500, a .* Ls)

# [0.03916, 0.04235, 0.04378, 0.04593, 0.04940, 0.05438, 0.05816, 0.06541, 0.07639, 0.08773, 0.09549, 0.11842, 0.13619, 0.17836, 0.22393, 0.27069]