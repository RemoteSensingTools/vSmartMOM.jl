using Revise
using Plots
using Pkg
# Pkg.activate(".")
using BenchmarkTools
using RadiativeTransfer
using RadiativeTransfer.Absorption
using RadiativeTransfer.Scattering
using RadiativeTransfer.vSmartMOM
using DelimitedFiles
using Interact
using Interpolations

# Load parameters from file
parameters = vSmartMOM.parameters_from_yaml("RadiativeTransfer/test/helper/ThreeBandsParameters.yaml")

# params_dict = YAML.load_file("RadiativeTransfer/src/vSmartMOM/ModelParameters/DefaultParameters.yaml")

# vSmartMOM.validate_yaml_parameters(params_dict)

# file = "RadiativeTransfer/src/vSmartMOM/ModelParameters/SampleProfile.yaml"
# profile_nc4 = vSmartMOM.read_atmos_profile(parameters.file, parameters.lat, parameters.lon, parameters.timeIndex);
# profile_yaml = vSmartMOM.read_atmos_profile(file);

# Sets all the "specific" parameters
# parameters = vSmartMOM.default_parameters();

# Generates all the derived attributes from above parameters
model = model_from_parameters(parameters);
#model.params.architecture = RadiativeTransfer.Architectures.GPU();

R = vSmartMOM.rt_run(model)


solar_R = readdlm("RadiativeTransfer/src/solar_merged_20160127_600_26316_100.out")

T = 5777

black_body = planck_spectrum(T, collect(1:1500), wavelength_flag=true)

R1 = planck_spectrum(T, collect((1e7/777):0.015:(1e7/757)))
R2 = planck_spectrum(T, reverse(1e7 ./ (collect((1e7/777):0.015:(1e7/757)))), wavelength_flag=true)


# h = 6.626070e-34 # J⋅Hz−1
# c = 3e8    # m / s
# k = 1.380649e-23 # J⋅K−1
# λs = collect(1:1500) .* 1e-9
# νs = 1 ./ λs 

# @manipulate for T in 3000:100:10000
#     Ls = map(ν-> (2*h*c^2) * (1/((1/ν)^5*(exp(h*c/((1/ν)*k*T)) - 1))) , νs)
#     plot(λs*1e6, Ls)
# end


peak = λs[argmax(Ls)]
title!("T=$(T), peak=$(peak*1e6)μm")

ν_min = νs[end] / 100
ν_max = νs[1] / 100

solar_idx_start = argmin(abs.(solar[:, 1] .- ν_min))
solar_idx_end   = argmin(abs.(solar[:, 1] .- ν_max))

itp = LinearInterpolation(solar[solar_idx_start:solar_idx_end, 1], 
                          solar[solar_idx_start:solar_idx_end, 2])

solar_interpolated = itp.(νs[argmin(abs.((νs ./ 100) .- solar[end,1])):(end-1)] / 100)

plot!(1e6 ./ νs[argmin(abs.((νs ./ 100) .- solar[end,1])):(end-1)], solar_interpolated * 1e13)

1 ./ νs[argmin(abs.((νs ./ 100) .- solar[end,1])):(end-1)]

a = ones(1500)
a[380:(end-1)] .= solar_interpolated

plot(1:1500, a .* Ls)