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
using DelimitedFiles
using Interact
using Interpolations
using InstrumentOperator
using ForwardDiff

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

#model.params.architecture = RadiativeTransfer.Architectures.GPU();

function run_auto(x)
    println(x)
    println(x[1] isa ForwardDiff.Dual)
    println(x[1])
    println(typeof(parameters.μ))
    parameters.μ  = [x[1]]
    parameters.σ  = [x[2]]
    parameters.nᵣ = [x[3]]
    parameters.nᵢ = [x[4]]
    model = model_from_parameters(parameters);
    # model.architecture = RadiativeTransfer.Architectures.CPU()
    model.params.architecture = RadiativeTransfer.Architectures.CPU()
    R = vSmartMOM.rt_run(model, i_band=1)
    return R[1,1,:]
end

x = [1.3, 2.0, 1.3, 0.000001];

# @time run_auto(x)
@time dfdx = ForwardDiff.jacobian(run_auto, x);

a = 10

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