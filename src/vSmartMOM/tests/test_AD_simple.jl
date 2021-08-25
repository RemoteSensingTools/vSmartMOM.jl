using Revise
using Plots
using Pkg
# Pkg.activate(".")
using RadiativeTransfer
using RadiativeTransfer.Absorption
using RadiativeTransfer.Scattering
using RadiativeTransfer.vSmartMOM
using ForwardDiff 
using Distributions

# Load parameters from file
parameters = vSmartMOM.parameters_from_yaml("RadiativeTransfer/test/helper/O2Parameters.yaml")
FT = Float64

# Runner is used to set AD fields as duals
function runner(x, parameters=parameters)

    parameters.scattering_params.rt_aerosols[1].τ_ref = x[1];
    parameters.scattering_params.rt_aerosols[1].p₀ = x[2];
    parameters.scattering_params.rt_aerosols[1].aerosol.nᵣ = x[3];
    parameters.scattering_params.rt_aerosols[1].aerosol.nᵢ = x[4];
    parameters.scattering_params.rt_aerosols[1].aerosol.size_distribution = LogNormal(log(x[5]), log(x[6]), check_args=false)

    model = model_from_parameters(parameters);
    R = vSmartMOM.rt_run(model, i_band=1);
    return R[1,1,:]
end

x = FT[0.1,90001.0,1.3, 1.0e-8, 1.3, 2.0]

# Run FW model:
# @time runner(x);
@time dfdx = ForwardDiff.jacobian(runner, x);