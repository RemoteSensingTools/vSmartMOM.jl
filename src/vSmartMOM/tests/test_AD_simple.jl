using Revise
using Plots
using Pkg
# Pkg.activate(".")
using RadiativeTransfer
using RadiativeTransfer.Absorption
using RadiativeTransfer.Scattering
using RadiativeTransfer.vSmartMOM
using ForwardDiff 

# Load parameters from file
parameters = vSmartMOM.parameters_from_yaml("RadiativeTransfer/test/helper/ThreeBandsParameters.yaml")

# default_parameters

# Sets all the "specific" parameters
# parameters = vSmartMOM.default_parameters();
function runner(x, parameters=parameters)
    parameters.τAer_ref = [x[1]];
    #@show parameters.p₀
    parameters.p₀ = [x[2]];
    parameters.nᵣ = [x[3]];
    @show parameters.p₀
    model = model_from_parameters(parameters);
    
    model.params.architecture = RadiativeTransfer.Architectures.GPU()
    R = vSmartMOM.rt_run(model, i_band=1);
    #@show J
    return R[1,1,:]#; R_SFI[1,1,:]
end

#x = [0.1,90001.0]
x = [0.1,90001.0,1.3]
# Run FW model:
# @time runner(x);
@time dfdx = ForwardDiff.jacobian(runner, x);