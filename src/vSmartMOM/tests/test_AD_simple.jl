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
parameters = vSmartMOM.parameters_from_yaml("test/helper/O2Parameters.yaml")
FT = Float64
# default_parameters

# Sets all the "specific" parameters
# parameters = vSmartMOM.default_parameters();
function runner(x, parameters=parameters)
    parameters.τAer_ref = [x[1]];
    parameters.p₀ = [x[2]];
    parameters.nᵣ = [x[3]];

    model = model_from_parameters(parameters);
    
    model.params.architecture = RadiativeTransfer.Architectures.GPU()
    R = vSmartMOM.rt_run(model, i_band=1);
    #@show J
    return R[1,1,:]#; R_SFI[1,1,:]
    #return R.τ_λ_all
end

#x = [0.1,90001.0]
x = FT[0.1,90001.0,1.3]
# Run FW model:
@time runner(x);
@time dfdx = ForwardDiff.jacobian(runner, x);