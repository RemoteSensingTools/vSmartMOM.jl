using Revise
using RadiativeTransfer
using RadiativeTransfer.Absorption
using RadiativeTransfer.Scattering
using RadiativeTransfer.vSmartMOM

# Load parameters from file
parameters = vSmartMOM.parameters_from_yaml("RadiativeTransfer/src/vSmartMOM/ModelParameters/DefaultParameters.yaml")

# Sets all the "specific" parameters
# parameters = vSmartMOM.default_parameters();

# Generates all the derived attributes from above parameters
 model = model_from_parameters(parameters);

# Run RT using the computed model
@time R_GPU, T_GPU = vSmartMOM.rt_run(model);

