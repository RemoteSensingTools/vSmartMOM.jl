using Revise
using Plots
using Pkg
# Pkg.activate(".")
using RadiativeTransfer
using RadiativeTransfer.Absorption
using RadiativeTransfer.Scattering
using RadiativeTransfer.vSmartMOM

# Load parameters from file
parameters = vSmartMOM.parameters_from_yaml("src/vSmartMOM/ModelParameters/DefaultParameters.yaml")

# default_parameters

# Sets all the "specific" parameters
# parameters = vSmartMOM.default_parameters();

# Generates all the derived attributes from above parameters
model = model_from_parameters(parameters);
#model.params.
model.params.architecture = RadiativeTransfer.Architectures.GPU();
@time R = vSmartMOM.rt_run(model);


