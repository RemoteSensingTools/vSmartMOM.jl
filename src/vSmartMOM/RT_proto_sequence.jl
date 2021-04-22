using Revise
# using Plots
using Pkg
# Pkg.activate(".")
using BenchmarkTools
using RadiativeTransfer
using RadiativeTransfer.Absorption
using RadiativeTransfer.Scattering
using RadiativeTransfer.vSmartMOM

# Load parameters from file
# parameters = vSmartMOM.parameters_from_yaml("RadiativeTransfer/src/vSmartMOM/ModelParameters/DefaultParameters.yaml")

# params_dict = YAML.load_file("RadiativeTransfer/src/vSmartMOM/ModelParameters/DefaultParameters.yaml")

# vSmartMOM.validate_yaml_parameters(params_dict)

# file = "RadiativeTransfer/src/vSmartMOM/ModelParameters/SampleProfile.yaml"
# profile_nc4 = vSmartMOM.read_atmos_profile(parameters.file, parameters.lat, parameters.lon, parameters.timeIndex);
# profile_yaml = vSmartMOM.read_atmos_profile(file);

# Sets all the "specific" parameters
parameters = vSmartMOM.default_parameters();

# Generates all the derived attributes from above parameters
model = model_from_parameters(parameters);
model.params.architecture = RadiativeTransfer.Architectures.GPU();

R = vSmartMOM.rt_run(model)
