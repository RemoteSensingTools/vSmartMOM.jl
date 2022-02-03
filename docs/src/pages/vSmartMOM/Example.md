# vSmartMOM Module Example

```julia
using vSmartMOM
using vSmartMOM.CoreRT

## 
## STEP 1: Load / Customize RT parameters
## 

# If you would like to load your own parameters from a YAML file
# (See required format at: https://github.com/remotesensingtools/vSmartMOM.jl/blob/main/src/CoreRT/ModelParameters/DefaultParameters.yaml)
parameters = parameters_from_yaml("vSmartMOM/src/CoreRT/DefaultParameters.yaml")

# OR if you would like to load a default set of parameters
parameters = default_parameters();

# You can then change any individual fields in parameters (parameters.field = ...)
# Please see descriptions of each field here: 
# https://remotesensingtools.github.io/vSmartMOM.jl/dev/pages/vSmartMOM/InputParametersGuide.html

## 
## STEP 2: Create a model from parameters
## 

# Creating this model will use the parameters object to calculate derived fields and attributes
model = model_from_parameters(parameters);

# Again, you can then change any derived field in model (model.field = ...)

## 
## STEP 3: Run the RT simulation
## 

R = rt_run(model)

```