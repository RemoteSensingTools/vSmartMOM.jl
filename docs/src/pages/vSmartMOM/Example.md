# vSmartMOM Module Example

## Forward RT

```julia
using vSmartMOM
using vSmartMOM.CoreRT

## STEP 1: Load / Customize RT parameters

# Load a default set of parameters
parameters = default_parameters();

# Or load from a YAML/TOML file (see format at:
# https://github.com/remotesensingtools/vSmartMOM.jl/blob/main/src/CoreRT/DefaultParameters.yaml)
# parameters = read_parameters("path/to/your/params.yaml")

# You can modify any field in parameters (it is a mutable struct):
# parameters.max_m = 5
# See the parameters guide for field descriptions:
# https://remotesensingtools.github.io/vSmartMOM.jl/dev/pages/vSmartMOM/InputParametersGuide/

## STEP 2: Create a model from parameters

# This computes all derived fields (absorption profiles, scattering phase functions,
# quadrature points, etc.) and returns an RTModel:
model = model_from_parameters(parameters);

# The RTModel is a hierarchical struct with sub-structs:
#   model.solver      — SolverConfig (polarization, quadrature, truncation)
#   model.geometry    — ObsGeometry (SZA, VZA, VAZ)
#   model.quad_points — QuadPoints (μ, weights)
#   model.atmosphere  — Atmosphere (profile, spec_bands)
#   model.optics      — Optics (rayleigh, aerosol, τ_abs, τ_rayl)
#   model.surfaces    — Vector{AbstractSurfaceType}

## STEP 3: Run the RT simulation

R, T = rt_run(model)
# R: TOA reflectance array (n_stokes × n_vza × n_spectral)
# T: BOA transmittance array (same shape)
```

## Linearized RT (Jacobians)

```julia
using vSmartMOM

parameters = default_parameters();

# Build both forward and linearized models:
model, lin_model = model_from_parameters(LinMode(), parameters);

# Determine Jacobian dimensions:
NAer  = length(parameters.scattering_params.rt_aerosols)
NGas  = size(lin_model.tau_dot_abs[1], 1)
NSurf = 1

# Run linearized RT:
R, T, dR, dT = rt_run(model, lin_model, NAer, NGas, NSurf)
# dR, dT contain analytic Jacobians with respect to aerosol, gas, and surface parameters
```

See the [Jacobians tutorial](../tutorials/Tutorial_Jacobians.md) for a detailed walkthrough of the Jacobian layout.
