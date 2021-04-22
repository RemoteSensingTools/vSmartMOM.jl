# RadiativeTransfer Module

## Introduction

The RadiativeTransfer module allows end-to-end simulation of radiative transfer (RT) throughout Earth's atmosphere and surface. Specifically, it:

1. Enables 1D vectorized plane-parallel RT modeling based on the Matrix Operator Method.
2. Incorporates fast, high fidelity simulations of scattering atmospheres containing haze and clouds – including pressure- and temperature-resolved absorption profiles of gaseous species in the atmosphere. 
3. Enables GPU-accelerated computations of the resulting hyperspectral reflectances/transmittances.

## Example

```julia
using RadiativeTransfer
using RadiativeTransfer.vSmartMOM

## 
## STEP 1: Load / Customize RT parameters
## 

# If you would like to load your own parameters from a YAML file
# (See required format at: https://github.com/RupeshJey/RadiativeTransfer.jl/blob/master/src/vSmartMOM/ModelParameters/DefaultParameters.yaml)
parameters = parameters_from_yaml("RadiativeTransfer/src/vSmartMOM/ModelParameters/DefaultParameters.yaml")

# OR if you would like to load a default set of parameters
parameters = default_parameters();

# You can then change any individual fields in parameters (parameters.field = ...)

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

## References 

```
1. Sanghavi, S., Davis, A.B. and Eldering, A., 2014. vSmartMOM: A vector matrix operator 
method-based radiative transfer model linearized with respect to aerosol properties. Journal of 
Quantitative Spectroscopy and Radiative Transfer, 133, pp.412-433.

2. Sanghavi, S. and Stephens, G., 2015. Adaptation of the delta-m and δ-fit truncation methods 
to vector radiative transfer: Effect of truncation on radiative transfer accuracy. Journal of 
Quantitative Spectroscopy and Radiative Transfer, 159, pp.53-68.

3. Grant, I.P. and Hunt, G.E., 1969. Discrete space theory of radiative transfer I. 
Fundamentals. Proceedings of the Royal Society of London. A. Mathematical and Physical Sciences, 
313(1513), pp. 183-197.

4. Plass, G.N., Kattawar, G.W. and Catchings, F.E., 1973. Matrix operator theory of radiative 
transfer. 1: Rayleigh scattering. Applied Optics, 12(2), pp.314-329.
```

## Defining Parameters for RT Simulation

```@docs
parameters_from_yaml
default_parameters
```

## Using the Parameters Object to Create a Model 

```@docs
model_from_parameters
```

## Performing the RT Simulation

```@docs
rt_run
```

## Types

### Parameters Type

```@docs
vSmartMOM.vSmartMOM_Parameters
```

### Model Type

```@docs
vSmartMOM.vSmartMOM_Model
```

### Surface Types

```@docs
vSmartMOM.AbstractSurfaceType
vSmartMOM.LambertianSurfaceScalar
vSmartMOM.LambertianSurfaceSpectrum
vSmartMOM.LambertianSurfacePolyFit

```

### Quadrature Types
```@docs
vSmartMOM.AbstractQuadratureType
vSmartMOM.RadauQuad
vSmartMOM.GaussQuadHemisphere
vSmartMOM.GaussQuadFullSphere

```

### Atmospheric Profile Type
```@docs
vSmartMOM.AtmosphericProfile
```