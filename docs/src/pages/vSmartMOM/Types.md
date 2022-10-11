# vSmartMOM Module Methods & Types

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
CoreRT.vSmartMOM_Parameters
```

### Model Type

```@docs
CoreRT.vSmartMOM_Model
```

### Surface Types

```@docs
CoreRT.AbstractSurfaceType
CoreRT.LambertianSurfaceScalar
CoreRT.LambertianSurfaceSpectrum

```

### Quadrature Types
```@docs
CoreRT.AbstractQuadratureType
CoreRT.RadauQuad
CoreRT.GaussQuadHemisphere
CoreRT.GaussQuadFullSphere

```

### Atmospheric Profile Type
```@docs
CoreRT.AtmosphericProfile
```