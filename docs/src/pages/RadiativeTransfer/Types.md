# RadiativeTransfer Module Methods & Types

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