# vSmartMOM Module Methods & Types

## Defining Parameters for RT Simulation

Use [`read_parameters`](@ref) as the general loader for YAML/TOML files, in-memory
dictionaries, and typed IO sources. Use [`default_parameters`](@ref) for the built-in
default scene.

The explicit helpers [`parameters_from_file`](@ref),
[`parameters_from_dict`](@ref), [`parameters_from_source`](@ref), and
[`parameters_from_yaml`](@ref) are also available when a call site should state the
input kind directly.

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
CoreRT.RTModel
```

### Surface Types

```@docs
CoreRT.AbstractSurfaceType
CoreRT.LambertianSurfaceScalar
CoreRT.LambertianSurfaceSpectrum
CoreRT.LambertianSurfaceLegendre
CoreRT.LambertianSurfaceSpline
CoreRT.rpvSurfaceScalar
CoreRT.RossLiSurfaceScalar
CoreRT.CoxMunkSurface
CoreRT.CanopySurface
CoreRT.CanopySurface_from_prospect
CoreRT.invalidate_canopy_cache!
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
