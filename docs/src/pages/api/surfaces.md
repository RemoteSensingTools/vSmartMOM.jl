# Surfaces API

Surface models are implemented in `CoreRT` because they form the lower boundary
of the adding-doubling column.

## Surface Types

```@docs
vSmartMOM.CoreRT.AbstractSurfaceType
vSmartMOM.CoreRT.LambertianSurfaceScalar
vSmartMOM.CoreRT.LambertianSurfaceSpectrum
vSmartMOM.CoreRT.LambertianSurfaceLegendre
vSmartMOM.CoreRT.LambertianSurfaceSpline
vSmartMOM.CoreRT.rpvSurfaceScalar
vSmartMOM.CoreRT.RossLiSurfaceScalar
vSmartMOM.CoreRT.CoxMunkSurface
vSmartMOM.CoreRT.CanopySurface
vSmartMOM.CoreRT.CanopySurface_from_prospect
```

## Surface Construction and Optics Helpers

```@docs
vSmartMOM.CoreRT.create_surface_layer!
vSmartMOM.CoreRT.water_refractive_index
vSmartMOM.CoreRT.fresnel_coefficients
vSmartMOM.CoreRT.fresnel_mueller
vSmartMOM.CoreRT.stokes_rotation_matrix
```
