# Surfaces

**For:** users selecting a lower-boundary BRDF and developers adding new surface models.

**Next:** [Add a Surface BRDF](../extending/surfaces.md), [Surfaces tutorial](../tutorials/Tutorial_Surfaces.md), [API Reference](../api_reference.md).

Surface models live in `CoreRT` because they are part of the lower-boundary interaction in the adding-doubling solver. Scene files refer to them through safe parser entries such as `LambertianSurfaceScalar(0.15)`.

## Scene File Surface Tags

| Tag | Arguments | Use case | Notes |
| --- | --- | --- | --- |
| `LambertianSurfaceScalar(albedo)` | scalar albedo | Band-wide isotropic surface | Analytical surface Jacobian |
| `LambertianSurfaceSpectrum([a1, a2, ...])` | spectral albedo vector | Per-grid-point isotropic surface | Vector length must match the spectral grid |
| `LambertianSurfaceLegendre([c0, c1, ...])` | Legendre coefficients | Smooth spectral albedo parameterization | Evaluated over the spectral grid |
| `rpvSurfaceScalar(rho0, rho_c, k, Theta)` | RPV empirical BRDF parameters | Scalar land BRDF with hotspot/angular shape | Scalar only |
| `RossLiSurfaceScalar(fvol, fgeo, fiso)` | Ross-Li kernel weights | MODIS/RAMI-style land BRDF | Scalar only |
| `CoxMunkSurface(wind_speed)` | 10 m wind speed in `m/s` | Ocean glint and polarization | Uses built-in water refractive-index lookup unless overridden |
| `CoxMunkSurface(wind_speed=U)` | keyword form | Ocean glint with optional keyword extensions | Supports `n_water`, `whitecap_albedo`, `include_whitecaps`, `shadowing` |

Example:

```yaml
radiative_transfer:
  surface:
    - LambertianSurfaceScalar(0.15)
```

## Canopy Wrapper

Vegetation canopies are configured through an optional top-level `canopy:` section. The parser wraps each band surface in [`CanopySurface`](@ref), so the existing surface becomes the soil BRDF unless a separate canopy soil tag is supplied.

```yaml
canopy:
  LAI: 3.0
  n_layers: 2
  leaf_reflectance: 0.4
  leaf_transmittance: 0.05
  soil: from_surface
```

Internally the canopy is treated as scattering sub-layers above the soil boundary and is combined with the same adding-doubling interaction machinery used by the atmospheric column. The current canopy BRDF path is elastic-only.

## Implementation Map

- Surface types are defined in `src/CoreRT/types.jl`.
- Surface layer construction lives under `src/CoreRT/Surfaces/`.
- Scene-file parsing is registered in `BRDF_MAP` in `src/IO/Parameters.jl`.
- New surface types should add a concrete `AbstractSurfaceType`, a `create_surface_layer!`/`reflectance` implementation as appropriate, and a parser registration.

## Useful APIs

- [`AbstractSurfaceType`](@ref)
- [`LambertianSurfaceScalar`](@ref), [`LambertianSurfaceSpectrum`](@ref), [`LambertianSurfaceLegendre`](@ref), [`LambertianSurfaceSpline`](@ref)
- [`rpvSurfaceScalar`](@ref), [`RossLiSurfaceScalar`](@ref), [`CoxMunkSurface`](@ref)
- [`CanopySurface`](@ref), [`CanopySurface_from_prospect`](@ref), [`create_surface_layer!`](@ref)
