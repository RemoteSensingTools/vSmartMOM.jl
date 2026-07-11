# `canopy` block (optional)

Wraps each spectral band's existing `radiative_transfer.surface` BRDF as
**soil** underneath a multi-layer vegetation canopy (`CoreRT.CanopySurface`),
solved internally via adding-doubling before interacting with the
atmospheric RT. The mere presence of this top-level block is what activates
the wrapping — see [`Schema/surface.md` § Composite (canopy + soil)](surface.md#composite-canopy--soil).

Parsed by `_parse_canopy_section` in
[`src/IO/Parameters.jl`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/src/IO/Parameters.jl).

## Note on schema coverage

The JSON Schema at
[`schemas/vsmartmom-parameters.schema.json`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/schemas/vsmartmom-parameters.schema.json)
currently declares `canopy` as `additionalProperties: true` with only a
block-level description — none of the per-field validation or autocomplete
below is wired into the schema yet. Field names and defaults on this page
come from the parser and the `CanopySurface` constructor, not from the JSON
Schema, so editor autocomplete (see [`Schema.md`](../Schema.md#editor-support--autocomplete--inline-validation))
will not catch a typo'd `canopy:` field name.

## Fields

| Field | Type | Default | Description |
|---|---|---|---|
| `LAI` | `Real` | `3.0` | Total leaf area index. |
| `n_layers` | `Integer` | `1` | Number of canopy sub-layers solved via adding-doubling internally (`1` = big-leaf). |
| `leaf_reflectance` | `Real` or `Vector{Real}` | `0.4` | Leaf reflectance. Scalar or spectral vector on `leaf_optics_grid`. |
| `leaf_transmittance` | `Real` or `Vector{Real}` | `0.05` | Leaf transmittance. Scalar or spectral vector on `leaf_optics_grid`. |
| `leaf_optics_grid` | `Vector{Real}` | `nothing` (scalar leaf optics) | Wavelength/wavenumber grid matching spectral `leaf_reflectance` / `leaf_transmittance`. In YAML this is a plain numeric vector paired with `grid_unit`; built in Julia it can also be a `Unitful` vector (`u"nm"` / `u"cm^-1"`). |
| `grid_unit` | `String` | `"nm"` | Unit of `leaf_optics_grid`: `"nm"` or `"cm_inv"`. |
| `clumping` | `Real`, `String`, or `Dict` | none (`NoClumping`) | Canopy clumping factor. A bare number `Ω` → `ConstantClumping(Ω)`; `"none"` → `NoClumping()`; `{type: "chen_leblanc", Ω₀: .., c: .., e: ..}` → `ChenLeblancClumping`. Affects the effective G-function used in propagation, not the Z-matrix normalization. |
| `include_atm` | `Bool` | `false` | Include within-canopy gas absorption between sub-layers. |
| `canopy_dp` | `Real` [hPa] | `nothing` | Pressure thickness of the within-canopy air column (e.g. `3.0` hPa ≈ 30 m forest). Only meaningful when `include_atm: true`. |
| `soil` | `String` or `"from_surface"` | `"from_surface"` | Soil BRDF underneath the canopy. `"from_surface"` reuses the band's existing `radiative_transfer.surface` entry as soil; otherwise a surface constructor string using the same vocabulary as [`Schema/surface.md`](surface.md). |
| `n_leaf_quadrature` | `Integer` | `CanopyOptics.CanopyQuadrature()` default | Leaf-angle quadrature point count (advanced; forwarded to `CanopyOptics.CanopyQuadrature`). |
| `n_azimuth_quadrature` | `Integer` | `CanopyOptics.CanopyQuadrature()` default | Azimuth quadrature point count (advanced; forwarded to `CanopyOptics.CanopyQuadrature`). |

Not YAML-configurable (Julia-only, set by constructing `CanopySurface`
directly): `LAD` (leaf-angle distribution; defaults to
`CanopyOptics.spherical_leaves()`), `canopy_scattering` (defaults to a
`BiLambertianCanopyScattering` built from the mean `leaf_reflectance` /
`leaf_transmittance`), and `lai_fractions` (per-sub-layer LAI split; defaults
to uniform). For PROSPECT-derived spectral leaf optics, build the surface in
Julia with `CanopySurface_from_prospect(leaf, λgrid; soil=..., LAI=...)`
instead of the YAML `canopy:` block (see the comment block in
`config/canopy_forest.yaml`).

## Example

Forest canopy over a dark Lambertian soil, with within-canopy absorption
(`config/canopy_forest.yaml`):

```yaml
radiative_transfer:
  spec_bands:
    - "[19417.0 19418.0]"
  surface:
    - LambertianSurfaceScalar(0.10)   # soil BRDF, wrapped by canopy below
  polarization_type: Stokes_IQUV()
  nstreams: 11
  truncation: NoTruncation()
  depol: -1
  float_type: Float64
  architecture: default_architecture

canopy:
  LAI: 5.0
  n_layers: 4
  leaf_reflectance: 0.45
  leaf_transmittance: 0.05
  include_atm: true
  canopy_dp: 3.0        # ~30 m forest ≈ 3 hPa
```

Full runnable scene:
[`config/canopy_forest.yaml`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/config/canopy_forest.yaml).
Minimal test scene:
[`test/test_parameters/CanopyTest.yaml`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/test/test_parameters/CanopyTest.yaml).

## See also

- [`Schema/surface.md`](surface.md#composite-canopy--soil) — the `surface:`
  field this block wraps, and the Phase C `component_m_max` trait for
  `CanopySurface`
- [`docs/src/pages/api/surfaces.md`](../../api/surfaces.md) —
  `CanopySurface` / `CanopySurface_from_prospect` docstrings
- [`docs/src/pages/tutorials/Tutorial_Canopy.md`](../../tutorials/Tutorial_Canopy.md)
  — runnable canopy tutorial
