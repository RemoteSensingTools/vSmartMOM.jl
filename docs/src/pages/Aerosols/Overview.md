# Aerosols

**For:** users reading aerosol inputs and developers stabilizing aerosol scheme support.

**Next:** [Configure a Scene](../IO/Overview.md), [Mie & Rayleigh (Concepts)](../concepts/03b_scattering.md), [Library](../api_reference.md).

The `Aerosols` module is a user-facing framework for aerosol input data, refractive-index lookup tables, and scheme-specific optical-property preparation. This API is still being stabilized after the `sanghavi-unified` merge, so treat it as available but more likely to evolve than the CoreRT, Absorption, and Scattering APIs.

## Supported Schemes

| Scheme | Type | Input style | Use case |
| --- | --- | --- | --- |
| TOMAS-15 | `TOMAS15Scheme` | Size-resolved concentrations in 15 logarithmic diameter bins | Microphysics-aware aerosol columns |
| Two-moment | `TwoMomentScheme` | AOD, effective radius, and fixed geometric width per species | Compact aerosol-state inputs from chemistry or retrieval systems |

Both schemes are configured from YAML and read from NetCDF through `read_aerosol_data`. The returned `AerosolData` stores the scheme, species data, coordinates, and NetCDF metadata.

## Refractive Indices

The refractive-index database is separate from the aerosol state:

```julia
using vSmartMOM

db = vSmartMOM.Aerosols.load_refractive_index_database("refractive_indices.yaml")
n = vSmartMOM.Aerosols.get_refractive_index(db, "sulfate_suso", 0.55)
```

Wavelengths are in micrometers. `get_refractive_index` interpolates the real and imaginary parts linearly and returns a complex refractive index.

## Relationship To Scattering

`Aerosols` handles data ingestion and scheme organization. The lower-level Mie and Greek-coefficient calculations are explained in [Mie & Rayleigh (Concepts)](../concepts/03b_scattering.md). In a complete workflow, aerosol input data are converted into optical properties that the CoreRT layer assembly can consume.

## Current Caveats

- The module is public but still WIP; names and configuration details may be cleaned up before registration.
- TOMAS-15 and two-moment support are present; additional schemes should be added through new `AerosolScheme` subtypes rather than one large parser branch.
- Heavy NetCDF fixtures are not part of the unit-test path.

## Useful APIs

- `AerosolScheme`, `TOMAS15Scheme`, `TwoMomentScheme`
- `AerosolData`, `AerosolSpeciesData`
- `RefractiveIndexDatabase`, `RefractiveIndexLUT`
- `read_aerosol_data`, `load_refractive_index_database`, `get_refractive_index`
- `list_species`, `wavelength_range`
