# IO API

The IO layer validates scene inputs and dispatches loading by source type:
ordinary files, dictionaries, NetCDF sources, and GEOS-Chem products.

## Atmospheric Profiles

```@docs
vSmartMOM.read_atmos_profile
vSmartMOM.read_atmos_profile_dict
```

## Typed Sources

```@docs
vSmartMOM.GeosChemSource
vSmartMOM.NetCDFGridSource
vSmartMOM.NetCDFSource
vSmartMOM.geoschem_to_dict
vSmartMOM.read_geoschem_profile
```
