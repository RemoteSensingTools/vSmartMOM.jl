# Aerosols API

The aerosol input layer is available, but its high-level API is still being
stabilized after the unified-branch merge.

## Schemes and Data Containers

```@docs
vSmartMOM.Aerosols.AerosolScheme
vSmartMOM.Aerosols.TOMASScheme
vSmartMOM.Aerosols.TOMAS15Scheme
vSmartMOM.Aerosols.TwoMomentScheme
vSmartMOM.Aerosols.AerosolSpeciesData
vSmartMOM.Aerosols.AerosolData
vSmartMOM.Aerosols.SectionalAerosolData
vSmartMOM.Aerosols.LogNormalFit
vSmartMOM.Aerosols.ConstantIntegrationPerBin
vSmartMOM.Aerosols.LinearIntegrationPerBin
vSmartMOM.Aerosols.DirectBinSum
```

## Refractive-Index Tables

```@docs
vSmartMOM.Aerosols.RefractiveIndexLUT
vSmartMOM.Aerosols.RefractiveIndexDatabase
vSmartMOM.Aerosols.read_aerosol_data
vSmartMOM.Aerosols.read_tomas
vSmartMOM.Aerosols.load_refractive_index_database
vSmartMOM.Aerosols.get_refractive_index
vSmartMOM.Aerosols.list_species
vSmartMOM.Aerosols.wavelength_range
vSmartMOM.Aerosols.compute_column_aod
```
