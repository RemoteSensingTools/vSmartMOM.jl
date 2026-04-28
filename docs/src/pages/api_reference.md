# API Reference

This page provides a grouped reference of key functions and types in vSmartMOM.jl.

## High-Level Functions

```@docs
vSmartMOM.read_parameters
vSmartMOM.parameters_from_file
vSmartMOM.parameters_from_source
vSmartMOM.parameters_from_dict
vSmartMOM.parameters_from_yaml
vSmartMOM.default_parameters
vSmartMOM.model_from_parameters
vSmartMOM.rt_run
```

## Linearized RT

```@docs
vSmartMOM.model_from_parameters_lin
vSmartMOM.rt_run_lin
vSmartMOM.CoreRT.ParameterLayout
vSmartMOM.CoreRT.n_total
vSmartMOM.CoreRT.aerosol_range
vSmartMOM.CoreRT.gas_range
vSmartMOM.CoreRT.surface_range
vSmartMOM.CoreRT.surface_index
vSmartMOM.CoreRT.n_layer_params
vSmartMOM.CoreRT.canopy_range
```

## HITRAN Data Management

```@docs
vSmartMOM.artifact
vSmartMOM.fetch_hitran
vSmartMOM.fetch_hitran_by_ids
vSmartMOM.set_hitran_edition!
vSmartMOM.get_hitran_edition
vSmartMOM.available_hitran_editions
vSmartMOM.hitran_info
vSmartMOM.hitran_is_cached
```

## Absorption

```@docs
vSmartMOM.Absorption.compute_absorption_cross_section
vSmartMOM.Absorption.absorption_cross_section
vSmartMOM.Absorption.read_hitran
vSmartMOM.Absorption.make_hitran_model
vSmartMOM.Absorption.make_interpolation_model
vSmartMOM.Absorption.AbstractBroadeningFunction
vSmartMOM.Absorption.Voigt
vSmartMOM.Absorption.Lorentz
vSmartMOM.Absorption.Doppler
```

## Scattering

```@docs
vSmartMOM.Scattering.make_mie_model
vSmartMOM.Scattering.compute_aerosol_optical_properties
vSmartMOM.Scattering.truncate_phase
vSmartMOM.Scattering.GreekCoefs
vSmartMOM.Scattering.Aerosol
vSmartMOM.Scattering.AerosolOptics
vSmartMOM.Scattering.NAI2
vSmartMOM.Scattering.PCW
```

## Inelastic Scattering

```@docs
vSmartMOM.InelasticScattering.AbstractRamanType
vSmartMOM.InelasticScattering.noRS
vSmartMOM.InelasticScattering.noRS_plus
vSmartMOM.InelasticScattering.RRS
vSmartMOM.InelasticScattering.RRS_plus
vSmartMOM.InelasticScattering.VS_0to1
vSmartMOM.InelasticScattering.VS_1to0
vSmartMOM.InelasticScattering.VS_0to1_plus
vSmartMOM.InelasticScattering.VS_1to0_plus
vSmartMOM.InelasticScattering.compute_γ_air_Cabannes!
vSmartMOM.InelasticScattering.compute_γ_air_Rayleigh!
vSmartMOM.InelasticScattering.compute_σ_Rayl_coeff!
vSmartMOM.InelasticScattering.compute_σ_Rayl_VibRaman_coeff_hires!
vSmartMOM.InelasticScattering.compute_σ_VibRaman_coeff!
vSmartMOM.InelasticScattering.apply_lineshape!
vSmartMOM.InelasticScattering.compute_energy_levels!
vSmartMOM.InelasticScattering.compute_stellar_Rayl
vSmartMOM.InelasticScattering.get_greek_raman
```

## Aerosols

```@docs
vSmartMOM.Aerosols.AerosolScheme
vSmartMOM.Aerosols.TOMAS15Scheme
vSmartMOM.Aerosols.TwoMomentScheme
vSmartMOM.Aerosols.AerosolSpeciesData
vSmartMOM.Aerosols.AerosolData
vSmartMOM.Aerosols.RefractiveIndexLUT
vSmartMOM.Aerosols.RefractiveIndexDatabase
vSmartMOM.Aerosols.read_aerosol_data
vSmartMOM.Aerosols.load_refractive_index_database
vSmartMOM.Aerosols.get_refractive_index
vSmartMOM.Aerosols.list_species
vSmartMOM.Aerosols.wavelength_range
```

## SolarModel

```@docs
vSmartMOM.SolarModel.planck_spectrum_wn
vSmartMOM.SolarModel.planck_spectrum_wl
vSmartMOM.SolarModel.solar_transmission_from_file
vSmartMOM.SolarModel.default_solar_transmission_path
vSmartMOM.SolarModel.default_solar_transmission
vSmartMOM.SolarModel.default_solar_spectrum_at_earth
```

## Surface Models

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
vSmartMOM.CoreRT.create_surface_layer!
vSmartMOM.CoreRT.water_refractive_index
vSmartMOM.CoreRT.fresnel_coefficients
vSmartMOM.CoreRT.fresnel_mueller
vSmartMOM.CoreRT.stokes_rotation_matrix
```

## Types

```@docs
vSmartMOM.CoreRT.vSmartMOM_Parameters
vSmartMOM.CoreRT.RTModel
vSmartMOM.CoreRT.AtmosphericProfile
vSmartMOM.CoreRT.ObsGeometry
vSmartMOM.CoreRT.QuadPoints
vSmartMOM.CoreRT.CompositeLayer
vSmartMOM.CoreRT.AddedLayer
```

## Architecture

```@docs
vSmartMOM.CPU
vSmartMOM.GPU
vSmartMOM.default_architecture
vSmartMOM.array_type
```
