# API Reference

This page provides a grouped reference of key functions and types in vSmartMOM.jl.

## High-Level Functions

```@docs
vSmartMOM.parameters_from_yaml
vSmartMOM.read_parameters
vSmartMOM.default_parameters
vSmartMOM.model_from_parameters
vSmartMOM.rt_run
vSmartMOM.rt_run_lin
vSmartMOM.model_from_parameters_lin
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
vSmartMOM.CoreRT.ParameterLayout
```

## Architecture

```@docs
vSmartMOM.CPU
vSmartMOM.GPU
vSmartMOM.default_architecture
vSmartMOM.array_type
```
