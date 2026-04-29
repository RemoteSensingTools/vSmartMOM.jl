# Scattering API

`Scattering` computes aerosol optical properties, Fourier/Greek coefficients,
and phase-matrix reconstructions used by the RT core.

## Mie and Aerosol Optics

```@docs
vSmartMOM.Scattering.make_mie_model
vSmartMOM.Scattering.compute_aerosol_optical_properties
vSmartMOM.Scattering.truncate_phase
```

## Output Types

```@docs
vSmartMOM.Scattering.GreekCoefs
vSmartMOM.Scattering.Aerosol
vSmartMOM.Scattering.AerosolOptics
```

## Fourier Decomposition Modes

```@docs
vSmartMOM.Scattering.NAI2
vSmartMOM.Scattering.PCW
```
