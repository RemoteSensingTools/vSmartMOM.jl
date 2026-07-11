# Scattering API

`Scattering` computes aerosol optical properties, Fourier/Greek coefficients,
and phase-matrix reconstructions used by the RT core.

## Mie and Aerosol Optics

```@docs
vSmartMOM.Scattering.make_mie_model
vSmartMOM.Scattering.compute_aerosol_optical_properties
vSmartMOM.Scattering.truncate_phase
```

## Phase Matrix Reconstruction

```@docs
vSmartMOM.Scattering.reconstruct_phase
```

## Convenience Computations

Useful for cross-sections or scalar phase-function outputs without running
the full RT-facing workflow.

```@docs
vSmartMOM.Scattering.phase_function
vSmartMOM.Scattering.compute_aerosol_XS
vSmartMOM.Scattering.compute_ref_aerosol_extinction
```

## Analytic Phase Functions

Lightweight scattering sources for idealized tests, sensitivity experiments,
and StandaloneSS validation. Not a separate RT implementation: each analytic
source is projected into `GreekCoefs` and then consumed by the same CoreRT
optical-property machinery as Mie-derived aerosols. `phase_matrix_first_column`
evaluates the single-scatter phase-matrix column needed for an unpolarized
direct solar beam at one exact sun-view geometry — the shared hook used by
StandaloneSS vector paths for the sun-to-atmosphere-to-sensor term.

```@docs
vSmartMOM.Scattering.AbstractAnalyticPhaseFunction
vSmartMOM.Scattering.HenyeyGreensteinPhaseFunction
vSmartMOM.Scattering.SyntheticPolarizedHenyeyGreensteinPhaseFunction
vSmartMOM.Scattering.greek_coefficients
vSmartMOM.Scattering.analytic_aerosol_optics
vSmartMOM.Scattering.phase_matrix_first_column
```

## Wigner Utilities (PCW)

```@docs
vSmartMOM.Scattering.compute_wigner_values
vSmartMOM.Scattering.save_wigner_values
vSmartMOM.Scattering.load_wigner_values
```

## Core Abstract Interfaces

```@docs
vSmartMOM.Scattering.AbstractAerosolType
vSmartMOM.Scattering.AbstractFourierDecompositionType
vSmartMOM.Scattering.AbstractPolarizationType
vSmartMOM.Scattering.AbstractTruncationType
```

## Output Types

```@docs
vSmartMOM.Scattering.GreekCoefs
vSmartMOM.Scattering.Aerosol
vSmartMOM.Scattering.AerosolOptics
vSmartMOM.Scattering.MieModel
```

## Linearized Output Types

```@docs
vSmartMOM.Scattering.linGreekCoefs
vSmartMOM.Scattering.linAerosolOptics
```

## Fourier Decomposition Modes

```@docs
vSmartMOM.Scattering.NAI2
vSmartMOM.Scattering.PCW
```

## Polarization Types

```@docs
vSmartMOM.Scattering.Stokes_IQUV
vSmartMOM.Scattering.Stokes_IQU
vSmartMOM.Scattering.Stokes_IQ
vSmartMOM.Scattering.Stokes_I
```

## Phase-function Truncation

```@docs
vSmartMOM.Scattering.AutoTruncation
vSmartMOM.Scattering.NoTruncation
vSmartMOM.Scattering.δBGE
```
