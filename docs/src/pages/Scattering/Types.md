# Scattering Module Methods & Types

This page groups the most commonly used APIs in a practical order.

## Theory Mapping (Quick)

Reference: [Sanghavi (2014)](https://doi.org/10.1016/j.jqsrt.2013.12.015).

- `compute_aerosol_optical_properties` (NAI2): Eq. (1) + Fourier framework around Eq. (17)
- `compute_aerosol_optical_properties` (PCW): Eq. (22) -> Eq. (24)
- `reconstruct_phase`: Greek/Fourier reconstruction around Eq. (17)
- `phase_function`, `compute_aerosol_XS`, `compute_ref_aerosol_extinction`: Eq. (1)-based cross-sections

For the expanded mapping table and implementation notes, see the [Scattering overview](Overview.md).

## High-level Workflow

### 1) Define model inputs

```@docs
make_mie_model
```

### 2) Compute optical properties

```@docs
compute_aerosol_optical_properties
```

### 3) Reconstruct phase matrix elements

```@docs
reconstruct_phase
```

## Convenience Computations

These are useful if you need cross-sections or scalar phase-function outputs without running the full RT-facing workflow.

```@docs
phase_function
compute_aerosol_XS
compute_ref_aerosol_extinction
truncate_phase
```

## Wigner Utilities (PCW)

```@docs
compute_wigner_values
save_wigner_values
load_wigner_values
```

## Types

### Core abstract interfaces

```@docs
AbstractAerosolType
AbstractFourierDecompositionType
AbstractPolarizationType
AbstractTruncationType
```

### Aerosol and model types

```@docs
Aerosol
MieModel
```

### Fourier decomposition types

```@docs
NAI2
PCW
```

### Polarization types

```@docs
Stokes_IQUV
Stokes_IQU
Stokes_I
```

### Truncation type

```@docs
δBGE
```

### Output types

```@docs
GreekCoefs
AerosolOptics
```

### Linearized output types

```@docs
linGreekCoefs
linAerosolOptics
```
