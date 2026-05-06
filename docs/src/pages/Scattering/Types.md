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

Use `make_mie_model` to combine aerosol microphysics, wavelength, polarization,
truncation, and quadrature settings into a `MieModel`.

### 2) Compute optical properties

Use `compute_aerosol_optical_properties` to produce RT-facing aerosol optical
properties: Fourier/Greek coefficients, single-scattering albedo, extinction,
and truncation information.

For analytic test or parameterized aerosols, use
`HenyeyGreensteinPhaseFunction` or
`SyntheticPolarizedHenyeyGreensteinPhaseFunction` with
`analytic_aerosol_optics`. This produces the same `AerosolOptics` container as
the Mie path and therefore feeds the normal CoreRT and StandaloneSS layer
optics.

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
```

## Analytic Phase Functions

Analytic phase-function aerosols are lightweight scattering sources for
idealized tests, sensitivity experiments, and StandaloneSS validation. They are
not a separate RT implementation: each analytic source is projected into
`GreekCoefs` and then consumed by the same CoreRT optical-property machinery as
Mie-derived aerosols.

```@docs
AbstractAnalyticPhaseFunction
HenyeyGreensteinPhaseFunction
SyntheticPolarizedHenyeyGreensteinPhaseFunction
greek_coefficients
analytic_aerosol_optics
phase_matrix_first_column
```

`phase_matrix_first_column` evaluates the single-scatter phase-matrix column
needed for an unpolarized direct solar beam at one exact sun-view geometry. It
is the shared hook used by StandaloneSS vector paths for the
sun-to-atmosphere-to-sensor term.

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
MieModel
```

`Aerosol` defines a size distribution and complex refractive index. See the
[Library](../api_reference.md) for its canonical docstring.

### Fourier decomposition types

The supported Fourier decomposition families are `NAI2` and `PCW`.

### Polarization types

```@docs
Stokes_IQUV
Stokes_IQU
Stokes_I
```

### Truncation type

```@docs
NoTruncation
δBGE
```

### Output types

`GreekCoefs` stores the Fourier/Greek coefficient representation consumed by
the RT kernels, and `AerosolOptics` groups those coefficients with extinction,
single-scattering albedo, and truncation outputs.

### Linearized output types

```@docs
linGreekCoefs
linAerosolOptics
```
