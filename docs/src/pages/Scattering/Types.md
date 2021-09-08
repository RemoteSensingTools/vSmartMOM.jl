# Scattering Module Methods & Types 

## Defining the Mie Computation Model

```@docs
make_mie_model
```

## Creating/Saving Wigner Matrices (if using PCW method)

```@docs
compute_wigner_values
save_wigner_values
load_wigner_values
```

## Computing Aerosol Optical Properties

```@docs
compute_aerosol_optical_properties
```

## Reconstructing Phase Function

```@docs
reconstruct_phase
```

## Types

### Aerosol Types

```@docs
Aerosol
```

### Fourier Decomposition Computation Types

```@docs
NAI2
PCW
```

### Polarization Types 

```@docs
Stokes_IQUV
Stokes_IQU
Stokes_I
```

### Truncation Types 

```@docs
δBGE
```

### Mie Computation Model Type

```@docs
MieModel
```

### Output Aerosol Optics Types 

```@docs
GreekCoefs
AerosolOptics
```