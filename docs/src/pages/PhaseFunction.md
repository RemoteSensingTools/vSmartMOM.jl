# Phase-Function Module

## Introduction

This module enables scattering phase-function calculation of atmospheric aerosols with different size distributions, incident wavelengths, and refractive indices. It can perform the calculation using either the Siewert NAI-2 or Domke PCW methods ([Suniti Sanghavi 2014](https://www.sciencedirect.com/science/article/pii/S0022407313004962)). 

You can calculate a phase function in three steps: 

1. Use [`make_univariate_aerosol`](@ref) to create an aerosol with selected distribution and properties
2. Use [`make_mie_model`](@ref) to set up all calculation parameters
3. Use [`compute_aerosol_optical_properties`](@ref) to perform the phase-function calculation using the defined model settings

## Example

```julia
using RadiativeTransfer
using RadiativeTransfer.PhaseFunction
using Distributions
using BenchmarkTools

## 
## STEP 1: Create the Aerosol
## 

# Aerosol particle distribution and properties 
μ  = 0.3                # Log mean radius
σ  = 6.82               # Log stddev of radius
r_max = 30.0            # Maximum radius
nquad_radius = 2500     # Number of quadrature points for integrating of size dist.
nᵣ = 1.3                # Real part of refractive index
nᵢ = 0.0                # Imag part of refractive index

size_distribution = LogNormal(log(μ), log(σ))

# Create the aerosol
aero = make_univariate_aerosol(size_distribution, r_max, nquad_radius, nᵣ, nᵢ)

## 
## STEP 2: Create the Mie Calculations model
## 

λ = 0.55   # Incident wavelength
polarization_type = Stokes_IQUV()
truncation_type = δBGE(10, 2)

# NAI2 Method
model_NAI2 = make_mie_model(NAI2(), aero, λ, polarization_type, truncation_type)

# PCW Method with saved/loaded Wigner tables
wigner_file_path = "PATH_TO_SAVED_WIGNER_MATRIX"
model_PCW = make_mie_model(PCW(), aero, λ, polarization_type, truncation_type, wigner_file_path)

# PCW Method with newly computed Wigner tables
wigner_A, wigner_B = compute_wigner_values(600) # Specify N_max
model_PCW_computed_wigner = make_mie_model(PCW(), aero, λ, polarization_type, truncation_type, wigner_A, wigner_B)

## 
## STEP 3: Perform the Mie Calculations
## 

aerosol_optics_NAI2 = compute_aerosol_optical_properties(model_NAI2);
aerosol_optics_PCW = compute_aerosol_optical_properties(model_PCW);

```

## Create the Aerosol

```@docs
make_univariate_aerosol
```

## Creating/Saving Wigner Matrices (if using PCW method)

```@docs
compute_wigner_values
save_wigner_values
load_wigner_values
```

## Defining the Mie Computation Model

```@docs
make_mie_model
```

## Computing Aerosol Optical Properties

```@docs
compute_aerosol_optical_properties
```

## Types

### Aerosol Types

```@docs
UnivariateAerosol
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