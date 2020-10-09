# Cross-Section Module

## Introduction

This module enables absorption cross-section calculations of atmospheric gases at different pressures, temperatures, and broadeners (Doppler, Lorentzian, Voigt). It uses the [HITRAN](https://hitran.org) energy transition database for calculations. While it enables lineshape calculations from scratch, it also allows users to create and save an interpolator object at specified wavelength, pressure, and temperature grids. It can perform these computations either on CPU or GPU (CUDA).

You can calculate an absorption cross-section in three steps: 

1. Use [`read_hitran`](@ref) to read in a HITRAN database file
2. Use [`make_hitran_model`](@ref) or [`make_interpolation_model`](@ref) to set up model parameters for the cross-section calculation
3. Use [`absorption_cross_section`](@ref) to perform the cross-section calculation using the defined model settings

## Example

```julia
using RadiativeTransfer
using RadiativeTransfer.CrossSection

## 
## STEP 1: Get the Hitran data
## 

# If you have a fixed-width file stored locally
hitran_data = read_hitran("path/to/file", mol=2, iso=1, ν_min=6000, ν_max=6400)

# If you would like to download the file from the hitran database
hitran_data = read_hitran(artifact"hitran_molec_id_2_CO2", mol=2, iso=1, ν_min=6000, ν_max=6400)

## 
## STEP 2: Create a model from parameters
## 

# These are some example models, but you can create/customize your model however you'd like. 
# Please see make_hitran_model documentation for optional arguments
model_doppler = make_hitran_model(hitran_data, Doppler())
model_lorentz = make_hitran_model(hitran_data, Lorentz())
model_voigt_CPU = make_hitran_model(hitran_data, Voigt())
model_voigt_GPU = make_hitran_model(hitran_data, Voigt(), architecture=Architectures.GPU())

# If you would prefer to create an interpolation, and then interpolate the cross-section 
# at other wavelengths/pressures/temperatures, you can use make_interpolation_model as such:

ν_grid = 6000:0.01:6400
pressures = 250:250:1250
temperatures = 100:75:400

model_interp = make_interpolation_model(hitran_data, Voigt(), ν_grid, pressures, temperatures)

# Note: please see make_interpolation_model documentation for optional parameters

## 
## STEP 3: Calculate the absorption cross section with the created model
## 

# You can specify the wavelength grid, pressure, and temperature. 
absorption_cross_section(model_*, 6000:0.01:6400, 1000.1, 296.1)

```

## Reading HITRAN Files

```@docs
read_hitran
```

## Defining Models with Cross-Section Parameters

```@docs
make_hitran_model
make_interpolation_model
```

## Computing Absorption Cross-Sections

```@docs
absorption_cross_section
```

## Types

### Hitran Data Structure Type

```@docs
HitranTable
```

### Broadening Function Types
```@docs
Doppler
Lorentz
Voigt
```

### Complex Error Function Types

```@docs
HumlicekErrorFunction
HumlicekWeidemann32VoigtErrorFunction
HumlicekWeidemann32SDErrorFunction
CPF12ErrorFunction
ErfcHumliErrorFunctionVoigt
ErfcHumliErrorFunctionSD
ErfcErrorFunction
```

### Cross-Section Model Types

```@docs
HitranModel
InterpolationModel
```

