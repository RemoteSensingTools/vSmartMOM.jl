# Absorption Module Example

```julia
using vSmartMOM
using vSmartMOM.Absorption

## 
## STEP 1: Get the Hitran data
## 

# If you have a fixed-width file stored locally
hitran_data = read_hitran("path/to/file", mol=2, iso=1, ν_min=6000, ν_max=6400)

# If you would like to download the file from the hitran database
hitran_data = read_hitran(artifact("CO2"), mol=2, iso=1, ν_min=6000, ν_max=6400)
# For a list of available molecules, use Absorption.show_molecules()

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
profile = absorption_cross_section(model_*, 6000:0.01:6400, 1000.1, 296.1)

# You can obtain the derivatives of the cross-section by setting the autodiff parameter to true: 
profile, derivs = absorption_cross_section(model, 6000:0.01:6400, 1000.1, 296.1, autodiff=true);

```