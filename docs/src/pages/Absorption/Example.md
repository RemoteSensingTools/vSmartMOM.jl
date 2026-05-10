# Absorption Module Example

```julia
using vSmartMOM
using vSmartMOM.Absorption

##
## STEP 1: Get the HITRAN data
##

# --- Option A: Legacy artifacts (HITRAN 2016, the default) ---
# The artifact is downloaded and cached automatically on first use.
hitran_data = read_hitran(artifact("CO2"), mol=2, ν_min=6000, ν_max=6400)

# --- Option B: Direct download from hitran.org (HITRAN 2024) ---
# Switch the active HITRAN edition, then use artifact() as before:
set_hitran_edition!("HITRAN2024")
hitran_data = read_hitran(artifact("CO2"), mol=2, ν_min=6000, ν_max=6400)

# Or download with a specific wavenumber range (avoids fetching the full database):
path = fetch_hitran("CO2"; numin=6000, numax=6400, edition="HITRAN2024")
hitran_data = read_hitran(path, mol=2, ν_min=6000, ν_max=6400)

# --- Option C: Load from a local .par file ---
# hitran_data = read_hitran("path/to/file.par", mol=2, ν_min=6000, ν_max=6400)

# By default, all isotopologues are included. To select only the main isotopologue:
# hitran_data = read_hitran(artifact("CO2"), mol=2, iso=1, ν_min=6000, ν_max=6400)

# For a list of available molecules, use:
# Absorption.show_molecules()
# Absorption.search_molecules("C")

##
## STEP 2: Create a model from parameters
##

# These are some example models, but you can create/customize your model however you'd like.
# Please see make_hitran_model documentation for optional arguments
model_doppler = make_hitran_model(hitran_data, Doppler())
model_lorentz = make_hitran_model(hitran_data, Lorentz())
model_voigt   = make_hitran_model(hitran_data, Voigt())

# GPU variant:
# model_voigt_GPU = make_hitran_model(hitran_data, Voigt(), architecture=Architectures.GPU())

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

# Compute cross-section at a given pressure (hPa) and temperature (K):
σ = absorption_cross_section(model_voigt, 6000:0.01:6400, 1013.0, 296.0)
```

See the [Absorption tutorial](../tutorials/Tutorial_Absorption.md) for a detailed walkthrough with plots.

See the [HITRAN Data Management](HITRAN_Data.md) page for details on edition switching, version tracking, and the download API.
