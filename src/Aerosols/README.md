# Aerosol Framework for vSmartMOM

A flexible aerosol module supporting multiple parameterization schemes for radiative transfer calculations.

## Overview

The aerosol framework provides a unified interface for working with different aerosol schemes:

- **TOMAS-15**: Size-resolved microphysics with 15 logarithmic bins (10 nm - 10 μm dry diameter)
- **Two-Moment**: Bulk aerosol optical properties (AOD + effective radius)
- **Extensible**: Easy to add new schemes (MAAM, CARMA, custom, etc.)

## Key Features

- **YAML-driven configuration**: Flexible, user-friendly scheme definitions
- **Wavelength-dependent refractive indices**: Comprehensive database with interpolation
- **Multiple aerosol species**: Dust, sea salt, sulfate, organic carbon, black carbon, water
- **Optical property calculations**: Extinction, scattering, absorption, SSA, phase functions
- **Integration with vSmartMOM RT**: Seamless coupling with radiative transfer core

## Quick Start

### 1. Reading TOMAS-15 Data

```julia
using vSmartMOM
using vSmartMOM.Aerosols

# Read size-resolved aerosol data
data = read_aerosol_data(
    "examples/aerosol_config_tomas15.yaml",
    "GEOSChem.Custom.20190702_0000z.nc4"
)

# Access species concentrations
dust_conc = data.species_data["DUST"].data["concentration"]  # [15 bins × 72 levels]
```

### 2. Reading Two-Moment Data

```julia
# Read bulk aerosol properties
data = read_aerosol_data(
    "examples/aerosol_config_two_moment.yaml",
    "GEOSChem.Aerosols.20190702_0000z.nc4"
)

# Access AOD and effective radius
so4_aod = data.species_data["so4"].data["aod"]
so4_radius = data.species_data["so4"].data["radius"]
```

### 3. Computing Optical Properties

```julia
# Load refractive index database
ri_db = load_refractive_index_database("data/refractive_indices_database.yaml")

# Compute optical properties at multiple wavelengths
wavelengths = [0.4, 0.55, 0.86, 1.6]  # μm
opt_props = compute_optical_properties(data, wavelengths, ri_db)

# Access results
extinction = opt_props["extinction"]  # [n_levels × n_wavelengths], km⁻¹
ssa = opt_props["ssa"]                # Single scattering albedo
g = opt_props["asymmetry_parameter"]  # Asymmetry parameter
```

### 4. Working with Refractive Indices

```julia
# Load database
db = load_refractive_index_database("data/refractive_indices_database.yaml")

# Get refractive index at specific wavelength
n_sulfate = get_refractive_index(db, "sulfate_suso", 0.55)  # At 550 nm
println("n = $(real(n_sulfate)) + $(imag(n_sulfate))i")

# Get at multiple wavelengths
wavelengths = [0.4, 0.55, 0.86]
n_vec = get_refractive_index(db, "organic_carbon", wavelengths)

# List available species
species = list_species(db)  # ["sulfate_suso", "organic_carbon", ...]

# Show database info
show_database_info(db)
```

## Configuration Files

### TOMAS-15 Configuration (`aerosol_config_tomas15.yaml`)

```yaml
aerosol_scheme:
  type: "TOMAS15"
  species:
    DUST:
      description: "Mineral dust"
      refractive_index: "dust_opac"
      density: 2600.0  # kg/m³
      molar_mass: 0.1  # kg/mol
    SF:
      description: "Sulfate (fine mode)"
      refractive_index: "sulfate_suso"
      density: 1770.0
      molar_mass: 0.098
    # ... more species
  
  size_bins:
    n_bins: 15
    diam_min_nm: 10.0
    diam_max_nm: 10000.0
    spacing: "logarithmic"

netcdf_mapping:
  concentration_pattern: "SpeciesConcVV_{species}{bin:02d}"
  meteorology:
    pressure: "Met_PMID"
    temperature: "Met_T"

processing_options:
  vertical_flip: true  # Convert BOA→TOA to TOA→BOA
```

### Two-Moment Configuration (`aerosol_config_two_moment.yaml`)

```yaml
aerosol_scheme:
  type: "TwoMoment"
  species:
    so4:
      description: "Sulfate"
      refractive_index: "sulfate_suso"
      sigma_g: 1.6
      aod_reference_wavelength: 0.55  # μm
      aod_variable: "AODHyg550nm_{species}"
      radius_variable: "Chem_AeroRadi{species}"
    # ... more species

processing_options:
  vertical_flip: false
```

## Data Structure

### Type Hierarchy

```
AerosolScheme (abstract)
├── TOMAS15Scheme
│   ├── species::Vector{String}
│   ├── n_bins::Int (15)
│   ├── bin_edges::Vector{Float64}
│   ├── bin_centers::Vector{Float64}
│   └── Physical properties (RI keys, densities, molar masses)
└── TwoMomentScheme
    ├── species::Vector{String}
    ├── sigma_g::Dict{String, Float64}
    ├── aod_wavelength::Dict{String, Float64}
    └── refractive_indices::Dict{String, String}

AerosolData{T<:AerosolScheme}
├── scheme::T
├── species_data::Dict{String, AerosolSpeciesData}
├── coordinates::Dict{String, Array}
└── metadata::Dict{String, Any}

AerosolSpeciesData
├── data::Dict{String, Array}
├── units::Dict{String, String}
└── description::String
```

### TOMAS-15 Data Format

Each species has a concentration array: `[n_bins=15 × n_levels]`

Size bins (dry diameter):
- Bin 1: 10.0 - 13.4 nm
- Bin 2: 13.4 - 17.9 nm
- ...
- Bin 15: 7464.5 - 10000.0 nm

### Two-Moment Data Format

Each species has:
- `aod`: Vector{Float64} [n_levels] - AOD at reference wavelength
- `radius`: Vector{Float64} [n_levels] - Effective radius (μm)
- Fixed σ_g (geometric standard deviation)

## Refractive Index Database

The database (`data/refractive_indices_database.yaml`) contains wavelength-dependent complex refractive indices for:

| Species | Description | Source | Wavelength Range |
|---------|-------------|--------|------------------|
| `sulfate_suso` | Sulfate (H₂SO₄/H₂O 75%) | OPAC | 0.3 - 3.75 μm |
| `organic_carbon` | Water-soluble organic carbon | OPAC | 0.3 - 4.0 μm |
| `black_carbon` | Black carbon (soot) | Bond & Bergstrom | 0.3 - 2.5 μm |
| `seasalt_sscm` | Sea salt (NaCl) | OPAC | 0.3 - 4.0 μm |
| `dust_opac` | Mineral dust | OPAC | 0.3 - 4.0 μm |
| `water` | Liquid water | refractiveindex.info | 0.3 - 4.0 μm |

Linear interpolation is used between wavelength points.

## Optical Property Calculations

### TOMAS-15 Method

For each size bin and species:

1. Get refractive index n(λ) from database
2. Compute size parameter x = 2πr/λ
3. Calculate Mie efficiencies Q_ext, Q_sca, Q_abs, g
4. Convert concentrations to number densities
5. Compute optical cross-sections
6. Sum over all bins and species

**Output**: Extinction, scattering, absorption coefficients (km⁻¹), SSA, g

### Two-Moment Method

For each species:

1. Scale AOD from reference to target wavelengths (Ångström law)
2. Get refractive index n(λ) from database
3. Compute SSA from refractive index
4. Partition extinction into scattering and absorption

**Output**: Same format as TOMAS-15

## Helper Functions

### Unit Conversions

```julia
# VMR to number concentration
n_conc = compute_number_concentration(vmr, pressure, temperature)  # #/cm³

# VMR to mass concentration
mass_conc = compute_mass_concentration(vmr, molar_mass, pressure, temperature)  # μg/m³
```

### Size Distributions

```julia
# Lognormal distribution dN/dr
dN_dr = lognormal_size_distribution(r, r_eff, σ_g)

# Effective radius from median radius
r_eff = effective_radius_from_moments(r_med, σ_g)

# Median radius from effective radius
r_med = median_radius_from_effective(r_eff, σ_g)
```

### Wavelength Scaling

```julia
# Scale AOD using Ångström exponent
aod_new = scale_aod_wavelength(aod_ref, λ_ref, λ_target, angstrom_exponent)
```

## Integration with vSmartMOM RT

The aerosol module is designed to integrate seamlessly with vSmartMOM's radiative transfer core:

```julia
# 1. Read aerosol data
aerosol_data = read_aerosol_data(config_file, netcdf_file)

# 2. Compute optical properties
ri_db = load_refractive_index_database(ri_database_file)
opt_props = compute_optical_properties(aerosol_data, wavelengths, ri_db)

# 3. Use in RT calculation
# (Integration with existing vSmartMOM RT setup)
# - Add aerosol extinction to gas absorption
# - Include aerosol scattering in phase function
# - Use aerosol SSA for radiative balance
```

## Testing

Run the test suite:

```julia
using Pkg
Pkg.test("vSmartMOM")
```

Or specifically for aerosols:

```julia
include("test/test_Aerosols.jl")
```

Tests cover:
- Refractive index database loading and interpolation
- TOMAS-15 and two-moment scheme construction
- Data reading from NetCDF files
- Helper function calculations
- Mie scattering approximations
- Optical property computations
- Physical constraint validation

## Future Extensions

The framework is designed to easily accommodate:

1. **MAAM (Modal Aerosol Model)**: Add `MAAMScheme` type
2. **CARMA**: Community Aerosol and Radiation Model for Atmospheres
3. **Custom schemes**: Define your own `AerosolScheme` subtype
4. **Full Mie implementation**: Replace placeholder with exact Mie calculation
5. **Phase matrix calculations**: Full scattering matrix for vector RT
6. **Aerosol mixing**: Internal/external mixtures, core-shell particles
7. **Hygroscopic growth**: Humidity-dependent size and refractive index

## References

- **TOMAS**: Adams, P. J., & Seinfeld, J. H. (2002). *Predicting global aerosol size distributions in general circulation models*. Journal of Geophysical Research, 107(D19).

- **OPAC**: Hess, M., Koepke, P., & Schult, I. (1998). *Optical properties of aerosols and clouds: The software package OPAC*. Bulletin of the American Meteorological Society, 79(5), 831-844.

- **Black Carbon**: Bond, T. C., & Bergstrom, R. W. (2006). *Light absorption by carbonaceous particles: An investigative review*. Aerosol Science and Technology, 40(1), 27-67.

- **GEOSChem**: Bey, I., et al. (2001). *Global modeling of tropospheric chemistry with assimilated meteorology: Model description and evaluation*. Journal of Geophysical Research, 106(D19), 23073-23095.

## Contact

For questions or issues, please open an issue on the vSmartMOM GitHub repository.
