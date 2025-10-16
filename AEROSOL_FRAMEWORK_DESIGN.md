# Flexible Aerosol Framework Design

## Overview

Design for a flexible aerosol system that supports multiple aerosol schemes (TOMAS-15, two-moment, future extensions) with wavelength-dependent refractive indices.

## Architecture

### 1. Type Hierarchy

```julia
# Abstract base type
abstract type AerosolScheme end

# Concrete implementations
struct TOMAS15Scheme <: AerosolScheme
    species::Vector{String}          # e.g., ["DUST", "SS", "SF", ...]
    n_bins::Int                       # 15
    diam_min::Float64                 # 10 nm
    diam_max::Float64                 # 10 μm (10000 nm)
    bin_edges::Vector{Float64}        # Diameter edges (nm)
    bin_centers::Vector{Float64}      # Diameter centers (nm)
end

struct TwoMomentScheme <: AerosolScheme
    species::Vector{String}           # e.g., ["so4", "ocpi", "sala", ...]
    sigma_g::Dict{String, Float64}    # Fixed geometric std dev per species
end

# Future extensibility
struct MAAMScheme <: AerosolScheme
    # Modal Aerosol Module (could add later)
end

struct CustomBinnedScheme <: AerosolScheme
    # Custom number of bins and ranges
    n_bins::Int
    bin_edges::Vector{Float64}
end
```

### 2. Aerosol Data Container

```julia
struct AerosolData{T<:AerosolScheme}
    scheme::T
    species_data::Dict{String, AerosolSpeciesData}
    metadata::Dict{String, Any}
end

# Species-specific data (flexible for different schemes)
struct AerosolSpeciesData
    # For TOMAS: concentration[bin, level]
    # For two-moment: aod[level], radius[level]
    data::Dict{String, Array}
    units::Dict{String, String}
end
```

### 3. Refractive Index System

```julia
# Refractive index LUT
struct RefractiveIndexLUT
    species::String
    wavelengths::Vector{Float64}      # μm
    n_real::Vector{Float64}
    n_imag::Vector{Float64}
    source::String                    # Reference/citation
end

# Container for all species
struct RefractiveIndexDatabase
    data::Dict{String, RefractiveIndexLUT}
end

# Interpolation function
function get_refractive_index(db::RefractiveIndexDatabase, 
                              species::String, 
                              λ::Float64)
    # Linear interpolation in wavelength
end
```

### 4. YAML Configuration Structure

#### For TOMAS-15:
```yaml
aerosol_scheme:
  type: "TOMAS15"
  file: "/path/to/GEOSChem.Custom.20190702_0000z.nc4"
  
  species:
    DUST:
      nc_prefix: "SpeciesConcVV_DUST"
      description: "Mineral dust"
      refractive_index: "dust_opac"  # Reference to LUT
    SS:
      nc_prefix: "SpeciesConcVV_SS"
      description: "Sea salt"
      refractive_index: "seasalt_sscm"
    SF:
      nc_prefix: "SpeciesConcVV_SF"
      description: "Sulfate"
      refractive_index: "sulfate_suso"
    # ... etc
  
  size_bins:
    n_bins: 15
    diam_min_nm: 10.0
    diam_max_nm: 10000.0
    type: "logarithmic"
  
  dimensions:
    lon: "Xdim"
    lat: "Ydim"
    lev: "lev"
    face: "nf"
    time: "time"

  meteorology:
    temperature: "Met_T"
    pressure_dry: "Met_DELPDRY"
    pressure_wet: "Met_DELP"
    surface_pressure: "Met_PS2WET"
```

#### For Two-Moment:
```yaml
aerosol_scheme:
  type: "TwoMoment"
  file: "/path/to/GEOSChem.Aerosols.20190601_0000z.nc4"
  
  species:
    so4:
      aod_var: "AODHyg550nm_SO4"
      radius_var: "Chem_AeroRadiSULF"
      sigma_g: 1.6
      description: "Sulfate aerosol"
      refractive_index: "sulfate_suso"
      aod_wavelength: 0.550  # μm (reference wavelength)
    ocpi:
      aod_var: "AODHyg550nm_OCPI"
      radius_var: "Chem_AeroRadiOC"
      sigma_g: 1.6
      description: "Organic carbon (hydrophilic)"
      refractive_index: "organic_carbon"
      aod_wavelength: 0.550
    # ... etc
  
  dimensions:
    lon: "lon"
    lat: "lat"
    lev: "lev"
    time: "time"
```

#### Refractive Index Database:
```yaml
refractive_indices:
  sulfate_suso:
    source: "GEOS-ESM/GEOSmie suso00"
    citation: "OPAC (Hess et al. 1998)"
    wavelengths: [0.300, 0.400, 0.550, 0.670, 1.000, 1.500, 2.000, 2.500]  # μm
    n_real:      [1.430, 1.430, 1.432, 1.432, 1.435, 1.440, 1.450, 1.870]
    n_imag:      [1e-8,  1e-8,  1e-8,  1e-8,  5e-3,  1e-2,  2e-2,  3.15e-2]
  
  organic_carbon:
    source: "GEOS-ESM/GEOSmie data"
    wavelengths: [0.440, 0.550, 0.670, 1.000, 1.600, 2.250]
    n_real:      [1.520, 1.520, 1.520, 1.510, 1.450, 1.384]
    n_imag:      [8e-3,  8e-3,  8e-3,  8e-3,  1e-2,  1.26e-3]
  
  dust_opac:
    source: "OPAC dust-like"
    wavelengths: [0.400, 0.550, 0.670, 1.000, 1.500, 2.250, 3.000]
    n_real:      [1.530, 1.530, 1.530, 1.530, 1.520, 1.500, 1.480]
    n_imag:      [8e-3,  8e-3,  8e-3,  8e-3,  1e-2,  1.5e-2, 2e-2]
  
  seasalt_sscm:
    source: "GEOS-ESM/GEOSmie sscm00"
    wavelengths: [0.400, 0.550, 0.670, 1.000, 1.500, 2.250]
    n_real:      [1.500, 1.500, 1.498, 1.490, 1.450, 1.430]
    n_imag:      [1e-8,  1e-8,  1e-8,  1e-3,  2e-3,  4e-3]
  
  black_carbon:
    source: "GEOS-ESM/GEOSmie soot00"
    wavelengths: [0.440, 0.550, 0.670, 1.000, 1.600, 2.250]
    n_real:      [1.750, 1.750, 1.780, 1.800, 1.810, 1.820]
    n_imag:      [0.440, 0.440, 0.460, 0.480, 0.500, 0.510]
```

### 5. Reader Interface

```julia
# Generic reader function
function read_aerosol_data(config_file::String)
    config = YAML.load_file(config_file)
    scheme_type = config["aerosol_scheme"]["type"]
    
    if scheme_type == "TOMAS15"
        return read_tomas15(config)
    elseif scheme_type == "TwoMoment"
        return read_two_moment(config)
    else
        error("Unknown aerosol scheme: $scheme_type")
    end
end

# TOMAS15 reader
function read_tomas15(config::Dict)
    # Read NetCDF, extract binned concentrations
    # Return AerosolData{TOMAS15Scheme}
end

# Two-moment reader
function read_two_moment(config::Dict)
    # Read NetCDF, extract AOD and radius
    # Return AerosolData{TwoMomentScheme}
end
```

### 6. Optical Property Calculator

```julia
# Generic interface
function compute_optical_properties(aerosol::AerosolData, 
                                    λ::Float64,
                                    ri_db::RefractiveIndexDatabase)
    # Dispatch based on scheme type
end

# TOMAS15 implementation
function compute_optical_properties(aerosol::AerosolData{TOMAS15Scheme},
                                    λ::Float64,
                                    ri_db::RefractiveIndexDatabase)
    # For each species and bin:
    #   1. Get refractive index at λ
    #   2. Compute Mie scattering for bin size
    #   3. Multiply by concentration
    #   4. Sum over bins
end

# Two-moment implementation
function compute_optical_properties(aerosol::AerosolData{TwoMomentScheme},
                                    λ::Float64,
                                    ri_db::RefractiveIndexDatabase)
    # For each species:
    #   1. Get refractive index at λ
    #   2. Scale AOD from reference wavelength to λ
    #   3. Use effective radius and σ_g for phase function
end
```

## Implementation Strategy

### Phase 1: Core Infrastructure
1. ✅ Create type definitions
2. ✅ Implement YAML parser for both schemes
3. ✅ Create refractive index database and LUT reader
4. ✅ Implement interpolation functions

### Phase 2: TOMAS15 Support
1. ✅ Create TOMAS15 reader (based on Python exploration)
2. ✅ Implement concentration → optical properties conversion
3. ✅ Add tests with GEOSChem.Custom file

### Phase 3: Two-Moment Support
1. ✅ Create two-moment reader (based on existing code)
2. ✅ Implement AOD scaling and optical properties
3. ✅ Add tests with GEOSChem.Aerosols file

### Phase 4: Integration
1. ✅ Create unified interface for vSmartMOM RT
2. ✅ Add wavelength-dependent calculations
3. ✅ Documentation and examples

## File Structure

```
src/
  Aerosols/
    Aerosols.jl                      # Main module
    types.jl                         # Type definitions
    schemes/
      tomas15.jl                     # TOMAS15 implementation
      two_moment.jl                  # Two-moment implementation
    refractive_index.jl              # RI database and interpolation
    readers.jl                       # NetCDF readers
    optical_properties.jl            # Mie/optical calculations
    
data/
  refractive_indices/
    default_ri_database.yaml         # Default RI values
    
examples/
  aerosol_tomas15_example.jl
  aerosol_two_moment_example.jl
  
test/
  test_aerosols.jl
  aerosol_config_tomas15.yaml        # Test config
  aerosol_config_two_moment.yaml     # Test config
```

## Benefits of This Design

1. **Extensible**: Easy to add new schemes (MAAM, custom bins, etc.)
2. **Flexible**: YAML configuration allows different NetCDF structures
3. **Wavelength-dependent**: Full spectral calculations supported
4. **Type-safe**: Julia dispatch ensures correct methods
5. **Testable**: Each component can be tested independently
6. **Documented**: Clear separation of concerns

## Next Steps

1. Implement core types and YAML reader
2. Create refractive index database
3. Port TOMAS15 exploration code to Julia
4. Integrate with existing two-moment code
5. Add to vSmartMOM RT workflow
