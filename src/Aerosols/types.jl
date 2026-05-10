"""
    Type definitions for flexible aerosol framework
"""

# ============================================================================
# Abstract base types
# ============================================================================

"""
    AerosolScheme

Abstract base type for different aerosol schemes.
Concrete implementations: TOMAS15Scheme, TwoMomentScheme, etc.
"""
abstract type AerosolScheme end

# ============================================================================
# Concrete scheme types
# ============================================================================

"""
    TOMAS15Scheme{FT} <: AerosolScheme

TOMAS microphysics with 15 size bins.
Size-resolved aerosol concentrations in logarithmically-spaced diameter bins.

# Fields
- `species::Vector{String}`: Aerosol species names (e.g., ["DUST", "SS", "SF", ...])
- `n_bins::Int`: Number of size bins (15)
- `diam_min::FT`: Minimum dry diameter (nm)
- `diam_max::FT`: Maximum dry diameter (nm)
- `bin_edges::Vector{FT}`: Bin edge diameters (nm), length n_bins+1
- `bin_centers::Vector{FT}`: Bin center diameters (nm), length n_bins
- `refractive_indices::Dict{String, String}`: Species → RI database key
- `densities::Dict{String, FT}`: Species → density (kg/m³)
- `molar_masses::Dict{String, FT}`: Species → molar mass (kg/mol)
"""
struct TOMAS15Scheme{FT} <: AerosolScheme
    species::Vector{String}
    n_bins::Int
    diam_min::FT
    diam_max::FT
    bin_edges::Vector{FT}
    bin_centers::Vector{FT}
    refractive_indices::Dict{String, String}
    densities::Dict{String, FT}
    molar_masses::Dict{String, FT}
end

"""
    TwoMomentScheme{FT} <: AerosolScheme

Two-moment aerosol scheme with lognormal distributions.
Each species characterized by AOD, effective radius, and fixed σ_g.

# Fields
- `species::Vector{String}`: Aerosol species names (e.g., ["so4", "ocpi", ...])
- `sigma_g::Dict{String, FT}`: Geometric standard deviation per species
- `aod_wavelength::Dict{String, FT}`: Reference wavelength for AOD (μm)
- `refractive_indices::Dict{String, String}`: Species → RI database key
"""
struct TwoMomentScheme{FT} <: AerosolScheme
    species::Vector{String}
    sigma_g::Dict{String, FT}
    aod_wavelength::Dict{String, FT}
    refractive_indices::Dict{String, String}
end

# ============================================================================
# Data containers
# ============================================================================

"""
    AerosolSpeciesData

Data for a single aerosol species.
Structure varies by scheme type.

# Fields
- `data::Dict{String, Array}`: Variable name → array data
  - TOMAS15: "concentration" → Array{Float64, 2} (n_bins × n_levels)
  - TOMAS15 NK: may also include nested species-fraction dictionaries
  - TwoMoment: "aod" → Vector{Float64}, "radius" → Vector{Float64}
- `units::Dict{String, String}`: Variable name → units
- `description::String`: Human-readable description
"""
struct AerosolSpeciesData
    data::Dict{String, Any}
    units::Dict{String, String}
    description::String
end

"""
    AerosolData{T<:AerosolScheme}

Container for all aerosol data from a given scheme.

# Fields
- `scheme::T`: The aerosol scheme specification
- `species_data::Dict{String, AerosolSpeciesData}`: Species name → data
- `coordinates::Dict{String, Array}`: Coordinate arrays (lon, lat, lev, time)
- `metadata::Dict{String, Any}`: Additional metadata from NetCDF
"""
struct AerosolData{T<:AerosolScheme}
    scheme::T
    species_data::Dict{String, AerosolSpeciesData}
    coordinates::Dict{String, Array}
    metadata::Dict{String, Any}
end

# ============================================================================
# Refractive index types
# ============================================================================

"""
    RefractiveIndexLUT{FT}

Lookup table for wavelength-dependent refractive index.

# Fields
- `species::String`: Species identifier
- `wavelengths::Vector{FT}`: Wavelength values (μm)
- `n_real::Vector{FT}`: Real part of refractive index
- `n_imag::Vector{FT}`: Imaginary part of refractive index
- `source::String`: Data source/reference
- `description::String`: Human-readable description
"""
struct RefractiveIndexLUT{FT}
    species::String
    wavelengths::Vector{FT}
    n_real::Vector{FT}
    n_imag::Vector{FT}
    source::String
    description::String
end

"""
    RefractiveIndexDatabase{FT}

Database of refractive indices for multiple aerosol species.

# Fields
- `data::Dict{String, RefractiveIndexLUT{FT}}`: Species key → LUT
"""
struct RefractiveIndexDatabase{FT}
    data::Dict{String, RefractiveIndexLUT{FT}}
end

# ============================================================================
# Helper constructors
# ============================================================================

"""
    TOMAS15Scheme(config::Dict, FT=Float64)

Construct TOMAS15Scheme from YAML configuration dictionary.
"""
function TOMAS15Scheme(config::Dict, FT=Float64)
    species_config = config["aerosol_scheme"]["species"]
    size_config = config["aerosol_scheme"]["size_bins"]
    
    # Extract species names
    species = collect(keys(species_config))
    
    # Size bin configuration
    n_bins = size_config["n_bins"]
    diam_min = FT(size_config["diam_min_nm"])
    diam_max = FT(size_config["diam_max_nm"])
    
    # Calculate logarithmic bin edges
    bin_edges = diam_min .* (diam_max / diam_min) .^ (FT.(collect(0:n_bins)) ./ FT(n_bins))
    bin_centers = sqrt.(bin_edges[1:end-1] .* bin_edges[2:end])
    
    # Extract species properties
    refractive_indices = Dict{String, String}()
    densities = Dict{String, FT}()
    molar_masses = Dict{String, FT}()
    
    for (sp, sp_config) in species_config
        refractive_indices[sp] = sp_config["refractive_index"]
        densities[sp] = FT(sp_config["density"])
        molar_masses[sp] = FT(sp_config["molar_mass"])
    end
    
    return TOMAS15Scheme{FT}(
        species, n_bins, diam_min, diam_max,
        bin_edges, bin_centers,
        refractive_indices, densities, molar_masses
    )
end

"""
    TwoMomentScheme(config::Dict, FT=Float64)

Construct TwoMomentScheme from YAML configuration dictionary.
"""
function TwoMomentScheme(config::Dict, FT=Float64)
    species_config = config["aerosol_scheme"]["species"]
    
    species = collect(keys(species_config))
    sigma_g = Dict{String, FT}()
    aod_wavelength = Dict{String, FT}()
    refractive_indices = Dict{String, String}()
    
    for (sp, sp_config) in species_config
        sigma_g[sp] = FT(sp_config["sigma_g"])
        aod_wavelength[sp] = FT(sp_config["aod_reference_wavelength"])
        refractive_indices[sp] = sp_config["refractive_index"]
    end
    
    return TwoMomentScheme{FT}(species, sigma_g, aod_wavelength, refractive_indices)
end
