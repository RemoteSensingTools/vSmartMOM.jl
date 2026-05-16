"""
    Type definitions for flexible aerosol framework
"""

# ============================================================================
# Abstract base types
# ============================================================================

"""
    AerosolScheme

Top-level abstract base type for aerosol-scheme metadata objects.
Kept un-parameterised for back-compat — code that wrote `x::AerosolScheme`
or `<: AerosolScheme` continues to work. New schemes should subtype the
parametric [`AbstractAerosolScheme{FT}`](@ref) below, which is itself
`<: AerosolScheme`.
"""
abstract type AerosolScheme end

"""
    AbstractAerosolScheme{FT} <: AerosolScheme

Parametric layer of the [`AerosolScheme`](@ref) hierarchy. Carries the
working float type so downstream sectional payloads (e.g.
[`SectionalAerosolData{FT, S}`](@ref)) can constrain `S` to a matching
`FT`.
"""
abstract type AbstractAerosolScheme{FT} <: AerosolScheme end

# ============================================================================
# Concrete scheme types
# ============================================================================

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
struct TwoMomentScheme{FT} <: AbstractAerosolScheme{FT}
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
