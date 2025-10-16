"""
    Refractive index database loading and interpolation
"""

"""
    load_refractive_index_database(yaml_file::String, FT=Float64)

Load refractive index database from YAML file.

# Arguments
- `yaml_file::String`: Path to YAML database file
- `FT`: Floating point type (default: Float64)

# Returns
- `RefractiveIndexDatabase{FT}`: Database with refractive indices for all species
"""
function load_refractive_index_database(yaml_file::String, FT=Float64)
    # Read YAML file
    yaml_data = YAML.load_file(yaml_file)
    ri_data = yaml_data["refractive_indices"]
    
    # Parse each species
    luts = Dict{String, RefractiveIndexLUT{FT}}()
    
    for (species_key, species_data) in ri_data
        wavelengths = FT.(species_data["wavelengths"])
        n_real = FT.(species_data["n_real"])
        n_imag = FT.(species_data["n_imag"])
        source = get(species_data, "source", "Unknown")
        description = get(species_data, "description", "")
        
        luts[species_key] = RefractiveIndexLUT{FT}(
            species_key, wavelengths, n_real, n_imag, source, description
        )
    end
    
    return RefractiveIndexDatabase{FT}(luts)
end

"""
    get_refractive_index(db::RefractiveIndexDatabase{FT}, species::String, λ::Real) where FT

Get complex refractive index at specified wavelength via linear interpolation.

# Arguments
- `db::RefractiveIndexDatabase{FT}`: Database with refractive indices
- `species::String`: Species key (e.g., "sulfate_suso", "organic_carbon")
- `λ::Real`: Wavelength in micrometers (μm)

# Returns
- `Complex{FT}`: Complex refractive index n = n_real + i*n_imag

# Throws
- `KeyError`: If species not found in database
- `ArgumentError`: If wavelength outside database range

# Example
```julia
db = load_refractive_index_database("data/refractive_indices_database.yaml")
n = get_refractive_index(db, "sulfate_suso", 0.55)  # At 550 nm
```
"""
function get_refractive_index(db::RefractiveIndexDatabase{FT}, species::String, λ::Real) where FT
    # Check if species exists
    if !haskey(db.data, species)
        error("Species '$species' not found in refractive index database. " *
              "Available species: $(keys(db.data))")
    end
    
    lut = db.data[species]
    
    # Check wavelength range
    λ_min, λ_max = extrema(lut.wavelengths)
    if λ < λ_min || λ > λ_max
        error("Wavelength $λ μm outside database range [$λ_min, $λ_max] μm for species '$species'")
    end
    
    # Linear interpolation for real and imaginary parts
    itp_real = LinearInterpolation(lut.wavelengths, lut.n_real)
    itp_imag = LinearInterpolation(lut.wavelengths, lut.n_imag)
    
    n_real_interp = itp_real(λ)
    n_imag_interp = itp_imag(λ)
    
    return Complex{FT}(n_real_interp, n_imag_interp)
end

"""
    get_refractive_index(db::RefractiveIndexDatabase{FT}, species::String, λ::AbstractVector) where FT

Get complex refractive indices at multiple wavelengths.

# Arguments
- `db::RefractiveIndexDatabase{FT}`: Database with refractive indices
- `species::String`: Species key
- `λ::AbstractVector`: Wavelengths in micrometers (μm)

# Returns
- `Vector{Complex{FT}}`: Complex refractive indices at each wavelength
"""
function get_refractive_index(db::RefractiveIndexDatabase{FT}, species::String, λ::AbstractVector) where FT
    return [get_refractive_index(db, species, λi) for λi in λ]
end

"""
    list_species(db::RefractiveIndexDatabase{FT}) where FT

List all available species in the refractive index database.

# Returns
- `Vector{String}`: Species keys
"""
function list_species(db::RefractiveIndexDatabase{FT}) where FT
    return collect(keys(db.data))
end

"""
    wavelength_range(db::RefractiveIndexDatabase{FT}, species::String) where FT

Get the wavelength range for a given species.

# Returns
- `Tuple{FT, FT}`: (λ_min, λ_max) in micrometers
"""
function wavelength_range(db::RefractiveIndexDatabase{FT}, species::String) where FT
    if !haskey(db.data, species)
        error("Species '$species' not found in database")
    end
    
    lut = db.data[species]
    return (minimum(lut.wavelengths), maximum(lut.wavelengths))
end

"""
    show_database_info(db::RefractiveIndexDatabase{FT}) where FT

Print summary information about the refractive index database.
"""
function show_database_info(db::RefractiveIndexDatabase{FT}) where FT
    println("Refractive Index Database")
    println("=" ^ 60)
    println("Number of species: ", length(db.data))
    println()
    
    for (species_key, lut) in sort(collect(db.data))
        λ_min, λ_max = extrema(lut.wavelengths)
        n_points = length(lut.wavelengths)
        
        println("Species: $species_key")
        println("  Description: $(lut.description)")
        println("  Source: $(lut.source)")
        println("  Wavelength range: $λ_min - $λ_max μm")
        println("  Number of data points: $n_points")
        println("  Sample n at $(lut.wavelengths[1]) μm: $(lut.n_real[1]) + $(lut.n_imag[1])i")
        println()
    end
end
