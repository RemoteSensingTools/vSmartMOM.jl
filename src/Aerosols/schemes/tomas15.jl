"""
    TOMAS-15 aerosol scheme reader
    
    Reads size-resolved aerosol concentrations from GEOSChem NetCDF files
    with 15 logarithmically-spaced size bins (10 nm - 10 μm dry diameter).
"""

"""
    read_tomas15(config::Dict, netcdf_file::String, FT=Float64)

Read TOMAS-15 aerosol data from NetCDF file.

# Arguments
- `config::Dict`: YAML configuration dictionary
- `netcdf_file::String`: Path to NetCDF file
- `FT`: Floating point type (default: Float64)

# Returns
- `AerosolData{TOMAS15Scheme{FT}}`: Container with size-resolved aerosol data

# Data Structure
Each species has concentrations in mol/mol dry air for 15 size bins.
NetCDF variables: SpeciesConcVV_[SPECIES][BIN] (e.g., SpeciesConcVV_DUST01)

# Configuration
The config dict must specify:
- aerosol_scheme/species: Dict of species with properties
- aerosol_scheme/size_bins: Size bin configuration
- netcdf_mapping: Variable name patterns
- processing_options/vertical_flip: Whether to flip BOA→TOA to TOA→BOA
"""
function read_tomas15(config::Dict, netcdf_file::String, FT=Float64)
    # Construct scheme from config
    scheme = TOMAS15Scheme(config, FT)
    
    # Open NetCDF file
    ds = NCDataset(netcdf_file)
    
    try
        # Extract coordinates
        coordinates = extract_coordinates(ds, config)
        n_levels = length(coordinates["lev"])
        
        # Get processing options
        vertical_flip = get(config["processing_options"], "vertical_flip", false)
        
        # Get meteorology for number density conversion (if needed)
        met_vars = config["netcdf_mapping"]["meteorology"]
        
        # Read each species
        species_data = Dict{String, AerosolSpeciesData}()
        
        for species_name in scheme.species
            # NetCDF variable pattern
            var_pattern = replace(
                config["netcdf_mapping"]["concentration_pattern"],
                "{species}" => species_name,
                "{bin:02d}" => ""  # Will be replaced in loop
            )
            
            # Read all bins for this species
            concentrations = zeros(FT, scheme.n_bins, n_levels)
            
            for ibin in 1:scheme.n_bins
                var_name = replace(var_pattern, "{bin:02d}" => @sprintf("%02d", ibin))
                
                if haskey(ds, var_name)
                    # Read data (assumes dimensions: nf, Xdim, Ydim, lev, time)
                    # Extract first time step and aggregate over cubed-sphere faces/horizontal
                    data_full = Array(ds[var_name])
                    
                    # Average over horizontal dimensions (faces, X, Y)
                    # Shape: (nf=6, Xdim=24, Ydim=24, lev=72, time=1)
                    data_profile = dropdims(
                        mean(data_full, dims=(1, 2, 3)),
                        dims=(1, 2, 3)
                    )  # Now shape: (lev, time)
                    
                    # Extract first time step
                    concentrations[ibin, :] = data_profile[:, 1]
                else
                    @warn "Variable $var_name not found for species $species_name, bin $ibin"
                    concentrations[ibin, :] .= 0.0
                end
            end
            
            # Apply vertical flip if requested
            if vertical_flip
                concentrations = reverse(concentrations, dims=2)
            end
            
            # Store species data
            data_dict = Dict{String, Array}(
                "concentration" => concentrations
            )
            
            units_dict = Dict{String, String}(
                "concentration" => "mol mol-1 dry air"
            )
            
            species_data[species_name] = AerosolSpeciesData(
                data_dict,
                units_dict,
                "TOMAS-15 size-resolved $(species_name) aerosol"
            )
        end
        
        # Extract metadata
        metadata = extract_metadata(ds)
        
        return AerosolData(scheme, species_data, coordinates, metadata)
        
    finally
        close(ds)
    end
end

"""
    compute_number_concentration(vmr::AbstractVector, pressure::AbstractVector, 
                                  temperature::AbstractVector)

Convert volume mixing ratio to number concentration.

# Arguments
- `vmr::AbstractVector`: Volume mixing ratio (mol/mol)
- `pressure::AbstractVector`: Pressure (Pa)
- `temperature::AbstractVector`: Temperature (K)

# Returns
- `Vector{Float64}`: Number concentration (#/cm³)

# Formula
n = VMR × (P / (k_B × T)) × 10^-6

where k_B = 1.380649e-23 J/K (Boltzmann constant)
"""
function compute_number_concentration(
    vmr::AbstractVector,
    pressure::AbstractVector,
    temperature::AbstractVector
)
    k_B = 1.380649e-23  # J/K
    
    # Number density of air (molecules/m³)
    n_air = pressure ./ (k_B .* temperature)
    
    # Number concentration (#/m³)
    n_aerosol = vmr .* n_air
    
    # Convert to #/cm³
    return n_aerosol .* 1e-6
end

"""
    compute_mass_concentration(vmr::AbstractVector, molar_mass::Float64,
                               pressure::AbstractVector, temperature::AbstractVector)

Convert volume mixing ratio to mass concentration.

# Arguments
- `vmr::AbstractVector`: Volume mixing ratio (mol/mol)
- `molar_mass::Float64`: Molar mass (kg/mol)
- `pressure::AbstractVector`: Pressure (Pa)
- `temperature::AbstractVector`: Temperature (K)

# Returns
- `Vector{Float64}`: Mass concentration (μg/m³)

# Formula
ρ = VMR × (P × M) / (R × T) × 10^9

where R = 8.314462618 J/(mol·K) (universal gas constant)
"""
function compute_mass_concentration(
    vmr::AbstractVector,
    molar_mass::Float64,
    pressure::AbstractVector,
    temperature::AbstractVector
)
    R = 8.314462618  # J/(mol·K)
    
    # Mass concentration (kg/m³)
    mass_conc = (vmr .* pressure .* molar_mass) ./ (R .* temperature)
    
    # Convert to μg/m³
    return mass_conc .* 1e9
end

"""
    bin_volume(diam_nm::Float64)

Calculate volume of a spherical particle from its diameter.

# Arguments
- `diam_nm::Float64`: Particle diameter (nm)

# Returns
- `Float64`: Particle volume (nm³)
"""
function bin_volume(diam_nm::Float64)
    radius_nm = diam_nm / 2.0
    return (4.0/3.0) * π * radius_nm^3
end
