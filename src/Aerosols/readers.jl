"""
    Generic aerosol data reader with YAML configuration support
"""

"""
    aerosol_processing_options(config::Dict)

Return aerosol reader options from either the documented top-level
`processing_options` table or the legacy `aerosol_scheme.options` table.
Missing options are represented by an empty dictionary.
"""
function aerosol_processing_options(config::Dict)
    if haskey(config, "processing_options")
        return config["processing_options"]
    end

    aerosol_scheme = get(config, "aerosol_scheme", Dict())
    return get(aerosol_scheme, "options", Dict())
end

"""
    read_aerosol_data(config_file::String, netcdf_file::String, FT=Float64)

Read aerosol data from NetCDF file using YAML configuration.

Dispatches to scheme-specific reader based on configuration.

# Arguments
- `config_file::String`: Path to YAML configuration file
- `netcdf_file::String`: Path to NetCDF data file
- `FT`: Floating point type (default: Float64)

# Returns
- `AerosolData{T}`: Aerosol data container with scheme-specific type

# Example
```julia
# TOMAS scheme
data = read_aerosol_data(
    "examples/aerosol_config_tomas15.yaml",
    "GEOSChem.Custom.20190702_0000z.nc4"
)

# Two-moment scheme
data = read_aerosol_data(
    "examples/aerosol_config_two_moment.yaml",
    "GEOSChem.Aerosols.20190702_0000z.nc4"
)
```
"""
function read_aerosol_data(config_file::String, netcdf_file::String, FT=Float64)
    # Load configuration
    config = YAML.load_file(config_file)
    
    # Determine scheme type
    scheme_name = config["aerosol_scheme"]["type"]
    
    if scheme_name in ("TOMAS", "TOMAS15")
        return read_tomas(TOMASScheme(config, FT), config, netcdf_file)
    elseif scheme_name == "TwoMoment"
        return read_two_moment(config, netcdf_file, FT)
    else
        error("Unknown aerosol scheme: $scheme_name. Supported: TOMAS/TOMAS15, TwoMoment")
    end
end

"""
    extract_coordinates(ds::NCDataset, config::Dict)

Extract coordinate arrays from NetCDF dataset.

# Returns
- `Dict{String, Array}`: Coordinate name → array
"""
function extract_coordinates(ds::NCDataset, config::Dict)
    coords = Dict{String, Array}()
    
    coord_names = ["lon", "lat", "lev", "time"]
    
    for coord in coord_names
        if haskey(ds, coord)
            coords[coord] = Array(ds[coord])
        end
    end
    
    return coords
end

"""
    extract_metadata(ds::NCDataset)

Extract global metadata from NetCDF dataset.

# Returns
- `Dict{String, Any}`: Metadata dictionary
"""
function extract_metadata(ds::NCDataset)
    metadata = Dict{String, Any}()
    
    # Global attributes
    for (key, val) in ds.attrib
        metadata[key] = val
    end
    
    # Dimensions
    metadata["dimensions"] = Dict(name => size for (name, size) in ds.dim)
    
    return metadata
end

"""
    flip_vertical(data::AbstractArray)

Flip vertical dimension from BOA→TOA to TOA→BOA ordering.

# Arguments
- `data::AbstractArray`: Data array with vertical dimension

# Returns
- `AbstractArray`: Flipped array
"""
function flip_vertical(data::AbstractArray)
    # Assumes vertical is last dimension
    ndims_data = ndims(data)
    return reverse(data, dims=ndims_data)
end
