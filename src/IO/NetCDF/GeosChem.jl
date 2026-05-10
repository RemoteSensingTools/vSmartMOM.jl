# GEOSChem-specific NetCDF reading and conversion
# Included directly into IO module (not a separate module)

using .Sources: GeosChemSource
using ..CoreRT: compute_atmos_profile_fields
# load_config is imported at top level of IO module

"""
    _geoschem_error(msg)

Raise a stable `ArgumentError` for invalid GEOS-Chem NetCDF input.
"""
@inline _geoschem_error(msg) = throw(ArgumentError(msg))

"""
    _require_geoschem(cond, msg)

Validate GEOS-Chem NetCDF metadata and raise `ArgumentError` when invalid.
"""
@inline _require_geoschem(cond, msg) = cond ? nothing : _geoschem_error(msg)

"""
    geoschem_to_dict(src::GeosChemSource) -> Dict

Read a GEOSChem NetCDF4 file at the specified grid location and convert
to a Dict compatible with `parameters_from_dict`.

# Arguments
- `src::GeosChemSource`: Source specification with file path and grid indices

# Returns
- `Dict`: Configuration dictionary with atmospheric_profile and absorption sections

# Notes
The GEOSChem file structure has dimensions (Xdim × Ydim × nf × lev × time):
- `nf`: Which of the 6 cubed-sphere faces (1-6)
- `Xdim, Ydim`: Location on that face
- `lev`: Model layers, indexing from BOA (lev=1) to TOA (lev=72)

Data is automatically flipped to match vSmartMOM convention (TOA→BOA indexing).

# Example
```julia
src = GeosChemSource("GEOSChem.Custom.20190101_0000z.nc4", 10, 20, 1)
config_dict = geoschem_to_dict(src)
# Now use with parameters_from_dict(config_dict)
```
"""
function geoschem_to_dict(src::GeosChemSource)
    idx, idy, idf = src.idx, src.idy, src.idf
    
    # Read the NetCDF file
    ds = NCDataset(src.path)
    
    # Extract location information
    lat = ds["lats"][idx, idy, idf]
    lon = ds["lons"][idx, idy, idf]
    time = ds["time"][1]
    
    @info "Reading GEOSChem scene at lat=$(lat)°, lon=$(lon)°, time=$(time)"
    
    # Build pressure grid (73 levels from surface pressure + layer thicknesses)
    # Met_DELP: pressure thickness of each layer [hPa]
    # Met_PS2WET: surface pressure [hPa]
    dp = ds["Met_DELP"].var[idx, idy, idf, :, 1]
    sp = ds["Met_PS2WET"].var[idx, idy, idf, 1]
    
    # Pressure at half-levels (boundaries): [TOA, ..., Surface]
    # GCHP stores bottom-to-top, we flip to top-to-bottom
    pressure_half = reverse([sp; sp .+ cumsum(-dp)])
    
    # Temperature profile [K] - flip from BOA→TOA to TOA→BOA
    _require_geoschem(ds["Met_T"].attrib["units"] == "K", "Expected temperature in Kelvin")
    temperature = reverse(ds["Met_T"].var[idx, idy, idf, :, 1])
    
    # Specific humidity [g/kg] - flip from BOA→TOA to TOA→BOA
    _require_geoschem(ds["Met_SPHU"].attrib["units"] == "g kg-1", "Expected specific humidity in g/kg")
    q = reverse(ds["Met_SPHU"].var[idx, idy, idf, :, 1])
    
    # Read volume mixing ratios for trace gases
    vmr = Dict{String, Vector{Float64}}()
    
    # Define molecules to extract (all available in GEOSChem)
    molecules_to_read = ["N2O", "CH4", "C2H6", "CO2", "CO", "H2O"]
    molecules_available = String[]
    
    for molecule in molecules_to_read
        var_name = "SpeciesConcVV_$(molecule)"
        if haskey(ds, var_name)
            _require_geoschem(ds[var_name].attrib["units"] == "mol mol-1 dry", "Expected VMR in mol/mol")
            # Flip from BOA→TOA to TOA→BOA
            vmr[molecule] = reverse(Float64.(ds[var_name].var[idx, idy, idf, :, 1]))
            push!(molecules_available, molecule)
        end
    end
    
    close(ds)
    
    # Build the configuration dictionary
    config = Dict{String, Any}()
    
    # Atmospheric profile section
    config["atmospheric_profile"] = Dict{String, Any}(
        "T" => Float64.(temperature),
        "p" => Float64.(pressure_half),
        "q" => Float64.(q),
        "vmr" => vmr,
        "profile_reduction" => nothing  # No reduction by default
    )
    
    # Absorption section (if molecules are available)
    if !isempty(molecules_available)
        # Always include O2 (not in GEOSChem output, will use default)
        all_molecules = [["O2", molecules_available...]]
        config["absorption"] = Dict{String, Any}(
            "molecules" => all_molecules,
            "vmr" => vmr,
            "broadening" => "Voigt()",
            "CEF" => "HumlicekWeidemann32SDErrorFunction()",
            "wing_cutoff" => 25
        )
    end
    
    # Add metadata for reference
    config["_metadata"] = Dict{String, Any}(
        "source_file" => src.path,
        "source_type" => "GEOSChem",
        "latitude" => lat,
        "longitude" => lon,
        "time" => time,
        "grid_indices" => (idx, idy, idf)
    )
    
    return config
end

"""
    read_geoschem_profile(file::String, idx::Int, idy::Int, idf::Int) -> Dict

Convenience function to read a GEOSChem file and return the configuration Dict.

# Arguments
- `file::String`: Path to GEOSChem NetCDF4 file
- `idx::Int`: X-dimension index
- `idy::Int`: Y-dimension index  
- `idf::Int`: Face index (1-6)

# Returns
- `Dict`: Configuration dictionary compatible with `parameters_from_dict`

# Example
```julia
config = read_geoschem_profile("GEOSChem.Custom.20190101_0000z.nc4", 10, 20, 1)
params = parameters_from_dict(config)
```
"""
function read_geoschem_profile(file::String, idx::Int, idy::Int, idf::Int)
    src = GeosChemSource(file, idx, idy, idf)
    return geoschem_to_dict(src)
end

# Extend Formats.load_config for GeosChemSource
function load_config(src::GeosChemSource)
    return geoschem_to_dict(src)
end
