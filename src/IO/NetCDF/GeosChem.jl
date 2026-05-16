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
    # Shim: forward to the new open-once scene reader. The dict now reflects
    # bug fixes vs. the historical reader (q in kg/kg, H2O excluded from
    # `absorption.molecules`, O2 omitted when no VMR is provided). The
    # `aerosol_scheme=nothing` opt-out skips TOMAS IO that the legacy dict
    # discarded anyway.
    scene = read_gchp_scene(src.path, src.idx, src.idy, src.idf;
                            aerosol_scheme = nothing)
    @info "Reading GEOSChem scene at lat=$(scene.lat)°, lon=$(scene.lon)°, time=$(scene.time)"
    return scene_to_dict(scene)
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
