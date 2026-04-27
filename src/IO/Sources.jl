# IO Source types for different data formats
# Extends the IOSource pattern from Formats.jl

module Sources

using ..Formats: IOSource

export GeosChemSource, NetCDFGridSource, NetCDFSource

"""
    abstract type NetCDFSource <: IOSource

Abstract type for all NetCDF-based data sources.
Provides a hierarchy for different NetCDF formats (GEOSChem, WRF, GCHP, etc.).
"""
abstract type NetCDFSource <: IOSource end

"""
    _source_error(msg)

Raise a stable `ArgumentError` for invalid IO source configuration.
"""
@inline _source_error(msg) = throw(ArgumentError(msg))

"""
    _require_source(cond, msg)

Validate IO source constructor input and raise `ArgumentError` when invalid.
"""
@inline _require_source(cond, msg) = cond ? nothing : _source_error(msg)

"""
    GeosChemSource(path::String, idx::Int, idy::Int, idf::Int)

IO source for GEOSChem output files (NetCDF4 format).

# Fields
- `path::String`: Path to the GEOSChem NetCDF4 file (e.g., GEOSChem.Custom.YYYYMMDD_HHMMz.nc4)
- `idx::Int`: X-dimension index on the cubed-sphere grid
- `idy::Int`: Y-dimension index on the cubed-sphere grid
- `idf::Int`: Face index (1-6) for the cubed-sphere face

# Usage
```julia
src = GeosChemSource("GEOSChem.Custom.20190101_0000z.nc4", 10, 20, 1)
params = parameters_from_yaml(src)  # Reads NetCDF, converts to vSmartMOM_Parameters
```

# Notes
GEOSChem files have dimensions (Xdim × Ydim × nf × lev × time), where:
- `nf`: Cubed-sphere face (1-6)
- `Xdim, Ydim`: Location on that face
- `lev`: Model layers (1=BOA, 72=TOA for typical setup)
- `time`: Time dimension

The data is automatically flipped from GCHP convention (BOA→TOA) to vSmartMOM convention (TOA→BOA).
"""
struct GeosChemSource <: NetCDFSource
    path::String
    idx::Int
    idy::Int
    idf::Int
    
    function GeosChemSource(path::AbstractString, idx::Integer, idy::Integer, idf::Integer)
        _require_source(isfile(path), "File not found: $path")
        _require_source(idx > 0, "idx must be positive")
        _require_source(idy > 0, "idy must be positive")
        _require_source(1 ≤ idf ≤ 6, "idf must be between 1 and 6 (cubed-sphere face)")
        new(String(path), Int(idx), Int(idy), Int(idf))
    end
end

"""
    NetCDFGridSource(path::String, lat_idx::Int, lon_idx::Int)

Generic IO source for gridded NetCDF atmospheric data (WRF, GCHP, etc.).

# Fields
- `path::String`: Path to the NetCDF file
- `lat_idx::Int`: Latitude grid index
- `lon_idx::Int`: Longitude grid index

# Usage
```julia
src = NetCDFGridSource("wrfout.nc", 50, 100)
params = parameters_from_yaml(src)
```

# Notes
This is a generic source for rectilinear gridded data. For specialized formats
(like GEOSChem with cubed-sphere), use the specific source type instead.
"""
struct NetCDFGridSource <: NetCDFSource
    path::String
    lat_idx::Int
    lon_idx::Int
    
    function NetCDFGridSource(path::AbstractString, lat_idx::Integer, lon_idx::Integer)
        _require_source(isfile(path), "File not found: $path")
        _require_source(lat_idx > 0, "lat_idx must be positive")
        _require_source(lon_idx > 0, "lon_idx must be positive")
        new(String(path), Int(lat_idx), Int(lon_idx))
    end
end

# Future extensibility examples (commented out for now):
# struct GCHPSource <: NetCDFSource ... end  # For GCHP variants
# struct WRFSource <: NetCDFSource ... end   # For WRF output
# struct CLMSource <: NetCDFSource ... end   # For CLM/CESM data

end # module Sources
