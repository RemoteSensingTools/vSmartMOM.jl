# GEOSChem Integration for vSmartMOM

This document describes the new clean API for integrating GEOSChem NetCDF output files with vSmartMOM radiative transfer calculations.

## Overview

vSmartMOM now supports reading atmospheric profiles and trace gas concentrations directly from GEOSChem NetCDF4 output files using a type-safe, extensible IOSource system. This eliminates the need for manual data extraction and conversion.

## Architecture

The implementation follows Julia best practices with multiple dispatch:

```
src/IO/
├── Sources.jl              # IOSource type hierarchy
├── NetCDF/
│   ├── NetCDF.jl          # NetCDF utilities module
│   └── GeosChem.jl        # GEOSChem-specific reader
└── ...
```

### Type Hierarchy

```julia
IOSource (from Formats.jl)
├── FileSource              # YAML/TOML files
├── DictSource              # Direct Dict
└── NetCDFSource            # Abstract for NetCDF data
    ├── GeosChemSource      # GEOSChem cubed-sphere
    └── NetCDFGridSource    # Generic gridded data
```

## Usage

### Quick Start

```julia
using vSmartMOM

# Create a data source
src = GeosChemSource("GEOSChem.Custom.20190101_0000z.nc4", 10, 20, 1)

# Option 1: Use with default RT parameters
params = read_parameters(src)
model = model_from_parameters(params)
R = rt_run(model)

# Option 2: Get Dict and customize
config = geoschem_to_dict(src)
config["radiative_transfer"]["spec_bands"] = [(1e7/777):0.015:(1e7/757)]
params = read_parameters(config)
```

### Complete Example

```julia
using vSmartMOM

# Specify GEOSChem grid location
file = "GEOSChem.Custom.20190101_0000z.nc4"
idx, idy, idf = 10, 20, 1  # X, Y indices and cubed-sphere face

# Create source
src = GeosChemSource(file, idx, idy, idf)

# Read and convert to configuration Dict
config = geoschem_to_dict(src)
# This automatically extracts:
#   - Temperature profile (T)
#   - Pressure levels (p)
#   - Specific humidity (q)
#   - Trace gas VMRs (CO2, CO, CH4, N2O, C2H6, H2O)
#   - Geographic location (lat, lon)

# Add radiative transfer parameters
config["radiative_transfer"] = Dict(
    "spec_bands" => [(1e7/777):0.015:(1e7/757)],
    "surface" => ["LambertianSurfaceScalar(0.15)"],
    "quadrature_type" => "GaussQuadFullSphere()",
    "polarization_type" => "Stokes_I()",
    "max_m" => 3,
    "Δ_angle" => 2.0,
    "l_trunc" => 20,
    "depol" => 0.0,
    "float_type" => "Float64",
    "architecture" => "default_architecture"
)

config["geometry"] = Dict(
    "sza" => 60.0,
    "vza" => [0.0, 30.0, 60.0],
    "vaz" => [0.0, 0.0, 0.0],
    "obs_alt" => 1000.0
)

# Convert to parameters and run
params = read_parameters(config)
model = model_from_parameters(params)
R = rt_run(model)
```

### Convenience Function (See examples/geoschem_integration.jl)

```julia
include("examples/geoschem_integration.jl")

# Simple usage with defaults
R = run_rt_with_geoschem_v2("GEOSChem.Custom.20190101_0000z.nc4", 10, 20, 1)

# With custom parameters
R = run_rt_with_geoschem_v2(
    "GEOSChem.Custom.20190101_0000z.nc4", 
    10, 20, 1;
    spec_bands = [(1e7/777):0.015:(1e7/757)],
    sza = 45.0,
    vza = [0.0, 30.0, 60.0],
    architecture = GPU()
)
```

## API Reference

### Types

#### `GeosChemSource`

```julia
GeosChemSource(path::String, idx::Int, idy::Int, idf::Int)
```

Represents a GEOSChem data source at a specific cubed-sphere grid location.

**Arguments:**
- `path`: Path to GEOSChem NetCDF4 file
- `idx`: X-dimension index (must be > 0)
- `idy`: Y-dimension index (must be > 0)
- `idf`: Cubed-sphere face (1-6)

**Validation:**
- File existence is checked
- Indices must be positive
- Face must be 1-6

### Functions

#### `geoschem_to_dict`

```julia
geoschem_to_dict(src::GeosChemSource) -> Dict
```

Read GEOSChem file and convert to vSmartMOM configuration Dict.

**Returns:** Dict with keys:
- `"atmospheric_profile"`: T, p, q, vmr, profile_reduction
- `"absorption"`: molecules, vmr, broadening, CEF, wing_cutoff
- `"_metadata"`: source_file, latitude, longitude, time, etc.

#### `read_geoschem_profile`

```julia
read_geoschem_profile(file::String, idx::Int, idy::Int, idf::Int) -> Dict
```

Convenience function that creates a `GeosChemSource` and calls `geoschem_to_dict`.

## Data Conventions

### GEOSChem File Structure

GEOSChem output files have dimensions: `(Xdim, Ydim, nf, lev, time)`

- **nf**: Cubed-sphere face (1-6)
- **lev**: Vertical levels (typically 72)
  - Level 1 = Bottom of Atmosphere (BOA)
  - Level 72 = Top of Atmosphere (TOA)
- **time**: Time dimension

### Variable Mapping

| GEOSChem Variable | vSmartMOM Field | Units | Notes |
|------------------|----------------|-------|-------|
| `Met_T` | `T` | K | Temperature profile |
| `Met_PS2WET` | Surface pressure | hPa | Used to compute pressure grid |
| `Met_DELP` | Layer thickness | hPa | Used to compute pressure grid |
| `Met_SPHU` | `q` | g/kg | Specific humidity |
| `SpeciesConcVV_*` | `vmr[molecule]` | mol/mol | Trace gas VMRs |
| `lats`, `lons` | Metadata | degrees | Geographic location |

### Automatic Processing

1. **Vertical Flipping**: GEOSChem uses BOA→TOA indexing; vSmartMOM uses TOA→BOA. Data is automatically reversed.

2. **Pressure Grid**: Computed from surface pressure (`Met_PS2WET`) and layer thicknesses (`Met_DELP`).

3. **Trace Gases**: Automatically extracts available molecules:
   - N2O, CH4, C2H6, CO2, CO, H2O
   - O2 is added to molecule list but uses vSmartMOM defaults (not in GEOSChem output)

4. **Units**: All conversions handled automatically to match vSmartMOM conventions.

## Extensibility

The IOSource pattern makes it easy to add support for other formats:

```julia
# Example: Add WRF support
struct WRFSource <: NetCDFSource
    path::String
    lat_idx::Int
    lon_idx::Int
    time_idx::Int
end

function geoschem_to_dict(src::WRFSource)
    # WRF-specific reading logic
    # ...
    return config_dict
end
```

## Comparison with Old Approach

### Old (Hacky) Approach
```julia
# Manual struct definition
struct GeosChemData
    data::Dict{String, Any}
    units::Dict{String, String}
end

# Manual reading and conversion
function read_gchp(file, idx, idy, idf)
    # ... lots of manual NetCDF reading ...
    return GeosChemData(data, units)
end

# Manual RT setup
geos = read_gchp(file, idx, idy, idf)
parameters = default_parameters()
parameters.T = geos.data["Met_T"]
parameters.p = geos.data["pressure"]
# ... manual field assignment ...
```

### New (Clean) Approach
```julia
# Type-safe source specification
src = GeosChemSource(file, idx, idy, idf)

# Automatic conversion
params = read_parameters(src)

# Or customize
config = geoschem_to_dict(src)
# ... modify config ...
params = read_parameters(config)
```

## Benefits

- **Type Safety**: IOSource types catch errors at compile time
- **Extensibility**: Easy to add WRF, GCHP, CLM, etc.
- **Clean Separation**: I/O logic separated from RT logic
- **No Duplication**: Reuses existing parameter system
- **Idiomatic Julia**: Multiple dispatch, clear abstractions
- **Backwards Compatible**: Old code still works
- **Well Documented**: Types, functions, examples
- **Testable**: Easy to mock sources for testing  

## Dependencies

- **NCDatasets.jl**: NetCDF file reading (already in Project.toml)

## Testing

See `examples/geoschem_integration.jl` for comprehensive examples and test patterns.

## Future Enhancements

Potential additions:
- [ ] WRF output support (`WRFSource`)
- [ ] GCHP variants (`GCHPSource`)
- [ ] CLM/CESM data (`CLMSource`)
- [ ] Caching for repeated reads
- [ ] Lazy loading for large files
- [ ] Interpolation to arbitrary lat/lon
- [ ] Time series support
