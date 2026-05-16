"""
    IO

Input/output helpers for loading YAML, dictionary, atmospheric-profile, and
NetCDF/GEOS-Chem configuration sources into vSmartMOM parameter objects.
"""
module IO

using YAML
using NCDatasets
using DocStringExtensions
using Printf
using ..vSmartMOM
using ..Absorption
using ..Scattering
using ..Architectures
using ..CoreRT
using ..Aerosols

# Load modules in dependency order
include("Formats.jl")
using .Formats
import .Formats: load_config  # Import for extension in GeosChem.jl

include("Parameters.jl")
include("AtmosProfile.jl")

include("Sources.jl")
using .Sources

# Include NetCDF readers directly (no nested module)
include("NetCDF/GeosChem.jl")
include("NetCDF/GCHPScene.jl")
include("Benchmark/netcdf_writer.jl")
include("Benchmark/aod_diagnostic.jl")

export read_parameters, parameters_from_file, parameters_from_source,
       parameters_from_yaml, parameters_from_dict,
       read_atmos_profile, read_atmos_profile_dict
# NetCDF/GEOSChem exports
export GeosChemSource, NetCDFGridSource, NetCDFSource
export geoschem_to_dict, read_geoschem_profile
# Open-once GCHP scene API
export GCHPFile, GCHPScene, scene_at, scenes,
       read_gchp_scene, scene_to_dict, parameters_from_scene,
       compute_scene_aod
# Benchmark writer + scene-loop driver
export write_scene_result, generate_benchmark,
       write_gchp_aod_diagnostic, write_gchp_aod_diagnostic_bulk

"""
    parameters_from_file(path::AbstractString) -> vSmartMOM_Parameters

Load a registered configuration file and return a fully populated
[`vSmartMOM_Parameters`](@ref) object.

The file format is selected by extension through the `Formats.load_config`
registry. Built-in extensions are `.yaml`, `.yml`, and `.toml`. The loaded
document must contain the same parameter schema accepted by
[`parameters_from_dict`](@ref).

# Arguments
- `path::AbstractString`: Path to a registered configuration file.

# Returns
- `params::vSmartMOM_Parameters{FT}`: Parsed parameter struct ready for
  [`model_from_parameters`](@ref).

# Example
```julia
using vSmartMOM
params = parameters_from_file("config/my_scene.toml")
params.architecture = vSmartMOM.Architectures.CPU()
model = model_from_parameters(params)
```

# See also
- [`parameters_from_source`](@ref) for typed sources such as GEOS-Chem NetCDF.
- [`parameters_from_dict`](@ref) for loading from an in-memory Dict.
- [`parameters_from_yaml`](@ref) for the legacy YAML-specific wrapper.
- [`model_from_parameters`](@ref) to build the RT model.
"""
function parameters_from_file(file_path::AbstractString)
    cfg = Formats.load_config(Formats.FileSource(String(file_path)))
    return parameters_from_dict(cfg)
end

"""
    parameters_from_source(src::Formats.IOSource) -> vSmartMOM_Parameters

Load parameters from a typed IO source.

This is the extension point for non-file or domain-specific inputs: implement
`Formats.load_config(src::MySource)` to return a configuration `Dict`, then
`parameters_from_source(src)` will reuse the shared parameter parser.
"""
function parameters_from_source(src::Formats.IOSource)
    cfg = Formats.load_config(src)
    return parameters_from_dict(cfg)
end

"""
    parameters_from_yaml(path::AbstractString) -> vSmartMOM_Parameters

Compatibility wrapper for YAML parameter files.

New code should prefer [`parameters_from_file`](@ref), which also supports TOML
and any future registered file formats.
"""
function parameters_from_yaml(file_path::AbstractString)
    ext = lowercase(splitext(String(file_path))[2])
    if ext ∉ (".yaml", ".yml")
        throw(ArgumentError("parameters_from_yaml only accepts .yaml/.yml files; use parameters_from_file for $(ext) input."))
    end
    return parameters_from_file(file_path)
end

parameters_from_yaml(src::Formats.IOSource) = parameters_from_source(src)
end # module
