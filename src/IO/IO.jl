# IO submodule: centralizes input/output for vSmartMOM
module IO

using YAML
using NCDatasets
using DocStringExtensions
using ..vSmartMOM
using ..Absorption
using ..Scattering
using ..Architectures
using ..CoreRT

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

export read_parameters, parameters_from_yaml, parameters_from_dict,
       read_atmos_profile, read_atmos_profile_dict
# NetCDF/GEOSChem exports
export GeosChemSource, NetCDFGridSource, NetCDFSource
export geoschem_to_dict, read_geoschem_profile

"""
    parameters_from_yaml(path::AbstractString) -> vSmartMOM_Parameters

Load a YAML configuration file and return a fully populated
[`vSmartMOM_Parameters`](@ref) object.

Parses the `radiative_transfer`, `geometry`, `atmospheric_profile`,
`absorption`, and `scattering` YAML sections.  Surface types are resolved
via `BRDF_MAP`; quadrature, polarization, and broadening via their
respective maps in `Parameters.jl`.

# Arguments
- `path::AbstractString`: Path to the YAML configuration file.

# Returns
- `params::vSmartMOM_Parameters{FT}`: Parsed parameter struct ready for
  [`model_from_parameters`](@ref).

# Example
```julia
using vSmartMOM
params = parameters_from_yaml("config/my_scene.yaml")
params.architecture = vSmartMOM.Architectures.CPU()
model = model_from_parameters(params)
```

# See also
- [`parameters_from_dict`](@ref) for loading from an in-memory Dict.
- [`model_from_parameters`](@ref) to build the RT model.
"""
function parameters_from_yaml(file_path::AbstractString)
    cfg = Formats.load_config(Formats.FileSource(String(file_path)))
    return parameters_from_dict(cfg)
end


"Load parameters from an IOSource using the formats registry"
function parameters_from_yaml(src::Formats.IOSource)
    cfg = Formats.load_config(src)
    return parameters_from_dict(cfg)
end
end # module
