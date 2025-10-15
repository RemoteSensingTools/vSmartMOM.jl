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

Load a YAML configuration file into a `vSmartMOM_Parameters` object.

- Dispatches to the formats registry based on file extension (currently YAML).
- Validates required keys and types; supports optional absorption and scattering blocks.
- See docs: vSmartMOM → User-Defined RT Parameters and IO Guide for schema and examples.
"""

"Load parameters from a file path using the formats registry"
function parameters_from_yaml(file_path::AbstractString)
    cfg = Formats.load_config(Formats.FileSource(String(file_path)))
    return parameters_from_dict(cfg)
end

end # module
