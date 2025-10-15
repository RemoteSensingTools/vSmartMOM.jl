# IO submodule: centralizes input/output for vSmartMOM
module IO

using YAML
using DocStringExtensions
using ..vSmartMOM
using ..Absorption
using ..Scattering
using ..Architectures
using ..CoreRT

include("Parameters.jl")
include("AtmosProfile.jl")
include("Formats.jl")
using .Formats

export read_parameters, parameters_from_yaml, parameters_from_dict,
       read_atmos_profile, read_atmos_profile_dict

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
