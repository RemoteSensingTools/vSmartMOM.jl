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

export read_parameters, parameters_from_yaml, parameters_from_dict,
       read_atmos_profile, read_atmos_profile_dict

end # module
