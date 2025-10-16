"""
    Aerosols Module

Flexible aerosol framework supporting multiple schemes:
- TOMAS-15: Size-resolved aerosol microphysics (15 bins)
- Two-Moment: Bulk aerosol properties (AOD + effective radius)
- Future: MAAM, custom schemes, etc.

Includes wavelength-dependent refractive index database and optical property calculations.
"""
module Aerosols

using NCDatasets
using YAML
using Interpolations
using Statistics
using Printf

# Export main types
export AerosolScheme, TOMAS15Scheme, TwoMomentScheme
export AerosolData, AerosolSpeciesData
export RefractiveIndexLUT, RefractiveIndexDatabase

# Export main functions
export read_aerosol_data
export load_refractive_index_database
export get_refractive_index
export compute_optical_properties

# Include submodules
include("types.jl")
include("refractive_index.jl")
include("readers.jl")
include("schemes/tomas15.jl")
include("schemes/two_moment.jl")
include("optical_properties.jl")

end # module
