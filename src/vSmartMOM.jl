#=
 
This file is the entry-point for the vSmartMOM module. 

It includes this module's source files and exports the relevant keywords.  
 
=#

module vSmartMOM
using Pkg.Artifacts
using LinearAlgebra
using Distributions
using Parameters
using DocStringExtensions
using UnPack
using UnicodePlots

# Export Architecture functions
export CPU, GPU, default_architecture, array_type

# Export the artifact convenience function
export artifact

# RT mode types (forward vs linearized)
abstract type RT_Mode end
struct FwdMode <: RT_Mode end
struct LinMode <: RT_Mode end
export FwdMode, LinMode

# GPU/CPU Architecture (from Oceanigans)
include("Architectures.jl")
using .Architectures

# Artifacts
include("Artifacts/artifact_helper.jl")

# Absorption Cross Section module:
include("Absorption/Absorption.jl")

# Mie Phase Function module:
include("Scattering/Scattering.jl")

# Inelastic Scattering module:
include("Inelastic/InelasticScattering.jl")

# CoreRT module:
include("CoreRT/CoreRT.jl")
using .CoreRT

# SolarModel module:
include("SolarModel/SolarModel.jl")

# IO submodule (must come after CoreRT types are defined)
include("IO/IO.jl")
using .IO

# Export some vSmartMOM functions
export default_parameters, parameters_from_yaml, model_from_parameters, rt_run, read_parameters, read_atmos_profile
# Export linearized RT functions
export rt_run_lin, model_from_parameters_lin
# Export new hierarchical model types
export RTModel, AbstractRTModel, SolverConfig, Atmosphere, RayleighScattering, AerosolState, Optics, OpticsLin
# Export GEOSChem/NetCDF integration
export GeosChemSource, NetCDFGridSource, NetCDFSource, geoschem_to_dict, read_geoschem_profile

using .Architectures
using .Absorption

# Module initialization - GPU setup happens in CUDAExt when CUDA is loaded
function __init__()
    # Nothing to do here; GPU detection is handled by the CUDA extension
end

end
