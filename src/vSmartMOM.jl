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
# Export GEOSChem/NetCDF integration
export GeosChemSource, NetCDFGridSource, NetCDFSource, geoschem_to_dict, read_geoschem_profile

using .Architectures
using .Absorption

# Module initialization - GPU setup happens in CUDAExt when CUDA is loaded
function __init__()
    # Small delay to let CUDA extension initialize first if it's loading
    sleep(0.05)
    
    # Only show CPU message if CUDA is truly not available
    # (CUDAExt will show its own message if GPU is available)
    if !has_cuda()
        # Pretty startup message for CPU-only mode
        println()
        println("┌─────────────────────────────────────────────────────────")
        println("│ ⚡ vSmartMOM Radiative Transfer")
        println("├─────────────────────────────────────────────────────────")
        println("│  Execution Mode: CPU")
        println("│  Backend:        KernelAbstractions.CPU()")
        println("│")
        println("│  💡 For GPU acceleration:")
        println("│     julia> using Pkg; Pkg.add(\"CUDA\")")
        println("│     (Requires NVIDIA GPU with compatible drivers)")
        println("└─────────────────────────────────────────────────────────")
        println()
    end
end

end
