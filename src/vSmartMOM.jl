#=
 
This file is the entry-point for the vSmartMOM module. 

It includes this module's source files and exports the relevant keywords.  
 
=#

module vSmartMOM
using Pkg.Artifacts
using LinearAlgebra
using Distributions
using CUDA
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

# Export some vSmartMOM functions
export default_parameters, parameters_from_yaml, model_from_parameters, rt_run
export lin_model_from_parameters
using .Architectures
using .Absorption

# Perform some GPU setup when the module is loaded
function __init__()
    @hascuda begin
        @info "CUDA-enabled GPU(s) detected"
        for (gpu, dev) in enumerate(CUDA.devices())
            @info "$dev: $(CUDA.name(dev))"
        end
	CUDA.allowscalar(false)
    end
end

end
