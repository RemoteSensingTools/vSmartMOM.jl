module RadiativeTransfer
using Pkg.Artifacts
using LinearAlgebra
using Distributions
using CUDA

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
include("CrossSection/CrossSection.jl")

# Mie Phase Function module:
include("PhaseFunction/PhaseFunction.jl")

# Perform some GPU setup when the module is loaded
function __init__()
    @hascuda begin
        @info "CUDA-enabled GPU(s) detected:"
        for (gpu, dev) in enumerate(CUDA.devices())
            @info "$dev: $(CUDA.name(dev))"
        end
	CUDA.allowscalar(false)
    end
end

end
