module RadiativeTransfer
using Pkg.Artifacts

export Architectures, CPU, GPU,device

# Export the Cross Section models
export HitranModel, InterpolationModel

using LinearAlgebra
using Distributions
using CUDA
using KernelAbstractions

# GPU/CPU Architecture (from Oceanigans)
include("Architectures.jl")

# Absorption Cross Section module:
include("CrossSection/CrossSection.jl")

# Mie Phase Function module:
include("PhaseFunction/PhaseFunction.jl")

using .Architectures
using .CrossSection

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
