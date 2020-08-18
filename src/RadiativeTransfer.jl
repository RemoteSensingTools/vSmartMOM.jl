module RadiativeTransfer

export Architectures, CPU, GPU,device

# Export the Cross Section models
export HitranModel, InterpolationModel


using CUDA
using KernelAbstractions

# GPU/CPU Architecture (from Oceanigans)
include("Architectures.jl")

# Absorption Cross Section stuff:
include("CrossSection/CrossSection.jl")

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
