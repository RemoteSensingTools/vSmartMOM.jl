"""
    vSmartMOMCUDAExt

Package extension that loads CUDA support when CUDA.jl is available.
This allows vSmartMOM to work without CUDA on systems without compatible GPUs.

Requires Julia 1.9+ with package extensions feature.
"""
module vSmartMOMCUDAExt

using vSmartMOM
using vSmartMOM.Architectures
using vSmartMOM.CoreRT
using CUDA
using KernelAbstractions

# Re-export CUDAKernels for code that needs CUDADevice
import CUDA.CUDAKernels

# Extend GPU device to use CUDA backend
Architectures.devi(::vSmartMOM.Architectures.GPU) = CUDA.CUDABackend(; always_inline=true)

# Extend array_type to return CuArray for GPU
Architectures.array_type(::vSmartMOM.Architectures.GPU) = CuArray

# Extend architecture detection for CuArrays
Architectures.architecture(::CuArray) = vSmartMOM.Architectures.GPU()

# No need to override synchronize_if_gpu - it checks has_cuda() and calls CUDA.synchronize()
# which will be available from this extension

# Include GPU-specific batched operations
include("gpu_batched_cuda.jl")

# Module initialization - called when extension is loaded
function __init__()
    # CUDA.functional() checks:
    # 1. CUDA libraries are available
    # 2. A compatible GPU is present
    # 3. The CUDA toolkit version matches Julia's requirements
    # 4. The driver version is compatible
    
    if CUDA.functional()
        try
            # Make CUBLAS available to CoreRT for batched GPU operations (make_added_layer, doubling)
            CoreRT.CUBLAS_ref[] = CUDA.CUBLAS
            # Additional safety check: try to actually use CUDA
            # This catches edge cases where CUDA.functional() returns true
            # but operations fail due to version mismatches
            test_arr = CUDA.CuArray([1.0f0])
            CUDA.synchronize()
            
            # If we got here, CUDA really works!
            Architectures._has_cuda[] = true
            Architectures._sync_gpu[] = CUDA.synchronize

            CUDA.allowscalar(false)
            
        catch e
            @warn "vSmartMOM GPU initialization failed, falling back to CPU" exception=e
            Architectures._has_cuda[] = false
        end
    else
        Architectures._has_cuda[] = false
    end
end

end # module vSmartMOMCUDAExt
