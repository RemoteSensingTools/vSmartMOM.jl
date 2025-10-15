"""
    vSmartMOMCUDAExt

Package extension that loads CUDA support when CUDA.jl is available.
This allows vSmartMOM to work without CUDA on systems without compatible GPUs.

Requires Julia 1.9+ with package extensions feature.
"""
module vSmartMOMCUDAExt

using vSmartMOM
using vSmartMOM.Architectures
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
            # Additional safety check: try to actually use CUDA
            # This catches edge cases where CUDA.functional() returns true
            # but operations fail due to version mismatches
            test_arr = CUDA.CuArray([1.0f0])
            CUDA.synchronize()
            
            # If we got here, CUDA really works!
            Architectures._has_cuda[] = true
            Architectures._sync_gpu[] = CUDA.synchronize
            
            # Pretty startup message
            println()
            println("┌─────────────────────────────────────────────────────────")
            println("│ 🚀 vSmartMOM GPU Acceleration: ENABLED")
            println("├─────────────────────────────────────────────────────────")
            println("│  CUDA Runtime:  v$(CUDA.runtime_version())")
            println("│  CUDA Driver:   v$(CUDA.driver_version())")
            println("│  Detected GPUs: $(length(CUDA.devices()))")
            for (gpu, dev) in enumerate(CUDA.devices())
                mem_gb = round(CUDA.totalmem(dev) / 1024^3, digits=1)
                println("│    ├─ GPU $gpu: $(CUDA.name(dev)) ($(mem_gb) GB)")
            end
            println("│")
            println("│  💡 Tip: Use CPU() to force CPU execution if needed")
            println("└─────────────────────────────────────────────────────────")
            println()
            
            CUDA.allowscalar(false)
            
        catch e
            # CUDA failed functional test
            println()
            println("┌─────────────────────────────────────────────────────────")
            println("│ ⚠️  vSmartMOM GPU Acceleration: UNAVAILABLE")
            println("├─────────────────────────────────────────────────────────")
            println("│  CUDA.jl is installed but GPU initialization failed.")
            println("│  ")
            println("│  Common causes:")
            println("│   • CUDA toolkit/driver version mismatch")
            println("│   • Outdated GPU drivers")
            println("│   • GPU in use by another process")
            println("│  ")
            println("│  Falling back to CPU execution.")
            println("│  Error: $(sprint(showerror, e))")
            println("│")
            println("│  💡 Try: CUDA.versioninfo() for diagnostics")
            println("└─────────────────────────────────────────────────────────")
            println()
            Architectures._has_cuda[] = false
        end
    else
        # CUDA.functional() returned false
        reason = if !CUDA.has_cuda()
            "CUDA runtime library not found on system"
        elseif !CUDA.has_cuda_gpu()
            "No CUDA-capable GPU detected"
        else
            "CUDA available but not functional (version mismatch?)"
        end
        
        println()
        println("┌─────────────────────────────────────────────────────────")
        println("│ ℹ️  vSmartMOM GPU Acceleration: NOT AVAILABLE")
        println("├─────────────────────────────────────────────────────────")
        println("│  Reason: $reason")
        println("│  Using:  CPU backend (KernelAbstractions.CPU())")
        println("│")
        println("│  💡 For GPU support, ensure:")
        println("│     • NVIDIA GPU with CUDA support")
        println("│     • Compatible CUDA drivers installed")
        println("│     • CUDA toolkit matches Julia requirements")
        println("└─────────────────────────────────────────────────────────")
        println()
        Architectures._has_cuda[] = false
    end
end

end # module vSmartMOMCUDAExt
