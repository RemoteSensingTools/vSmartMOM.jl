# CUDA/GPU Setup Guide for vSmartMOM.jl

## Overview

vSmartMOM.jl supports both **CPU** and **GPU** (CUDA) execution. CUDA support is **completely optional** - the package works perfectly fine on CPU-only systems.

## Architecture Selection

vSmartMOM uses `KernelAbstractions.jl` for hardware abstraction, allowing the same code to run on both CPU and GPU.

### Automatic Architecture Detection

By default, vSmartMOM automatically detects the best available architecture:

```julia
using vSmartMOM

# Automatically selects GPU if available, otherwise CPU
arch = default_architecture()  # Returns GPU() or CPU()
```

### Manual Architecture Selection

You can explicitly specify the architecture:

```julia
# Force CPU execution (always available)
using vSmartMOM
params = default_parameters()
# ... configure params ...
model = model_from_parameters(params, arch=CPU())

# Request GPU execution (requires CUDA)
model = model_from_parameters(params, arch=GPU())
```

## GPU/CUDA Requirements

### What Gets Checked

vSmartMOM performs comprehensive CUDA availability checks:

1. ✅ **CUDA.jl installed** - Is the Julia CUDA package available?
2. ✅ **CUDA runtime libraries** - Are CUDA libraries installed on the system?
3. ✅ **GPU hardware present** - Is there a CUDA-capable NVIDIA GPU?
4. ✅ **Driver compatibility** - Is the GPU driver version compatible?
5. ✅ **Toolkit version** - Does the CUDA toolkit match Julia's requirements?
6. ✅ **Functional test** - Can we actually allocate and use GPU memory?

If ANY of these checks fail, vSmartMOM automatically falls back to CPU mode.

### Installing CUDA Support

```julia
using Pkg
Pkg.add("CUDA")
```

Then restart Julia and load vSmartMOM:

```julia
using vSmartMOM
# You should see: "vSmartMOM CUDA extension loaded successfully"
# Plus info about detected GPUs
```

## Common Issues and Solutions

### Issue 1: "CUDA.jl reported functional GPU but encountered an error"

**Cause**: CUDA toolkit version mismatch with your GPU driver or Julia's CUDA.jl version.

**Solution**:
```julia
using CUDA
CUDA.versioninfo()  # Check versions

# Update CUDA.jl to latest
using Pkg
Pkg.update("CUDA")

# Or install a specific version if needed
Pkg.add(name="CUDA", version="5.5")
```

### Issue 2: "No CUDA-capable GPU detected"

**Causes**:
- No NVIDIA GPU in the system
- GPU not recognized by drivers
- GPU in use by another process

**Solution**:
```bash
# Check if GPU is visible to system
nvidia-smi

# If nvidia-smi works, check Julia can see it:
```
```julia
using CUDA
CUDA.has_cuda_gpu()  # Should return true
```

### Issue 3: Package loads but GPU not being used

**Check**:
```julia
using vSmartMOM

# Check what architecture is being used
arch = default_architecture()
println("Using: ", arch)  # Should show GPU() if CUDA works

# Check if CUDA extension loaded
println("CUDA available: ", vSmartMOM.Architectures.has_cuda())
```

### Issue 4: LD_LIBRARY_PATH warnings

If you see warnings about CUDA libraries loaded from system paths:

```
WARNING: CUDA runtime library `libcublasLt.so.12` was loaded from a system path
```

**This is usually not a problem**, but to clean it up:

```bash
# Unset LD_LIBRARY_PATH for CUDA libraries
unset LD_LIBRARY_PATH

# Or remove CUDA paths specifically
export LD_LIBRARY_PATH=$(echo $LD_LIBRARY_PATH | tr ':' '\n' | grep -v cuda | tr '\n' ':')
```

## Performance Considerations

### When to Use GPU

✅ **GPU is faster for**:
- Large spectral bands (>100 wavelengths)
- High angular resolution (many viewing angles)
- Multiple aerosol layers
- Batch processing multiple scenarios

⚠️ **CPU might be comparable for**:
- Small problems (<50 wavelengths, few angles)
- Single scenarios
- Systems with high-end CPUs

### Benchmarking

```julia
using vSmartMOM, BenchmarkTools

params = default_parameters()
# ... configure params ...

# CPU
model_cpu = model_from_parameters(params, arch=CPU())
@btime rt_run($model_cpu)

# GPU (if available)
model_gpu = model_from_parameters(params, arch=GPU())
@btime rt_run($model_gpu)
```

## Technical Details

### Package Extension Architecture

vSmartMOM uses Julia 1.9+ package extensions (weakdeps) for optional CUDA support:

- **Without CUDA.jl**: Package loads normally, only CPU support available
- **With CUDA.jl**: Extension `vSmartMOMCUDAExt` loads automatically, enabling GPU support

This means:
- No CUDA dependency for CPU-only users
- No import errors on systems without GPUs
- CUDA code only compiled when actually needed

### Float Precision

GPU performance benefits from using `Float32`:

```julia
# Float32 (faster on GPU, lower memory)
params = default_parameters(FT=Float32)

# Float64 (default, more accurate, slower on GPU)
params = default_parameters(FT=Float64)
```

### Memory Management

GPU memory is limited. For large problems:

```julia
using CUDA

# Check available GPU memory
CUDA.available_memory()

# Monitor memory usage
@time model = model_from_parameters(params, arch=GPU())
CUDA.memory_status()
```

## Troubleshooting Command Summary

```julia
# Check CUDA availability
using CUDA
CUDA.functional()        # Overall status
CUDA.has_cuda()          # Runtime available?
CUDA.has_cuda_gpu()      # GPU detected?
CUDA.versioninfo()       # Detailed version info

# Check vSmartMOM CUDA status
using vSmartMOM
vSmartMOM.Architectures.has_cuda()  # Is CUDA support loaded?
default_architecture()               # What will be used?

# Test GPU functionality
using CUDA
test = CUDA.CuArray([1.0f0, 2.0f0, 3.0f0])
CUDA.@sync sum(test)  # Should return 6.0 if GPU works
```

## Getting Help

If you continue to have CUDA-related issues:

1. **Collect diagnostic information**:
   ```julia
   using CUDA, vSmartMOM
   CUDA.versioninfo()
   println("vSmartMOM CUDA available: ", vSmartMOM.Architectures.has_cuda())
   ```

2. **Check system CUDA installation**:
   ```bash
   nvidia-smi
   nvcc --version  # If CUDA toolkit installed
   ```

3. **Open an issue** at https://github.com/RemoteSensingTools/vSmartMOM.jl/issues
   - Include the diagnostic information above
   - Mention your OS, Julia version, and GPU model

## Summary

- ✅ **CUDA is optional** - vSmartMOM works great on CPU
- ✅ **Automatic detection** - No configuration needed in most cases
- ✅ **Robust fallback** - Failures automatically switch to CPU
- ✅ **Clear error messages** - Know exactly what's needed for GPU
- ✅ **No dependencies** - CUDA.jl only loaded when installed and working
