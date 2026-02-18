# # GPU Acceleration Tutorial
#
# ### Introduction
# vSmartMOM can run the RT solver on NVIDIA GPUs via CUDA.jl, with
# no changes to your input YAML files. This tutorial covers:
#
# 1. How to switch between CPU and GPU
# 2. Architecture-aware array handling
# 3. Performance considerations
# 4. Common pitfalls
#
# ---
#
# ### Load packages

using vSmartMOM

# Check whether CUDA is available at runtime
const CAN_USE_GPU = try
    @eval using CUDA
    CUDA.functional()
catch
    false
end

println("CUDA available: $CAN_USE_GPU")

# ---
#
# ### 1. Switching between CPU and GPU
#
# The architecture is set via `params.architecture`. All array allocations
# and kernel launches automatically adapt.

params = parameters_from_yaml("test/test_parameters/PureRayleighParameters.yaml")

## CPU run
params.architecture = vSmartMOM.Architectures.CPU()
model_cpu = model_from_parameters(params)
R_cpu, T_cpu, _, _, _, _, _ = rt_run(model_cpu)
println("CPU result shape: ", size(R_cpu), "  eltype: ", eltype(R_cpu))

## GPU run (only if CUDA is available)
if CAN_USE_GPU
    params.architecture = vSmartMOM.Architectures.GPU()
    model_gpu = model_from_parameters(params)
    R_gpu, T_gpu, _, _, _, _, _ = rt_run(model_gpu)
    println("GPU result shape: ", size(R_gpu), "  eltype: ", eltype(R_gpu))

    ## Compare CPU vs GPU (should agree to machine precision)
    max_diff = maximum(abs.(R_cpu .- R_gpu))
    println("Max CPU-GPU difference: $max_diff")
else
    println("Skipping GPU run — CUDA not available on this system.")
end

# ---
#
# ### 2. Architecture-aware arrays
#
# Internally, vSmartMOM uses `array_type(architecture)` to select the
# array constructor:
# - `CPU()` → `Array`
# - `GPU()` → `CuArray`
#
# When writing code that works on both, use these patterns:
#
# ```julia
# arr_type = array_type(architecture)
# x = arr_type(zeros(FT, n))          # allocate on the right device
# y = collect(gpu_array)              # safely copy GPU→CPU (no-op on CPU)
# z = Matrix(Diagonal(v))            # dense from diagonal (always CPU)
# ```
#
# **Avoid** calling `Array(gpu_array)` directly — use `collect()` instead.
# This makes intent explicit and works identically on CPU and GPU.

# ---
#
# ### 3. Performance considerations
#
# **Batch size matters.** The RT solver uses batched matrix operations
# (`batched_mul` from NNlib.jl). GPU performance scales with the number
# of spectral points (`nSpec`) — a batch of 100+ wavelengths is needed
# to saturate modern GPUs.
#
# **Memory.** Each atmospheric layer allocates `r, t, j` matrices of shape
# `(nμ, nμ, nSpec)`. With full Stokes polarization (`nμ = 4 * Nquad`)
# and many wavelengths, GPU memory can become a bottleneck.
#
# **First call overhead.** Julia compiles GPU kernels on first use.
# The second call is much faster. Use a small warm-up run.

# ---
#
# ### 4. Common pitfalls
#
# **Scalar indexing.** Accessing individual elements of a `CuArray`
# triggers slow GPU-to-CPU transfers. The code avoids this by using
# broadcasted operations and KernelAbstractions.jl kernels for
# element-wise work.
#
# **Non-contiguous slices.** `batched_mul` requires contiguous memory.
# If you slice a 4D array as `A[i,:,:,:]`, the result may not be
# contiguous. The code uses `copy()` in these cases to ensure
# contiguity.
#
# **Type stability.** GPU compilation is sensitive to type instabilities.
# The RT kernels enforce `FT<:Real` as the element type, so
# `ForwardDiff.Dual` numbers never reach CUDA kernels.

if CAN_USE_GPU
    println("\n--- GPU device info ---")
    println("Device:  ", CUDA.name(CUDA.device()))
    println("Memory:  ", round(CUDA.total_memory() / 1e9, digits=1), " GB")
end
