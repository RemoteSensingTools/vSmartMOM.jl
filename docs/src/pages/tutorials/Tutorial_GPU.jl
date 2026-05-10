# # GPU Acceleration
#
# vSmartMOM can run the RT solver on optional GPU backends. CUDA.jl is the
# mature NVIDIA path. Metal.jl support is experimental for Apple Silicon
# Float32 runs and starts with the core batched multiply/inverse path.
#
# ## 1) Checking GPU availability

using vSmartMOM

if vSmartMOM.Architectures.has_cuda()
    println("CUDA is available — GPU mode enabled.")
elseif vSmartMOM.Architectures.has_metal()
    println("Metal is available — Apple GPU mode enabled.")
else
    println("CUDA not available — running CPU-only examples.")
end

# The default architecture auto-detects:
arch = vSmartMOM.Architectures.default_architecture()
println("Default architecture: ", arch)

# ## 2) Selecting architecture
#
# You can explicitly choose CPU or GPU:
#
# ```julia
# params.architecture = vSmartMOM.Architectures.CPU()   # force CPU
# params.architecture = vSmartMOM.Architectures.GPU()   # force CUDA GPU
# params.architecture = vSmartMOM.Architectures.MetalGPU()  # force Metal GPU
# ```
#
# When loading from YAML, set `architecture: GPU` for CUDA or
# `architecture: MetalGPU` for Apple Silicon Metal, or override after loading:

params = read_parameters(
    joinpath(pkgdir(vSmartMOM),
             "test", "test_parameters", "PureRayleighParameters.yaml"))

params.architecture = vSmartMOM.Architectures.CPU()
model_cpu = model_from_parameters(params)
R_cpu, = rt_run(model_cpu)
println("CPU result shape: ", size(R_cpu))

# ## 3) GPU execution
#
# When a GPU is available, simply switch the architecture:
#
# ```julia
# params.architecture = vSmartMOM.Architectures.GPU()
# model_gpu = model_from_parameters(params)
# R_gpu, = rt_run(model_gpu)
# ```
#
# The same code paths are used — arrays are automatically moved to GPU
# via backend arrays (`CuArray` or `MtlArray`), and the RT kernels use
# `KernelAbstractions.jl` for device-portable kernel dispatch.

# ## 4) What runs on GPU
#
# The following operations are fully GPU-accelerated:
#
# | Component                | GPU support |
# |:-------------------------|:------------|
# | Elemental layer (`elemental!`) | ✅ Full |
# | Doubling (`doubling!`)         | ✅ Full |
# | Interaction (`interaction!`)   | ✅ Full |
# | Batched matrix inverse         | ✅ CUBLAS on CUDA; portable KA kernel on Metal, with a local-memory guard |
# | Linearized/Jacobian runs       | ✅ CUDA/CPU; Metal not yet validated |
# | Phase function computation     | ❌ CPU only |
# | Absorption cross-sections      | ✅ Full |
# | Postprocessing (VZA interp.)   | ✅ Full |
#
# Phase function computation (Mie/NAI2) is always done on CPU and then
# transferred to GPU.  This is typically a one-time cost per aerosol type.
#
# ## 5) GPU-specific YAML parameters
#
# For GPU runs, a smaller spectral grid often makes sense to fit in
# GPU memory.  The `test/test_parameters/O2Parameters_GPU.yaml` config
# demonstrates a reduced grid (nSpec=60 vs 6837 for CPU):
#
# ```yaml
# radiative_transfer:
#   architecture: GPU
#   spec_bands:
#     - "(1e7/771):0.2:(1e7/759)"   # coarser grid
# ```
#
# On Apple Silicon, use Float32 and the Metal architecture:
#
# ```yaml
# radiative_transfer:
#   float_type: Float32
#   architecture: MetalGPU
# ```
#
# Metal support currently targets modest stream/Stokes dimensions. The batched
# inverse path checks the backend local-memory budget before launch and errors
# clearly when a scene is too large for the shared-memory LU kernel. With the
# current 32 KiB guard, Float32 matrices with `N = Nquad * nStokes >= 64` are
# rejected.

# ## 6) Benchmarking CPU vs GPU
#
# The benchmark script at `test/benchmark_rt.jl` provides timing
# comparisons.  A typical A100 GPU speedup for the O₂-A band:
#
# ```
# Forward noRS (nSpec=17):  CPU ~0.3s, GPU ~0.05s  (6× speedup)
# Forward noRS (nSpec=60):  CPU ~1.0s, GPU ~0.08s  (12× speedup)
# ```
#
# Speedup increases with spectral resolution because the GPU excels at
# the batched matrix operations that scale with `nSpec`.

# ## 7) Memory considerations
#
# GPU memory usage scales as `O(nμ² × nSpec)` for the layer matrices.
# With `nμ = 18` (9 streams, Stokes-I) and `nSpec = 6000`, each 3D
# matrix is about 15 MB in Float64.  The solver needs ~20 such matrices
# simultaneously, totaling ~300 MB — well within a 40 GB A100.
#
# For very high spectral resolution, consider splitting into sub-bands
# or using Float32 (supported but with reduced accuracy).
