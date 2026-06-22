module Architectures
# Shout out to https://github.com/CliMA/Oceananigans.jl/blob/master/src/Architectures.jl

export
    @hascuda,
    AbstractArchitecture, CPU, GPU, MetalGPU,
    devi, array_type,
    architecture,
    default_architecture,
    synchronize_if_gpu,
    has_cuda, has_metal,
    ka_backend, default_mie_precision_policy, has_gpu_mie

using KernelAbstractions

# Check if CUDA is available and functional
# This will be true when the CUDAExt extension is loaded
const _has_cuda = Ref(false)
has_cuda() = _has_cuda[]

# Check if Metal is available and functional.
# This will be true when the Metal extension is loaded on a supported Mac.
const _has_metal = Ref(false)
has_metal() = _has_metal[]

# Synchronization function reference - set by loaded GPU backend extensions.
const _sync_gpu = Ref{Function}(() -> nothing)
sync_device() = _sync_gpu[]()

"""
    AbstractArchitecture
Abstract supertype for compute architectures supported by vSmartMOM.
"""
abstract type AbstractArchitecture end

"""
    CPU <: AbstractArchitecture
Run on a single-core of a CPU.
"""
struct CPU <: AbstractArchitecture end

"""
    GPU <: AbstractArchitecture
Run on a single NVIDIA CUDA GPU.
"""
struct GPU <: AbstractArchitecture end

"""
    MetalGPU <: AbstractArchitecture

Run on an Apple Silicon GPU through Metal.jl. Metal support is loaded through
the optional `vSmartMOMMetalExt` package extension and requires Float32 arrays.
The current Metal batched inverse path uses a conservative 32 KiB threadgroup
memory guard, so large stream/Stokes matrices should use CPU or CUDA until a
global-memory fallback lands.
"""
struct MetalGPU <: AbstractArchitecture end

"""
    @hascuda expr
A macro to compile and execute `expr` only if CUDA is installed and available. Generally used to
wrap expressions that can only be compiled if `CuArrays` and `CUDAnative` can be loaded.
"""
macro hascuda(expr)
    return has_cuda() ? :($(esc(expr))) : :(nothing)
end

# CPU device always available
@inline devi(::CPU) = KernelAbstractions.CPU()

# GPU device - will be defined in CUDAExt when CUDA is loaded
# No fallback method defined to avoid precompilation errors

@inline architecture(::Array) = CPU()

"""
    ka_backend(arch::AbstractArchitecture)

Return the KernelAbstractions backend object for `arch`, suitable for launching
KA kernels (`@kernel`/`ndrange=`) and `KernelAbstractions.allocate`.

- `CPU()` → `KernelAbstractions.CPU()` (always available).
- `GPU()` → `CUDA.CUDABackend()`, but **only** when the `vSmartMOMCUDAExt`
  extension is loaded (i.e. `using CUDA`). Without it, the default method below
  throws an informative error rather than failing cryptically inside a kernel
  launch. This mirrors the weak-dependency pattern used for `devi`/`array_type`.

This is the architecture→backend mapping consumed by the GPU Mie pipeline
(`compute_aerosol_optical_properties` on a `GPU()`-architecture `MieModel`).
"""
@inline ka_backend(::CPU) = KernelAbstractions.CPU()

# GPU/Metal ka_backend is provided by the matching backend extension as a MORE
# SPECIFIC method (e.g. `ka_backend(::GPU)` in vSmartMOMCUDAExt). We define only
# a generic fallback on the abstract type here — adding a concrete `::GPU`
# method later is a refinement, not a method overwrite, so this stays
# precompile-safe (the same pattern the codebase uses for `devi`/`array_type`).
# The fallback throws an actionable "load the backend" error.
ka_backend(arch::AbstractArchitecture) = error(
    "Architectures.ka_backend($(typeof(arch))) is unavailable: load the " *
    "matching GPU backend (`using CUDA` for GPU, `using Metal` for MetalGPU) " *
    "to enable its KernelAbstractions backend.")

"""
    default_mie_precision_policy(arch::AbstractArchitecture, ::Type{FT})
    default_mie_precision_policy(arch::AbstractArchitecture)            # forwarding shim → Float64

Return the default GPU Mie `Dₙ`-recursion precision policy for `arch` and output
float type `FT`, used when a `MieModel`'s `precision_policy` is `nothing` (auto).
The concrete policy types (`NativeFloat64`, `DSEmulated`) live in the Scattering
module, so the returned value is supplied by the loaded GPU backend extension:

- `CPU()` → `nothing` (the CPU path ignores the precision policy entirely).
- `GPU()` (set by `vSmartMOMCUDAExt`) is **FT-aware**:
    * `Float64` → `Scattering.NativeFloat64()`, the full-FP64 `Dₙ` recursion
      appropriate for A100/V100-class hardware.
    * `Float32` → `Scattering.DSEmulated()`, the Float32-native path. `Float32`
      models cannot use `NativeFloat64` (its kernel asserts `FT === Float64`), so
      auto-selecting it would break the several existing GPU YAMLs built with
      `float_type: Float32`. `DSEmulated` runs the device `Dₙ` recursion in
      Float32 double-single pairs and the host-side bulk reduction is already
      Float64-widened + Neumaier-compensated, so accuracy lands at the Float32
      representational floor — the best a Float32 model can do.

Calling this on a `GPU()` without CUDA.jl loaded throws, matching `ka_backend`.
The 1-arg method is a back-compat forwarding shim that defaults `FT = Float64`.
"""
@inline default_mie_precision_policy(::CPU, ::Type) = nothing
@inline default_mie_precision_policy(::CPU)          = nothing

# Concrete GPU/Metal policies are supplied by the backend extension as more
# specific methods (precompile-safe, same rationale as `ka_backend` above).
default_mie_precision_policy(arch::AbstractArchitecture, ::Type) = error(
    "Architectures.default_mie_precision_policy($(typeof(arch)), FT) is " *
    "unavailable: load the matching GPU backend (`using CUDA` / `using Metal`) " *
    "to enable its Mie precision policy.")

# 1-arg forwarding shim: anything still calling the old single-arg form gets the
# Float64 default policy (the historical behavior).
@inline default_mie_precision_policy(arch::AbstractArchitecture) =
    default_mie_precision_policy(arch, Float64)

"""
    has_gpu_mie(arch::AbstractArchitecture) -> Bool

Capability trait: does `arch` have a GPU Mie pipeline (`ka_backend` +
`default_mie_precision_policy` + the KernelAbstractions Mie kernels)?

- Defaults to `false` for **every** architecture, including `GPU()` and
  `MetalGPU()`. A backend extension opts in by adding a more-specific method
  once it actually wires up the GPU Mie path. `vSmartMOMCUDAExt` sets
  `has_gpu_mie(::GPU) = true`; the Metal extension defines no Mie path, so
  `MetalGPU()` keeps the default `false` and Mie is computed on the CPU (RT
  arrays still run on Metal, exactly as before this trait existed).

The default lives on the abstract supertype and the CUDA ext adds the concrete
`::GPU` refinement — the same precompile-safe pattern as `ka_backend` (a more
specific method is a refinement, not a method overwrite).

Consumed by the aerosol-optics router: `true` → GPU Mie pipeline; `false` →
CPU Mie compute (with a one-time `@warn` when `arch` is not `CPU()`).
"""
@inline has_gpu_mie(::AbstractArchitecture) = false
@inline has_gpu_mie(::CPU) = false

"""
    array_type(arch::AbstractArchitecture)

Return the array constructor for the given architecture (`Array` for CPU,
`CuArray` for CUDA `GPU`, and `MtlArray` for `MetalGPU` when the matching
optional backend extension is loaded).
"""
@inline array_type(::CPU) = Array

# GPU array_type - will be defined in CUDAExt/MetalExt when loaded
# No fallback method defined to avoid precompilation errors

"""
    default_architecture()

Return `GPU()` if CUDA is available, `MetalGPU()` if Metal is available, and
otherwise `CPU()`.
"""
default_architecture() = has_cuda() ? GPU() : has_metal() ? MetalGPU() : CPU()

# Synchronization - calls the sync function set by the loaded GPU backend.
@inline synchronize_if_gpu() = sync_device()

end
