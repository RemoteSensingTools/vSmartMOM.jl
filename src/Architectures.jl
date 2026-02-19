module Architectures
# Shout out to https://github.com/CliMA/Oceananigans.jl/blob/master/src/Architectures.jl

export
    @hascuda,
    AbstractArchitecture, CPU, GPU,
    devi, array_type, 
    architecture,
    default_architecture,
    synchronize_if_gpu,
    has_cuda

using KernelAbstractions

# Check if CUDA is available and functional
# This will be true when the CUDAExt extension is loaded
const _has_cuda = Ref(false)
has_cuda() = _has_cuda[]

# Synchronization function reference - will be set by CUDAExt
const _sync_gpu = Ref{Function}(() -> nothing)
sync_device() = _sync_gpu[]()

"""
    AbstractArchitecture
Abstract supertype for architectures supported by Oceananigans.
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
    array_type(arch::AbstractArchitecture)

Return the array constructor for the given architecture (`Array` for CPU, `CuArray` for GPU).
"""
@inline array_type(::CPU) = Array

# GPU array_type - will be defined in CUDAExt when CUDA is loaded
# No fallback method defined to avoid precompilation errors

"""
    default_architecture()

Return `GPU()` if CUDA is available, otherwise `CPU()`.
"""
default_architecture() = has_cuda() ? GPU() : CPU()

# Synchronization - calls the sync function (set by CUDAExt if CUDA is loaded)
@inline synchronize_if_gpu() = sync_device()

end