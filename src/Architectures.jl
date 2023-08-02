module Architectures
# Shout out to https://github.com/CliMA/Oceananigans.jl/blob/master/src/Architectures.jl

export
    @hascuda,
    AbstractArchitecture, CPU, GPU,
    devi, array_type, 
    architecture,
    default_architecture,
    synchronize_if_gpu

using CUDA

using KernelAbstractions
using CUDAKernels

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

devi(::CPU) = KernelAbstractions.CPU()
devi(::GPU) = CUDAKernels.CUDADevice()

         architecture(::Array)   = CPU()
#@hascuda 
        architecture(::CuArray) = GPU()

         array_type(::CPU) = Array
#@hascuda 
        array_type(::GPU) = CuArray

default_architecture = has_cuda() ? GPU() : CPU()

synchronize_if_gpu() = has_cuda() ? synchronize() : nothing

end