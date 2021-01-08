module RTM

using LinearAlgebra
using NCDatasets
using ProgressMeter
using Distributions
using Interpolations
using Polynomials
using DelimitedFiles
using ..PhaseFunction
using FastGaussQuadrature 
using CUDA
using TimerOutputs
using StatsBase
using KernelAbstractions
using KernelAbstractions.Extras
using StaticArrays
using TensorOperations
using NNlib
using Parameters


include("types.jl") 
include("rt_streams.jl")
include("constants.jl")
include("atmo_prof.jl")
include("rt_elemental.jl")
include("rt_doubling.jl")
include("rt_interaction.jl")
include("gpu_batched.jl")
include("rt_tools.jl")

export rt_set_streams

end