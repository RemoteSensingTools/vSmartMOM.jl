module vSmartMOM

using LinearAlgebra
using NCDatasets
using ProgressMeter
using Distributions
using Interpolations
using Polynomials
using DelimitedFiles
using ..Scattering
using ..Absorption
using ...RadiativeTransfer
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


include("Types/types.jl") 
include("Utils/rt_streams.jl")
include("Constants/constants.jl")
include("Utils/atmo_prof.jl")
include("Solvers/elemental.jl")
include("Solvers/doubling.jl")
include("Solvers/interaction.jl")
include("GPU/gpu_batched.jl")
include("Utils/rt_utils.jl")
include("Solvers/rt_run.jl")
include("GPU/CUDA_getri.jl")

export rt_set_streams, compute_absorption_profile!, ObsGeometry

end