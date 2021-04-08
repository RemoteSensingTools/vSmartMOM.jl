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
using ...Architectures

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
using DocStringExtensions
using JSON
using CSV
using YAML
using ForwardDiff
import NNlib.batched_mul

# More threads in LA wasn't really helpful, turn that off now and use Julia threads!
LinearAlgebra.BLAS.set_num_threads(1)

include("Types/types.jl") 
include("Utils/rt_streams.jl")
include("Constants/constants.jl")
include("Utils/atmo_prof.jl")
include("Solvers/elemental.jl")
include("Solvers/doubling.jl")
include("Solvers/interaction.jl")
include("GPU/gpu_batched.jl")
include("Utils/rt_utils.jl")
include("Utils/forwardDiff_tools.jl")
include("Solvers/rt_run.jl")
include("GPU/CUDA_getri.jl")
include("Utils/default_model.jl")
include("Surface/lambertian_surface.jl")

export rt_set_streams, compute_absorption_profile!, ObsGeometry, default_model, 
       model_from_parameters, parameters_from_json, parameters_from_yaml

end