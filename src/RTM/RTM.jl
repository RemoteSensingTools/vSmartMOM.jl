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

include("types.jl") 
include("constants.jl")
include("RTM_main.jl")
include("atmo_prof.jl")
include("rt_elemental.jl")
include("rt_doubling.jl")
include("rt_interaction.jl")
include("rt_tools.jl")

export rt_set_streams

end