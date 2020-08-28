module PhaseFunction

using Parameters                # For constructing HitranTable with keywords
using DocStringExtensions       # For simplifying docstring
using FastGaussQuadrature       # For fast Gauss-Legendre quadrature points
using JLD2                      # For saving and loading the interpolator
using ProgressMeter             # For showing progress, especially in creating the interpolator
using KernelAbstractions        # For heterogeneous (GPU+CPU) programming
using CUDA                      # For GPU programming
using Distributions             # Distributions from Julia 
using ..Architectures: device

include("types.jl")              # All types used in this module
include("mie_functions.jl")      # Mie file-related functions
include("wigner3j_recursive.jl") # Recursive Wigner 3j calculations

# Export the mie models
export compute_ab, UnivariateAerosol, comp_ab

end
