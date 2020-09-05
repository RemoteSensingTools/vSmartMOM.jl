module PhaseFunction

using Parameters                # 
using DocStringExtensions       # For simplifying docstring
using FastGaussQuadrature       # For fast Gauss-Legendre quadrature points
using JLD2                      # For saving and loading 
using ProgressMeter             # For showing progress
using KernelAbstractions        # For heterogeneous (GPU+CPU) programming
using CUDA                      # For GPU programming
using Distributions             # Distributions from Julia 
using ForwardDiff, DiffResults  # Automatic Differentiation tools
using ..Architectures: device

include("types.jl")              # All types used in this module
include("mie_functions.jl")      # Mie file-related functions
include("mie_bulk_methods.jl")   # Functions for Mie calculations over size distribution
include("wigner3j_recursive.jl") # Recursive Wigner 3j calculations

# Export the mie models
export compute_ab, UnivariateAerosol, comp_ab, wigner!, compute_mie_π_τ!

end
