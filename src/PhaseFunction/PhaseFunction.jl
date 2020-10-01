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
using LinearAlgebra
using BenchmarkTools
using StaticArrays
using IndexedTables

using ..Architectures: device

include("types.jl")                # All types used in this module
include("mie_helper_functions.jl") # Mie file-related functions
include("legendre_functions.jl")   # Recursions for associated Legendre Polynomials
include("mie_computations.jl")     # Functions for Mie calculations over size distribution
include("phase_truncation.jl")     # Functions for truncation
include("wigner3j_recursive.jl")   # Recursive Wigner 3j calculations
include("compute_Sl.jl")           # Sl_νν functions
include("aerosol.jl")
include("mie_model.jl")
include("compute_B.jl")


# Export make functions/types
export make_univariate_aerosol, make_mie_model

# Export types
export NAI2, PCW, UnivariateAerosol, MieModel, Stokes_IQUV, Stokes_I, Stokes_IQU, 
       δBGE, GreekCoefs, AerosolOptics

export compute_B, compute_ab, GreekCoefs, comp_ab, compute_mie_π_τ!, 
       compute_wigner_values, save_wigner_values, load_wigner_values, 
       compute_Sl, gausslegendre, compute_aerosol_optical_properties

end
