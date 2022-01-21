#=
 
This file is the entry-point for the Scattering module. 

It includes this module's source files and exports the relevant keywords.  
 
=#

module Scattering

using Parameters                # For default values in structs
using DocStringExtensions       # For simplifying docstring
using FastGaussQuadrature       # For fast Gauss-Legendre quadrature points
using JLD2                      # For saving and loading 
using ProgressMeter             # For showing progress
using KernelAbstractions        # For heterogeneous (GPU+CPU) programming
using CUDA                      # For GPU programming
using Distributions             # Distributions from Julia 
using ForwardDiff, DiffResults  # Automatic Differentiation tools
using LinearAlgebra             # For calculations
using StatsBase                 # Fit statistics for truncation
using StaticArrays              # 
# using UnicodePlots              # For plotting aerosol distribution

using ..Architectures: device

include("types.jl")                       # All types used in this module
include("mie_helper_functions.jl")        # Mie file-related functions
include("legendre_functions.jl")          # Recursions for associated Legendre Polynomials
include("truncate_phase.jl")              # Functions for truncation
include("compute_wigner_values.jl")       # Recursive Wigner 3j calculations
include("make_mie_model.jl")              # Convenience functions to create Mie model
include("compute_NAI2.jl")                # Compute phase function w/ NAI2
include("compute_PCW.jl")                 # Compute phase function w/ PCW
include("phase_function_autodiff.jl")     # Auto-differentiation
include("show_utils.jl")                  # Pretty-print


# Export make functions/types
export make_mie_model, reconstruct_phase

# Export types
export NAI2, PCW, Aerosol, MieModel, Stokes_IQUV, Stokes_I, Stokes_IQU, 
       δBGE, GreekCoefs, AerosolOptics, AbstractFourierDecompositionType

export compute_B, compute_ab, GreekCoefs, comp_ab, compute_mie_π_τ!, 
       compute_wigner_values, save_wigner_values, load_wigner_values, 
       compute_Sl, gausslegendre, compute_aerosol_optical_properties,
       compute_ref_aerosol_extinction,
       ConjugateTransposePairs, AbstractPolarizationType, 
       AbstractAerosolType, AbstractAerosolType, MieModel, 
       AbstractTruncationType, phase_function

end
