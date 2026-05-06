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
using Distributions             # Distributions from Julia 
using ForwardDiff, DiffResults  # Automatic Differentiation tools
using LinearAlgebra             # For calculations
using StatsBase                 # Fit statistics for truncation
using StaticArrays              # 

# Architectures used for CPU/GPU dispatch
using ..Architectures
using ..vSmartMOM: RT_Mode, FwdMode, LinMode

include("types.jl")                       # All types used in this module
include("types_lin.jl")                   # Derivative types for linearization
include("gpu_precision.jl")              # DS arithmetic, Neumaier accumulators (GPU-compatible)
include("mie_helper_functions.jl")        # Mie file-related functions
include("analytic_phase_functions.jl")    # Analytic phase functions -> Greek coefficients
include("mie_helper_functions_lin.jl")    # Linearized Mie functions
include("legendre_functions.jl")          # Recursions for associated Legendre Polynomials
include("truncate_phase.jl")              # Functions for truncation
include("truncate_phase_lin.jl")          # Linearized truncation
include("compute_wigner_values.jl")       # Recursive Wigner 3j calculations
include("make_mie_model.jl")              # Convenience functions to create Mie model
include("compute_NAI2.jl")                # Compute phase function w/ NAI2
include("compute_NAI2_lin.jl")            # Linearized NAI2 computation
include("gpu_mie_kernels.jl")            # KernelAbstractions Mie kernels (CPU+GPU)
include("compute_NAI2_gpu.jl")           # GPU-accelerated NAI2 entry point
include("compute_PCW.jl")                 # Compute phase function w/ PCW
include("phase_function_autodiff.jl")     # Auto-differentiation
include("show_utils.jl")                  # Pretty-print
include("compute_Z_matrices.jl")

# Export make functions/types
export make_mie_model, reconstruct_phase, greek_coefficients,
       greek_coefficients_from_scattering_matrix, analytic_aerosol_optics,
       phase_matrix_first_column, scattering_matrix

# Export types
export NAI2, PCW, Aerosol, MieModel, Stokes_IQUV, Stokes_I, Stokes_IQ,
       Stokes_IQU,
       δBGE, NoTruncation, GreekCoefs, AerosolOptics,
       AbstractFourierDecompositionType, AbstractAnalyticPhaseFunction,
       HenyeyGreensteinPhaseFunction,
       SyntheticPolarizedHenyeyGreensteinPhaseFunction

# Export linearized types
export linGreekCoefs, linAerosolOptics

export GreekCoefs, compute_mie_π_τ,
       compute_wigner_values, save_wigner_values, load_wigner_values,
       compute_Sl, gausslegendre, compute_aerosol_optical_properties,
       compute_ref_aerosol_extinction, truncate_phase,
       AbstractPolarizationType,
       AbstractAerosolType, AbstractAerosolType, MieModel,
       AbstractTruncationType, phase_function, compute_aerosol_XS

# GPU Mie exports
export compute_aerosol_optical_properties_gpu,
       MiePrecisionPolicy, NativeFloat64, DSEmulated,
       DoubleSingle, ComplexDS, NeumaierAccum,
       neumaier_add, neumaier_sum

end
