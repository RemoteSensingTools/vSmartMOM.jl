#=
 
This file is the entry-point for the Absorption module. 

It includes this module's source files and exports the relevant keywords.  
 
=#

module Absorption

using SpecialFunctions          # For complex error functions
using Parameters                # For constructing HitranTable with keywords
using DocStringExtensions       # For simplifying docstring
using Interpolations            # For interpolating in lookup tables and interpolating CS
using JLD2                      # For saving and loading the interpolator
using ProgressMeter             # For showing progress, especially in creating interpolator
using KernelAbstractions        # For heterogeneous (GPU+CPU) programming
using CUDA.CUDAKernels               # Access to CUDADevice
using CUDA                      # For GPU programming
using ForwardDiff, DiffResults  # For auto-differentiation
using NetCDF                    # For loading NetCDF files with constants
using ..Architectures           # For GPU/CPU convenience
using ..Architectures: CPU, GPU # Again for GPU/CPU convenience
import DataInterpolations: CubicSpline as DI_CS # For use in qoft

include("constants/constants.jl")                   # Scientific and mathematical constants
include("constants/mol_weights.jl")                 # Molecular weights
include("constants/TIPS_2017.jl")                   # Partition sums data 
include("types.jl")                                 # All types used in this module
include("read_hitran.jl")                                # HITRAN file-related functions
include("complex_error_functions.jl")               # CEFs used in line broadening
include("make_model_helpers.jl")                  # CS interpolator functions
include("compute_absorption_cross_section.jl")      # Cross-section from HITRAN
include("autodiff_helper.jl")                       # Auto-differentiation

# Export the Cross Section models
export AbstractCrossSectionModel, HitranModel, InterpolationModel

# Export the absorption_cross_section functions
export compute_absorption_cross_section, absorption_cross_section

# Export the hitran functions from hitran.jl
export read_hitran, make_hitran_model

# Export the broadening function types
export Doppler, Lorentz, Voigt, AbstractBroadeningFunction

# Export the complex error functions
export AbstractComplexErrorFunction, HumlicekErrorFunction, 
       HumlicekWeidemann32VoigtErrorFunction, HumlicekWeidemann32SDErrorFunction, 
       CPF12ErrorFunction, ErfcHumliErrorFunctionVoigt, ErfcHumliErrorFunctionSD, 
       ErfcErrorFunction

# Export the hitran table struct type
export HitranTable

# Export the interpolator functions
export make_interpolation_model, save_interpolation_model, load_interpolation_model

end
