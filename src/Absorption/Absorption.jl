module Absorption

using Parameters                # For constructing HitranTable with keywords
using DocStringExtensions       # For simplifying docstring
using Interpolations            # For interpolating in both lookup tables and qoft!
using JLD2                      # For saving and loading the interpolator
using ProgressMeter             # For showing progress, especially in creating interpolator
using KernelAbstractions        # For heterogeneous (GPU+CPU) programming
using CUDA                      # For GPU programming
using ForwardDiff, DiffResults  # For auto-differentiation
using ..Architectures           # For GPU/CPU convenience
using NetCDF                    # For loading NetCDF files with constants

include("Constants/constants.jl")                   # Scientific and mathematical constants
include("Types/types.jl")                           # All types used in this module
include("Hitran/hitran.jl")                         # HITRAN file-related functions
include("LineShapes/complex_error_functions.jl")    # CEFs used in line broadening
include("Constants/mol_weights.jl")                 # Molecular weights
include("Constants/TIPS_2017.jl")                   # Partition sums data 
include("TableInterpolation/partition_sums.jl")     # Partition sums interpolator 
include("TableInterpolation/cross_section_interpolator.jl") # CS interpolator functions
include("Hitran/absorption_cross_section.jl")       # Cross-section from HITRAN
include("Hitran/cross_section_autodiff.jl")         # Auto-differentiation

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
