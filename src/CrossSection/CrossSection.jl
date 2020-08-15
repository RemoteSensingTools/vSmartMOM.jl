module CrossSection

using Parameters                # For constructing HitranTable with keywords
using DocStringExtensions       # For simplifying docstring
using Interpolations            # For interpolating in both lookup tables and qoft!
using JLD2                      # For saving and loading the interpolator
using ProgressMeter             # For showing progress, especially in creating the interpolator
using KernelAbstractions        # For heterogeneous (GPU+CPU) programming
using CUDA                      # For GPU programming
using ..Architectures: device

include("constants.jl")         # Scientific and mathematical constants
include("types.jl")             # All types used in this module
include("hitran.jl")            # HITRAN file-related functions
include("complex_error_functions.jl")   # Complex error functions used in line broadening
include("mol_weights.jl")       # Molecular weights (TODO: replace with netCDF file)
include("TIPS_2017.jl")         # Partition sums data (TODO: replace with netCDF file)
include("partition_sums.jl")    # Partition sums interpolator (TODO: replace with LinearInterpolation)
include("cross_section_interpolator.jl") # Cross-section interpolator functions
include("absorption_cross_section.jl")   # Cross-section from HITRAN

# Export the Cross Section models
export HitranModel, InterpolationModel

# Export the absorption_cross_section functions
export compute_absorption_cross_section, absorption_cross_section

# Export the hitran functions from hitran.jl
export read_hitran, make_hitran_model

# Export the broadening function types
export Doppler, Lorentz, Voigt

# Export the complex error functions
export HumlicekErrorFunction, HumlicekWeidemann32VoigtErrorFunction, HumlicekWeidemann32SDErrorFunction, CPF12ErrorFunction, ErfcHumliErrorFunctionVoigt, ErfcHumliErrorFunctionSD, ErfcErrorFunction

# Export the hitran table struct type
export HitranTable

# Export the interpolator functions
export make_interpolation_model, save_interpolation_model, load_interpolation_model



end
