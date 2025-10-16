#=
 
This file is the entry-point for the vSmartMOM module. 

It includes this module's source files and exports the relevant keywords.  
 
=#

module CoreRT

using Interpolations               # For interpolating the vmr's
using LinearAlgebra                # For linear algebra routines
using ProgressMeter                # Showing progress in for loops
using Distributions                # Distributions of aerosols
using Parameters                   # For keyword arguments in structs
using ..Scattering                 # Use scattering module
using ..Absorption                 # Use absorption module
using ..InelasticScattering        # Use Inelastic Scattering module
using ...vSmartMOM                 # Use parent RadiativeTransfer module
using ...Architectures             # Use Architectures module


using CUDA                         # GPU CuArrays and functions
using KernelAbstractions           # Abstracting code for CPU/GPU
using KernelAbstractions.Extras
using CUDA.CUDAKernels

using Unitful                      # For parsing 
using UnitfulEquivalences          # For converting between wavenumber / wavelength
using FastGaussQuadrature          # Computes quadrature points (Gauss, legendre, Radau,...)
using TimerOutputs                 # For timing sections of the code
using DocStringExtensions          # For documenting
using YAML                         # For reading properties files 
using ForwardDiff                  # Automatic Differentiation
using NNlib                        # For batched multiplications
import NNlib.batched_mul           # Required to overwrite batched_mul for Duals
using NCDatasets                   # For loading absco lookup tables

using JLD2
import Base.show                   # For overloading show for custom types

#using InelasticScattering

# More threads in LA wasn't really helpful, can be turned off here:
# LinearAlgebra.BLAS.set_num_threads(1)

# Constants and Types
include("constants.jl")                        # Scientific constants
include("types.jl")  
include("types_lin.jl")                    # Types for linearized RT


# All custom types for this module
# Raman additions
#include("Inelastic/types.jl")
#include("Inelastic/inelastic_helper.jl")
#include("Inelastic/raman_atmo_prop.jl")

# Solvers
include("CoreKernel/elemental.jl")             # Elemental 
include("CoreKernel/elemental_lin.jl") 
include("CoreKernel/elemental_inelastic.jl")   # Elemental for inelastic scattering
include("CoreKernel/elemental_inelastic_plus.jl")   # Elemental for inelastic scattering
include("CoreKernel/doubling.jl")              # Doubling
include("CoreKernel/doubling_lin.jl") 
include("CoreKernel/doubling_inelastic.jl")    # Doubling for elastic + inelastic scattering 
include("CoreKernel/interaction.jl")           # Interaction
include("CoreKernel/interaction_lin.jl") 
include("CoreKernel/interaction_inelastic.jl") # Interaction for elastic + inelastic scattering 
include("CoreKernel/interaction_multisensor.jl") # Suniti: ms
include("CoreKernel/interaction_ss.jl") #for single scattering contribution only
include("CoreKernel/interlayer_flux.jl")       # Suniti: ms
include("CoreKernel/rt_kernel.jl")             # Handle Core RT (Elemental/Doubling/Interaction)
include("CoreKernel/rt_kernel_lin.jl")
include("CoreKernel/rt_kernel_ss.jl")          # Single scattering only: Handle Core RT (Elemental/Doubling/Interaction)
include("CoreKernel/rt_kernel_multisensor.jl") # Suniti: ms
include("CoreKernel/lin_added_layer_all_params.jl") # 3 parameters -> all parameters
include("postprocessing_vza.jl")               # Postprocess (Azimuthal Weighting)
include("postprocessing_vza_lin.jl") 
include("postprocessing_vza_ms.jl")
include("rt_run.jl")                           # Starting point for RT 
include("rt_run_lin.jl")                        
include("rt_run_multisensor.jl")  
# GPU
include("gpu_batched.jl")                   # Batched operations

# Utilities / Helper Functions
include("atmo_prof.jl")                     # Helper Functions for Handling Atmospheric Profiles
include("atmo_prof_lin.jl")
include("rt_helper_functions.jl")           # Miscellaneous Utility Functions
include("rt_helper_functions_lin.jl") 
include("rt_set_streams.jl")                # Set streams before RT
include("parameters_from_yaml.jl")          # Loading in parameters from YAML file
include("model_from_parameters.jl")         # Converting parameters to derived model attributes
include("lin_model_from_parameters.jl")
include("show_utils.jl")                    # Pretty-printing objects
include("LayerOpticalProperties/compEffectiveLayerProperties.jl")
include("LayerOpticalProperties/compEffectiveLayerProperties_lin.jl")

# Surfaces
include("lambertian_surface.jl")            # Lambertian Surface 
include("lambertian_surface_lin.jl") 


# Functions to export
export parameters_from_yaml,                # Getting parameters from a file
       model_from_parameters,               # Converting the parameters to model 
       rt_run, rt_run_ss,                             # Run the RT code
       default_parameters                   # Set of default parameters
export lin_added_layer_all_params  # 3 parameters -> all parameters 
                   # vSmartMOM Modes
# Export types to show easily
export GaussQuadFullSphere, LambertianSurfaceScalar, LambertianSurfaceSpectrum

end