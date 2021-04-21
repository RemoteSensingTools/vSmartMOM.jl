module vSmartMOM

using LinearAlgebra                # For linear algebra routines
using NCDatasets                   # For reading netcdf atmospheric profiles
using ProgressMeter                # Showing progress in for loops
using Distributions                # Distributions of aerosols
using Interpolations               # Interpolations? <<Christian>>
using Polynomials                  # Polynomials? <<Christian>>
using DelimitedFiles               # For reading ASCII files 
using ..Scattering                 # Use scattering module
using ..Absorption                 # Use absorption module
using ...RadiativeTransfer         # Use parent RadiativeTransfer module
using ...Architectures             # Use Architectures module


using CUDA                         # GPU CuArrays and functions
using KernelAbstractions           # Abstracting code for CPU/GPU
using KernelAbstractions.Extras

using FastGaussQuadrature          # FastGaussQuadrature? <<Christian>>
using TimerOutputs                 # For timing sections of the code
using StatsBase                    # StatsBase? <<Christian>>
using StaticArrays                 # StaticArrays? <<Christian>>
using TensorOperations             # TensorOperations? <<Christian>>
using Parameters                   # For keyword arguments in structs
using DocStringExtensions          # For documenting
using YAML                         # For reading properties files 
using ForwardDiff                  # Automatic Differentiation
using NNlib                        # For batched multiplications
import NNlib.batched_mul           # 

# More threads in LA wasn't really helpful, turn that off now and use Julia threads!
LinearAlgebra.BLAS.set_num_threads(1)

# Constants and Types
include("Constants/constants.jl")         # Scientific constants
include("Types/types.jl")                 # All custom types for this module

# Solvers
include("Solvers/elemental.jl")           # Elemental 
include("Solvers/doubling.jl")            # Doubling
include("Solvers/interaction.jl")         # Interaction
include("Solvers/rt_kernel.jl")           # Handle Core RT (Elemental/Doubling/Interaction)
include("Solvers/postprocessing_vza.jl")  # Postprocess (Azimuthal Weighting)
include("Solvers/rt_run.jl")              # Starting point for RT 

# GPU
include("GPU/gpu_batched.jl")             # Batched operations
include("GPU/CUDA_getri.jl")              # Custom getri_strided_batched!

# Utilities / Helper Functions
include("Utils/atmo_prof.jl")             # Helper Functions for Hanling Atmospheric Profiles
include("Utils/rt_utils.jl")              # Miscellaneous Utility Functions
include("Utils/rt_streams.jl")            # Set streams before RT
include("Utils/forwardDiff_tools.jl")     # Helpers for Forward Differentiation
include("Utils/model_parameters.jl")      # Handling Model Parameters 

# Surfaces
include("Surface/lambertian_surface.jl")  # Lambertian Surface 

# Functions to export
export parameters_from_yaml,              # Getting parameters from a file
       model_from_parameters,             # Converting the parameters to model 
       rt_run,                            # Run the RT code
       default_parameters                 # Set of default parameters

end