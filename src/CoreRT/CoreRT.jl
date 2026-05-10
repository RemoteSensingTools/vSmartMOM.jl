#=
 
This file is the entry-point for the vSmartMOM module. 

It includes this module's source files and exports the relevant keywords.  
 
=#

"""
    CoreRT

Core radiative-transfer implementation: model types, layer optical-property
construction, adding-doubling kernels, surface coupling, and solver entry
points.
"""
module CoreRT

using Interpolations               # For interpolating the vmr's
using LinearAlgebra                # For linear algebra routines
using ProgressMeter                # Showing progress in for loops
using Distributions                # Distributions of aerosols
using Parameters                   # For keyword arguments in structs
using Unitful                      # For unit-aware conversions
using UnitfulEquivalences          # For spectral equivalences (nm ↔ cm⁻¹)
using ..Scattering                 # Use scattering module
using ..Absorption                 # Use absorption module
using ..InelasticScattering        # Use Inelastic Scattering module
using ...vSmartMOM                 # Use parent RadiativeTransfer module
using ...Architectures             # Use Architectures module

using KernelAbstractions           # Abstracting code for CPU/GPU
using KernelAbstractions.Extras

using FastGaussQuadrature          # Computes quadrature points (Gauss, legendre, Radau,...)
using TimerOutputs                 # For timing sections of the code
using DocStringExtensions          # For documenting
using YAML                         # For reading properties files 
using ForwardDiff                  # Automatic Differentiation
using NNlib                        # For batched multiplications
import NNlib.batched_mul           # Required to overwrite batched_mul for Duals
using NCDatasets                   # For loading absco lookup tables
using QuadGK
using CanopyOptics
using SpecialFunctions: erfc         # For Smith (1967) shadowing in Cox-Munk
using StaticArrays                   # For Fresnel/Cox-Munk Mueller matrices

import Base.show                   # For overloading show for custom types

#using InelasticScattering

# More threads in LA wasn't really helpful, can be turned off here:
# LinearAlgebra.BLAS.set_num_threads(1)

# Constants and Types
include("constants.jl")                        # Scientific constants
include("Sources/types.jl")                    # v0.6 AbstractSource / SourceSet vocabulary (Phase 1) — must precede types.jl since RTModel.sources::AbstractSource
include("types.jl")
include("parameter_layout.jl")                 # ParameterLayout struct for Jacobian indexing
include("types_lin.jl")                        # Types for linearized RT
include("Sources/solar_beam.jl")               # v0.6 SolarBeam + PreparedSolarBeam (Phase 2)
include("Sources/surface_sif.jl")              # v0.6 SurfaceSIF + PreparedSurfaceSIF (Phase 5)

# Note: Raman types (AbstractRamanType, noRS, RRS, etc.) come from 
# InelasticScattering module via `using ..InelasticScattering` above.
# Inelastic layer types (CompositeLayerRS, AddedLayerRS) are in types.jl.

# Shared helpers
include("CoreKernel/rt_helpers.jl")

# Solvers -- Elemental
include("CoreKernel/elemental.jl")             # Elemental (elastic)
include("CoreKernel/elemental_lin.jl")         # Elemental (linearized)
include("CoreKernel/elemental_inelastic.jl")   # Elemental for inelastic scattering
include("CoreKernel/elemental_inelastic_plus.jl")   # Elemental for inelastic scattering (VS/RVRS)
include("CoreKernel/elemental_canopy.jl")

# Solvers -- Doubling
include("CoreKernel/doubling.jl")              # Doubling (elastic)
include("CoreKernel/doubling_lin.jl")          # Doubling (linearized)
include("CoreKernel/doubling_inelastic.jl")    # Doubling for elastic + inelastic scattering 

# Solvers -- Interaction
include("CoreKernel/interaction.jl")           # Interaction (elastic)
include("CoreKernel/interaction_lin.jl")       # Interaction (linearized)
include("CoreKernel/interaction_hdrf.jl")      # Addl. surface interaction for RAMI output
include("CoreKernel/interaction_inelastic.jl") # Interaction for elastic + inelastic scattering 
include("CoreKernel/interaction_multisensor.jl") # Multi-sensor interaction
include("CoreKernel/interaction_ss.jl")        # Single scattering interaction
include("CoreKernel/interlayer_flux.jl")       # Interlayer flux

# Solvers -- RT Kernels
include("CoreKernel/rt_kernel.jl")             # Handle Core RT (Elemental/Doubling/Interaction)
include("CoreKernel/rt_kernel_lin.jl")         # Linearized RT kernel
include("CoreKernel/rt_kernel_ss.jl")          # Single scattering RT kernel
include("CoreKernel/rt_kernel_multisensor.jl") # Multi-sensor RT kernel
include("CoreKernel/lin_added_layer_all_params.jl") # 3 core params -> all state params

# Postprocessing
include("tools/postprocessing_vza.jl")               # Postprocess (Azimuthal Weighting)
include("tools/postprocessing_vza_lin.jl")           # Postprocess linearized
include("tools/postprocessing_vza_ms.jl")

# RT Run entry points
include("rt_run.jl")                           # Starting point for RT 
include("rt_run_lin.jl")                       # Linearized RT run
include("rt_run_multisensor.jl")

# CPU batched operations (always available)
include("tools/ka_batched_kernels.jl")         # Portable KA batched kernels
include("tools/cpu_batched.jl")                   # CPU batched linear algebra operations

# Utilities / Helper Functions
include("tools/atmo_prof.jl")                     # Helper Functions for Handling Atmospheric Profiles
include("tools/atmo_prof_lin.jl")                 # Linearized atmosphere profile
include("tools/rt_helper_functions.jl")           # Miscellaneous Utility Functions
include("tools/rt_helper_functions_lin.jl")       # Linearized helper functions
include("tools/rt_set_streams.jl")                # Set streams before RT
include("tools/model_from_parameters.jl")         # Converting parameters to derived model attributes
include("tools/lin_model_from_parameters.jl")     # Linearized model from parameters
include("tools/show_utils.jl")                    # Pretty-printing objects
include("LayerOpticalProperties/compEffectiveLayerProperties.jl")
include("LayerOpticalProperties/delta_m_truncation.jl")         # δ-M truncation + chain rule
include("LayerOpticalProperties/compEffectiveLayerProperties_lin.jl")

# Surfaces
include("Surfaces/lambertian_surface.jl")            # Lambertian Surface 
include("Surfaces/lambertian_surface_lin.jl")        # Linearized Lambertian Surface
include("Surfaces/rpv_surface.jl")                   # RPV Surface 
include("Surfaces/rossli_surface.jl")                # Ross-Li Surface
include("Surfaces/canopy_surface.jl")                # Canopy + soil composite surface
include("Surfaces/fresnel.jl")                       # Fresnel reflection utilities
include("Surfaces/water_refraction.jl")              # Built-in water refractive index
include("Surfaces/coxmunk_surface.jl")               # Cox-Munk ocean surface
include("Surfaces/coxmunk_surface_lin.jl")           # Linearized Cox-Munk (Jacobians)

# Phase C of the v0.7 Fourier/Stream Resolution refactor — per-component
# trait dispatch. Comes after all surface/source/scatterer types are
# defined (types.jl + Sources/ + Surfaces/) but before
# tools/model_from_parameters.jl which consumes the aggregator.
include("component_m_max.jl")

# Functions to export
export model_from_parameters,               # Converting the parameters to model
       model_from_parameters_lin,           # Convenience alias for linearized model
       rt_run, rt_run_lin, rt_run_ss,       # Run the RT code (forward, linearized, single scatter)
       default_parameters                   # Set of default parameters
export lin_added_layer_all_params!,           # 3 params -> all params chain rule
       OpticalPropertyJacobian,               # AD boundary struct alias
       RawAerosolJacobian,                    # AD boundary for raw aerosol derivatives
       delta_m_forward, delta_m_truncation_lin, # δ-M truncation functions
       ParameterLayout,                       # Jacobian parameter layout descriptor
       n_total, aerosol_range, gas_range,     # ParameterLayout accessors
       surface_range, surface_index,
       canopy_range, n_layer_params

# Export new hierarchical model types
export AbstractRTModel, RTModel,
       SolverConfig, RTNumericalParameters,
       Atmosphere, RayleighScattering, AerosolState, Optics, OpticsLin,
       RTModelLin

# v0.6 source-term vocabulary (Phase 1: types only; behaviour lands in Phase 2+)
export AbstractSource, AbstractPreparedSource,
       NoSource, SourceSet,
       AbstractSourceADMode, AnalyticSourceJacobian, ForwardDiffSourceJacobian,
       NoSourceJacobian, source_ad_mode

# Phase 2: SolarBeam + the prepare_source seam
# Phase 4: BlackbodySource constructor sugar (returns a SolarBeam with Planck F₀)
# Phase 5: SurfaceSIF + double-dispatch surface_source_contribute!
export SolarBeam, PreparedSolarBeam, prepare_source, prepare_sources, BlackbodySource,
       SurfaceSIF, PreparedSurfaceSIF, surface_source_contribute!

# Export types to show easily
export RadauQuad, GaussLegQuad,
       LambertianSurfaceScalar, LambertianSurfaceSpectrum,
       CanopySurface, CanopySurface_from_prospect, invalidate_canopy_cache!,
       CoxMunkSurface, water_refractive_index

end
