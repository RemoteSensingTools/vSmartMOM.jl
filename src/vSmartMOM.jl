#=
 
This file is the entry-point for the vSmartMOM module. 

It includes this module's source files and exports the relevant keywords.  
 
=#

"""
    vSmartMOM

Top-level package module for configuring and running vector adding-doubling
radiative-transfer simulations.
"""
module vSmartMOM
using Pkg.Artifacts
using LinearAlgebra
using Distributions
using Parameters
using DocStringExtensions
using UnPack
using UnicodePlots
using Dates
using Downloads: Downloads
using SHA: sha256
using Scratch: @get_scratch!

# Export Architecture functions
export CPU, GPU, MetalGPU, default_architecture, array_type

# Export the artifact convenience function and HITRAN management
export artifact
export fetch_hitran, fetch_hitran_by_ids
export set_hitran_edition!, get_hitran_edition, available_hitran_editions, hitran_info
export hitran_is_cached

# RT mode types (forward vs linearized)
abstract type RT_Mode end

"""
    FwdMode()

Marker selecting the forward radiative-transfer model construction path.
"""
struct FwdMode <: RT_Mode end

"""
    LinMode()

Marker selecting the linearized model construction path used for Jacobians.
"""
struct LinMode <: RT_Mode end
export FwdMode, LinMode

# GPU/CPU Architecture (from Oceanigans)
include("Architectures.jl")
using .Architectures

# Absorption Cross Section module (loaded before Artifacts so mol_names/global_ids are available):
include("Absorption/Absorption.jl")

# Artifacts — HITRAN data management (preferences, API client, artifact dispatch)
include("Artifacts/hitran_preferences.jl")
include("Artifacts/hitran_api.jl")
include("Artifacts/artifact_helper.jl")

# Mie Phase Function module:
include("Scattering/Scattering.jl")

# Inelastic Scattering module:
include("Inelastic/InelasticScattering.jl")

# CoreRT module:
include("CoreRT/CoreRT.jl")
using .CoreRT

# Standalone exact single-scattering module:
include("StandaloneSS/StandaloneSS.jl")
using .StandaloneSS

# SolarModel module:
include("SolarModel/SolarModel.jl")

# IO submodule (must come after CoreRT types are defined)
include("IO/IO.jl")
import .IO: parameters_from_file, parameters_from_source,
            parameters_from_yaml, parameters_from_dict, read_parameters,
            read_atmos_profile, read_atmos_profile_dict,
            GeosChemSource, NetCDFGridSource, NetCDFSource,
            geoschem_to_dict, read_geoschem_profile

# Aerosols module — user-facing flexible aerosol framework (TOMAS-15 +
# two-moment schemes). WIP header in Aerosols.jl documents follow-up cleanup.
include("Aerosols/Aerosols.jl")
using .Aerosols

# SIF emission data + loaders
import DataInterpolations: ExtrapolationType
include("SIF_emission/sif_loader.jl")
export load_sif_spectrum, load_ficus_reflectance, sif_data_path, build_sif_source

# Export some vSmartMOM functions
export default_parameters, parameters_from_file, parameters_from_source,
       parameters_from_yaml, parameters_from_dict,
       model_from_parameters, rt_run, read_parameters,
       read_atmos_profile, read_atmos_profile_dict
# Export linearized RT functions
export rt_run_lin, model_from_parameters_lin
# Export standalone exact single-scattering API
export StandaloneSS, run_exact_ss, ExactSSConfig, SSGeometry,
       LambertianSSSurface, RayleighSSContributor, HGAerosolSSContributor,
       AbsorptionSSContributor
# Export new hierarchical model types
export RTModel, AbstractRTModel, SolverConfig, Atmosphere, RayleighScattering, AerosolState, Optics, OpticsLin
# Export GEOSChem/NetCDF integration
export GeosChemSource, NetCDFGridSource, NetCDFSource, geoschem_to_dict, read_geoschem_profile

using .Architectures
using .Absorption

# Module initialization - GPU setup happens in optional backend extensions.
function __init__()
    # Nothing to do here; GPU detection is handled by backend extensions.
end

end
