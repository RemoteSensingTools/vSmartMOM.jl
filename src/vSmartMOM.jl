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

export HenyeyGreensteinPhaseFunction,
       SyntheticPolarizedHenyeyGreensteinPhaseFunction,
       greek_coefficients, analytic_aerosol_optics,
       phase_matrix_first_column

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
using .Scattering: HenyeyGreensteinPhaseFunction,
                   SyntheticPolarizedHenyeyGreensteinPhaseFunction,
                   greek_coefficients, analytic_aerosol_optics,
                   phase_matrix_first_column

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

# Aerosols module — user-facing flexible aerosol framework (TOMAS-15 +
# two-moment schemes). Loaded before IO so IO/NetCDF/GCHPScene.jl can type
# its `aerosols` field against AbstractAerosolBinData. Aerosols itself has
# no IO module dependency (verified: no `using .IO` in src/Aerosols/).
include("Aerosols/Aerosols.jl")
using .Aerosols

# IO submodule (must come after CoreRT types and Aerosols are defined)
include("IO/IO.jl")
import .IO: parameters_from_file, parameters_from_source,
            parameters_from_yaml, parameters_from_dict, read_parameters,
            read_atmos_profile, read_atmos_profile_dict,
            GeosChemSource, NetCDFGridSource, NetCDFSource,
            geoschem_to_dict, read_geoschem_profile,
            GCHPFile, GCHPScene, scene_at, scenes,
            read_gchp_scene, scene_to_dict, parameters_from_scene,
            compute_scene_aod,
            write_scene_result, generate_benchmark,
            write_gchp_aod_diagnostic, write_gchp_aod_diagnostic_bulk

# SIF emission data + loaders
import DataInterpolations: ExtrapolationType
include("SIF_emission/sif_loader.jl")
export load_sif_spectrum, load_ficus_reflectance, sif_data_path, build_sif_source

# Export some vSmartMOM functions
export default_parameters, parameters_from_file, parameters_from_source,
       parameters_from_yaml, parameters_from_dict,
       model_from_parameters, rt_run, read_parameters,
       read_atmos_profile, read_atmos_profile_dict
# GEOS-Chem / GCHP scene IO
export GeosChemSource, geoschem_to_dict, read_geoschem_profile
export GCHPFile, GCHPScene, scene_at, scenes,
       read_gchp_scene, scene_to_dict, parameters_from_scene,
       compute_scene_aod
export write_scene_result, generate_benchmark,
       write_gchp_aod_diagnostic, write_gchp_aod_diagnostic_bulk
# Export linearized RT functions
export rt_run_lin, model_from_parameters_lin
# gchp-io: atmosphere/surface split (Phase C) + scenario sweep (Phase D)
export rt_run_atmosphere, rt_run_surface, rt_run_multi_surface, AtmosphereRTCache
export ScenarioSweep, SweepResult, run_sweep,
       SceneOptics, scene_optics, model_for_sza
# Export standalone exact single-scattering API
export StandaloneSS, run_exact_ss, ExactSSConfig, SSGeometry,
       LambertianSSSurface, CoxMunkSSSurface, RayleighSSContributor,
       SSMeasurementSelector,
       HGAerosolSSContributor, GreekCoefsSSContributor,
       AbsorptionSSContributor,
       determine_required_l_aerosol,
       determine_required_l_from_moments, determine_required_nbrdf,
       determine_required_nbrdf_coxmunk, determine_required_nquad,
       determine_required_nquad_inner, determine_required_nstreams,
       exact_ss_config_from_model,
       run_exact_ss_with_jacobians, chain_rule_combine_dP,
       chain_rule_combine_dτ, chain_rule_combine_dϖ,
       chain_rule_combine_surface_brdf, selected_measurement_jacobian,
       selected_measurements, surface_brdf_wind_jacobian,
       truncated_ss_path1, truncated_ss_path2, apply_back_correction!
# Export new hierarchical model types
export RTModel, AbstractRTModel, SolverConfig, Atmosphere, RayleighScattering,
       AerosolState, LayerResolvedAerosolOptics, Optics, OpticsLin
# Export v0.6 source-term abstraction so `using vSmartMOM` is enough for
# `SolarBeam`, `SurfaceSIF`, `BlackbodySource`, etc. — re-exports the names
# the CoreRT submodule already exports.
export AbstractSource, AbstractPreparedSource,
       NoSource, SourceSet,
       AbstractSourceADMode, AnalyticSourceJacobian, ForwardDiffSourceJacobian,
       NoSourceJacobian, source_ad_mode,
       SolarBeam, PreparedSolarBeam, BlackbodySource,
       SurfaceSIF, PreparedSurfaceSIF,
       ThermalEmission, PreparedThermalEmission, has_thermal_emission, contribute!,
       prepare_source, prepare_sources, surface_source_contribute!
# Export GEOSChem/NetCDF integration
export GeosChemSource, NetCDFGridSource, NetCDFSource, geoschem_to_dict, read_geoschem_profile

using .Architectures
using .Absorption

# Module initialization - GPU setup happens in optional backend extensions.
function __init__()
    # Nothing to do here; GPU detection is handled by backend extensions.
end

end
