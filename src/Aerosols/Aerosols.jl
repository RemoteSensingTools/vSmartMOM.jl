# WIP: This module is user-facing but the API may evolve in follow-up PRs.
# Landed as part of the sanghavi-unified merge to close the
# `using vSmartMOM.Aerosols` export gap; further cleanup is a known
# follow-up workstream.
"""
    Aerosols Module

Flexible aerosol framework supporting multiple schemes:
- TOMAS: Size-resolved aerosol microphysics (12/15/30/40 bins)
- Two-Moment: Bulk aerosol properties (AOD + effective radius)
- Future: MAAM, custom schemes, etc.

Includes wavelength-dependent refractive index database and optical property calculations.
"""
module Aerosols

using NCDatasets
using YAML
using Interpolations
using Statistics
using Printf
using FastGaussQuadrature: gausslegendre
using ..Scattering

# Export main types
export AerosolScheme, TOMASScheme, TOMAS15Scheme, TwoMomentScheme
export AerosolData, AerosolSpeciesData
export RefractiveIndexLUT, RefractiveIndexDatabase

# Export abstract-type hierarchy (gchp-io PR — see abstract_types.jl)
export AbstractSizeDistribution, LogNormalSizeDistribution,
       GammaSizeDistribution, DiscreteBinSizeDistribution
export AbstractBinIntegration, LogNormalFit,
       ConstantIntegrationPerBin, LinearIntegrationPerBin,
       DirectBinSum
export AbstractRefractiveIndex, ConstantRI, WavelengthDependentRI
export SpeciesComposition
export AbstractMixingRule, ExternalMixing, VolumeWeightedMixing,
       MaxwellGarnettMixing, BruggemanMixing
export AbstractAerosolScheme
export AerosolSizeGrid
export AbstractAerosolBinData, SectionalAerosolData
export nbins, nlayers, species_list, bin_composition, densities
export refractive_index_key, read_aerosol_cell
export effective_ri, external_ri_components, compute_aerosol_optics
export compute_column_aod

# Export main functions
export read_aerosol_data, read_tomas, read_tomas15
export load_refractive_index_database
export get_refractive_index
export list_species
export wavelength_range
export compute_optical_properties

# Include submodules
include("types.jl")
include("refractive_index.jl")
include("abstract_types.jl")
include("readers.jl")
include("schemes/tomas15.jl")
include("schemes/two_moment.jl")
include("aod_diagnostics.jl")
include("optical_properties.jl")

end # module
