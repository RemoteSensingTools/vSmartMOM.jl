using Distributions
using DocStringExtensions

"""
    type AbstractAerosolType
Abstract aerosol type 
"""
abstract type AbstractAerosolType end

"""
    type AbstractGreekType
Abstract aerosol greeo coefficient computation type 
"""
abstract type AbstractGreekType end

struct NAI2  <: AbstractGreekType end

struct Domke <: AbstractGreekType end

"""
    struct UnivariateAerosol{FT}

A struct which provides all model parameters needed for cross-section 
calculations using HITRAN data

# Fields
$(DocStringExtensions.FIELDS)
"""
struct UnivariateAerosol{FT,FT2} <: AbstractAerosolType
    "Univariate size distribution"
    size_distribution::ContinuousUnivariateDistribution
    "Maximum radius `[μm]`"
    r_max::FT
    "Number of quadrature points for integration over size distribution"
    nquad_radius::Int
    "Real part of refractive Index"
    nᵣ::FT2
    "Imaginary part of refractive Index"
    nᵢ::FT2
end