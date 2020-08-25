using Distributions
using DocStringExtensions

"""
    type AbstractAerosolType
Abstract aerosol type 
"""
abstract type AbstractAerosolType end

"""
    struct UnivariateAerosol{FT}

A struct which provides all model parameters needed for cross-section 
calculations using HITRAN data

# Fields
$(DocStringExtensions.FIELDS)
"""
struct UnivariateAerosol{FT} <: AbstractAerosolType
    "Univariate size distribution"
    size_distribution::ContinuousUnivariateDistribution
    "Maximum radius `[μm]`"
    r_max::FT
    "Number of quadrature points for integration over size distribution"
    nquad_radius::Int
    "Real part of refractive Index"
    nᵣ::FT
    "Imaginary part of refractive Index"
    nᵢ::FT
end