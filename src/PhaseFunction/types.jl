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
    type AbstractPolarizationType
Abstract Polarization type 
"""
abstract type AbstractPolarizationType end

# Use full Stokes Vector:
struct FullStokes <: AbstractPolarizationType end

# Use scalar only:
struct Scalar <: AbstractPolarizationType end

"""
    type AbstractTruncationType
Abstract greek coefficient truncation type 
"""
abstract type AbstractTruncationType end

struct δBGE{FT} <: AbstractTruncationType
    "Trunction length for legendre terms"
    l_max::Int
    "Exclusion angle for forward peak (in fitting procedure) `[degrees]`"
    Δ_angle::FT
end


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