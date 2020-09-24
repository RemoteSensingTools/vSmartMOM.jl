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

struct Domke{FT} <: AbstractGreekType end

"""
    type AbstractQuadratureType
Abstract Quadrature type 
"""
abstract type AbstractQuadratureType end


"""
    type AbstractPolarizationType
Abstract Polarization type 
"""
abstract type AbstractPolarizationType end

"Use full Stokes Vector ([I,Q,U,V])"
struct Stokes_IQUV <: AbstractPolarizationType end

"Use part of Stokes Vector ([I,Q,U])"
struct Stokes_IQU <: AbstractPolarizationType end

"Use scalar only ([I]):"
struct Stokes_I <: AbstractPolarizationType end

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

A struct which provides all univariate aerosol parameters needed for Mie 
computation

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

"""
    struct UnivariateAerosol{FT}

A struct which holds all Greek coefficient lists (over l) in one object

# Fields
$(DocStringExtensions.FIELDS)
"""
struct GreekCoefs{FT}
    α::Array{FT,1}
    β::Array{FT,1}
    γ::Array{FT,1}
    δ::Array{FT,1}
    ϵ::Array{FT,1}
    ζ::Array{FT,1}
end