#####
##### Types for describing an aerosol's composition for Mie computation
#####

"""
    type AbstractAerosolType
Abstract aerosol type 
"""
abstract type AbstractAerosolType end

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

# TODO: struct MultivariateAerosol{FT,FT2} <: AbstractAerosolType

#####
##### Types of Fourier Decomposition (NAI2 or PCW)
#####

"""
    type AbstractFourierDecompositionType
Abstract aerosol Fourier Decomposition computation type (NAI2 and PCW)
"""
abstract type AbstractFourierDecompositionType end

"""
    type NAI2

Perform Siewart's numerical integration method, NAI-2, to compute aerosol phase function 
decomposition. See: http://adsabs.harvard.edu/full/1982A%26A...109..195S
"""
struct NAI2  <: AbstractFourierDecompositionType end

"""
    type PCW

Perform Domke's Precomputed Wigner Symbols method, PCW, to compute aerosol phase function 
decomposition. See: http://adsabs.harvard.edu/full/1984A%26A...131..237D
"""
struct PCW <: AbstractFourierDecompositionType end

#####
##### Types of Polarization (which Stokes vector to use)
#####

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

#####
##### Types of Truncation (for Legendre terms)
#####

"""
    type AbstractTruncationType

Abstract greek coefficient truncation type 
"""
abstract type AbstractTruncationType end

"""
    type δBGE{FT} <: AbstractTruncationType
        
# Fields
$(DocStringExtensions.FIELDS)
"""
struct δBGE{FT} <: AbstractTruncationType
    "Trunction length for legendre terms"
    l_max::Int
    "Exclusion angle for forward peak (in fitting procedure) `[degrees]`"
    Δ_angle::FT
end

#####
##### Model that specifies the Mie computation details 
#####

"""
    type MieModel

Model to hold all Mie computation details for NAI2 and PCW

# Fields
$(DocStringExtensions.FIELDS)
"""
@with_kw struct MieModel{FDT<:AbstractFourierDecompositionType}

    computation_type::FDT
    aerosol::AbstractAerosolType
    λ
    polarization_type::AbstractPolarizationType
    truncation_type::AbstractTruncationType

    wigner_A = zeros(1, 1, 1)
    wigner_B = zeros(1, 1, 1)

end

#####
##### Types that are needed for the output of the Fourier decomposition
#####

"""
    struct GreekCoefs{FT}

A struct which holds all Greek coefficient lists (over l) in one object. 
See eq 16 in Sanghavi 2014 for details. 

# Fields
$(DocStringExtensions.FIELDS)
"""
struct GreekCoefs{FT}
    "Greek matrix coefficient α, is in B[2,2]"
    α::Array{FT,1}
    "Greek matrix coefficient β, is in B[1,1] (only important one for scalar!)"
    β::Array{FT,1}
    "Greek matrix coefficient γ, is in B[2,1],B[1,2]"
    γ::Array{FT,1}
    "Greek matrix coefficient δ, is in B[4,4]"
    δ::Array{FT,1}
    "Greek matrix coefficient ϵ, is in B[3,4] and - in B[4,3]"
    ϵ::Array{FT,1}
    "Greek matrix coefficient ζ, is in B[3,3]"
    ζ::Array{FT,1}
end

""" Extend Base.isapprox (≈) to compare two GreekCoefs """
function Base.:isapprox(greek_coefs_a::GreekCoefs, greek_coefs_b::GreekCoefs) 
    field_names = fieldnames(GreekCoefs)
    return all([getproperty(greek_coefs_a, field) ≈ getproperty(greek_coefs_b, field) for field in field_names])
end

"""
    struct AerosolOptics

A struct which holds all computed aerosol optics

# Fields
$(DocStringExtensions.FIELDS)
"""
@with_kw struct AerosolOptics{FT}
    "Greek matrix"
    greek_coefs::GreekCoefs
    "Single Scattering Albedo"
    ω̃::FT
    "Extinction coefficient"
    k::FT
    "Derivatives"
    derivs = zeros(1)
end

""" Extend Base.isapprox (≈) to compare two AerosolOptics """
function Base.:isapprox(aerosol_optics_a::AerosolOptics, aerosol_optics_b::AerosolOptics) 
    field_names = fieldnames(AerosolOptics)
    return all([getproperty(aerosol_optics_a, field) ≈ getproperty(aerosol_optics_b, field) for field in field_names])
end



# """
#     type AbstractQuadratureType
# Abstract Quadrature type 
# """
# abstract type AbstractQuadratureType end