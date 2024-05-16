#=
 
This file contains all types that are used in the Scattering module:

- `AbstractAerosolTypes` specify aerosol properties
- `AbstractFourierDecompositionTypes` specify the decomposition method to use (NAI2 vs PCW)
- `AbstractPolarizationTypes` specify the polarization type (I/IQU/IQUV)
- `AbstractTruncationTypes` specify the type of truncation for legendre terms
- `GreekCoefs` holds all greek coefficients 
- `ScatteringMatrix` holds all computed phase function elements
- `AerosolOptics` holds all computed aerosol optical properties

=#

"""
    type AbstractAerosolType
Abstract aerosol type 
"""
abstract type AbstractAerosolType end

"Aerosol type with its properties (size distribution and refractive index)"
mutable struct Aerosol{}
    "Univariate size distribution"
    size_distribution::ContinuousUnivariateDistribution
    "Real part of refractive index"
    nᵣ
    "Imag part of refractive index"
    nᵢ
end

# TODO: struct MultivariateAerosol{FT,FT2} <: AbstractAerosolType

#=

Types of Fourier Decomposition (NAI2 or PCW)

=#

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

#=

Types of Polarization (which Stokes vector to use)

=#

"""
    type AbstractPolarizationType

Abstract Polarization type 
"""
abstract type AbstractPolarizationType  end

"""
    struct Stokes_IQUV{FT<:AbstractFloat}

A struct which defines full Stokes Vector ([I,Q,U,V]) RT code

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct Stokes_IQUV{FT<:AbstractFloat} <: AbstractPolarizationType
    "Number of Stokes components (int)"
    n::Int = 4
    "Vector of length `n` for ... (see eq in Sanghavi )"
    D::Array{FT}  = [1.0, 1.0, -1.0, -1.0]
    "Incoming Stokes vector for scalar only"
    I₀::Array{FT} = [1.0, 0.0, 0.0, 0.0] #assuming completely unpolarized incident stellar radiation
end

"""
    struct Stokes_IQU{FT<:AbstractFloat}

A struct which defines Stokes Vector ([I,Q,U]) RT code

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct Stokes_IQU{FT<:AbstractFloat} <: AbstractPolarizationType
    "Number of Stokes components (int)" 
    n::Int = 3
    "Vector of length `n` for ... (see eq in Sanghavi )"
    D::Array{FT}  = [1.0, 1.0, -1.0]
    "Incoming Stokes vector for scalar only"
    I₀::Array{FT} = [1.0, 0.0, 0.0] #assuming linearly unpolarized incident stellar radiation
end

"""
    struct Stokes_I{FT<:AbstractFloat}

A struct which define scalar I only RT code

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct Stokes_I{FT<:AbstractFloat} <: AbstractPolarizationType 
    "Number of Stokes components (int)"
    n::Int = 1
    "Vector of length `n` for ... (see eq in Sanghavi )"
    D::Array{FT} = [1.0]
    "Incoming Stokes vector for scalar only"
    I₀::Array{FT} = [1.0]
end

#=

Types of Truncation (for Legendre terms)

=#

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

#=

Model that specifies the Mie computation details

=#

"""
    type MieModel

Model to hold all Mie computation details for NAI2 and PCW

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct MieModel{FDT<:AbstractFourierDecompositionType, FT}

    computation_type::FDT
    aerosol::Aerosol
    λ::Real
    polarization_type::AbstractPolarizationType
    truncation_type::AbstractTruncationType

    "Maximum radius `[μm]`"
    r_max::FT
    "Number of quadrature points for integration over size distribution"
    nquad_radius::Int

    wigner_A = zeros(1, 1, 1)
    wigner_B = zeros(1, 1, 1)

end

#=

Types that are needed for the output of the Fourier decomposition

=#

"""
    struct GreekCoefs{FT}

A struct which holds all Greek coefficient lists (over l) in one object. 
See eq 16 in Sanghavi 2014 for details. 

# Fields
$(DocStringExtensions.FIELDS)
"""
struct GreekCoefs{FT<:Union{AbstractFloat, ForwardDiff.Dual}}
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
    struct ScatteringMatrix

A struct which holds all computed phase function elements. 
f₁₁ represents the phase function p for the Intensity (first Stokes Vector element) and is normalized as follows:
1/4π ∫₀²⁽ᵖⁱ⁾ dϕ ∫₋₁¹ p(μ) dμ  = 1
    
# Fields
$(DocStringExtensions.FIELDS)
""" 
struct ScatteringMatrix{FT}
    f₁₁::FT
    f₁₂::FT
    f₂₂::FT
    f₃₃::FT
    f₃₄::FT
    f₄₄::FT
end

"""
    struct AerosolOptics

A struct which holds all computed aerosol optics

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct AerosolOptics{FT<:Union{AbstractFloat, ForwardDiff.Dual}}
    "Greek matrix"
    greek_coefs::GreekCoefs
    "Single Scattering Albedo"
    ω̃::Union{FT, AbstractArray{FT}}
    "Extinction cross-section"
    k::Union{FT, AbstractArray{FT}}
    "Truncation factor" 
    fᵗ::FT
    "Derivatives"
    derivs = zeros(1)
end

""" Extend Base.isapprox (≈) to compare two AerosolOptics """
function Base.:isapprox(aerosol_optics_a::AerosolOptics, aerosol_optics_b::AerosolOptics) 
    field_names = fieldnames(AerosolOptics)
    return all([getproperty(aerosol_optics_a, field) ≈ getproperty(aerosol_optics_b, field) for field in field_names])
end