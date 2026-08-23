#=
 
This file contains all types that are used in the Scattering module:

- `AbstractAerosolTypes` specify aerosol properties
- `AbstractFourierDecompositionTypes` specify the decomposition method to use (NAI2 vs PCW)
- `AbstractPolarizationTypes` specify the polarization type (I/IQ/IQU/IQUV)
- `AbstractTruncationTypes` specify the type of truncation for legendre terms
- `GreekCoefs` holds all greek coefficients 
- `ScatteringMatrix` holds all computed phase function elements
- `AerosolOptics` holds all computed aerosol optical properties

=#
# TODO: struct MultivariateAerosol{FT,FT2} <: AbstractAerosolType


"""
    type AbstractFourierDecompositionType
Abstract aerosol Fourier Decomposition computation type (NAI2 and PCW)
"""
# abstract type AbstractFourierDecompositionType end

"""
    type NAI2

Perform Siewart's numerical integration method, NAI-2, to compute aerosol phase function 
decomposition. See: http://adsabs.harvard.edu/full/1982A%26A...109..195S
"""
# struct NAI2  <: AbstractFourierDecompositionType end

"""
    type PCW

Perform Domke's Precomputed Wigner Symbols method, PCW, to compute aerosol phase function 
decomposition. See: http://adsabs.harvard.edu/full/1984A%26A...131..237D
"""
# struct PCW <: AbstractFourierDecompositionType end


"""
    struct GreekCoefs{FT}

A struct which holds all Greek coefficient lists (over l) in one object. 
See eq 16 in Sanghavi 2014 for details. 

# Fields
$(DocStringExtensions.FIELDS)
"""
mutable struct linGreekCoefs{FT<:AbstractFloat}
    "Greek matrix coefficient α, is in B[2,2]"
    α̇::AbstractArray{FT}
    "Greek matrix coefficient β, is in `B[1,1]` (only important one for scalar!)"
    β̇::AbstractArray{FT}
    "Greek matrix coefficient γ, is in B[2,1],B[1,2]"
    γ̇::AbstractArray{FT}
    "Greek matrix coefficient δ, is in B[4,4]"
    δ̇::AbstractArray{FT}
    "Greek matrix coefficient ϵ, is in B[3,4] and - in B[4,3]"
    ϵ̇::AbstractArray{FT}
    "Greek matrix coefficient ζ, is in B[3,3]"
    ζ̇::AbstractArray{FT}
end
""" 
    struct ScatteringMatrix

A struct which holds all computed phase function elements. 
f₁₁ represents the phase function p for the Intensity (first Stokes Vector element) and is normalized as follows:
1/4π ∫₀²⁽ᵖⁱ⁾ dϕ ∫₋₁¹ p(μ) dμ  = 1
    
# Fields
$(DocStringExtensions.FIELDS)
""" 
struct linScatteringMatrix{FT<:AbstractFloat}
    ḟ₁₁::AbstractArray{FT}
    ḟ₁₂::AbstractArray{FT}
    ḟ₂₂::AbstractArray{FT}
    ḟ₃₃::AbstractArray{FT}
    ḟ₃₄::AbstractArray{FT}
    ḟ₄₄::AbstractArray{FT}
end

"""
    struct AerosolOptics

A struct which holds all computed aerosol optics

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct linAerosolOptics{FT<:AbstractFloat}
    "Greek matrix"
    lin_greek_coefs::linGreekCoefs
    "Single Scattering Albedo"
    ω̃̇::AbstractArray{FT}
    "Extinction cross-section"
    k̇::AbstractArray{FT}
    "Truncation factor" 
    ḟᵗ::AbstractArray{FT}
    "Wavenumber nodes corresponding to `phase_lin_greek`"
    phase_ν::Union{Nothing,Vector{FT}} = nothing
    "Post-truncation Greek-coefficient tangents retained at phase nodes"
    phase_lin_greek::Union{Nothing,Vector{linGreekCoefs{FT}}} = nothing
    #"Derivatives"
    #derivs = zeros(1)
end
