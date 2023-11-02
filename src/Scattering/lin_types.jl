#=
 
This file contains all types that are used in the Scattering module:

- `GreekCoefs` holds all greek coefficients 
- `ScatteringMatrix` holds all computed phase function elements
- `AerosolOptics` holds all computed aerosol optical properties

=#
#=

Types that are needed for the output of the Fourier decomposition

=#
#=
"""
    struct GreekCoefs{FT}

A struct which holds all Greek coefficient lists (over l) in one object. 
See eq 16 in Sanghavi 2014 for details. 

# Fields
$(DocStringExtensions.FIELDS)
"""
struct dGreekCoefs{FT<:Union{AbstractFloat, ForwardDiff.Dual}}
    "Greek matrix coefficient α, is in B[2,2]"
    dα::Array{FT,2} 
    "Greek matrix coefficient β, is in B[1,1] (only important one for scalar!)"
    dβ::Array{FT,2}
    "Greek matrix coefficient γ, is in B[2,1],B[1,2]"
    dγ::Array{FT,2}
    "Greek matrix coefficient δ, is in B[4,4]"
    dδ::Array{FT,2}
    "Greek matrix coefficient ϵ, is in B[3,4] and - in B[4,3]"
    dϵ::Array{FT,2}
    "Greek matrix coefficient ζ, is in B[3,3]"
    dζ::Array{FT,2}
end


""" Extend Base.isapprox (≈) to compare two GreekCoefs """
function Base.:isapprox(d_greek_coefs_a::dGreekCoefs, d_greek_coefs_b::dGreekCoefs) 
    field_names = fieldnames(dGreekCoefs)
    return all([getproperty(dgreek_coefs_a, field) ≈ getproperty(dgreek_coefs_b, field) for field in field_names])
end

 
=#
"""
    struct ScatteringMatrix

A struct which holds all computed phase function elements. 
f₁₁ represents the phase function p for the Intensity (first Stokes Vector element) and is normalized as follows:
1/4π ∫₀²⁽ᵖⁱ⁾ dϕ ∫₋₁¹ p(μ) dμ  = 1
    
# Fields
$(DocStringExtensions.FIELDS) 
"""

struct dScatteringMatrix{FT}
    df₁₁::Array{FT,2}
    df₁₂::Array{FT,2}
    df₂₂::Array{FT,2}
    df₃₃::Array{FT,2}
    df₃₄::Array{FT,2}
    df₄₄::Array{FT,2}
end

"""
    struct AerosolOptics

A struct which holds all computed aerosol optics

# Fields
$(DocStringExtensions.FIELDS)
"""
#Base.@kwdef 
Base.@kwdef mutable struct dAerosolOptics{FT<:AbstractFloat}
    "derivatives of Greek matrix w.r.t nᵣ, nᵢ, r₀ and σ₀"
    d_greek_coefs::Vector{GreekCoefs{FT}}
    "derivatives of Single Scattering Albedo w.r.t nᵣ, nᵢ, r₀ and σ₀"
    dω̃::Vector{FT}
    "derivatives of Extinction cross-section w.r.t nᵣ, nᵢ, r₀ and σ₀"
    dk::Vector{FT}
    "derivatives of Extinction cross-section at reference wavelength w.r.t nᵣ, nᵢ, r₀ and σ₀"
    dk_ref::Vector{FT}
    "derivatives of Truncation factor w.r.t nᵣ, nᵢ, r₀ and σ₀" 
    dfᵗ::Vector{FT}
    #d_greek_coefs=d_greek_coefs, dω̃=dω̃, dk=d_bulk_C_ext, dk_ref=zeros(FT, 4), dfᵗ=zeros(FT, 4)

end

""" Extend Base.isapprox (≈) to compare two AerosolOptics """
function Base.:isapprox(daerosol_optics_a::dAerosolOptics, daerosol_optics_b::dAerosolOptics) 
    field_names = fieldnames(dAerosolOptics)
    return all([getproperty(daerosol_optics_a, field) ≈ getproperty(daerosol_optics_b, field) for field in field_names])
end