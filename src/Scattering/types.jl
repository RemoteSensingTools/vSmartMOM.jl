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

"""
    AbstractAnalyticPhaseFunction

Analytic phase/scattering matrix source that can be converted to Greek
coefficients and then used by the standard MOM optical-property path.
"""
abstract type AbstractAnalyticPhaseFunction end

"""
    HenyeyGreensteinPhaseFunction(g)

Scalar Henyey-Greenstein phase function,
`(1 - g^2) / (1 + g^2 - 2g cosΘ)^(3/2)`, normalized so its sphere average is
one.
"""
Base.@kwdef struct HenyeyGreensteinPhaseFunction{FT<:Real} <: AbstractAnalyticPhaseFunction
    "Henyey-Greenstein asymmetry parameter; must satisfy `abs(g) < 1`."
    g::FT
end

"""
    SyntheticPolarizedHenyeyGreensteinPhaseFunction(; g, polarization_fraction)

Diagnostic polarizing Henyey-Greenstein-like scattering matrix. The `f11`
element is standard Henyey-Greenstein; `f12/f11` follows the bounded toy law
`polarization_fraction * (1 - cosΘ^2) / (1 + cosΘ^2)`. This is intended for
tests and sensitivity experiments, not as a Mie substitute.
"""
Base.@kwdef struct SyntheticPolarizedHenyeyGreensteinPhaseFunction{
    FT<:Real, PF<:Real
} <: AbstractAnalyticPhaseFunction
    "Henyey-Greenstein asymmetry parameter; must satisfy `abs(g) < 1`."
    g::FT
    "Maximum synthetic fractional linear polarization; must satisfy `abs(p) <= 1`."
    polarization_fraction::PF
end

"""
    Aerosol

Aerosol microphysical properties: particle size distribution and complex
refractive index.  Used as input to [`MieModel`](@ref) for computing
single-scattering optical properties via Lorenz-Mie theory.

The refractive index convention is `n = nᵣ - i·nᵢ`, where positive `nᵢ`
indicates absorption.

# Fields
- `size_distribution::ContinuousUnivariateDistribution`: Particle radius distribution (e.g., `LogNormal`). Units: μm.
- `nᵣ`: Real part of the refractive index (relative to air).
- `nᵢ`: Imaginary part of the refractive index (absorption).

# Example
```julia
using Distributions
aer = Aerosol(LogNormal(log(0.3), 0.4), 1.3, 0.01)
```
"""
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
    AbstractTruncationType

Abstract supertype for phase-function truncation methods. All concrete
methods are dispatched through [`truncate_phase`](@ref) and supply
`l_max(t)` (the per-band Legendre cutoff that the RT pipeline allocates
for). Subtypes:

* [`NoTruncation`](@ref) — identity. Use when the phase function has
  no sharp forward peak (canopy, isotropic scattering, smooth Rayleigh).
* [`δBGE`](@ref) — δ-BGE-fit (Sanghavi & Stephens 2015, JQSRT 159
  §3); recommended for hyperspectral retrievals.

The atmospheric `Δ_angle` (forward exclusion half-angle) lives inside
the truncation type that needs it, not as a free parameter on
[`vSmartMOM_Parameters`](@ref) — different methods have different
hyper-parameters and `NoTruncation` has none.
"""
abstract type AbstractTruncationType end

"""
    NoTruncation(; l_max=typemax(Int))

Identity truncation — phase functions are passed through unchanged.

This is the correct choice for radiative transfer through media whose
phase function has no sharp forward peak: canopy bi-Lambertian
scattering (the `f_tr → 0` limit of Sanghavi & Stephens 2015 Eq. 8 is
exactly the identity), isotropic scattering, and smooth Rayleigh.
For Mie aerosol or ice-cloud forward peaks use [`δBGE`](@ref) instead.
"""
Base.@kwdef struct NoTruncation <: AbstractTruncationType
    "Per-band Legendre cutoff used by the RT pipeline. Defaults to
    `typemax(Int)`, which is interpreted downstream as 'use the full
    Greek-coefficient length'."
    l_max::Int = typemax(Int)
end

"""
    δBGE{FT}(l_max, Δ_angle)

δ-BGE-fit truncation, vector form (Sanghavi & Stephens 2015, JQSRT 159,
§3 — extension of Hu et al. 2000 to vector RT). Fits truncated Legendre
coefficients outside the forward exclusion cone of half-angle
`Δ_angle` and renormalises by the retained scattering fraction
``c_0 = 1 - f^t``.

Recommended over plain δ-m for hyperspectral retrievals because δ-m
has known DSE and PTE errors near exact backscatter (Sanghavi &
Stephens 2015 §2.4).

# Fields
$(DocStringExtensions.FIELDS)
"""
struct δBGE{FT} <: AbstractTruncationType
    "Truncation length for Legendre terms"
    l_max::Int
    "Exclusion angle for forward peak (in fitting procedure) `[degrees]`"
    Δ_angle::FT
end

"""
    l_max(t::AbstractTruncationType) -> Int

Per-band Legendre cutoff supplied by the truncation type. RT pipeline
code that needs a finite cutoff with `NoTruncation` should clamp to the
actual Greek-coefficient length, e.g.
`min(length(greek.β), l_max(truncation))`.
"""
@inline l_max(t::AbstractTruncationType) = t.l_max

#=

Model that specifies the Mie computation details

=#

"""
    MieModel{FDT<:AbstractFourierDecompositionType, FT}

Configuration for a Lorenz–Mie scattering computation.  Specifies the aerosol
(size distribution + refractive index), wavelength, polarization type,
truncation strategy, and integration parameters.  The `computation_type`
selects between NAI-2 (Siewert) and PCW (Domke) Fourier decomposition
algorithms.

Pre-computed Wigner symbol tables (`wigner_A`, `wigner_B`) can be supplied
for PCW; they default to trivial placeholders when unused.

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
    GreekCoefs{FT<:Real}

Expansion coefficients of the 4×4 scattering (phase) matrix in generalised
spherical functions (the "Greek" coefficients).  Six independent coefficient
vectors (`α, β, γ, δ, ϵ, ζ`) fully describe the azimuthal Fourier
decomposition of the scattering matrix **B** for a given particle or mixture.
See Eq. 16 in Sanghavi (2014) for the mapping to **B** elements.

For scalar (intensity-only) RT, only `β` (the phase-function expansion) is
used.

# Fields
$(DocStringExtensions.FIELDS)
"""
mutable struct GreekCoefs{FT<:Real}
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
    AerosolOptics{FT<:Real}

Computed aerosol single-scattering optical properties for one aerosol type
at one (or more) wavelengths.  Produced by integrating the Mie solution over
the particle size distribution.

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct AerosolOptics{FT<:Real}
    "Greek matrix"
    greek_coefs::GreekCoefs
    "Single Scattering Albedo"
    ω̃::Union{FT, AbstractArray{FT}}
    "Extinction cross-section"
    k::Union{FT, AbstractArray{FT}}
    "Truncation factor" 
    fᵗ::Union{FT, AbstractArray{FT}}
    "Derivatives"
    derivs = zeros(1)
end

""" Extend Base.isapprox (≈) to compare two AerosolOptics """
function Base.:isapprox(aerosol_optics_a::AerosolOptics, aerosol_optics_b::AerosolOptics) 
    field_names = fieldnames(AerosolOptics)
    return all([getproperty(aerosol_optics_a, field) ≈ getproperty(aerosol_optics_b, field) for field in field_names])
end
