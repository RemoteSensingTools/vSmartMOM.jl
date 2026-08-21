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

The struct is parameterized by the refractive-index element type `FT`. `FT` is
typically a plain `AbstractFloat` (`Float64`/`Float32`) but may also be a
`ForwardDiff.Dual` so the autodiff path can track derivatives with respect to
`nᵣ`/`nᵢ`; for that reason `FT` is intentionally **not** constrained to
`<:AbstractFloat`. The outer convenience constructor `Aerosol(dist, nᵣ, nᵢ)`
promotes `nᵣ` and `nᵢ` to a common type, so existing call sites continue to
work unchanged.

# Fields
- `size_distribution::ContinuousUnivariateDistribution`: Particle radius distribution (e.g., `LogNormal`). Units: μm.
- `nᵣ::FT`: Real part of the refractive index (relative to air).
- `nᵢ::FT`: Imaginary part of the refractive index (absorption).

# Example
```julia
using Distributions
aer = Aerosol(LogNormal(log(0.3), 0.4), 1.3, 0.01)
```
"""
mutable struct Aerosol{FT}
    "Univariate size distribution"
    size_distribution::ContinuousUnivariateDistribution
    "Real part of refractive index"
    nᵣ::FT
    "Imag part of refractive index"
    nᵢ::FT
end

# Outer convenience constructor: promote nᵣ/nᵢ to a common element type so that
# all existing `Aerosol(dist, nᵣ, nᵢ)` call sites keep working unchanged. We use
# `promote` (not `promote_type ∘ AbstractFloat`) so that ForwardDiff `Dual`
# arguments are preserved for the autodiff path.
function Aerosol(size_distribution::ContinuousUnivariateDistribution, nᵣ, nᵢ)
    nᵣp, nᵢp = promote(nᵣ, nᵢ)
    return Aerosol{typeof(nᵣp)}(size_distribution, nᵣp, nᵢp)
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
    struct Stokes_IQ{FT<:AbstractFloat}

A struct which defines Stokes Vector ([I,Q]) RT code.

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct Stokes_IQ{FT<:AbstractFloat} <: AbstractPolarizationType
    "Number of Stokes components (int)"
    n::Int = 2
    "Vector of length `n` for ... (see eq in Sanghavi )"
    D::Array{FT}  = [1.0, 1.0]
    "Incoming Stokes vector for scalar only"
    I₀::Array{FT} = [1.0, 0.0] #assuming linearly unpolarized incident stellar radiation
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
[`CoreRT.vSmartMOM_Parameters`](@ref vSmartMOM.CoreRT.vSmartMOM_Parameters)
— different methods have different hyper-parameters and `NoTruncation` has
none.
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
    AutoTruncation()

Phase D — deferred-decision marker for `truncation: auto` in YAML
(mirrors VLIDORT's `DO_DELTAM_SCALING` philosophy). At
`model_from_parameters` time it resolves deterministically:

- No aerosol scattering, or all aerosols' `length(greek.β) - 1`
  fits within `stream_l_cap` ⇒ `NoTruncation()` (typical for
  Rayleigh-only or smooth-aerosol scenes).
- Aerosols with `phase_lmax > stream_l_cap` ⇒
  `δBGE(stream_l_cap, Δ_angle)`.

The resolver emits an `@info` line stating which branch was taken
so the user can always see what was applied. `AutoTruncation` is
**deliberately not threaded through the Mie/RT kernels** — it is a
build-time placeholder that gets replaced before any kernel sees it.

User-facing knobs in YAML:

| Value                     | Meaning                                            |
|---------------------------|----------------------------------------------------|
| `truncation: auto`        | This deferred-decision mode (Phase D recommended)  |
| `truncation: NoTruncation()` / `null` | Exactly no transform; errors if coefs exceed cap |
| `truncation: δBGE(N, Δ)`  | Explicit; used for benchmarks / cross-validation   |
"""
struct AutoTruncation <: AbstractTruncationType end

"""
    l_max(t::AbstractTruncationType) -> Int

Per-band Legendre cutoff supplied by the truncation type. RT pipeline
code that needs a finite cutoff with `NoTruncation` should clamp to the
actual Greek-coefficient length, e.g.
`min(length(greek.β), l_max(truncation))`. `AutoTruncation()` reports
`typemax(Int)` because it is unresolved — the model builder replaces it
before any kernel runs.
"""
@inline l_max(t::AbstractTruncationType) = t.l_max
@inline l_max(::AutoTruncation) = typemax(Int)

#=

Model that specifies the Mie computation details

=#

"""
    MieModel{FDT<:AbstractFourierDecompositionType, FT, ARCH}

Configuration for a Lorenz–Mie scattering computation.  Specifies the aerosol
(size distribution + refractive index), wavelength, polarization type,
truncation strategy, integration parameters, and the compute architecture.
The `computation_type` selects between NAI-2 (Siewert) and PCW (Domke) Fourier
decomposition algorithms.

Pre-computed Wigner symbol tables (`wigner_A`, `wigner_B`) can be supplied
for PCW; they default to trivial 1×1×1 placeholders (of element type `FT`)
when unused (NAI2 path).

# FT semantics (three orthogonal precision axes)

`FT` is the **output** float type of the returned Greek coefficients and
optical scalars.  It is carried by `λ`, `r_max`, and `aerosol` (which must all
share the same `FT`) so that the entire `MieModel` is FT-consistent. This is
independent of two other precision choices:

1. `FT` (output type) — what the user consumes downstream in the RT pipeline.
2. The internal `Dₙ` continued-fraction recursion: on the CPU path the
   recursion is always promoted to `Float64` for numerical stability regardless
   of `FT` (see `_mie_dn_recursion!` in `mie_helper_functions.jl`);
   `Dual` inputs keep native arithmetic for AD compatibility.
3. `precision_policy` (GPU only) — `NativeFloat64` (default on CUDA) runs the
   GPU `Dₙ` recursion in hardware FP64; `DSEmulated` uses Float32
   double-single pairs (for Metal/L40S). This axis is inert on the CPU path.

# Architecture dispatch

`architecture::ARCH` is either `Architectures.CPU()` (default) or
`Architectures.GPU()`. `compute_aerosol_optical_properties(mie_model)`
dispatches on it: NAI2+CPU → analytic CPU path; NAI2+GPU → the
KernelAbstractions GPU pipeline; PCW+GPU → CPU fallback (no GPU PCW kernel).

# Fields
$(DocStringExtensions.FIELDS)
"""
struct MieModel{FDT<:AbstractFourierDecompositionType, FT, ARCH}

    computation_type::FDT
    "Aerosol microphysics; element type must match the model's `FT`"
    aerosol::Aerosol{FT}
    "Wavelength `[μm]`; must have the same float type as `r_max`"
    λ::FT
    polarization_type::AbstractPolarizationType
    truncation_type::AbstractTruncationType

    "Maximum radius `[μm]`"
    r_max::FT
    "Number of quadrature points for integration over size distribution"
    nquad_radius::Int

    "Precomputed Wigner ν=0 table (PCW only; trivial 1×1×1 placeholder for NAI2)"
    wigner_A::Array{FT,3}
    "Precomputed Wigner ν=2 table (PCW only; trivial 1×1×1 placeholder for NAI2)"
    wigner_B::Array{FT,3}

    "Compute architecture (`CPU()` or `GPU()`); selects CPU vs GPU Mie path"
    architecture::ARCH
    "GPU precision policy (`NativeFloat64`/`DSEmulated`); `nothing` = auto-select on GPU, ignored on CPU"
    precision_policy

end

# Keyword constructor — preserves the public `MieModel(; …)` API that the former
# `Base.@kwdef` provided, while enforcing the FT-consistent field types. `λ`, `r_max`
# and the (placeholder) Wigner tables are promoted to the aerosol's float type `FT`,
# exactly as `make_mie_model` does, so an `Aerosol{FT}` always yields a `MieModel{…,FT,…}`.
# (FT is read from the aerosol inside the body: a `where {FT}` param cannot be bound from
# a keyword argument — only positional args do that — so the type annotation route fails.)
function MieModel(; computation_type, aerosol, λ, polarization_type, truncation_type,
                  r_max, nquad_radius, wigner_A = nothing, wigner_B = nothing,
                  architecture = Architectures.CPU(), precision_policy = nothing)
    FT = typeof(aerosol.nᵣ)   # aerosol::Aerosol{FT}
    wA = wigner_A === nothing ? zeros(FT, 1, 1, 1) : convert(Array{FT,3}, wigner_A)
    wB = wigner_B === nothing ? zeros(FT, 1, 1, 1) : convert(Array{FT,3}, wigner_B)
    return MieModel(computation_type, aerosol, FT(λ), polarization_type, truncation_type,
                    FT(r_max), nquad_radius, wA, wB, architecture, precision_policy)
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
    α::AbstractArray{FT}
    "Greek matrix coefficient β, is in `B[1,1]` (only important one for scalar!)"
    β::AbstractArray{FT}
    "Greek matrix coefficient γ, is in B[2,1],B[1,2]"
    γ::AbstractArray{FT}
    "Greek matrix coefficient δ, is in B[4,4]"
    δ::AbstractArray{FT}
    "Greek matrix coefficient ϵ, is in B[3,4] and - in B[4,3]"
    ϵ::AbstractArray{FT}
    "Greek matrix coefficient ζ, is in B[3,3]"
    ζ::AbstractArray{FT}
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
    "Wavenumber nodes (cm⁻¹) at which the post-truncation phase coefficients were evaluated"
    phase_ν::Union{Nothing,Vector{FT}} = nothing
    "Post-truncation Greek coefficients at `phase_ν`; retained so Z is evaluated only at these nodes"
    phase_greek::Union{Nothing,Vector{GreekCoefs{FT}}} = nothing
end

""" Extend Base.isapprox (≈) to compare two AerosolOptics """
function Base.:isapprox(aerosol_optics_a::AerosolOptics, aerosol_optics_b::AerosolOptics) 
    field_names = fieldnames(AerosolOptics)
    return all([getproperty(aerosol_optics_a, field) ≈ getproperty(aerosol_optics_b, field) for field in field_names])
end
