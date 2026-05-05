abstract type AbstractSSContributor end
abstract type AbstractSSSurface end

"""
    SSGeometry(μ₀, μv, Δϕ)

Solar and viewing geometry for the standalone single-scatter solver.
`μ₀` is the positive solar cosine, `μv` are positive upward viewing cosines,
and `Δϕ` is the view-sun relative azimuth in radians for each view.
"""
Base.@kwdef struct SSGeometry{FT<:Real, V<:AbstractVector{FT}}
    μ₀::FT
    μv::V
    Δϕ::V
end

"""
    LambertianSSSurface(albedo)

Lambertian surface used by Phase 1 path 2. `albedo` may be a scalar or a
per-spectral vector.
"""
Base.@kwdef struct LambertianSSSurface{A} <: AbstractSSSurface
    albedo::A
end

Base.@kwdef struct RayleighSSContributor{A<:AbstractMatrix} <: AbstractSSContributor
    "Optical depth by layer and spectral point: (nLayer, nSpec)."
    τ::A
end

Base.@kwdef struct HGAerosolSSContributor{FT<:Real, A<:AbstractMatrix} <: AbstractSSContributor
    "Henyey-Greenstein asymmetry parameter."
    g::FT
    "Single-scattering albedo."
    ϖ::FT
    "Optical depth by layer and spectral point: (nLayer, nSpec)."
    τ::A
end

Base.@kwdef struct AbsorptionSSContributor{A<:AbstractMatrix} <: AbstractSSContributor
    "Absorption optical depth by layer and spectral point: (nLayer, nSpec)."
    τ::A
end

"""
    TruncatedAndExactScatteringOpticalProperties(truncated, exact, truncfac)

Dual-form optical-property container used by the broader exact-SFI design.
Phase 1 defines the type for API alignment but does not wire it into `rt_run`.
"""
struct TruncatedAndExactScatteringOpticalProperties{T, E, F}
    truncated::T
    exact::E
    truncfac::F
end

"""
    ExactSSConfig(; geometry, surface, contributors, I0, n_stokes=1)

Configuration for the Phase 1 standalone exact single-scatter solver.
Outputs are scalar Stokes-I arrays with shape `(nGeom, 1, nSpec)`.
"""
Base.@kwdef struct ExactSSConfig{GEO, SUR, CTRB, I}
    geometry::GEO
    surface::SUR
    contributors::CTRB
    I0::I
    n_stokes::Int = 1
    inner_nquad::Int = 16
    azimuth_nquad::Int = 64
end

τ_matrix(c::AbstractSSContributor) = c.τ

single_scattering_albedo(::RayleighSSContributor) = 1
single_scattering_albedo(c::HGAerosolSSContributor) = c.ϖ
single_scattering_albedo(::AbsorptionSSContributor) = 0

exact_phase_function(::RayleighSSContributor, cosΘ) =
    oftype(cosΘ, 0.75) * (one(cosΘ) + cosΘ^2)

function exact_phase_function(c::HGAerosolSSContributor, cosΘ)
    g = convert(typeof(cosΘ), c.g)
    return (one(g) - g^2) / (one(g) + g^2 - 2g * cosΘ)^(oftype(g, 1.5))
end

exact_phase_function(::AbsorptionSSContributor, cosΘ) = zero(cosΘ)
