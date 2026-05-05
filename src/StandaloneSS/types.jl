abstract type AbstractSSContributor end
abstract type AbstractSSSurface end

"""
    SSGeometry(μ₀, μv, Δϕ)

Solar and viewing geometry for the standalone single-scatter solver.
`μ₀` is the positive solar cosine, `μv` are positive upward viewing cosines,
and `Δϕ` is the view-sun relative azimuth in radians for each view.
"""
Base.@kwdef struct SSGeometry{M<:Real, V<:AbstractVector{<:Real}, A<:AbstractVector{<:Real}}
    μ₀::M
    μv::V
    Δϕ::A
end

"""
    LambertianSSSurface(albedo)

Lambertian surface used by Phase 1 path 2. `albedo` may be a scalar or a
per-spectral vector.
"""
Base.@kwdef struct LambertianSSSurface{A} <: AbstractSSSurface
    albedo::A
end

"""
    CoxMunkSSSurface(; wind_speed, n_water=nothing, whitecap_albedo=0.22,
                       include_whitecaps=true, shadowing=true)

Scalar Cox-Munk ocean surface for standalone exact path 2. `wind_speed` is the
10-m wind speed in m/s. `n_water` may be `nothing`, a scalar complex refractive
index, or a per-spectral vector of complex refractive indices. The `nothing`
fallback uses a fixed visible-water refractive index for this standalone path.
"""
struct CoxMunkSSSurface{W<:Real, N, A<:Real} <: AbstractSSSurface
    wind_speed::W
    n_water::N
    whitecap_albedo::A
    include_whitecaps::Bool
    shadowing::Bool
end

function CoxMunkSSSurface(; wind_speed::W, n_water=nothing,
                          whitecap_albedo = nothing,
                          include_whitecaps::Bool = true,
                          shadowing::Bool = true) where {W<:Real}
    albedo = whitecap_albedo === nothing ? W(0.22) : whitecap_albedo
    albedo isa Real ||
        throw(ArgumentError("CoxMunkSSSurface.whitecap_albedo must be real"))
    return CoxMunkSSSurface(wind_speed, n_water, albedo, include_whitecaps,
                            shadowing)
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
    ExactSSConfig(; geometry, surface, contributors, I0, n_stokes=1,
                    architecture=CPU())

Configuration for the Phase 1 standalone exact single-scatter solver.
Outputs are scalar Stokes-I arrays with shape `(nGeom, 1, nSpec)`.
"""
Base.@kwdef struct ExactSSConfig{GEO, SUR, CTRB, I, ARCH<:AbstractArchitecture}
    geometry::GEO
    surface::SUR
    contributors::CTRB
    I0::I
    n_stokes::Int = 1
    inner_nquad::Int = 16
    azimuth_nquad::Int = 64
    architecture::ARCH = CPU()
end

τ_matrix(c::AbstractSSContributor) = c.τ

single_scattering_albedo(::RayleighSSContributor) = 1
single_scattering_albedo(c::HGAerosolSSContributor) = c.ϖ
single_scattering_albedo(::AbsorptionSSContributor) = 0

exact_phase_function(::RayleighSSContributor, cosΘ) =
    oftype(cosΘ, 0.75) * (one(cosΘ) + cosΘ^2)

function exact_phase_function(c::HGAerosolSSContributor, cosΘ)
    FT = promote_type(typeof(c.g), typeof(cosΘ))
    g = convert(FT, c.g)
    x = convert(FT, cosΘ)
    return (one(g) - g^2) / (one(g) + g^2 - 2g * x)^(FT(1.5))
end

exact_phase_function(::AbsorptionSSContributor, cosΘ) = zero(cosΘ)
