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

Base.@kwdef struct RayleighSSContributor{A<:AbstractMatrix, D<:Real} <: AbstractSSContributor
    "Optical depth by layer and spectral point: (nLayer, nSpec)."
    τ::A
    "Rayleigh depolarization factor. `0` gives ideal molecular Rayleigh."
    depol::D = 0
end

Base.@kwdef struct HGAerosolSSContributor{FT<:Real, A<:AbstractMatrix} <: AbstractSSContributor
    "Henyey-Greenstein asymmetry parameter."
    g::FT
    "Single-scattering albedo."
    ϖ::FT
    "Optical depth by layer and spectral point: (nLayer, nSpec)."
    τ::A
end

"""
    GreekCoefsSSContributor(; greek_coefs, ϖ, τ)

Standalone exact-angle scattering contributor backed by the same
`Scattering.GreekCoefs` used by the full MOM solver. This is the preferred
bridge for vector aerosol single scattering in StandaloneSS.
"""
Base.@kwdef struct GreekCoefsSSContributor{
    G<:Scattering.GreekCoefs, W<:Real, A<:AbstractMatrix
} <: AbstractSSContributor
    "Greek phase-matrix coefficients from `Scattering`."
    greek_coefs::G
    "Single-scattering albedo."
    ϖ::W
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
                    polarization_type=nothing, architecture=CPU())

Configuration for the Phase 1 standalone exact single-scatter solver.
Outputs have shape `(nGeom, nStokes, nSpec)`. `polarization_type`, when
provided, determines the Stokes count for the currently supported vector
paths: Rayleigh atmospheric single scatter and direct-beam surface reflection.
"""
Base.@kwdef struct ExactSSConfig{GEO, SUR, CTRB, I, ARCH<:AbstractArchitecture, POL}
    geometry::GEO
    surface::SUR
    contributors::CTRB
    I0::I
    n_stokes::Int = 1
    polarization_type::POL = nothing
    inner_nquad::Int = 16
    azimuth_nquad::Int = 64
    architecture::ARCH = CPU()
end

_config_n_stokes(config::ExactSSConfig) =
    _config_n_stokes(config.polarization_type, config.n_stokes)

_config_n_stokes(::Nothing, n_stokes::Int) = n_stokes

function _config_n_stokes(polarization_type, n_stokes::Int)
    hasproperty(polarization_type, :n) ||
        throw(ArgumentError("ExactSSConfig.polarization_type must expose an `n` field"))
    return Int(getproperty(polarization_type, :n))
end

τ_matrix(c::AbstractSSContributor) = c.τ

single_scattering_albedo(::RayleighSSContributor) = 1
single_scattering_albedo(c::HGAerosolSSContributor) = c.ϖ
single_scattering_albedo(c::GreekCoefsSSContributor) = c.ϖ
single_scattering_albedo(::AbsorptionSSContributor) = 0

@inline function _rayleigh_depolarization_scale(depol::FT) where {FT}
    return (one(FT) - depol) / (one(FT) + depol / FT(2))
end

function exact_phase_function(c::RayleighSSContributor, cosΘ)
    FT = promote_type(typeof(c.depol), typeof(cosΘ))
    x = convert(FT, cosΘ)
    dpl_p = _rayleigh_depolarization_scale(convert(FT, c.depol))
    P₂ = (FT(3) * x^2 - one(FT)) / FT(2)
    return one(FT) + FT(0.5) * dpl_p * P₂
end

function exact_phase_function(c::HGAerosolSSContributor, cosΘ)
    FT = promote_type(typeof(c.g), typeof(cosΘ))
    g = convert(FT, c.g)
    x = convert(FT, cosΘ)
    return (one(g) - g^2) / (one(g) + g^2 - 2g * x)^(FT(1.5))
end

function exact_phase_function(c::GreekCoefsSSContributor, cosΘ)
    FT = promote_type(eltype(c.greek_coefs.β), typeof(cosΘ))
    x = convert(FT, cosΘ)
    return Scattering.reconstruct_phase(c.greek_coefs, [x]).f₁₁[1]
end

exact_phase_function(::AbsorptionSSContributor, cosΘ) = zero(cosΘ)
