#=

`SurfaceSIF` — first-class surface-emission source for the v0.6 source-term
refactor. Represents Solar-Induced (chlorophyll) Fluorescence as an
isotropic emission at the lower boundary that adds to the surface layer's
upwelling source vector `j₀⁻`.

This is the "surfaces produce R/T only; sources contribute to surface
j₀± via dispatch" pattern (Christian's design directive). The old
`RS_type.SIF₀` surface-emission path is intentionally not used by `rt_run`;
SIF now enters through clean double-dispatch:

  surface_source_contribute!(prepared_sources, surface, surface_added_layer, m, pol_type, arch)
    └─ dispatch on (source-type, surface-type) pairs to the right kernel.

# Source / surface dispatch table (Phase 5a — initial population)

| Source                | Surface                                        | Body                                                 |
|-----------------------|------------------------------------------------|------------------------------------------------------|
| `PreparedNoSource`    | any                                            | no-op                                                |
| `PreparedSourceSet`   | any                                            | iterate over members                                 |
| `PreparedSolarBeam`   | any                                            | no-op (solar surface reflection still in `create_surface_layer!`; moves out in 5c) |
| `PreparedSurfaceSIF`  | `LambertianSurfaceScalar/Legendre/Spline`      | factor 2 · SIF₀ broadcast to all Nquad streams (m=0) |
| `PreparedSurfaceSIF`  | other (`rpv`, `RossLi`, `CoxMunk`, `Canopy`)   | no-op (no SIF emission model for non-Lambertian)     |

# Units convention

`SurfaceSIF.SIF₀` is a Stokes vector (pol_n × nSpec) of hemispheric
isotropic emission with units **mW · m⁻² · cm⁻¹** (matching the source-
term unit decision in [solar_beam.jl](solar_beam.jl)). The factor of 2
absorbed at injection is `(1/π) · 2π`: `1/π` converts SIF irradiance to
Lambertian radiance, `2π` undoes the `weight = 0.5/π` azimuthal weighting
applied downstream in `postprocessing_vza!`.

=#

"""
    SurfaceSIF(; SIF₀=nothing, SIF760=nothing, mSIF=0,
               wavenumber_cm1=nothing,
               SIF755=nothing, slope=0, wavelength_nm=nothing) <: AbstractSource

User-facing surface fluorescence source. Carries an isotropic emission
spectrum and (optionally) is composed with other sources via `+`:

```julia
sources = SolarBeam() + SurfaceSIF(SIF₀ = sif_spec)
```

`SIF₀` is a Stokes vector (`pol_type.n` × `nSpec`) of surface-emitted
hemispheric irradiance (mW · m⁻² · cm⁻¹). The unpolarized component is
typically `SIF₀[1, :]`; higher Stokes components are zero unless the
canopy emission is polarized.

Alternatively, `SIF760` and `mSIF` define the retrievable unpolarized
spectral radiance

``L_SIF(ν) = SIF760 + mSIF * (ν - ν₇₆₀)``.

Both coefficients are appended to the Jacobian state vector, in that order.
`SIF760` is in mW m⁻² sr⁻¹ (cm⁻¹)⁻¹ and `mSIF` is in those units per cm⁻¹.
`wavenumber_cm1` must contain the model spectral grid. The legacy
`SIF755`/wavelength-domain `slope` interface remains supported.
For nonzero retrievable SIF, the source set must also contain a `SolarBeam` with an
explicit non-unit Fraunhofer `F₀`; unit-normalized radiances and physical SIF
must not be mixed.

When `SIF₀ === nothing`, [`prepare_source`](@ref) materialises a zero
matrix — the source is a no-op (useful as a placeholder).

# Fields
- `SIF₀ :: Union{Nothing, AbstractMatrix}`: surface emission Stokes
  vector or `nothing` (zero-default). Stored without an `eltype`
  constraint — the model's `FT` drives precision via `prepare_source`,
  matching [`SolarBeam`](@ref)'s FT-deferred design.
"""
struct SurfaceSIF <: AbstractSource
    SIF₀ :: Union{Nothing, AbstractMatrix}
    SIF755 :: Union{Nothing, Real}
    slope :: Real
    wavelength_nm :: Union{Nothing, AbstractVector}
    SIF760 :: Union{Nothing, Real}
    mSIF :: Real
    wavenumber_cm1 :: Union{Nothing, AbstractVector}
    ν_ref :: Real
end

function SurfaceSIF(; SIF₀=nothing, SIF760=nothing, SIF_ref=nothing, mSIF=0,
                    wavenumber_cm1=nothing, ν_ref=1e7/760,
                    SIF755=nothing, slope=0, wavelength_nm=nothing)
    SIF760 !== nothing && SIF_ref !== nothing && throw(ArgumentError(
        "Use SIF760, not both SIF760 and the deprecated SIF_ref alias."))
    SIF760 === nothing && (SIF760 = SIF_ref)
    SIF₀ !== nothing && (SIF755 !== nothing || SIF760 !== nothing) && throw(ArgumentError(
        "SurfaceSIF accepts either prescribed SIF₀ or retrievable SIF, not both."))
    SIF755 !== nothing && SIF760 !== nothing && throw(ArgumentError(
        "Use either legacy SIF755 or wavenumber-based SIF760, not both."))
    SIF755 !== nothing && wavelength_nm === nothing && throw(ArgumentError(
        "SurfaceSIF(SIF755=...) requires wavelength_nm on the model spectral grid."))
    SIF760 !== nothing && wavenumber_cm1 === nothing && throw(ArgumentError(
        "SurfaceSIF(SIF760=...) requires wavenumber_cm1 on the model spectral grid."))
    return SurfaceSIF(SIF₀, SIF755, slope, wavelength_nm,
                      SIF760, mSIF, wavenumber_cm1, ν_ref)
end

# Preserve the original positional constructor.
SurfaceSIF(SIF₀::Union{Nothing,AbstractMatrix}) = SurfaceSIF(SIF₀=SIF₀)

source_ad_mode(::SurfaceSIF) = AnalyticSourceJacobian()

Base.show(io::IO, s::SurfaceSIF) =
    s.SIF755 === nothing && s.SIF760 === nothing ?
        print(io, "SurfaceSIF(SIF₀=", s.SIF₀ === nothing ? "zeros" : summary(s.SIF₀), ")") :
    s.SIF760 !== nothing ?
        print(io, "SurfaceSIF(SIF760=", s.SIF760, ", mSIF=", s.mSIF,
              ", ν₇₆₀=", s.ν_ref, " cm⁻¹)") :
        print(io, "SurfaceSIF(SIF755=", s.SIF755, ", slope=", s.slope, ", λref=755 nm)")

"""
    PreparedSurfaceSIF{FT, AT} <: AbstractPreparedSource

Kernel-ready surface-fluorescence payload. `SIF₀` is materialised on the
model's array type at the right `(pol_type.n, nSpec)` shape and `FT`
precision.
"""
struct PreparedSurfaceSIF{FT<:AbstractFloat, AT<:AbstractMatrix, JT<:AbstractArray} <: AbstractPreparedSource
    SIF₀ :: AT
    SIḞ₀ :: JT
    n_parameters :: Int
end

source_ad_mode(::PreparedSurfaceSIF) = AnalyticSourceJacobian()

Base.show(io::IO, p::PreparedSurfaceSIF) =
    print(io, "PreparedSurfaceSIF(SIF₀=", summary(p.SIF₀), ")")

"""
    prepare_source(s::SurfaceSIF, FT, pol_n, nSpec, arr_type) -> PreparedSurfaceSIF

Resolve a [`SurfaceSIF`](@ref) into a kernel-ready
[`PreparedSurfaceSIF`](@ref). The default (`SIF₀ === nothing`)
materialises a zero matrix on the active architecture; a user-supplied
`SIF₀` is precision-converted and shape-checked.
"""
function prepare_source(s::SurfaceSIF, FT::Type{<:AbstractFloat},
                        pol_n::Integer, nSpec::Integer, arr_type)
    if s.SIF760 !== nothing
        length(s.wavenumber_cm1) == nSpec || throw(ArgumentError(
            "SurfaceSIF: wavenumber_cm1 length $(length(s.wavenumber_cm1)) does not match nSpec=$nSpec."))
        Δν = FT.(s.wavenumber_cm1) .- FT(s.ν_ref)
        L = FT(s.SIF760) .+ FT(s.mSIF) .* Δν
        SIF₀ = zeros(FT, pol_n, nSpec)
        SIḞ₀ = zeros(FT, pol_n, nSpec, 2)
        @views SIF₀[1, :] .= FT(π) .* L
        @views SIḞ₀[1, :, 1] .= FT(π)
        @views SIḞ₀[1, :, 2] .= FT(π) .* Δν
        dev, dotdev = arr_type(SIF₀), arr_type(SIḞ₀)
        return PreparedSurfaceSIF{FT, typeof(dev), typeof(dotdev)}(dev, dotdev, 2)
    elseif s.SIF755 !== nothing
        length(s.wavelength_nm) == nSpec || throw(ArgumentError(
            "SurfaceSIF: wavelength_nm length $(length(s.wavelength_nm)) does not match nSpec=$nSpec."))
        Δλ = FT.(s.wavelength_nm) .- FT(755)
        L = FT(s.SIF755) .+ FT(s.slope) .* Δλ
        # The surface kernel consumes hemispheric irradiance.  A Lambertian
        # radiance L corresponds to πL, hence this conversion at the seam.
        SIF₀ = zeros(FT, pol_n, nSpec)
        SIḞ₀ = zeros(FT, pol_n, nSpec, 2)
        @views SIF₀[1, :] .= FT(π) .* L
        @views SIḞ₀[1, :, 1] .= FT(π)
        @views SIḞ₀[1, :, 2] .= FT(π) .* Δλ
        dev = arr_type(SIF₀)
        dotdev = arr_type(SIḞ₀)
        return PreparedSurfaceSIF{FT, typeof(dev), typeof(dotdev)}(dev, dotdev, 2)
    elseif s.SIF₀ !== nothing
        size(s.SIF₀) == (pol_n, nSpec) || error(
            "SurfaceSIF: SIF₀ shape $(size(s.SIF₀)) does not match required " *
            "(pol_type.n, nSpec) = ($pol_n, $nSpec). Reshape your SIF spectrum " *
            "to match the model's polarization and spectral grid.")
        SIF₀_dev = arr_type(convert(Array{FT,2}, s.SIF₀))
        dotdev = arr_type(zeros(FT, pol_n, nSpec, 0))
        return PreparedSurfaceSIF{FT, typeof(SIF₀_dev), typeof(dotdev)}(SIF₀_dev, dotdev, 0)
    else
        SIF₀ = zeros(FT, pol_n, nSpec)
        dev = arr_type(SIF₀)
        dotdev = arr_type(zeros(FT, pol_n, nSpec, 0))
        return PreparedSurfaceSIF{FT, typeof(dev), typeof(dotdev)}(dev, dotdev, 0)
    end
end

surface_sif_parameter_count(::NoSource) = 0
surface_sif_parameter_count(s::SurfaceSIF) =
    s.SIF755 === nothing && s.SIF760 === nothing ? 0 : 2
surface_sif_parameter_count(::AbstractSource) = 0
surface_sif_parameter_count(s::SourceSet) = sum(surface_sif_parameter_count, s.sources)

function validate_sif_solar_spectrum(srcs::AbstractSource)
    surface_sif_parameter_count(srcs) == 0 && return nothing
    sif_nonzero = any(s -> s isa SurfaceSIF &&
                           ((s.SIF755 !== nothing && (!iszero(s.SIF755) || !iszero(s.slope))) ||
                            (s.SIF760 !== nothing && (!iszero(s.SIF760) || !iszero(s.mSIF)))),
                      srcs isa SourceSet ? srcs.sources : (srcs,))
    sif_nonzero || return nothing
    sources = srcs isa SourceSet ? srcs.sources : (srcs,)
    beam = findfirst(s -> s isa SolarBeam, sources)
    beam === nothing && throw(ArgumentError(
        "Nonzero SIF requires SolarBeam(F₀=...) with an explicitly supplied, non-unit Fraunhofer spectrum."))
    F₀ = sources[beam].F₀
    (F₀ === nothing || all(isone, @view(F₀[1, :]))) && throw(ArgumentError(
        "Nonzero SIF requires an explicitly supplied Fraunhofer F₀; the default/unit F₀=1 is not allowed."))
    return nothing
end

# ============================================================================
# Surface-side multiple dispatch — `surface_source_contribute!`
#
# Surfaces produce R/T matrices (in `create_surface_layer!`); sources
# contribute to the surface-layer source vector `j₀⁻` via this dispatch.
# Wired into `rt_run` as a single-line call after `create_surface_layer!`
# in a follow-up sub-phase.
# ============================================================================

"""
    surface_source_contribute!(prepared_sources::AbstractSource,
                                surface::AbstractSurfaceType,
                                surface_added_layer,
                                m::Int, pol_type, architecture)

Iterate `prepared_sources` and call the per-source per-surface kernel for
each member. Empty / NoSource → no-op; SourceSet → tuple unroll.
"""
surface_source_contribute!(::NoSource, ::AbstractSurfaceType,
                           _surface_added_layer, _m, _pol_type, _architecture) = nothing

function surface_source_contribute!(s::SourceSet, surface::AbstractSurfaceType,
                                    surface_added_layer, m::Integer, pol_type, architecture)
    @inbounds for src in s.sources
        surface_source_contribute!(src, surface, surface_added_layer, m, pol_type, architecture)
    end
    return nothing
end

# Default per-source surface contribution: no-op. Per-source-per-surface
# methods below provide the actual physics.
surface_source_contribute!(::AbstractPreparedSource, ::AbstractSurfaceType,
                           _surface_added_layer, _m, _pol_type, _architecture) = nothing

# Solar-beam contribution to the surface-layer source vector is currently
# computed inside `create_surface_layer!` from the resolved spectral F₀ for
# backward compatibility. A later sub-phase will move it out into a dedicated
# dispatch:
#   surface_source_contribute!(::PreparedSolarBeam, ::LambertianSurfaceScalar, …)
# The default (above) is a no-op; today's create_surface_layer! handles it.

"""
    surface_source_contribute!(prep::PreparedSurfaceSIF,
                                surface::Union{LambertianSurfaceScalar,
                                               LambertianSurfaceSpectrum,
                                               LambertianSurfaceLegendre,
                                               LambertianSurfaceSpline},
                                surface_added_layer, m, pol_type, architecture)

Inject hemispheric SIF emission into the Lambertian surface added-layer's
upwelling source vector at the m=0 Fourier moment. Bit-equal to today's
`inject_surface_SIF!(brdf, surface_added_layer, m, pol_type, SIF₀, arch)`.

The factor of 2 is `(1/π) · 2π` — `1/π` converts SIF irradiance to
Lambertian radiance, `2π` undoes the `weight = 0.5/π` azimuthal weighting
applied downstream in `postprocessing_vza!` so the isotropic SIF
contribution survives unweighted.
"""
function surface_source_contribute!(prep::PreparedSurfaceSIF,
        ::Union{LambertianSurfaceScalar, LambertianSurfaceSpectrum,
                LambertianSurfaceLegendre, LambertianSurfaceSpline},
        surface_added_layer,
        m::Integer, pol_type, architecture)
    m == 0 || return nothing
    iszero(prep.SIF₀) && return nothing
    FT = eltype(surface_added_layer.j₀⁻)
    arr_type = array_type(architecture)
    Nquad = size(surface_added_layer.j₀⁻, 1) ÷ pol_type.n
    surface_added_layer.j₀⁻[:, 1, :] .+= FT(2) .* arr_type(repeat(FT.(prep.SIF₀), Nquad))
    return nothing
end


function surface_source_contribute_lin!(::NoSource, ::AbstractSurfaceType,
        _layer, _layer_lin, _m, _pol, _arch, _sif_range)
    return nothing
end

function surface_source_contribute_lin!(sources::SourceSet, surface::AbstractSurfaceType,
        layer, layer_lin, m, pol, arch, indices)
    offset = 0
    for src in sources.sources
        n = src isa PreparedSurfaceSIF ? src.n_parameters : 0
        local_indices = n == 0 ? (1:0) : (first(indices) + offset):(first(indices) + offset + n - 1)
        surface_source_contribute_lin!(src, surface, layer, layer_lin, m, pol, arch, local_indices)
        offset += n
    end
    return nothing
end

surface_source_contribute_lin!(::AbstractPreparedSource, ::AbstractSurfaceType,
        _layer, _layer_lin, _m, _pol, _arch, _indices) = nothing

function surface_source_contribute_lin!(prep::PreparedSurfaceSIF,
        ::Union{LambertianSurfaceScalar, LambertianSurfaceSpectrum,
                LambertianSurfaceLegendre, LambertianSurfaceSpline},
        layer, layer_lin, m::Integer, pol_type, architecture, indices)
    m == 0 || return nothing
    prep.n_parameters == 0 && return nothing
    FT = eltype(layer.j₀⁻)
    Nquad = size(layer.j₀⁻, 1) ÷ pol_type.n
    arr_type = array_type(architecture)
    for (k, p) in enumerate(indices)
        @views layer_lin.ap_J̇₀⁻[:, 1, :, p] .+=
            FT(2) .* arr_type(repeat(FT.(prep.SIḞ₀[:, :, k]), Nquad))
    end
    return nothing
end

# Predicate for "any source set carries a SurfaceSIF" (mirrors `has_solar_beam`).
has_surface_sif(::NoSource) = false
has_surface_sif(::PreparedSurfaceSIF) = true
has_surface_sif(::AbstractPreparedSource) = false

function has_surface_sif(s::SourceSet)
    @inbounds for src in s.sources
        src isa PreparedSurfaceSIF && return true
    end
    return false
end
