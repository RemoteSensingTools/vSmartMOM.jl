#=

`SurfaceSIF` — first-class surface-emission source for the v0.6 source-term
refactor. Represents Solar-Induced (chlorophyll) Fluorescence as an
isotropic emission at the lower boundary that adds to the surface layer's
upwelling source vector `j₀⁻`.

This is the "surfaces produce R/T only; sources contribute to surface
j₀± via dispatch" pattern (Christian's design directive). Today's
`inject_surface_SIF!` shotgun-injects SIF₀ from `RS_type.SIF₀` into the
Lambertian surface added-layer. Phase 5 replaces that with a clean
double-dispatch:

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
    SurfaceSIF(; SIF₀ = nothing) <: AbstractSource

User-facing surface fluorescence source. Carries an isotropic emission
spectrum and (optionally) is composed with other sources via `+`:

```julia
sources = SolarBeam() + SurfaceSIF(SIF₀ = sif_spec)
```

`SIF₀` is a Stokes vector (`pol_type.n` × `nSpec`) of surface-emitted
hemispheric irradiance (mW · m⁻² · cm⁻¹). The unpolarized component is
typically `SIF₀[1, :]`; higher Stokes components are zero unless the
canopy emission is polarized.

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
end

SurfaceSIF(; SIF₀=nothing) = SurfaceSIF(SIF₀)

source_ad_mode(::SurfaceSIF) = AnalyticSourceJacobian()

Base.show(io::IO, s::SurfaceSIF) =
    print(io, "SurfaceSIF(SIF₀=", s.SIF₀ === nothing ? "zeros" : summary(s.SIF₀), ")")

"""
    PreparedSurfaceSIF{FT, AT} <: AbstractPreparedSource

Kernel-ready surface-fluorescence payload. `SIF₀` is materialised on the
model's array type at the right `(pol_type.n, nSpec)` shape and `FT`
precision.
"""
struct PreparedSurfaceSIF{FT<:AbstractFloat, AT<:AbstractMatrix} <: AbstractPreparedSource
    SIF₀ :: AT
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
    if s.SIF₀ === nothing
        SIF₀ = zeros(FT, pol_n, nSpec)
        return PreparedSurfaceSIF{FT, typeof(arr_type(SIF₀))}(arr_type(SIF₀))
    else
        size(s.SIF₀) == (pol_n, nSpec) || error(
            "SurfaceSIF: SIF₀ shape $(size(s.SIF₀)) does not match required " *
            "(pol_type.n, nSpec) = ($pol_n, $nSpec). Reshape your SIF spectrum " *
            "to match the model's polarization and spectral grid.")
        SIF₀_dev = arr_type(convert(Array{FT,2}, s.SIF₀))
        return PreparedSurfaceSIF{FT, typeof(SIF₀_dev)}(SIF₀_dev)
    end
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
# computed inside `create_surface_layer!` for backward compatibility. A
# later sub-phase will move it out into a dedicated dispatch:
#   surface_source_contribute!(::PreparedSolarBeam, ::LambertianSurfaceScalar, …)
# The default (above) is a no-op; today's create_surface_layer! handles it.

"""
    surface_source_contribute!(prep::PreparedSurfaceSIF,
                                surface::Union{LambertianSurfaceScalar,
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
        ::Union{LambertianSurfaceScalar, LambertianSurfaceLegendre, LambertianSurfaceSpline},
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
