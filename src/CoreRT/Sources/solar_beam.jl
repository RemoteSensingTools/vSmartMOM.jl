#=

`SolarBeam` — first-class solar / stellar / lunar collimated source for the
v0.6 source-term refactor. Phase 2 lands the user-facing type plus
`prepare_source` so users can write

```julia
R, T = rt_run(model; sources = SolarBeam(F₀ = solar_spec))
```

The SFI math itself is **unchanged in Phase 2**: F₀ still flows through
`RS_type.F₀` via the existing `get_elem_rt_SFI!` kernel. Phase 3 relocates
the analytic Sanghavi 2014 tangent into `source_tangent!`, and Phase 5
removes `RS_type.F₀` ownership entirely. See the user-facing source-term
guide at `docs/src/pages/extending/sources.md` for the design summary.

=#

"""
    SolarBeam(; F₀ = nothing, sza = nothing)

User-facing collimated direct-beam source. Carries the solar Stokes
irradiance spectrum and an optional zenith angle — **without committing
to a floating-point precision**, so the same `SolarBeam` works in any
model regardless of `FT`. [`prepare_source`](@ref) does the
precision/architecture conversion against the model's `FT` at solve time
(and that's also the AD seam — see Pillar D in the v0.6 plan).

`F₀`, when supplied, must be a `(pol_type.n, nSpec)` matrix:
- `F₀[1, :]` is unpolarized Stokes-I irradiance per spectral point.
- `F₀[2:end, :]` is incident polarization (Q/U/V), zero by default.

When `F₀ === nothing`, `prepare_source` materialises a unit Stokes-I
vector at the model's `FT` — bit-equal to today's `RS_type.F₀ = ones`
default.

# Fields
- `F₀ :: Union{Nothing, AbstractMatrix}`: solar irradiance Stokes vector
  or `nothing` for the unit default. Stored without an `eltype`
  constraint so users can pass a `Vector{Float32}` matrix into a
  `Float64` model (or vice versa) — the conversion happens once in
  `prepare_source`.
- `sza :: Union{Nothing, Real}`: **advisory only in v0.6.** vSmartMOM
  currently always reads SZA from `RTModel.obs_geom.sza` (which is fixed
  at model construction by `parameters_from_yaml(...).sza`). The
  `SolarBeam.sza` field is reserved for a future per-source-geometry
  override (rt_run will rebuild quad_points from this when set), but
  setting it today is a no-op — the model's geometry wins. To change
  SZA today, set `params.sza` before `model_from_parameters(params)`.

# AD mode
[`source_ad_mode`](@ref) returns [`AnalyticSourceJacobian`](@ref); Phase 3
provides the analytic `source_tangent!` body (relocated from
`get_elem_rt_SFI_fused!`).

# Examples

```julia
sb = SolarBeam()                                   # default unit Stokes I
sb = SolarBeam(; sza = 35.0)                       # override SZA
sb = SolarBeam(; F₀ = my_solar_irradiance)         # custom spectrum
```
"""
struct SolarBeam <: AbstractSource
    F₀  :: Union{Nothing, AbstractMatrix}
    sza :: Union{Nothing, Real}
end

SolarBeam(; F₀=nothing, sza=nothing) = SolarBeam(F₀, sza)

source_ad_mode(::SolarBeam) = AnalyticSourceJacobian()

Base.show(io::IO, sb::SolarBeam) =
    print(io, "SolarBeam(F₀=", sb.F₀ === nothing ? "default" : summary(sb.F₀),
          ", sza=", sb.sza === nothing ? "model.obs_geom" : sb.sza, ")")

"""
    PreparedSolarBeam{FT, AT<:AbstractMatrix} <: AbstractPreparedSource

Kernel-ready solar-beam payload. `F₀` is materialised on the model's
array type (CPU `Array` or `CuArray`) at the right shape
`(pol_type.n, nSpec)`. The geometry indices live in the
[`QuadPoints`](@ref) struct already on the model — `PreparedSolarBeam`
holds neither μ₀ nor iμ₀ to avoid duplicating mutable state.

# Fields
- `F₀ :: AT`: solar irradiance Stokes vector on the active architecture.
"""
struct PreparedSolarBeam{FT<:AbstractFloat, AT<:AbstractMatrix} <: AbstractPreparedSource
    F₀ :: AT
end

source_ad_mode(::PreparedSolarBeam) = AnalyticSourceJacobian()

Base.show(io::IO, p::PreparedSolarBeam) =
    print(io, "PreparedSolarBeam(F₀=", summary(p.F₀), ")")

"""
    prepare_source(sb::SolarBeam, FT::Type, pol_n::Int, nSpec::Int, arr_type) -> PreparedSolarBeam

Resolve a [`SolarBeam`](@ref) into a kernel-ready
[`PreparedSolarBeam`](@ref).

The default (`F₀ === nothing`) materialises a unit Stokes-I irradiance
matching today's `rt_run` allocation — `F₀[1, :] .= 1`, all other Stokes
components zero. A user-provided `F₀` is reshaped/converted to the
requested `(pol_n, nSpec)` shape and `FT` precision; when the user-side
shape doesn't match, this is an error rather than a silent broadcast.

This is the AD seam (constraint 3): everything above this call can use
ForwardDiff `Dual` numbers; the returned `PreparedSolarBeam.F₀` is plain
`FT`. A future `prepare_source_with_tangent` will mirror this signature
and return the source-parameter Jacobian alongside.
"""
function prepare_source(sb::SolarBeam, FT::Type{<:AbstractFloat},
                        pol_n::Integer, nSpec::Integer, arr_type)
    if sb.F₀ === nothing
        F₀ = zeros(FT, pol_n, nSpec)
        @inbounds @views F₀[1, :] .= one(FT)
        return PreparedSolarBeam{FT, typeof(arr_type(F₀))}(arr_type(F₀))
    else
        size(sb.F₀) == (pol_n, nSpec) || error(
            "SolarBeam: F₀ shape $(size(sb.F₀)) does not match required " *
            "(pol_type.n, nSpec) = ($pol_n, $nSpec). Reshape your spectrum " *
            "to match the model's polarization and spectral grid.")
        F₀_dev = arr_type(convert(Array{FT,2}, sb.F₀))
        return PreparedSolarBeam{FT, typeof(F₀_dev)}(F₀_dev)
    end
end

"""
    prepare_sources(srcs::AbstractSource, FT, pol_n, nSpec, arr_type) -> AbstractPreparedSource

Walk the [`SourceSet`](@ref) tuple and call [`prepare_source`](@ref) on
each member, returning a structurally-parallel [`SourceSet`](@ref) of
[`AbstractPreparedSource`](@ref). [`NoSource`](@ref) round-trips
unchanged. A bare single source is wrapped in a one-element `SourceSet`
so the rest of the orchestration only has to iterate.
"""
prepare_sources(::NoSource, ::Type{<:AbstractFloat}, ::Integer, ::Integer, _arr_type) =
    NoSource()

function prepare_sources(s::SourceSet, FT::Type{<:AbstractFloat},
                         pol_n::Integer, nSpec::Integer, arr_type)
    return SourceSet(map(src -> prepare_source(src, FT, pol_n, nSpec, arr_type), s.sources))
end

function prepare_sources(s::AbstractSource, FT::Type{<:AbstractFloat},
                         pol_n::Integer, nSpec::Integer, arr_type)
    return SourceSet((prepare_source(s, FT, pol_n, nSpec, arr_type),))
end

"""
    extract_solar_F₀(prepared, FT, pol_n, nSpec, arr_type) -> AbstractMatrix

Return the `F₀` matrix carried by the first
[`PreparedSolarBeam`](@ref) in `prepared` (a [`PreparedSolarBeam`](@ref),
a `SourceSet` of prepared sources, or [`NoSource`](@ref)). When no
`PreparedSolarBeam` is present, return the unit Stokes-I default — the
same matrix today's `rt_run` would have allocated when
`size(RS_type.F₀) != (pol_type.n, nSpec)`.

This helper is the Phase-2 bridge between the new source vocabulary and
the legacy `RS_type.F₀` channel still consumed by the kernels. Phase 5
removes the bridge by giving the kernels direct access to
`PreparedSolarBeam.F₀`.
"""
function extract_solar_F₀(::NoSource, FT::Type{<:AbstractFloat},
                          pol_n::Integer, nSpec::Integer, arr_type)
    F₀ = zeros(FT, pol_n, nSpec)
    return arr_type(F₀)
end

function extract_solar_F₀(p::PreparedSolarBeam, ::Type{<:AbstractFloat},
                          ::Integer, ::Integer, _arr_type)
    return p.F₀
end

function extract_solar_F₀(s::SourceSet, FT::Type{<:AbstractFloat},
                          pol_n::Integer, nSpec::Integer, arr_type)
    @inbounds for src in s.sources
        if src isa PreparedSolarBeam
            return src.F₀
        end
    end
    # User supplied an explicit `SourceSet` that does not contain a
    # PreparedSolarBeam — they meant "no solar contribution". Return zeros
    # so the SFI kernel produces zero source vectors. Today's bit-equal
    # default (unit Stokes-I) only applies when the user did NOT specify
    # `sources=` at all, in which case `model.sources` defaults to a
    # `SolarBeam()` and this branch is never reached.
    F₀ = zeros(FT, pol_n, nSpec)
    return arr_type(F₀)
end

"""
    has_solar_beam(prepared) -> Bool

Predicate used by `rt_run` to derive the legacy `SFI::Bool` flag from a
prepared source set. Returns `true` whenever a [`PreparedSolarBeam`](@ref)
is present, mirroring today's "SFI is on by default" behaviour.
"""
has_solar_beam(::NoSource) = false
has_solar_beam(p::PreparedSolarBeam) = true
has_solar_beam(::AbstractPreparedSource) = false

function has_solar_beam(s::SourceSet)
    @inbounds for src in s.sources
        src isa PreparedSolarBeam && return true
    end
    return false
end

# ============================================================================
# BlackbodySource — Planck-derived collimated source (Phase 4)
#
# Units convention (Phase 4 source-unit decision):
#   F₀  [mW · m⁻² · cm⁻¹]  spectral irradiance (Stokes vector by component)
#   B(ν,T) returned by `vSmartMOM.SolarModel.planck_spectrum_wn` is the
#   Planck SPECTRAL RADIANCE in mW · m⁻² · sr⁻¹ · cm⁻¹. Multiplying by
#   `factor=π` converts radiance from a Lambertian disk to incident
#   hemispheric irradiance (mW · m⁻² · cm⁻¹) — directly comparable to
#   `SolarBeam` F₀ = 1 (which we adopt as 1 mW · m⁻² · cm⁻¹ for unit
#   normalization). All sources in a `SourceSet` should share these units;
#   `R = rt_run(model; sources=...)` is then in mW · m⁻² · sr⁻¹ · cm⁻¹.
#
# `BlackbodySource` is a constructor function (not a struct) — it
# materialises the Planck spectrum at user-construction time and returns a
# `SolarBeam` with the right `F₀`. The kernel sees a `SolarBeam`; the
# "this is thermal" semantics live entirely above the `prepare_source`
# seam. Once we add full ThermalEmission (atmospheric volume Planck
# integral) in a later phase, the two concepts will be cleanly separated:
#   - `BlackbodySource(T, spec_band)`  → boundary/external Planck beam
#   - `ThermalEmission(T_layers)`      → in-atmosphere volume emission
# ============================================================================

"""
    BlackbodySource(T, spec_band; pol_component=1, pol_n=3, factor=π, scale=1) -> SolarBeam

Construct a [`SolarBeam`](@ref) whose `F₀` is the Planck spectrum at
temperature `T` (K) on spectral grid `spec_band` (cm⁻¹). Useful for
modelling thermal-IR scenes such as a Carbon-I-like 2-2.4 µm setup with a
hot lab source illuminating CO₂/CH₄/H₂O absorption.

`F₀[pol_component, :] = factor · scale · B(ν, T)` where
`B(ν,T)` is the Planck radiance from `planck_spectrum_wn` in
mW · m⁻² · sr⁻¹ · cm⁻¹. The default `factor = π` converts radiance to
hemispheric irradiance for a Lambertian-disk source — set `factor = 1`
for a head-on collimated lab beam where `F₀` is the radiance directly.

# Arguments
- `T :: Real`: blackbody temperature in K.
- `spec_band :: AbstractVector{<:Real}`: spectral grid in cm⁻¹.
- `pol_component :: Integer = 1`: which Stokes component carries the
  source (1=I, 2=Q, 3=U, 4=V); only `1` is physically meaningful for an
  unpolarized blackbody.
- `pol_n :: Integer = 3`: number of Stokes components in the model
  (1=Stokes_I, 3=Stokes_IQU, 4=Stokes_IQUV). Match the model's
  polarization to avoid `prepare_source` shape errors.
- `factor :: Real = π`: geometric factor (π for Lambertian-disk → hemisphere irradiance).
- `scale :: Real = 1`: extra normalization multiplier.

# Units
`F₀` has units **mW · m⁻² · cm⁻¹** (irradiance). All sources in a
`SourceSet` must agree on units, so a default `SolarBeam(F₀=ones(...))`
should be interpreted as `1 mW · m⁻² · cm⁻¹`. The radiance returned by
`rt_run` is in **mW · m⁻² · sr⁻¹ · cm⁻¹**.

# Example: Carbon-I-like 2-2.4 µm with a 1500 K source
```julia
spec_band = collect(4167:0.1:5000)            # cm⁻¹ for 2-2.4 µm
sources   = BlackbodySource(1500, spec_band)  # SolarBeam with Planck F₀
R, T      = rt_run(model; sources = sources)  # R in mW·m⁻²·sr⁻¹·cm⁻¹
```
"""
function BlackbodySource(T::Real, spec_band::AbstractVector{<:Real};
                         pol_component::Integer = 1,
                         pol_n::Integer = 3,
                         factor::Real = π,
                         scale::Real = 1.0)
    pol_component <= pol_n ||
        error("BlackbodySource: pol_component=$pol_component exceeds pol_n=$pol_n.")
    ν_grid = collect(float.(spec_band))
    radiance = vSmartMOM.SolarModel.planck_spectrum_wn(T, ν_grid)  # mW/m²/sr/cm⁻¹
    nSpec = length(ν_grid)
    F₀ = zeros(eltype(radiance), pol_n, nSpec)
    coeff = float(factor) * float(scale)
    @inbounds @views F₀[pol_component, :] .= coeff .* radiance
    return SolarBeam(F₀ = F₀)
end

# ============================================================================
# Phase 3 — source_tangent! seam
#
# `source_tangent!(::PreparedSolarBeam, ...)` is the named hand-written
# linearization entry point for the solar beam. It satisfies constraint 2
# (hand-written linearization keeps working — bit-equal numerics) and
# constraint 3 (future AD ↔ hand-written boundary is named).
#
# In Phase 3 the body is a thin wrapper around the existing fused kernel
# `CoreRT.CoreKernel.get_elem_rt_SFI_fused!` defined in
# `src/CoreRT/CoreKernel/elemental_lin.jl`. Math is unchanged; the
# Sanghavi 2014 App. C tangents (∂J/∂τ, ∂J/∂ϖ, ∂J/∂Z) plus the Bug-22
# per-parameter beam-attenuation derivative remain in that kernel.
#
# A future sub-phase can either:
#   (a) relocate the kernel body into this file (literal source-locality
#       refactor — zero numeric change), or
#   (b) replace the wrapper with a `LayerContext` / `LayerLinContext`-based
#       dispatch (Pillar A in the v0.6 plan) once the context structs land.
# Either evolution preserves bit-equality.
# ============================================================================

"""
    source_tangent!(prep::PreparedSolarBeam,
                    j₀⁺, j₀⁻, J̇₀⁺, J̇₀⁻, ap_J̇₀⁺, ap_J̇₀⁻,
                    ϖ_λ, dτ_λ, τ_sum, τ̇_sum, Z⁻⁺, Z⁺⁺,
                    dτ̇, ϖ̇, Ż⁻⁺, Ż⁺⁺,
                    qp_μN, ndoubl, wct02, nStokes, I₀, iμ0, D, nparams,
                    architecture)

Hand-written analytic source tangent for a [`PreparedSolarBeam`](@ref).
Writes forward source vectors (`j₀⁺`, `j₀⁻`), 3-core derivatives
(`J̇₀±[i,1,n,1:3]` w.r.t. `τ, ϖ, Z`), and per-physical-parameter chain-rule
slabs (`ap_J̇₀±[i,1,n,iparam]`) — all in one fused GPU pass — by delegating
to `CoreRT.CoreKernel.get_elem_rt_SFI_fused!`. Bit-equal to today's
in-line dispatch in `elemental_lin.jl`.

Implements Sanghavi 2014 App. C Eqs. (C.8)–(C.10) plus the Bug-22
beam-attenuation chain-rule fix `J̇₀± += J₀± · (-τ̇_sum/μ₀)`.

# AD seam
This is the **below-the-seam** path: arguments are plain `FT`. Above the
seam, [`prepare_source`](@ref) materialises the F₀ from a user-facing
[`SolarBeam`](@ref) (potentially under ForwardDiff `Dual` numbers in
v0.7+). The seam is named here so future AD work lands without renaming.

# Mode trait
[`source_ad_mode`](@ref) returns [`AnalyticSourceJacobian`](@ref) — this
function is the authoritative path for solar-beam Jacobians.
"""
function source_tangent!(prep::PreparedSolarBeam,
                         j₀⁺, j₀⁻, J̇₀⁺, J̇₀⁻, ap_J̇₀⁺, ap_J̇₀⁻,
                         ϖ_λ, dτ_λ, τ_sum, τ̇_sum,
                         Z⁻⁺, Z⁺⁺,
                         dτ̇, ϖ̇, Ż⁻⁺, Ż⁺⁺,
                         qp_μN, ndoubl, wct02, nStokes, I₀, iμ0, D, nparams,
                         architecture)
    # Reuse the existing fused SFI + chain-rule kernel. `get_elem_rt_SFI_fused!`
    # is forward-referenced here — its definition lives in
    # `src/CoreRT/CoreKernel/elemental_lin.jl` which is included after this file.
    device = devi(architecture)
    kernel! = get_elem_rt_SFI_fused!(device)
    kernel!(j₀⁺, j₀⁻,
            J̇₀⁺, J̇₀⁻,
            ap_J̇₀⁺, ap_J̇₀⁻,
            ϖ_λ, dτ_λ, τ_sum, τ̇_sum,
            Z⁻⁺, Z⁺⁺, prep.F₀,
            dτ̇, ϖ̇, Ż⁻⁺, Ż⁺⁺,
            qp_μN, ndoubl, wct02,
            nStokes, I₀, iμ0, D, nparams,
            ndrange = size(j₀⁺))
    return nothing
end

# ============================================================================
# Multiple-dispatch entry points for NoSource and SourceSet (Phase 3.5)
#
# These dispatchers replace `if SFI` branching across the elastic kernels.
# `NoSource` and an empty `SourceSet` route to a no-op; a populated
# `SourceSet` iterates over its tuple at compile time and dispatches each
# member to the correct concrete `contribute!` / `source_tangent!`.
#
# The signature here matches the per-source forms above so a single call
# site `contribute!(prepared_sources, ...)` works regardless of whether
# `prepared_sources` is a `NoSource`, a single `PreparedSolarBeam`, or a
# `SourceSet` of mixed prepared sources.
# ============================================================================

# Forward `contribute!` — to be defined per concrete prepared source. SolarBeam
# is added in a follow-up sub-phase that wires the elemental forward kernel
# through this dispatch. NoSource and SourceSet dispatchers are universal:
contribute!(::NoSource, args...; kwargs...) = nothing

function contribute!(s::SourceSet, args...; kwargs...)
    @inbounds for src in s.sources
        contribute!(src, args...; kwargs...)
    end
    return nothing
end

# Linearized `source_tangent!` — NoSource and SourceSet dispatchers.
source_tangent!(::NoSource, args...; kwargs...) = nothing

function source_tangent!(s::SourceSet, args...; kwargs...)
    @inbounds for src in s.sources
        source_tangent!(src, args...; kwargs...)
    end
    return nothing
end
