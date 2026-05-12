#=

This file declares the v0.6 source-term type vocabulary.

The MOM solver is mathematically affine: optical properties define the
operator A, source terms define the additive RHS b. Adding/doubling/interaction
propagate (A, b) generically without knowing whether b came from a solar
beam, thermal emission, surface fluorescence, or a lidar pulse.

`AbstractSource` makes that affine structure explicit at the type level.
A source is a user-facing description of an additive RHS contribution; a
`PreparedSource` (declared by Phase 2 implementations) is its kernel-ready,
isbits-friendly counterpart, mirroring the `Aerosol` ↔ `AerosolOptics`
two-tier idiom already in vSmartMOM.

# Composition contract

```
SolarBeam + SurfaceSIF             # ⇒ SourceSet((SolarBeam, SurfaceSIF))
SourceSet + ThermalEmission        # ⇒ SourceSet((..., ThermalEmission))
NoSource  + s                      # ⇒ s   (NoSource is the additive identity)
```

Internally, a `SourceSet` carries a concrete `Tuple` of sources so the hot
loop can iterate them without dynamic dispatch (CPU and GPU).

# Source contract (declared in Phase 2)

Each concrete source type implements:

```
prepare_source(::AbstractSource, model, iBand)  -> AbstractPreparedSource
contribute!(j₀⁺, j₀⁻, ::AbstractPreparedSource, layer_ctx)
source_tangent!(j̇₀⁺, j̇₀⁻, ::AbstractPreparedSource, layer_ctx, layer_lin_ctx)   # if AD-mode = Analytic
```

`prepare_source` is the seam between user-parameter space (where
ForwardDiff or analytic upstream linearization can live) and kernel space
(plain `FT`, hand-written `contribute!` / `source_tangent!`). The seam is
named here in v0.6; the AD half lands in v0.7+ without renaming the API.

# AD differentiation modes

A future `prepare_source_with_tangent` will mirror `prepare_source` for
parameter-space AD. Each source declares its differentiation mode via the
`source_ad_mode` trait so dispatch can route to analytic tangents,
ForwardDiff over `prepare_source`, or skip Jacobian computation entirely.

=#

"""
    AbstractSource

Top-level abstract type for v0.6 first-class source terms (solar, thermal,
surface fluorescence, lidar, …). Distinct from the legacy
`AbstractSourceType` (which has unused subtypes `DNI` / `SFI` and is the
type slot reserved by the v0.5 dispatch design for future thermal RT —
left untouched in v0.6 to avoid breaking external code).

A concrete source must satisfy the contract documented in
`src/CoreRT/Sources/types.jl`: user-facing configuration that round-trips
through `prepare_source` into a [`AbstractPreparedSource`](@ref) before
reaching kernel code.
"""
abstract type AbstractSource end

"""
    AbstractPreparedSource

Kernel-ready, isbits-friendly counterpart of an [`AbstractSource`](@ref).
Concrete prepared sources (e.g. `PreparedSolarBeam`) hold values already
materialised on the active architecture — μ₀-index for the solar stream,
broadcast-shaped F₀ on the device array type, attenuation buffers, etc.

Prepared sources flow through the elemental hot loop via `contribute!`
(forward) and `source_tangent!` (linearized). They do not own their
buffers' lifetime; allocation is pinned at model build time.
"""
abstract type AbstractPreparedSource end

"""
    NoSource()

Additive identity for source composition. Useful as an explicit dispatch
target for "no active source" and as the default when `rt_run` is called
without a `sources=` kwarg in legacy paths.

```julia
NoSource() + s == s   # for any AbstractSource s
```
"""
struct NoSource <: AbstractSource end

"""
    SourceSet(sources::Tuple) <: AbstractSource

Type-stable composite of source contributions. The internal tuple holds
one `AbstractSource` per affine RHS contributor; the elemental hot loop
unrolls over it at compile time.

Construct via the `+` operator, never directly:

```julia
src = SolarBeam(F₀=F₀_spec, sza=35) + SurfaceSIF(strength=sif)
# ⇒ SourceSet((SolarBeam(...), SurfaceSIF(...)))
```

Iteration, length, and indexing forward to the underlying tuple.
"""
struct SourceSet{S<:Tuple} <: AbstractSource
    sources::S
end

SourceSet() = SourceSet{Tuple{}}(())
SourceSet(srcs::AbstractSource...) = SourceSet(srcs)

# Composition (`+`) ===========================================================
#
# The composition rules are explicit per pair to make NoSource the additive
# identity and to flatten SourceSet inputs (no nested SourceSets). Methods
# are written so the most-specific match wins under Julia's dispatch — see
# the comment-block at the end of this file for the dispatch table.

# NoSource × NoSource → NoSource (identity preserved)
Base.:+(::NoSource, ::NoSource) = NoSource()

# NoSource × concrete (and concrete × NoSource): drop the NoSource side
Base.:+(::NoSource, b::AbstractSource) = b
Base.:+(a::AbstractSource, ::NoSource) = a

# Disambiguation: NoSource × SourceSet — pick the SourceSet
Base.:+(::NoSource, b::SourceSet) = b
Base.:+(a::SourceSet, ::NoSource) = a

# SourceSet × SourceSet — flatten without nesting
Base.:+(a::SourceSet, b::SourceSet) = SourceSet((a.sources..., b.sources...))

# SourceSet × concrete and concrete × SourceSet — flatten
Base.:+(a::SourceSet, b::AbstractSource) = SourceSet((a.sources..., b))
Base.:+(a::AbstractSource, b::SourceSet) = SourceSet((a, b.sources...))

# concrete × concrete — wrap into a fresh SourceSet
Base.:+(a::AbstractSource, b::AbstractSource) = SourceSet((a, b))

# Tuple-like access ==========================================================

Base.iterate(s::SourceSet) = iterate(s.sources)
Base.iterate(s::SourceSet, state) = iterate(s.sources, state)
Base.length(s::SourceSet) = length(s.sources)
Base.eltype(::Type{SourceSet{S}}) where {S} = eltype(S)
Base.getindex(s::SourceSet, i::Integer) = s.sources[i]
Base.firstindex(::SourceSet) = 1
Base.lastindex(s::SourceSet) = length(s)

Base.show(io::IO, ::NoSource) = print(io, "NoSource()")
function Base.show(io::IO, s::SourceSet)
    print(io, "SourceSet(")
    for (i, src) in enumerate(s.sources)
        i > 1 && print(io, ", ")
        show(io, src)
    end
    print(io, ")")
end

# AD differentiation modes ===================================================

"""
    AbstractSourceADMode

Trait hierarchy describing how a source's parameters participate in
linearized RT. Each source declares its mode via [`source_ad_mode`](@ref);
the linearization driver dispatches to the appropriate path.

Concrete modes:
- [`AnalyticSourceJacobian`](@ref): hand-written `source_tangent!` is the
  authoritative path. Today's solar SFI tangents (Sanghavi 2014 App. C)
  belong here.
- [`ForwardDiffSourceJacobian`](@ref): declared but not implemented in
  v0.6 — reserved for v0.7+ source-parameter Jacobians via ForwardDiff
  through `prepare_source`.
- [`NoSourceJacobian`](@ref): the source contributes no parameters to the
  Jacobian (e.g. `NoSource`).
"""
abstract type AbstractSourceADMode end

"""
    AnalyticSourceJacobian()

Source ships a hand-written `source_tangent!` that fills the analytic
core derivatives `(∂j/∂τ, ∂j/∂ϖ, ∂j/∂Z)` plus any source-parameter
column claims in the linearized layer's `ap_J̇₀±` slabs.
"""
struct AnalyticSourceJacobian <: AbstractSourceADMode end

"""
    ForwardDiffSourceJacobian()

Reserved for v0.7+. Source-parameter Jacobians come from ForwardDiff
through `prepare_source_with_tangent`; optical-parameter tangents
remain analytic via `source_tangent!`. v0.6 declares this trait so
future code lands without renaming the API.
"""
struct ForwardDiffSourceJacobian <: AbstractSourceADMode end

"""
    NoSourceJacobian()

Source contributes no parameters to the Jacobian and has no
`source_tangent!` body. `NoSource` uses this mode.
"""
struct NoSourceJacobian <: AbstractSourceADMode end

"""
    source_ad_mode(::AbstractSource) -> AbstractSourceADMode

Trait. Concrete sources override this; the default is
[`AnalyticSourceJacobian`](@ref) so adding a new analytic source requires
no extra ceremony.
"""
source_ad_mode(::AbstractSource) = AnalyticSourceJacobian()
source_ad_mode(::NoSource) = NoSourceJacobian()
# A SourceSet's mode is informational only — the linearization driver
# iterates the tuple and dispatches per-source via `source_ad_mode` on
# each member. Default to Analytic so callers asking the set as a whole
# get the same behaviour as a homogeneous analytic set.
source_ad_mode(::SourceSet) = AnalyticSourceJacobian()

# ============================================================================
# v0.7 Phase A.2a — Per-source slot traits
# ============================================================================
#
# `source_key(::AbstractPreparedSource) :: Union{Nothing, Symbol}`
#
#   The NamedTuple key under which this source's j₀± / J₀± live in the
#   `AddedLayer.j₀_by_src` / `CompositeLayer.J₀_by_src` containers.
#
#   - `:solar`      → SolarBeam (special: lives in legacy `j₀⁺/j₀⁻` fields,
#                     not in the by_src NamedTuple, so this returns `nothing`
#                     to skip slot allocation. The legacy field handles it.)
#   - `:thermal`    → ThermalEmission (volume Planck)
#   - `:surface_sif` → SurfaceSIF (lives in legacy slot too — surface
#                     injection only, no doubling — so returns `nothing`)
#   - `nothing`     → no slot needed (legacy path / surface-only sources)
#
# `source_expk_init(::AbstractPreparedSource, dτ_λ, μ₀, arr_type) -> AbstractArray{FT,1}`
#
#   The initial per-spectral expk vector for this source's doubling
#   recurrence. Squared each doubling step in the noRS doubling loop.
#   For thermal this is `ones(nSpec)`; the doubling math then reduces to
#   the legacy Fortran TIR recipe (`rt_doubling.f90:191-197`).
#
# Concrete sources override these via methods in their respective files
# (solar_beam.jl returns `nothing` so the legacy path stays; thermal_emission.jl
# returns `:thermal` and `ones`).
# ============================================================================

"""
    source_key(::AbstractPreparedSource) -> Union{Nothing, Symbol}

Trait: which NamedTuple key (under `AddedLayer.j₀_by_src` /
`CompositeLayer.J₀_by_src`) carries this source's per-source j₀± / J₀±.

Returns `nothing` for sources that live in the legacy `j₀⁺/j₀⁻` slot
(solar, surface SIF) — `make_added_layer` will skip slot allocation for
them. Concrete v0.7 source types that want their own slot override this
to return a `Symbol` (e.g. `:thermal`).
"""
source_key(::AbstractPreparedSource) = nothing
source_key(::NoSource) = nothing

"""
    source_expk_init(prep, dτ_λ, μ₀, arr_type) -> AbstractArray{FT,1}

Trait: the initial per-spectral `expk` vector for this source's doubling
recurrence. Each doubling step squares it (so after `n` doublings the
factor reflects the full layer thickness from the source's reference
frame).

- Solar: `exp(-dτ_λ/μ₀)` — the direct-beam attenuation across the lower
  sub-layer. *Lives in the legacy path*; not consulted via this trait.
- Thermal: `ones(eltype(dτ_λ), length(dτ_λ))` — the bottom sub-layer's
  Planck emission is **not** attenuated relative to the top's (both
  sub-layers self-emit). The squaring is then a no-op (1² = 1), so the
  doubling recurrence collapses to the Fortran TIR formula
  `tj⁺ = j⁺ + tt_gp · (j⁺ + r⁻⁺ · j⁻)` (rt_doubling.f90:196).

Default implementation returns a `ones` vector matching the architecture
of `dτ_λ` (covers thermal and any future isotropic-self-generated source).
"""
function source_expk_init(::AbstractPreparedSource, dτ_λ::AbstractArray, _μ₀, arr_type)
    FT = eltype(dτ_λ)
    return arr_type(ones(FT, length(dτ_λ)))
end

# Dispatch table for `+` (debugging aid) =====================================
#
# Rules in order of decreasing specificity (Julia picks the most specific
# match):
#
#   (NoSource, NoSource)             → NoSource()
#   (NoSource, SourceSet)            → b
#   (SourceSet, NoSource)            → a
#   (NoSource, AbstractSource)       → b
#   (AbstractSource, NoSource)       → a
#   (SourceSet, SourceSet)           → flatten
#   (SourceSet, AbstractSource)      → append
#   (AbstractSource, SourceSet)      → prepend
#   (AbstractSource, AbstractSource) → SourceSet((a, b))
#
# All NoSource × SourceSet and SourceSet × NoSource paths are explicit so
# the (NoSource, AbstractSource) and (AbstractSource, SourceSet) methods
# never collide.
