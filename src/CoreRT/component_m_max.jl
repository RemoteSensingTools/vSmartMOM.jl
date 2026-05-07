# =========================================================================
# component_m_max — per-component Fourier-support traits
# =========================================================================
#
# Phase C of the Fourier/Stream Resolution refactor (plan:
# ~/.claude/plans/gpt-also-had-some-velvety-whale.md). Each component
# (surface, source, scatterer) declares the maximum Fourier order it
# contributes to via `component_m_max(c, ctx) :: Int`. The model-build
# aggregator takes `maximum(component_m_max(...))` across the active
# components and clamps against the user's resolving cap.
#
# Phase C lands the dispatch *infrastructure* behind a flag on
# `SolverConfig.use_component_traits` (default `false`). When `false`,
# the historical aggregator at `model_from_parameters.jl` keeps running
# and behavior is bit-equal to Phase B. When flipped to `true` (a
# follow-up commit, after Sanghavi sign-off on the Cox-Munk runtime
# impact), Cox-Munk / RossLi / RPV / canopy stop being silently
# half-truncated — they get their full `user_l_cap` Fourier resolution.
#
# Trait values:
# - LambertianSurface*  → 0   (m=0 is exact for any Lambertian)
# - CoxMunk/rpv/RossLi/Canopy → ctx.user_l_cap (no scheme-imposed cap)
# - RayleighScattering  → 2   (β₀, β₁, β₂ exhaust the phase function)
# - AerosolOptics       → length(greek.β) - 1 (under NoTruncation),
#                          clamped to user_l_cap (under δBGE)
# - SolarBeam           → 0   (neutral — see note below)
# - SurfaceSIF          → 0   (isotropic)
# - NoSource            → 0
# - SourceSet           → maximum across its sources
#
# `SolarBeam → 0` is deliberate. The naïve choice `typemax(Int)` would
# pin every run to `ctx.user_l_cap` because the aggregator takes a max
# — exactly the opposite of "the solar beam doesn't dictate the loop
# bound." With `0`, the aggregator falls through to the surface and
# scatterer contributions, which is the correct behavior.

"""
    component_m_max(component, ctx) -> Int

Maximum Fourier order this component contributes to. Used by
[`_aggregate_m_max`](@ref) when `SolverConfig.use_component_traits == true`.
The `ctx` NamedTuple carries `user_l_cap`, `stream_l_cap`, and
`truncation` so each component can decide whether to defer to the
projection cap or to declare a tighter bound of its own.
"""
function component_m_max end

# Default: no specific opinion. Falls through to `ctx.user_l_cap` —
# this is a conservative, "I cover the full stream cap" stance for
# components that haven't opted into a tighter trait.
component_m_max(::Any, ctx) = ctx.user_l_cap

# ---- Sources ------------------------------------------------------------
# `0` is neutral against a max-only aggregator — it doesn't pin the
# loop, it lets surface/scatterer drive.
component_m_max(::SolarBeam, ctx) = 0
component_m_max(::SurfaceSIF, ctx) = 0
component_m_max(::NoSource, ctx) = 0
component_m_max(s::SourceSet, ctx) = isempty(s.sources) ? 0 :
    maximum(component_m_max(src, ctx) for src in s.sources)

# ---- Surfaces -----------------------------------------------------------
component_m_max(::LambertianSurfaceScalar, ctx) = 0
component_m_max(::LambertianSurfaceSpectrum, ctx) = 0
component_m_max(::LambertianSurfaceLegendre, ctx) = 0
component_m_max(::LambertianSurfaceSpline, ctx) = 0

# Cox-Munk, RPV, Ross-Li, canopy: the scheme has no built-in finite
# Fourier support, so we let the user's projection cap drive. Phase
# C-flip turns this into the full Cox-Munk Fourier resolution that
# the previous aggregator was silently halving.
component_m_max(::CoxMunkSurface, ctx) = ctx.user_l_cap
component_m_max(::rpvSurfaceScalar, ctx) = ctx.user_l_cap
component_m_max(::RossLiSurfaceScalar, ctx) = ctx.user_l_cap
component_m_max(::CanopySurface, ctx) = ctx.user_l_cap

# ---- Scatterers ---------------------------------------------------------

# Rayleigh has only β₀, β₁, β₂ non-zero — exact m_max = 2. We support
# both an instance dispatch (when a constructed RayleighScattering is
# in scope) and a Type{...} dispatch (used at model-build time when
# the full Rayleigh state hasn't been assembled yet — see
# `_band_components` in `model_from_parameters.jl`).
component_m_max(::RayleighScattering, ctx) = 2
component_m_max(::Type{RayleighScattering}, ctx) = 2

# AerosolOptics: rely on the truncated greek-coef array length.
# Under δBGE the array is already clamped to `user_l_cap`. Under
# NoTruncation the length reflects the Mie series; we let it through
# (a coarse aerosol with too many coefs will hit the parser's
# validate-against-`user_l_cap` guard at config time, not here).
function component_m_max(a::Scattering.AerosolOptics, ctx)
    β = a.greek_coefs.β
    return max(length(β) - 1, 0)
end

"""
    _aggregate_m_max(active_components, ctx) -> Int

Aggregate the per-component Fourier supports into a single per-band
loop bound (order). Internal helper used by `model_from_parameters`
when `SolverConfig.use_component_traits == true`.

`active_components` is an iterable of components participating in the
band — typically `(rayleigh, aerosols..., surface, sources)`. The
aggregator takes the maximum Fourier support, then clamps to
`ctx.user_l_cap` and (if non-nothing) `ctx.m_max_override`.
"""
function _aggregate_m_max(active_components, ctx)
    isempty(active_components) && return 0
    m_max = maximum(component_m_max(c, ctx) for c in active_components)
    m_max = min(m_max, ctx.user_l_cap)
    if haskey(ctx, :m_max_override) && ctx.m_max_override !== nothing
        m_max = min(m_max, ctx.m_max_override)
    end
    return max(m_max, 0)
end
