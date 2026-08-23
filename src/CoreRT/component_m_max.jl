# =========================================================================
# component_m_max — per-component Fourier-support traits
# =========================================================================
#
# Phase C of the Fourier/Stream Resolution refactor. Each component
# (surface, source, scatterer) declares the maximum Fourier order it
# contributes to via `component_m_max(c, ctx) :: Int`. The model-build
# aggregator takes `maximum(component_m_max(...))` across the active
# components and clamps against the user's resolving cap.
#
# Phase C lands the dispatch *infrastructure* behind the flag
# `SolverConfig.use_component_traits` — **default `true`** in v2.1.
#
# Relation to the runtime Fourier-convergence test: these traits fix the
# STATIC per-band ceiling `m_max` at model build (what the scene could
# require); `RTNumericalParameters.fourier_convergence` (see
# `AbstractFourierConvergence`) may then terminate the loop BELOW that
# ceiling at run time, once the accumulated Stokes-I contribution at the
# user view angles stabilises. The two are independent and compose:
# traits bound the loop, convergence exits it early.
# Cox-Munk / RossLi / RPV / canopy now run to their full
# `user_l_cap` Fourier resolution instead of being silently
# half-truncated by the legacy aggregator.  Flip to `false` to fall
# back to the historical `min(ceil((l_max+1)/2), params.max_m)`
# aggregator at `model_from_parameters.jl` (bit-equal to Phase B —
# kept as an escape hatch for byte-equality regression).
#
# Trait values:
# - LambertianSurface*  → 0   (m=0 is exact for any Lambertian)
# - CoxMunk/rpv/RossLi/Canopy → ctx.user_l_cap (no scheme-imposed cap)
# - RayleighScattering  → 2   (β₀, β₁, β₂ exhaust the phase function)
# - AerosolOptics       → largest retained β degree. By default this is
#                          length(greek.β)-1; when `greek_beta_cutoff` is set,
#                          it is the last l with |β_l| >= cutoff at any band
#                          wavelength, subsequently capped by user_l_cap.
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
The `ctx` NamedTuple carries `user_l_cap`, `stream_l_cap`, `truncation`, and
optionally `greek_beta_cutoff`, so each component can declare a tighter bound.
The aerosol cutoff applies only to `abs(β_l)`; the other Greek families are
not Fourier-support tests.
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

# AerosolOptics uses a two-stage support estimate. Mie integration first
# allocates the conservative series implied by the maximum size parameter.
# The optional post-integration filter then removes a numerically negligible
# high-l tail by finding the last |β_l| above `greek_beta_cutoff`.
# Under δBGE the array is already clamped to `user_l_cap`. Under
# NoTruncation the length reflects the Mie series; we let it through
# (a coarse aerosol with too many coefs will hit the parser's
# validate-against-`user_l_cap` guard at config time, not here).
function component_m_max(a::Scattering.AerosolOptics, ctx)
    β = a.greek_coefs.β
    full_support = max(size(β, 1) - 1, 0)
    cutoff = haskey(ctx, :greek_beta_cutoff) ? ctx.greek_beta_cutoff : nothing
    cutoff === nothing && return full_support

    # β can be `(l,)` or `(l,nSpec)`. Retain l when at least one spectral point
    # reaches the threshold, making the decision conservative over the band.
    # Use magnitude because β_l can be signed. Do not use α/γ/δ/ϵ/ζ: those
    # families may legitimately be zero or change sign and are not the scalar
    # phase-function support benchmark. β₀ is normalized to one, so a valid
    # positive cutoff always leaves at least m=0.
    βhost = Array(β)
    retained = if ndims(βhost) == 1
        findlast(x -> abs(x) >= cutoff, βhost)
    else
        findlast(il -> maximum(abs, @view(βhost[il, :])) >= cutoff,
                 axes(βhost, 1))
    end
    return retained === nothing ? 0 : retained - 1
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

Note the naming collision with `rt_run`'s `m_max_override` keyword, which
has the OPPOSITE sense: here (model build) the override NARROWS the trait
ceiling; there (runtime, used by `rt_run_atmosphere` to size surface-replay
caches) it only WIDENS the loop bound. They are independent knobs that
happen to share a name.
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
