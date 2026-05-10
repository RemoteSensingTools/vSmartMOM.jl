# Sources (v0.6 source-term refactor)

The radiative-transfer source term — what drives the affine RHS of the
MOM operator equation. v0.6 introduced an `AbstractSource` hierarchy
so multiple sources (solar + SIF, solar + thermal, etc.) can compose
cleanly. v0.7 surfaces the source list through the trait dispatch in
`component_m_max(...)`.

## Default

If no source is configured, `model_from_parameters` defaults to
`SolarBeam()` — equivalent to the pre-v0.6 fixed-solar-only behavior.

## Programmatic API (recommended)

Sources are first-class types under `vSmartMOM.CoreRT`. The cleanest
way to set them is via the `sources=` kwarg on `model_from_parameters`
or `rt_run`:

```julia
using vSmartMOM, vSmartMOM.CoreRT

# Solar-only (default)
model = model_from_parameters(params)

# Solar + surface SIF
model = model_from_parameters(params; sources = SolarBeam() + SurfaceSIF())

# Per-call override
R, T = rt_run(model; sources = BlackbodySource(temperature=300.0,
                                                spec_band=1))
```

## YAML configuration

A YAML-side `sources:` list is **forthcoming** (Phase 7-ish, not
shipping in v0.7). For now, configure sources programmatically after
`parameters_from_yaml`. The reason is that some source types (e.g.
`BlackbodySource`) need a model-resolved spectral band that doesn't
exist until after `model_from_parameters` runs.

## Source vocabulary

| Type | Purpose | `component_m_max(::T, ctx)` |
|------|---------|-----------------------------|
| `SolarBeam(; F₀=nothing, sza=nothing)` | Direct solar beam at the model's `sza` (or `obs_geom.sza` if not overridden). `F₀` is the per-spectral-point Stokes-vector solar irradiance (default Stokes-I unit if not set). | `0` (neutral — doesn't pin the loop) |
| `SurfaceSIF(; SIF₀=nothing)` | Sun-induced fluorescence at the lower boundary. Adds isotropically to `j₀⁻[:, 1, :]`. | `0` (m=0 only) |
| `BlackbodySource(temperature, spec_band; factor=π)` | Boundary/surface Planck emission at `T`. The `factor=π` default makes its `F₀` directly comparable to a unit `SolarBeam(F₀=1)`. | `0` (isotropic) |
| `NoSource()` | Identity / no-op source. Use as an explicit "I disabled the source" marker. | `0` |
| `SourceSet((s₁, s₂, ...))` | Tuple-typed composition. Built by `s₁ + s₂` operator. | `maximum(component_m_max.(sources, Ref(ctx)))` |
| `ThermalEmission` | **Type stub only in v0.6**, not yet wired. The atmospheric volume Planck integral. Reserved name. | — |
| `DiffuseBoundary` | **Type stub only**. Reserved for diffuse top-of-atmosphere or surface boundary sources. | — |

`SolarBeam → 0` is intentional — the trait aggregator takes a `max` over
components, so `0` is "neutral, let surfaces/scatterers drive". An
earlier draft used `typemax(Int)` and would have pinned every run to
`user_l_cap`; that was a Codex-flagged bug fixed before Phase C landed.

## Source-unit convention (v0.6)

All source intensities are in **mW · m⁻² · cm⁻¹**. `rt_run` output
(reflectance / transmittance) is in **mW · m⁻² · sr⁻¹ · cm⁻¹**.
`BlackbodySource` defaults `factor=π` so a unit Planck source is
comparable to `SolarBeam(F₀=1)`.

## Composability

```julia
src1 = SolarBeam()
src2 = SurfaceSIF()
src3 = BlackbodySource(temperature=300.0, spec_band=1)

# `+` flattens; no nested SourceSet
combined = src1 + src2 + src3   # SourceSet((src1, src2, src3))

# NoSource is identity
combined + NoSource() == combined   # true

# Iteration order is deterministic (a Tuple, not a Set)
for s in combined.sources
    println(typeof(s).name.name)
end
```

`SourceSet` is materialized to a concrete `Tuple` at
`model_from_parameters` time so the hot RT loop sees a fully
type-stable iterable.

## See also

- [`docs/src/pages/extending/sources.md`](../../extending/sources.md) —
  full architecture writeup (how to add a new source type)
- [`docs/src/pages/release_notes.md`](../../release_notes.md) — v0.6
  source-term refactor section
- [`Schema/radiative_transfer.md`](radiative_transfer.md) — Fourier
  loop bound interaction with source traits
