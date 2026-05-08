# Source terms

`vSmartMOM` v0.6 introduces a first-class **source-term** abstraction.
Solar beams, surface fluorescence, and (in the near future) thermal
emission and lidar pulses are all concrete subtypes of
`AbstractSource`, composed via `+`, and dispatched per layer through
multiple dispatch — no `if SFI` branching, no `RS_type.F₀`/`RS_type.SIF₀`
ownership leak.

The MOM solver is mathematically affine:

```
optical properties define the operator A
sources         define the additive RHS b
```

Doubling and adding propagate `(A, b)` generically. The new design makes
that explicit at the type level.

## User API

```julia
using vSmartMOM
using vSmartMOM.CoreRT

params = parameters_from_yaml("config/o2_a_band.yaml")
model  = model_from_parameters(params)

# Default — RTModel.sources defaults to a SolarBeam (unit Stokes I irradiance).
R, T = rt_run(model)

# Explicit solar beam with custom F₀ (irradiance, mW · m⁻² · cm⁻¹).
R, T = rt_run(model; sources = SolarBeam(F₀ = solar_irradiance))

# Solar + surface fluorescence (Lambertian-only emission at m=0).
R, T = rt_run(model; sources = SolarBeam() + SurfaceSIF(SIF₀ = sif_spec))

# Thermal-IR / Carbon-I-style scene (1500 K blackbody source at 2-2.4 µm).
spec_band = collect(model.atmosphere.spec_bands[1])   # cm⁻¹
sources   = BlackbodySource(1500, spec_band)          # SolarBeam with Planck F₀
R, T = rt_run(model; sources = sources)
```

`sources` can also be set at model build time:

```julia
model = model_from_parameters(params; sources = SolarBeam() + SurfaceSIF())
```

The `sources=` kwarg on `rt_run` overrides `model.sources` for that
specific solve; both flow through the same dispatch.

## Composition

Sources compose via `+` and are stored as a type-stable `SourceSet` of
concrete prepared types in the hot loop:

| Expression                                         | Result                                        |
|----------------------------------------------------|-----------------------------------------------|
| `NoSource() + s`                                   | `s`                                           |
| `SolarBeam() + SurfaceSIF()`                       | `SourceSet((SolarBeam, SurfaceSIF))`          |
| `(SolarBeam() + SurfaceSIF()) + ThermalEmission()` | `SourceSet((SolarBeam, SurfaceSIF, ...))`     |
| `SourceSet + SourceSet`                            | flattened — no nesting                        |

`SourceSet` iteration is unrolled at compile time, so dispatching to per-
source kernels is type-stable on both CPU and GPU.

## Built-in source types (v0.6)

| Type              | Where it contributes                       | Math                                                  |
|-------------------|--------------------------------------------|-------------------------------------------------------|
| `NoSource`        | nowhere                                    | identity for source composition                       |
| `SolarBeam`       | atmospheric layer `j₀±` per Fourier moment | exact finite-δ single-scatter (Fell 1997 Eqs. 1.52-1.54) |
| `BlackbodySource` | atmospheric layer `j₀±`                    | sugar around `SolarBeam` with `F₀ = factor · π · B(ν, T)` |
| `SurfaceSIF`      | surface layer `j₀⁻` at m=0                 | factor-2 broadcast across Nquad streams (Lambertian only) |

`BlackbodySource(T, spec_band)` is a constructor that returns a
`SolarBeam` with F₀ filled from the Planck spectrum at temperature `T`.
Use it for Carbon-I-like scenes (hot lab source illuminating CO₂/CH₄/H₂O
absorption in the 2-2.4 µm range).

Future versions will add `ThermalEmission` (atmospheric volume Planck
integral), `DiffuseBoundary`, and (later) `LidarPulse`.

## Units convention

| Quantity                   | Unit                              |
|----------------------------|-----------------------------------|
| `SolarBeam.F₀`             | mW · m⁻² · cm⁻¹                   |
| `SurfaceSIF.SIF₀`          | mW · m⁻² · cm⁻¹                   |
| `rt_run` output (R, T, J)  | mW · m⁻² · sr⁻¹ · cm⁻¹            |
| `SolarBeam(F₀ = ones(...))`| 1 mW · m⁻² · cm⁻¹ per spectral pt |

All sources in a `SourceSet` should share these units so additive
composition makes physical sense. `BlackbodySource` defaults to `factor =
π` (Lambertian-disk → hemispheric irradiance) so its `F₀` is comparable
to a `SolarBeam(F₀ = 1)` baseline.

## Differentiation

`AbstractSourceADMode` traits declare how each source's parameters
participate in linearized RT:

```julia
abstract type AbstractSourceADMode end
struct AnalyticSourceJacobian   <: AbstractSourceADMode end
struct ForwardDiffSourceJacobian <: AbstractSourceADMode end   # reserved for v0.7+
struct NoSourceJacobian          <: AbstractSourceADMode end
```

The `AD seam` is `prepare_source(::AbstractSource, FT, pol_n, nSpec,
arr_type)`. Above the seam, a future `prepare_source_with_tangent` can
trace through user parameters with `ForwardDiff`. Below, the kernels run
on plain `FT` arrays and the analytic Sanghavi 2014 App. C tangents (for
`SolarBeam`) — bit-equal to today's hand-written linearization.

## Dispatch architecture (for software folks)

```
                      AbstractSource
                            ▲
              ┌─────────────┼──────────────────────┐
              │             │                      │
        NoSource    SourceSet{S<:Tuple}     concrete sources
                                              ▲
                              ┌───────────────┼─────────────────┐
                          SolarBeam    SurfaceSIF     (ThermalEmission, ...)
```

`prepare_source` lifts each user-facing source to a `Prepared*` form on
the model's `arr_type` and `FT`:

```
                  AbstractPreparedSource
                            ▲
              ┌─────────────┼──────────────────────┐
              │             │                      │
        NoSource    SourceSet           concrete prepared sources
                                              ▲
                              ┌───────────────┼─────────────────┐
                     PreparedSolarBeam   PreparedSurfaceSIF      ...
```

The hot loop calls one of two dispatchers:

```
contribute!(prepared_sources, j₀⁺, j₀⁻, layer_ctx...)         # forward
source_tangent!(prepared_sources, j₀⁺, j₀⁻, J̇₀⁺, J̇₀⁻, ap..., layer_ctx...)  # linearized
```

`NoSource` → no-op; `SourceSet` → unrolled tuple iteration with per-
source dispatch; concrete prepared → kernel call.

For surface contributions, the same pattern applies but with double-
dispatch on (source-type, surface-type):

```
surface_source_contribute!(prepared_sources, surface, surface_added_layer, m, pol_type, arch)
```

The dispatch table below shows how each (source × surface) pair is
handled:

| Source                | Surface                 | Body                                   |
|-----------------------|-------------------------|----------------------------------------|
| `PreparedNoSource`    | any                     | no-op                                  |
| `PreparedSourceSet`   | any                     | iterate                                |
| `PreparedSolarBeam`   | any                     | (currently in `create_surface_layer!`; will move to dispatch in a later sub-phase) |
| `PreparedSurfaceSIF`  | `LambertianSurface*`    | factor-2 SIF₀ broadcast (m=0 only)     |
| `PreparedSurfaceSIF`  | non-Lambertian          | no-op                                  |

## Architecture invariants

1. **Sources only contribute on the elemental / surface step.** Doubling
   and interaction propagate `j₀±` and `J̇₀±` affinely without knowing
   which source produced them. This is why the elastic linearized path
   is fully `if SFI`-free in v0.6: the math handles `j₀±=0` as a natural
   no-op.
2. **Per-Fourier-moment cleanliness.** Sources are prepared once before
   the Fourier loop; each `m` runs through the same operator → source-
   contribution → propagation pipeline.
3. **AD seam at `prepare_source`.** User-parameter space (potentially
   `Dual`) lives above; kernel-space (plain `FT`) lives below.
4. **`RTModel.sources` is the source of truth.** Defaults to
   `SolarBeam()` to preserve historical unit-Stokes-I behavior; users
   override at construction or per-`rt_run` call.

## See also

- [`src/CoreRT/Sources/types.jl`](https://github.com/cfranken/vSmartMOM.jl/blob/main/src/CoreRT/Sources/types.jl) — `AbstractSource`, `SourceSet`, `NoSource`, AD-mode traits.
- [`src/CoreRT/Sources/solar_beam.jl`](https://github.com/cfranken/vSmartMOM.jl/blob/main/src/CoreRT/Sources/solar_beam.jl) — `SolarBeam`, `PreparedSolarBeam`, `BlackbodySource`, `source_tangent!`.
- [`src/CoreRT/Sources/surface_sif.jl`](https://github.com/cfranken/vSmartMOM.jl/blob/main/src/CoreRT/Sources/surface_sif.jl) — `SurfaceSIF`, `surface_source_contribute!`.
- [`test/test_sources.jl`](https://github.com/cfranken/vSmartMOM.jl/blob/main/test/test_sources.jl) — per-phase regression tests with end-to-end bit-equality assertions.
- The original v0.6 design note: [`dev_notes/source_terms_architecture_v0_6.md`](https://github.com/cfranken/vSmartMOM.jl/blob/main/dev_notes/source_terms_architecture_v0_6.md).
- The execution plan: `~/.claude/plans/gpt-also-had-some-velvety-whale.md`.

## API reference

### Source vocabulary

```@docs
vSmartMOM.CoreRT.AbstractSource
vSmartMOM.CoreRT.AbstractPreparedSource
vSmartMOM.CoreRT.NoSource
vSmartMOM.CoreRT.SourceSet
```

### Concrete sources

```@docs
vSmartMOM.CoreRT.SolarBeam
vSmartMOM.CoreRT.PreparedSolarBeam
vSmartMOM.CoreRT.BlackbodySource(::Real, ::AbstractVector{<:Real})
vSmartMOM.CoreRT.SurfaceSIF
vSmartMOM.CoreRT.PreparedSurfaceSIF
```

### AD-mode traits

```@docs
vSmartMOM.CoreRT.AbstractSourceADMode
vSmartMOM.CoreRT.AnalyticSourceJacobian
vSmartMOM.CoreRT.ForwardDiffSourceJacobian
vSmartMOM.CoreRT.NoSourceJacobian
vSmartMOM.CoreRT.source_ad_mode(::vSmartMOM.CoreRT.AbstractSource)
```

### Dispatch entry points

```@docs
vSmartMOM.CoreRT.prepare_source(::vSmartMOM.CoreRT.SolarBeam, ::Type{<:AbstractFloat}, ::Integer, ::Integer, ::Any)
vSmartMOM.CoreRT.prepare_source(::vSmartMOM.CoreRT.SurfaceSIF, ::Type{<:AbstractFloat}, ::Integer, ::Integer, ::Any)
vSmartMOM.CoreRT.prepare_sources(::vSmartMOM.CoreRT.NoSource, ::Type{<:AbstractFloat}, ::Integer, ::Integer, ::Any)
vSmartMOM.CoreRT.surface_source_contribute!(::vSmartMOM.CoreRT.PreparedSurfaceSIF, ::Union{vSmartMOM.CoreRT.LambertianSurfaceLegendre, vSmartMOM.CoreRT.LambertianSurfaceScalar, vSmartMOM.CoreRT.LambertianSurfaceSpline}, ::Any, ::Integer, ::Any, ::Any)
vSmartMOM.CoreRT.surface_source_contribute!(::vSmartMOM.CoreRT.NoSource, ::vSmartMOM.CoreRT.AbstractSurfaceType, ::Vararg{Any, 4})
```
