# vSmartMOM.jl v0.6 Source-Term Architecture

**Status:** working design note for v0.6 discussion.
**Scope:** no implementation here; this is the intended architecture for cleanly composing solar SFI, thermal emission, SIF, lidar/lunar beams, and future sources.
**Related notes:** `dev_notes/vsmartmom_dispatch_design_v0_6.md`, `dev_notes/standalone_ss_solver_plan.md`, `docs/dev_notes/LINEARIZATION_BUGS.md`.

## Motivation

The current source handling is physically correct in important places, but the ownership boundaries are blurry:

- The direct solar source-function integration path is wired through `RS_type.F₀`, `SFI::Bool`, and `get_elem_rt_SFI!`.
- Surface SIF is injected inside Lambertian surface machinery via `inject_surface_SIF!`, even though SIF is a source term, not a BRDF property.
- `AbstractRamanType` currently mixes medium physics, inelastic redistribution state, and external source vectors (`F₀`, `SIF₀`).
- Adding thermal emission, lidar, lunar illumination, or multiple concurrent sources would require more cross-cutting edits in `rt_run`, elemental kernels, surfaces, and linearization code.

The v0.6 goal should be to make source terms first-class, composable objects while keeping the MOM operator path clean.

## Core Rule

For fixed optical properties, the MOM RT problem is affine:

```text
Optical properties define the operator.
Sources define additive right-hand-side terms.
```

The forward elemental step should conceptually become:

```julia
r, t = layer_operator(optics)
j    = layer_source(sources, optics, state)
```

The adding/doubling code should not know whether `j` came from solar scattering, thermal emission, SIF, lidar, or another source. It should only propagate `J`.

## User-Facing Shape

Define a source abstraction:

```julia
abstract type AbstractSourceType end

struct NoSource{FT} <: AbstractSourceType end

struct SourceSet{S<:Tuple}
    sources::S
end
```

Users should be able to write:

```julia
sources = SourceSet((
    SolarBeamSource(...),
    ThermalEmissionSource(...),
    SurfaceSIFSource(...),
))
```

The source tuple should be concrete and type-stable. Avoid `Vector{AbstractSourceType}` in hot paths.

`SourceSet` should support the same composition idiom as optical-property mixing:

```julia
src = SolarBeamSource(F₀=solar_spec, sza=35) + ThermalEmissionSource(T=profile_T)
src = src + SurfaceSIFSource(strength=sif, spectrum=sif_spectrum)
```

This should produce:

```julia
SourceSet((SolarBeamSource(...), ThermalEmissionSource(...), SurfaceSIFSource(...)))
```

`Base.iterate(::SourceSet)` should expose the tuple so internal code can write:

```julia
for src in sources
    contribute!(j₀⁺, j₀⁻, src, layer_ctx)
end
```

`NoSource()` is useful as an explicit, dispatchable default for no active sources; `SourceSet(())` can remain the structural representation of an empty list.

## Source Contract

Each source type should implement a small interface:

```julia
prepare_source(source, model, iBand, arch) -> prepared_source

source_lmax_requirement(source, model, iBand) -> Int
source_mmax_requirement(source, model, iBand) -> Int

add_layer_source!(j₀⁺, j₀⁻, prepared_source, layer_ctx)
add_surface_source!(j₀⁺, j₀⁻, prepared_source, surface_ctx)
add_boundary_source!(J₀⁺, J₀⁻, prepared_source, boundary_ctx)
```

Most sources implement only one of the `add_*_source!` methods; the others are no-ops.

The minimal forward contract can be summarized as one method:

```julia
contribute!(j₀⁺, j₀⁻, source, layer_ctx)
```

`contribute!` adds this source's current layer contribution in place. Self-gating belongs inside the source implementation: a surface-only source returns immediately for atmospheric layers; a thermal source returns immediately for unsupported Fourier moments; a disabled source returns immediately for zero amplitude.

The elemental source path should become:

```julia
fill!(j₀⁺, 0)
fill!(j₀⁻, 0)
for src in sources
    contribute!(j₀⁺, j₀⁻, src, layer_ctx)
end
```

Doubling and interaction stay source-agnostic.

Optional location traits can support fast paths without becoming correctness requirements:

```julia
abstract type SourceLocation end
struct ExternalSourceLocation     <: SourceLocation end
struct InAtmosphereSourceLocation <: SourceLocation end
struct AtSurfaceSourceLocation    <: SourceLocation end

source_location(::SolarBeamSource) = ExternalSourceLocation()
source_location(::ThermalEmissionSource) = InAtmosphereSourceLocation()
source_location(::SurfaceSIFSource) = AtSurfaceSourceLocation()
```

These traits let a kernel or driver skip whole classes of no-op calls when profitable, while retaining the simple `contribute!` contract.

## Source Categories

### `CollimatedBeamSource`

External directed illumination: solar, lunar, stellar, or active lidar.

Owns:

- incident direction and boundary side
- Stokes spectrum `F₀`
- azimuth convention
- optional beam geometry metadata

This replaces the current hardwired solar SFI `F₀` path. It contributes atmospheric layer source terms through exact finite-δ single-scatter formulas. Solar and lunar are just different instances; lidar may later add path/gating metadata.

The initial implementation should likely be named `SolarBeam` or `SolarBeamSource` and preserve today's unpolarized solar convention by default, while allowing a full incident Stokes spectrum later:

```julia
SolarBeamSource(F₀=solar_spec, sza=35)
SolarBeamSource(F₀=solar_stokes_spec, direction=..., polarization=...)
```

### `ThermalEmissionSource`

Atmospheric volume emission. It contributes directly to `j₀⁺/j₀⁻`, typically only for `m=0`.

The intended numerical implementation is the exact layer integral for a source function linear in optical depth, using the within-layer Planck/temperature-gradient method. This source should own the thermal emission model, while the optical operator remains unchanged.

The source API should decide whether user-facing temperatures are specified at layer centers or layer boundaries. The exact linear-in-τ method naturally wants boundary values or an internally constructed `dB/dτ`; if users provide center temperatures, conversion belongs in `prepare_source`.

### `SurfaceEmissionSource`

Surface-origin emission: SIF, surface thermal emission, fluorescence from canopy, etc.

This should replace `inject_surface_SIF!` living inside Lambertian surface handling. The surface BRDF remains responsible for reflection matrices. Surface sources are responsible for emission vectors.

For SIF, the source is isotropic and only contributes to `m=0`; any normalization needed to compensate downstream Fourier weighting belongs in the source implementation, documented there.

The design preference is that SIF lives as a source, not as an extension of `AbstractSurfaceType`. The surface BRDF owns reflection; `SurfaceSIFSource` owns emitted radiance. A surface may still provide metadata or convenience constructors for source creation.

### `DiffuseBoundarySource`

User-specified diffuse radiance/irradiance incident at TOA or BOA on the quadrature streams. This is useful for validation, thermal/MW boundary conditions, and non-solar backgrounds.

### Inelastic Processes Are Not Just Sources

Raman/VS should be represented as medium redistribution processes, not merely as `AbstractSourceType`s. A Raman-lidar scene is:

```text
collimated source beam + inelastic redistribution process
```

This keeps source ownership separate from wavelength-redistribution physics. In the long run, `RS_type` should lose `F₀/SIF₀` ownership and become an inelastic/redistribution configuration.

## Multiple Source Channels

Default output can sum source contributions into the existing one-column `J`:

```julia
J[nμ, 1, nSpec]
```

But the architecture should allow source-separated propagation:

```julia
J[nμ, nSourceChannels, nSpec]
```

The adding/doubling equations already support this mathematically because source propagation is matrix times RHS. A multi-column `J` lets us propagate solar, thermal, SIF, and diagnostic components in one RT solve, then either sum channels or expose them separately.

This is also useful for Jacobians and debugging because source-separated outputs make it clear which physical term caused a residual.

## Differentiation Boundary

This design should simplify the AD/hand-written interface.

The operator and source derivatives are separate:

```julia
ṙ, ṫ = layer_operator_tangent(optics, opticṡ)
j̇     = layer_source_tangent(sources, optics, opticṡ, source_statė)
```

Adding/doubling then uses a generic affine chain rule:

```text
J_out  = A(optics) * J_in + j(source, optics)
J̇_out = Ȧ * J_in + A * J̇_in + j̇
```

The key rule:

```text
source_tangent! returns only source RHS derivatives.
It does not know about doubling, interaction, or postprocessing.
```

Each source can declare its differentiation mode:

```julia
abstract type AbstractSourceADMode end
struct AnalyticSourceJacobian <: AbstractSourceADMode end
struct ForwardDiffSourceJacobian <: AbstractSourceADMode end
struct NoSourceJacobian <: AbstractSourceADMode end
```

Examples:

- `SolarBeamSource`: analytic tangent for `τ`, `ϖ`, `Z`, `F₀`, and possibly geometry.
- `ThermalEmissionSource`: analytic tangent for temperature, Planck function, optical depth, and layer-gradient terms.
- `SurfaceSIFSource`: analytic or AD tangent depending on whether SIF spectrum/scale is a retrieval parameter.
- `LidarSource`: analytic tangent for amplitude/pointing later; no Jacobian initially is acceptable.

The AD boundary should be:

```text
AD upstream:
  user parameters -> prepared optical/source scalar fields

Hand-written downstream:
  MOM operator propagation + affine source propagation
```

This avoids editing the adding/doubling linearized machinery every time a new source is added.

There are two differentiated surfaces:

| Parameter family | Preferred path | Reason |
|---|---|---|
| Source parameters (`F₀`, thermal profile, SIF strength/spectrum, lidar amplitude) | `ForwardDiff` through type-generic `contribute!` | Few parameters per source; source-local AD is acceptable and keeps source physics self-contained. |
| Core optical properties (`τ`, `ϖ`, `Z`) | hand-written analytic `source_tangent!` / `linearize_core_OD!` per source | Many layer/spectral parameters; closed forms are short and performance-sensitive. |

Source superposition is linear:

```text
dJ_total/dx = Σ_source dJ_source/dx
```

So source Jacobians do not produce cross-terms between sources. The combined linearized source path mirrors the forward path:

```julia
for src in sources
    source_tangent!(j̇₀⁺, j̇₀⁻, src, layer_ctx, layer_lin_ctx)
end
```

Today's linearized solar elemental source becomes the `SolarBeamSource` analytic tangent rather than a special case in the global linearized elemental kernel.

## Fourier And Stream Requirements

Sources should participate in the same band-level angular-resolution logic as surfaces and optical properties.

Each source reports requirements:

```julia
source_lmax_requirement(source, model, iBand)
source_mmax_requirement(source, model, iBand)
```

Then per band:

```text
lmax_band = min(l_trunc, max(optics_lmax, surface_lmax, source_lmax...))
mmax_band = compatible_mmax(lmax_band)
nstreams_band = streams_for_lmax(lmax_band)
```

`l_trunc` is only a ceiling. It should not force every band/source/surface to use the same angular order. Rayleigh over Lambertian with thermal emission should be cheap; coarse aerosol or wind-broadened Cox-Munk may request more.

## Dispatch Rule

Use multiple dispatch for physical source categories and workflow boundaries.

Use plain `if/else` inside numerical kernels for local decisions:

- `m == 0`
- singular `μᵢ == μ₀` or `μᵢ == μⱼ` limits
- zero-weight user VZA nodes
- thermal layer-gradient degeneracies
- optional source disabled by zero amplitude

This matches the v0.6 dispatch principle: dispatch for physics, branches for scalar numerics.

## Documentation Contract For Sources

Every source type must document:

- where it injects: atmospheric layer, surface, boundary, or redistribution process
- direction convention: `+` down, `-` up
- Fourier normalization and `m=0` treatment
- whether it is collimated, diffuse, isotropic, or angularly resolved
- whether it changes angular-resolution requirements
- supported differentiation mode
- source-channel behavior: summed-only or separable
- expected units of the source spectrum

Users should be able to understand exactly how a source contributes without reading surface or Raman internals.

Every concrete source should have a docstring with:

```julia
SolarBeamSource(...)
ThermalEmissionSource(...)
SurfaceSIFSource(...)
```

including one short `rt_run(model; sources=...)` example. Concentrating source physics in source docstrings is part of the user-facing cleanup.

## Migration Plan

1. Introduce `AbstractSourceType`, `SourceSet`, and prepared source states.
2. Add `SolarBeamSource` and `contribute!(::SolarBeamSource, ...)`; refactor `get_elem_rt_SFI!` around this implementation with zero behavior change. `rt_run` accepts `sources=` and defaults to the current solar beam.
3. Treat `SFI::Bool` as legacy/internal: `true` maps to a default `SolarBeamSource`; `false` maps to `NoSource()` or `SourceSet(())`. Deprecate the flag after source tests are stable.
4. Move `inject_surface_SIF!` behavior into `SurfaceSIFSource`; keep a compatibility shim during transition.
5. Add `ThermalEmissionSource` as the first genuinely new source. This validates the layer-source interface and the thermal/MW direction.
6. Split `RS_type` responsibilities: keep inelastic redistribution in Raman/VS process types, move external illumination and SIF into sources.
7. Add source-separated RHS columns only after the one-column summed `SourceSet` path is stable.
8. Extend linearization by implementing `source_tangent!` per source rather than editing adding/doubling chain-rule code.
9. Defer `LidarPulseSource` and time-resolved source axes to v0.7+; the source abstraction should not force time-bin support into v0.6.

## Example User API

Backward-compatible default:

```julia
R, T = rt_run(model)
```

Explicit solar:

```julia
src = SolarBeamSource(F₀=solar_spec, sza=35)
R, T = rt_run(model; sources=src)
```

Multi-source:

```julia
src = SolarBeamSource(F₀=solar_spec, sza=35) +
      ThermalEmissionSource(T=atmosphere.T) +
      SurfaceSIFSource(strength=sif, spectrum=sif_spectrum)

R, T = rt_run(model; sources=src)
```

Source-parameter AD:

```julia
dR_dF₀ = ForwardDiff.derivative(1.0) do F0
    rt_run(model; sources=SolarBeamSource(F₀=F0, sza=35))[1]
end
```

Core optical-property analytic linearization:

```julia
R, T, dR, dT = rt_run(model, lin_model; sources=src)
```

## Open Questions

1. **Tuple vs vector.** The favored internal representation is `SourceSet{<:Tuple}` for type stability and inlining. Do we also need a runtime-flexible vector wrapper for user-generated source lists, converting to a tuple at model construction?
2. **Incident polarization API.** Should every external source expose `incident_stokes(source, iBand)` or should `SolarBeamSource` preserve the current unpolarized convention and add full Stokes support later?
3. **Spectral source representation.** Should source spectra be precomputed arrays per band, or callables evaluated in `prepare_source`? Current vSmartMOM patterns favor precomputed arrays.
4. **Thermal profile granularity.** Should users supply layer-center temperatures, boundary temperatures, or Planck source values directly?
5. **Surface-source coupling.** Can `SurfaceSIFSource` remain fully independent of surface BRDFs, or should surfaces provide source-construction helpers for physically coupled cases like canopy fluorescence?
6. **No-source default.** Is `NoSource()` useful enough as a public singleton, or should the public default be `sources=SourceSet(())`?
7. **Flag deprecation.** How long should `SFI::Bool` remain as a compatibility alias?

## Desired End State

Adding a new source should require:

```text
define source type
prepare source state
implement forward add_*_source!
optionally implement source_tangent!
document normalization and units
```

It should not require changing:

```text
elemental operator construction
doubling
interaction
surface BRDF reflection matrices
generic linearized adding/doubling propagation
```

That is the main payoff: sources become clean affine RHS contributors, while optical operators, surface reflection, inelastic redistribution, and linearization stay better delineated.
