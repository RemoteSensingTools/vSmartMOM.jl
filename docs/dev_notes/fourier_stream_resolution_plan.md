# Fourier / Stream Resolution Setup Plan

Review target: Sanghavi and Claude.

Status: proposal for an implementable near-term refactor. The goal is to make
the RT setup easier for users while preserving the important numerical
relationships between spherical harmonic order, Fourier order, and quadrature
streams.

Sanghavi review update: use `nstreams` as the primary user choice. These are
streams **per hemisphere**. They are not "2-stream solver" terminology:
classical two-stream means one stream per hemisphere, whereas
`nstreams = 3` here means three upward and three downward discrete ordinates
before any zero-weight output nodes are appended.

## Problem

The current public setup asks users for `max_m` and `l_trunc`, while `Nquad`
is derived indirectly by `rt_set_streams`. This creates three problems:

1. `max_m` currently behaves like a **count** in code (`for m = 0:max_m-1`),
   but the name reads like a maximum Fourier order.
2. `l_trunc` is often treated as the actual order for all cases, when it
   should be a **ceiling**. Rayleigh over Lambertian only needs `l_max = 2`;
   a smooth fine-mode aerosol may need less than the ceiling; a coarse aerosol
   may need more but must be truncated to the ceiling.
3. `Nquad` includes zero-weight SZA/VZA output nodes after quadrature setup, so
   it is not the right public resolution knob. The true solver resolution is
   the base number of half-space streams before output-node augmentation.

## Guiding Rules

Use Sanghavi's convention:

```julia
m_max == l_max
n_fourier_moments = m_max + 1
Nstreams = (l_max + 2) ÷ 2
```

Consequences:

| `l_max` | `Nstreams` | `m_max` | Fourier loop |
|---:|---:|---:|---|
| 0 | 1 | 0 | `m = 0:0` |
| 2 | 2 | 2 | `m = 0:2` |
| 15 | 8 | 15 | `m = 0:15` |
| 25 | 13 | 25 | `m = 0:25` |
| 30 | 16 | 30 | `m = 0:30` |

For even `l_max`, the formula gives `l_max = 2*Nstreams - 2`; for odd
`l_max`, it gives `l_max = 2*Nstreams - 1`. This avoids the `l_max = 0`
corner case without special casing.

## Public API Direction

Users should not specify `Nquad`; it is an augmented kernel dimension after
zero-weight SZA/VZA output nodes are inserted.

Make `nstreams` the main user-facing resolution knob. This is more intuitive
for VLIDORT/DISORT users and maps directly to solver cost. The harmonic-order
ceiling implied by streams is:

```julia
stream_l_cap = 2*nstreams - 1
```

Recommended user-facing fields:

```yaml
radiative_transfer:
  nstreams: 13         # streams per hemisphere; implies stream_l_cap = 25
  m_max: null          # optional lower hard cap; null means auto up to stream_l_cap
  truncation: null     # no phase-function forward-peak modification
  convergence: null    # optional Fourier convergence; null means fixed loop
```

TOML has no `null`, so TOML users should omit optional auto fields:

```toml
[radiative_transfer]
nstreams = 13
# m_max omitted means auto
# truncation omitted means no phase-function forward-peak modification
# convergence omitted means fixed loop
```

`truncation` is not the same thing as `nstreams`. `nstreams` controls the
solver's angular resolving power. `truncation` controls whether the phase
function is modified before projection into that finite resolving power:

```yaml
radiative_transfer:
  nstreams: 13
  truncation: "δBGE(25, 2.0)"      # truncate/compress a sharp forward peak
```

`truncation: null` / omitted should mean `NoTruncation()`: use the phase
coefficients as-is when they fit within `stream_l_cap`. If an aerosol phase
function has more coefficients than the stream cap and no truncation method is
specified, the setup should warn or error rather than silently applying a plain
coefficient cutoff. Users can then either increase `nstreams`, choose `δBGE`,
or explicitly request a plain cutoff if we decide to expose that advanced mode.

Keep `l_trunc` as a legacy/expert alias for an exact harmonic cap. If `l_trunc`
is supplied without `nstreams`, derive:

```julia
nstreams = (l_trunc + 2) ÷ 2
m_max = l_trunc
```

If both `nstreams` and `l_trunc` are supplied, require
`l_trunc <= 2*nstreams - 1`; otherwise the requested harmonic cap is not
resolvable by the requested streams.

Validation rule:

- Solar / scattering scenes should require `nstreams >= 3` because Rayleigh
  needs `m = 0, 1, 2`.
- Thermal-only or source-only scenes may allow `nstreams = 1` later, but that
  should be tied to an explicit source / RT-mode decision rather than accepted
  silently for solar-scattering runs.
- Error messages should say "streams per hemisphere" so users do not confuse
  `nstreams = 2` here with a classical two-stream solver.

## File Format / Reader Migration

This change touches the full config path:

```text
YAML/TOML/Dict
  -> parameters_from_file / read_parameters
  -> parameters_from_dict
  -> vSmartMOM_Parameters
  -> model_from_parameters
  -> RTModel / SolverConfig / QuadPoints
```

The migration should happen at the parser boundary first, so the rest of the
code receives normalized fields.

### Schema

Current schema requires `radiative_transfer.max_m`. New schema should:

1. Remove `max_m` from required fields.
2. Add `nstreams` as the preferred resolution field. It may be optional in the
   parser with default `13`, but new examples should write it explicitly.
3. Add optional `m_max` as either integer or null for YAML/Dict.
4. Add optional `truncation`, allowing null / omission for no phase-function
   forward-peak truncation.
5. Keep optional `l_trunc` and `max_m` for backward compatibility.
6. Validate `nstreams >= 3` for solar/scattering scenes; reserve
   `nstreams < 3` for future explicit thermal-only / source-only modes.
7. Add optional `convergence`.

TOML-specific rule: because TOML has no null, absence means auto. This means
`parameters_from_dict` must use `haskey`, not only `get(..., nothing)`, because
YAML `null` and missing TOML keys both need to normalize to `nothing`.

Suggested normalized parse:

```julia
rt = params_dict["radiative_transfer"]

legacy_l_trunc =
    haskey(rt, "l_trunc") ? Int(rt["l_trunc"]) : nothing

nstreams =
    if haskey(rt, "nstreams")
        Int(rt["nstreams"])
    elseif legacy_l_trunc !== nothing
        (legacy_l_trunc + 2) ÷ 2
    else
        13
    end

stream_l_cap = 2*nstreams - 1

m_max_override =
    if haskey(rt, "m_max")
        isnothing(rt["m_max"]) ? nothing : Int(rt["m_max"])
    elseif haskey(rt, "max_m")
        # Legacy `max_m` was a count. Convert count -> maximum order.
        max(Int(rt["max_m"]) - 1, 0)
    else
        nothing
    end

legacy_l_cap_override = legacy_l_trunc
user_l_cap = something(m_max_override, legacy_l_cap_override, stream_l_cap)
```

Validation:

```julia
if user_l_cap > stream_l_cap
    throw(ArgumentError(
        "requested m/l cap $user_l_cap exceeds nstreams=$nstreams " *
        "resolving power $stream_l_cap; increase nstreams or lower the cap"))
end

if nstreams < 3 && !allows_sub_rayleigh_streams(rt)
    throw(ArgumentError(
        "nstreams=$nstreams is below the solar/scattering minimum. " *
        "Use at least 3 streams per hemisphere; classical two-stream " *
        "terminology corresponds to one stream per hemisphere."))
end
```

Phase truncation parse:

```julia
truncation =
    if haskey(rt, "truncation")
        rt["truncation"] === nothing ?
            Scattering.NoTruncation() :
            _truncation_from_string(rt["truncation"], FT)
    elseif haskey(rt, "nstreams")
        # New-style config: omission means no phase-function modification.
        Scattering.NoTruncation()
    else
        # Legacy compatibility only: old configs omitted `truncation` but
        # expected δBGE(l_trunc, Δ_angle). Deprecate this behavior.
        Scattering.δBGE{FT}(user_l_cap, Δ_angle)
    end
```

Important distinction:

- `stream_l_cap` / `user_l_cap` is the **projection cap**: coefficients above
  this order cannot be represented by the chosen streams.
- `truncation` is a **phase-function transform**: `δBGE` or delta-M modifies
  the phase function and optical properties so the unresolved forward peak is
  represented consistently.
- `NoTruncation()` means no phase-function transform. It does not magically
  make a 1000-term coarse aerosol fit into 13 streams.

If `NoTruncation()` is active and `length(greek.β)-1 > user_l_cap`, do not
silently discard the tail. The near-term behavior should be at least a warning,
and preferably an error for aerosol/coarse-particle cases:

```julia
if truncation isa Scattering.NoTruncation && phase_lmax > user_l_cap
    throw(ArgumentError(
        "phase coefficients require l=$phase_lmax but nstreams=$nstreams " *
        "can only resolve l=$user_l_cap. Increase nstreams or set truncation."))
end
```

For compatibility, `max_m` should warn once:

```julia
Base.depwarn(
    "`radiative_transfer.max_m` is a historical Fourier-count field. " *
    "Use `m_max` for maximum Fourier order or omit it for auto.",
    :parameters_from_dict)
```

### `vSmartMOM_Parameters`

Change the input parameter struct from:

```julia
max_m::Integer
l_trunc::Integer
```

to:

```julia
m_max_override::Union{Nothing, Int}
nstreams::Int
stream_l_cap::Int
legacy_l_cap_override::Union{Nothing, Int}
truncation::AbstractTruncationType
convergence::FourierConvergenceConfig
```

To reduce churn, an intermediate implementation can keep a deprecated
`max_m::Integer` field and add new fields beside it. The parser would fill
both:

```julia
legacy_max_m_count = isnothing(m_max_override) ? 0 : m_max_override + 1
```

but new code should consume `m_max_override`.

### Model Creation

`model_from_parameters` should be responsible for deriving all band-specific
resolution quantities after optical properties and surfaces are known:

```julia
resolution = derive_resolution(params, optics, surfaces, sources)
```

where:

```julia
struct BandResolution
    l_max::Vector{Int}
    m_max::Vector{Int}
    n_fourier_moments::Vector{Int}
    nstreams::Vector{Int}
end
```

Near-term model creation flow:

1. Parse config into normalized `nstreams`, `stream_l_cap`, `m_max_override`,
   `truncation`, and convergence settings.
2. Build geometry and preliminary quadrature using `stream_l_cap` to preserve
   current one-global-quadrature behavior.
3. Compute Rayleigh/aerosol/surface/source objects.
4. Derive per-band `m_max` with component traits.
5. Store `m_max_bands`, `l_max_bands`, and `nstreams_bands` in `SolverConfig`.
6. Keep one global `quad_points` initially.

Follow-up model creation flow:

1. Derive per-band resolution before allocating RT workspaces.
2. Build `quad_points_bands::Vector{QuadPoints}`.
3. Let each band use its own `NquadN` and workspace dimensions.

This staged approach makes today's change mostly parser/model metadata plus
loop bounds. It avoids forcing a full per-band workspace refactor immediately.

## Internal Model Fields

Separate the quantities that are currently mixed together:

```julia
struct QuadPoints{FT}
    # existing fields...
    Nstreams::Int  # base weighted half-space streams
    Nquad::Int     # augmented nodes, including zero-weight VZA/SZA
end

struct SolverConfig
    nstreams::Int
    stream_l_cap::Int
    m_max_override::Union{Nothing, Int}
    truncation::AbstractTruncationType
    m_max_bands::Vector{Int}
    l_max_bands::Vector{Int}
    nstreams_bands::Vector{Int}
end
```

`Nquad` remains a kernel sizing detail. `Nstreams` is the solver-resolution
quantity used for caps, summaries, and VLIDORT comparisons.

## Component `m_max` Traits

Define exact or effective Fourier support per component:

```julia
component_m_max(::RayleighScattering, ctx) = 2
component_m_max(::LambertianSurfaceScalar, ctx) = 0
component_m_max(::LambertianSurfaceSpectrum, ctx) = 0
component_m_max(::SurfaceSIF, ctx) = 0

component_m_max(a::AerosolOptics, ctx) =
    min(effective_lmax(a.greek_coefs, ctx), ctx.user_l_cap)

component_m_max(::CoxMunkSurface, ctx) = ctx.user_l_cap  # near-term cap
component_m_max(::rpvSurfaceScalar, ctx) = ctx.user_l_cap
component_m_max(::RossLiSurfaceScalar, ctx) = ctx.user_l_cap
component_m_max(::CanopySurface, ctx) = ctx.user_l_cap
```

For Rayleigh, this makes the current implicit behavior explicit. Rayleigh
Greek coefficients have length 3, so terms above `m = 2` are exactly zero.

For aerosols, validate before capping:

```julia
phase_lmax(greek) = length(greek.β) - 1

function validate_phase_resolution!(greek, ctx)
    if ctx.truncation isa Scattering.NoTruncation &&
       phase_lmax(greek) > ctx.user_l_cap
        throw(ArgumentError(
            "NoTruncation cannot fit l=$(phase_lmax(greek)) into " *
            "user_l_cap=$(ctx.user_l_cap); increase nstreams or use δBGE."))
    end
end

effective_lmax(greek, ctx) =
    ctx.truncation isa Scattering.NoTruncation ?
        phase_lmax(greek) :
        min(phase_lmax(greek), ctx.user_l_cap)
```

Later, `effective_lmax` can use coefficient-tail criteria or phase-function
error estimates. Today, length-and-ceiling is enough.

## Per-Band Aggregation

For each spectral band:

```julia
stream_l_cap = 2*params.nstreams - 1
user_l_cap = something(params.m_max_override,
                       params.legacy_l_cap_override,
                       stream_l_cap)

ctx = (; stream_l_cap,
       user_l_cap,
       truncation = params.truncation)

l_required = maximum(component_m_max(component, ctx)
                     for component in active_components_for_band)

m_max = min(l_required, user_l_cap, stream_l_cap)

if params.m_max_override !== nothing
    m_max = min(m_max, params.m_max_override)
end

l_max = m_max
Nstreams = (l_max + 2) ÷ 2
n_fourier_moments = m_max + 1
```

Then build per-band quadrature from `l_max`, not from the global ceiling when
band-dependent savings are enabled.

Near-term conservative path:

1. Compute per-band `m_max_bands`.
2. Keep one global `QuadPoints` for the maximum required `l_max` across bands.
3. Use per-band Fourier loop bounds immediately.

Second step:

1. Store per-band `QuadPoints`.
2. Use smaller `NquadN` workspaces for bands with lower `l_max`.
3. This is where most compute savings appear.

## Fourier Loop Semantics

Rename the current public/config concept:

```julia
# old
for m = 0:max_m - 1

# new
for m = 0:m_max
```

Internally, use names literally:

- `m_max`: maximum Fourier order included.
- `n_fourier_moments`: count, equal to `m_max + 1`.

Compatibility bridge:

- Keep reading old `max_m` for one release.
- Interpret old `max_m` as the historical count.
- Convert with `m_max_override = max_m - 1`.
- Emit a deprecation warning recommending `m_max`.

## Per-Component Skipping Inside the Loop

The outer loop should run to the band maximum. Individual components should
skip after their own support ends:

```julia
for m = 0:m_max_band
    include_rayleigh = m <= component_m_max(rayleigh, ctx)
    include_aerosol  = m <= component_m_max(aerosol, ctx)
    include_surface  = m <= component_m_max(surface, ctx)
end
```

This matters when, for example, Cox-Munk needs high Fourier orders but
Rayleigh is exactly zero above `m = 2`.

Near-term implementation can simply pass zero `Z` for skipped components or
avoid calling `compute_Z_moments` for them. The latter is preferred because it
avoids unnecessary Legendre setup and allocation.

## Optional Fourier Convergence

Convergence should be an optional layer on top of deterministic component
caps, not a replacement for them.

Suggested config:

```yaml
radiative_transfer:
  convergence:
    enabled: false
    rtol: 1.0e-5
    atol: 1.0e-8
    min_m: 2
    consecutive: 2
```

Suggested criterion:

```julia
converged =
    m >= min_m &&
    maxabs(ΔR_m) <= atol + rtol * max(maxabs(R_sum_I), intensity_floor) &&
    satisfied_for_consecutive_moments
```

Notes:

- Do not use Q/U/V relative error against their own near-zero values.
- Scale polarized convergence by intensity plus absolute tolerances.
- Keep convergence disabled by default for retrieval/Jacobian workflows until
  reproducibility and derivative behavior are validated.
- Enable it first for exploratory forward-only runs.

## Source Composition

The source abstraction should allow multiple sources. The current source
design already supports this through `SourceSet`:

```julia
SolarBeam() + SurfaceSIF()
```

Config should eventually expose a list of sources:

```yaml
radiative_transfer:
  sources:
    - SolarBeam()
    - SurfaceSIF()
```

Default remains `SolarBeam()` for backward compatibility. Source support in
the Fourier plan is just another `component_m_max` contributor:

- direct solar beam follows atmospheric/surface scattering needs;
- SIF isotropic surface source contributes only `m = 0`;
- thermal/lunar/lidar sources can add their own traits later.

## Implementation Plan We Can Start Today

### Step 1: Rename Semantics Internally

- Add `m_max_bands` to `SolverConfig`.
- Keep `max_m_bands` temporarily as a compatibility alias if needed.
- Change forward and linearized loops to `for m = 0:m_max_bands[iBand]`.
- Update summaries to print `m_max`, `n_fourier_moments`, `l_max`, and
  `Nstreams`.

### Step 2: Add `Nstreams`

- Add `Nstreams` to `QuadPoints`.
- In `GaussLegQuad`, compute `Nstreams = (l_cap + 2) ÷ 2` before appending
  zero-weight SZA/VZA nodes.
- Audit `RadauQuad`; align it with the same convention or document why it
  intentionally differs.
- Keep `Nquad` as the augmented node count.

### Step 3: Add Component Traits

- Implement `component_m_max` methods for Rayleigh, aerosols, Lambertian
  surfaces, Cox-Munk, RPV, Ross-Li, canopy, and sources.
- Aggregate per-band `m_max_bands`.
- Replace scalar `get_max_m(model)` consumption with per-band resolution.

### Step 4: Separate Projection Cap From Phase Truncation

- For Rayleigh-only/Lambertian cases, derive `m_max = 2` even when
  `nstreams = 13`.
- For Lambertian + no scattering, derive `m_max = 0`.
- For aerosols with `NoTruncation()`, require `length(β)-1 <= user_l_cap` or
  tell the user to increase `nstreams` / choose `δBGE`.
- For aerosols with `δBGE`, transform the phase function to the chosen cap;
  do not treat this as equivalent to dropping high-order Greek coefficients.
- For Cox-Munk/canopy/hotspot BRDFs, use `user_l_cap` as the near-term cap.

### Step 5: Config Migration

- Add `nstreams` as the preferred schema field and examples default.
- Document that `nstreams` means streams per hemisphere, not total streams and
  not classical two-stream terminology.
- Validate `nstreams >= 3` unless the scene is explicitly thermal-only or
  otherwise source-only.
- Add `m_max: null` and `truncation: null` to YAML examples.
- Keep old `max_m` accepted for one release as historical count.
- Keep old `l_trunc` accepted as a legacy exact cap; convert it to
  `nstreams = (l_trunc + 2) ÷ 2` and `m_max = l_trunc`.
- Update default YAML/TOML examples to omit `max_m` and use `nstreams`.
- Update VLIDORT fixtures explicitly:
  - `NSTREAMS = 8`
  - `l_max = m_max = 15`
  - `n_fourier_moments = 16`

### Step 6: Optional Convergence

- Add a disabled-by-default convergence config struct.
- Implement convergence only for forward RT first.
- Add tests that fixed-loop and convergence-disabled results are unchanged.
- Add one forward-only smoke test where convergence exits early for a simple
  Lambertian/Rayleigh scene.

### Step 7: Per-Band Quadrature

- First use per-band Fourier bounds with the existing global quadrature.
- Then add per-band `QuadPoints` and workspaces.
- This is the larger performance change and can be staged after semantics are
  stable.

## Acceptance Criteria

1. Old configs still run with a deprecation warning for `max_m`.
2. VLIDORT solar tester maps cleanly:
   `NSTREAMS=8 -> l_max=m_max=15 -> m=0:15`.
3. Rayleigh over Lambertian uses only `m=0:2` under default `nstreams=13`.
4. Lambertian-only/no-scattering runs use only `m=0`.
5. Coarse aerosol with many coefficients either uses an explicit phase
   truncation method (`δBGE`) or errors under `NoTruncation()` instead of
   silently dropping Greek coefficients.
6. Solar/scattering configs reject `nstreams < 3` with an error that says
   streams are per hemisphere.
7. `Nstreams` and `Nquad` are visibly distinct in summaries and diagnostics.
8. Forward and linearized paths use the same per-band `m_max` decision.

## Open Questions For Review

1. Should `truncation: null` error or warn when aerosol coefficients exceed
   stream resolving power?
2. Should we expose an explicit unsafe/plain cutoff mode for expert debugging,
   distinct from `NoTruncation()` and `δBGE`?
3. For Cox-Munk, is `user_l_cap = 2*nstreams-1` acceptable as the temporary
   cap, or should wind speed set a larger/smaller effective cap?
4. What is the right source / RT-mode flag for allowing `nstreams < 3`
   (`nstreams = 1` included) in thermal-only / source-only scenes?
5. Should convergence ever be allowed in linearized retrieval runs, or should
   it remain forward-only?
6. For per-band quadrature, do we want separate workspaces per band
   immediately, or first land the semantic cleanup with global quadrature?
7. Should old `max_m` be deprecated immediately, or kept indefinitely as
   `n_fourier_moments` for YAML compatibility?
