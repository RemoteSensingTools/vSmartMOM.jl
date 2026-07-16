# Release Notes and Migration

**For:** users moving to the registration-ready vSmartMOM 2.0 interface.

**Next:** [Quick Start](quickstart.md), [Configure a Scene](IO/Overview.md), [Compute Jacobians](jacobians.md).

This page summarizes the user-visible changes in the 2.0 line. It is written as
a migration guide, not as a complete git history.

## Unreleased — observer height selection and multi-level radiances

`geometry.obs_alt` is now a geometric height selection in **km above BOA**,
not a pressure. Scalar and vector forms distinguish a single observing level
from the historical full-column output:

- `0` returns BOA downwelling only; an interior scalar `H` returns both
  directions at `H`; a scalar at or above the model TOA returns TOA
  upwelling only.
- `[0]` returns TOA upwelling and BOA downwelling.
- `[H1, H2]` returns both directions at each strict-interior height, and
  `[0, H]` adds the TOA and BOA endpoints.

Existing full-column scenes that used scalar `0` or a pressure-like value such
as `1000.0` should migrate to `obs_alt: [0]`; those scalars now mean BOA-only
and a height above TOA, respectively.

An interior height becomes an exact atmospheric interface. vSmartMOM
interpolates layer-centred input profiles in log pressure, recomputes derived
vertical fields, and builds the optical state on the reframed grid. Profile
reduction preserves all requested interfaces; if `K` distinct interior
heights require more than the requested reduced layering, the solver raises
the count to `K + 1` and warns.

Humidity and vector VMR inputs may now use either `N` layer-center values or
`N+1` pressure-interface values; interface fields are normalized in log
pressure before height insertion. The GEOS-Chem reader also converts
`Met_SPHU` from its declared g/kg units to the internal kg/kg mass fraction.

Forward `rt_run` now returns an `ObserverRTResult` with named `toa`, `boa`, and
`levels` fields. Historical tuple iteration remains available, so scenes with
`obs_alt: [0]` continue to support `R, T = rt_run(model)`. See the
[`geometry` schema](IO/Schema/geometry.md) for the complete convention and
output fields. Strict-interior outputs currently use the full multiple-scatter
forward driver; the analytic-linearized and single-scatter drivers reject
interior requests explicitly while continuing to honor endpoint-only choices.

## v2.1.0 — Fourier / stream resolution + source-term refactor

> **vSmartMOM** = **V**ector **S**imulated **M**easurements of the **A**tmosphere using **R**adiative **T**ransfer based on the **M**atrix **O**perator **M**ethod.

This is the package version that ships the schema generation we
internally called `v0.7` (Fourier / Stream Resolution refactor) plus
the `v0.6` source-term abstraction. *Treat it as a breaking schema
migration*: `max_m` / `l_trunc` still parse for tactical reasons but
are no longer the recommended idiom in any in-tree YAML — every
example uses `nstreams` + `truncation: auto` now.

The new public input is **`nstreams`** (weighted streams per
hemisphere). The schema-generation v0.7 label persists in some docs
(`stream_l_cap = 2·nstreams - 1` and the per-band trait dispatch);
the package version on the registry is v2.1.0.

### What changed for users

```yaml
# Before (still works — legacy schema)
radiative_transfer:
  quadrature_type: GaussLegQuad()
  max_m: 14
  l_trunc: 25

# After (new schema, recommended)
radiative_transfer:
  nstreams: 8            # public contract: stream_l_cap = 2·N - 1
  truncation: auto       # NoTruncation if phase fits, δBGE otherwise
  # quadrature_type omitted ⇒ GaussLegQuad() (Sanghavi: cheaper +
  # 5–50× more accurate per stream than RadauQuad on Rayleigh)
```

### Highlights

- **`nstreams ≥ 3`** is enforced — Rayleigh contributes through m=2.
- **`truncation: auto`** mirrors VLIDORT's `DO_DELTAM_SCALING`: the
  resolver picks `NoTruncation()` when the phase function fits within
  `stream_l_cap`, `δBGE(N, Δ_angle)` otherwise. Logs the choice via
  `@info` so the user always sees what was applied.
- **`GaussLegQuad()`** is the default quadrature when `quadrature_type`
  is omitted. `RadauQuad()` remains an option but is documented as
  expert/legacy. (See the natraj benchmark in
  [`docs/src/pages/benchmarks.md`](benchmarks.md) for the per-scheme
  accuracy comparison.)
- **JSON Schema** at
  [`schemas/vsmartmom-parameters.schema.json`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/schemas/vsmartmom-parameters.schema.json)
  describes the v0.7 schema. Wired to TOML via the repo's
  `.taplo.toml`; YAML editors honour it via
  `# yaml-language-server: $schema=...`. Setup recipe in
  [`Schema.md`](IO/Schema.md#editor-support--autocomplete--inline-validation).
- **Per-block schema docs** under
  [`docs/src/pages/IO/Schema/`](IO/Schema.md) — one page per top-level
  YAML/TOML block.
- **Gas absorption engine swap (externalAbsorption branch).**
  `model_from_parameters` now computes LBL gas absorption via
  [AtmosphericAbsorption.jl](https://github.com/RemoteSensingTools/AtmosphericAbsorption.jl)
  (TIPS-2021 partition sums) rather than the internal `Absorption` module (TIPS-2017).
  Users with high-temperature retrievals (> ~500 K) may see small numerical
  differences in the resulting optical depths; results at standard atmospheric
  temperatures are practically unchanged.  The standalone `Absorption` module
  API (`read_hitran` / `make_hitran_model` / `absorption_cross_section`)
  remains available unchanged for direct σ(ν, T, p) queries.

### Internal cleanup

- `Nstreams` and `m_max_bands` are first-class fields on `QuadPoints`
  and `SolverConfig`. Forward and lin RT paths share a single
  `_derive_m_max_bands` helper, fixing a latent precedence-bug in the
  lin path at even `l_max`.
- Per-component `component_m_max(c, ctx)` traits drive the per-band
  Fourier loop bound. Cox-Munk / RPV / Ross-Li / canopy now run to
  their full `user_l_cap` instead of the silently-half-truncated
  count-aggregator output.

### Migration guidance

- **YAMLs that explicitly set `max_m` / `l_trunc`** keep working
  byte-identically — they hit the legacy parsing branch.
- **YAMLs that hand-set `max_m` to upper-bound a Cox-Munk loop**
  (e.g. `config/ocean_coxmunk.yaml`) get a numeric change: the
  Fourier loop now runs to the user-set order rather than the
  half-truncated aggregator output. The new behavior matches the
  user's explicit request.
- **lin Jacobians** may shift by one Fourier moment for configs with
  even `l_max` (intentional fix to a latent precedence bug;
  documented in the Phase B commit).

See [`docs/dev_notes/fourier_stream_resolution_plan.md`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/docs/dev_notes/fourier_stream_resolution_plan.md)
for the design rationale.

## Platform Support

vSmartMOM 2.0 supports Julia 1.10 and newer Julia 1.x releases listed in
`Project.toml`. CUDA is an optional weak dependency, and the package includes
compatibility with CUDA.jl 6.

CPU remains the default portable path:

```julia
params = read_parameters(joinpath(pkgdir(vSmartMOM), "config", "quickstart.yaml"))
params.architecture = CPU()
model = model_from_parameters(params)
R, T = rt_run(model)
```

GPU runs use the same high-level flow, but set `params.architecture = GPU()` and
load CUDA in the environment that runs the model.

## Model Hierarchy

`model_from_parameters(params)` now returns an `RTModel{ARCH, FT}`. The model
organizes solver state into named sub-objects:

- `solver`: polarization, quadrature, truncation, and numerical settings.
- `geometry`: solar/viewing geometry and resolved observer-output interfaces.
- `quad_points`: quadrature abscissae and weights.
- `atmosphere`: atmospheric profile and spectral bands.
- `optics`: Rayleigh, aerosol, gas, and linearized optical-property state.
- `surfaces`: one surface model per spectral band.
- `sources`: source-term scene (v0.6 — see below). Default
  `SolarBeam()` preserves historical unit-Stokes-I behavior.

The old flat `vSmartMOM_Model` container has been replaced. Compatibility
getters still expose many legacy fields, but new code should prefer the
hierarchical fields and accessor functions documented in the
[Library](api_reference.md).

## v0.6 Source-Term Refactor

Source terms (solar, surface fluorescence, blackbody, future thermal /
lidar / lunar) are now first-class composable types under
`AbstractSource`. The MOM solver is mathematically affine — optical
properties define the operator `A`, sources define the additive RHS `b` —
and v0.6 makes that explicit at the type level.

```julia
# Default (unchanged behavior — unit Stokes-I solar beam):
R, T = rt_run(model)

# Explicit solar with custom irradiance (mW · m⁻² · cm⁻¹):
R, T = rt_run(model; sources = SolarBeam(F₀ = solar_spec))

# Solar + surface fluorescence:
R, T = rt_run(model; sources = SolarBeam() + SurfaceSIF(SIF₀ = sif_spec))

# Thermal-IR / Carbon-I-style scene (1500 K blackbody at 2-2.4 µm):
spec_band = collect(model.atmosphere.spec_bands[1])
R, T = rt_run(model; sources = BlackbodySource(1500, spec_band))
```

Highlights:

- `RTModel.sources::AbstractSource` is a new top-level field. Set it
  with `model_from_parameters(params; sources = …)` or override per-call
  with the `sources=` kwarg on `rt_run` / `rt_run_lin` / `rt_run_ss`.
- `SolarBeam` (FT-deferred — the model's float type drives precision via
  `prepare_source`), `SurfaceSIF`, `BlackbodySource` (constructor sugar
  that returns a `SolarBeam` with Planck-derived F₀), `NoSource`,
  `SourceSet` (compile-time-unrolled tuple) ship as built-ins.
- Composition via `+`: `SolarBeam() + SurfaceSIF()` ⇒
  `SourceSet((SolarBeam, SurfaceSIF))`. `NoSource` is the identity.
- Multiple-dispatch entry points `contribute!` (forward) and
  `source_tangent!` (linearized) replace the legacy `if SFI` branching
  across the elastic linearized path. Surfaces use a parallel
  `surface_source_contribute!(prepared_sources, surface, …)` double-
  dispatch table for surface-side contributions (SIF on Lambertian
  surfaces today; non-Lambertian dispatches to a no-op).
- The Sanghavi 2014 App. C analytic tangents and the Bug-22 beam-
  attenuation chain-rule fix are unchanged inside
  `get_elem_rt_SFI_fused!`; they're now reached via
  `source_tangent!(::PreparedSolarBeam, …)` — the named hand-written
  linearization seam.
- AD seam at `prepare_source(::AbstractSource, FT, pol_n, nSpec, arr_type)`:
  user-parameter space (potentially `Dual`) above; kernel space (plain
  `FT`) below. A future `prepare_source_with_tangent` will mirror this
  signature and return source-parameter Jacobians; v0.6 reserves the
  `ForwardDiffSourceJacobian` AD-mode trait for that path.
- Source-unit convention adopted: `F₀` and `SIF₀` in
  **mW · m⁻² · cm⁻¹**, `rt_run` outputs in
  **mW · m⁻² · sr⁻¹ · cm⁻¹**. `BlackbodySource` defaults to
  `factor = π` (Lambertian-disk → hemisphere irradiance) so its `F₀` is
  comparable to a unit `SolarBeam(F₀ = 1)` baseline.

Legacy compatibility: the `inject_surface_SIF!(brdf, …, SIF₀, arch)` path
that consumes `RS_type.SIF₀` still runs in parallel for back-compat with
existing `rs.SIF₀ = …; rt_run_test_ss(rs, model, 1)` test patterns. If a
user supplies BOTH `model.sources` containing a `SurfaceSIF` AND sets
`RS_type.SIF₀`, the surface SIF is double-counted; choose one API at a
time. A future PR will retire `RS_type.F₀` / `RS_type.SIF₀` ownership in
favour of the dispatch system.

Full architecture writeup: [Sources](extending/sources.md). Test coverage
is in `test/test_sources.jl` (Phase 1 → 5b regression, 86/86 pass).

## Parameter Loading

Use `read_parameters` as the public entry point:

```julia
params = read_parameters("scene.yaml")
params = read_parameters(Dict("architecture" => "CPU"))
params = read_parameters(GeosChemSource("GEOSChem.Custom.nc4", 1, 1, 1))
```

The explicit lower-level aliases remain available:

- `parameters_from_file(path)` for extension-dispatched YAML/TOML files.
- `parameters_from_dict(dict)` for in-memory dictionaries.
- `parameters_from_source(src)` for `IOSource` values.
- `parameters_from_yaml(path)` for YAML/YML files only.

`parameters_from_yaml("scene.toml")` now raises `ArgumentError`; use
`read_parameters("scene.toml")` or `parameters_from_file("scene.toml")`.

## Jacobians

Linearized RT is exposed through `LinMode`, `model_from_parameters_lin`, and
`rt_run_lin`:

```julia
params = read_parameters(joinpath(pkgdir(vSmartMOM), "test", "test_parameters", "JacobianTestFast.yaml"))
model, lin_model = model_from_parameters_lin(params)
NAer  = length(params.scattering_params.rt_aerosols)
NGas  = size(lin_model.τ̇_abs[1], 1)
NSurf = 1
R, T, R_jac, T_jac = rt_run_lin(model, lin_model, NAer, NGas, NSurf)
```

Use `ParameterLayout` to map Jacobian columns to aerosol, gas, surface, and
canopy parameter blocks. Raman/inelastic linearization is not complete in this
release line; keep linearized scenes on the elastic path unless a feature branch
or later release states otherwise.

## Test and Data Policy

The registration-oriented test suite avoids machine-local data assumptions.
RAMI remains useful as an external canopy validation workflow, but it is not a
unit-test fixture in the default test environment.

The EMIT/MODTRAN comparison scripts remain in the repository as benchmark and
analysis utilities, but their execution tests are disabled in CI because they
depend on external LUTs and machine-local scene inputs. The Phase 6 test suite
still parses those scripts so ordinary import and syntax drift is caught.

SIF helper functions are exported because the data loaders exist, but the
end-to-end SIF workflow is kept under [Experimental Helpers](api/experimental.md)
until the data-fixture and product policy are settled.

## Documentation and Citation

The manual has been reorganized around task pages first and a module-grouped
[Library](api_reference.md) second. The docs build runs with
`checkdocs = :exports` and without `warnonly`, so exported public symbols and
cross references are checked as part of the release gate.

The repository now ships a root `CITATION.bib` with the JOSS software paper and
the core method references. See [References](vSmartMOM/References.md) for the
human-readable citation guide.

## Correctness Fixes

The following bugs were silently wrong in earlier releases; all are fixed in
the current codebase.

- **δBGE forward-cone exclusion (`Δ_angle`) now applied in production.**
  The `Δ_angle` parameter of `δBGE(N, Δ_angle)` controls how much of the
  forward-scattering peak is excluded before fitting the truncation coefficient.
  Previously, `Δ_angle` was accepted in the constructor but silently ignored
  during the production truncation fit inside `model_from_parameters`.
  It is now applied correctly.  The **default `Δ_angle = 0` is unchanged**,
  so most users see no numerical change.  Users who explicitly set `Δ_angle > 0`
  in their configs should re-validate their results against this release.

- **`LambertianSurfaceSpectrum` surface layer now works.**
  Calling `rt_run` on a model with a `LambertianSurfaceSpectrum` surface
  previously raised a `MethodError` because `create_surface_layer!` was not
  implemented for that type.  The missing method has been added.

- **Canopy soil-albedo Jacobian axis fix.**
  The canopy linearization (`create_surface_layer_lin!`) wrote the soil-albedo
  derivative along the wrong array axis, producing incorrect Jacobians for
  canopy scenes with soil Jacobians enabled.  The axis index is now correct.

- **Wavelength-path wing-cutoff window fix (standalone Absorption API).**
  When computing absorption cross-sections with the standalone `Absorption`
  module, the per-line wing-cutoff window was computed in wavelength space
  rather than wavenumber space, producing asymmetric cutoffs for lines near
  band edges at longer wavelengths.  The window is now applied consistently
  in wavenumber space.

- **VS Raman call-chain and atomic-mass fix.**
  The vibrational Raman scattering (VS) call chain contained a unit-inconsistency
  in the atomic-mass factor used to compute the Raman shift.  Corrected values
  will produce numerically different VS Raman spectra.  Rotational Raman (RRS)
  is unaffected.

- **Mie Dₙ recurrence: Dual-number-safe initialization.**
  The downward Dₙ recurrence used in Mie coefficient computation started from
  a hard-coded `Float64` initial value, causing a type-promotion error when
  the refractive index was a `ForwardDiff.Dual` number (AD path).  The
  initializer is now type-generic, restoring correct analytic derivatives for
  Mie-sensitive parameters via `compute_aerosol_optical_properties(LinMode(), …)`.

## Dependencies

- **LogExpFunctions** compat widened to `0.3` (was `~0.2`) to resolve downstream
  conflicts when using vSmartMOM alongside other packages.
- **Parameters.jl** compat widened to allow `0.13` in addition to `0.12`, reducing
  resolver friction for users that also depend on Turing.jl or Optim.jl.

## Performance

- **GPU Mie dispatch** — CUDA architecture now runs Mie coefficient computation
  on-GPU via `make_mie_model(...; architecture = GPU())`.  Measured speedup:
  **10.7–12.9×** vs CPU for typical aerosol grids (1 500 quadrature points,
  NAI2 decomposition).  `NativeFloat64` is the default precision policy;
  `DSEmulated` is available for extended-precision Float32 accumulation.
  Metal and other architectures fall back to CPU Mie automatically.

- **Flat-Z zero-copy path** — Rayleigh-only and analytic-phase-function layers
  share a single pre-computed Z matrix across all Fourier moments rather than
  re-copying it per moment, reducing memory traffic in the Fourier loop.

- **m-invariant Fourier-loop cache** — Optical properties that do not vary
  with Fourier moment `m` (Rayleigh, absorption) are computed once per layer
  and cached across moments, cutting per-moment setup cost for scenes with
  many moments.

- **GPU sync stripping** — Redundant device-synchronization calls between
  kernel launches have been removed.  Measured improvement: **~19.5%** faster
  GPU forward RT on the O₂-A band test case.

- **`BatchContext` batch API** — New `BatchContext` / `update_model!` /
  `update_aerosol_loading!` / `update_aerosol_microphysics!` API enables
  efficient multi-scene loops (ensemble retrievals, parameter sweeps) by
  caching Mie, Fourier decomposition, and HITRAN parsing across calls.
  See the [CoreRT API](api/core_rt.md#batch-scene-processing) for details.

## Known WIP Areas

- Aerosol scene input is documented, but the high-level API is still being
  stabilized.
- Internal export cleanup is intentionally deferred; this release does not
  tighten public API surface area at the same time as the registration cutover.
- Unified offline source-function integration and thermal emission are design
  topics, not implemented user-facing features in this release line.

## Future Developments

- Mixed automatic differentiation and hand-coded linearization.
- Coupled ocean-atmosphere radiative transfer.
- More line-shape and line-mixing functionality, including expansion to more
  spectral databases.
- Raman performance and memory-pressure improvements.
- Correlated-k and other spectral dimension-reduction methods.
