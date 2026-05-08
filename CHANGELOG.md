# Changelog

## v2.1.0 — Fourier-stream resolution + source-term refactor

**Schema migration release.** Treat as breaking for downstream code
that hard-codes `params.max_m` / `params.l_trunc` field reads, raw
access to `SolarBeam`'s `F₀` on `RS_type`, or `Δ_angle = 2.0` as the
implicit default.  The legacy YAML parser branch (`max_m`/`l_trunc`)
is preserved for tactical reuse, but every in-tree YAML has been
migrated to `nstreams` + `truncation: auto`; users following examples
will land on the v0.7 schema.

### What's new since v2.0.0 — at-a-glance

If you last looked at `vSmartMOM.jl` on `main`, this is what changed.
Most items have their own subsection further down, or a dedicated
docs page in `docs/src/pages/`.

| Theme | Summary |
|-------|---------|
| **Source-term abstraction (v0.6)** | First-class `AbstractSource` types — `SolarBeam`, `BlackbodySource`, `SurfaceSIF`, `NoSource`, composed via `+` into `SourceSet`. Replaces the implicit `RS_type.F₀` ownership. Forward/linearized RT now flow through `prepare_source` → `contribute!` / `source_tangent!` dispatchers. AD seam at `prepare_source` (user-parameter space ⇄ kernel space). See `docs/src/pages/extending/sources.md`. |
| **Fourier / stream resolution schema (v0.7)** | `nstreams` is the primary user input; `stream_l_cap = 2·nstreams - 1` regardless of quadrature. `truncation: auto` (conservative in v2.1, see below) and `truncation: "δBGE(L, Δ)"` are the recommended idioms. Per-component `component_m_max(c, ctx)` traits — Cox-Munk / RPV / RossLi / canopy now run to their full `user_l_cap` instead of being silently half-truncated. |
| **VLIDORT 2.8.3 baseline suite** | New `test/vlidort_baseline/` with Siewert 2000 Problem IIA (Stokes-I), and `solar_tester` scalar (Task 1) + vector (Task 1, IQU) reference data committed in tree. No PyVLIDORT / Fortran runtime needed for regression. |
| **Stokes IQ polarization** | New `Stokes_IQ()` polarization type — 2-channel I/Q, between scalar I and IQU, useful for instruments with linear polarization but no Q/U separation. Forward, single-scatter, surface BRDFs (incl. Cox-Munk SS wind Jacobian), and Greek single-scatter all covered. |
| **Standalone exact single-scatter** | `StandaloneSS` driver (`run_exact_ss`, `ExactSSConfig`, `SSGeometry`) for high-precision SS-only computations independent of the adding-doubling solver. RTModel SS adapter included. |
| **Aerosols module** | Multi-moment + TOMAS-15 aerosol size-distribution support; analytic phase-function aerosols (Henyey-Greenstein); `compute_aerosol_optical_properties` exported; vector phase + surface BRDF chain rules wired through. |
| **Canopy surface (PROSAIL/4SAIL)** | `CanopySurface` BRDF integrated into the standard `rt_run` pipeline (legacy `rt_run_canopy` removed). Canopy block in YAML schema; per-band soil + leaf reflectance/transmittance. |
| **HITRAN edition selection** | `set_hitran_edition!`, `get_hitran_edition`, `available_hitran_editions`, `fetch_hitran`, `fetch_hitran_by_ids` — fetches from hitran.org with SHA-256 provenance. Defaults to the legacy 2016 artifacts; switch with one call. |
| **Numerics tuning surface** | `RTNumericalParameters` carries `dτ_max_threshold`, `dτ_min_floor`, `blas_threads` knobs that previously lived as hard-coded `0.001` / `1024·eps(FT)`. YAML key `radiative_transfer.numerics` exposes them. |
| **Documentation rewrite** | DocumenterVitepress migration, narrative-thread restructure (Concepts arc), per-block schema docs under `docs/src/pages/IO/Schema/`, JSON-Schema autocomplete for TOML/YAML, auto-generated benchmarks page (Siewert IIA et al). README + landing page now expand the **vSmartMOM** acronym (Vector Simulated Measurements of the Atmosphere using Radiative Transfer based on the Matrix Operator Method). |
| **`RTModel` hierarchy (v2.0 finalised)** | `model_from_parameters` returns a hierarchical `RTModel{ARCH, FT}` with physics-based sub-structs (`SolverConfig`, `Atmosphere`, `Optics`, `Surfaces`, `Sources`). `Base.getproperty` shim provides backward-compat field access. |
| **CI hardening** | Test matrix now Julia 1.10 / 1.11 / 1.12; `actions/checkout@v4`; `julia-actions/cache@v2` step; Aqua quality gate (stale deps, ambiguities, piracy, undocumented exports) included in `runtests.jl`. |

### Breaking changes vs v2.0.0

These will require user action if you have downstream code:

- `params.max_m` and `params.l_trunc` are still parsed but no longer
  the canonical knobs.  In-tree YAMLs all use `nstreams` now.
- `SolverConfig.max_m` / `get_max_m(model)` have been removed in
  favour of `m_max_bands` (per-band order). If you read these,
  replace with `model.solver.m_max_bands[iBand]`.
- `RS_type.F₀` is no longer the public interface for solar irradiance
  — pass `sources = SolarBeam(F₀ = …)` to `model_from_parameters` or
  `rt_run` instead.  Existing code that wrote to `RS_type.F₀`
  directly still works in the legacy code path but the field is
  scheduled for removal in a follow-up.
- `Δ_angle` is now an optional YAML field, default `0.0` (was a
  required field with implicit `2.0` examples).  Set explicitly when
  reproducing VLIDORT `DO_DELTAM_SCALING = 2°` benchmarks.
- `rt_run_canopy` (canopy-only entry point) was removed — call
  `rt_run` with `CanopySurface` as the BRDF.
- `parameters_from_yaml` is now YAML-only; use `parameters_from_file`
  / `read_parameters` for TOML or registry-dispatched inputs.

### Highlights

First-class `AbstractSource` type vocabulary (`SolarBeam`,
`BlackbodySource`, `SurfaceSIF`); user-facing `nstreams` /
`truncation: auto` Fourier-resolution schema; per-model BLAS thread
cap; VLIDORT 2.8.3 baseline validation suite (Siewert 2000 Problem IIA,
solar_tester scalar/vector); `Δ_angle` is now optional and defaults to
`0.0` (legacy `2.0` matches VLIDORT `DO_DELTAM_SCALING` and can still
be set explicitly).

### New user-facing schema (Phase D)

- **`nstreams` is the primary resolution knob.** Public contract:
  `stream_l_cap = 2·nstreams - 1` regardless of `quadrature_type`.
  Minimum 3 (Rayleigh contributes through m=2). Default 8 when
  legacy `max_m`/`l_trunc` are absent — gives `stream_l_cap = 15`,
  which fits typical aerosol δ-fit truncations and gives ~4-decimal
  Rayleigh convergence at standard atmospheric depth.
- **`truncation: auto`** — VLIDORT-`DO_DELTAM_SCALING`-style.
  `NoTruncation()` for Rayleigh-only scenes; `δBGE(stream_l_cap,
  Δ_angle)` once aerosols are configured. The per-band Mie loop
  applies δBGE only when `length(greek.β) > stream_l_cap`, so a
  small-particle Mie series shorter than the cap automatically
  falls back to `NoTruncation()` — no crash, no user action
  needed. Logs the chosen branch via `@info`.
- **`quadrature_type` defaults to `GaussLegQuad()`** when omitted
  (Sanghavi: 5–50× more accurate per stream than `RadauQuad` on
  Rayleigh-only). `RadauQuad` remains supported but is documented
  as expert/legacy.
- **JSON Schema updated** —
  `schemas/vsmartmom-parameters.schema.json` covers the new fields.
  The repo's `.taplo.toml` already wires it to TOML configs; YAML
  editors honour it via `# yaml-language-server: $schema=...`. See
  `docs/src/pages/IO/Schema.md` for the per-block reference and
  editor-setup recipes.

### Internal architecture (Phases A–C)

- **`Nstreams` field on `QuadPoints`** (Phase A) — count of nonzero
  weights, distinct from the augmented `Nquad` total. Surfaces in
  `Base.show` and `conventions.md` §6.
- **`m_max_bands` order semantics on `SolverConfig`** (Phase B) —
  replaces the scalar `max_m` and the count-form `max_m_bands`.
  Forward and lin RT paths now flow through `_derive_m_max_bands`,
  fixing a latent precedence bug at even `l_max` in the lin path
  (`Int(ceil(l_max+1)/2)` → `Int(ceil((l_max+1)/2))`).
- **`component_m_max(c, ctx)` traits** (Phase C) — per-component
  Fourier support dispatch (Rayleigh → 2, Lambertian → 0,
  CoxMunk/RPV/RossLi/Canopy → `user_l_cap`, SolarBeam → 0,
  SurfaceSIF → 0). Resolves Phase B's regression for Cox-Munk
  forward; the Phase C P2 follow-up clamps trait-derived `m_max`
  against the actual quadrature `Nstreams` cap and includes active
  sources in the aggregation.

### Bit-equality

- All 64 in-tree YAMLs run byte-equal through Phase D — they carry
  explicit `max_m` / `l_trunc` so the parser hits the legacy branch.
- Phase B fixed the lin-only precedence bug at even `l_max`; lin
  Jacobians shift by one Fourier moment for affected configs
  (intentional silent-bug fix).
- Phase C+P2 unified forward and lin Fourier-loop bounds. Cox-Munk
  forward goes from the half-truncated count-aggregator output
  back to its full `user_l_cap` resolution.

### See also

- Plan: `docs/dev_notes/fourier_stream_resolution_plan.md`
- Schema index: `docs/src/pages/IO/Schema.md`

## v2.0.0

### Breaking Changes

- **`RTModel` replaces `vSmartMOM_Model`**: `model_from_parameters()` now returns an `RTModel{ARCH, FT}` with physics-based sub-structs (`SolverConfig`, `Atmosphere`, `Optics`, surfaces). The old flat `vSmartMOM_Model` struct has been removed. A `Base.getproperty` shim provides backward-compatible field access (e.g. `model.τ_abs`, `model.profile`).

- **`RTModelLin` replaces `vSmartMOM_Lin`**: The linearized model container returned by `model_from_parameters(LinMode(), params)` is now named `RTModelLin`.

- **`rt_run_canopy` removed**: Canopy RT is now handled through the standard `rt_run()` pipeline using `CanopySurface` as the surface type. The standalone `rt_run_canopy()` function has been removed.

- **`parameters_from_yaml` is YAML-only**: Use `parameters_from_file` or `read_parameters` for TOML or registry-dispatched inputs. Passing a non-`.yaml`/`.yml` path to `parameters_from_yaml` now raises `ArgumentError`.

### New Features

- **Hierarchical model architecture**: `RTModel` organizes state into `solver`, `geometry`, `quad_points`, `atmosphere`, `optics`, and `surfaces` sub-structs, clearly separating fixed configuration from differentiable state.

- **`ParameterLayout`**: New struct for Jacobian index arithmetic — provides `aerosol_range()`, `gas_range()`, `surface_range()`, `n_total()` accessors.

- **Task-based documentation**: The manual now starts from runnable tasks (quick start, scene configuration, Jacobians, GPU, and extension guides) and keeps module pages as reference material.

- **Module-grouped Library**: Public docstrings are organized by module in a dedicated Library section, with a generated function index for lookup.

- **Repository citation file**: `CITATION.bib` now ships at the repository root with the JOSS software paper and core method references.

- **Cox-Munk ocean surface**: Full polarized BRDF with Fresnel reflection, wind-speed-dependent roughness, whitecap correction, and Smith (1967) shadowing. Linearized variant included.

- **Accessor functions**: `architecture(model)`, `CoreRT.polarization_type(model)`, `CoreRT.float_type(model)`, `CoreRT.n_aerosols(model)`, `get_surface(model, iBand)`, `get_spec_bands(model)`.

### Improvements

- Julia compat is declared for Julia 1.10, 1.11, and 1.12. CUDA remains a weak dependency and now allows CUDA.jl 6.
- `read_parameters` is the preferred public loader for file paths, dictionaries, and `IOSource` values. `parameters_from_file`, `parameters_from_dict`, and `parameters_from_source` remain available as explicit aliases.
- Parser validation now throws `ArgumentError` instead of relying on `@assert`, so configuration checks are not elided by optimized Julia runs.
- The docs include a shipped CPU quick-start scene and a small docs smoke runner for the examples most likely to drift.
- The docs build runs with `checkdocs = :exports` and no `warnonly` suppressions.
- Canopy surface integration now targets CanopyOptics.jl 0.2, including the
  analytic Stokes-aware canopy Z-matrix API used by `CanopySurface`.
- Replaced `eval(Meta.parse(...))` in YAML spec_bands parsing with a safe literal parser.
- Fixed `Vector{Any}` fields in Raman types to use proper parametric types.
- Cleaned up dead code: removed unused prototype files, plot scripts, and old test utilities.
- Cleaned up dated and stale TODO comments across inelastic scattering and surface modules.

### Known Limitations

- Linearized RT (Jacobians) does not yet support Raman scattering (`RS_type = noRS()` only).
- SIF helper exports exist, but the end-to-end SIF workflow remains experimental until fixture/data policy is settled.
- EMIT/MODTRAN comparison scripts are parsed in CI but not executed by default because they require external LUTs and machine-local scene inputs.
- Aerosol scene input is documented, but the higher-level aerosol API is still being stabilized.
- Wigner 3j symbol validation test remains disabled by default (~60s runtime).
