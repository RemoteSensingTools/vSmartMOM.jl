# Changelog

## v2.1.0 — Fourier-stream resolution + source-term refactor

Additive minor release over v2.0.0. Highlights: first-class `AbstractSource`
type vocabulary (`SolarBeam`, `BlackbodySource`, `SurfaceSIF`); user-facing
`nstreams` / `truncation: auto` Fourier-resolution schema; per-model BLAS
thread cap; VLIDORT 2.8.3 baseline validation suite (Siewert 2000 Problem IIA,
solar_tester scalar/vector). All v2.0.0 YAMLs run byte-equal — `max_m` /
`l_trunc` keep working through the legacy parser branch.

### New user-facing schema (Phase D)

- **`nstreams` is the primary resolution knob.** Public contract:
  `stream_l_cap = 2·nstreams - 1` regardless of `quadrature_type`.
  Minimum 3 (Rayleigh contributes through m=2). Default 13 when
  legacy `max_m`/`l_trunc` are absent.
- **`truncation: auto`** — VLIDORT-`DO_DELTAM_SCALING`-style mode.
  Resolves at model build time: `NoTruncation()` when phase fits
  within `stream_l_cap`, `δBGE(N, Δ_angle)` otherwise. Logs the
  choice via `@info`.
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
