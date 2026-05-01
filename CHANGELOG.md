# Changelog

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
