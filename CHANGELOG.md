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

- **Cox-Munk ocean surface**: Full polarized BRDF with Fresnel reflection, wind-speed-dependent roughness, whitecap correction, and Smith (1967) shadowing. Linearized variant included.

- **Accessor functions**: `architecture(model)`, `CoreRT.polarization_type(model)`, `CoreRT.float_type(model)`, `CoreRT.n_aerosols(model)`, `get_surface(model, iBand)`, `get_spec_bands(model)`.

### Improvements

- Replaced `eval(Meta.parse(...))` in YAML spec_bands parsing with a safe literal parser.
- Fixed `Vector{Any}` fields in Raman types to use proper parametric types.
- Cleaned up dead code: removed unused prototype files, plot scripts, and old test utilities.
- Cleaned up dated and stale TODO comments across inelastic scattering and surface modules.

### Known Limitations

- Linearized RT (Jacobians) does not yet support Raman scattering (`RS_type = noRS()` only).
- Wigner 3j symbol validation test remains disabled by default (~60s runtime).
