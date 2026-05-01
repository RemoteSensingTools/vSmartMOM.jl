# Release Notes and Migration

**For:** users moving to the registration-ready vSmartMOM 2.0 interface.

**Next:** [Quick Start](quickstart.md), [Configure a Scene](IO/Overview.md), [Compute Jacobians](jacobians.md).

This page summarizes the user-visible changes in the 2.0 line. It is written as
a migration guide, not as a complete git history.

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
- `geometry`: solar/viewing geometry and observer altitude.
- `quad_points`: quadrature abscissae and weights.
- `atmosphere`: atmospheric profile and spectral bands.
- `optics`: Rayleigh, aerosol, gas, and linearized optical-property state.
- `surfaces`: one surface model per spectral band.

The old flat `vSmartMOM_Model` container has been replaced. Compatibility
getters still expose many legacy fields, but new code should prefer the
hierarchical fields and accessor functions documented in the
[Library](api_reference.md).

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
R, T, R_jac, T_jac = rt_run_lin(model, lin_model)
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

## Known WIP Areas

- Aerosol scene input is documented, but the high-level API is still being
  stabilized.
- Internal export cleanup is intentionally deferred until after v2.0.0 so this
  release does not tighten public API surface area at the same time as the
  registration cutover.
- Unified offline source-function integration and thermal emission are design
  topics, not implemented user-facing features in this release line.

## Future Developments

- Mixed automatic differentiation and hand-coded linearization.
- Coupled ocean-atmosphere radiative transfer.
- More line-shape and line-mixing functionality, including expansion to more
  spectral databases.
- Raman performance and memory-pressure improvements.
- Correlated-k and other spectral dimension-reduction methods.
