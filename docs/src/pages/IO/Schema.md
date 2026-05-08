# Configuration schema (IO)

vSmartMOM accepts configuration as **YAML**, **TOML**, or in-memory
`Dict`. All three are normalized through the same
`parameters_from_dict` pipeline (see [`Parameters.jl`](https://github.com/cfranken/vSmartMOM.jl/blob/main/src/IO/Parameters.jl)).
This page is the top-level index for the configuration schema; each
top-level block has its own detail page.

## Top-level blocks

| Block                                         | Required | Purpose |
|-----------------------------------------------|:--------:|---------|
| [`radiative_transfer`](Schema/radiative_transfer.md) | ✓ | RT solver, quadrature, polarization, Fourier resolution (v0.7 `nstreams` knob), truncation |
| [`geometry`](Schema/geometry.md)              | ✓        | SZA, VZA(s), VAZ(s), observer altitude |
| [`atmospheric_profile`](Schema/atmospheric_profile.md) | ✓ | T, p, q, profile reduction |
| [`absorption`](Schema/absorption.md)          | ◯        | Gas absorption: VMRs, broadening, HITRAN edition pointer |
| [`scattering`](Schema/aerosols.md)            | ◯        | Mie aerosol size distribution, phase-function, decomposition |
| `canopy` *(forthcoming)*                      | ◯        | Vegetation canopy + soil composite surface — page not yet authored |
| `sources` (programmatic only)                 | ◯        | See [`Schema/sources.md`](Schema/sources.md) — set via `model_from_parameters(...; sources=...)` for now |

The `surface:` field inside `radiative_transfer` has its own
detailed page at [`Schema/surface.md`](Schema/surface.md) covering
the BRDF vocabulary (`LambertianSurfaceScalar`, `CoxMunkSurface`,
`rpvSurfaceScalar`, etc.) and the Phase C trait dispatch table.

## v0.7 schema migration (Phase D)

vSmartMOM v0.7 promotes **`nstreams`** to the primary user-facing
resolution knob:

```yaml
radiative_transfer:
  nstreams: 13          # weighted streams per hemisphere
                        # public contract: stream_l_cap = 2·nstreams - 1
  truncation: auto      # NoTruncation if phase fits, δBGE otherwise
```

Legacy `max_m` / `l_trunc` configs **continue to work**; the parser
detects the legacy schema by the presence of either field and applies
the historical aggregator + `δBGE` default. New configs should prefer
`nstreams`. See [`radiative_transfer.md`](Schema/radiative_transfer.md)
for the precedence rules and [`docs/dev_notes/fourier_stream_resolution_plan.md`](../../dev_notes/fourier_stream_resolution_plan.md)
for the design rationale.

## Editor support — autocomplete + inline validation

vSmartMOM ships a JSON Schema at
[`schemas/vsmartmom-parameters.schema.json`](https://github.com/cfranken/vSmartMOM.jl/blob/main/schemas/vsmartmom-parameters.schema.json),
plus a `.taplo.toml` at the repo root that wires it to all
configuration TOML files in `config/`, `test/test_parameters/`,
`test/benchmarks/`, `test/vlidort_baseline/configs/`, `examples/`,
and `sandbox/`.

### TOML — Taplo (VS Code or CLI)

Install the [Even Better TOML extension](https://marketplace.visualstudio.com/items?itemName=tamasfe.even-better-toml)
or the standalone `taplo` CLI:

```bash
cargo install taplo-cli --locked
taplo lint config/lambertian_land.toml      # validates against the schema
taplo format config/*.toml                   # formats per .taplo.toml
```

### YAML — yaml-language-server

VS Code's [YAML extension](https://marketplace.visualstudio.com/items?itemName=redhat.vscode-yaml)
or any editor running `yaml-language-server` honors a `$schema`
directive at the top of a YAML file:

```yaml
# yaml-language-server: $schema=https://raw.githubusercontent.com/cfranken/vSmartMOM.jl/main/schemas/vsmartmom-parameters.schema.json

radiative_transfer:
  nstreams: 13
  ...
```

For local development without internet, point `$schema` at the local
file:

```yaml
# yaml-language-server: $schema=../../../schemas/vsmartmom-parameters.schema.json
```

Once wired, hovering over any field shows the inline description, and
typos / out-of-range values surface as squiggles before runtime.

## Minimal new-schema YAML example

```yaml
radiative_transfer:
  spec_bands:
    - "[18867.92 18868.92]"
  surface:
    - LambertianSurfaceScalar(0.0)
  polarization_type: Stokes_IQU()
  Δ_angle: 2.0
  depol: 0.0
  float_type: Float64
  architecture: Architectures.CPU()
  # nstreams omitted ⇒ default 13
  # truncation omitted ⇒ "auto"
  # quadrature_type omitted ⇒ GaussLegQuad()

geometry:
  sza: 30
  vza: [10, 20, 40]
  vaz: [0, 0, 0]
  obs_alt: 1000.0

atmospheric_profile:
  T: [260, 262, 264]
  p: [0.1, 100.0, 800.0, 1005.0]
  profile_reduction: -1
```

## Reading from external sources

```julia
using vSmartMOM

# YAML / TOML / Dict — all routed through parameters_from_dict
params = parameters_from_yaml("config/quickstart.yaml")
params = parameters_from_file("config/quickstart.toml")
params = parameters_from_dict(my_dict)
```

## See also

- [`docs/src/pages/IO/Examples.md`](Examples.md) — end-to-end usage
  examples
- [`docs/src/pages/conventions.md`](../conventions.md) — sign conventions
  for VLIDORT cross-validation
- [`docs/dev_notes/fourier_stream_resolution_plan.md`](../../dev_notes/fourier_stream_resolution_plan.md)
  — v0.7 Fourier/stream resolution design memo
