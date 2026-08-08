# CLAUDE.md - vSmartMOM.jl

> **First read for any agent: [AGENTS.md](AGENTS.md)**
> AGENTS.md gives a 3-minute orientation — the narrative spine (problem → MOM →
> layer optics → solver → architecture), the two non-obvious tricks that
> compound (exact finite-δ elemental + scattering-only `N_doubl` sizing), the
> code-anchor table, and pointers into the Concepts arc and the verified
> equation crib sheet at [docs/dev_notes/theory_references.md](docs/dev_notes/theory_references.md).
> Read it before doing anything substantive.

> **⚠ CONVENTIONS — read before any cross-validation or data-import work:**
> [docs/src/pages/conventions.md](docs/src/pages/conventions.md) documents
> vSmartMOM's Hovenier-style γ sign, Stokes Q/U/V signs, and azimuth (Δφ)
> definition, plus exactly which sign-flips are needed when comparing to /
> importing from VLIDORT-format data (Siewert 2000 PROBLEM_IIA, PROBLEMIII,
> solar_tester). Most "vSmartMOM doesn't match VLIDORT" bugs are convention
> mismatches caught in this page.

## Project Overview

vSmartMOM.jl is a Julia package for vectorized atmospheric radiative transfer using the Matrix Operator Method (MOM). It computes polarized radiances and analytic Jacobians across spectral bands, supporting elastic/inelastic (Raman) scattering, multiple surface BRDF types, and CPU/GPU execution.

## Build & Test

```bash
# Run full test suite
julia --project=test test/runtests.jl

# Run a single test file
julia --project=test -e 'using vSmartMOM; include("test/test_Scattering.jl")'

# Run via Pkg
julia --project=test -e 'using Pkg; Pkg.test()'
```

No linter or formatter configured. Julia compatibility: 1.9-1.12.

## Architecture

### Module Load Order (`src/vSmartMOM.jl`)

1. **Architectures** (`src/Architectures.jl`) — CPU/GPU abstraction via KernelAbstractions
2. **Absorption** (`src/Absorption/Absorption.jl`) — Line-by-line cross-sections (HITRAN, Voigt/Doppler/Lorentz). Note: since the `externalAbsorption` merge, the RT production pipeline (`model_from_parameters`) computes LBL absorption via the external [AtmosphericAbsorption.jl](https://github.com/RemoteSensingTools/AtmosphericAbsorption.jl) package (TIPS-2021). The internal `Absorption` module remains the standalone API for direct σ(ν, T, p) queries plus CIA/MT-CKD/LUT paths.
3. **Artifacts** (`src/Artifacts/`) — HITRAN data management (edition preferences, hitran.org API client, artifact dispatch)
4. **Scattering** (`src/Scattering/Scattering.jl`) — Mie scattering (NAI2 and PCW decomposition)
5. **InelasticScattering** (`src/Inelastic/InelasticScattering.jl`) — Raman scattering (RRS, VS)
6. **CoreRT** (`src/CoreRT/CoreRT.jl`) — Core RT solver (adding-doubling, surfaces, layer optics)
7. **SolarModel** (`src/SolarModel/SolarModel.jl`) — Solar irradiance spectrum
8. **IO** (`src/IO/IO.jl`) — YAML config parsing, NetCDF readers, GEOSChem integration

GPU is a weak dependency via `ext/vSmartMOMCUDAExt.jl` (loads when CUDA.jl is present).

### Main Pipeline

```
YAML/TOML/Dict
    → read_parameters()        → vSmartMOM_Parameters   (unified entry point)
    → model_from_parameters()  → RTModel
    → rt_run(model)            → ObserverRTResult (named endpoint/level radiances)
    → rt_run_toa(model)        → TOA upwelling only (opt-in external-solar SFI)
```

`parameters_from_yaml(path)` is the YAML-specific alias and still works; use `read_parameters` for TOML or `Dict` inputs.

With the default `obs_alt: [0]`, `ObserverRTResult` iterates as the historical
forward tuple, so `R, T = rt_run(model)` still binds TOA upwelling and BOA
downwelling. Interior-height radiances are available through
`result.levels`.

`rt_run_toa` requires `model.quad_points.external_solar == true`. In this
opt-in Gauss/SFI representation, scalar `μ₀` is excluded from the diffuse
operator and retained on `phase_qp_μ` as the exact direct-beam source column.
The path currently supports forward elastic `noRS` with Lambertian surfaces;
it does not allocate or postprocess BOA, HDR, or BHR, and it does not support
linearized, Raman/VRS, `rt_run_ss`, non-Lambertian, or interior-sensor runs.
The embedded-`μ₀` representation remains the default; unsupported paths
reject external-solar models rather than falling back silently.

Linearized variant: `model_from_parameters(LinMode(), params)` then
`rt_run(model, lin_model, NAer, NGas, NSurf)` returns an
`ObserverRTResultLin`. It remains iterable as `(R, T, dR, dT)` and exposes
strict-interior radiances/Jacobians through `result.levels`.

### RTModel Hierarchy (Oceananigans-style)

`model_from_parameters()` returns an `RTModel{ARCH, FT}` with physics-based sub-structs:

```
RTModel{ARCH, FT} <: AbstractRTModel{ARCH, FT}
├── architecture :: ARCH                    # CPU() or GPU()
├── solver       :: SolverConfig{FT}        # polarization, quadrature, truncation, m_max_bands
├── geometry     :: ObsGeometry{FT}         # angles + resolved observer interfaces
├── quad_points  :: QuadPoints{FT}          # diffuse qp_μ/wt_μ + phase_qp_μ; optional external μ₀
├── atmosphere   :: Atmosphere{FT}          # profile + spec_bands
├── optics       :: Optics{FT}             # ALL optical properties
│   ├── rayleigh :: RayleighScattering{FT}  # greek_rayleigh, greek_cabannes, ϖ_Cabannes
│   ├── aerosols :: AerosolState{FT}        # aerosol_optics, τ_aer
│   ├── τ_abs    :: Vector{Matrix{FT}}      # absorption optical depth
│   └── τ_rayl   :: Vector{Matrix{FT}}      # Rayleigh optical depth
└── surfaces     :: Vector{AbstractSurfaceType}  # per-band BRDF
```

**Accessor functions** (work on RTModel):
- `architecture(model)`, `array_type(model)` — from Architectures
- `CoreRT.polarization_type(model)`, `CoreRT.float_type(model)`, `CoreRT.n_aerosols(model)`
- `get_surface(model, iBand)`, `get_surfaces(model)`, `get_spec_bands(model)`

**Convenience forwarding**: `model.τ_abs`, `model.profile`, `model.obs_geom`, etc. work via `Base.getproperty` override.

**AD boundary**: Differentiable state lives in `optics` (τ_abs, τ_aer, aerosol_optics) and `surfaces`. Fixed config lives in `solver`, `geometry`, `quad_points`.

### CoreRT Solver Flow (Adding-Doubling)

For each Fourier moment m = 0..m_max_bands[iBand] (v0.7 — order-semantics; trait-derived per-component bound):
1. **Elemental** — single-scattering layer → AddedLayer (r, t, j)
2. **Doubling** — double thin layers ndoubl times to full optical depth
3. **Interaction** — combine layers top-to-bottom: CompositeLayer (R, T, J) + AddedLayer
4. **Surface** — create surface layer from BRDF, interact with composite atmosphere
5. **Postprocessing** — azimuthal weighting, VZA interpolation

### Key Types

| Type | Location | Purpose |
|------|----------|---------|
| `vSmartMOM_Parameters` | `src/CoreRT/types.jl` | User config (from YAML) |
| `RTModel` | `src/CoreRT/types.jl` | Hierarchical RT model |
| `SolverConfig` | `src/CoreRT/types.jl` | RT solver settings (polarization, quadrature, truncation) |
| `Atmosphere` | `src/CoreRT/types.jl` | AtmosphericProfile + spec_bands |
| `Optics` | `src/CoreRT/types.jl` | All optical properties (rayleigh, aerosol, abs, rayl) |
| `AtmosphericProfile` | `src/CoreRT/types.jl` | T, p, q, VMR profiles |
| `ObsGeometry` | `src/CoreRT/types.jl` | SZA, VZA, VAZ, requested/resolved observer heights |
| `ObserverRTResult` | `src/CoreRT/types.jl` | Named TOA, BOA, and interior-height forward outputs |
| `LevelRadiance` | `src/CoreRT/types.jl` | Co-located up/down radiances at one interior interface |
| `ObserverRTResultLin` | `src/CoreRT/types_lin.jl` | Named endpoint and interior-height linearized outputs |
| `LevelRadianceLin` | `src/CoreRT/types_lin.jl` | Co-located radiances and analytic Jacobians at one interior interface |
| `QuadPoints` | `src/CoreRT/types.jl` | Diffuse quadrature, phase-evaluation grid, scalar μ₀, and external-solar flag |
| `CompositeLayer` | `src/CoreRT/types.jl` | Accumulated R, T, J matrices (uppercase) |
| `AddedLayer` | `src/CoreRT/types.jl` | Single-layer r, t, j matrices (lowercase) |
| `Aerosol` | `src/Scattering/types.jl` | Size distribution + refractive index (nr, ni) |
| `MieModel` | `src/Scattering/types.jl` | Mie computation config |
| `GreekCoefs` | `src/Scattering/types.jl` | Phase matrix expansion coefficients |
| `AerosolOptics` | `src/Scattering/types.jl` | Computed SSA, extinction, greek coefs |
| `HitranTable` | `src/Absorption/types.jl` | HITRAN line parameters |
| `CoreScatteringOpticalProperties` | `src/CoreRT/types.jl` | Per-layer tau, omega, Z+/Z- (supports + and *) |
| `ParameterLayout` | `src/CoreRT/parameter_layout.jl` | Jacobian index arithmetic |

### Surfaces (`src/CoreRT/Surfaces/`)

| Type | File | Key Parameters |
|------|------|----------------|
| `LambertianSurfaceScalar` | `lambertian_surface.jl` | albedo |
| `LambertianSurfaceSpectrum` | `lambertian_surface.jl` | albedo vector |
| `rpvSurfaceScalar` | `rpv_surface.jl` | rho0, rho_c, k, Theta |
| `RossLiSurfaceScalar` | `rossli_surface.jl` | fvol, fgeo, fiso |
| `CoxMunkSurface` | `coxmunk_surface.jl` | wind_speed, n_water, whitecaps, shadowing |
| `CanopySurface` | `canopy_surface.jl` | soil, LAI, n_layers, leaf R/T |

All surfaces implement `create_surface_layer!()`. Linearized variants have `_lin` suffix files.

### Polarization Types (`src/Scattering/types.jl`)

- `Stokes_I` — scalar intensity only (n=1)
- `Stokes_IQU` — linear polarization (n=3)
- `Stokes_IQUV` — full Stokes polarization (n=4)

## Conventions

- **FT type parameter**: Most types/functions parameterized by float type (Float32/Float64)
- **`_lin` suffix**: Linearized (Jacobian) variant of functions and types
- **Unicode variables**: `τ` (optical depth), `ϖ` (SSA), `μ` (cosine zenith) used directly
- **Sign convention**: `+` = incoming/downward, `-` = outgoing/upward
- **Layer naming**: CompositeLayer uses uppercase (R, T, J), AddedLayer uses lowercase (r, t, j)
- **3D matrices**: diffuse RT matrices are `(NquadN, NquadN, nSpec)` where
  `NquadN = Nquad * n_stokes`. External-solar phase matrices may carry one
  additional exact μ₀ row/column, while the diffuse operators remain square
  on `qp_μ`.
- **External solar is not a stream**: with `external_solar=true`, scalar μ₀ and
  `iμ₀_phase` provide `Zₘ(μᵢ,μ₀)` source coupling; legacy operator indices
  `iμ₀`/`iμ₀Nstart` are zero sentinels. Five weighted streams plus one
  distinct VZA in IQU gives an 18×18 diffuse operator.
- **Spectral units**: wavenumber in cm⁻¹ internally; wavelength in micrometers for Mie
- **Vertical units**: profile pressure is hPa; `obs_alt` is geometric km above BOA
- **Profile direction**: TOA to BOA (top of atmosphere to bottom)

## File Structure

```
src/
  vSmartMOM.jl                # Entry point, module definition
  Architectures.jl            # CPU/GPU abstraction
  Artifacts/                  # HITRAN data management
    hitran_preferences.jl     # Edition selection (set/get_hitran_edition!, available_hitran_editions)
    hitran_api.jl             # Download client (fetch_hitran, fetch_hitran_by_ids)
    artifact_helper.jl        # Unified artifact() dispatch (legacy Artifacts vs scratch cache)
    download_hitran.jl        # Offline artifact generation script (developer-only)
  Absorption/
    Absorption.jl             # Module entry
    types.jl                  # HitranTable, broadening types, CEF types
    read_hitran.jl            # HITRAN file parser
    compute_absorption_cross_section.jl
    complex_error_functions.jl
    constants/                # Molecular weights, TIPS_2017, physical constants
  Scattering/
    Scattering.jl             # Module entry
    types.jl                  # Aerosol, MieModel, GreekCoefs, polarization types
    types_lin.jl              # Linearized variants
    compute_NAI2.jl           # NAI2 Fourier decomposition
    compute_PCW.jl            # Precomputed Wigner decomposition
    compute_Z_matrices.jl     # Z scattering matrices
    mie_helper_functions.jl   # Mie a,b coefficients, Q factors
  Inelastic/
    InelasticScattering.jl    # Module entry
    types.jl                  # noRS, RRS, VS_0to1, VS_1to0
  CoreRT/
    CoreRT.jl                 # Module entry
    types.jl                  # Parameters, Model, Layers, Surfaces, QuadPoints
    types_lin.jl              # Linearized layer types
    parameter_layout.jl       # Jacobian index layout
    rt_run.jl                 # Forward RT entry
    rt_run_lin.jl             # Linearized RT entry
    DefaultParameters.yaml    # Default config template
    CoreKernel/               # elemental, doubling, interaction solvers (+_lin, +_inelastic)
    Surfaces/                 # Surface BRDF implementations
    LayerOpticalProperties/   # Effective layer optics, delta-M truncation
    tools/                    # model_from_parameters, atmo_prof, postprocessing, helpers
  SolarModel/                 # Solar irradiance
  IO/
    IO.jl                     # Module entry
    Parameters.jl             # YAML -> vSmartMOM_Parameters (BRDF_MAP, QUAD_MAP, etc.)
    AtmosProfile.jl           # Atmospheric profile readers
    Sources.jl                # GeosChemSource, NetCDFGridSource
    NetCDF/                   # GEOSChem NetCDF reader
ext/
  vSmartMOMCUDAExt.jl         # CUDA weak dependency extension
config/                       # Example YAML configurations
test/
  runtests.jl                 # Test orchestrator (11 test sets)
  test_helpers.jl             # run_lin_rt, rel_errors, fd_jacobian_R
  test_parameters/            # Test YAML configs
```

## Common Workflows

### Adding a New Surface Model

1. Define struct subtyping `AbstractSurfaceType` in `src/CoreRT/types.jl`
2. Create `src/CoreRT/Surfaces/mysurf_surface.jl` implementing `create_surface_layer!()`
3. Include the file in `src/CoreRT/CoreRT.jl`
4. Add type name to `BRDF_MAP` in `src/IO/Parameters.jl`
5. Export the type from `CoreRT.jl`
6. Add YAML config example in `config/`
7. Add test in `test/`

### Adding a New Test

1. Create `test/test_myfeature.jl` with `@testset` blocks
2. Add `@testset` entry in `test/runtests.jl`
3. Test YAML files go in `test/test_parameters/`
4. Use helpers from `test/test_helpers.jl` (e.g., `run_lin_rt`, `rel_errors`, `fd_jacobian_R`)

### Running Forward RT

```julia
using vSmartMOM
params = read_parameters("config/lambertian_land.yaml")   # preferred; also accepts TOML/Dict
model  = model_from_parameters(params)
R, T   = rt_run(model)
```

### Running Linearized RT (Jacobians)

```julia
using vSmartMOM
params = read_parameters("config/ocean_coxmunk.yaml")
model, lin_model = model_from_parameters(LinMode(), params)
NAer = length(params.scattering_params.rt_aerosols)
NGas = size(lin_model.τ̇_abs[1], 1)
NSurf = 1
R, T, dR, dT = rt_run(model, lin_model, NAer, NGas, NSurf)
```

### Writing or modifying YAML/TOML configs

**Schema reference**: [`docs/src/pages/IO/Schema/`](docs/src/pages/IO/Schema/)
has one page per top-level block (`radiative_transfer.md`,
`geometry.md`, `surface.md`, `aerosols.md`, etc.) with full field
reference and examples.

**v0.7 schema (Phase D)**: prefer `nstreams` over the legacy
`max_m` / `l_trunc`. Public contract: `stream_l_cap = 2·nstreams - 1`
regardless of `quadrature_type`. `nstreams ≥ 3` (Rayleigh m=2
minimum). `truncation: auto` is the recommended default — resolves
to `NoTruncation()` if the phase function fits, `δBGE(N, Δ_angle)`
otherwise. Legacy YAMLs with `max_m`/`l_trunc` keep working unchanged.

**Editor support**: TOML configs get autocomplete + validation via
the repo's `.taplo.toml` (wires `schemas/vsmartmom-parameters.schema.json`).
For YAML, add `# yaml-language-server: $schema=...` at the top —
recipe in [`Schema.md`](docs/src/pages/IO/Schema.md#editor-support--autocomplete--inline-validation).

### Modifying YAML Config Parsing

All string-to-type maps are in `src/IO/Parameters.jl`:
- `BRDF_MAP` — surface types
- `QUAD_MAP` — quadrature types
- `POLARIZATION_MAP` — polarization types
- `BROADENING_MAP` — line broadening functions
- `DECOMP_MAP` — Fourier decomposition (NAI2/PCW)
- `ARCH_MAP` — CPU/GPU architecture
- `FLOAT_MAP` — Float32/Float64
