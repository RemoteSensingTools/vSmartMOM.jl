# CLAUDE.md - vSmartMOM.jl

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
2. **Artifacts** (`src/Artifacts/`) — HITRAN data download helpers
3. **Absorption** (`src/Absorption/Absorption.jl`) — Line-by-line cross-sections (HITRAN, Voigt/Doppler/Lorentz)
4. **Scattering** (`src/Scattering/Scattering.jl`) — Mie scattering (NAI2 and PCW decomposition)
5. **InelasticScattering** (`src/Inelastic/InelasticScattering.jl`) — Raman scattering (RRS, VS)
6. **CoreRT** (`src/CoreRT/CoreRT.jl`) — Core RT solver (adding-doubling, surfaces, layer optics)
7. **SolarModel** (`src/SolarModel/SolarModel.jl`) — Solar irradiance spectrum
8. **IO** (`src/IO/IO.jl`) — YAML config parsing, NetCDF readers, GEOSChem integration

GPU is a weak dependency via `ext/vSmartMOMCUDAExt.jl` (loads when CUDA.jl is present).

### Main Pipeline

```
YAML config
    → parameters_from_yaml()   → vSmartMOM_Parameters
    → model_from_parameters()  → RTModel
    → rt_run(model)            → (R, T) reflectance/transmittance
```

Linearized variant: `model_from_parameters(LinMode(), params)` then `rt_run(model, lin_model, NAer, NGas, NSurf)` returns `(R, T, dR, dT)`.

### RTModel Hierarchy (Oceananigans-style)

`model_from_parameters()` returns an `RTModel{ARCH, FT}` with physics-based sub-structs:

```
RTModel{ARCH, FT} <: AbstractRTModel{ARCH, FT}
├── architecture :: ARCH                    # CPU() or GPU()
├── solver       :: SolverConfig{FT}        # polarization, quadrature, truncation, max_m
├── geometry     :: ObsGeometry{FT}         # sza, vza, vaz, obs_alt
├── quad_points  :: QuadPoints{FT}          # μ₀, qp_μ, wt_μ, Nquad
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

For each Fourier moment m = 0..max_m:
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
| `ObsGeometry` | `src/CoreRT/types.jl` | SZA, VZA, VAZ, observer altitude |
| `QuadPoints` | `src/CoreRT/types.jl` | Quadrature points, weights, mu0 |
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
- **3D matrices**: RT matrices are `(NquadN, NquadN, nSpec)` where `NquadN = Nquad * n_stokes`
- **Spectral units**: wavenumber in cm⁻¹ internally; wavelength in micrometers for Mie
- **Pressure**: hPa in YAML configs; Pa for obs_alt
- **Profile direction**: TOA to BOA (top of atmosphere to bottom)

## File Structure

```
src/
  vSmartMOM.jl                # Entry point, module definition
  Architectures.jl            # CPU/GPU abstraction
  Artifacts/                  # HITRAN data artifact helpers
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
params = parameters_from_yaml("config/lambertian_land.yaml")
model  = model_from_parameters(params)
R, T   = rt_run(model)
```

### Running Linearized RT (Jacobians)

```julia
using vSmartMOM
params = parameters_from_yaml("config/ocean_coxmunk.yaml")
model, lin_model = model_from_parameters(LinMode(), params)
NAer = length(params.scattering_params.rt_aerosols)
NGas = size(lin_model.tau_dot_abs[1], 1)
NSurf = 1
R, T, dR, dT = rt_run(model, lin_model, NAer, NGas, NSurf)
```

### Modifying YAML Config Parsing

All string-to-type maps are in `src/IO/Parameters.jl`:
- `BRDF_MAP` — surface types
- `QUAD_MAP` — quadrature types
- `POLARIZATION_MAP` — polarization types
- `BROADENING_MAP` — line broadening functions
- `DECOMP_MAP` — Fourier decomposition (NAI2/PCW)
- `ARCH_MAP` — CPU/GPU architecture
- `FLOAT_MAP` — Float32/Float64
