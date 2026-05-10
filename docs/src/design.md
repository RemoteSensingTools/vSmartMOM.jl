# Architecture and Design

This document describes the internal architecture of vSmartMOM.jl, the data flow from parameters to radiance, and the key type hierarchy.

## Module Dependency Graph

```
vSmartMOM (root module)
│
├── Architectures          CPU/GPU abstraction (KernelAbstractions)
│
├── Absorption             HITRAN line-by-line + ABSCO lookup tables
│   └── uses: Architectures
│
├── Scattering             Mie theory (NAI2 / PCW), phase functions, Greek coefficients
│   └── uses: Architectures
│
├── InelasticScattering    Rotational Raman (RRS), Vibrational (VS) scattering
│   └── uses: Scattering
│
├── Aerosols               Size distributions, refractive indices, TOMAS-15, two-moment
│
├── CoreRT                 Radiative transfer solver (adding-doubling method)
│   └── uses: Absorption, Scattering, InelasticScattering, Architectures
│
├── SolarModel             Solar/stellar irradiance
│   └── uses: CoreRT
│
└── IO                     YAML/TOML/Dict parameter files, NetCDF, atmospheric profiles
    └── uses: CoreRT (for types)
```

## Data Flow: Parameters to Radiance

The complete forward model pipeline is:

### 1. Parameter Loading

```
YAML/TOML file, Dict, or IOSource  ──>  read_parameters()  ──>  vSmartMOM_Parameters
```

`read_parameters` is the user-facing loader. It dispatches on file paths, in-memory dictionaries, and typed IO sources such as GEOS-Chem NetCDF. The explicit aliases `parameters_from_file`, `parameters_from_dict`, and `parameters_from_source` call the same parser after loading the configuration.

`vSmartMOM_Parameters` holds all configuration: spectral bands, geometry, atmospheric profile, gas absorption settings, aerosol properties, surface type, and architecture (CPU/GPU).

### 2. Model Construction

```
vSmartMOM_Parameters  ──>  model_from_parameters()  ──>  RTModel
```

This step:
- Computes atmospheric profile (pressure, temperature, layer boundaries)
- Computes Rayleigh scattering optical depths per layer
- Runs Mie calculations for each aerosol type (greek coefficients, single-scattering albedo, extinction)
- Computes gas absorption cross-sections (HITRAN or ABSCO lookup)
- Sets up quadrature points for the discrete ordinate method

For the linearized model, `model_from_parameters(LinMode(), params)` additionally returns an `RTModelLin` struct with analytic derivatives of optical properties.

### 3. Optical Property Assembly

```
RTModel  ──>  constructCoreOpticalProperties()  ──>  CoreScatteringOpticalProperties[]
```

Per Fourier moment `m`, per layer `iz`, this combines Rayleigh + aerosol + gas absorption into:
- `τ` — total optical depth
- `ϖ` — single-scattering albedo  
- `Z⁺⁺`, `Z⁻⁺` — scattering phase matrix (truncated Legendre expansion)

These are the **3 core optical parameters** that the RT solver differentiates with respect to.

### 4. RT Solver (Adding-Doubling)

```
CoreScatteringOpticalProperties  ──>  rt_kernel!()  ──>  R, T, J (per layer)
```

The RT kernel processes each layer in three steps:

**Elemental** (`elemental!`): Computes the thin-layer single-scattering solution for reflection `r`, transmission `t`, and source `j` matrices.

**Doubling** (`doubling!`): Doubles the elemental layer `ndoubl` times to get the full homogeneous layer. Uses the geometric series `(I - RR)⁻¹` and batched matrix operations.

**Interaction** (`interaction!`): Adds the current layer to the composite (accumulated) layer using the adding method. Four cases depending on the scattering interface type:
- `ScatteringInterface_00`: No scattering in either layer
- `ScatteringInterface_01`: Scattering only in new layer
- `ScatteringInterface_10`: Scattering only in existing composite
- `ScatteringInterface_11`: Scattering in both

### 5. Surface and Postprocessing

```
CompositeLayer  ──>  create_surface_layer!()  ──>  interaction!()  ──>  postprocessing_vza!()
```

After processing all atmospheric layers, the surface BRDF is applied (Lambertian, RPV, or Ross-Li), a final interaction step couples the atmosphere to the surface, and the result is interpolated to the requested viewing zenith angles.

## Three Parallel RT Paths

The codebase supports three RT modes, each with dedicated kernel files:

| Mode | Files | Description |
|------|-------|-------------|
| **Forward** | `elemental.jl`, `doubling.jl`, `interaction.jl` | Computes R, T radiance. Supports elastic + inelastic (Raman). |
| **Linearized** | `elemental_lin.jl`, `doubling_lin.jl`, `interaction_lin.jl` | Computes R, T plus analytic Jacobians ∂R/∂(τ,ϖ,Z) per layer. Elastic only. |
| **Inelastic** | `elemental_inelastic.jl`, `doubling_inelastic.jl`, `interaction_inelastic.jl` | Handles cross-wavelength coupling for Raman scattering. 4D matrices. |

## Key Type Hierarchy

### Layer Types

```
AbstractLayer
├── AddedLayer{FT}           Forward RT: r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻, j₀⁺, j₀⁻  (3D arrays)
├── AddedLayerRS{FT}         Inelastic RT: same + ie* fields              (3D + 4D arrays)
├── CompositeLayer{FT}       Accumulated atmosphere: R⁻⁺, T⁺⁺, ...       (3D arrays)
└── CompositeLayerRS{FT}     Inelastic composite: same + ie* fields       (3D + 4D arrays)

AbstractLayerLin
├── AddedLayerLin{FT}        Linearized layer: ṙ⁻⁺, ṫ⁺⁺, ... + ap_ṙ⁻⁺, ap_ṫ⁺⁺, ...
└── CompositeLayerLin{FT}    Linearized composite: Ṙ⁻⁺, Ṫ⁺⁺, ...
```

- `AddedLayer` holds single-layer RT matrices (reflection, transmission, source).
- `CompositeLayer` holds the accumulated atmosphere from TOA down to the current level.
- The `_lin` variants carry derivatives with respect to the 3 core parameters (τ, ϖ, Z) and optionally all state parameters (`ap_*` fields).
- The `_RS` variants add inelastic (Raman) cross-wavelength coupling matrices as 4D arrays.

### Optical Property Types

```
AbstractOpticalProperties
├── CoreScatteringOpticalProperties{FT,FT2,FT3}    τ, ϖ, Z⁺⁺, Z⁻⁺
├── CoreAbsorptionOpticalProperties{FT}             τ_abs
└── CoreDirectionalScatteringOpticalProperties      τ, ϖ, Z⁺⁺, Z⁻⁺, G (directional)

UmbrellaCoreScatteringOpticalProperties{FWD,LIN}    Bundles forward + linearized
```

These support `+` (combining layers) and `*` (scaling) operations.

### Scattering Types

```
AbstractRamanType{FT}
├── noRS{FT}                No Raman scattering (elastic only)
├── RRS{FT}                 Rotational Raman Scattering
├── VS_0to1{FT}             Vibrational Stokes (0→1)
├── VS_1to0{FT}             Vibrational Anti-Stokes (1→0)
└── *_plus variants          Concatenated mode for multi-band
```

These dispatch the RT solver to the appropriate elemental/doubling/interaction implementations.

### Architecture Abstraction

```
AbstractArchitecture
├── CPU()      Arrays are Array, kernels run on CPU
└── GPU()      Arrays are CuArray, kernels run via KernelAbstractions + CUDA
```

Key functions:
- `devi(arch)` — returns the KernelAbstractions device
- `array_type(arch)` — returns `Array` or `CuArray`
- `architecture(arr)` — infers architecture from array type
- `synchronize_if_gpu()` — no-op on CPU, `CUDA.synchronize()` on GPU

## Linearized RT and Jacobian Workflow

### The AD Boundary

The codebase enforces a clean boundary between two worlds:

```
                     AD Zone                           Pure Float Zone
           ┌───────────────────────────┐    ┌──────────────────────────────────┐
           │  ForwardDiff.Dual allowed │    │  FT ∈ {Float32, Float64} only   │
           │                           │    │                                  │
           │  Mie cross-sections       │    │  elemental!, doubling!,          │
           │  Gas absorption           │    │  interaction!, rt_kernel!        │
           │  Atmospheric profiles     │    │                                  │
           │  Surface BRDF params      │    │  Analytic ∂R/∂(τ,ϖ,Z)           │
           └──────────┬────────────────┘    └──────────────┬───────────────────┘
                      │                                    │
                      ▼                                    ▼
           OpticalPropertyJacobian            lin_added_layer_all_params!
           ∂(τ,ϖ,Z⁺⁺,Z⁻⁺)/∂xⱼ      ──────>  Chain rule: ∂R/∂xⱼ
```

The `OpticalPropertyJacobian` (alias for `CoreScatteringOpticalPropertiesLin`) is the handoff struct. It contains `(τ̇, ϖ̇, Ż⁺⁺, Ż⁻⁺)` — the derivatives of the four core optical properties with respect to each retrieval parameter.

### Chain Rule

The chain rule `lin_added_layer_all_params!` maps per-layer optical property derivatives to radiance derivatives:

```
∂R/∂x_j = Σ_layers [ ∂R/∂τ · ∂τ/∂x_j + ∂R/∂ϖ · ∂ϖ/∂x_j + ∂R/∂Z · ∂Z/∂x_j ]
```

The RT kernels compute `∂R/∂τ`, `∂R/∂ϖ`, `∂R/∂Z` analytically (stored in `AddedLayerLin.ṫ⁺⁺[1:3,...]`). The optical property code provides `∂τ/∂x_j`, `∂ϖ/∂x_j`, `∂Z/∂x_j` either analytically or via ForwardDiff.

### Parameter Layout via `ParameterLayout`

Instead of hardcoding `7*NAer + NGas + NSurf`, the codebase uses a `ParameterLayout` struct
(defined in `src/CoreRT/parameter_layout.jl`):

```julia
layout = ParameterLayout(aerosol_params=7, n_aerosols=NAer,
                         n_gases=NGas, n_surface=NSurf)
n_total(layout)              # total Nparams
aerosol_range(layout, iaer)  # indices for aerosol iaer
gas_range(layout)            # indices for gas VMRs
surface_index(layout, iBand) # surface param for band iBand
```

| Index range               | Parameter                                       | Method   |
|---------------------------|------------------------------------------------ |----------|
| `aerosol_range(layout, i)` | `τ_ref, nᵣ, nᵢ, rₘ, σ_g, p₀, σ_p` per aerosol | Analytic |
| `gas_range(layout)`        | Gas VMR scaling factors                         | Analytic |
| `surface_range(layout)`    | Surface albedo / BRDF parameters                | Analytic |

### Adding New Parameters

To add derivatives with respect to a new parameter `θ`:

1. **Analytic path**: Compute `∂τ/∂θ`, `∂ϖ/∂θ`, `∂Z/∂θ` analytically in `constructCoreOpticalProperties` and add to the `OpticalPropertyJacobian` arrays. Extend `ParameterLayout` to include the new parameter block.

2. **ForwardDiff path**: Wrap the optical property generation in `ForwardDiff.jacobian`:
   ```julia
   function optical_props_wrt_theta(θ_vec)
       # ... compute τ, ϖ as functions of θ_vec ...
       return [τ; ϖ]
   end
   J = ForwardDiff.jacobian(optical_props_wrt_theta, θ₀)
   # Extract ∂τ/∂θ, ∂ϖ/∂θ from J, fill OpticalPropertyJacobian
   ```

3. **Hybrid path** (recommended for new parameters): Use `ForwardDiff.Dual` for the
   "outer" derivatives (parameter → raw optical properties) and retain the analytic
   RT kernel for the "inner" derivatives (optical properties → radiance). The boundary
   is the `OpticalPropertyJacobian` struct.

### Strategy by Parameter Type

| Parameter          | Recommended Method | Reason                                    |
|--------------------|--------------------|-------------------------------------------|
| Albedo             | Analytic           | Trivial derivative, already implemented   |
| RPV/RossLi BRDF    | ForwardDiff        | Low-dimensional, simple surface code      |
| τ_ref (aerosol OD) | Analytic           | Trivial: `∂τ/∂τ_ref = τ/τ_ref`           |
| p₀, σ_p (height)   | Analytic           | Already in `atmo_prof_lin.jl`             |
| nᵣ, nᵢ (refr. idx) | Analytic Mie       | Mie series is AD-hostile                  |
| rₘ, σ_g (size dist) | Analytic Mie      | Same reason                               |
| Gas VMR            | Analytic           | `∂τ_abs/∂VMR = cross_section`             |
| Surface pressure   | ForwardDiff        | Affects many paths (Rayleigh, absorption) |
| Temperature        | ForwardDiff        | Future; affects absorption cross-sections |

### Planned: 4 Core Variables with δ-M Truncation

Currently, the RT kernel differentiates with respect to 3 core variables `(τ, ϖ, Z)`,
and the δ-M truncation factor `fᵗ` is folded into these variables by `createAero` before
they enter the kernel. This means `createAero` must manually expand the chain rule for
δ-M truncation with respect to each aerosol sub-parameter.

The planned design promotes `fᵗ` to a **4th core variable**:

```
Current:   createAero → (τ_mod, ϖ_mod, Z) → RT kernel(3 core vars)
Planned:   Mie/profile → (τ_raw, ω̃, fᵗ, Z_raw) → RT kernel(4 core vars, δ-M inside)
```

This change:
- Moves δ-M truncation from the AD boundary into the analytic RT kernel
- Simplifies the chain rule: no manual δ-M derivative expansion in `createAero`
- Defines a `RawOpticalPropertyJacobian` at the AD boundary holding `∂(τ_raw, ω̃, fᵗ, Z_raw)/∂x`
- Makes it trivial to add new aerosol parameters or switch between analytic and ForwardDiff
  for the outer derivatives

## Directory Layout

```
src/
├── vSmartMOM.jl              Main module entry point
├── Architectures.jl          CPU/GPU abstraction
├── CoreRT/
│   ├── CoreRT.jl             Module entry point
│   ├── types.jl              Forward RT types
│   ├── types_lin.jl          Linearized RT types
│   ├── parameter_layout.jl   ParameterLayout for Jacobian indexing
│   ├── constants.jl          Scientific constants
│   ├── rt_run.jl             Forward RT entry point
│   ├── rt_run_lin.jl         Linearized RT entry point
│   ├── CoreKernel/           Elemental, doubling, interaction solvers
│   ├── LayerOpticalProperties/  Optical property assembly
│   ├── Surfaces/             Surface BRDF models
│   └── tools/                Helpers, postprocessing, model construction
├── Scattering/               Mie theory and phase functions
├── Absorption/               Gas absorption (HITRAN, ABSCO)
├── Inelastic/                Raman scattering
├── Aerosols/                 Aerosol microphysics
├── SolarModel/               Solar irradiance
└── IO/                       File I/O (YAML, NetCDF)

ext/
└── vSmartMOMCUDAExt.jl       CUDA extension (optional)

test/                         Unit tests (@testset)
sandbox/                      Development scripts, benchmarks, prototyping
docs/                         Documenter.jl documentation
```
