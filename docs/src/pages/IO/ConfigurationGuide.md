# Configuration Guide — Building a Scene, Step by Step

This guide builds a vSmartMOM scene up from the simplest possible setup,
adding one capability at a time. **Every example is shown in two equivalent
forms** — a config file (TOML or YAML) and the equivalent pure-Julia
construction with no file at all — because they map onto exactly the same
parser.

For the per-field reference tables, see the [Input Schema](Schema.md) pages
(one per block). This page is the *narrative*: what to reach for, in what
order, and what has to line up.

## Two ways to configure — they are the same thing

```
  config file (TOML/YAML)  ─┐
                            ├─► read_parameters ─► vSmartMOM_Parameters ─► model_from_parameters ─► rt_run
  in-script Dict           ─┘
```

`read_parameters` dispatches on its argument:

| Input | Call |
|-------|------|
| a `.toml` / `.yaml` path | `read_parameters("scene.toml")` |
| an in-memory `Dict` | `read_parameters(cfg)` |
| an existing `vSmartMOM_Parameters` | returned unchanged |

**The in-script `Dict` mirrors the file 1:1.** Each top-level block becomes a
`Dict` entry; the typed fields (surfaces, polarization type, architecture, …)
are given as *the same strings* the file uses — they are parsed identically.
So anything you can write in TOML you can build in a script, and vice versa,
without ever touching disk.

```julia
using vSmartMOM
params = read_parameters("scene.toml")     # from a file
# — or, identically, with no file: —
params = read_parameters(cfg)              # cfg::Dict, see every example below
model  = model_from_parameters(params)
R, T   = rt_run(model)
```

---

## Level 0 — Minimal pure-Rayleigh scalar scene

The smallest valid scene: one spectral point, a Lambertian surface, a
two-layer Rayleigh atmosphere, scalar (intensity-only) radiation, no
absorption, no aerosols.

**File (TOML):**
```toml
[radiative_transfer]
spec_bands        = ["[12987.0]"]            # one wavenumber, cm⁻¹
surface           = ["LambertianSurfaceScalar(0.15)"]
nstreams          = 3                          # ≥3 (Rayleigh reaches m=2)
polarization_type = "Stokes_I()"
truncation        = "NoTruncation()"
depol             = -1                          # -1 ⇒ standard Rayleigh depol
float_type        = "Float64"
architecture      = "CPU()"

[geometry]
sza     = 60.0
vza     = [60.0]
vaz     = [180.0]
obs_alt = [0.0]                                 # TOA + BOA radiances

[atmospheric_profile]
T = [250.0, 275.0]                              # K, full levels (layer centers)
p = [100.0, 500.0, 1000.0]                      # hPa, half levels (length = #layers + 1)
profile_reduction = -1                          # -1 ⇒ keep all layers
```

**In-script (Dict), exactly equivalent:**
```julia
cfg = Dict(
  "radiative_transfer" => Dict(
    "spec_bands"        => ["[12987.0]"],
    "surface"           => ["LambertianSurfaceScalar(0.15)"],
    "nstreams"          => 3,
    "polarization_type" => "Stokes_I()",
    "truncation"        => "NoTruncation()",
    "depol"             => -1,
    "float_type"        => "Float64",
    "architecture"      => "CPU()"),
  "geometry" => Dict(
    "sza" => 60.0, "vza" => [60.0], "vaz" => [180.0], "obs_alt" => [0.0]),
  "atmospheric_profile" => Dict(
    "T" => [250.0, 275.0], "p" => [100.0, 500.0, 1000.0],
    "profile_reduction" => -1))

params = read_parameters(cfg)
model  = model_from_parameters(params)
R, T   = rt_run(model)
```

**What each field means / what must line up**

- `spec_bands` — one entry per band. Two syntaxes: a single point/list
  `"[12987.0]"`, or a range `"(1e7/775):0.05:(1e7/755)"` (start:step:stop, in
  cm⁻¹). Wavelength↔wavenumber: `ν[cm⁻¹] = 1e7 / λ[nm]`. See the
  [`radiative_transfer` schema](Schema/radiative_transfer.md) (the `spec_bands` field).
- `surface` — **one BRDF per band**. `length(surface) == length(spec_bands)`.
- `nstreams` — angular resolution. Public contract: `stream_l_cap =
  2·nstreams − 1`. Minimum 3.
- `T` has length `#layers`; `p` (half levels) has length `#layers + 1`. They
  define the layering — get this off-by-one right.
- `obs_alt` is a height selection in **km above BOA**. `[0]` requests the
  standard TOA-upwelling and BOA-downwelling pair. Scalar `0` requests BOA
  only; a scalar interior `H` requests both directions at `H`; `[0, H]`
  requests TOA, BOA, and `H`. See the
  [`geometry` schema](Schema/geometry.md) for the complete scalar/vector
  convention and the named `ObserverRTResult` output.

---

## Level 1 — Add polarization

Switch the Stokes vector from intensity-only to full linear (or full)
polarization. Only `polarization_type` changes.

```toml
polarization_type = "Stokes_IQU()"     # I, Q, U  (n = 3)
# or "Stokes_IQUV()" for full Stokes (n = 4); "Stokes_I()" is scalar (n = 1)
```
```julia
cfg["radiative_transfer"]["polarization_type"] = "Stokes_IQU()"
```

The output `R`/`T` gain a Stokes dimension: `R[ivza, istokes, ispec]`. Nothing
else in the config changes; the solver propagates the larger phase matrix.

---

## Level 2 — Add trace-gas absorption

Introduce a new top-level `absorption` block. This is where line-by-line gas
optical depth enters.

```toml
[absorption]
fixed_molecules    = [["O2"]]          # one list per band; no Jacobian
variable_molecules = [[]]              # one list per band; Jacobian computed
vmr                = { O2 = 0.21 }     # volume mixing ratio per molecule
broadening         = "Voigt()"
CEF                = "HumlicekWeidemann32SDErrorFunction()"
wing_cutoff        = 40                # cm⁻¹ from line center
```
```julia
cfg["absorption"] = Dict(
  "fixed_molecules"    => [["O2"]],
  "variable_molecules" => [[]],
  "vmr"                => Dict("O2" => 0.21),
  "broadening"         => "Voigt()",
  "CEF"                => "HumlicekWeidemann32SDErrorFunction()",
  "wing_cutoff"        => 40)
```

### Deep dive — the `absorption` block

**Per-band molecule lists.** `fixed_molecules[ib]` and
`variable_molecules[ib]` each hold the molecules active in band `ib`. There is
**one list per spectral band** — band 1 might carry `O2`, band 2 `CO2`, etc.

- **fixed** — abundance held constant, *no* Jacobian computed.
- **variable** — abundance is a state-vector element, Jacobian computed (the
  linearized path differentiates these).

**What must match**
- `length(fixed_molecules) == length(variable_molecules) == length(spec_bands)`.
- Every molecule that appears in either list **must have a `vmr` entry**, or
  parsing errors.
- Each molecule must actually have lines inside its band's wavenumber range —
  a molecule with no lines in-band contributes nothing.

**VMR sources.** A `vmr` value is either
- a **scalar** (e.g. `O2 = 0.21`) — constant with height, or
- a **vector** with either one value per layer (`N`) or one value per pressure
  interface (`N+1`). Interface values are mapped to layer centers in log
  pressure before any observer-height interfaces are inserted.

**Line shape.** `broadening` ∈ `Voigt()` (default; Doppler+pressure, the right
choice through the troposphere/stratosphere), `Lorentz()` (pressure limit),
`Doppler()` (thermal limit). `CEF` is the complex-error-function evaluator for
Voigt (`HumlicekWeidemann32SDErrorFunction()`). `wing_cutoff` (cm⁻¹) drops line
contributions farther than that from each spectral point.

**Where cross-sections come from.** With no `LUTfiles`, lines are pulled from
HITRAN on the fly. With `LUTfiles` (one vector per band, parallel to the
molecule lists) the precomputed look-up tables are used instead — much faster,
and the production pipeline path. `cia_files` / `mtckd_file` add
collision-induced and continuum absorption.

> **Legacy form.** Older configs use a single `molecules = [[O2], [CO2]]`
> list (all treated as fixed) instead of `fixed_molecules`/`variable_molecules`.
> It still parses for non-H₂O species; prefer the explicit split for new work.

---

## Level 3 — Add water vapor via specific humidity `q`

H₂O is special: in the modern `fixed_molecules`/`variable_molecules` form it is
**never listed**. Instead it is derived from the specific-humidity profile `q`
on the `atmospheric_profile` block:

```toml
[atmospheric_profile]
T = [250.0, 275.0]
p = [100.0, 500.0, 1000.0]
q = [0.001, 0.005]        # specific humidity, per full level (mass fraction)
```
```julia
cfg["atmospheric_profile"]["q"] = [0.001, 0.005]
```

vSmartMOM converts `q` to an H₂O volume mixing ratio per layer:

```
vmr_h2o = q / (1 − q) · (M_dry / M_H₂O)
```

and adds the H₂O line absorption automatically (HITRAN on the fly, or the H₂O
LUT if supplied). `q` may contain `N` layer-center values (same as `T`) or
`N+1` pressure-interface values; interface values are mapped to layer centers
in log pressure. All-zero `q` (the default) means a dry atmosphere.

> **H₂O is driven only by `q`.** Do not list `H2O` in `fixed_molecules`,
> `variable_molecules`, or the legacy `molecules` form; the parser rejects it
> to prevent double-specifying water vapor.

---

## Level 4 — Multiple bands

Bands are independent: each gets its own spectral range, surface, molecule
lists, and aerosol Z-matrices. Just make every per-band list the same length.

```toml
[radiative_transfer]
spec_bands = ["(1e7/775):0.05:(1e7/755)", "(1e7/1640):0.1:(1e7/1590)"]
surface    = ["LambertianSurfaceScalar(0.15)", "LambertianSurfaceScalar(0.30)"]

[absorption]
fixed_molecules    = [["O2"], ["CO2"]]
variable_molecules = [[],     ["CO2"]]
vmr                = { O2 = 0.21, CO2 = 400e-6 }
```
```julia
cfg["radiative_transfer"]["spec_bands"] =
  ["(1e7/775):0.05:(1e7/755)", "(1e7/1640):0.1:(1e7/1590)"]
cfg["radiative_transfer"]["surface"] =
  ["LambertianSurfaceScalar(0.15)", "LambertianSurfaceScalar(0.30)"]
cfg["absorption"]["fixed_molecules"]    = [["O2"], ["CO2"]]
cfg["absorption"]["variable_molecules"] = [[], ["CO2"]]
cfg["absorption"]["vmr"]                = Dict("O2" => 0.21, "CO2" => 400e-6)
```

**What must match:** every per-band list — `surface`, `fixed_molecules`,
`variable_molecules`, and (if used) `LUTfiles` — must have one entry per band.

---

## Level 5 — Add aerosols (Mie scattering)

Add a `scattering` block. Each aerosol is a size distribution + refractive
index + a vertical placement + a reference optical depth.

```toml
[scattering]
r_max        = 50.0       # µm, max radius for the size-distribution quadrature
nquad_radius = 1000       # radius quadrature points
λ_ref        = 0.770      # µm, wavelength at which τ_ref is defined
decomp_type  = "NAI2()"

[[scattering.aerosols]]
τ_ref = 0.001             # aerosol optical depth at λ_ref
μ     = 1.0               # geometric mean radius (µm), lognormal
σ     = 1.5               # geometric std dev (≥ 1)
nᵣ    = 1.3               # real refractive index
nᵢ    = 1.0e-8            # imaginary (absorbing) part
p₀    = 700.0             # vertical peak pressure (hPa)
σp    = 50.0              # vertical width (hPa)
```
```julia
cfg["scattering"] = Dict(
  "r_max" => 50.0, "nquad_radius" => 1000, "λ_ref" => 0.770,
  "decomp_type" => "NAI2()",
  "aerosols" => [Dict(
    "τ_ref" => 0.001, "μ" => 1.0, "σ" => 1.5,
    "nᵣ" => 1.3, "nᵢ" => 1.0e-8,
    "p₀" => 700.0, "σp" => 50.0)])
```

### Deep dive — the `scattering` block

**Per-aerosol fields**

| Field | Meaning |
|-------|---------|
| `τ_ref` | optical depth of this aerosol at `λ_ref` (the normalization) |
| `μ`, `σ` | lognormal size distribution — geometric mean radius (µm) and geometric std dev (`σ ≥ 1`) |
| `nᵣ`, `nᵢ` | complex refractive index `nᵣ − i·nᵢ`; `nᵢ ≥ 0` is absorption |
| `p₀`, `σp` | **pressure-form** vertical profile: a normal in pressure, peak `p₀` (hPa), width `σp` |
| `z₀`, `σ₀` | **altitude-form** vertical profile (preferred): lognormal in altitude (km) |

**Vertical placement — exactly one form.** Give *either* `(p₀, σp)` *or*
`(z₀, σ₀)`, never both and never neither. The profile shape multiplies `τ_ref`
into a per-layer optical depth.

**Band-level scattering knobs**

- `r_max`, `nquad_radius` — the radius quadrature `Mie` integrates over. Larger
  `nquad_radius` resolves the size distribution better at more cost.
- `λ_ref` — the wavelength `τ_ref` is anchored to.
- `decomp_type` — phase-matrix Fourier decomposition: `NAI2()` (cheaper, GPU-
  capable, the default) or `PCW()` (precomputed Wigner; faster for very large
  `nquad_radius`, CPU-only).

**Where the Mie work happens.** For each band and aerosol, Lorenz–Mie theory is
integrated over the radius quadrature to produce the extinction, single-
scattering albedo `ω̃`, and the `GreekCoefs` (phase-matrix moments). These feed
the layer optical properties (the per-layer `τ`, `ϖ`, `Z⁺⁺`, `Z⁻⁺`).
Truncation (`truncation = "auto"`) applies δBGE forward-peak removal only when
the Greek series is longer than `stream_l_cap`.

### Spectral resolution — what varies within a band, and what doesn't

Not every optical property is resolved at every spectral point. The forward and
linearized paths are **consistent** here:

| Property | Resolution (forward = linearized) |
|----------|-----------------------------------|
| Gas absorption `τ_abs(ν)` | **per wavelength** (line-by-line) |
| Rayleigh `τ_rayl(ν)` (σ ∝ ν⁴) | **per wavelength** |
| Aerosol extinction `k` & AOD `τ_aer` | **per wavelength** (endpoint interpolation, see below) |
| Aerosol single-scattering albedo `ω̃` | **per wavelength** (endpoint interpolation) |
| Aerosol phase matrix `Z⁺⁺/Z⁻⁺` (Greek coefs) | **per band** (endpoint average — *not* per wavelength) |

**Absorption and Rayleigh are fully spectral** — each is stored as an
`(nSpec × nLayers)` matrix, recomputed at every wavenumber.

**Aerosol AOD is anchored at `λ_ref` and re-scaled per wavelength.** `τ_ref` is
the optical depth at `λ_ref`; the model rescales it by the Mie extinction ratio:

```
τ_aer(λ, z) = τ_ref · ( k(λ) / k(λ_ref) ) · vertical_profile(z)
```

To get `k(λ)` and `ω̃(λ)` without a full Mie calculation at every spectral point,
the aerosol Mie is evaluated at the two **band endpoints** and `k`/`ω̃` are
**linearly interpolated** in wavenumber onto each spectral point — so AOD does
vary across the band. The phase matrix (Greek coefficients) is the **average of
the two endpoints** — a single phase matrix per band, *not* recomputed per
wavelength. A single-wavelength band uses the band-center value directly.

> **Approximation & future work.** Linear endpoint interpolation assumes the
> aerosol optics vary smoothly and near-linearly across the band — true for
> bands narrow compared to the aerosol's spectral structure, which is **a key
> reason to keep spectral bands small**. A future refinement could add the
> band-center Mie point and use **quadratic (3-point) interpolation**, and could
> resolve the phase matrix per wavelength for wide bands. (See the code comment
> at the aerosol Mie call site in `model_from_parameters.jl`.)

### Refractive index: in-band `nᵣ/nᵢ` vs. reference `n_ref`

Two refractive-index roles exist in the aerosol block:

| Symbol | Config key | Role |
|--------|-----------|------|
| `nᵣ − i·nᵢ` | per-aerosol `nᵣ`, `nᵢ` | **In-band index** — used for all Mie calculations within the spectral band (`k(λ)`, `ω̃(λ)`, Greek coefs). |
| `n_ref` | `[scattering].n_ref` (optional) | **Reference index** — used *only* to compute `k_ref = k(λ_ref, n_ref)`, the extinction at the normalisation wavelength that anchors `τ_ref`. |

**If `n_ref` is not set, it defaults to the aerosol's own `nᵣ/nᵢ`**, so in-band n == reference n and the distinction is irrelevant. Setting `n_ref` is only needed when you want `τ_ref` to be normalised at a refractive index different from the aerosol's spectral value (e.g. a standard atmospheric convention).

> **No in-band n(λ) dispersion.** The refractive index `nᵣ/nᵢ` is constant
> across the band — wavelength-dependent dispersion of the refractive index is
> not yet supported. Each aerosol uses a single (nᵣ, nᵢ) pair for all Mie
> calculations within the band.

---

## Level 6 — Analytic phase-function aerosols (fast prototyping)

Skip Mie entirely with an analytic phase function — handy for synthetic scenes
and quick sensitivity tests.

```toml
[[scattering.aerosols]]
τ_ref          = 0.1
phase_function = "HenyeyGreensteinPhaseFunction(g=0.7)"
ssa            = 0.95
p₀ = 700.0
σp = 50.0
```

`ssa` (single-scattering albedo) is required here because there is no Mie
computation to derive it from. **Caveat:** analytic phase-function aerosols are
not supported in the linearized (Jacobian) path.

---

## Level 7 — Non-Lambertian surfaces

Swap the per-band `surface` string for any supported BRDF; the parser builds
the matching surface type.

```toml
surface = ["rpvSurfaceScalar(0.1, 0.05, 0.6, -0.1)"]   # RPV land
# "RossLiSurfaceScalar(fiso, fvol, fgeo)"               — Ross-Li (MODIS-style)
# "CoxMunkSurface(wind_speed, n_water)"                 — ocean glint
```

Each BRDF has its own parameters (see [surface](Schema/surface.md)). They are
per-band like everything else.

---

## Level 8 — Canopy

Wrap the band surface(s) in a structural canopy (leaves over soil) by adding a
`canopy` block.

```toml
[canopy]
LAI                = 3.0
n_layers           = 1
leaf_reflectance   = 0.4
leaf_transmittance = 0.05
soil               = "from_surface"   # reuse the band BRDF as understory
```
```julia
cfg["canopy"] = Dict(
  "LAI" => 3.0, "n_layers" => 1,
  "leaf_reflectance" => 0.4, "leaf_transmittance" => 0.05,
  "soil" => "from_surface")
```

The canopy wraps each band's BRDF in a `CanopySurface`; the soil below is the
true lower boundary. See the [`surface` schema](Schema/surface.md) (the
"Composite (canopy + soil)" section) for `CanopySurface` configuration.

---

## Level 9 — Raman (inelastic) scattering

Raman is selected through the `RS_type` passed to the run (the elastic config
above is unchanged). The pure-elastic default is `noRS()`; rotational Raman is
`RRS(...)`. See the [Inelastic concepts](../concepts/08_inelastic.md) page for
the source-term machinery and how the Cabannes/Rayleigh split feeds it.

---

## Level 10 — Jacobians (linearized RT)

Same config; build the model in linearized mode and run the linearized solver.
Molecules listed under `variable_molecules`, the aerosol parameters, and the
surface parameters become the differentiable state.

```julia
params           = read_parameters(cfg)
model, lin_model = model_from_parameters(LinMode(), params)
NAer  = length(params.scattering_params.rt_aerosols)
NGas  = size(lin_model.τ̇_abs[1], 1)
NSurf = 1
R, T, dR, dT = rt_run(model, lin_model, NAer, NGas, NSurf)
```

(Note: the linearized path is pure-elastic — Raman is not linearized.)

---

## Level 11 — Precision, GPU, and solver knobs

```toml
[radiative_transfer]
float_type   = "Float32"        # or "Float64"
architecture = "GPU()"          # or "CPU()"

[radiative_transfer.numerics]
blas_threads = 4                # cap BLAS threads for big spectral batches
verbose      = false
```
```julia
cfg["radiative_transfer"]["float_type"]   = "Float32"
cfg["radiative_transfer"]["architecture"] = "GPU()"
cfg["radiative_transfer"]["numerics"] =
  Dict("blas_threads" => 4, "verbose" => false)
```

`Float32` roughly halves memory and speeds the kernels at the cost of per-pixel
precision; `GPU()` dispatches the adding-doubling kernels to CUDA when
available. `numerics.blas_threads` caps BLAS threads (4–8, or 1 for large
spectral batches).

---

## Where to go next

- Per-block field reference: [Input Schema](Schema.md).
- Editor autocomplete + validation for TOML/YAML: [Schema editor support](Schema.md#editor-support--autocomplete--inline-validation).
- Worked end-to-end runs: [IO Examples](Examples.md) and the
  [Tutorials](../tutorials/Tutorial_IO.md).
