# `scattering` block (optional)

Mie aerosol configuration: per-aerosol size distribution, refractive
index, vertical distribution, and the Fourier decomposition method.
Skip this block entirely for clear-sky / Rayleigh-only / canopy-only
runs.

## Required when present

- **`aerosols`** — `Vector{Dict}`. One dictionary per aerosol species
  (see [Per-aerosol fields](#per-aerosol-fields) below).
- **`r_max`** — `Real` [µm]. Maximum particle radius for the size
  distribution quadrature. 30 µm covers most tropospheric Mie cases.
- **`nquad_radius`** — `Integer`. Number of radius quadrature points.
  Typical: 250–500.
- **`λ_ref`** — `Real` [µm]. Reference wavelength for the aerosol
  optical depth normalization.
- **`decomp_type`** — `String`. Fourier decomposition algorithm:
  - `"NAI2()"` — Siewert NAI-2 (default for most cases; cheaper).
  - `"PCW()"` — Domke PCW (precomputed Wigner symbols; uses more
    memory but cheaper for very large `nquad_radius`).

## Optional fields

- **`n_ref`** — `Complex`. Reference refractive index for the AOD
  normalization. Defaults to the first aerosol's `(nᵣ, nᵢ)` if
  omitted.

## Per-aerosol fields

Each entry in `aerosols:` is a dict. Two forms are supported.

### Mie/lognormal aerosol (full Mie computation)

| Field | Type | Description |
|-------|------|-------------|
| `τ_ref` | `Real` | Aerosol optical depth at `λ_ref`. |
| `μ` | `Real` | Geometric mean radius [µm]. |
| `σ` | `Real` | Geometric standard deviation. |
| `nᵣ` | `Real` | Real part of refractive index. |
| `nᵢ` | `Real` | **Negative** imaginary part (vSmartMOM uses `n = nᵣ - i·nᵢ`; see commit `b8227c2`). |

### Analytic phase-function aerosol (skip Mie computation)

| Field | Type | Description |
|-------|------|-------------|
| `τ_ref` | `Real` | Aerosol optical depth at `λ_ref`. |
| `phase_function` | `String` or `Dict` | Analytic phase function constructor. Allowed: `HenyeyGreensteinPhaseFunction(g=0.7)`, `SyntheticPolarizedHenyeyGreensteinPhaseFunction(g=0.45, polarization_fraction=0.6)`. |
| `ssa` *(optional)* | `Real` | Single-scattering albedo. Default `1.0`. Alias: `ϖ`. |
| `μ`, `σ`, `nᵣ`, `nᵢ` *(optional, all-or-none)* | `Real` | Mie microphysical metadata; retained but the analytic `phase_function` overrides the optics. |

The analytic path is intended for fast prototyping and synthetic
test scenes; linearized Jacobians are currently defined only for
Mie/lognormal aerosols (analytic phase functions throw in `LinMode`).

### Vertical distribution (required, exactly one form)

| Form | Fields | Use case |
|------|--------|----------|
| Altitude-form (Sanghavi convention, **preferred**) | `z₀`, `σ₀` | Log-normal in altitude. Generalizes cleanly across atmospheric profiles. |
| Pressure-form (legacy unified convention) | `p₀`, `σp` | Normal in pressure. The historical default before the sanghavi-unified merge. |

Both forms with both pairs is rejected.

## Examples

### Single Mie aerosol with altitude-form profile

```yaml
scattering:
  aerosols:
    - τ_ref: 0.05
      μ:    0.20      # geometric mean radius (µm)
      σ:    1.6
      nᵣ:   1.45
      nᵢ:   0.001
      z₀:   3.0       # altitude-form: log-normal centered at 3 km
      σ₀:   1.5
  r_max: 30.0
  nquad_radius: 300
  λ_ref: 0.55
  decomp_type: "NAI2()"
```

### Two-aerosol setup (fine + coarse mode)

```yaml
scattering:
  aerosols:
    - τ_ref: 0.10
      μ: 0.15
      σ: 1.5
      nᵣ: 1.45
      nᵢ: 0.001
      z₀: 2.0
      σ₀: 1.2
    - τ_ref: 0.02
      μ: 1.5
      σ: 2.0
      nᵣ: 1.50
      nᵢ: 0.005
      z₀: 5.0
      σ₀: 2.0
  r_max: 50.0
  nquad_radius: 500
  λ_ref: 0.55
  decomp_type: "PCW()"
```

### Henyey-Greenstein analytic phase function (no Mie)

```yaml
scattering:
  aerosols:
    - τ_ref: 0.5
      phase_function: "HenyeyGreensteinPhaseFunction(g=0.8)"
      ssa: 0.95
      p₀: 800.0     # pressure-form profile
      σp: 100.0
  r_max: 30.0
  nquad_radius: 300
  λ_ref: 0.55
  decomp_type: "NAI2()"
```

## Phase D interaction with `truncation: auto`

When `radiative_transfer.truncation = auto` (the v0.7 default):

- No `scattering` block ⇒ `NoTruncation()` (Rayleigh only fits any
  `stream_l_cap`).
- `scattering` block present ⇒ `δBGE(stream_l_cap, Δ_angle)` (Mie
  phase functions typically have hundreds of Greek moments).

The choice is logged at `@info` level so you always see what was
applied. To force `NoTruncation()` even with aerosols (e.g. for a
benchmark that requires no transform), set
`truncation: "NoTruncation()"` explicitly.

## See also

- [`Schema/radiative_transfer.md`](radiative_transfer.md) — `truncation`
  modes (`auto`, `NoTruncation`, `δBGE`)
- [`Schema/atmospheric_profile.md`](atmospheric_profile.md) — vertical
  grid that pressure-form profiles bind to
- Aerosol module reference: [`docs/src/pages/api/scattering.md`](../../api/scattering.md)
