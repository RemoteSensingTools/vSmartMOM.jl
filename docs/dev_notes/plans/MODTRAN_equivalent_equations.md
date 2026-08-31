# MODTRAN-equivalent fields from paired vSmartMOM runs

This note spells out how the six MODTRAN-style atmospheric correction fields
written by `test/benchmarks/modtran_equivalent_fields.jl` are derived from the
two vSmartMOM noRS runs produced by:

- `test/benchmarks/emit_modtran_noRS_scenarios.jl`      — Lambertian α = 0.0
- `test/benchmarks/emit_modtran_noRS_scenarios_alb02.jl` — Lambertian α = 0.2

Both drivers write the same per-scenario spectra:

| variable    | meaning                                                     |
|-------------|-------------------------------------------------------------|
| `R_tot`     | TOA reflected radiance at sensor VZA (Stokes I, F₀ = 1)     |
| `T_tot`     | BOA transmitted radiance at sensor VZA (Stokes I, F₀ = 1)   |
| `tau_total` | Total vertical atmospheric optical depth (abs + Ray + aer)  |
| `hemR_tot`  | Bi-hemispheric reflectance `bhr_uw` at TOA (Stokes-I flux)  |
| `hemT_tot`  | Bi-hemispheric transmittance `bhr_dw` at BOA (Stokes-I flux)|

The vSmartMOM solar source is normalised so that the downwelling TOA flux is
`F_TOA↓ = μ_s`; all radiances below are Stokes I with `F₀ = 1`. `bhr_dw` and
`bhr_uw` are the hemispheric flux integrals defined in
`src/CoreRT/CoreKernel/interaction_hdrf.jl:23`.

## Notation

| symbol            | meaning                                            |
|-------------------|----------------------------------------------------|
| `μ_s = cos(θ_s)`  | cosine of solar zenith angle                       |
| `μ_v = cos(θ_v)`  | cosine of view zenith angle (= `cos(180° − TSZ)`) |
| `τ`               | total vertical atmospheric optical depth           |
| `ρ_s`             | Lambertian surface albedo                          |
| `ρ_s2`            | second-albedo value (= 0.2 in this pipeline)       |
| `L(ρ_s)`          | `R_tot` (TOA radiance) at the given surface albedo |
| `bhr_dw(ρ_s)`     | `hemT_tot` (BOA hemispheric flux) at albedo `ρ_s`  |

## Governing relation

For a plane-parallel atmosphere above a Lambertian surface of albedo `ρ_s` the
TOA reflectance at cosine `μ_v`, solar cosine `μ_s`, is

```
ρ_TOA(ρ_s, μ_s, μ_v) = ρ_atm(μ_s, μ_v)
                      + T↓(μ_s) · ρ_s · T↑(μ_v) / (1 − ρ_s · S),         (1)
```

where

- `ρ_atm` — path reflectance (atmosphere over a black surface)
- `T↓(μ_s) = T↓_dir(μ_s) + T↓_dif(μ_s)` — total downwelling transmittance
- `T↑(μ_v) = T↑_dir(μ_v) + T↑_dif(μ_v)` — total upwelling transmittance
- `S` — atmospheric spherical albedo (bi-hemispheric reflectance from below)

The dimensionless reflectance corresponding to a radiance L (Stokes I, F₀ = 1) is
`ρ = π · L / μ_s`.

## Field-by-field equations

Let `0`-subscripts denote the α = 0 run and `2`-subscripts the α = 0.2 run. The
driver-produced spectra enter as follows.

### 1. Path reflectance `rhoatm`

At `ρ_s = 0`, equation (1) collapses to `ρ_TOA = ρ_atm`, so

```
rhoatm(λ) = π · L0(λ) / μ_s                                              (2)
         = π · R_tot(α=0, λ) / μ_s
```

Source fields: `alb0.R_tot`.

### 2. Direct-beam transmittances `transm_down_dir`, `transm_up_dir`

Beer–Lambert on the vertical column:

```
transm_down_dir(λ) = exp(−τ(λ) / μ_s)                                    (3a)
transm_up_dir  (λ) = exp(−τ(λ) / μ_v)                                    (3b)
```

Source field: `alb0.τ_total` (either run's `τ_total` is equivalent; the α = 0
JLD2 is used).

### 3. Atmospheric spherical albedo `sphalb`

Start from the BOA downward flux. With solar input, the series of surface
multiple-reflections gives

```
bhr_dw(ρ_s) = bhr_dw(0) · Σ_{k=0}^∞ (ρ_s · S)^k
            = bhr_dw(0) / (1 − ρ_s · S).                                  (4)
```

Solving (4) for `S` with the two available albedos:

```
sphalb(λ) = (1 − bhr_dw(0, λ) / bhr_dw(ρ_s2, λ)) / ρ_s2                  (5)
         = (1 − hemT_tot(α=0, λ) / hemT_tot(α=0.2, λ)) / 0.2
```

This is the multiple-scattering spherical albedo; no single-scatter proxy is
required. Source fields: `alb0.hemT_tot`, `alb02.hemT_tot`.

### 4. Total downwelling transmittance and its diffuse part

From the solar source normalisation, `bhr_dw(0) = μ_s · T↓(μ_s)`, so

```
T↓(μ_s, λ) = hemT_tot(α=0, λ) / μ_s                                      (6)
transm_down_dif(λ) = T↓(μ_s, λ) − transm_down_dir(λ)                     (7)
```

Source field: `alb0.hemT_tot`.

### 5. Total upwelling transmittance and its diffuse part

Use (1) with `ρ_s = ρ_s2`:

```
π · (L2 − L0) / μ_s = ρ_s2 · T↓(μ_s) · T↑(μ_v) / (1 − ρ_s2 · S).          (8)
```

Solving for `T↑(μ_v)` and substituting (5), (6):

```
T↑(μ_v, λ) = π · (L2(λ) − L0(λ)) · (1 − ρ_s2 · S(λ)) / (ρ_s2 · bhr_dw(0, λ))   (9)
transm_up_dif(λ) = T↑(μ_v, λ) − transm_up_dir(λ)                         (10)
```

Source fields: `alb0.R_tot`, `alb02.R_tot`, `alb0.hemT_tot`, plus `S` from (5).

## Exact code → equation mapping

Locations in `test/benchmarks/modtran_equivalent_fields.jl` (after the
2026-04-24 rewrite):

| file line                                  | equation                  |
|--------------------------------------------|---------------------------|
| `@. rho = π * L0 / μ_s`                    | (2) — `rhoatm`            |
| `@. S = (1 - hT0 / hT2) / RHO_S2`          | (5) — `sphalb`            |
| `@. tdw_dir = exp(-τ / μ_s)`               | (3a) — `transm_down_dir`  |
| `@. tup_dir = exp(-τ / μ_v)`               | (3b) — `transm_up_dir`    |
| `T_down_tot = hT0 ./ μ_s`                  | (6)                       |
| `@. tdw_dif = T_down_tot - tdw_dir`        | (7) — `transm_down_dif`   |
| `α_ρ = @. π * (L2 - L0) / μ_s`             | LHS of (8)                |
| `T_up_tot = @. α_ρ * (1 - RHO_S2*S) / (RHO_S2*T_down_tot)` | (9)     |
| `@. tup_dif = T_up_tot - tup_dir`          | (10) — `transm_up_dif`    |

The JLD2 keys consumed by the driver:

| code name in `modtran_equivalent_fields.jl` | JLD2 key from driver file     |
|----------------------------------------------|-------------------------------|
| `alb0.R_tot`                                 | `R_tot`     (α = 0 JLD2)     |
| `alb0.τ_total`                               | `τ_total`   (α = 0 JLD2)     |
| `alb0.hemT_tot`                              | `hemT_tot`  (α = 0 JLD2)     |
| `alb02.R_tot`                                | `R_tot`     (α = 0.2 JLD2)   |
| `alb02.hemT_tot`                             | `hemT_tot`  (α = 0.2 JLD2)   |
| `alb02.τ_total`                              | `τ_total`   (α = 0.2 JLD2, consistency check only) |

## Sanity bounds

- `sphalb` should lie in [0, 1]. Large `τ · μ_s^{-1}` leaves `bhr_dw(0) ≈
  bhr_dw(ρ_s2)` and `sphalb → 0`; small-τ clear atmospheres give `sphalb ~
  Rayleigh bihemispheric reflectance` (~0.08 at 550 nm).
- `transm_down_dir ≤ T↓(μ_s) ≤ 1`, so `transm_down_dif ≥ 0`. The driver still
  clamps to zero to guard against tiny negative residuals from finite-precision
  interpolation / matrix inversion.
- Beer-Lambert identity `ln(transm_up_dir)·μ_v = ln(transm_down_dir)·μ_s =
  −τ_total` holds exactly in vSmartMOM output by construction (see Part I of
  `plans/HANDOFF_LUT_REBUILD_2026-04-23.md`). MODTRAN breaks this by ~0.2 in
  AOT because of its 6S-style direct/diffuse partitioning; that difference is
  expected.

## Zero-flux guard

In deep line cores (e.g. O₂ A-band line centers, H₂O 1.94 μm overtone cores),
both `bhr_dw(0)` and `bhr_dw(ρ_s2)` collapse to machine zero and equations (5)
and (9) reduce to 0/0. The driver applies the guards

```
S(λ)        = 0                     if bhr_dw(ρ_s2, λ) ≤ eps(Float64)   (5′)
T↑(μ_v, λ)  = exp(−τ(λ)/μ_v)        if bhr_dw(0, λ)    ≤ eps(Float64)   (9′)
```

The `T↑` fallback makes `transm_up_dif = 0` inside the same saturated cells.
Both substitutions are physically meaningful: no photon reaches the surface,
so no multi-reflection channel is resolvable and any upwelling radiance must
be pure direct-beam.
