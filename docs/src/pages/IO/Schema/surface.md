# `surface` field (in `radiative_transfer`)

The bidirectional reflectance (BRDF) at the lower boundary, one per
spectral band. Each entry is a constructor-call-shaped string parsed
safely (no `eval` — a small whitelist of constructors).

## Field signature

```yaml
radiative_transfer:
  surface:
    - "<surface_constructor>"   # one entry per spec_band
```

The `surface` vector must have **the same length as `spec_bands`**.

## Allowed constructors

### Lambertian (m=0 only — diffuse, no Fourier dependence)

| Constructor | Use case |
|---|---|
| `LambertianSurfaceScalar(albedo::Real)` | Single albedo for all spectral points in the band. The simplest BRDF. Phase C trait: `m_max = 0`. |
| `LambertianSurfaceSpectrum([albedos...])` | Vector of albedos, one per spectral point in the band. Use when the surface albedo varies appreciably across the band (e.g. snow in NIR). |
| `LambertianSurfaceLegendre([a₀, a₁, a₂, ...])` | Legendre expansion in the spectral band's normalized coordinate. |
| `LambertianSurfaceSpline(λ, ρ)` | Spline through tabulated `(λ, ρ)` pairs. For full spectral coverage use the `LambertianSurfaceSpline` helper directly in Julia rather than YAML.|

### Vegetation BRDFs (full Fourier; trait: `m_max = user_l_cap`)

| Constructor | Source | Use case |
|---|---|---|
| `rpvSurfaceScalar(ρ₀, ρ_c, k, Θ)` | Rahman-Pinty-Verstraete | Vegetation BRDF; common in canopy retrievals. |
| `RossLiSurfaceScalar(fvol, fgeo, fiso)` | Roujean / Wanner | MODIS BRDF kernel decomposition. Three-parameter linear model.|

### Ocean (full Fourier + sun-glint; trait: `m_max = user_l_cap`)

| Constructor | Use case |
|---|---|
| `CoxMunkSurface(wind_speed=U, n_water=...)` | Cox-Munk wave-slope statistics + Fresnel + (optional) whitecaps + Smith shadowing. Full polarization. |

`CoxMunkSurface` keyword arguments:

- `wind_speed` *(required, m/s)* — 10 m wind. Drives the Cox-Munk
  variance via Cox & Munk (1954) `σ² = 0.003 + 0.00512·U`.
- `n_water` — complex refractive index. Default `nothing` triggers the
  built-in Segelstein 1981 spectral lookup.
- `whitecap_albedo` *(default `0.22`)* — Lambertian whitecap albedo
  (Koepke 1984).
- `include_whitecaps` *(default `true`)* — whether to include the
  whitecap-fraction Lambertian addend.
- `shadowing` *(default `true`)* — Smith (1967) shadow masking.

### Composite (canopy + soil)

`CanopySurface` is configured via a top-level `canopy:` block, not as
a string in `surface:`. The parser wraps each band's surface entry as
the soil BRDF inside the canopy when `canopy:` is present. See
`Schema/canopy.md` (forthcoming).

## Examples

### Lambertian land

```yaml
radiative_transfer:
  spec_bands: ["[18867 18870]"]
  surface:    ["LambertianSurfaceScalar(0.05)"]
```

### Cox-Munk ocean, default Segelstein water

```yaml
radiative_transfer:
  spec_bands: ["[18867 18870]"]
  surface:    ["CoxMunkSurface(wind_speed=5.0)"]
```

### Multiband (one BRDF per band)

```yaml
radiative_transfer:
  spec_bands:
    - "[20000 20003]"
    - "[12500 12503]"
  surface:
    - "LambertianSurfaceScalar(0.05)"
    - "LambertianSurfaceScalar(0.40)"   # higher NIR albedo for soil
```

## Phase C trait dispatch

Each surface declares its `component_m_max`:

| Surface | `component_m_max(::T, ctx)` |
|---|---|
| `LambertianSurface*` | `0` (m=0 is exact) |
| `CoxMunkSurface` | `ctx.user_l_cap` (no scheme-imposed cap) |
| `rpvSurfaceScalar` / `RossLiSurfaceScalar` | `ctx.user_l_cap` |
| `CanopySurface` | `ctx.user_l_cap` |

A Cox-Munk + Lambertian band loops to the user-set `stream_l_cap`; a
Lambertian-only Rayleigh band loops only to `m=0:2` (Rayleigh trait
`= 2`). See [`docs/dev_notes/fourier_stream_resolution_plan.md`](https://github.com/cfranken/vSmartMOM.jl/blob/main/docs/dev_notes/fourier_stream_resolution_plan.md)
for the full trait table.

## See also

- [`Schema/radiative_transfer.md`](radiative_transfer.md) — `nstreams`,
  `truncation`, `polarization_type` (the surface choice often drives
  these knobs)
- [`docs/src/pages/conventions.md`](../../conventions.md) — sign
  conventions for VLIDORT cross-validation
