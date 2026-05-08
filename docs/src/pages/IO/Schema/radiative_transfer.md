# `radiative_transfer` block

Top-level configuration for the RT solver — quadrature, polarization,
truncation, and the projection-cap knobs introduced in v0.7. The full
file at `src/CoreRT/DefaultParameters.yaml` shows a representative
config; use this page to look up what each field does, what the
defaults are, and which legacy aliases still work.

## Required fields

- **`spec_bands`** — `Vector{String}`. Each string is a Julia range
  expression in cm⁻¹, e.g. `"(1e7/777):0.015:(1e7/757)"`. Multiple
  bands run as a band-batched solve.
- **`surface`** — `Vector{String}`. One BRDF constructor per spectral
  band. See [`Schema/surface.md`](surface.md) for the surface vocabulary
  (`LambertianSurfaceScalar`, `CoxMunkSurface`, `rpvSurfaceScalar`, etc.).
- **`polarization_type`** — `String`. One of `"Stokes_I()"`,
  `"Stokes_IQ()"`, `"Stokes_IQU()"`, `"Stokes_IQUV()"`. Sets the Stokes
  vector dimension (1/2/3/4 channels respectively).
- **`Δ_angle`** — `Real`, **optional** (default `0.0`). Forward-peak
  exclusion angle in degrees, consumed by `δBGE` truncation. Ignored
  when `truncation: NoTruncation()` or when `truncation: auto` resolves
  to `NoTruncation()`. Legacy YAMLs commonly set `2.0` to reproduce
  VLIDORT's `DO_DELTAM_SCALING` default — set explicitly when matching
  a benchmark that pins this value.
- **`depol`** — `Real`. Rayleigh depolarization factor (0 = no
  depolarization, the Natraj-2009 idealization).
- **`float_type`** — `String`. `"Float64"` or `"Float32"`. Float32 is
  fine for most retrievals and gives a ~2× speedup on large spectral
  windows; Float64 is required for reproducing published benchmarks at
  full precision.
- **`architecture`** — `String`. `"Architectures.CPU()"`,
  `"Architectures.GPU()"`, or `"Architectures.MetalGPU()"`. The CUDA
  GPU path is loaded via the `vSmartMOMCUDAExt` weak-dependency
  extension (`using CUDA` activates it).

## Resolution knobs (Phase D — v0.7)

Pick **one** of the schemas below. The new `nstreams` knob is the
recommended path; legacy `max_m` + `l_trunc` configs keep working.

### New schema (recommended)

- **`nstreams`** — *optional, default `8`*. **Weighted streams per
  hemisphere**, the primary user-facing resolution knob. Public
  contract: `stream_l_cap = 2·nstreams - 1` regardless of
  `quadrature_type`. For VLIDORT users: this is the Julia equivalent
  of VLIDORT's `NSTREAMS`.
- **`m_max`** — *optional, default `null`*. Explicit Fourier loop
  bound (order). When set, clamps the trait-aggregator output below
  `stream_l_cap`. Useful when you want to truncate the loop more
  aggressively than the per-component traits would.
- **`truncation`** — *optional, default `auto`*. Controls how the
  aerosol phase function is treated relative to the projection cap:
  - `"auto"` *(default)*: VLIDORT-`DO_DELTAM_SCALING`-style. With no
    aerosols, resolves to `NoTruncation()`. With aerosols, picks
    `δBGE(stream_l_cap, Δ_angle)` and the per-band Mie loop applies
    it only when `length(greek.β) > stream_l_cap` — short phase
    functions automatically fall back to `NoTruncation()` so the
    truncation step can never crash on a fits-the-cap series.
    Logs the choice via `@info`.
  - `"NoTruncation()"` or `null`: exactly no transform. Coefficients
    above `stream_l_cap` are silently dropped from the Fourier
    projection (per Phase C trait clamps).
  - `"δBGE(L, Δ)"`: explicit δ-BGE-fit truncation (Sanghavi & Stephens
    2015). Used for benchmark reproducibility when you want exact
    parity with a published numerical convention.

### Legacy schema (still supported, deprecated v0.7)

- **`max_m`** — `Integer`. Number of Fourier moments (count, loop runs
  `m = 0:max_m-1`). Triggers legacy parsing — `nstreams` will be left
  unset and the old `min(ceil((l_max+1)/2), max_m)` aggregator
  semantics apply.
- **`l_trunc`** — `Integer`. Legendre cutoff. Plays the role of
  `stream_l_cap` under the legacy schema.

Setting **both** `nstreams` and `l_trunc` produces a `@warn` — the
legacy `l_trunc` is ignored and `nstreams` wins.

## Other optional fields

- **`quadrature_type`** — *optional, default `"GaussLegQuad()"`*. The
  recommended default. `"RadauQuad()"` is supported but documented as
  expert/legacy: it includes the SZA as a quadrature node (handy for
  DNI-on-node configurations) at the cost of being 5–50× less accurate
  per stream than Gauss-Legendre on Rayleigh-only setups (see the
  Natraj-2009 commit `f9403eb` and the v0.7 release notes). Sanghavi's
  guidance: keep Radau as an option, but reach for Gauss by default.
- **`numerics`** — *optional*. Sub-block with rarely-touched solver
  knobs:
  - `dτ_max_threshold` *(default ≈ 0.0125)* — doubling-step cap.
  - `dτ_min_floor`     *(default ≈ 1e-7)*  — minimum layer thickness.
  - `blas_threads`     *(default `null`)*  — per-model BLAS thread
    cap. `null` leaves `BLAS.get_num_threads()` alone; an integer
    pins it for the duration of `rt_run`. See commit `37390d2`.

## Example — minimal new schema

A new-schema config can omit `nstreams` (default `8`),
`quadrature_type` (default `GaussLegQuad`), and `truncation` (default
`auto`). The leanest valid `radiative_transfer` block is:

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
```

This expands to: `nstreams = 8` (so `stream_l_cap = 15`),
`quadrature_type = GaussLegQuad()`, `truncation = auto`. For a
Rayleigh-only Lambertian scene this resolves to `NoTruncation()` and
the Fourier loop runs `m = 0:2` (Phase C trait dispatch).

## Example — Cox-Munk ocean retrieval

```yaml
radiative_transfer:
  spec_bands: ["[18867 18870]"]
  surface: ["CoxMunkSurface(wind_speed=5.0)"]
  polarization_type: Stokes_IQU()
  nstreams: 16          # 16 weighted streams ⇒ stream_l_cap = 31
  truncation: auto      # auto → δBGE(31, 2.0) once Mie aerosol is added
  Δ_angle: 2.0
  depol: 0.0
  float_type: Float64
  architecture: Architectures.CPU()
```

## See also

- [`Schema/surface.md`](surface.md) — BRDF vocabulary
- [`Schema/aerosols.md`](aerosols.md) — Mie size distribution + refractive index
- [`Schema/sources.md`](sources.md) — `SolarBeam`, `SurfaceSIF`, `BlackbodySource`
- [`docs/src/pages/conventions.md`](../../conventions.md) §6 — `Nstreams` vs `Nquad` distinction
- [`docs/dev_notes/fourier_stream_resolution_plan.md`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/docs/dev_notes/fourier_stream_resolution_plan.md) — design memo for the v0.7 schema migration
