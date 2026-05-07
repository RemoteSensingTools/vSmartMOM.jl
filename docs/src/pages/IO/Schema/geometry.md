# `geometry` block

The observation geometry — solar zenith, viewing zenith(s) and
azimuth(s), and the observer altitude. All angles in degrees.

## Required fields

- **`sza`** — `Real`. Solar zenith angle [°]. Single scalar (no
  multi-SZA support yet — for that, run multiple `rt_run` calls). The
  solar beam is a single direction; multi-SZA is a Phase E.5 follow-up
  on top of the v0.6 source-term refactor.
- **`vza`** — `Vector{Real}`. Viewing zenith angles [°]. One entry per
  output line of sight. Each entry is appended to the quadrature as a
  zero-weight output node so the solver can interpolate exactly to the
  requested viewing direction (see [`conventions.md` §6](../../conventions.md#6-quadrature-streams-nstreams-vs-nquad)).
- **`vaz`** — `Vector{Real}`. Viewing azimuth angles [°], one per
  `vza`. Convention: `vaz = 0°` is the principal plane on the sun's
  side; `vaz = 180°` is the anti-solar principal plane. See
  [`conventions.md` §3](../../conventions.md#3-azimuth-angle-φ) for the
  Hovenier vs Mishchenko sign-difference notes.
- **`obs_alt`** — `Real`. Observer altitude in **Pa** (not hPa). For
  TOA observers use the top-of-atmosphere pressure (typically 1 Pa or
  smaller). For surface-based observers use a value larger than the
  atmospheric profile's surface pressure.

## Examples

### Single nadir TOA observation

```yaml
geometry:
  sza: 45.0
  vza: [0.0]
  vaz: [0.0]
  obs_alt: 1000.0     # Pa, near-TOA
```

### Multi-VZA, principal-plane sweep

```yaml
geometry:
  sza: 78.46          # acos(0.2), the Natraj-2009 setup
  vza: [10.0, 20.0, 40.0, 60.0]
  vaz: [0.0, 0.0, 0.0, 0.0]   # all in the principal plane
  obs_alt: 1000.0
```

### Out-of-plane multi-azimuth

```yaml
geometry:
  sza: 35.0
  vza: [10.0, 20.0, 40.0]
  vaz: [10.0, 90.0, 170.0]    # solar_tester vector setup
  obs_alt: 1000.0
```

## See also

- [`Schema/atmospheric_profile.md`](atmospheric_profile.md) — the
  pressure grid `obs_alt` is interpreted against
- [`conventions.md`](../../conventions.md) — Q/U/V signs and azimuth
  conventions vs VLIDORT / Mishchenko / Hovenier
