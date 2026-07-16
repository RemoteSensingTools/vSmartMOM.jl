# `geometry` block

The observation geometry: solar zenith, viewing zenith(s), viewing
azimuth(s), and the vertical levels at which radiances are requested. Angles
are in degrees; observer heights are in kilometres above the model bottom of
atmosphere (BOA).

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
- **`obs_alt`** — `Union{Real, Vector{Real}}`. A scalar height or a
  non-empty vector of heights, in **km above BOA**. Values must be finite
  and nonnegative. Scalar and vector forms deliberately have different
  endpoint semantics, as shown below.

## Observer-height convention

Let `H_TOA` be the geometric altitude of the model's top interface. It is
derived hydrostatically from the supplied pressure, temperature, and humidity
profiles.

| Input | Radiances returned |
|---|---|
| `0` | BOA downwelling only |
| `H`, where `0 < H < H_TOA` | Upwelling and downwelling at height `H` only |
| `H`, where `H >= H_TOA` | TOA upwelling only |
| `[0]` | TOA upwelling and BOA downwelling |
| `[H1, H2, ...]` | Upwelling and downwelling at every strict-interior height; any value at or above `H_TOA` also selects TOA |
| `[0, H1, H2, ...]` | TOA, BOA, and every strict-interior requested height |

Thus `[0]` is the normal full-column setting and preserves the historical
two-output `R, T = rt_run(model)` workflow. A zero inside a vector is an
endpoint sentinel: it selects **both** TOA and BOA. By contrast, scalar `0`
selects BOA only.

When migrating a scene that previously used `obs_alt: 0` or a pressure-like
value such as `obs_alt: 1000.0` merely to obtain the standard endpoint pair,
change it to `obs_alt: [0]`. Values are no longer interpreted as pressure.

Duplicate interior heights are collapsed. Interior results are ordered from
the highest requested altitude to the lowest, independently of input order.

## How interior heights affect the atmosphere

Every strict-interior height `H` becomes an exact atmospheric half-level
interface. If needed, vSmartMOM splits the layer containing `H`, maps the
height to pressure by interpolation in log pressure, and reframes the complete
vertical state onto the new grid. Temperature and the normalized layer-centred
humidity/VMR fields are interpolated in log pressure; humidity and vector VMRs
supplied on pressure interfaces are first mapped to layer centers. Full-level
pressures, water VMR, dry/wet column densities, and layer thicknesses are recomputed.
Rayleigh, gas, aerosol, and other derived layer properties are subsequently
built on this final grid.

Profile reduction retains every requested interior height as an interface. If
`K` distinct interior heights were requested but `profile_reduction` asks for
fewer than `K + 1` layers, vSmartMOM raises the effective layer count to
`K + 1` and emits a warning containing the requested heights and effective
count. Endpoint selections (`0` and values at or above TOA) do not consume an
interior interface.

## Examples

### Default TOA + BOA outputs

```yaml
geometry:
  sza: 45.0
  vza: [0.0]
  vaz: [0.0]
  obs_alt: [0]
```

### BOA only

```yaml
geometry:
  sza: 45.0
  vza: [0.0]
  vaz: [0.0]
  obs_alt: 0
```

### One interior level only

```yaml
geometry:
  sza: 45.0
  vza: [0.0]
  vaz: [0.0]
  obs_alt: 5.0       # upwelling + downwelling at 5 km above BOA
```

### TOA + BOA + two interior levels

```yaml
geometry:
  sza: 35.0
  vza: [10.0, 20.0, 40.0]
  vaz: [10.0, 90.0, 170.0]
  obs_alt: [0, 2.0, 8.0]
```

## Reading height-aware results

Forward `rt_run` returns an `ObserverRTResult` with named endpoint and
interior-level outputs:

```julia
result = rt_run(model)

result.toa                 # TOA upwelling, or nothing when not requested
result.boa                 # BOA downwelling, or nothing when not requested
result.toa_altitude_km

for level in result.levels
    @show level.height_km
    up = level.upwelling
    diffuse_down = level.downwelling
    direct_down = level.unscattered_downwelling
    total_down = total_downwelling(level)
end
```

Each `LevelRadiance` also contains `boundary_index`,
`inelastic_upwelling`, and `inelastic_downwelling`. The unscattered field is
nonzero only when the requested VZA resolves to the solar ordinate; it uses
the same VAZ-independent m=0 carrier and `1/(2π)` normalization as the
historical BOA output. The diffuse and direct fields remain separate so the
collimated beam is never silently folded into the MOM diffuse solution.

`ObserverRTResult`
retains the historical seven-slot iteration/indexing order
`(toa, boa, inelastic_toa, inelastic_boa, hdr, bhr_uw, bhr_dw)`. Therefore,
when `obs_alt: [0]`, existing code can continue to write:

```julia
R, T = rt_run(model)
```

Use the named `levels` records for interior-height outputs. At present,
forward interior-height runs support `SolarBeam` and `NoSource`; thermal and
surface-emission sources require multisensor source-slot propagation and
produce a clear error. `CanopySurface` also produces a clear error because its
spectral setup and atmosphere-canopy interleaving are not represented in the
legacy forward multisensor kernel.

### Reading linearized height-aware results

The analytic full-multiple-scatter path uses the same height convention and
returns an `ObserverRTResultLin`:

```julia
model, lin_model = model_from_parameters(LinMode(), params)
NAer = isnothing(params.scattering_params) ? 0 :
       length(params.scattering_params.rt_aerosols)
NGas = size(lin_model.τ̇_abs[1], 1)
NSurf = 1

result_lin = rt_run_lin(model, lin_model, NAer, NGas, NSurf)

result_lin.toa
result_lin.boa
result_lin.toa_jacobian
result_lin.boa_jacobian
result_lin.layout

for level in result_lin.levels
    up = level.upwelling
    down = level.downwelling
    direct = level.unscattered_downwelling

    dup = level.upwelling_jacobian
    ddown = level.downwelling_jacobian
    ddirect = level.unscattered_downwelling_jacobian
    dtotal = total_downwelling_jacobian(level)
end
```

Each `LevelRadianceLin` Jacobian appends the parameter dimension to the
corresponding forward array, giving shape
`nVZA × nStokes × nSpec × nParams`. Its last axis follows
`result_lin.layout`. The unscattered-beam tangent is reported separately from
the diffuse tangent and obeys

```math
\frac{\partial L_{\mathrm{direct}}}{\partial x_j}
= -\frac{L_{\mathrm{direct}}}{\mu_0}
  \frac{\partial \tau_{\mathrm{above}}}{\partial x_j}.
```

It is nonzero only for atmospheric parameter columns and solar-aligned output
ordinates; surface-parameter columns are zero. Use
`total_downwelling(level)` and `total_downwelling_jacobian(level)` when a
diffuse-plus-direct total is required. The BOA endpoint fields `boa` and
`boa_jacobian` retain the historical total-radiance convention and already
include this carrier; only strict-interior records keep it separate.

`ObserverRTResultLin` remains iterable as the historical
`(toa, boa, toa_jacobian, boa_jacobian)` tuple, so endpoint-only callers can
continue to write `R, T, dR, dT = rt_run_lin(...)`.

Linearized interior-height radiances are currently supported on the elastic
`noRS` path with `SolarBeam` or `NoSource`. Thermal and surface-emission
sources require linearized multisensor source-slot propagation and are
rejected explicitly, as is `CanopySurface`. Raman-active linearization remains
unsupported. The single-scatter driver `rt_run_ss` is still full-column-only
and rejects strict-interior requests explicitly; its endpoint selection
continues to follow the scalar/vector convention above.

## See also

- [`Schema/atmospheric_profile.md`](atmospheric_profile.md) — the input
  profile that defines BOA, TOA, and the vertical interpolation grid
- [`conventions.md`](../../conventions.md) — Q/U/V signs and azimuth
  conventions vs VLIDORT / Mishchenko / Hovenier
