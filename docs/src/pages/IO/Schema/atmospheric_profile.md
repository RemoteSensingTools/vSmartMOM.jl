# `atmospheric_profile` block

Vertical profile of the model atmosphere: temperature, pressure,
specific humidity, and (optionally) a level-reduction request. Profile
direction is **TOA → BOA** (top of atmosphere first); GEOSChem-style
data that arrives as BOA → TOA is flipped at parse time by
`geoschem_to_dict`, the reader behind `GeosChemSource`.

## Required fields

- **`T`** — `Vector{Real}`. Temperature [K] at full levels (layer
  centers). Length = number of layers.
- **`p`** — `Vector{Real}`. Pressure [hPa] at half levels (layer
  boundaries). Length = number of layers + 1. Must be monotonic in
  pressure (TOA = lowest p first).
- **`profile_reduction`** — `Union{Integer, Nothing}`. When an integer
  `N > 0`, model construction reframes the atmosphere onto at most `N` layers,
  interpolating intensive profiles and recomputing integrated layer fields.
  `null` (YAML) or `-1` keeps the full layering; TOML has no null, so use
  `-1` there. Observer interfaces requested
  through `geometry.obs_alt` are always retained. If `K` distinct interior
  observer heights require more interfaces than the requested reduction
  permits, the effective count is raised to `K + 1` layers and a warning
  reports the adjustment.

## Optional fields

- **`q`** — `Vector{Real}`. Specific humidity (water-vapor mass fraction,
  kg kg⁻¹). It may contain `N` layer-centred values (the preferred form) or
  `N+1` values on the pressure interfaces; interface values are interpolated
  to layer centers in log pressure. Default: zeros (dry atmosphere). Used by
  the absorption pipeline to derive the H₂O VMR automatically; do not also
  list H₂O in the absorption molecule lists. A legacy H₂O key in the `vmr`
  dictionary is ignored.

## Examples

### Single-layer Rayleigh scene (Natraj test fixtures)

```yaml
atmospheric_profile:
  T: [231.62]
  p: [0.14, 0.22]            # 1-layer column, ~0.08 hPa thick
  profile_reduction: null
```

### Multi-layer with profile reduction

```yaml
atmospheric_profile:
  T: [220.0, 230.0, 250.0, 270.0, 280.0, 290.0]
  p: [0.01, 1.0, 10.0, 100.0, 500.0, 800.0, 1013.25]
  q: [1e-6, 1e-5, 1e-4, 1e-3, 5e-3, 1e-2]
  profile_reduction: 5       # collapse 6 layers → 5
```

## Reading from external sources

For real GEOSChem outputs, use the high-level reader rather than YAML:

```julia
using vSmartMOM
src = vSmartMOM.GeosChemSource("/path/to/geoschem.nc4", 10, 20, 1) # x, y, face indices
config = vSmartMOM.geoschem_to_dict(src)
ap = config["atmospheric_profile"]
profile = vSmartMOM.read_atmos_profile_dict(Dict(
    "T" => ap["T"], "p_half" => ap["p"], "q" => ap["q"]))
```

`geoschem_to_dict` returns a complete configuration dictionary.
`read_atmos_profile_dict` converts a standalone profile mapping to an
`AtmosphericProfile`. That mapping names the pressure interfaces `p_half`
(or gives hybrid `ak`, `bk` and `p_surf`), whereas the scene block above
calls them `p`, hence the renaming. Use `read_parameters(src)` when building
a complete RT scene.

## See also

- [`Schema/absorption.md`](absorption.md) — VMR profiles bind to the
  same vertical grid as `T`
- [`Schema/geometry.md`](geometry.md) — `obs_alt` selects geometric heights
  in km above BOA; interior heights are inserted into this pressure grid
