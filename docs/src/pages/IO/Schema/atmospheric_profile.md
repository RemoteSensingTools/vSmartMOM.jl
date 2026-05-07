# `atmospheric_profile` block

Vertical profile of the model atmosphere: temperature, pressure,
specific humidity, and (optionally) a level-reduction request. Profile
direction is **TOA → BOA** (top of atmosphere first); GEOSChem-style
data that arrives as BOA → TOA is flipped at parse time inside the
`Sources.GeosChemSource` reader.

## Required fields

- **`T`** — `Vector{Real}`. Temperature [K] at full levels (layer
  centers). Length = number of layers.
- **`p`** — `Vector{Real}`. Pressure [hPa] at half levels (layer
  boundaries). Length = number of layers + 1. Must be monotonic in
  pressure (TOA = lowest p first).
- **`profile_reduction`** — `Union{Integer, Nothing}`. When an integer
  `N`, the parser collapses adjacent layers down to ≤ `N` layers
  using mass-weighted averaging (preserves total optical depth, makes
  the solver faster). `nothing` keeps the full layering.

## Optional fields

- **`q`** — `Vector{Real}`. Specific humidity (water vapor mass mixing
  ratio) at the same full levels as `T`. Default: zeros (dry
  atmosphere). Used by the absorption pipeline when H₂O is in the
  active molecule list.

## Examples

### Single-layer Rayleigh scene (Natraj test fixtures)

```yaml
atmospheric_profile:
  T: [231.62]
  p: [0.14, 0.22]            # 1-layer column, ~0.08 hPa thick
  profile_reduction: nothing
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
src = vSmartMOM.GeosChemSource("/path/to/geoschem.nc4")
profile = vSmartMOM.IO.read_atmosphere(src; level_lon=120, level_lat=30)
```

The reader produces an `AtmosphericProfile` directly; bypass the YAML
`atmospheric_profile` block when using it.

## See also

- [`Schema/absorption.md`](absorption.md) — VMR profiles bind to the
  same vertical grid as `T`
- [`Schema/geometry.md`](geometry.md) — `obs_alt` is interpreted in Pa
  against this pressure grid
