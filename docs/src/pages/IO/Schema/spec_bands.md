# `spec_bands` field (in `radiative_transfer`)

The spectral grid(s) the RT solve runs on, in **wavenumber (cm⁻¹)**
internally. `spec_bands` is a `Vector{String}` — one string per band — and
each string is parsed **without `eval`** by `_safe_parse_spec_band` in
[`src/IO/Parameters.jl`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/src/IO/Parameters.jl).
Multiple entries run as a **band-batched solve**: the whole pipeline
(absorption, Mie, RT) loops over bands, one 3-D `(NquadN, NquadN, nSpec)`
array set per band.

## Field signature

```yaml
radiative_transfer:
  spec_bands:
    - "<band_expression>"   # one entry per band
```

`surface:` (and, for per-band absorption settings, `fixed_molecules` /
`variable_molecules`) must have the same length as `spec_bands`.

## Supported expression forms

### 1. Explicit list

```yaml
spec_bands:
  - "[19417.0 19418.0]"             # space-separated
  - "[19417.0, 19417.5, 19418.0]"   # comma-separated
```

A bracketed list of numbers, each a full-precision spectral point in cm⁻¹.
Used for tiny single-point or hand-picked-point scenes
(`config/quickstart.yaml`, `config/canopy_forest.yaml`).

### 2. Range: `start:step:stop`

```yaml
spec_bands:
  - "(1e7/780):0.1:(1e7/750)"    # wavenumber range from wavelength math
  - "12950.0:0.02:13050.0"       # plain wavenumber range
```

`start`, `step`, and `stop` are each evaluated by a small safe
recursive-descent arithmetic parser (`_safe_parse_number`) supporting
`+ - * /`, unary sign, and parentheses — e.g.
`"(1e7/775):0.05:(1e7/755)"` is a common way to express a wavelength window
in wavenumber terms, since ``\nu_{\mathrm{cm^{-1}}} = 10^7 / \lambda_{\mathrm{nm}}``.
No `eval` is used, so only numeric literals and these four operators are
accepted — arbitrary Julia expressions are rejected. A 2-part `start:stop`
(unit step) form is also accepted.

### 3. Wavelength input via `collect(...)u"nm"`

```yaml
spec_bands:
  - "collect(400.0:1.0:700.0)u\"nm\""     # range, given in nm
  - "collect([400.0, 700.0])u\"nm\""      # explicit list, given in nm
```

Wrapping either the range or the bracketed-list form in
`collect(...)u"nm"` interprets the enclosed values as **wavelengths in
nanometers**; they are converted to wavenumber via
``\nu = 10^7 / \lambda_{nm}`` and the result is sorted ascending
(wavelength-ascending inputs invert to wavenumber-descending, so the sort
keeps the internal grid monotonically increasing). `collect(...)u"cm_inv"`
is also accepted as a no-op tag (values already in the internal unit); a
bare `collect(...)` with no unit suffix is equivalent to omitting `collect`
entirely.

## Units and ordering

- Internal unit is always **cm⁻¹**; wavelength-tagged (`u"nm"`) input is
  converted at parse time, before the value ever reaches the solver.
- A raw wavenumber `start:step:stop` range should already be ascending —
  `collect` follows the sign of `step` and the parser does not re-sort a
  plain (non-`u"nm"`) range or list.

## Examples

### Single band, single point (quickstart)

```yaml
radiative_transfer:
  spec_bands:
    - "[12987.0]"
```

### O₂ A-band, wavelength-derived wavenumber range

```yaml
radiative_transfer:
  spec_bands:
    - "(1e7/780):0.1:(1e7/750)"
```

### Multiple bands (one BRDF and one absorption VMR list per band)

```yaml
radiative_transfer:
  spec_bands:
    - "[20000 20003]"
    - "[12500 12503]"
  surface:
    - "LambertianSurfaceScalar(0.05)"
    - "LambertianSurfaceScalar(0.40)"
```

See [`Schema/absorption.md`](absorption.md) for how `fixed_molecules` /
`variable_molecules` line up per band, and [`Schema/surface.md`](surface.md)
for the per-band `surface:` vocabulary.

## See also

- [`Schema/radiative_transfer.md`](radiative_transfer.md) — the rest of the
  `radiative_transfer` block this field lives in
- [`Schema/absorption.md`](absorption.md) — per-band molecule lists that
  must line up in length with `spec_bands`
- [`docs/src/pages/IO/ConfigurationGuide.md`](../ConfigurationGuide.md) —
  worked examples building up a scene from Level 0
