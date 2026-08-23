# `absorption` block (optional)

Gas absorption configuration: which molecules are active, where to
load HITRAN line data from, broadening / line-shape choices, and the
volume mixing ratios (VMRs) per molecule. Skip this block entirely
for clear-sky / Rayleigh-only / canopy-only runs.

## Required when present

- **`vmr`** — `Dict{String, Union{Real, Vector{Real}}}`. Dry-air molar
  ratios (mol species per mol dry air) for non-H₂O molecules. Scalar =
  uniform VMR; a vector may use the `N`
  layer centers or the `N+1` pressure interfaces. Interface-defined vectors
  are interpolated to layer centers in log pressure before observer-height
  reframing.
  Example: `{ O2: 0.21, CO2: [3.8e-4, 3.9e-4, 4.0e-4, 4.1e-4] }`.
  H₂O is derived from the atmospheric-profile `q` field; a legacy H₂O key in
  this dictionary is ignored.
- **`broadening`** — `String`. Line shape function. One of
  `"Voigt()"`, `"Doppler()"`, `"Lorentz()"`. Voigt is correct across
  the troposphere and stratosphere; Doppler/Lorentz are diagnostic
  modes for the high-altitude / surface limits respectively.
- **`CEF`** — `String`. Complex error function for Voigt evaluation.
  One of `"HumlicekWeidemann32SDErrorFunction()"`. Speed/accuracy
  trade-off rarely worth touching.
- **`wing_cutoff`** — `Integer`. Absorption-line wing cutoff in cm⁻¹
  (typically 25). Lines outside this distance from each spectral
  point are skipped.

### Molecule-list fields (one of)

- **`fixed_molecules`** — `Vector{Vector{String}}`, one inner list per
  spec_band. Molecules whose VMR is *not* a state variable in the
  linearization. Preferred form.
- **`variable_molecules`** — `Vector{Vector{String}}`, one inner list
  per spec_band. Molecules whose VMR *is* a state variable
  (Jacobians computed via the analytic `lin_*` path).
- **`molecules`** — *legacy*. Single molecule list applied to all
  bands; treated as fully-fixed. Don't combine with
  `variable_molecules` in new configs.

H₂O is auto-handled when `q` is set in the profile — don't list it
manually. Its internal `vmr_h2o` value is also a dry-air molar ratio,
`N_H2O/N_dry`, rather than a moist-air mole fraction.

## Optional CIA and continuum fields

- **`cia_files`** — list of HITRAN CIA inputs. An entry may be a legacy path
  string, or a mapping with:
  - **`path`** — required CIA file path.
  - **`reference_codes`** — optional reference-code string or nonempty list
    of strings, such as `"54"` or `["54"]`.
  - **`negative_policy`** — optional `error` (default) or `clamp_zero`.
- **`mtckd_file`** — optional AER MT_CKD H₂O-continuum NetCDF path.

A bare CIA path is accepted only when different HITRAN reference families do
not overlap on the requested model grid. Ambiguous overlap fails closed and
must be resolved with `reference_codes`. Negative cross sections likewise
fail by default; `clamp_zero` must be requested explicitly.

The surface-pressure tangent includes the hydrostatic-column and midpoint-
pressure response of CIA and MT_CKD at fixed temperature, humidity, and VMR.
Their abundance Jacobians are not yet implemented. Likewise, the analytic
line-absorption pressure tangent currently holds the pressure-dependent line
cross section fixed. Forward calculations at independently specified pressure
states recompute every opacity and are unaffected by these derivative limits.

```yaml
cia_files:
  - path: "/data/HITRAN_CIA/O2-O2_2024.cia"
    reference_codes: ["54"]
    negative_policy: error
  - "/data/HITRAN_CIA/O2-N2_2024.cia"
```

## ABSCO LUT files

`LUTfiles` is an array of file arrays, one per spectral band. Each inner array must contain the
non-H₂O absorbers in the same order as `fixed_molecules[band]` followed by
`variable_molecules[band]`. An optional H₂O file can occur anywhere in the inner array; it is
recognized from its ABSCO molecule ID and is driven by `atmospheric_profile.q`.

Three formats are accepted by extension:

- `.hdf`, `.h5`, `.hdf5`: original ABSCO tables, read by AtmosphericAbsorption. vSmartMOM infers
  O₂/WCO₂/SCO₂ from `spec_bands`, reads only that continuous slab, retains all H₂O foreign-broadener
  nodes, and places the Float32 cube on the configured CPU/GPU architecture.
- `.absco`: a portable table produced by `AtmosphericAbsorption.save_absco_lut`; it is restored on
  the configured architecture.
- all other extensions: the legacy vSmartMOM `InterpolationModel`/JLD2 loader.

For a three-band Float32 CUDA run, the relevant TOML fragment is:

```toml
[radiative_transfer]
spec_bands = ["12950:0.01:13200", "6170:0.01:6270", "4800:0.01:4900"]
float_type = "Float32"
architecture = "GPU()"
# ...the other radiative_transfer fields...

[absorption]
fixed_molecules = [["O2"], ["CO2"], ["CO2"]]
variable_molecules = [[], [], []]
LUTfiles = [
  ["${ENV:ABSCO_ROOT}/o2_v52.hdf",  "${ENV:ABSCO_ROOT}/h2o_v52.hdf"],
  ["${ENV:ABSCO_ROOT}/co2_v52.hdf", "${ENV:ABSCO_ROOT}/h2o_v52.hdf"],
  ["${ENV:ABSCO_ROOT}/co2_v52.hdf", "${ENV:ABSCO_ROOT}/h2o_v52.hdf"],
]
vmr = { O2 = 0.2095, CO2 = 0.00042 }
broadening = "Voigt()" # retained for HITRAN fallback; ABSCO interpolation is linear
CEF = "HumlicekWeidemann32SDErrorFunction()"
wing_cutoff = 25
```

Set `ENV["ABSCO_ROOT"]` before reading the file. The O₂/CO₂ absorber abundance still multiplies
optical depth separately. For every ABSCO gas LUT, `profile.vmr_h2o` is automatically passed to the
table's broadener axis; the H₂O absorber LUT is additionally multiplied by the profile's H₂O VMR.

## HITRAN edition

Set the HITRAN edition once per session via:

```julia
using vSmartMOM
vSmartMOM.set_hitran_edition!(2020)   # default; or 2016
```

This is a per-session preference, not a YAML field. The Artifacts
infrastructure caches HITRAN snapshots in `~/.julia/scratchspaces/`
and downloads on-demand from hitran.org if needed. See the
[HITRAN Data Management](../../Absorption/HITRAN_Data.md) page for the
dispatch flow and edition-selection workflow.

## Examples

### Single-band O₂ A-band retrieval

```yaml
absorption:
  vmr:
    O2: 0.21
    N2: 0.78
    Ar: 0.009
  variable_molecules:
    - ["O2"]                # one inner list per spec_band
  fixed_molecules:
    - ["N2", "Ar"]
  broadening: "Voigt()"
  CEF: "HumlicekWeidemann32SDErrorFunction()"
  wing_cutoff: 25
```

### Multi-band with H₂O profile

```yaml
atmospheric_profile:
  T: [220.0, 250.0, 280.0]
  p: [0.1, 100.0, 500.0, 1013.25]
  q: [1e-6, 1e-3, 5e-3]      # H₂O auto-included as variable
  profile_reduction: nothing

absorption:
  vmr:
    O2:  0.21
    N2:  0.78
    CO2: 4.0e-4
  variable_molecules:
    - ["O2"]
    - ["CO2"]
  fixed_molecules:
    - ["N2"]
    - ["N2", "O2"]            # cross-band different fixed sets allowed
  broadening: "Voigt()"
  CEF: "HumlicekWeidemann32SDErrorFunction()"
  wing_cutoff: 25
```

## See also

- [`Schema/atmospheric_profile.md`](atmospheric_profile.md) — `T`, `p`,
  `q` profile that VMRs bind to
- [`Schema/spec_bands.md`](spec_bands.md) — spectral grid the absorption
  cross-sections are evaluated on
