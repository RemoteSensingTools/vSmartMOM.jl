# `absorption` block (optional)

Gas absorption configuration: which molecules are active, where to
load HITRAN line data from, broadening / line-shape choices, and the
volume mixing ratios (VMRs) per molecule. Skip this block entirely
for clear-sky / Rayleigh-only / canopy-only runs.

## Required when present

- **`vmr`** — `Dict{String, Union{Real, Vector{Real}}}`. Volume mixing
  ratios per molecule, on the same vertical grid as
  `atmospheric_profile.T`. Scalar = uniform VMR; vector = profile.
  Example: `{ O2: 0.21, H2O: [1e-6, 1e-5, 1e-4, 1e-3] }`.
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
manually.

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
- `Schema/spec_bands.md` *(forthcoming)* — spectral grid the absorption
  cross-sections are evaluated on
