# Aerosol polarization in the synthetic OCO measurement

## Finding recorded 2026-08-27

The current truth-map experiment shows that aerosol-induced polarization has
a large effect on an OCO-like single-analyzer measurement in the O2 A-band.
This is a useful scientific result of the truth-map construction and should
be retained alongside the aerosol assumptions rather than treated as a
plotting detail.

The quantitative example below compares state 009 with state 001: urban
surface, SIF off, XCO2 = 380 ppm, SZA = 30 degrees, nadir viewing, and aerosol
AOD760 = 0.28 versus no aerosol. The representative OCO analyzer row is

```text
(M11, M12, M13) = (0.5, -0.470727354288, -0.168569713831).
```

After applying the plot-only division by `M11 = 0.5`, its response is

```text
I_OCO / M11 = I + 0.941454709 Q - 0.337139428 U.
```

At the current principal-plane geometry, `U = 0`. The difference between the
aerosol and clear elastic/noRS calculations therefore retains a substantial
`Q` contribution:

| Aerosol-minus-clear elastic/noRS diagnostic | Value |
|---|---:|
| Continuum mean of convolved Stokes-I difference | 0.0686634 |
| Continuum mean analyzer-Q contribution | -0.0126833 |
| Continuum mean normalized OCO response | 0.0559801 |
| Continuum reduction from polarization | 18.47% |
| Analyzer-minus-I relative RMS over the band | 23.49% |

Thus it is reasonable to describe the band-wide effect as approximately 25%
for this case, while retaining the more precise distinction that the selected
continuum reduction is approximately 18.5%. Cabannes alone shows the same
behavior: its analyzer-minus-I relative RMS is 23.73%, and its continuum
polarization reduction is 18.66%. The common Cabannes and Rayleigh/noRS
polarization responses cancel when those two simulations are differenced.

## Interpretation and polarimetric opportunity

The result demonstrates a strong aerosol-induced change in the measured
polarization state. An OCO-like instrument with enough independent
polarization measurements to separate `I`, `Q`, and possibly `U` could exploit
this sensitivity to detect aerosols and constrain their optical properties.
A single fixed analyzer already responds to the effect, but cannot by itself
separate polarization from intensity; a useful polarimetric design requires
multiple analyzer orientations, modulation states, or equivalent independent
polarization information.

The approximately 25% value is not universal. It applies to the current
geometry, surface, aerosol mixture, AOD, and representative analyzer row. Its
dependence on viewing geometry, relative azimuth, surface type, AOD, aerosol
height, size, and refractive index should be evaluated before using it as a
detection threshold. The present evidence is best described as strong
aerosol-driven depolarization for this controlled case.

## Why the pure RRS component behaves differently

The pure rotational-Raman component is already strongly depolarized, so the
same analyzer has very little additional effect on it:

| Pure RRS diagnostic | Clear | AOD760 = 0.28 |
|---|---:|---:|
| Continuum analyzer contribution relative to I | 0.798% | 0.784% |
| Analyzer-minus-I relative RMS over the band | 0.868% | 0.849% |

This sub-percent response is negligible compared with the approximately
23.5% elastic aerosol effect, which explains why the pure RRS curves show
almost no polarization-induced separation in the diagnostic plots.

The RRS intensity is not mathematically independent of aerosol radiative
transport: for state 009 minus state 001, its continuum intensity increases by
approximately 6.04% (and its band RMS change is approximately 6.03%). The
appropriate conclusion is therefore that pure RRS is nearly insensitive to
the **polarimetric analyzer effect** at the considered AOD, not that aerosol
scattering has exactly zero influence on its intensity.

## Measurement-vector convention

The `M11 = 0.5` division is for diagnostic plots only. Production variables in
`truth_map/OCO_radiances/OCO2sims_NNN.nc` retain the raw `oco_gain.jl`
analyzer response. The measurement vector and its corresponding measurement
covariance must use that same unnormalized convention.

Relevant products and scripts:

- `truth_map/state009_minus001_aerosol_difference_o2_components.png`;
- `truth_map/rayleigh_three_bands/state009_rayleigh_three_bands.png`;
- `scripts/plot_truth_state_o2_components.py`;
- `scripts/plot_aerosol_effect_o2_components.py`;
- `inversion/instrument/representative_stokes_coefficients.nc`.
