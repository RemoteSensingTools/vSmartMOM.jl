# Aerosols for the OCO truth map

## Purpose and scope

This directory documents the vertical aerosol basis chosen for the 64-scene
OCO experiment. The experiment targets tropospheric and upper-troposphere /
lower-stratosphere (UTLS) aerosol sensitivity. It does **not** attempt to model
the optically much thinner, quiescent background Junge layer.

The profiles are normalized altitude-lognormal distributions. For species
`a`, its 760 nm column optical depth is distributed as

```text
p_a(z) = exp[-(ln(z)-ln(z_med,a))^2/(2 sigma_a^2)]
         / (sqrt(2*pi) sigma_a z)

Delta tau(a,k) = tau760(a) * integral[p_a(z) dz, z_bottom(k), z_top(k)]
                 / integral[p_a(z) dz, 0, z_model_top].
```

The layer integrals are evaluated as exact lognormal-CDF differences, rather
than by evaluating the profile at layer centers. A single column normalization
at the model top ensures that the layer AODs sum exactly to the prescribed
column AOD; it removes only the tiny mathematical tail above the atmosphere.

## Why altitude-lognormal profiles?

Sanghavi et al. (2012) used this same normalized functional family for aerosol
height retrievals from the O2 A and B bands. A broad mode close to the surface
can approximate an exponentially decreasing boundary-layer profile, whereas a
finite-altitude mode can describe aerosol introduced or produced aloft and
subsequently mixed in an atmosphere of decreasing density.

Hollstein and Filipitsch (2014) tested sums of these modes against the LIVAS
climatology, which is based primarily on CALIPSO/CALIOP 532 nm extinction
profiles. One lognormal mode reconstructed 59% of their realistic global
profile proxy, two modes 74%, and three modes 79%. This supports lognormal
modes as a compact empirical basis, while also showing that they cannot encode
all multilayer and irregular profile structure. Their reconstruction accuracy
varied strongly by region and more weakly by season.

The OMI aerosol ATBD used aerosol-type-dependent exponential profiles for
several near-surface types and elevated Gaussian layers for smoke and dust.
Those assumptions are qualitatively representable within this lognormal
family: a broad low mode approximates exponential decay, and a narrower
finite-height mode represents an elevated transported layer.

This differs from the ACOS/OCO Gaussian-in-relative-pressure parameterization.
That parameterization remains useful as a constrained retrieval state, but it
is not imposed on this controlled truth experiment.

## Selected modes

| Species | AOD at 760 nm | True mode (km) | `sigma_ln_z` | Median (km) | Interpretation |
|---|---:|---:|---:|---:|---|
| Tropospheric sulfate | 0.1935471100 | 1.2 | 0.49 | 1.525 | Lower-tropospheric mode |
| Organic carbon | 0.0807084200 | 1.8 | 0.40 | 2.112 | Somewhat more elevated pollution/transport mode |
| UTLS sulfate | 0.0057444777 | 12.0 | 0.10 | 12.121 | Narrow deliberately enhanced UTLS layer |

The optical depths preserve the current 0.28 total at 760 nm. The UTLS share
is therefore an experiment-specific enhanced layer, not an estimate of the
background Junge-layer AOD.

## Composition and particle-size distributions

All three components are represented by homogeneous spherical particles and
Mie scattering. Their complex refractive indices are fixed scenario values,
`m = n_real + i*n_imag`; they are not functions of wavelength, humidity,
temperature, or chemical aging in the present truth calculation.

| Species | Composition represented | `n_real` | `n_imag` | Median radius (micron) | Geometric sigma |
|---|---|---:|---:|---:|---:|
| Tropospheric sulfate | Weakly absorbing sulfate accumulation-mode proxy | 1.43 | 1.0e-8 | 0.10 | 1.60 |
| Organic carbon | Moderately absorbing organic-carbon accumulation-mode proxy | 1.53 | 1.0e-2 | 0.20 | 1.80 |
| UTLS sulfate | Weakly absorbing sulfate-droplet proxy | 1.45 | 1.0e-8 | 0.15 | 1.50 |

For each species, particle radius follows the normalized lognormal number
distribution

```text
n(r) = exp[-(ln(r)-ln(r_med))^2/(2 ln(sigma_g)^2)]
       / (sqrt(2*pi) ln(sigma_g) r).
```

Here `r_med` is the median particle radius and `sigma_g` is the dimensionless
geometric standard deviation. These particle-size parameters are independent
of the altitude-lognormal parameters above. Mie integration extends to a
particle radius of 10 micron using 300 radius-quadrature points. Optical
properties are referenced at 0.55 micron and spectrally evaluated/interpolated
by the truth-map forward model.

The labels describe optical proxies rather than complete chemical mixtures.
In particular, real atmospheric sulfate is normally an aqueous sulfuric-acid
solution whose refractive index and size respond to composition and ambient
conditions; organic carbon spans a wide range of absorptivity. The fixed
indices deliberately isolate vertical-distribution and radiative-transfer
effects, but should not be interpreted as a humidity-resolved aerosol model.

## Fourier-support evaluation

The untruncated Mie allocation at the shortest O2 A-band wavelength
(757.001655 nm) previously implied `m_max = 220` from the maximum particle
size parameter. That ceiling is conservative because the lognormal radius
distributions place very little integrated weight in their largest-radius
tails.

The truth-map configuration therefore sets:

```yaml
radiative_transfer:
  greek_beta_cutoff: 1.0e-5
```

After integrating each complete size distribution into its Greek
coefficients, the model finds the last degree satisfying
`abs(β_l) >= 1e-5`. The evaluation is beta-only: α, γ, δ, ϵ, and ζ are not
used because they can legitimately be non-positive and are not the scalar
phase-function support benchmark. The resulting supports at this wavelength
are:

| Species | Last retained degree |
|---|---:|
| Tropospheric sulfate | 24 |
| Organic carbon | 99 |
| UTLS sulfate | 22 |

The maximum across species gives the aerosol effective-support candidate
`m_max = 99`, instead of 220. With at least 50 streams and no tighter explicit
cap, the benchmark therefore evaluates 100 Fourier moments (`m = 0:99`)
instead of 221 (`m = 0:220`). For a spectral band, a degree is retained if it
passes at any wavelength, after which the maximum across all active components
and the configured stream/user caps determine the final band value. Thus the
production truth-map configuration's current `nstreams: 5` still imposes its
smaller `stream_l_cap = 9`; the value 99 describes the untruncated-resolution
investigation, not that capped production run.

This cutoff changes the Fourier-loop length only. It does not automatically
reduce the configured quadrature streams; resolving order 99 requires at
least 50 streams under `stream_l_cap = 2*nstreams - 1`.

### Median versus mode

For the parameterization above,

```text
z_mode = z_median * exp(-sigma_ln_z^2).
```

The input accepted by `Distributions.LogNormal` is the median parameter, not
the height at which extinction is maximal. The table deliberately specifies a
physical modal height first and calculates
`z_median = z_mode * exp(sigma_ln_z^2)`. This prevents the earlier ambiguity in
which `z0=12 km` was described as a peak although it was used as a median.

The lower-mode numerical values are controlled scenario choices informed by
the functional behavior in Sanghavi et al.; Hollstein and Filipitsch support
the model family and number of modes, not these three universal parameters.
The tropospheric widths were additionally constrained so that the continuous
extinction density at 10 km is no more than `1e-4` of its value at the mode.
Solving this condition gives maximum widths of approximately 0.494 for the
1.2 km sulfate mode and 0.400 for the 1.8 km organic-carbon mode; the adopted
values are 0.49 and 0.40, respectively.

## Files and reproduction

- `aerosol_vertical_profiles.dat`: exact boundaries, parameters, CDF layer
  fractions, and layer AOD for both the 12- and 16-layer grids.
- `continuous_vertical_profiles.png`: continuous extinction/AOD-density
  profiles, with their modes and medians distinguished.
- `layer_integrated_profiles_12_vs_16.png`: exact layer AOD alongside model
  boundaries for both vertical grids.

Regenerate with:

```bash
julia --project=. RRS_XCO2/scripts/build_truth_map_aerosols.jl
python3 RRS_XCO2/scripts/plot_truth_map_aerosols.py
```

## References and observational limitations

1. Sanghavi, S., Martonchik, J. V., Landgraf, J., and Platt, U. (2012),
   *Retrieval of the optical depth and vertical distribution of particulate
   scatterers in the atmosphere using O2 A- and B-band SCIAMACHY observations
   over Kanpur: a case study*, Atmospheric Measurement Techniques 5,
   1099–1119. <https://doi.org/10.5194/amt-5-1099-2012>
2. Hollstein, A. and Filipitsch, F. (2014), *Global representation of aerosol
   vertical profiles by sums of lognormal modes: Consequences for the passive
   remote sensing of aerosol heights*, Journal of Geophysical Research:
   Atmospheres 119, 8899–8907. <https://doi.org/10.1002/2014JD021472>
3. OMI Algorithm Theoretical Basis Document, Volume III, aerosol vertical
   distribution assumptions. <https://omi.fmi.fi/docs/OMI_ATBD_Volume_3_V2.pdf>

The LIVAS/CALIPSO climatology is not error-free: optically thin layers can be
missed, extinction depends on an assumed lidar ratio, and climatological
averaging smooths individual multilayer scenes. The plots here must therefore
be understood as reproducible basis profiles, not unique climatological truth.
