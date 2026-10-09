# Future bottom-layer XCO2 retrieval campaign

Status: truth-production scripts prepared; no bottom-layer simulations have
been launched. Do not launch a partition on a host until its current
full-column retrieval worker has stopped. Do not move active full-column
products until all current jobs have completed and been archived.

## Purpose

The current campaign changes CO2 uniformly through the atmospheric column.
The future campaign will instead keep a 400 ppm background in every layer and
change CO2 only in the bottom model layer. This creates an idealized test of
how corrected and uncorrected OCO-like retrievals respond to a near-surface
source or sink.

This construction represents a horizontally homogeneous, instantaneous CO2
burden in the lowest layer. It is not itself an emission flux. Converting a
surface flux into a mole-fraction profile would additionally require the
boundary-layer depth, mixing time, wind field, horizontal footprint, and
atmospheric transport.

## Truth profiles

The current atmosphere has 16 layers in TOA-to-BOA order. Retrieval layer 16
is the bottom layer. At 1000 hPa surface pressure it spans approximately
937.509--1000 hPa and 0--0.541 km and contains

```text
w_bottom = 0.06229780736
```

of the dry-air column. With all other layers held at 400 ppm,

```text
Delta XCO2 = w_bottom * (CO2_bottom - 400 ppm).
```

The agreed truth grid is:

| Case | Layers 1:15 | Layer 16 | Resulting XCO2 | Delta XCO2 |
|---|---:|---:|---:|---:|
| stronger sink | 400 ppm | 360 ppm | 397.508088 ppm | -2.491912 ppm |
| weaker sink | 400 ppm | 380 ppm | 398.754044 ppm | -1.245956 ppm |
| control | 400 ppm | 400 ppm | 400.000000 ppm | 0 ppm |
| weaker source | 400 ppm | 420 ppm | 401.245956 ppm | +1.245956 ppm |
| stronger source | 400 ppm | 440 ppm | 402.491912 ppm | +2.491912 ppm |

The 400 ppm control is required even though it is not one of the four
perturbed cases. It supplies the matched reference for spectral differences,
noise sensitivity, and retrieval bias.

If the atmospheric pressure grid or surface pressure changes, these XCO2
values must be recomputed from dry-air column weights; the factor above must
not be copied to a different grid.

## Truth-map indexing and component reuse

The 80 states retain the original nesting order
`surface -> aerosol -> SIF -> CO2`, with the five bottom-layer CO2 values as
the fastest-changing index. The exact index is

```text
20*(surface_index - 1) + 10*(aerosol_index - 1)
    + 5*(sif_index - 1) + bottom_co2_index.
```

The accepted 400 ppm full-column state is physically identical to the new
400 ppm control for every surface/aerosol/SIF combination, but its index
changes because the innermost dimension now has five values. The complete
mapping is generated in `truth/control_reuse_map.dat`.

No new O2 A-band solve is required: CO2 is absent from the A-band absorber
list, so all three O2 arrays are reused from the matching accepted 400 ppm
scene. Weak and strong CO2 contain no SIF, so their newly computed no-SIF
arrays are reused in the paired SIF-on scenes. Only the four non-control
bottom-layer profiles are solved for each surface/aerosol combination. This
reduces each machine's RT work to 16 profiles times two CO2 bands.

The dedicated producer, launcher, exact machine partition, atomic publication
rules, and recovery precautions are documented in
[`../../bottom_layer_XCO2_retrievals/README.md`](../../bottom_layer_XCO2_retrievals/README.md).

## Literature basis

The state spacing is based on anomalies relative to the local atmospheric
background, rather than on historical absolute CO2 values. This distinction
matters because the background concentration differs among measurement years.

Broad lower-atmospheric observations support natural source and sink
perturbations of roughly 10--40 ppm:

- Photosynthetic uptake over Iowa produced approximately 15--17 ppm depletion
  through an approximately 2 km planetary boundary layer
  ([Ramanathan et al., 2015](https://doi.org/10.1002/2014GL062749);
  [Menzies et al., 2014](https://doi.org/10.1175/JTECH-D-13-00128.1)).
- Amazon ecosystem respiration produced early-morning enhancements of roughly
  10--20 ppm averaged through about 750 m
  ([Lloyd et al., 2007](https://doi.org/10.5194/bg-4-759-2007)).
- Siberian forest and bog measurements found non-combustion enhancements of
  about 15--22 ppm through the lowest 300 m, with more localized diurnal
  gradients approaching 40 ppm
  ([Kozlova et al., 2008](https://doi.org/10.1029/2008GB003209)).
- A West Siberian taiga record found a 29 ppm seasonal amplitude in the
  planetary boundary layer
  ([Sasakawa et al., 2013](https://doi.org/10.1002/jgrd.50755)).

Much larger concentrations have been measured within forest canopies and very
shallow stable nocturnal layers. For example, early-morning values reached
about 480 ppm over forest and 540 ppm over a deforested Amazon site
([Acevedo et al., 2008](https://doi.org/10.1029/2007JG000612)). Those local,
mostly nocturnal maxima do not justify filling the complete 541 m model layer
with the same concentration, particularly for a reflected-sunlight
observation. They are therefore excluded from the primary grid.

The implied column anomalies also occupy a useful OCO-like detection range:

- Monthly modeled biospheric signals had global mean absolute values of about
  0.5 ppm in February and 1.3 ppm in July; the July 90th percentile was
  approximately 3.8 ppm
  ([Miller et al., 2018](https://doi.org/10.5194/acp-18-6785-2018)).
- Bias-corrected OCO-2 v11.1 comparisons with TCCON retained approximately
  0.72--0.85 ppm scatter, depending on observing mode
  ([Das et al., 2025](https://doi.org/10.1029/2024EA003935)).
- OCO-2 detected natural volcanic enhancements of approximately 1--2 ppm at
  Kilauea and about 3.4 ppm at Yasur
  ([Johnson et al., 2020](https://doi.org/10.1029/2020GL090507);
  [Schwandner et al., 2017](https://doi.org/10.1126/science.aam5782)).

Consequently, the +/-1.246 ppm cases probe the practical detection boundary,
whereas the +/-2.492 ppm cases should be clearly detectable under favorable
conditions. Detection of the column perturbation does not guarantee that the
retrieval can localize it to layer 16: the retrieved change may spread among
correlated lower-tropospheric CO2 layers.

## Retrieval prior

The future campaign will retain the current 400 ppm CO2 prior mean and the
current mapped ACOS layer-CO2 prior covariance, including its off-diagonal
terms. The covariance provenance and mapping are documented in
[`CO2_PRIOR_COVARIANCE_AUDIT.md`](CO2_PRIOR_COVARIANCE_AUDIT.md). No special
constraint based on the known bottom-layer truth may be introduced.

The four fixed upper-atmosphere CO2 layers remain at their truth value, which
is 400 ppm for every profile in this campaign. The active lower-layer prior
means also start at 400 ppm.

Before launching the complete ensemble, the noiseless control and four
perturbed cases must be used to check:

1. the measurement-space signal relative to the frozen measurement covariance;
2. the retrieved XCO2 and layer-16 averaging-kernel response;
3. whether the signal is detected but vertically redistributed by the prior;
4. closure of the forward and analytical-Jacobian convolution/resampling paths.

The general a priori state definition remains in [`README.md`](README.md), and
the required readiness gates remain in
[`TRUTH_FORWARD_LINEARIZED_CONSISTENCY.md`](TRUTH_FORWARD_LINEARIZED_CONSISTENCY.md).

### Round-5/round-6 separation

Round 5 retains this mapped tapered CO2 mean and covariance without any
modification, allowing a direct comparison with round 4. Its only prior/state
changes are fixing both SIF coefficients and reducing the UTLS sulfate
`ln(AOD760)` and `ln(z0/km)` standard deviations by a factor of ten. See
[`ROUND5_FIXED_SIF_PRIOR.md`](ROUND5_FIXED_SIF_PRIOR.md).

Round 6 instead fixes both SIF coefficients while restoring the original
UTLS uncertainties; all non-SIF means and covariances remain unchanged. See
[`ROUND6_FIXED_SIF_PRIOR.md`](ROUND6_FIXED_SIF_PRIOR.md). An altitude-dependent
CO2 coupling relaxation remains an unnamed future experiment and must be
validated separately; it is not included in either round 5 or round 6.

## Surface-curvature prior for the future campaign

The present prior has `sigma(P2) = 1e-4` in every band. The weak-CO2 desert
truth has

```text
P2_desert_weak = 5.177341280743e-3,
```

which is 51.8 current prior standard deviations from the zero-curvature prior
mean. This can force a continuum-shape mismatch into CO2 or other retrieved
parameters.

For the bottom-layer campaign, use the provisional slope and curvature
settings in all three bands

```text
sigma(P1_O2_A) = sigma(P1_weak_CO2) = sigma(P1_strong_CO2) = 2.0e-3,
sigma(P2_O2_A) = sigma(P2_weak_CO2) = sigma(P2_strong_CO2) = 2.0e-3.
```

This places the desert weak-band truth at 2.59 standard deviations while
avoiding either a desert-specific or band-selective relaxation. Keep the `P0`
prior definition unchanged. Apply the same `P1` and `P2` priors to corrected
and uncorrected retrievals and freeze them across all noise perturbations.
Confirm the setting with the noiseless desert 400 ppm control (state 043)
before starting the full ensemble. Its corrected retrieval must converge,
pass the per-band fit-quality test, and finish with OE outcome 1. The paired
uncorrected retrieval contains the intentional RRS/noRS model mismatch, so it
is a structural/finite-result check rather than a scientific-closure gate.

The no-SIF-first bottom-layer ensemble also centers the absolute
wavelength-space SIF-slope prior at zero. Its one-sigma width remains
three times the former value, `2.625e-3 mW m-2 sr-1 nm-2`. The earlier mean
of `-3.5e-3` placed a zero slope four prior standard deviations away under
the old width and would have imposed an avoidable nonzero-SIF pull on every
truth case. The SIF reference-radiance prior remains
`0.1 +/- 0.25 mW m-2 sr-1 nm-1`.

The truth coefficients and their ECOSTRESS provenance are recorded in
[`../../surface_albedos/lambertian_legendre_inputs.dat`](../../surface_albedos/lambertian_legendre_inputs.dat)
and [`../../surface_albedos/PROVENANCE.md`](../../surface_albedos/PROVENANCE.md).

This prior change belongs only to the future bottom-layer campaign. It must
not alter or retrospectively reinterpret the active full-column campaign.

## Campaign-product separation

After the active computations finish and pass validation, archive campaign
products under separate top-level roots beneath `RRS_XCO2`:

```text
RRS_XCO2/
  full_column_XCO2_retrievals/
    truth/
    retrievals/
      corrected/
      uncorrected/
    manifests/
    logs/

  bottom_layer_XCO2_retrievals/
    truth/
    retrievals/
      corrected/
      uncorrected/
    manifests/
    logs/
```

Shared source code remains in `RRS_XCO2/inversion`; it should not be copied
into each product tree. Every manifest must identify the campaign, truth
profile definition, prior/covariance version, source revision, instrument
processing configuration, and random-noise seed.

Do not move files while active workers still use the present `truth_map`,
`corrected`, or `uncorrected` paths. At campaign completion, first stop or
verify completion of all workers, validate product counts and files, then move
the products and update every manifest and script path in one controlled
transition. The current worker state is tracked in
[`../ACTIVE_CAMPAIGN_STATUS.md`](../ACTIVE_CAMPAIGN_STATUS.md).
