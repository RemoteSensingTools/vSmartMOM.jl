# RRS–XCO2 Optimal-Estimation Inversions

This directory is the central guideline for the inversions built from the
`RRS_XCO2/truth_map` simulations. It should be updated whenever the state
definition, measurement processing, covariances, or retrieval assumptions
change.

The authoritative truth/forward/linearized consistency and retrieval-readiness
contract is
[`retrieval_setup/TRUTH_FORWARD_LINEARIZED_CONSISTENCY.md`](retrieval_setup/TRUTH_FORWARD_LINEARIZED_CONSISTENCY.md).
It must contain no unresolved gate before a retrieval campaign is started.

## Directory layout

- `corrected/`: class-specific products and run documentation using a
  Rayleigh-based measurement vector.
- `uncorrected/`: class-specific products and run documentation using an RRS
  + Cabannes measurement vector.
- `old/`: archived products, manifests, logs, plots, and prior files from
  superseded retrieval campaigns. Nothing under this directory is an active
  retrieval input.
- `instrument/`: the shared Stokes-analyzer, Gaussian-convolution, and
  synthetic-sampling operator used to construct both measurement classes.
- `retrieval_setup/`: the a priori state/covariance definition, fixed-state
  mask convention, compact retrieval-state layout, and documented OCO-2
  convergence and fit-quality criteria.

Shared Julia implementation lives at this directory's top level so both
classes execute the same solver and forward adapter. Python is reserved for
plotting.

## Retrieval method

Both classes will use optimal estimation (OE). For state vector `x`,
measurement vector `y`, forward model `F(x)`, measurement-error covariance
`S_e`, prior state `x_a`, and prior covariance `S_a`, the retrieval minimizes

```text
J(x) = (y - F(x))' S_e^-1 (y - F(x))
     + (x - x_a)' S_a^-1 (x - x_a).
```

The OCO-2 detector-noise contribution and completed a priori are defined
below. The first suite deliberately includes detector noise only; any future
forward-model, spectroscopic, or calibration-error terms must be documented
before being added.

### Forward model used by both retrieval classes

Both the corrected and uncorrected inversions use vSmartMOM's **linearized
elastic forward model**, with the Raman-scattering type set to `noRS` (the
source-code spelling is `RS_type::noRS`, instantiated as
`InelasticScattering.noRS{FT}()`). The analytic linearization supplies the
measurement Jacobian with respect to the complete retrieval state.

RRS is therefore not included in the retrieval forward model or its
Jacobians. The controlled distinction between the two inversion classes is in
the measurement vector:

- the corrected measurement is based on the truth-map Rayleigh simulation and
  is consistent with the `noRS` retrieval physics;
- the uncorrected measurement contains truth-map Cabannes + RRS radiances, but
  is fitted with the same `noRS` linearized forward model.

This common forward model is essential: differences between the two retrieval
classes can then be attributed to the uncorrected inelastic-scattering signal
in the measurements rather than to different retrieval implementations.

The shared OCO configuration selects componentwise Fourier convergence. At
each order, I, Q, and U must independently satisfy
`abs(delta_S) <= 1e-5*abs(S_partial)` at every wavelength and view. The
package default requests two successive passing moments; this validated
exact-nadir workflow explicitly requests one. The guard retains orders through
`min(2,m_max)`: this
prevents scalar Rayleigh's structural `beta_1=0` from hiding its nonzero
`beta_2`, while an exact lower-support problem is evaluated only through its
precomputed `m_max`. The combined
forward/analytic-linearization calculation
uses that single forward-field decision: all Jacobian columns are accumulated
through the accepted order and stop with the forward field. The package-wide
default remains the complete Fourier sum for workflows that have not validated
an early-exit criterion for their geometry.

Both classes use the retrieval-selective `OCO_RRS_synth` linearization:

```julia
model, lin_model = model_from_parameters(OCO_RRS_synth(), params)
result = rt_run_lin(model, lin_model; i_band=iband, sources=band_sources)
K_band = globalize_jacobian(result.toa_jacobian, result.layout)
```

This propagates only the active band-local physical columns through the MOM
kernels (12 in O2 A, 22 in each CO2 band) and maps them into one 30-column
cross-band basis. The aerosol columns at this stage are derivatives with
respect to native `tau_ref` and `z0`; the AOD760-reference conversion and
log-coordinate chain rule are part of retrieval-state assembly, as detailed
in `retrieval_setup/README.md`.

Selection occurs before Fourier-dependent phase mixing and before allocation
of the elemental/doubling/interaction tangent operators. It is therefore a
computational mask, not an output-only slice. The implementation and every
finite-difference validation are documented in
`docs/dev_notes/selective_jacobian_plans.md`.

The direct solar source is identical to the truth-map source: the 5777 K
Planck spectrum is multiplied by the high-resolution transmission spectrum in
`solar.out` before constructing `SolarBeam`. A smooth Planck-only source is
not a valid retrieval approximation because it removes the Fraunhofer lines.
The source path is recorded in every completed retrieval file.

## Common retrieved state

The corrected and uncorrected retrievals use the same physical state:

1. CO2 volume-mixing ratio in every atmospheric layer, `VMR_CO2[i]` for
   `i = 1:N_layers`.
2. Surface pressure, `p_surface`.
3. Three band-local Legendre surface coefficients in each OCO-2 band:
   `surface_P0[band]`, `surface_P1[band]`, and `surface_P2[band]`, for the
   O2 A-band, weak CO2 band, and strong CO2 band (nine coefficients total).
4. Two SIF parameters, `SIF760` and `mSIF`. SIF applies only to the O2 A-band;
   its contribution and derivatives in both CO2 bands are zero.
5. The natural logarithm of three aerosol optical depths, one each for the
   sulfate, organic-carbon, and UTLS/stratospheric aerosol components. The
   physical AODs use the documented 760 nm reference convention.
6. The natural logarithm of the three positive aerosol characteristic heights,
   one for each aerosol component.

For the current 16-layer atmosphere this gives a 34-element retrieval-facing
state (with aerosol AOD and height represented in log coordinates):

```text
16 CO2 VMRs + 1 surface pressure + 9 surface coefficients
+ 2 SIF parameters + 3 aerosol AODs + 3 aerosol heights = 34.
```

The current prior fixes the four CO2 layers centered above 10 km by assigning
zero variance, so the numerical OE solve has 30 active variables. In this
synthetic suite those fixed layers are set to each scene's uniform truth CO2
value before the solve; this keeps all four 380/400/420/440 ppm truth families
representable without adding retrieval coordinates. The full 34-element state
is retained for forward-model bookkeeping and output.

For SIF-on truth scenes, `0.5` denotes the version-2 reference-wavelength
normalization `2pi * L_lambda(760 nm) = 0.5 mW m^-2 nm^-1`, or
`L_lambda(760 nm) = 0.5/(2pi) mW m^-2 sr^-1 nm^-1` in each isotropic
upwelling BOA stream. It is not a wavelength-integrated SIF total. `SIF760`
and `mSIF` remain the native per-wavenumber retrieval coordinates.

The compact retrieval layout, units, fixed-layer mask, and its mapping to the
larger native vSmartMOM Jacobian are defined in `retrieval_setup/README.md`.
Atmospheric layers are listed TOA-to-BOA, matching vSmartMOM's internal
profile direction.

### CO2 a priori covariance

The active CO2 prior mean remains uniformly 400 ppm. Layers 1:4 are fixed as
described above. The covariance for active layers 5:16 is adapted from the
fixed 20-pressure-level CO2 covariance stored in all four OCO-2 aggregate
orbit files selected by `oco_gain.jl`. If `H` linearly maps the ACOS pressure
levels to the RRS-XCO2 layer centers in normalized pressure, the 16-layer
matrix is

```text
S_CO2,16 = H * S_CO2,ACOS * transpose(H).
```

The retrieval uses the marginal `5:16,5:16` block, retaining its off-diagonal
terms. It is not conditioned on the scene-truth values assigned to the fixed
upper layers: doing so would leak synthetic truth information into the active
prior. CO2 remains uncorrelated a priori with all non-CO2 state elements,
matching the inspected OCO-2 files.

At 1000 hPa the old diagonal prior implied `sigma(XCO2)=3.689 ppm`. The mapped
correlated covariance gives `sigma(XCO2)=13.716 ppm`; using the mapped
variances without their correlations would give only `5.074 ppm`. The exact
old and new covariance entries are stored in
`retrieval_setup/co2_prior_covariances.dat`. Their extraction, provenance,
mapping, conditioning decision, and orbit-by-orbit audit are documented in
`retrieval_setup/CO2_PRIOR_COVARIANCE_AUDIT.md`.

Aerosol AODs and profile heights are exponentiated before each forward-model
evaluation. Native physical Jacobian columns are converted by
`K_logq = K_q*q`; the corresponding prior and reported posterior covariances
use the local inverse/forward Jacobian transformations documented in the
retrieval-setup guide.

### Fixed aerosol properties

The following aerosol properties are not retrieved:

- complex refractive index/composition;
- particle size-distribution parameters;
- width of each aerosol vertical distribution.

Their values and spectral dependence remain those documented for the truth
map. Only each component's AOD and characteristic height vary in the
retrieval.

## Measurement-vector classes

### Corrected retrievals

The measurement vector is defined from the Rayleigh simulations in the truth
map. This represents the RRS-corrected case and is retrieved with the
linearized `noRS` forward model described above.

### Uncorrected retrievals

The measurement vector is defined from the truth-map RRS + Cabannes
simulations. Before retrieval it undergoes the shared analyzer projection,
Gaussian convolution, and synthetic OCO-2 resampling described in
`instrument/README.md`. The retrieval itself still uses the same linearized
`noRS` forward model; it does not model RRS explicitly.

For a controlled corrected-versus-uncorrected comparison, both measurement
classes should ultimately use the same finalized instrument convolution,
sampling grid, spectral masks, and noise convention; only the source
radiative-transfer component should differ. This requirement should be
confirmed when the measurement operator is finalized.

## Shared instrument processing

Both measurement classes use the same fixed operator, in this order:

1. apply a representative nadir OCO L1B analyzer row using the exact
   `oco_gain.jl` convention `M11*I - M12*Q + M13*U`, with no further
   normalization;
2. convert radiance density from per cm-1 to per nm;
3. convolve in wavelength with a Gaussian whose FWHM is 0.04, 0.08, and
   0.10 nm in the O2 A, weak CO2, and strong CO2 bands, respectively;
4. sample at `758:0.015:772`, `1594:0.031:1619`, and
   `2042:0.04:2082` nm.

The identical fixed operator is applied column-by-column to the linearized
`noRS` forward-model Jacobian, yielding a sampled `(measurement, parameter)`
matrix consistent with each measurement vector.

This is enforced by `validate_truth_forward_closure.jl`: for an aerosol-on,
SIF-off state, all 30 analytic columns are passed through the same analyzer,
spectral-density conversion, Gaussian convolution, and resampling routine and
compared with central perturbations of the high-resolution Stokes spectrum.
The 2026-08-29 ABSCO closure passed in all bands; worst absolute differences
were `5.15e-12`, `2.16e-12`, and `9.38e-13` for O2 A, weak CO2, and strong
CO2, respectively. The complete report is
`retrieval_setup/ABSCO_CLOSURE.md`.

The analyzer row is selected by a sign-preserving circular-median procedure
over all available nadir soundings. Gaussian integration includes the
nonuniform wavelength weights implied by the uniform 0.1 cm-1 truth grid.
Full +/-6 sigma source support is mandatory. A separate post-truth script
generated the eight missing short-wavelength strong-band points; these points
were validated and merged into the completed truth files. See
`instrument/README.md` for provenance, exact coefficients, commands, output
schema, and validation.

The appended strong-band shoulder does not redefine the surface state. The
three Legendre coefficients retain the truth map's 2042--2084 nm base-band
coordinate; an exact affine coefficient transform is used on the expanded
solve grid, and its three-by-three chain rule is applied to the analytic
surface Jacobian.

## Measurement-noise covariance

For every truth-map state, the instrument-noise covariance follows Eq. (3-8)
of the OCO-2 L1B ATBD. Each energy-radiance measurement is converted to photon
radiance, evaluated with the Table 3-5 MaxMS value and representative
wavelength-dependent `InstrumentHeader/snr_coef` arrays, and converted back to
energy-radiance NEN. Table 3-6 MinMS is used only to flag samples outside the
specified dynamic range.

The requested covariance assumes independent spectral-sample noise:

```text
S_e,instrument = Diagonal(NEN(lambda)^2).
```

It is signal dependent. Corrected and uncorrected O2 measurements therefore
receive separately evaluated covariance diagonals. Their weak- and strong-CO2
measurements are identical and have identical covariance entries. The compact
products are stored under `truth_map/OCO_radiances/noise_covariances/`; see its
README and `instrument/README.md` for the complete equation, coefficient
estimator, units, schema, and validation.

This diagonal covariance currently represents detector noise only. It does
not include spectral correlations, spectroscopy or forward-model error,
radiometric systematics, scene-inhomogeneity noise, or uncertainty caused by
using representative rather than sounding-specific instrument coefficients.

For each truth case and perturbation index, a standardized vector is drawn
from `Uniform(-sqrt(3),sqrt(3))`. It has zero mean and unit variance, so
`noise_std .* u` has covariance `S_e`. The corrected and uncorrected members
of a pair use the identical `u`, scaled by their separately stored noise
standard deviations. Distinct truth states—including matched aerosol and
no-aerosol scenes—use distinct deterministic seeds and therefore independent
standardized draws. Their `S_e` products are also evaluated independently from
their respective simulated radiances. `S_e` is frozen for the full retrieval
and is never recomputed from an iterated radiance. Every retrieval file stores
the noiseless measurement, perturbed measurement, normalized draw, noise
standard deviation, and exact physical `injected_measurement_noise` vector.

Perturbations 01:10 use these random draws. Perturbation 11 is the unperturbed
case: its normalized draw and injected-noise vector are identically zero, and
its retrieval measurement equals the noiseless truth measurement exactly.

## First retrieval suite

The current suite selects the 32 no-SIF truth states: four surfaces, four CO2
VMRs, and aerosol/no-aerosol cases. Ten noise perturbations plus one
unperturbed experiment and two measurement classes give 704 retrieval solves
(352 paired experiments). Both SIF
parameters remain active; only the truth SIF is zero. The linearized forward
model uses nine streams and `OCO_RRS_synth`.

Products are named `retrieval_stateNNN_perturbationPP.nc`, with `PP=01:11`,
under `corrected/` and `uncorrected/`. Each file records every accepted or
rejected LM trial, state, timing, per-band reduced chi-square, convergence
outcome, terminal Jacobian, gain, posterior covariance, and averaging kernel.
The shared `retrieval_manifest.dat` gives all 704 paths, paired seeds, and an
explicit `noise_injected` flag.

Perturbation 11 (the noiseless measurement) is computed first within each
truth state, followed by perturbations 01--10. Corrected/uncorrected members
remain adjacent, and this scheduling change does not alter perturbation
indices, canonical pair/retrieval IDs, filenames, or deterministic noise seeds.

Each completed product also records the dry-air-column-weighted CO2 VMR in
ppm. `XCO2` is the terminal value, `a_priori_XCO2` is the value at the prior,
and `XCO2_at_trial` follows every evaluated LM trial. The weighting is
`sum(VMR_CO2[z] * VCD_dry[z]) / sum(VCD_dry[z])`; dry-air columns are
recomputed at the trial surface pressure using the fixed humidity profile.
The four non-active upper layers remain part of the column and use the
scene-specific fixed CO2 value recorded in the file metadata.

## Corrected-versus-uncorrected gain analysis

`plot_corrected_vs_uncorrected_errors.py` compares matched retrievals in the
same 3-by-7 compact-state layout used throughout this study. Its default mode
uses filled markers for the terminal retrieval error relative to truth,

```text
x_retrieved - x_truth,
```

and open markers for the corresponding linearized retrieval error predicted
from the noiseless retrieval and that realization's terminal gain matrix,

```text
x_retrieved,0 - x_truth + G_terminal * Delta y_injected.
```

Here `x_retrieved,0` is the terminal state retrieved from perturbation 11, the
exactly noiseless measurement. Corrected and uncorrected retrievals use their
own class-specific noiseless terminal states. Each noisy realization still
uses its own terminal gain matrix rather than reusing the perturbation-11 gain.

The horizontal and vertical coordinates are the corrected and uncorrected
values, respectively. Marker color identifies the paired perturbation; index
11, the unperturbed measurement, is black and uses a square instead of the
circles used for perturbations 01:10. When the gain overlay is enabled, its
open square coincides exactly with its filled square because its injected
noise is zero. The companion whitespace table records the raw `G*Delta y`
increments and the recentered predictions numerically.

The gain product is formed in the native 30-coordinate state and then mapped
to the physical panel units. AOD and aerosol-height increments use the tangent
of the exponential transformation at the terminal state. SIF760 uses the
wavenumber-to-wavelength density conversion. The twelve CO2 increments are
condensed to XCO2 using the retrieval's dry-column weights, including the
first-order contribution from the surface-pressure increment. The script
checks that this independent XCO2 mapping reproduces the saved terminal XCO2
before accepting a retrieval.

The difference between an open and filled marker isolates departures from the
linearized noise propagation about the noiseless retrieval. Such departures
can arise because the terminal gain changes with the realization, because the
retrieval is nonlinear, or because convergence and prior regularization cause
the full retrieval path to differ from the local gain prediction. The common
noiseless offset retains any class-specific Raman model mismatch and other
systematic retrieval error already present without detector noise.

Usage:

```bash
# Truth error with the recentered linear gain/noise prediction (default)
python3 plot_corrected_vs_uncorrected_errors.py 49

# Same truth-error plot without the gain/noise overlay
python3 plot_corrected_vs_uncorrected_errors.py 49 --no-gain-noise

# Optional prior-referenced displacement; the truth-centered overlay is off
python3 plot_corrected_vs_uncorrected_errors.py 49 \
    --reference prior --no-gain-prediction
```

The `--gain-prediction` and `--no-gain-prediction` switches explicitly control
the overlay; `--gain-noise` and `--no-gain-noise` remain accepted aliases. The
overlay is truth-referenced by definition and therefore cannot be combined
with `--reference prior`. Both overlay modes use the same default PNG and
table names, so rerunning with the opposite switch replaces the prior
rendering. Use `--output` and `--table-output` when both renderings should be
retained side by side.

## Active bottom-layer campaign

The source/sink experiment documented in
[`retrieval_setup/BOTTOM_LAYER_XCO2_RETRIEVAL_PLAN.md`](retrieval_setup/BOTTOM_LAYER_XCO2_RETRIEVAL_PLAN.md)
is now active. Its five-state CO2 order is `360, 380, 400, 420, 440 ppm` in
the bottom layer, with a fixed 400 ppm background above it. The production
log-AOD prior is `sigma[ln(AOD760)] = 0.75`. All no-SIF retrievals must pass a
shared Curry/Wurst barrier before any SIF-on retrieval is released; see
`../bottom_layer_XCO2_retrievals/RETRIEVAL_RESTART_SIGMA_LNAOD_0P75.md` for
the frozen prior hash and restart record.

## Deferred extensions

- possible numerical bounds if the tightened log-aerosol prior still
  misbehaves;
- land-glint geometries after the nadir suite;
- any future correlated or forward-model contribution to `S_e`.

## Reproducibility rule

Every retrieval product should retain, directly or through a manifest, the
truth-map scene index, measurement class, input simulation path, wavelength
grid, instrument-processing configuration, state-vector layout, covariance
versions, vSmartMOM revision, numeric precision, and convergence status.
