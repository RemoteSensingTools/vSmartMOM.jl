# Retrieval setup: a priori state and covariance

The planned near-surface source/sink experiment, its literature basis, truth
profiles, future surface-curvature prior, and campaign-product separation are
documented in
[`BOTTOM_LAYER_XCO2_RETRIEVAL_PLAN.md`](BOTTOM_LAYER_XCO2_RETRIEVAL_PLAN.md).
It is a future campaign and does not modify the active full-column retrievals.

The nonlinear stopping rule, bandwise spectral-fit classification, and their
OCO-2 provenance are documented separately in
[`OCO2_CONVERGENCE_CRITERIA.md`](OCO2_CONVERGENCE_CRITERIA.md).

The comparison of the archived diagonal CO2 prior with the vertically
correlated 20-level prior found in the four OCO-2 orbit files is recorded in
[`CO2_PRIOR_COVARIANCE_AUDIT.md`](CO2_PRIOR_COVARIANCE_AUDIT.md). Its mapped
covariance was approved for the replacement retrieval campaign. The exact old
and mapped matrices are tabulated in
[`co2_prior_covariances.dat`](co2_prior_covariances.dat).

Round 5 fixes both SIF coefficients and tightens only the UTLS sulfate AOD and
height priors. It deliberately keeps the complete round-4 CO2 prior unchanged;
the exact definition and comparison guard are in
[`ROUND5_FIXED_SIF_PRIOR.md`](ROUND5_FIXED_SIF_PRIOR.md). Round 6 fixes SIF but
restores the original UTLS uncertainties, with every non-SIF prior unchanged:
[`ROUND6_FIXED_SIF_PRIOR.md`](ROUND6_FIXED_SIF_PRIOR.md). The proposed CO2
vertical-covariance change remains a separate, unnamed future experiment.

The compatibility and validation audit for Christian's analytical-Jacobian
performance work is in
[`ANALYTICAL_JACOBIAN_SPEEDUP_INTEGRATION.md`](ANALYTICAL_JACOBIAN_SPEEDUP_INTEGRATION.md).

## Status

This directory records the a priori specification for the 16-layer RRS-XCO2
retrieval experiment. The physical state has 34 entries. Four upper-atmosphere
CO2 entries have exactly zero prescribed uncertainty and are fixed, leaving 30
active retrieval variables. For this synthetic experiment, their runtime values
are fixed to the uniform CO2 truth value of the current scene (380, 400, 420,
or 440 ppm). They are not unconditionally fixed to the 400 ppm active-state
prior, because that would make three of the four truth families impossible to
represent.

All requested quantities are numerically defined. The agreed small nonzero
SIF reference prior is `0.1 mW m^-2 sr^-1 nm^-1`; it also normalizes the
fractional slope. `build_apriori.jl` has generated the authoritative
`apriori_states.nc` and human-readable `apriori_states.dat` products.

## Statistical convention

The `+/-` values supplied for this experiment are one-standard-deviation
(`1 sigma`) uncertainties, not hard bounds or 95% intervals. The CO2 layers
use the mapped ACOS covariance, including its off-diagonal terms. The two
wavelength-space SIF coefficients also acquire covariance when transformed
into native wavenumber coefficients. Other parameter blocks are independent.

The full covariance is positive semidefinite but singular because four CO2
layers have zero variance. The OE implementation must not invert that full
matrix. It should remove the fixed entries and invert the resulting 30 by 30
active-state covariance, or enforce those entries with an explicit fixed-state
mask. A tiny replacement variance would change the stated prior and must not be
introduced silently.

## Compact state-vector layout

The retrieval-facing order is:

```text
 1       surface pressure
 2:17    CO2 VMR, layer 1 (TOA) through layer 16 (BOA)
18:20    ln(AOD760) for sulfate, organic carbon, and UTLS sulfate
21:23    ln(z0/km) for the three aerosol profile locations
24:26    O2 A-band surface P0, P1, P2
27:29    weak-CO2 surface P0, P1, P2
30:32    strong-CO2 surface P0, P1, P2
33:34    SIF760 and mSIF
```

This compact order is not the native full vSmartMOM Jacobian order. The native
linearized model carries seven columns per aerosol,
`[tau_ref, n_real, n_imag, r_median, sigma_g, z0, sigma0]`. Retrieval assembly
must select only columns 1 and 6 for each aerosol, hold the other five aerosol
properties fixed, and map all selected columns into the compact order above.

This selection is now implemented by the `OCO_RRS_synth` Jacobian flavour:

```julia
model, lin_model = model_from_parameters(OCO_RRS_synth(), params)
result = rt_run_lin(model, lin_model; i_band=iband, sources=band_sources)
K_band_global = globalize_jacobian(result.toa_jacobian, result.layout)
```

Wall-clock savings relative to the historical full physical Jacobian can be
reproduced with [`benchmark_reduced_jacobian.jl`](benchmark_reduced_jacobian.jl).
It reports warmed model-construction and per-band adding-doubling times
separately, and rejects a timing run unless every retained compact column
agrees with its full-Jacobian counterpart. The default uses 256 representative
0.1 cm^-1 samples per band; set `JACOBIAN_BENCH_NSPEC=0` for the complete
fine-resolution grids.

It propagates 12 local columns in the O2 A band and 22 local columns in each
CO2 band, then scatters them into the shared 30-column active state. Fixed
aerosol microphysics/H2O tangents are not computed. The native aerosol AOD
column is still `tau_ref` at the scattering configuration's reference
wavelength; conversion to the retrieval's AOD760 and subsequent logarithmic
chain rule must be applied by the retrieval-state mapping.

The exact local kernel layouts are:

```text
O2 A:       psurf, (tau_ref,z0)x3, surface(P0:P2), SIF760, mSIF
weak CO2:   psurf, (tau_ref,z0)x3, CO2(layer 5:16), surface(P0:P2)
strong CO2: psurf, (tau_ref,z0)x3, CO2(layer 5:16), surface(P0:P2)
```

The globalizer reorders these into the compact retrieval-facing order above
and inserts exact zeros for band-inactive parameters. It does not change
physical coordinates. At fixed aerosol microphysics,
`K_log(AOD760) = tau_ref*K_tau_ref` and `K_log(z0) = z0*K_z0`; CO2 columns
must be multiplied by `1e-6` when the OE state is expressed in ppm rather than
unit VMR.

All selected parameter classes have been checked against two-sided finite
differences. Maximum relative discrepancies were `8.8e-6` or smaller for
pressure, aerosol loading/height, selected layer absorption, all three
surface coefficients, and both SIF parameters. The complete allocation,
mapping, transformation, and validation record is in
`docs/dev_notes/selective_jacobian_plans.md`.

The physical aerosol height is `z0 = exp(mu)` in the model's
`LogNormal(log(z0), sigma0)` distribution. It is the input median, not the
height of maximum extinction. The compact retrieval coordinate is the natural
logarithm of this positive median. This distinction is essential when forming
and perturbing the height Jacobian.

## Atmospheric prior

Surface pressure is `1000 hPa` with `sigma = 50 hPa`.

The active CO2 prior is uniformly `400 ppm` (`400e-6` VMR). The fixed ACOS
20-pressure-level covariance is linearly mapped in normalized pressure to the
actual 16-layer centers. Layers 5:16 use the marginal mapped block, including
all off-diagonal entries, in internal TOA-to-BOA order:

| Layer | Center height (km) | `x_a` (ppm) | `sigma` (ppm) | Status |
|---:|---:|---:|---:|---|
| 1 | 40.792618 | 400 | 0 | fixed |
| 2 | 16.678894 | 400 | 0 | fixed |
| 3 | 13.178297 | 400 | 0 | fixed |
| 4 | 10.971210 | 400 | 0 | fixed |
| 5 | 9.325729 | 400 | 6.839 | active |
| 6 | 7.991297 | 400 | 7.795 | active |
| 7 | 6.854076 | 400 | 8.851 | active |
| 8 | 5.850269 | 400 | 9.954 | active |
| 9 | 4.944914 | 400 | 12.340 | active |
| 10 | 4.126884 | 400 | 14.877 | active |
| 11 | 3.370160 | 400 | 18.304 | active |
| 12 | 2.663361 | 400 | 22.375 | active |
| 13 | 2.015736 | 400 | 27.033 | active |
| 14 | 1.404734 | 400 | 32.809 | active |
| 15 | 0.820332 | 400 | 37.758 | active |
| 16 | 0.268311 | 400 | 43.527 | active |

The `400 ppm` entries shown for fixed layers 1:4 document the scene-independent
prior product. At runtime, `run_retrievals.jl` replaces those four non-active
values with `truth.xco2_ppm` before every solve. The output NetCDF records this
choice in `fixed_upper_co2_layers`, `fixed_upper_co2_ppm`, and
`fixed_upper_co2_source`. This is a controlled synthetic-closure convention;
an observational retrieval would instead require an external upper-atmosphere
CO2 prior.

At 1000 hPa this correlated active block gives a CO2-only
`sigma(XCO2)=13.715633 ppm`. The archived altitude-binned diagonal covariance
gave only `3.689129 ppm`. The covariance is deliberately marginal rather than
conditioned on the fixed scene-truth upper layers, because conditioning would
leak synthetic truth information into the active prior. The exact covariance
comparison and provenance are in `CO2_PRIOR_COVARIANCE_AUDIT.md`.

## Aerosol prior

All three physical AODs have `q_a = 0.02`. AOD means column optical depth at
760 nm. The retrieval coordinate is `u=ln(q)`. The production standard
deviation is `sigma_u=0.75`. This replaces the legacy direct choices
`sigma_u=5` and then `sigma_u=2`, both of which admitted overly large aerosol
steps during LM iteration:

```text
u_a = ln(0.02) = -3.912023005428
sigma_u = 0.75
variance(u) = 0.5625
```

The final tightening was tested with noiseless perturbation 11 for bottom-layer
states 001 and 013. State 001 is the urban, clear, no-SIF case with a 360 ppm
bottom CO2 layer; state 013 is the urban, aerosol, no-SIF 400 ppm control. This
pair deliberately exercises both the clear-scene boundary, where the retrieved
AODs should remain small, and the `AOD760=0.28` aerosol case that exposed the
large coupled trial steps. The comparison changes only the three independent
`ln(AOD760)` variances, leaving the prior means and every other covariance
entry fixed.

| Component | Physical `z0` (km) | `ln(z0/km)` prior | Extinction mode (km) | `sigma[ln(z0/km)]` |
|---|---:|---:|---:|---:|
| Sulfate | 1.525651538 | 0.422421556852 | 1.2 | 0.1 |
| Organic carbon | 2.112319568 | 0.747786665004 | 1.8 | 0.1 |
| UTLS sulfate | 12.120602005 | 2.494906649787 | 12.0 | 0.1 |

The 10% height uncertainty is applied to the physical `z0` median and converts
to `sigma[ln(z0/km)] = 0.1` at first order. The profile widths remain fixed at
0.49, 0.40, and 0.10, respectively.

For either AOD or height, retrieval iteration evaluates `q=exp(u)` before
updating the forward model. If `K_q = dF/dq` is the native physical Jacobian,
the compact-state column and covariance transformations are

```text
K_u = dF/du = K_q*q
S_u = D^-1*S_q*D^-1,       D = Diagonal(q_a)
S_q ~= D*S_u*D             (local physical-space reporting).
```

This guarantees positive forward-model AODs and heights without changing the
native vSmartMOM optical-property parameterization. `sigma_u=0.75` is a direct
log-coordinate prior choice, not a first-order transformation of a symmetric
physical-space AOD error bar. `AEROSOL_LN_AOD_SIGMA` remains available for
isolated sensitivity studies without changing the production default.

## Surface prior

For each scene and band, the surface prior is `(true P0, 0, 0)`. Its standard
deviations are `(0.10*P0, 1e-3, 1e-4)`. Consequently the prior is
surface-class-dependent, but is identical for all 16 truth states that share a
surface class.

Those values describe the full-column campaign. The bottom-layer CO2
campaign uses its own generated prior with `sigma(P1)=sigma(P2)=2e-3` in all
three bands; `P0` is unchanged. `build_apriori.jl` exposes common and
per-band environment controls so this does not alter the full-column product.

| Surface | Band | `P0` prior | `sigma(P0)` |
|---|---|---:|---:|
| urban | O2 A | 0.2715374708583 | 0.0271537470858 |
| urban | weak CO2 | 0.2522486169619 | 0.0252248616962 |
| urban | strong CO2 | 0.2176911056292 | 0.0217691105629 |
| rural | O2 A | 0.4316337681355 | 0.0431633768136 |
| rural | weak CO2 | 0.2649756999501 | 0.0264975699950 |
| rural | strong CO2 | 0.1335417267332 | 0.0133541726733 |
| desert | O2 A | 0.4186071552534 | 0.0418607155253 |
| desert | weak CO2 | 0.4972482231007 | 0.0497248223101 |
| desert | strong CO2 | 0.4821355300133 | 0.0482135530013 |
| forest | O2 A | 0.4682630263673 | 0.0468263026367 |
| forest | weak CO2 | 0.2557722421504 | 0.0255772242150 |
| forest | strong CO2 | 0.1100470188428 | 0.0110047018843 |

This is a truth-informed synthetic prior: only `P0` is initialized from the
truth, while the true nonzero `P1` and `P2` values are deliberately not used.
The source values are in `surface_albedos/lambertian_legendre_inputs.dat`.

## SIF parameterization

The requested wavelength-domain specification at 760 nm is

```text
L_lambda(760 nm) = 0.1 mW m^-2 sr^-1 nm^-1
sigma[L_lambda(760 nm)] = 0.25 mW m^-2 sr^-1 nm^-1
fractional slope g = -0.035 nm^-1
sigma(g) = 0.25*abs(g) = 0.00875 nm^-1.
```

The vSmartMOM linearized source instead retrieves the additive wavenumber
coefficients

```text
L_nu(nu) = SIF760 + mSIF*(nu - nu760),
```

where `SIF760` has units `mW m^-2 sr^-1 (cm^-1)^-1` and `mSIF` has units
`mW m^-2 sr^-1 (cm^-1)^-2`. Spectral densities obey

```text
L_nu = L_lambda * lambda_nm^2 / 1e7.
```

Thus the reference-radiance uncertainty converts unambiguously to

```text
sigma(SIF760) = 0.25 * 760^2 / 1e7 = 0.01444
                mW m^-2 sr^-1 (cm^-1)^-1.
```

The agreed reference scale gives

```text
dL_lambda/dlambda = -0.035*0.1 = -0.0035
sigma[dL_lambda/dlambda] = 0.00875*0.1 = 0.000875
```

The exact local coefficient conversion is

```text
a = L_lambda(760)
b = dL_lambda/dlambda = g*L_scale
SIF760 = a*lambda^2/1e7
mSIF = -(b*lambda^4 + 2*a*lambda^3)/1e14,  lambda = 760 nm.
```

`build_apriori.jl` treats `(a,b)` as independent before conversion. The linear transformation also
creates a nonzero covariance between `SIF760` and `mSIF`; this should be
retained rather than forcing the native SIF block to be diagonal.

For the no-SIF-first bottom-layer campaign, `b` is centered at zero while its
sigma is loosened to three times the historical value, `2.625e-3`, in
wavelength-space units. This removes the displacement of the no-SIF truth
created by the historical `b=-0.0035` mean and broadens the independent slope
dimension. The 0.1 reference-radiance mean and its 0.25 uncertainty remain
unchanged.

## Covariance construction

The retrieval covariance is block diagonal by parameter class. Its CO2 block
is the mapped, vertically correlated ACOS covariance; its native SIF block is
correlated through the exact wavelength-to-wavenumber coefficient transform.
All remaining blocks are diagonal in their stated retrieval coordinates. CO2
has zero a priori cross-covariance with every non-CO2 state element, matching
the four inspected OCO-2 files.

The machine-readable products are generated by `build_apriori.jl`, whose
default SIF reference is the agreed value `0.1`. It writes
`apriori_states.nc` and `apriori_states.dat` with:

- one 34-element prior for each of the four surface classes (or 64 scene
  records pointing to those four unique priors);
- the full 34 by 34 covariance and fixed-state mask;
- the reduced 30-element active state, covariance, and full-to-active index
  map;
- the CO2 covariance model, pressure nodes, mapping provenance, and implied
  prior XCO2 uncertainty;
- the exact compact-to-native vSmartMOM Jacobian column map.

## Corrected-versus-uncorrected reporting convention

All future corrected-versus-uncorrected truth-state comparisons use the same
compact terminal-state table.  The twelve retrieved layer CO2 VMRs are not
listed individually; they are replaced by the derived dry-air-column `XCO2`.
Every other retrieval quantity is retained in this order:

1. `XCO2` and surface pressure;
2. the three physical AOD760 values and, optionally, their derived sum;
3. the three physical aerosol profile-center heights;
4. `P0`, `P1`, and `P2` for each of the O2 A, weak-CO2, and strong-CO2 bands;
5. `SIF760` and `mSIF` in their native retrieval units.

The required columns are:

```text
parameter
truth
corrected mean +/- sample SD
uncorrected mean +/- sample SD
paired (uncorrected - corrected) mean +/- sample SD
```

Means and sample standard deviations are evaluated across the matched noise
perturbations, and the number of perturbations must be stated.  Corrected and
uncorrected members of each pair must use the identical normalized noise
draw.  These ensemble standard deviations must not be described as posterior
uncertainties.  Convergence counts and per-band reduced chi-squared statistics
are reported separately from the state table.

For an aerosol-free truth scene, the table still lists the configured nominal
profile-center heights, but marks them as spectrally non-identifiable because
their corresponding true AODs are zero.  Comparison plots follow the same
compact parameter ordering, show truth and prior reference lines, distinguish
corrected from uncorrected members, and show the ensemble mean +/- one sample
standard deviation.
