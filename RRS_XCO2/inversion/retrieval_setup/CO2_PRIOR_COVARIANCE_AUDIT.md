# OCO-2 CO2 a priori covariance audit

Status: **approved and implemented on 2026-08-31**. The diagonal-prior
retrieval products were moved to `../old/`, the generated RRS-XCO2 prior now
uses the mapped covariance, and replacement no-aerosol retrievals were
launched/queued on Curry and Wurst. Aerosol retrievals remain gated on their
final OCO-radiance and noise-covariance inputs.

## Inputs inspected

The four aggregate orbit files selected by
`/home/sanghavi/code/github/OCORaman/OCORaman/src/OCOPlots/oco_gain.jl` were
inspected directly:

1. `Nadir_37214_TallVeg.h5` (6,842 soundings)
2. `Nadir_42610_Himalayas.h5` (4,698 soundings)
3. `Nadir_41856_Sahara.h5` (13,377 soundings)
4. `Nadir_45442_Andes.h5` (7,910 soundings)

The relevant datasets are:

- `RetrievedStateVector/state_vector_names`
- `RetrievedStateVector/state_vector_apriori`
- `RetrievalResults/apriori_covariance_matrix`
- `RetrievalResults/vector_pressure_levels_apriori`
- `RetrievalResults/retrieved_dry_air_column_layer_thickness`

The first 20 state-vector elements are `CO2 VMR for Press Lvl 1` through
`CO2 VMR for Press Lvl 20`, ordered from TOA to the surface. The pressure
levels are essentially equally spaced in normalized pressure.

## What the orbit files establish

The 20 by 20 CO2 block of `apriori_covariance_matrix` is bit-for-bit identical
for every sounding in all four files. The maximum absolute difference from the
first TallVeg sounding is zero. CO2 covariance with all 44 non-CO2 state
elements is also exactly zero in every sounding. Thus ACOS uses:

- a sounding-dependent CO2 prior mean profile; and
- one fixed, vertically correlated CO2 covariance matrix.

The fixed 1-sigma values, from TOA to surface, are:

```text
1.437  2.515  3.233  3.982  5.552  6.594  7.447  8.199  9.174  10.056
12.188 14.150 17.149 19.941 23.972 27.758 32.716 37.406 42.677 47.690 ppm
```

The covariance is strongly correlated through most of the troposphere. For
example, the adjacent-level correlations from ACOS pressure levels 6 through
18 are approximately 0.95--0.97. The complete 20 by 20 block is positive
definite. Its eigenvalues span approximately 0.945 to 7,589 ppm2, giving a
condition number of about 8,031 before state scaling.

Using each sounding's dry-air subcolumns to form pressure-level XCO2 weights,
this covariance gives a typical CO2-only a priori XCO2 uncertainty of about
13.9 ppm. Retaining its diagonal variances but discarding its correlations
would give only about 4.6 ppm. The off-diagonal terms therefore do not merely
smooth the profile: they preserve uncertainty in vertically coherent profile
shifts that affect XCO2.

The climate-model provenance needs a precise qualification. O'Dell et al.
(2012) state that the ACOS **mean prior profiles** came from an LMDz forward
model run, aggregated by month, latitude, and land/ocean class and adjusted to
GLOBALVIEW surface CO2. They describe the nonzero covariance off-diagonals as
an imposed smoothness constraint, with diagonal variability decreasing with
altitude and the matrix scaled to give an approximately 12 ppm XCO2 prior
uncertainty. The publication does not say that this fixed covariance was
itself calculated directly from a climate-model ensemble. O'Dell et al.
(2018) likewise state that CO2 is the exception to ACOS's otherwise diagonal
a priori covariance.

References:

- [O'Dell et al. (2012), ACOS algorithm description](https://doi.org/10.5194/amt-5-99-2012), especially Sect. 2.1 and Fig. 2.
- [O'Dell et al. (2018), OCO-2 ACOS version 8](https://doi.org/10.5194/amt-11-6539-2018), especially Sect. 3 and Table 2.
- [Kulawik et al. (2019), OCO-2 error-analysis validation](https://doi.org/10.5194/amt-12-5317-2019), Sect. 2.1.

## Sounding-dependent ACOS prior means

The covariance is fixed, but the mean profiles are not. Dry-column-weighted
prior XCO2 statistics in the four files are:

| Orbit/scene | Mean (ppm) | Standard deviation (ppm) | Range (ppm) |
|---|---:|---:|---:|
| TallVeg | 414.802 | 0.996 | 413.408--417.810 |
| Himalayas | 417.788 | 1.001 | 415.524--419.687 |
| Sahara | 417.659 | 2.471 | 414.874--422.391 |
| Andes | 416.850 | 2.042 | 414.279--421.958 |

The current synthetic experiment deliberately uses a common active-layer
prior mean of 400 ppm. Changing both the mean and covariance at the same time
would alter the scientific experiment and confound their effects. The
implemented change therefore retains the 400 ppm mean and changes only the CO2
covariance.

## Comparison with the current RRS-XCO2 prior

The current prior in `build_apriori.jl` is diagonal in CO2 and assigns 0, 4, 8,
or 40 ppm according to layer-center altitude. The upper four layers are fixed;
the lower 12 are active. Its dry-column-weighted XCO2 uncertainty is only
**3.689 ppm**.

To compare like with like, the fixed ACOS pressure-level covariance was mapped
to the RRS-XCO2 layer centers. Let `r` be normalized pressure and let `H`
linearly interpolate the 20 ACOS pressure-level VMRs to the 16 RRS-XCO2 layer
centers. The mapped covariance is

```text
S16 = H * S20_ACOS * transpose(H).
```

This is also the layer-average mapping for a profile represented linearly in
pressure, apart from the very small within-layer water/gravity variation. The
comparison is:

| RRS layer | Center (km) | Center (hPa) | Current sigma (ppm) | Mapped ACOS sigma (ppm) |
|---:|---:|---:|---:|---:|
| 5 | 9.326 | 281.351 | 4 | 6.839 |
| 6 | 7.991 | 343.842 | 4 | 7.795 |
| 7 | 6.854 | 406.333 | 4 | 8.851 |
| 8 | 5.850 | 468.824 | 4 | 9.954 |
| 9 | 4.945 | 531.316 | 4 | 12.340 |
| 10 | 4.127 | 593.807 | 4 | 14.877 |
| 11 | 3.370 | 656.298 | 4 | 18.304 |
| 12 | 2.663 | 718.789 | 8 | 22.375 |
| 13 | 2.016 | 781.281 | 8 | 27.033 |
| 14 | 1.405 | 843.772 | 8 | 32.809 |
| 15 | 0.820 | 906.263 | 40 | 37.758 |
| 16 | 0.268 | 968.754 | 40 | 43.527 |

The existing setup is therefore substantially tighter than ACOS from about
1--9 km. It is comparable only in the two lowest layers. More importantly,
its zero cross-layer covariances suppress vertically coherent CO2 changes.

For the active 12 layers:

| Prior construction | CO2-only sigma(XCO2) |
|---|---:|
| Current diagonal RRS-XCO2 prior | 3.689 ppm |
| Mapped ACOS diagonal only | 5.074 ppm |
| Mapped ACOS covariance, all correlations retained | 13.716 ppm |
| Mapped ACOS conditional covariance given exact upper four layers | 11.011 ppm |

The recommended correlated prior is 3.72 times looser in XCO2 standard
deviation (13.82 times larger in variance) than the current one.

## Implemented construction

1. Keep the active-layer mean at 400 ppm so this campaign isolates the effect
   of the covariance change.
2. Map the fixed 20-level ACOS covariance onto the 16 RRS-XCO2 pressure grid
   using `S16 = H*S20*H'` in normalized pressure.
3. Continue fixing layers 1:4 to the current scene truth for synthetic closure.
4. Use the **marginal** active block `S16[5:16,5:16]`, not a covariance
   conditioned on the fixed upper values. Conditioning would use knowledge of
   the synthetic truth to tighten and shift the active prior, which is not
   information an observational retrieval would possess.
5. Preserve zero a priori covariance between CO2 and non-CO2 state elements,
   matching the four orbit files.
6. Store the source grid, interpolation matrix, mapped covariance, and its
   provenance in the generated prior NetCDF so every retrieval is auditable.

The mapped 12 by 12 active block is positive definite. Its raw condition
number is about 11,470 and the corresponding correlation-matrix condition
number is about 1,392. The retrieval already solves in prior-standard-deviation
scaled coordinates and uses Float64, so this is numerically workable, but the
regenerated prior and one representative retrieval should still receive an
explicit Cholesky/closure test before launching the suite.

## Implementation record

1. `build_apriori.jl` now reads the approved `acos_mapped` lower triangle from
   `co2_prior_covariances.dat`, converts ppm2 to VMR2, and rejects a missing,
   duplicate, non-positive-definite, or incorrectly fixed matrix.
2. `apriori_states.nc` and `apriori_states.dat` were regenerated. The NetCDF
   stores the covariance model, source, mapping, pressure nodes, and implied
   XCO2 uncertainty.
3. Regression tests verify positive definiteness, representative marginal
   sigmas, nonzero off-diagonals, and the exact 13.715633 ppm column
   uncertainty. Retrieval state mapping and the truth solar/spectral operator
   tests also pass.
4. The obsolete Curry worker was stopped without touching the truth-map
   workers. Its 53 corrected and 52 uncorrected products, old generated prior,
   manifests, logs, and analyses are preserved in `../old/`.
5. Replacement manifests contain 704 planned solves. The 352 no-aerosol
   solves have ready inputs and were launched/queued in a disjoint 3:1
   Wurst:Curry state allocation. The aerosol half is marked `inputs_ready=0`
   until its instrument products exist.
