# Synthetic OCO-2 measurement operator

This directory contains the common instrument-processing stage used by both
optimal-estimation retrieval classes. It converts each high-resolution
truth-map Stokes spectrum into a scalar, synthetic OCO-2-like measurement.
The corrected and uncorrected retrievals must use this identical operator.

## Synthetic bands

The wavelength limits, sampling intervals, and Gaussian FWHM values are the
OCO-2 entries in Table 1 of
`RRS_XCO2/resources/The_Cross-Calibration_of_Spectral_Radiances_and_Cr.pdf`.
The real, sounding-dependent OCO-2 ILS is deliberately replaced by a Gaussian.

| Band | Synthetic grid (nm) | FWHM (nm) | Sigma (nm) | Samples |
|---|---:|---:|---:|---:|
| O2 A-band | `758.0:0.015:772.0` | 0.04 | 0.0169864360 | 934 |
| Weak CO2 | `1594.0:0.031:1619.0` | 0.08 | 0.0339728720 | 807 |
| Strong CO2 | `2042.0:0.04:2082.0` | 0.10 | 0.0424660900 | 1001 |

Julia's range semantics retain the last point not exceeding the stop. Thus the
O2 grid ends at 771.995 nm and the weak-CO2 grid at 1618.986 nm; the strong
grid lands exactly on 2082.0 nm.

## OCO Stokes analyzer

OCO L1B sounding files store
`FootprintGeometry/footprint_stokes_coefficients` with four weighting factors
for each band, footprint, and frame. The dataset description calls them
"weighting factors applied to the Stokes parameters calculated by the
radiative transfer code to compute the radiance." They are an analyzer row,
not a complete invertible 4x4 Mueller matrix.

Following
`/home/sanghavi/code/github/OCORaman/OCORaman/src/OCOPlots/oco_gain.jl`
exactly, the scalar response is

```text
y_analyzer(lambda) = M11 I(lambda) - M12 Q(lambda) + M13 U(lambda).
```

There is no subsequent division by `M11`: the OCO coefficients already carry
the analyzer normalization (`M11 = 0.5`). The current nadir truth spectra have
three Stokes components, and the L1B `M14` values used here are zero.

### Representative-row estimator

The coefficients rotate with sounding geometry. A component-wise signed mean
would shrink `M12` and `M13` toward zero even though every individual analyzer
has `hypot(M12,M13)/M11` essentially equal to one. To avoid this cancellation:

1. pool every valid footprint/frame from the four available nadir L1B files;
2. compute the analyzer direction `atan(M13/M11,M12/M11)`;
3. take its circular median;
4. select the actual observed coefficient row nearest that median.

This gives every sounding equal weight, preserves signs, and preserves the
physical analyzer magnitude. The generated NetCDF and human-readable table
record source filenames, estimator metadata, sample counts, the raw signed
mean (diagnostic only), and the selected rows:

| Band | M11 | M12 | M13 | M14 | Rows pooled |
|---|---:|---:|---:|---:|---:|
| O2 A-band | 0.5 | -0.470727354288 | -0.168569713831 | 0 | 254024 |
| Weak CO2 | 0.5 | -0.478100210428 | -0.146356388927 | 0 | 254024 |
| Strong CO2 | 0.5 | -0.478100597858 | -0.146355077624 | 0 | 254024 |

Regenerate these resources with:

```bash
julia --project=. \
  RRS_XCO2/inversion/instrument/derive_representative_stokes_coefficients.jl
```

Set `OCO_L1B_ROOT` if the L1B files move. `OCO_L1B_MODE=all` is available for
diagnostics, but `nadir` is the production default because the truth scenes
are nadir observations.

## Processing order and units

For each independently saved truth component, the operator performs:

1. the signed OCO analyzer projection above;
2. the spectral-density conversion
   `L_per_nm = L_per_cm-1 * 1e7/lambda_nm^2`;
3. convolution with a normalized wavelength-space Gaussian;
4. direct evaluation at the synthetic sample centers.

The truth basis is uniform in wavenumber and consequently nonuniform in
wavelength. `SyntheticOCO2.jl` sorts it in wavelength, uses trapezoidal
wavelength integration weights in the Gaussian numerator and denominator,
and therefore preserves a constant spectrum. It does not assume a uniform
wavelength grid or first interpolate onto an arbitrary dense grid.

The Gaussian is evaluated through +/-6 sigma by default. Production processing
requires full source coverage over every kernel; clipped edge kernels are not
silently renormalized.

The retrieval forward model and its analytic Jacobians use this same fixed
operator. `process_stokes_jacobian` accepts the canonical high-resolution
layout `(stokes, wavelength, parameter)` and returns `(sample, parameter)`.
Since the representative analyzer, Gaussian ILS, and sample grid are fixed,
the instrument-stage Jacobian is simply the instrument operator applied to
each forward-model derivative column.

## Supplemental convolution shoulders

The current source ranges cover the +/-6 sigma support for the O2 A-band and
weak CO2 band. The strong-band basis starts at 2042.0396945 nm, while the first
synthetic sample is 2042.0 nm and requires support to 2041.7452035 nm.
Therefore eight additional 0.1 cm^-1 samples are required on the strong-band
short-wavelength side, spanning approximately 2041.7062--2041.9980 nm.

The missing points were first generated as a separately validated staging
dataset with:

```bash
CUDA_DEVICE=1 RRS_XCO2_FLOAT_TYPE=Float32 \
  julia --project=. \
  RRS_XCO2/scripts/generate_truth_map_convolution_shoulders.jl
```

The generator refuses to run before all 32 aerosol truth files are marked
complete. Surface Legendre polynomials retain the complete base band's
coordinate when evaluated at the supplemental points. Once validated, merge
the staging files atomically into all 64 canonical scenes and both wavelength
files with:

```bash
julia --project=. RRS_XCO2/scripts/merge_convolution_shoulders.jl
```

The resulting `strong_co2` dimension contains 995 points spanning about
2041.7062--2084.0001 nm. The first 987 values are the unchanged base spectrum;
the last eight are the appended convolution shoulder. The staging directory
may be removed only after the integrated files pass the full validation.

## Measurement classes and outputs

`process_truth_map.jl` writes `OCO2sims_NNN.nc` files under
`RRS_XCO2/truth_map/OCO_radiances/` by default.

Processed variables use the explicit form
`I_OCO_<component>_<band>`. For example, an O2 A-band scene contains
`I_OCO_rayleigh_o2a`, `I_OCO_cabannes_o2a`, `I_OCO_rrs_o2a`,
`I_OCO_corrected_o2a`, and `I_OCO_uncorrected_o2a`. Each wavelength
coordinate records the sampling interval, Gaussian FWHM, and Gaussian sigma.

- O2 corrected measurement: processed Rayleigh.
- O2 uncorrected measurement: processed Cabannes + processed RRS.
- Weak/strong CO2 corrected and uncorrected measurements: the same processed
  Rayleigh simulation, because RRS is excluded from those truth bands.

Processed Rayleigh, Cabannes, and RRS O2 components are retained separately
for quality control. Output radiance units are `mW m-2 sr-1 nm-1`.

The aerosol truth-map comparison reveals an approximately 23.5% band-RMS
polarization effect on the elastic aerosol signal for the current AOD760 =
0.28 case, while the already depolarized pure RRS component changes by less
than 1% under the analyzer. The interpretation, numerical decomposition,
polarimetric opportunity, and geometry-specific caveats are recorded in
`RRS_XCO2/truth_map_aerosols/MEASUREMENT_FINDINGS.md`.

Production command, after truth and shoulder completion:

```bash
julia --project=. RRS_XCO2/inversion/instrument/process_truth_map.jl
```

The script requires completed truth scenes, finite data, and full Gaussian
support. Completion is recorded with the general `simulation_complete=1`
attribute; existing chunked aerosol files remain compatible through their
`chunked_simulation_complete=1` attribute. Development overrides exist but are
labeled unsafe in the script help and must not be used to generate retrieval
measurements.

## Validation

Run the focused tests with:

```bash
julia --project=. RRS_XCO2/inversion/instrument/test_synthetic_oco2.jl
```

After production processing, validate all 64 NetCDF products with:

```bash
julia --project=. RRS_XCO2/inversion/instrument/validate_oco_radiances.jl
```

Compare their absolute radiance scale with the four nadir OCO-2 L1B files
used by `oco_gain.jl` with:

```bash
python3 RRS_XCO2/inversion/instrument/compare_oco2_observed_radiances.py
```

The method, unit audit, numerical ranges, and one bright-desert strong-band
flag are documented in
`RRS_XCO2/truth_map/OCO_radiances/validation_against_OCO2/README.md`.

The tests cover grid endpoints/counts, the exact `oco_gain.jl` signed analyzer
formula, rejection of determinant normalization, the per-wavenumber to per-
wavelength Jacobian, constant/affine preservation on a reversed nonuniform
source grid, finite-difference verification of instrument-processed Jacobian
columns, and rejection of the currently under-supported strong-band edge.

The 0.1 cm^-1 strong-band basis spacing is about 0.0425 nm, slightly coarser
than the 0.04 nm output sampling and about 2.35 basis points per 0.10 nm FWHM.
This follows the selected truth-basis resolution, but the convolved strong-band
spectra should receive an explicit convergence comparison before retrievals
are treated as final.

## OCO-2 measurement-noise covariance

The state-dependent diagonal measurement-noise covariance follows Eq. (3-8)
of `RRS_XCO2/resources/OCO_L1B_ATBD.pdf`. The equation is evaluated in photon
radiance with the Table 3-5 MaxMS values and the wavelength-dependent photon
and background coefficients from `InstrumentHeader/snr_coef`. Table 3-6 MinMS
is retained as a dynamic-range check, not inserted into the equation as a
floor.

Because the truth map has no unique footprint or acquisition date, the
representative coefficient resource is a pointwise median of 32 coefficient
spectra: four nadir L1B files used by `oco_gain.jl` times eight footprints.
Generate it, test the equation, generate all scene covariances, and validate
them with:

```bash
julia --project=. RRS_XCO2/inversion/instrument/derive_representative_snr_coefficients.jl
julia --project=. RRS_XCO2/inversion/instrument/test_oco2_noise.jl
FORCE=1 julia --project=. RRS_XCO2/inversion/instrument/generate_noise_covariances.jl
julia --project=. RRS_XCO2/inversion/instrument/validate_noise_covariances.jl
```

Each `truth_map/OCO_radiances/noise_covariances/OCO2noise_NNN.nc` stores the
2742-element corrected and uncorrected measurement vectors, NEN vectors, SNR,
dynamic-range flags, and `Se_diagonal_corrected` / `Se_diagonal_uncorrected`.
The dense matrix is reconstructed only when needed as `Diagonal(Se_diagonal)`.
Full provenance, schema, numerical ranges, and limitations are documented in
`truth_map/OCO_radiances/noise_covariances/README.md`.
