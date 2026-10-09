# Round 7: fixed SIF, tight UTLS, imperfect Raman correction

Approved experiment: retain the Round-5 retrieval model and priors, replacing
the ideal corrected O2 observations by an independently LUT-corrected version
of the existing uncorrected observations. Round 5 and Round 6 are read-only.

## Scientific definition

For each of the 80 scenes, average the three named `o2a_surface_P0/P1/P2`
coordinates over **all ten valid Round-5 uncorrected noisy members 01:10**.
Do not include noiseless member 11, average across scenes, or condition this
selection on the success of the corresponding ideal-corrected retrieval.
All 800 required noisy uncorrected members passed completeness, convergence,
and fit-quality checks when this release was prepared.

Reconstruct the albedo polynomial in its original complete-band wavenumber
coordinate. Surface pressure and geometry are independent scene metadata:
1000 hPa, SZA 30 degrees, VZA 0 degrees, relative azimuth 0 degrees. Neither
retrieved pressure nor truth albedo is used to estimate the correction.

At each native high-resolution LUT wavenumber, evaluate

```
Delta_stokes = LUT_Cabannes + LUT_RRS - LUT_Rayleigh
correction_observation = H(Delta_stokes)
y_round7 = y_round5_uncorrected - correction_observation
```

The **same aerosol-free, SIF-off LUT** is used for every combination of
aerosol and SIF truth. Interpolation is linear in surface pressure, cos(SZA),
and the reconstructed albedo at that wavenumber. Exact pressure nodes are
used without consulting adjacent slices. The current 1000-hPa slice is
complete; no extrapolation or clipping is permitted.

The SZA brackets are approximately 26.002838 and 31.988062 degrees, with
cos(SZA) interpolation weights 0.35300449 and 0.64699551, respectively.
The albedo coordinate remains based on the original 757--773 nm truth band,
not the LUT's wider domain. The correction is evaluated on native LUT nodes
within that band, avoiding an unnecessary interpolation between two offset
0.1-cm^-1 grids. Wavelength is computed as `1e7/Float64(wn)`.

`H` is the unchanged Round-5 synthetic OCO operator: signed analyzer
`M11*I - M12*Q + M13*U`, per-wavenumber to per-nm density conversion,
wavelength-space Gaussian convolution with FWHM 0.04 nm and six-sigma support,
and sampling at the 934 O2 observation centers. The Python implementation
used during generation is independently compared against the pinned Julia
`SyntheticOCO2.process_stokes_spectrum` for all 80 generated corrections.
There is no additional M11 or pi normalization. Both CO2 bands are unchanged.

One correction is frozen per scene and applied to all eleven realizations.
The exact uncorrected noise draw, injected noise, standard deviation, and
diagonal covariance are retained. They are **not** taken from the old
ideal-corrected measurement class, and no new noise is generated. Retrievals
treat this estimated correction as fixed; uncertainty in its coefficient
estimate and resulting inter-member correlations are not propagated into Se.

## Baseline and known approximations

The forward/OE implementation is the clean accelerated Round-5 checkout at
`7acab57000dae259207a6760faae156cfa1734f6`, Julia 1.12.5, Float32 GPU, nine
streams, 16 layers, and 28 active coordinates. The exact original Round-5
priors are reused: fixed SIF759 and mSIF, UTLS log-AOD sigma 0.075 and
log-height sigma 0.01. CO2 and other priors are unchanged.

The LUT has a 12-layer AFGL atmosphere and historical numerical settings;
it is not solver/profile matched to the 16-layer truth. Its parent contains
promoted chunks from multiple commits. The **whole-file SHA256**, rather than
the parent's commit attribute alone, identifies the LUT:

```
502a3070ff6db5041285222172a2584e89c7845dca8841adf21538097b1fe1b2
```

Using spectral albedo to query a scalar-albedo LUT is itself approximate,
because Raman scattering couples different wavelengths. Residual errors
therefore cannot be attributed exclusively to albedo retrieval error, SIF,
or aerosols. No LUT regeneration or atmospheric retuning is part of Round 7.

## Saved inputs and controls

The isolated campaign directory is
`RRS_XCO2/bottom_layer_XCO2_retrievals/round7_fixed_sif_imperfect_correction/`.

`observations/OCO2round7_NNN.nc` retains:

- uncorrected, ideal-corrected reference, and imperfectly corrected noiseless
  observations;
- all eleven original uncorrected and new imperfectly corrected observations;
- exact noise draws, injected noise, standard deviation, covariance and seeds;
- correction before and after the instrument, interpolated LUT components,
  albedo spectrum, coefficient mean, and the ten input coefficient vectors;
- metadata, SIF-v2 provenance, source retrieval paths/checksums, and LUT inputs.

`observations/provenance.json` records the 80-scene audit and all source hashes.
`observations/SHA256SUMS` seals the 80 NetCDF files and that audit. Existing
observation releases are never overwritten.

Some old local SIF-on truth/noise products still use the obsolete `total_0p5`
definition. **Round 7 does not read those products.** It reconstructs inputs
from the self-contained, completed, definition-v2 Round-5 retrieval outputs,
which contain the actual measurements and noise used on Gattaca. The new
checksummed input release and SIF-v2 validation replace the old raw-truth
publication gate; this is not an override allowing legacy SIF inputs.

New outputs use `retrievals_nosif/corrected/` and
`retrievals_sif/corrected/` for filename compatibility, with explicit Round-7
and `imperfectly_corrected` provenance. `corrected` no longer means ideal
Rayleigh in this campaign. Existing Round-5 uncorrected outputs are reused
as controls, never relabeled as newly executed Round-7 retrievals.

## Initial observation audit

Differences from the ideal corrected reference, in band-RMS instrument-noise
units (no retrieval fit is involved in these numbers):

| Truth subset | Before correction | After LUT correction |
|---|---:|---:|
| No aerosol, no SIF | 0.2771--0.3438 | 0.00928--0.01142 |
| No aerosol, with SIF | 0.2770--0.3438 | 0.00933--0.01147 |
| Aerosol, no SIF | 0.2963--0.3554 | 0.01496--0.02425 |
| Aerosol, with SIF | 0.2962--0.3553 | 0.01510--0.02434 |

These are diagnostics of this particular LUT/campaign, not claims of general
retrieval accuracy or solver equivalence.

## Execution and safety

There are 880 new imperfectly corrected solves (80 scenes x 11 realizations).
The 880 uncorrected controls need no rerun. Curry physical GPU 1 receives all
no-aerosol scenes; Wurst physical GPU 1 receives all aerosol scenes. Each
worker runs no-SIF and SIF blocks sequentially, 440 solves per machine.

The launcher verifies the pinned source, sealed release/input manifests,
prior identities, host, GPU UUID and GPU vacancy. It binds physical GPU 1 by
UUID, with process-local CUDA ordinal 0, and holds a per-worker lock. A busy
GPU causes refusal; another user's process must never be stopped or displaced.
Before full launch, CPU preflight and small noiseless retrieval tests must
pass. Resume verifies existing outputs rather than overwriting them.
