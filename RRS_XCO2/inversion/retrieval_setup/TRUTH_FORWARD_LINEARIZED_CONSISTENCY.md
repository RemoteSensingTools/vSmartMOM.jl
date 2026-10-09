# Truth, forward-model, and analytic-Jacobian consistency contract

This document is the retrieval-readiness guide for the RRS–XCO2 experiment.
It records which quantities must be identical between the truth generator and
the retrieval forward/linearized model, how each equality is tested, the
accepted numerical tolerances, and every intentional difference. A retrieval
campaign is not ready merely because all NetCDF files exist: every required
gate in the final section must pass.

The production retrieval Jacobian is the analytic Jacobian from linearized
vSmartMOM. Finite differences appear below only as independent regression
tests; they are never used to generate retrieval Jacobians.

## Current status

Status on 2026-08-29: **not yet ready for retrieval production**.

- The canonical truth builder and retrieval model agree closely in all three
  bands for an aerosol-on, SIF-off closure state.
- The new weak- and strong-CO2 truth spectra use the common ABSCO/profile
  configuration and have passed completion/provenance validation.
- The archived aerosol A-band spectra predate exact retained-core grid
  construction. Their known Float32 node displacement causes the old active
  truth-product closure to fail. Exact-grid A-band regeneration is in progress
  with `scripts/regenerate_o2_preserve_co2.jl`; an interpolation diagnostic is
  not an approved replacement.
- The current closure report is therefore correctly marked **FAIL**. Do not
  start retrieval production until a newly generated A-band makes it pass.

The authoritative numerical report is
[`ABSCO_CLOSURE.md`](ABSCO_CLOSURE.md). The report must say `Overall result:
PASS` before retrievals are released.

## Shared physical and numerical inputs

Except for the intentional differences listed later, truth and retrieval must
use the following common configuration.

| Quantity | Required common value or construction | Check |
|---|---|---|
| Geometry | SZA 30 deg, VZA 0 deg, relative azimuth 0 deg | NetCDF provenance and model parameters |
| Solar representation | `external_solar=true`; scalar solar direction and dedicated Z0/R0/T0 source columns | closure model construction |
| Precision | Float32 for production truth and retrieval | model parameters and output provenance |
| Polarization | Stokes I,Q,U | array shape and instrument test |
| Surface pressure | 1000 hPa truth; retrieved from a 1000-hPa prior | state mapping and exact profile checks |
| Vertical grid | one common materialized 16-layer p/T/q profile | element-for-element profile comparison |
| Humidity | identical specific humidity `q`, H2O VMR, wet/dry columns | element-for-element profile comparison |
| O2 VMR | 0.21 | model profile and NetCDF provenance |
| CO2 truth | uniform 380, 400, 420, or 440 ppm | truth-state table |
| CO2 retrieval state | layers 5:16 active; layers 1:4 fixed to the current scene's truth VMR | state-mapping regression |
| O2 spectroscopy | ABSCO v5.2 `o2_v52_v2.jld2` | table path and optical-depth closure |
| A-band H2O | rebuilt HITRAN `H2O.jld2`, driven by the common `q` profile | table path and optical-depth closure |
| Weak CO2 spectroscopy | ABSCO v5.2 weak-band CO2 and H2O tables | table paths and optical-depth closure |
| Strong CO2 spectroscopy | ABSCO v5.2 strong-band CO2 and H2O tables | table paths and optical-depth closure |
| Solar spectrum | the same high-resolution `solar.out` interpolator multiplied by the same Planck spectrum and scale | exact source-array regression |
| Aerosol microphysics | the same size distributions, refractive indices, fixed widths, and three vertical-profile families | shared YAML and optical closure |
| Aerosol state mapping | retrieval `ln(AOD760)` converted through the tabulated fixed `tau_ref(550)/AOD760` factors; `ln(z0)` converted to physical height | mapping and finite-difference regressions |
| Aerosol optics | full-band endpoint/reference-node interpolation, with canonical band anchors independent of Raman chunks | optical closure and shoulder invariance test |
| Truncation | delta-BGE with `l_trunc=16`, `stream_l_cap=17`, and `greek_beta_cutoff=1e-5` for aerosol production | shared configuration |
| Surface | three Legendre coefficients per band, defined on each canonical physical band | surface-coordinate regression |
| Source | identical physical solar source; SIF handled as described under intentional differences | source regression |

The shared definitions live primarily in
[`scripts/common.jl`](../../scripts/common.jl). Truth construction is in
[`scripts/generate_truth_map.jl`](../../scripts/generate_truth_map.jl) and
[`scripts/generate_truth_map_aerosol_chunked.jl`](../../scripts/generate_truth_map_aerosol_chunked.jl).
The A-band-only correction that preserves the accepted CO2 arrays is
[`scripts/regenerate_o2_preserve_co2.jl`](../../scripts/regenerate_o2_preserve_co2.jl).
The retrieval adapter is
[`VSmartMOMForward.jl`](../VSmartMOMForward.jl).

## Angular-grid contract

`external_solar=true` keeps SZA outside the square diffuse operator. The
current nadir VZA is not a nine-stream Gauss-Legendre node and is appended as
one zero-weight output node.

- Aerosol truth and retrieval: 9 weighted streams + nadir VZA = 10 diffuse
  nodes, hence a 30 x 30 IQU operator.
- Aerosol-free truth: 5 weighted streams + nadir VZA = 6 diffuse nodes, hence
  an 18 x 18 IQU operator.
- Aerosol-free retrieval: 9 weighted streams + nadir VZA = 10 diffuse nodes.

The five-versus-nine-stream aerosol-free difference is deliberate. Rayleigh
and Lambertian components have sufficiently low angular support, but a
five-stream truth versus nine-stream retrieval comparison must still pass the
final radiance tolerance after A-band regeneration.

Changing SZA alone cannot remove the tenth diffuse node. The extra node is the
nadir VZA, not the external solar direction.

## Spectral-grid and Raman-shoulder contract

There is one canonical retained O2 A-band grid. Truth and retrieval must never
construct independent nominally equivalent Float32 ranges.

1. noRS retrieval calculations use the 2,735-point canonical retained grid
   directly.
2. RRS truth calculations add 234 cm-1 shoulders on both sides of each core
   chunk.
3. Only the shoulders are constructed. The canonical core array is inserted
   verbatim into the expanded solve grid.
4. Cabannes and RRS use the same expanded grid and retained index range. Pure
   monochromatic Rayleigh/noRS uses those exact retained core nodes directly,
   without evaluating the discarded shoulders.
5. The required invariant is exact, not approximate:

   ```julia
   rrs_solve_grid[keep] == noRS_grid
   ```

6. Retrieval noRS does not evaluate the solve-only Raman shoulders.

The shared constructor is `RRSXCO2Common.raman_solve_grid`. The boundary
regression in [`test_forward_state_mapping.jl`](../test_forward_state_mapping.jl)
passes for both Float32 and Float64 (6/6 assertions).

The historical construction started one long Float32 range at the left Raman
shoulder. Its changed range origin displaced some internal retained nodes by
one Float32 ULP, at most 9.765625e-4 cm-1. The effect recurs across many
64-point chunks and is strongest on steep O2 and Fraunhofer line slopes; it is
not confined to one wavelength.

The weak CO2 band uses its canonical 1,281-point grid. The strong CO2 solve
grid contains the 987-point physical surface-definition band plus eight
short-wavelength convolution-support nodes, for 995 points total. The surface
coefficient transform preserves the original three-term polynomial on that
extended grid.

## Forward versus linearized-forward closure

The linearized run returns both its forward Stokes spectrum and the analytic
Jacobian. Its forward spectrum is compared with a separate nonlinear noRS
run before any instrument processing. These arrays are not expected to be
bit-for-bit identical in Float32 because the fused linearized kernels have a
different floating-point evaluation order. They must satisfy the numerical
closure tolerance.

For aerosol state 009, the 2026-08-29 closure produced:

| Band | canonical truth vs retrieval noRS, relative L2 | retrieval noRS vs linearized forward, relative L2 | OCO truth vs retrieval, relative L2 | OCO retrieval vs linearized forward, relative L2 |
|---|---:|---:|---:|---:|
| O2 A | 1.2582e-6 | 7.7538e-6 | 4.3962e-7 | 7.3108e-6 |
| Weak CO2 | 1.8669e-7 | 5.4412e-7 | 9.6338e-8 | 5.2831e-7 |
| Strong CO2 | 2.3154e-7 | 1.5256e-7 | 1.3138e-7 | 1.4013e-7 |

The current acceptance limits are:

- canonical truth versus retrieval: relative L2 <= 2e-5;
- retrieval nonlinear versus linearized forward: relative L2 <= 1e-5;
- the same limits after OCO processing.

The optical-property comparisons for the same closure were:

| Band | tau_abs relative L2 | tau_Rayleigh relative L2 | tau_aerosol relative L2 |
|---|---:|---:|---:|
| O2 A | 3.0368e-8 | 0 | 2.9271e-8 |
| Weak CO2 | 3.6467e-8 | 0 | 1.4519e-7 |
| Strong CO2 | 4.2064e-8 | 0 | 1.3350e-7 |

All T, p_half, p_full, q, H2O VMR, dry-column, O2 VMR, and CO2 VMR arrays
matched element-for-element in all three bands.

## Analytic-Jacobian checks

Production Jacobians are computed analytically by linearized vSmartMOM. The
analytic derivatives pass through elemental, doubling, interaction, source,
and surface operations. `OCO_RRS_synth` selects only the retrieval parameters;
inactive derivatives are not propagated and selected compact columns agree
with the corresponding full-model columns.

The current selective-Jacobian regression passed on 2026-08-29:

- overflow-safe Float32 pressure-column tangent: 5/5;
- aerosol/gas analytic derivatives versus central finite differences: 3/3;
- pressure/surface/SIF analytic derivatives versus central finite differences:
  7/7;
- OCO retrieval plan layout and mapping: 18/18;
- compact analytic propagation versus selected full columns: 9/9.

Recorded Float64 CPU finite-difference errors are:

| Analytic parameter class | Maximum absolute error | Maximum relative error |
|---|---:|---:|
| aerosol `tau_ref` | 1.84e-10 | 2.77e-9 |
| aerosol physical `z0` | 3.39e-12 | 1.06e-6 |
| selected layer absorption | 1.08e-11 | 2.50e-8 |
| surface pressure | 8.84e-12 | 8.79e-6 |
| surface P0 | 1.17e-12 | 1.78e-7 |
| surface P1 | 1.17e-12 | 1.78e-7 |
| surface P2 | 1.17e-12 | 1.78e-7 |
| SIF760 | 1.69e-11 | 8.74e-7 |
| mSIF | 1.32e-9 | 4.58e-7 |

These finite differences validate each derivative implementation; they are
not a replacement for the analytic retrieval Jacobian.

## Mueller processing, convolution, and resampling

Truth spectra and every analytic Jacobian column use the same functions in
[`SyntheticOCO2.jl`](../instrument/SyntheticOCO2.jl), in the same order:

1. OCO analyzer projection, `M11*I - M12*Q + M13*U`;
2. no production normalization of the analyzer response;
3. conversion from radiance per cm-1 to radiance per nm using
   `1e7/lambda_nm^2`;
4. normalized wavelength-space Gaussian convolution with trapezoidal source
   weights and six-sigma support;
5. evaluation at the synthetic OCO sample centers.

The Gaussian FWHM values are 0.04, 0.08, and 0.10 nm for O2 A, weak CO2, and
strong CO2. The corresponding synthetic grids contain 934, 807, and 1,001
samples.

`process_truth_map.jl` applies `process_stokes_spectrum` to truth radiances.
`VSmartMOMForward.evaluate_oco_forward` applies that same function to the
linearized forward radiance and `process_stokes_jacobian` to every analytic
column. The instrument operator is linear and its Jacobian test passed all 25
current assertions. In the full closure, the worst processed-column error was
4.06e-12 absolute and 2.82e-15 relative in the A band; the weak and strong
bands were similarly small.

The display-only division by 0.5 used in diagnostic plots is never applied to
the measurement vector, analytic Jacobian, or measurement covariance.

## Raman shoulder checks

The closure suite separately requires:

- exact retained canonical nodes in the corrected grid constructor;
- zero change in OCO convolution when solve-only shoulder samples are present
  but removed by `keep` before instrument processing;
- canonical full-band versus a production-size 256-point noRS core solved
  with +/-234 cm-1 shoulders: relative L2 <= 2e-5.

The production Rayleigh/noRS component is consequently evaluated on the
retained core alone. Cabannes and RRS remain on the shoulder-expanded solve;
the canonical full-band aerosol and surface anchors are identical in both
paths.

The most recent corrected-constructor results were:

- retained-node maximum difference: 0 cm-1;
- convolution difference caused by discarded solve-only shoulders: 0;
- full-band versus shouldered-core elastic relative L2: 1.4042e-6.

These checks concern the current constructor. They do not retroactively repair
the archived aerosol A-band product.

## Intentional differences

The following are part of the experiment and must not be misdiagnosed as
implementation inconsistencies.

1. **RRS measurement versus noRS retrieval.** Corrected measurements use the
   truth Rayleigh spectrum. Uncorrected measurements use truth Cabannes+RRS.
   Both retrieval classes deliberately fit these measurements with the same
   linearized noRS forward model. The uncorrected model discrepancy is the
   signal under study.
2. **SIF representation.** Truth SIF retains the supplied high-resolution
   spectral shape and is normalized at 760 nm by
   `2pi * L_lambda(760 nm) = 0.5 mW m^-2 nm^-1`. This is an unweighted
   upward-solid-angle integral at one wavelength, not a wavelength-integrated
   SIF area. Retrieval SIF is the two-parameter line
   `SIF760 + mSIF*(nu-nu760)`. For a SIF-off suite both are exactly zero;
   SIF-on retrievals intentionally include the line-shape approximation.
3. **Aerosol-free stream count.** Aerosol-free truth uses five weighted
   streams; all retrievals use nine. Aerosol truth and retrieval both use
   nine.
4. **Retrieved state.** Surface pressure, lower 12 CO2 layers, aerosol AODs
   and heights, surface coefficients, SIF760, and mSIF change during retrieval.
   Those state changes are the purpose of the forward model, not a truth-model
   inconsistency.
5. **Upper CO2 layers.** The top four CO2 layers are excluded from the active
   state because their prior variance is zero. They are set to the current
   truth scene's uniform CO2 value before every retrieval.

## Known failing product check

The copied archived aerosol A-band currently differs from the aligned
retrieval noRS spectrum by:

- high resolution: relative L2 6.5613e-4, maximum absolute difference
  5.1967e-2;
- after OCO processing: relative L2 3.6182e-4, maximum absolute difference
  6.8451e-2.

This fails the closure thresholds. The differences are distributed around
many steep spectral features and survive convolution; they are not a single
outlier. The diagnostic plot is
`truth_map/aerosol_chunked/state009_o2_components_before_after_grid_correction.png`.
The interpolated curve shown there is diagnostic only.

After production, `validate_regenerated_truth.jl` requires all 64 scenes to
carry exact-grid regenerated-O2 provenance, verifies finite A-band and CO2
arrays, and proves bit-identical A-band results across the four XCO2-only
states in each physical group. It is followed—not replaced—by
`validate_truth_forward_closure.jl`.

## Required retrieval-readiness gates

Run these gates after regenerating the A-band truth:

1. Validate every high-resolution scene for finite values, completion flags,
   state identity, spectroscopy paths, profile provenance, and wavelength
   arrays.
2. Confirm the exact A-band retained-grid invariant for every truth mode.
3. Run one aerosol-on SIF-off truth/retrieval/linearized closure through all
   three bands.
4. Run one aerosol-free SIF-off closure to quantify the deliberate five-
   versus-nine-stream difference.
5. Run the selective analytic-Jacobian regression.
6. Run the state/profile/source/surface mapping regression.
7. Run the instrument forward/Jacobian regression.
8. Regenerate all Mueller-processed, convolved, and resampled OCO radiances
   from the accepted truth files.
9. Validate all OCO radiances and regenerate their frozen noise covariances.
10. Require `ABSCO_CLOSURE.md` to report `Overall result: PASS`; no failed or
    stale readiness sentinel may remain.

Representative commands are:

```bash
# GPU truth/retrieval/linearized closure
CUDA_VISIBLE_DEVICES=1 CUDA_DEVICE=0 \
  julia --project=. RRS_XCO2/inversion/validate_truth_forward_closure.jl

# Analytic-Jacobian regressions; run from test/ because paths are relative
cd test
julia --project=. -e '
  using Test, vSmartMOM
  include("test_selective_jacobians.jl")
'

# Return to the repository root for retrieval-specific tests
cd ..
julia --project=. RRS_XCO2/inversion/test_forward_state_mapping.jl
julia --project=. RRS_XCO2/inversion/instrument/test_synthetic_oco2.jl
```

Any future change to spectroscopy, atmospheric-profile preparation, angular
quadrature, aerosol interpolation/truncation, surface coordinates, source
spectra, wavelength grids, analyzer coefficients, convolution, resampling, or
Jacobian state mapping must update this document and rerun the affected gate.
