# Noise normalization, systematic differences, and retrieval impact

## What noise σ means

For detector sample `i`, `σᵢ = sqrt(Se_diagonal[i])` is the assumed standard
deviation of measurement noise in **convolved detector radiance**, in
`mW m⁻² sr⁻¹ nm⁻¹`. A reported spectral difference of 0.05 noise σ means
`|y₂ᵢ − y₁ᵢ| / σᵢ = 0.05`. It is 5% of that sample's one-standard-deviation
noise, not 5% of its radiance and not an XCO2 uncertainty. It is also distinct
from the Gaussian convolution kernel's wavelength width.

The archived case uses a diagonal measurement covariance, frozen throughout
retrieval. `verify_noise.py` confirms exact equality with the archived noise
product and reconstructs its photon-plus-background formula from the stored
signal, photon/background coefficients, dynamic-range constants, and photon
energy. No new noise estimate is fitted to the model differences.

For O2 case 035, σ ranges from **0.010861 to 0.115816**, with median **0.092924**
in the radiance units above. Median archived signal-to-noise ratio is **630.36**.
At that SNR, 0.05 noise σ corresponds to about **0.0079% of radiance**.

## Are the differences systematic?

They are deterministic numerical/model differences with a coherent spectral
component. For the reference terminal state, defining the difference as
Float64 minus Float32 after convolution:

| Comparison | Mean difference / noise σ | Negative samples | Adjacent-sample correlation |
|---|---:|---:|---:|
| RT arithmetic only, Float32 optical inputs fixed | −0.019264 | 98.07% | 0.9025 |
| Full precision change, spectral grid matched | −0.019750 | 95.40% | 0.8361 |
| Full precision change, native grids | −0.039020 | 77.62% | 0.8277 |

The nearby optimized terminal state gives essentially the same pattern.
These differences should not be modeled as independent zero-mean measurement
noise. A coherent pattern can shift retrieved parameters even when each
individual difference is smaller than one detector noise standard deviation.
The shift depends on its projection onto the Jacobians and on the prior;
the maximum spectral difference alone cannot determine it.

With the Jacobian held fixed, a useful local diagnostic for switching the
forward model while keeping the observation fixed is

```math
\delta x \approx -G\,\delta y,\qquad
G=(K^T S_e^{-1}K+S_a^{-1})^{-1}K^T S_e^{-1}.
```

This approximation does not include changed Jacobians, nonlinear trajectory,
or finite stopping effects. The full inversions below measure their combined
effect under the campaign's actual stopping rule.

## Measured full-inversion shifts

`precision_retrieval.jl` runs all three bands for the corrected, no-SIF case 035
and archived noise realization 10. It holds the observation, noise covariance,
prior, initial state, surface/source parameterization, and stopping rule fixed.
Both Float64 runs retain the Float32 elemental floor and threshold. The second
Float64 run also uses the exact promoted Float32 solve nodes **in all bands**.
Every model evaluation includes the full instrument convolution and sampling.

All three runs converge with the same five-trial decision sequence:
**accepted, rejected, accepted, accepted, accepted**. The stopping threshold
is the archived `d_sigma_sq_scaled < 2`; it was not tightened or relaxed. The
Float32 replay reproduces the earlier optimized retrieval's states, costs,
measurement, Jacobian, and posterior exactly.

| Configuration | Retrieved XCO2 (ppm) | Change from Float32 (ppm) | Surface-pressure change (hPa) | Largest active-layer CO2 change (ppm) |
|---|---:|---:|---:|---:|
| Float32 baseline | 402.826427 | 0 | 0 | 0 |
| Float64, native grids | 402.808031 | **−0.018396** | −0.067639 | 0.068088 |
| Float64, exact Float32 grids | 402.828921 | **+0.002494** | −0.003739 | 0.012073 |

All columns are evaluated with one common Float64 dry-air-column diagnostic,
including each retrieved surface pressure, so reporting precision does not
contaminate the comparison. The baseline's original Float32 diagnostic reports
402.8264255 ppm; changing only that diagnostic changes it by 0.0000011 ppm.

The baseline posterior XCO2 uncertainty is **0.314147 ppm**. The measured
precision-induced shifts are approximately **0.0586** and **0.00794** times
that posterior standard deviation. Both are small for this case, but they
are systematic shifts, not reductions in random measurement noise.
These are differences between numerical workflows, not estimates of either
workflow's absolute bias against truth or an ensemble-average bias.

The new complete solves took **34.56 s** for Float32, **40.80 s** for Float64
on native grids, and **40.36 s** for Float64 on matched grids. Each configuration
was evaluated at the saved state before timing its inversion; these are single
inversion timings, not medians of repeated performance trials. Even the full
Float64 runs are about **17×** shorter than the archived 698.15 s evaluation
total for this case, subject to the archived-versus-isolated timing caveat in
the [study speed summary](../../suniti_inversions.md).

The earlier maximum difference of **0.000449 ppm** concerned physical/matrix
versus local/source Jacobian propagation, with both using Float32. It must
not be used as the Float32-versus-Float64 precision effect.

## Relation to the isolated O2 diagnostics

The local gain predicts an XCO2 change of **+0.0000427 ppm** for the isolated
O2 RT-only precision difference, with the other two bands' model differences
set to zero. This is a linear estimate; an isolated RT-only full inversion
was not run. It must not be confused with promoting preparation and RT in
all three bands.

For the complete three-band fixed-state differences, the local predictions
are **−0.017877 ppm** on native grids and **+0.003028 ppm** on matched grids.
The full inversions give −0.018396 and +0.002494 ppm respectively. This
comparison illustrates why signed, Jacobian-projected spectral differences
are more informative about retrieval shifts than the largest pixel residual.

## Grid locations versus fitted wavelength shifts

This 30-coordinate retrieval has no wavelength-shift or stretch parameter.
Its active state contains pressure, twelve CO2 layers, aerosol loadings/heights,
surface polynomial coefficients, and SIF. The independently constructed
Float32 and Float64 grids therefore differ before any fitted spectral shift.

`verify_grid_shift.jl` tests whether the discrepancy could be absorbed by one.
The O2 grids have a mean wavelength offset of 0.000015007 nm, but an affine
coordinate fit leaves 0.000016164 nm RMS of irregular pointwise displacement.
The maximum original displacement is 0.000038665 nm. Thus there is both an
offset and nonuniform quantization of the sample coordinates.

An offline weighted fit of the convolved spectral residual to a wavelength
shift removes only **0.119% of its squared noise-weighted norm** for the
grid-only Float64 comparison. Adding a stretch removes only **0.134%**.
RMS falls from 0.124020 to 0.123946 noise σ for a shift alone; the fitted shift
is about 0.000000699 nm. The derivative used for this projection is checked
by halving its finite-difference step. This is not a new retrieval with
calibration parameters, and those parameters remain absent from the study.

Changing where a spectrum is sampled is different from shifting all its
physical features. Uneven sample locations and integration of narrow features
cannot generally be repaired by a global shift. A promising implementation
direction is one canonical Float64 spectral coordinate grid, shared by optical
lookup, solar/source evaluation, and convolution, with operator arithmetic
precision controlled separately. Merely relabeling already computed Float32
spectra as Float64 cannot restore lost coordinate information. Matching the
current Float32 nodes demonstrates consistency, not absolute accuracy of those
nodes; a canonical-grid change still needs discretization and retrieval tests.

## Validation and scope

`verify_precision_retrieval.py` checks identical observations, covariance,
priors, initial states, exact Float32 replay, final costs, state differences,
and artifact hashes. `noise-verified.json` records the noise provenance and
signed-residual statistics. All 154 tracked live study source/config hashes
remain unchanged.

This is one scene and one noise realization, using the campaign's convergence
criterion. It establishes the numerical impact for that case, not a bound
across all scenes or a fully minimized mathematical optimum. The original
Jacobian-implementation comparison's spectral gate failures remain recorded.
Future precision validation should include convolved maximum/RMS residuals,
their signed structure, state-space projections, and actual paired retrievals.
Broader scenes and noise realizations are needed to characterize an ensemble
bias, and independent references or discretization convergence are still
needed for absolute forward accuracy.

From `test/`, use the same study/GPU environment as `replay.jl`:

```bash
julia --project=. ../docs/dev_notes/jacobian_batched/evidence/convergence/precision_retrieval.jl
python3 ../docs/dev_notes/jacobian_batched/evidence/convergence/verify_precision_retrieval.py "$REPLAY_OUTPUT"
python3 ../docs/dev_notes/jacobian_batched/evidence/convergence/verify_noise.py "$REPLAY_OUTPUT" "$STUDY_ROOT"
```

Raw JLD2 files remain in `/tmp/vsmartmom-convergence-replay`. Committed
`precision-retrieval.toml`, `precision-retrieval-verified.json`, and
`precision-retrieval.log` record the outcomes and hashes.
