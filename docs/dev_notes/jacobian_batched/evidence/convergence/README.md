# Full retrieval convergence comparison

This isolated driver uses Sanghavi's unchanged `OptimalEstimation.jl` and
`VSmartMOMForward.jl` with the performance branch's elastic RT implementation.
It compares:

- **Reference:** `jacobian_basis=:physical`, `jacobian_adding=:matrix`, ordinary
  deep copies of trial parameters.
- **Optimized:** `jacobian_basis=:local`, `jacobian_adding=:source`, independent
  trial copies sharing read-only absorption LUT storage.

Both paths use the same current optical physics and precision fixes. This
isolates the propagation/copy optimizations from changes relative to older
archived package versions. The first pair includes compilation and should
not be interpreted as a warmed benchmark.

## Cases and reproducibility

Four cases replay archived observations: clear urban state 001 and aerosol
rural state 035, each corrected and uncorrected, perturbation 10, from
`bottom_layer_XCO2_retrievals/retrievals_acos_mapped_tapered_vertical_correlation_nosif`.
The driver reads the archived measurement, frozen measurement variance, full
prior covariance/state, first trial state, and stored OE stopping settings.
Fixed upper CO2 remains 400 ppm. No truth generation or campaign runner is
included; both inversions use elastic `noRS`.

Two additional cases add a controlled SIF signal to the state-035 observations.
At the archived final atmosphere, the direct physical reference computes the
measurement difference between zero SIF and `campaign_sif_state()` coordinates.
That difference is added to the existing observation. These are **synthetic
linear-SIF extension tests**, not replays of released full-template SIF/Raman
truth. The prior and frozen noise remain those of the original archive.

Run from `test/`, using the data environment in
[the original replay](../suniti_replay/README.md):

```bash
REPLAY_OUTPUT=/tmp/vsmartmom-convergence-replay julia --project=. \
  ../docs/dev_notes/jacobian_batched/evidence/convergence/replay.jl
python3 ../docs/dev_notes/jacobian_batched/evidence/convergence/verify.py \
  /tmp/vsmartmom-convergence-replay
```

Host: A100 GPU 0, Float32 IQU, 9 weighted streams, 16 layers, 5,011 solve
wavelengths, 2,742 measurements, 30 retrieval coordinates; Julia 1.12.6,
`JULIA_NUM_THREADS=4`, `OPENBLAS_NUM_THREADS=1`, `CUDA_VISIBLE_DEVICES=0`.
The host requires the actual Julia executable at
`/home/cfranken/.julia/juliaup/julia-1.12.6+0.x64.linux.gnu/bin/julia`.
The verifier uses Python 3.6's HDF5/NumPy and Python 3.11's TOML parser.

The driver saves each final state, measurement/Jacobian, posterior covariance,
averaging kernel, full trial-state/cost/acceptance history, and input arrays.
It checks the template state and all six LUT coefficient hashes after every
pair. `verify.py` independently reads HDF5 arrays, checks finite values and
identical inversion inputs, and verifies input-file hashes.

Before inspecting paired results, the comparison limits were set to:
0.01 ppm XCO2, 0.01 prior standard deviations for every state coordinate,
0.01 measurement-noise standard deviations for every sample, 0.1% final cost,
and 0.1% relative Frobenius difference of the prior-scaled posterior covariance.
Outcomes and full accepted/rejected decision sequences must agree. These
are implementation-comparison criteria; they do not establish the scientific
adequacy of the underlying retrieval.

The verifier separately reports drift from archived XCO2 and archived stopping
outcomes for the four original observations. Those differences include earlier
intentional optical/precision changes and must not be attributed to the new
copy policy or local/source propagation alone.

## Measured retrieval results

All twelve inversions converged with outcome 1 (all bands pass the study's
fit-quality criterion). Every pair has the same full accepted/rejected
sequence: three accepted trials for clear scenes; five trials, including
one rejection, for aerosol scenes.

The predeclared **full comparison gate does not pass**: three of six pairs
exceed the 0.01-noise spectral limit. All state, XCO2, cost, posterior, and
iteration-decision criteria pass. The limits were not relaxed after the run.

| Scene / observation / SIF | Absolute XCO2 difference (ppm) | Max spectral difference (noise σ) | Full gate |
|---|---:|---:|---|
| state001_corrected_siffalse | 3.75265e-06 | 0.000297489 | pass |
| state001_uncorrected_siffalse | 1.87626e-06 | 0.000246398 | pass |
| state035_corrected_siffalse | 0.000448989 | 0.0172382 | FAIL: spectral |
| state035_corrected_siftrue | 0.000110839 | 0.0207552 | FAIL: spectral |
| state035_uncorrected_siffalse | 5.44618e-05 | 0.00871695 | pass |
| state035_uncorrected_siftrue | 0.000142444 | 0.0133061 | FAIL: spectral |

Across the six pairs, maximum state-coordinate disagreement is 0.000236 prior
standard deviations, maximum relative cost change is 0.0108%, and maximum
relative change in the prior-scaled posterior covariance is 0.00243%.
The four original observations also retain their archived convergence outcomes,
trial counts, and rejection counts. Their reference XCO2 shifts from the
older archive by at most 0.001012 ppm.

`verification.json` contains the complete metrics and raw-array hashes.
`verify.py` intentionally exits **1** for this evidence because of the
unresolved spectral gate. Raw JLD2 arrays remain under
`/tmp/vsmartmom-convergence-replay`; compact results and logs are retained here.

The extension guards passed 29 focused checks; the existing selective
Jacobian/finite-difference tests passed 51 checks and source regressions passed
86 checks (166 distinct checks total). The strict documentation build passed.
