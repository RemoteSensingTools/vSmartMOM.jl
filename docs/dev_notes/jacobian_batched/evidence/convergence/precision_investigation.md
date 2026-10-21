# Precision isolation: O2 retrieval case 035

Investigated 2026-09-07 against production code at `3d5bd106`. This extends the
matched-control end-to-end probe in [README.md](README.md). Both saved terminal
states of `state035_corrected_siffalse` are used; the O2 band has 2,735 solve
nodes and 934 instrument samples. All results below are CUDA results. Noise
normalization uses the original retrieval observation variance. All
noise-normalized comparisons below are **after the study's Gaussian instrument
convolution and detector sampling**, including the Jacobian predictions. The
[explicit convolution audit](convolution_accuracy.md) reconstructs those
measurements from the saved high-resolution spectra and compares each stage.

## Frozen boundary experiment

`precision_frozen.jl` extracts the production linearized driver into a separately
named, in-memory diagnostic method. It freezes each Fourier moment's complete
core optics and optical tangents, quadrature, angles, surface coefficients,
prepared solar/SIF payloads, and Fourier weights. Promoting Float32 payloads to
Float64 preserves their numerical values exactly. Demoting Float64 payloads
includes their representational rounding. The underlying production elemental,
doubling, adding, surface, and postprocessing routines are used unchanged.

This separates preparation precision from RT arithmetic without reimplementing
the MOM equations. The instrument operator and its coordinate grid stay fixed
within each preparation pair. The two numerical thresholds are held at the
Float32 campaign values. Every run uses Fourier orders 0–3 and doubling counts
`[4,4,5,6,6,4,4,5,6,7,8,9,9,10,9,6]`.

Four identity checks reproduce the public API's radiance and Jacobian arrays
exactly. `verify_frozen.py` additionally checks exact equality of the native
measurement/Jacobian arrays against the independently saved earlier Float32
and matched-floor Float64 probes. All eight runs pass these controls, finite
output checks, Fourier/doubling checks, and output hash verification.

### Absolute differences at fixed state

Maximum absolute differences, in observation noise standard deviations:

| Change | Reference terminal state | Optimized terminal state |
|---|---:|---:|
| RT Float32 → Float64, preparation frozen at Float32 | 0.047837 | 0.050631 |
| RT Float32 → Float64, preparation frozen at Float64 | 0.051716 | 0.050581 |
| Preparation and instrument grid Float32 → Float64, RT fixed at Float64 | 0.890241 | 0.890223 |
| Complete workflow Float32 → Float64 | 0.908749 | 0.906982 |

These maxima occur at potentially different samples and must not be added as
a scalar error budget. Changing RT precision accounts for a smaller spectral
difference than changing preparation and instrument coordinates. The largest
RT-only Jacobian-column relative L2 difference is 0.00181.

### Local derivative consistency

For the displacement between the same two saved states, compare `Δy` with
the initial state's `K Δx`:

| Preparation | RT | max \|Δy\| / noise | max \|K Δx\| / noise | max remainder / noise |
|---|---|---:|---:|---:|
| Float32 | Float32 | 0.017234 | 0.000424401 | 0.017277 |
| Float32 | Float64 | 0.000421523 | 0.000424393 | 0.000038476 |
| Float64 | Float32 | 0.013744 | 0.000424208 | 0.013935 |
| Float64 | Float64 | 0.000425881 | 0.000424197 | 0.000004795 |

Promoting only RT arithmetic reduces the Float32-preparation remainder about
449-fold. Retaining Float32 RT with Float64 preparation still produces a large
remainder. Thus RT precision is a major source of the local inconsistency in
this case. Preparation precision also matters: its remaining local remainder
at Float64 RT is about eight times the full-Float64 result. This experiment
does not yet distinguish elemental initialization, repeated doubling, adding,
or reconstruction as the responsible RT stage.

## Instrument-coordinate control

The live study's `scripts/common.jl::surface_basis_grids` constructs Float32 and
Float64 ranges independently. The package parameter and atmosphere types also
store spectral coordinates as `Vector{Vector{FT}}`. Changing `FT` therefore
changes sample locations as well as arithmetic.

`verify_frozen_instrument.jl` compares the two instrument coordinate mappings
on each **identical** high-resolution Stokes array. This diagnostic changes
spectral labels without recomputing the physical spectrum; it is not a
physically consistent alternative forward calculation. Its maximum change is
0.4359–0.4369 noise σ. Holding one of these instrument maps common while changing
preparation at Float64 RT leaves a maximum difference of 0.9963–0.9971 σ.
Preparation of the spectrum itself therefore also contributes substantially.

The physically consistent grid control, `precision_grid.jl`, rebuilds Float64
optics and sources at the exact promoted Float32 solve nodes and processes them
at those same nodes. The original grids differ by at most 0.000647083 cm⁻¹
(0.000038665 nm). Both grid-control solves preserve the same Fourier orders
and doubling counts.

| Comparison | Reference state, max noise σ | Optimized state, max noise σ |
|---|---:|---:|
| Change grid only within full Float64 preparation/RT | 0.890616 | 0.890617 |
| Preparation Float32 → Float64, grid matched and RT Float64 | 0.144835 | 0.144816 |
| Complete Float32 → Float64, grid matched | 0.158414 | 0.159316 |

Matching grid coordinates reduces the reference state's complete-workflow
maximum from 0.908749 to 0.158414 noise σ, and its RMS from 0.129320 to 0.026650.
The preparation-only RMS at fixed Float64 RT falls from 0.124684 to 0.014014.
The grid construction is therefore a major contributor to the earlier absolute
gap. Remaining preparation differences include sources, profile/input casts,
optical calculations, quadrature, and surface coefficients; they have not yet
been isolated individually.

The full Float64 calculation on matched nodes retains the small local
remainder: 0.000004799 noise σ. Changing the grid did not reintroduce the
Float32 RT sensitivity in this diagnostic.

## Scope and next decisions

The original three strict spectral-gate failures remain recorded. No gate,
truth data, priors, production kernels, or live study sources were changed.
Float64 is an internal precision control here, not an independent accuracy
reference. The subsequent [full three-band retrieval test](precision_retrieval_impact.md)
converges for this case with XCO2 shifts of −0.018396 ppm for native Float64
grids and +0.002494 ppm for exact promoted Float32 grids. Broader scenes/noise
realizations and individual preparation/RT stages still require investigation.

Future precision controls should distinguish spectral/input coordinates,
optical preparation, operator arithmetic, and numerical discretization.
Changing a single `FT` currently changes all four. Preserve supplied-tangent
boundaries and explicit numerical controls when extending the implementation;
new optical/source providers need precision and finite-difference checks at
that boundary as well as end-to-end retrieval checks. A mixed-precision policy
needs stage-specific accuracy and performance measurements before adoption.

The next bounded experiments are:

1. At the same grid and Float64 RT, swap prepared sources, then scalar/phase
   optics separately to locate the remaining 0.145 σ preparation difference.
2. Trace one frozen layer through elemental initialization and each doubling,
   then atmospheric adding, to locate the Float32 RT sensitivity. Test promoted
   stages only after establishing where the discrepancy grows.
3. Measure any candidate policy's cost and perturbation-size convergence,
   then rerun full retrievals with a declared, unchanged accuracy budget.

## Reproduction

Use the study data environment documented in the parent evidence directory
and run from `test/`, with `REPLAY_OUTPUT` pointing to the saved convergence
outputs. The scripts read the live adapter; they never modify it.

```bash
julia --project=. ../docs/dev_notes/jacobian_batched/evidence/convergence/precision_frozen.jl
python3 ../docs/dev_notes/jacobian_batched/evidence/convergence/verify_frozen.py "$REPLAY_OUTPUT"
julia --project=. ../docs/dev_notes/jacobian_batched/evidence/convergence/verify_frozen_instrument.jl "$REPLAY_OUTPUT"
julia --project=. ../docs/dev_notes/jacobian_batched/evidence/convergence/precision_grid.jl
python3 ../docs/dev_notes/jacobian_batched/evidence/convergence/verify_grid.py "$REPLAY_OUTPUT"
```

This machine used the Julia 1.12.6 binary directly because its launcher is
broken. Raw JLD2 files stay in `/tmp/vsmartmom-convergence-replay`; committed
JSON/TOML summaries, hashes, and logs record the evidence.
