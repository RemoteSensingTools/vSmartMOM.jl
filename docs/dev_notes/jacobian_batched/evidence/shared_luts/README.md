# Sharing read-only absorption tables across retrieval trials

The package exposes `copy_parameters(params; share_luts=true)` for independent
trial-state copies that reuse loaded gas/H₂O table storage. Ordinary `deepcopy`
and the helper's default still copy the tables. Atmospheric state, surface and
aerosol parameters, and the LUT container lists remain independent; aliases
within the copied state are preserved. Shared tables, including coefficient
arrays, grids and mutable metadata, must remain read-only.

The implementation uses a private deepcopy wrapper and one memo for the whole
parameter graph. It marks mutable table storage as shared before traversing the
parameters. Marking only the outer LUT would be insufficient: Julia reconstructs
immutable LUT/interpolator wrappers and would otherwise still copy their arrays.
It does not walk coefficient-array elements, change Base's policy for existing
types, or cache state-dependent optical depths.

## Measured results

Same-process A100 GPU 0 comparison on 2026-09-07, Julia 1.12.6, four Julia
threads, one OpenBLAS thread; medians of three warmed evaluations:

| Copy policy | Complete evaluation | RT/Jacobian | Host allocations |
|---|---:|---:|---:|
| Deep-copy tables | 8.305 s | 5.041 s | 4.670 GB |
| Share read-only table storage | **5.368 s** | 4.873 s | **0.954 GB** |

This saves **35% of evaluation time** (1.55× faster) and **80% of host
allocations**. The isolated copy probe falls from 3.190 s / 3.716 GB to
0.328 ms / 44,832 bytes. RT is unchanged; the small difference between its
timer medians is run-time variability rather than a new solver optimization.

The original and perturbed states both produce **bitwise-identical measurements
and all 30 Jacobian columns** across copy policies. Returning from the perturbed
state to the original state reproduces its outputs bitwise. The parameter
template and all six loaded coefficient-array SHA-256 hashes remain unchanged.
All 154 study source/configuration hashes match the prior inspection before
and after the replay. The 53 portable checks and strict Documenter/Vitepress
build passed. Raw timings, array hashes, isolation results and logs are stored
alongside this README.

## Validation and reproduction

Run the portable checks from the performance worktree's `test/` directory:

```bash
julia --project=. test_copy_parameters.jl
```

The 53 checks cover shared coefficient storage, independent profiles/VMRs,
aerosols/surfaces/grids, independent LUT/H₂O lists, preserved state aliases,
unchanged ordinary deepcopy, absent tables and H₂O sentinels. Repeated model
construction and full source-adding RT produce exact radiance/Jacobian agreement
between default and shared copies while preserving the original coefficients.

`replay.jl` loads the unchanged study adapter, then defines two separately named
methods in memory. Both use the five-direction local basis and source adding;
only the shared method replaces `deepcopy(evaluator.base_parameters)` with
`copy_parameters(evaluator.base_parameters; share_luts=true)`. It includes no
campaign runner and writes only under `REPLAY_OUTPUT`.

Use the study data environment from [the original replay](../suniti_replay/README.md),
with `CUDA_VISIBLE_DEVICES=0`, four Julia threads and one OpenBLAS thread:

```bash
# From test/, with STUDY_ROOT, REPLAY_STATE and the data variables already set:
REPLAY_OUTPUT=/tmp/vsmartmom-shared-luts-replay julia --project=. \
  ../docs/dev_notes/jacobian_batched/evidence/shared_luts/replay.jl
python3 ../docs/dev_notes/jacobian_batched/evidence/shared_luts/verify.py \
  /tmp/vsmartmom-shared-luts-replay
```

On the benchmark host the direct Julia binary is
`/home/cfranken/.julia/juliaup/julia-1.12.6+0.x64.linux.gnu/bin/julia`; the shell's
Juliaup launcher is incomplete. Python 3.6 supplies HDF5/NumPy and Python 3.11
parses TOML, through the neighboring selected-basis verification helper.

The full A100 replay covers all 5,011 wavelengths and 30 instrument-level
Jacobian columns. Each copy policy is warmed once and timed three times with
CUDA synchronization and garbage collection between samples. The isolated copy
probe is also warmed and repeated three times. After timing, a perturbed state
changes pressure, CO₂, aerosol loading/height, surface and SIF before returning
to the original state. Both states must agree exactly between copy policies,
the repeated original state must reproduce its output, and the template and
SHA-256 hashes of all loaded coefficient arrays must remain unchanged.
The Python check independently compares the saved arrays bitwise.

`adapter.patch` is a concrete integration patch against the inspected live
study adapter. It enables both source adding and shared LUT copies and requires
the package version containing these APIs. `git apply --check` passed against
the inspected checkout; the patch was not applied to it. Full retrieval
convergence and active campaign migration remain separate validation steps.
