# Retrieval-selected local basis: validation and replay

This follow-up removes unselected microphysical phase directions before local
MOM workspace allocation. The three fixed-microphysics aerosols in the inspected
Sanghavi retrieval use five local directions instead of seventeen. See the
[local-basis derivation](../../local_basis.md) and
[retrieval audit](../../suniti_inversions.md).

The fresh baseline is commit `6cbaccc1`, run from the isolated worktree
`/tmp/vsmartmom-jacobian-basis-baseline`. The modified package is in
`/home/cfranken/code/gitHub/vSmartMOM-jacobian-perf`. Both use the same test
dependency manifest, Julia 1.12.6, A100 PCIe 40 GB GPU 0, four Julia threads,
and one OpenBLAS thread. GPU runs are sequential. All study source/configuration
hashes in [the earlier manifest](../suniti_replay/study-source-sha256.json)
were checked against the live checkout before and after the comparison.

## Replay results

Medians of three warmed complete evaluations:

| Configuration | Evaluation | RT/Jacobian |
|---|---:|---:|
| Baseline matrix | 15.911 s | 12.293 s |
| Selected matrix | 11.254 s | 7.977 s |
| Baseline source | 14.081 s | 11.273 s |
| Selected source | 8.194 s | 4.905 s |

Source adding is 1.72× faster end to end and 2.30× faster in RT. Its complete
measurement and 30-column Jacobian arrays are bitwise identical before and
after pruning. Automatic matrix adding preserves the measurement bitwise and
changes its most affected Jacobian column by 1.03e-5 relative L2; O₂ switches
from physical columns to the local basis. The selected source/matrix difference
is below 5.0e-5 detector-noise units and 5.83e-6 per-column relative L2.

`baseline-timings.toml` and `selected-timings.toml` contain all samples;
`verification.json` contains independent comparisons and output hashes.
`baseline-replay.log` and `selected-replay.log` include warm-up and copy probes.
The prior archive discrepancy is unchanged by pruning; it remains a separate
consequence of the preceding optical-precision repair.

## Regression results

- CPU: 5,164 local optical checks, 51 retrieval-selection checks, and 376
  solar/SIF source-adding checks: **5,591 passes**.
- CUDA with scalar indexing disabled: 2,187 local optical checks and 376
  solar/SIF source-adding checks: **2,563 passes**.
- The strict Documenter/Vitepress build passed. Metal was not hardware-tested.

The optical checks require identical forward optical arrays in physical and
local coordinates, compare contracted phase/scalar derivatives, assert exact
gas-phase zeros, and compare complete matrix/source RT outputs. Coverage includes
fixed and partially selected microphysics, Float32/Float64, and embedded/external solar;
independent finite differences cover gas propagation and the existing selected
pressure/aerosol/surface/SIF parameter classes. Logs are retained alongside this
README.

## Reproduce

Use the data environment variables from the
[original replay instructions](../suniti_replay/README.md), and run from each
worktree's `test/` directory. Set `REPLAY_FAST=true` and use distinct output
directories:

```bash
# Baseline worktree, commit 6cbaccc1:
REPLAY_OUTPUT=/tmp/vsmartmom-basis-baseline-replay julia --project=. \
  ../docs/dev_notes/jacobian_batched/evidence/suniti_replay/replay.jl

# Modified worktree:
REPLAY_OUTPUT=/tmp/vsmartmom-pruned-basis-replay julia --project=. \
  ../docs/dev_notes/jacobian_batched/evidence/suniti_replay/replay.jl
```

On this host the shell's `julia` launcher points to an incomplete Juliaup backup;
the runs used `/home/cfranken/.julia/juliaup/julia-1.12.6+0.x64.linux.gnu/bin/julia`
directly. No Juliaup installation was changed.

`matrix` means automatic basis selection with matrix adding. At the baseline
it uses physical O₂ columns and a 17-direction local basis in CO₂; the modified
version uses a five-direction local basis in every band. `source` explicitly
selects local/source adding in every band, including O₂ SIF. The physical
retrieval coordinates, state, spectral grids, forward optics and instrument
operator are preserved. Timings include parameter deepcopy and exclude one
complete warm-up per mode; each median has three synchronized samples.

Portable CPU checks, from `test/`:

```bash
julia --project=. -e 'include("test_local_jacobian.jl")'
julia --project=. -e 'include("test_selective_jacobians.jl"); include("test_source_adding_sif.jl"); include("test_source_adding.jl")'
```

The GPU driver reuses the optical-check definitions without repeating the CPU
suite, adds a gas-only selection with no pressure/aerosol columns, and runs
the SIF/solar source-adding regressions:

```bash
CUDA_VISIBLE_DEVICES=0 VSMARTMOM_SOURCE_GPU_TEST=true julia --project=. \
  ../docs/dev_notes/jacobian_batched/evidence/selected_basis/cuda-tests.jl
```

Independent verification reads the JLD2/NetCDF arrays using HDF5 and recomputes
per-column norms using NumPy. It checks every replay state against the archive,
finite output arrays, measurement differences below 0.001 detector-noise units,
and each Jacobian-column relative L2 difference below 3e-4. These checks apply
between current implementations and before/after pruning. Comparison with the
historical archive is reported separately because the preceding precision
repair intentionally changed the old mixed Float32/Float64 calculation.

```bash
python3 ../docs/dev_notes/jacobian_batched/evidence/selected_basis/verify.py \
  /tmp/vsmartmom-basis-baseline-replay /tmp/vsmartmom-pruned-basis-replay
```

The verification script uses the host's Python 3.6 HDF5/NumPy installation
and invokes Python 3.11 for TOML parsing. It records hashes of the raw JLD2
outputs and timing files. This is a complete single-state forward/Jacobian
replay, not a rerun of inversion convergence or an active retrieval campaign.
