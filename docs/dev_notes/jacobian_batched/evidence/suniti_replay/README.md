# Isolated RRS_XCO2 replay evidence

See [the audit](../../suniti_inversions.md) for study configuration, timing
boundaries, and limits. This benchmark requires the study checkout and its
external absorption/solar data; it is not part of the portable test suite.
It only includes the optimal-estimation definitions and forward adapter, never
the campaign runners. The adapter's solver call is replaced in a separately
named in-memory method for the source-adding comparison. The original method
and study files remain unchanged.

From this optimization worktree's `test/` directory:

```bash
export CUDA_VISIBLE_DEVICES=0 CUDA_DEVICE=0
export JULIA_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1
export STUDY_ROOT=/home/sanghavi/code/github/uni_vSmartMOM/RRS_XCO2
export VSMARTMOM_HITRAN_LUT_DIR=/home/sanghavi/data/HITRAN_LUTs
export VSMARTMOM_SIF_DATA_DIR=/home/sanghavi/code/github/uni_vSmartMOM/src/SIF_emission
export SOLAR_OUT=/home/sanghavi/Raman_misc/workflows/worktrees/uni_vSmartMOM_sanghavi_2025-03-18/src/SolarModel/solar.out
export REPLAY_STATE="$STUDY_ROOT/bottom_layer_XCO2_retrievals/retrievals_acos_mapped_tapered_vertical_correlation_nosif/corrected/retrieval_state035_perturbation10.nc"
export REPLAY_OUTPUT=/tmp/vsmartmom-suniti-replay
export REPLAY_FAST=true
julia --project=. ../docs/dev_notes/jacobian_batched/evidence/suniti_replay/replay.jl
```

`matrix` uses the unmodified adapter with matrix adding and automatic basis
selection: physical columns for O₂, a local basis for the two CO₂ bands.
`source` forces a local optical basis and equivalent-source adding for **all three
bands**, including SIF in O₂. The original pre-SIF replay used physical matrix
adding for O₂ and local/source adding for the two CO₂ bands; its timings are
retained separately as `before-sif-timings.toml`. The older files label the
default matrix run `physical` and the alternative `fast`; `physical` there did
not force physical coordinates in the CO₂ bands. These mode names describe the
RT implementation; both compute a complete forward model and Jacobian.

The replay outputs measurements and Jacobians in JLD2 and per-trial timings in
TOML. It reports absolute differences from the archived measurement/Jacobian
and measurement differences in detector-noise units. Large absolute Jacobian
entries use heterogeneous coordinate units; compare column-relative norms as
well. Compilation is warmed once per mode. Each of three samples includes
parameter deepcopy, construction, all bands, and instrument processing.

`study-source-sha256.json` records 154 source/configuration hashes relative to
`/home/sanghavi/code/github/uni_vSmartMOM`, including the locally modified live
adapter. These were verified unchanged during the audit. The recorded source
version and the benchmark environment matter: the archived campaign used a
different dependency environment, so its terminal timing is historical context,
not a fresh controlled comparison.

Portable regressions, also run from `test/`:

```bash
julia --project=. -e 'include("test_source_adding_sif.jl"); include("test_source_adding.jl")'
CUDA_VISIBLE_DEVICES=0 VSMARTMOM_SOURCE_GPU_TEST=true julia --project=. -e 'include("test_source_adding_sif.jl"); include("test_source_adding.jl")'
```

The SIF regression covers Float32/Float64, scalar/IQU, embedded/external solar,
scalar/Legendre surfaces, zero/nonzero retrievable SIF, prescribed SIF, and
multiple sources with distinct SIF coordinate systems. It compares complete
matrix/source Jacobians and independent forward radiances, with central finite
differences for emission amplitude/slope, gas absorption, and albedo. The
separate preexisting regression exercises solar-only behavior.

`cpu-sif-tests.log` and `cuda-tests.log` contain the final 284-check SIF runs.
The latter also contains the 92 solar-only CUDA checks. `cpu-solar-regression.log`
contains the 92 solar-only CPU checks and an earlier 258-check SIF run, before
the extra albedo and multiple-source tests were added. `docs-build.log` records
the successful strict Documenter/Vitepress build for the SIF extension.
The subsequent optical-precision checks and study replay are recorded separately;
see the audit for the mixture and absorption-phase repairs they identified.

`final-cpu-tests.log` and `final-cuda-tests.log` cover the final combined runs:
4,100 CPU passes and 1,397 CUDA passes. `final-docs-build.log` records the strict
build after the repairs. The CUDA driver includes the optical-check definitions
without repeating the CPU suite, then runs the three CUDA optical cases and
the SIF/solar source suites. `phase-roundoff.log` records the cancelled entries
used to set the per-scatterer Float32 absolute allowance; the forward optical
arrays and gas phase zeros are still checked exactly.

`campaign-trials.csv` contains the 289 recorded trial timings from 64 archived
retrievals, with filenames relative to `STUDY_ROOT`; `campaign-summary.json`
aggregates these by aerosol condition. These include cold campaign calls and
different worker conditions and should not be treated as controlled timings.
