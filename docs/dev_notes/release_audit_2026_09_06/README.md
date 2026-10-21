# Reproduction bundle

These standalone audit scripts accompany the [release evaluation](../release_readiness_2026-09-06.md) and [Jacobian profile](../jacobian_bottlenecks_2026-09-06.md) of commit `b524a0d85ae36ab0eecf8b23a939f4d6ce62c456`. They do not modify solver code and are not included in the default test suite.

Use the candidate's test environment, with `vSmartMOM` developed from its parent directory, and run from `test/`:

```bash
cd test
julia --project=. -e 'using Pkg; Pkg.develop(path=".."); Pkg.instantiate()'
julia --project=. ../docs/dev_notes/release_audit_2026_09_06/probes.jl
julia --project=. ../docs/dev_notes/release_audit_2026_09_06/extra_probes.jl
julia --project=. ../docs/dev_notes/release_audit_2026_09_06/coxmunk_parity.jl
julia --project=. ../docs/dev_notes/release_audit_2026_09_06/nnlib_compat.jl
```

The probes print measured outcomes and caught errors. The reports explain expected invariants and which outcomes demonstrate defects. `nnlib_compat.jl` must run in a fresh process to compare behavior before and after importing vSmartMOM.

`gpu_runner.jl` executes the five shipped CUDA test files using corrected absolute include paths and requires functional CUDA hardware. It does not relax assertions. `jacobian_profile.jl` takes workload/backend settings through environment variables; see its companion report for commands and limitations.

The full CPU suite was invoked with:

```bash
CUDA_VISIBLE_DEVICES='' PHASE1B_CPU=1 JULIA_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1 \
  julia --startup-file=no --project=. runtests.jl
```

The documentation environment was developed against the same candidate and built from `docs/` with `CI=false`, CUDA devices hidden and `julia --startup-file=no --project=. make.jl`.

Selected output and profile files are retained in [evidence](evidence/). Raw output may contain terminal progress-control characters and audit-machine paths. Host allocation totals are cumulative allocations, not peak RAM or VRAM. Sampled allocation profiles report only sampled bytes. Dependency versions resolved during this audit are preserved separately as `resolved-environments.json`; generated manifests remain outside the tracked release artifact.
