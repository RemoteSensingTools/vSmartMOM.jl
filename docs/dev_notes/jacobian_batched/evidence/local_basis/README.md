# Local-basis evidence

See [the investigation](../../local_basis.md) for scene definitions, equations,
interpretation and reproduction commands. These are staged development
measurements on an A100 PCIe 40 GB, Julia 1.12.6, Float64 unless indicated.
Timing scripts warm each path and synchronize CUDA; TOML files preserve all
three samples, cumulative allocations and the maximum parity error.
CPU timings are not performance claims because the machine was shared.
Log copies normalize terminal line endings and trailing whitespace; numerical
values and test outcomes are unchanged.

| Files | What they establish |
|---|---|
| `local-final-{iqu,scalar}-10000.toml` | Five layers, absorption-free fourteen-column layout (five zero gas columns); physical vs local propagation after the elemental limit fix. |
| `local-real-iqu-10000.toml` | Twenty layers, 69 columns with 60 active gas columns; ten independently prepared spectral chunks. Precedes the exact-zero contraction guard and elemental limit fix. |
| `medium-streams{8,16}.toml` | Previous per-column fallback vs blocked products vs local basis, 512 wavelengths. |
| `vendor-streams{8,16}.toml` | Same scenes with solve-owned cuBLAS pointer batches. The blocked reference is untimed and has empty sample arrays. Includes source fingerprint and configuration. |
| `external-iqu-10000.toml` | Actual external-solar TOA solve, N=15; untiled/tiled physical propagation, local basis and forward, with source fingerprint. |
| `small-tile-scalar-10000.toml` | Embedded-solar scalar solve, N=6, after expanding the shared-memory gate; same-run untiled reference. |
| `local-core-*.log` | Supplied τ/ϖ/β₂ directions, black boundary, repeated spectral optics; CPU all-column endpoint finite differences and 10,000-point CUDA state-size scaling. |
| `source-response-*.log`, `layer-response-cuda512.log` | Experimental equivalent-source adding; CPU finite differences, CUDA parity at 512 points, then 10,000-point timings without allocating the dense comparison workspace. |
| `response-streams16-*.log` | Independent-source responses at sixteen streams: CPU finite differences, CUDA parity for 1/18 columns, then 1/60-column timings. |
| `source-sweep-*.log`, `response-modes-final-cpu.log` | Equivalent-source sweeps without a suffix cache; all-column CPU finite differences and CUDA comparisons. The mode-comparison run also checks the extracted source/matrix helper. |
| `elemental-boundary-before.log`, `elemental-boundary-cuda.log` | Independent forward-kernel finite differences at zero τ, ϖ and projected phase: 28 pre-fix failures; corrected CPU/CUDA checks pass. |
| `local-actual-external-*.log` | Corrected constructor-keyword fixtures: 1,634 CPU and 311 CUDA checks, including actual external-solar aerosol phase columns and TOA output. |
| `vendor-products.log`, `tile-crossover*.log` | Isolated product-rule crossover experiments, not complete-solve speedups. At 512 points N=6 is approximately neutral; larger tested operators benefit from shared memory. |
| `final-vendor-tiles.log`, `vendor-workspace-cpu.log` | Pointer ownership, active views, backend product parity and CPU workspace regression. |
| `final-local-docs-build.log` | Strict documentation build with deployment disabled, after the elemental limit correction. |
| `final-regression-suite.log` | Full CPU suite: 8,360 passed, 15 skipped/broken; GPU-specific cases are unavailable in this run. |
| `final-vendor-docs-build.log` | Strict documentation build after backend-workspace additions, deployment disabled. |
| `shared-core-fixture-cpu.log` | Parity and all-column endpoint finite differences after extracting the common core fixture. |
| `cached-inverse-{cpu,cuda}.log`, `cached-inverse-streams16.toml` | Rejected inverse-cache experiment: checks pass, complete-solve timing shows no gain. These are not the production implementation. |

The older nominal external-solar local-basis fixtures used an ignored YAML
entry. Their counts are not evidence of external-solar aerosol coverage;
the `local-actual-external-*` runs supersede them. External-solar runs expose
TOA only; embedded-solar runs compare both endpoints.

Allocation totals measure traffic, not peak resident memory. Historical files
without a source fingerprint are tied to the implementation stages above;
do not describe every file as a measurement of the final branch state.

The rejected inverse cache is preserved as `cached_inverse_experiment.patch`,
which applies to the production source at `4b5a2408`. In an isolated checkout
apply it with `git apply --unidiff-zero` and the patch path. With the benchmark
files available, run from `test/`:

```bash
VSMARTMOM_JACOBIAN_GPU_TEST=true julia --project=. test_jacobian_vendor.jl
AUDIT_BACKEND=cuda AUDIT_NSPEC=512 AUDIT_STREAMS=16 julia --project=. \
  ../docs/dev_notes/jacobian_batched/evidence/local_basis/cached_inverse_benchmark.jl
```

The reproduction driver writes separate reference/cached TOML files; the
recorded experiment used a physical reference with uncached inverses followed by
cached physical/local/forward measurements in one benchmark loop. All paths
were warmed. The later cached-sweep follow-up was terminated after rejection
of the production inverse-cache change and is excluded from this record.
