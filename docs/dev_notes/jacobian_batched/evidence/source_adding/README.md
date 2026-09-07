# Equivalent-source adding evidence

See [the implementation report](../../source_adding.md) for the algorithm,
supported configurations, results and limitations.

The two `source-real-*.toml` files contain warmed A100 timings, allocation
measurements, environment configuration, and matrix/source parity errors for
69-column IQU retrievals with real CO₂/CH₄/H₂O absorption. Each processes
10,000 wavelengths in ten independently prepared 1,000-point chunks. Their
scalar Lambertian benchmark path is unchanged by the final Legendre-only
direct-beam repair. The source fingerprints describe files on disk at benchmark
completion; final correctness tests are recorded separately.

`source-end-to-end-10000.toml` uses independent model construction inside the
timed boundary: 32.4561 s combined versus 28.1427 s forward-only (1.153×).
It records synchronized construction and solve subtimers. This full-rebuild
ratio has a different cost balance from repeated solves on prepared optics.

`source-streams16.toml` and `.log` contain the prepared-model 57×57-operator
comparison and CUDA/host profiles: 512 wavelengths, 20 layers, 69 columns,
29.3663 s combined versus 11.0652 s forward (2.654×). This larger scene checks
forward-radiance parity; it does not compare all Jacobian columns with matrix
adding at that size. Reproduce with the real-absorption command below, replacing
the source-comparison flag by `AUDIT_PROFILE_ONLY=true AUDIT_PROFILE_MODE=source`,
and setting `AUDIT_STREAMS=16 AUDIT_NSPEC=512 AUDIT_CHUNKS=1`.

`phase-storage-10000.toml` profiles 10,000 wavelengths in one batch, with
synthetic absorption and one aerosol. `phase_storage_benchmark.jl` records the
separate angular-node, interpolation, basis, and layer-mixture timings and
asserts that the local basis is shared across layers.

The CPU suite was run in contiguous segments while fixing failures:

- `source-adding-full-suite.log`: original full run, all testsets through aerosol
  reference passed. It stopped at three local-phase comparisons with cancellation
  residuals of approximately 2.4×10⁻¹⁵ in analytically zero gas columns. The
  Float64 absolute tolerance was adjusted from 10eps to 32eps; exact gas-phase
  zero checks were retained.
- `source-adding-suite-resumed.log`: the complete local-basis testset passed
  (1,634 checks), then the new solar-view source test exposed the Legendre
  surface's duplicate BOA direct carrier. This is a historical failing run,
  not final validation.
- `source-adding-suite-final.log`: all testsets from source adding through the
  final VLIDORT baseline passed (4,514 checks), including selected layouts,
  existing surface Jacobians, model updates, Aqua/JET and the quality gates.
  Combined with 2,353 passes before local optics in the original run and
  1,634 local-basis passes in the resumed run, this gives **8,501 CPU passes
  and 15 existing skips/broken tests** across the contiguous suite segments.

`source-adding-cuda-validated.log` records the final source-adding regression:
**92/92 passed** on CUDA. Coverage includes scalar and IQU Stokes layouts,
Float32/Float64, embedded/external solar, two aerosol modes, scalar/Legendre
surfaces, independent forward parity, gas finite differences, and zero-albedo
finite differences. `source-adding-docs-final.log` records the successful
strict local Documenter/Vitepress build (`CI=false`, no deployment).

All Julia tests were launched from `test/` with Julia 1.12.6, four Julia threads
and one BLAS thread. CPU validation disables CUDA visibility; GPU validation
uses GPU 0 with scalar CUDA indexing disabled. Metal was not hardware-tested.
Copied logs have terminal control codes and trailing whitespace removed.

Reproduce the real-absorption benchmark from `test/`:

```sh
CUDA_VISIBLE_DEVICES=0 JULIA_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1 \
  AUDIT_BACKEND=cuda AUDIT_COMPARE_SOURCE=true \
  AUDIT_EXTERNAL_SOLAR=false AUDIT_NSPEC=1000 AUDIT_CHUNKS=10 \
  AUDIT_LAYERS=20 AUDIT_GASES=CO2,CH4,H2O \
  AUDIT_NU_MIN=6150 AUDIT_NU_MAX=6250 AUDIT_OUTPUT=/tmp/source-real.toml \
  julia --project=. ../docs/dev_notes/jacobian_batched/local_basis_benchmark.jl
```

Repeat with `AUDIT_EXTERNAL_SOLAR=true` for the TOA-only row. Run the phase
profile with the same Julia/thread/GPU settings:

```sh
AUDIT_OUTPUT=/tmp/phase-storage.toml \
  julia --project=. ../docs/dev_notes/jacobian_batched/phase_storage_benchmark.jl
```

Add `AUDIT_END_TO_END=true` to time parameter parsing, independent forward or
LinMode model construction, and the complete solve. This mode compares source
adding with an independently built forward-only model and records synchronized
construction/solve subtimers. JIT compilation and artifact downloads are warmed;
Mie, spectroscopy and upstream optical derivatives are inside the timed call.


## Subsequent construction optimization

See [construction report](../../construction_cost.md) for the shared HITRAN
cache, Mie function boundary and warmed LUT/direct-HITRAN comparisons.
`construction-profile.log` is the baseline before these changes.
`mie-before.log` and `mie-after.log` isolate CPU NAI-2 values and derivatives;
`mie-parity.log` records bitwise comparisons with the implementation at
`8bc60941` for Float32/Float64, wavelengths 0.76/1.626 μm and imaginary indices
0/0.01, including all six Greek arrays, extinction, SSA and their tangents.

`lut-end-to-end.toml` and `.log` record the cached legacy LUT run and one-time
load metadata. This is the updated full-construction timing boundary with LUT
startup excluded, not the earlier repeated-HITRAN-parsing boundary.
`regressions.log` contains 161 successful targeted CPU checks (initial cache,
aerosol reference finite differences, source adding, H₂O self broadening).
`construction-cache-final.log` records the extended 13-check cache regression,
including BatchContext object reuse. `construction-docs.log` records the
successful strict local docs build. No deployment was performed.

`hitran-cached-end-to-end.toml` and `.log` record the corresponding direct-HITRAN
run with shared parsed data: 3.174 s forward / 6.473 s with Jacobians (2.039×).
Both current benchmark tables record source fingerprints and all three samples.
