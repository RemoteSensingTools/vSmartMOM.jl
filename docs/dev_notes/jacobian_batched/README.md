# Batched Jacobian propagation

Latest implementation: [equivalent-source adding and spectral phase costs](source_adding.md).
Constructor optimization: [spectroscopy reuse and Mie derivative specialization](construction_cost.md).
Local optical basis: [factorization and remaining adding cost](local_basis.md).
Previous investigation: [IQU bottlenecks and core-optics cost](iqu_followup.md).
Scientific review: [core derivatives, truncation and state-vector scaling](paper_review.md).
Historical implementation: [Fortran ordering and ideas retained for Julia](fortran_ordering.md).

Development branch: `perf/jacobian-batched-propagation`, based on integration
commit `b524a0d85ae36ab0eecf8b23a939f4d6ce62c456`.

The implementation applies the forward solver's recent performance approach to
physical-parameter tangents: reduce launches, keep operands on the device, and
reuse scratch. It changes doubling and `ScatteringInterface_11` interaction,
which dominated the audited Jacobian workload. The analytic equations and public
Jacobian layouts stay the same.

## Design

`src/CoreRT/CoreKernel/jacobian_batched.jl` provides two operations:

- `_jmul!`: an in-place product where either operand can carry the parameter axis.
- `_jprod!`: the fused product rule `dA*B + A*dB` for every requested column.

On GPU, wavelength and parameter form one launch axis. A forward 3D matrix is
shared directly across parameters; a derivative operand stays 4D. There are no
per-column CuArray views to materialize and no duplicated forward matrices or
host-generated pointer arrays for these products. For CUDA square operators of
size 6–32 and at least 512 wavelengths, a workgroup stages operands in shared
memory for reuse across all output entries. Sizes 33–64 use solve-owned cuBLAS
pointer batches generated on the device, with portable blocked kernels as the
layout fallback. Smaller operators and source-vector products use the simpler
element kernel. CPU uses in-place BLAS slices, with a serial path for
small work and a coarser wavelength × parameter partition for larger work.

The propagation workspace lives with `AddedLayerLin` and is reused across layers
and Fourier moments. Doubling saves old sources, matrices and tangents until all
terms consuming them are complete. General interaction computes both directions
before committing the new composite; the two directions share one implementation
of the product-rule algebra. The forward inverse is computed once per direction
and reused for every parameter.

The path is enabled for Float32/Float64 CPU arrays and GPU operators up to
32×32. CUDA additionally uses 16×16 blocked products for operators 33–64
with at least 512 wavelengths, validated in the 8/16-stream IQU follow-up.
Larger operators, smaller batches above 32, other numeric types and manually
constructed layers without the workspace use the existing implementation.
These are conservative measured coverage limits, not universal crossover
points. Metal uses the portable kernel machinery but was not hardware-tested here.

`CoreRT._BATCHED_JACOBIANS_ENABLED[]` is a private diagnostic switch for A/B tests.
Do not change it concurrently with RT solves. The reference equations remain in
their original files. Special interaction cases 00/01/10 retain their existing implementations. Phase-tangent construction now shares
angular tables and uses static Stokes-block accumulators; zero tangent buffers
are allocated directly on their backend.

The workspace trades persistent scratch for much lower cumulative allocation.
The reference A/B path also constructs this scratch, so measurements compare
propagation strategies within the new branch. At the 64-point scalar case this
adds roughly 3.7 MB to the reference solve versus the original audit; it is small
relative to its 1.55 GB of cumulative allocation.

The 10,000-wavelength evaluation exposed another launch bottleneck in the
embedded-solar scalar Lambertian source tangents. The original surface builder
loops over wavelength and atmospheric parameter, performing a separate
matrix-vector product per pair, plus one albedo product per wavelength. For
10,000 wavelengths and fourteen total columns, that is 140,000 tiny products
at m=0. The new path flattens wavelength × atmospheric parameter into a single
matrix multiplication and evaluates all albedo columns with a second product.
It also keeps the incident work arrays in `FT` rather than accidentally using
Float64 for Float32 scenes. The private A/B switch covers this surface change
as well; external-solar surface construction keeps its existing implementation.

## Results after surface batching (before the IQU follow-up)

With both propagation and Lambertian source batching, commit `8aa6079e`:

| Scene | Reference full | Final full | Speedup | Reference selected | Final selected | Speedup |
|---|---:|---:|---:|---:|---:|---:|
| CUDA, scalar, 10,000 wavelengths | 16.716 s | 1.816 s | 9.20× | 4.263 s | 0.513 s | 8.31× |
| CUDA, IQU, 10,000 wavelengths | 38.268 s | 16.441 s | 2.33× | 8.164 s | 4.119 s | 1.98× |

These use the same Float64, five-layer, absorption-free scenes and hardware
specified below. Final timings are three warmed medians; reference medians are
reused from the first-stage runs. Both final runs also evaluate a fresh reference
solve and assert full R/T/dR/dT parity and selected dR column parity. Maximum
absolute Jacobian disagreement is `4.44e-16`; radiance disagreement is at most
`1.39e-17`. These checks preserve the reference implementation's behavior; they
do not settle the separate scientific issues in the release audit.

Full cumulative device allocation falls from 223.74 to **14.40 GB** for scalar
and from 1,911.63 to **121.94 GB** for IQU (15.5–15.7× reductions). Sampled
full-phase process memory, including CUDA's reserved pool and warmup, is
**25.2 → 8.2 GiB** for scalar and **39.5 → 32.3 GiB** for IQU. These sampled
pool peaks depend on allocation history and are not exact live-array peaks.
Both cases completed on the 40 GB A100 without spectral chunking.

Final full solves still allocate **3.13 GB** on the host for scalar and
**25.88 GB** for IQU. The preceding profile identifies phase-tangent construction
and mixing as material host work. Reducing those allocations and moving more of
that assembly to the device is the next candidate; it needs its own profile and
validation, especially for larger polarized operators. Current warmed forward
baselines are 0.136 s and 0.636 s, respectively, so Jacobians remain substantially
more expensive than a forward solve.

## Propagation-only measurements

These first-stage measurements are from commit `7cb8ae98`, before batching the
Lambertian source tangents. Final 10,000-point results are recorded separately
below; the first-stage table and profile explain how the second bottleneck was
found.

Median of three warmed, synchronized solves on an AMD EPYC 7H12 / NVIDIA A100
PCIe 40 GB, Julia 1.12.6, CUDA.jl 5.11.3, Float64, BLAS one thread. CPU measurements
use one Julia thread; CUDA uses four host threads. The server was shared with
other work. These are workload-specific measurements, not dedicated-machine
performance guarantees.

All scenes use five layers, one aerosol, Lambertian surface, embedded solar and
no gas absorption. There are fourteen native Jacobian columns, including five
unused gas slots retained by the legacy layout. The selected case requests four
columns: pressure, aerosol optical depth, aerosol profile location and albedo.
The IQU case has an 18×18 operator; scalar cases have 6×6 operators.

| Backend / scene | Reference full | Batched full | Speedup | Reference selected | Batched selected |
|---|---:|---:|---:|---:|---:|
| CPU, scalar, 64 wavelengths | 1.530 s | 1.034 s | 1.48× | 0.410 s | 0.266 s |
| CPU, IQU, 64 wavelengths | 7.697 s | 3.932 s | 1.96× | 1.451 s | 1.036 s |
| CUDA, scalar, 64 wavelengths | 2.918 s | 0.429 s | 6.81× | 0.811 s | 0.250 s |
| CUDA, scalar, 512 wavelengths | 3.565 s | 0.763 s | 4.67× | 0.893 s | 0.342 s |
| CUDA, IQU, 64 wavelengths | 3.330 s | 0.593 s | 5.61× | 0.834 s | 0.256 s |
| CUDA, scalar, 10,000 wavelengths | 16.716 s | 10.462 s | 1.60× | 4.263 s | 2.735 s |
| CUDA, IQU, 10,000 wavelengths | 38.268 s | 25.768 s | 1.49× | 8.164 s | 6.488 s |

Full CPU cumulative allocations fall from 1.553 GB to 0.181 GB in the scalar
case, and from 12.592 GB to 1.118 GB in the IQU case. Host allocation counters
are not peak RAM or device-memory measurements.

At 10,000 wavelengths, the scalar full solve's cumulative GPU allocation
falls from **223.74 GB to 14.41 GB**; the selected solve falls from **59.21 GB
to 3.42 GB**. The full-phase sampled process memory peak falls from **25.2 GiB
to 7.2 GiB**, including the reserved pool and warmup. This is a sampled peak,
not exact live memory. Separate `CUDA.@timed` allocation counters measure
cumulative device allocation; they do not imply that hundreds of GB are live
simultaneously. The 10,000-point runs explicitly reclaim the pool before each
phase and then warm the exact path before timing.

Reference/batched radiances agree exactly in the CPU benchmarks. Maximum absolute
differences across CUDA benchmark radiances are about `1e-17`; derivative
differences across the measured cases are at most `4.5e-16`.

The CUDA 64-point full-solve CUPTI profile records:

| Activity | Original audit | Batched propagation |
|---|---:|---:|
| Kernel launch calls | 72,002 | 14,672 |
| Host-to-device copy calls | 114,192 | 5,992 |
| Device-to-device copy calls | 4,676 | 386 |
| Dominant small GEMM kernel calls | 35,215 | 1,075 |

The original profile used the same scene on the integration commit. The batched
profile is retained with this note. Counts measure calls, not transfer volume.
The instrumented trace took 0.795 s versus 6.03 s in the original audit; use the
uninstrumented timing table for speedups because profiling changes runtime.

## Correct forward baseline

A separate commit fixes the angular cache's use of `length(β)` for matrix-valued
Greek coefficients. The angular axis is `size(β,1)`; wavelength count must not
inflate the angular tables. The new test constructs a spectral aerosol model,
checks the actual table extent, and checks bitwise forward parity with the cache
disabled. Both A/B propagation modes in the timing table include this fix.

The forward/Jacobian cost is still workload-dependent. Reusing inverses does not
make all requested physical derivatives free. The public concepts page and agent
onboarding now describe the algorithm without an unconditional `<2× forward`
promise; the inverse-tangent formula's erroneous minus sign in that explanation
is also corrected.

## Validation

- 1,536 additional surface comparisons passed on each of CPU and CUDA,
  covering Float32/Float64, scalar/IQU, one/multiple wavelengths, dark and
  spectral illumination, zero/nonzero albedo, zero/three atmospheric columns,
  an unused state-vector slot, and source clearing between m=0 and m=1.
- 416 CPU operator comparisons passed, covering Float32/Float64, scalar and IQU,
  one/multiple wavelengths, zero/multiple doublings, an inactive parameter prefix,
  and general interaction with nonsymmetric matrices and nonzero tangents.
- The same 416 comparisons passed on CUDA with scalar indexing disabled.
- Existing end-to-end checks passed: 36 Jacobian finite-difference assertions,
  42 selective Jacobian assertions, 73 external-solar assertions, and 227
  multisensor assertions.
- After the final surface storage-layout adjustment, the 36 Jacobian
  finite-difference, 42 selective-column, and 7 Float32 assertions were rerun
  and passed. The broader surface-change regression run also passed the 73
  external-solar and 86 source-routing assertions.
- Angular-table tests passed, including the five new spectral-cache assertions.
- All **20,291** assertions in the broader CUDA suite passed through the audit
  harness, which corrects the pre-existing test-runner include paths.
- Before the final surface optimization, the complete CPU suite passed **4,107** assertions across 58 top-level
  testsets, with the same 14 Broken/skipped entries as the integration audit
  (process exit 0). The optional Raman reference comparison was enabled. JET
  continues to exceed its advisory baseline, as it did on the integration commit.

This work preserves the reference propagation's numerical behavior. It does not
resolve the separate aerosol reference-normalization or Cox–Munk correctness
findings from the [release audit](../release_readiness_2026-09-06.md). Those still need their own physics corrections
and regression coverage before release.

Selected raw outputs and the CUPTI summary are retained in [evidence](evidence/),
and timing medians in [results.json](results.json).

## Reproduce

Run from `test/` after instantiating its environment against this checkout:

```bash
julia --project=. test_jacobian_batched.jl
julia --project=. test_lambertian_jacobian_batched.jl
VSMARTMOM_JACOBIAN_GPU_TEST=true julia --project=. test_jacobian_batched.jl
VSMARTMOM_JACOBIAN_GPU_TEST=true julia --project=. test_lambertian_jacobian_batched.jl
julia --project=. ../docs/dev_notes/jacobian_batched/gpu_runner.jl

AUDIT_BACKEND=cuda AUDIT_NSPEC=64 JULIA_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1 \
  julia --project=. ../docs/dev_notes/jacobian_batched/benchmark.jl
```

Use `AUDIT_BACKEND=cpu` and `JULIA_NUM_THREADS=1` for CPU results,
`AUDIT_NSPEC=10000` for the large batch, `AUDIT_POL='Stokes_IQU()'` for polarization,
`AUDIT_TIMERS=true` for the section timings, and `AUDIT_CUPTI=true` to collect
the separate CUDA profile. `AUDIT_TIMER_ONLY=true` runs just a warmup and one
batched full solve with section timers. Add `AUDIT_HOST_PROFILE=true` to sample
that solve with Julia Profile instead. Model construction
and compilation are excluded from the warmed solve timings. The synthetic
four-column plan is not a complete OCO retrieval benchmark.
`AUDIT_FAST_REFERENCE=true` runs the reference once for parity, then warms and
measures only the optimized path; use the default for a fresh complete A/B
timing comparison. The final large-batch runs used this option and reuse the
first-stage reference timing medians; Float64 reference arithmetic is unchanged.

For the sampled GPU memory record, run from `test/` with an existing output
folder and choose an available GPU index:

```bash
AUDIT_OUTPUT=/tmp/jacobian-results python3 \
  ../docs/dev_notes/jacobian_batched/monitor_gpu.py /path/to/julia 0 scalar-10000
```

An optional final argument `'Stokes_IQU()'` selects polarization. This wrapper
sets 10,000 wavelengths, four Julia threads and one BLAS thread. It samples
process memory approximately every 250 ms, including the reserved CUDA pool.
The full phases include compilation/warmup and the three timed solves. The
selected phases also include a full solve to retain outputs for parity, so
those memory peaks must not be interpreted as selected-only memory needs.

## Large-batch profiling and next targets

At 10,000 wavelengths, the first-stage scalar TimerOutput covers only 22% of
10.7 s total wall time. Optical properties are 57% of those timed sections;
that does **not** establish their share of total runtime. A subsequent warmed
host profile finds surface creation in 3,842 of 6,181 samples containing the
RT driver (about 62% inclusive), with 3,572 samples at the atmospheric-source
matrix-vector loop itself. This identified the source batching change above.
See `evidence/timer-new-cuda-n10000.log` and
`evidence/host-profile-cuda-n10000.log`; both precede the surface fix.

A separate warmed post-change CUDA timer sample assigns about 49% of measured
section time to optical-property/Jacobian assembly, 38% to doubling, 8% to
interaction and 4% to elemental work. These are host timer sections, not a
partition of device execution or total wall time. The full-solve bottleneck has
shifted toward the upstream derivative handoff; retaining and reusing its phase
and mixing work is now a higher priority. See `evidence/timer-new-cuda.log`.

The remaining launches include source attenuation, broadcasts and writeback,
forward products/inverses, and the special interaction cases. Further fusion
should target these measured costs. For operators beyond the measured CUDA
range, compare tiled product-rule kernels with vendor batched BLAS before
extending the current gate.
Reuse of upstream Mie/phase tangents and more general active layouts remain
separate opportunities. Each extension needs parity and finite-difference checks
as well as warmed runtime, allocation and launch-count measurements.

The opt-in [equivalent-source adding path](source_adding.md) integrates local
forcing vectors with Lambertian surface derivatives in the normal RT driver.
