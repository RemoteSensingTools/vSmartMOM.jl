# Jacobian bottleneck investigation — 2026-09-06

Candidate: `integration/surface-split-multisensor`, `b524a0d85ae36ab0eecf8b23a939f4d6ce62c456`. Companion to the [release readiness audit](release_readiness_2026-09-06.md). No solver changes were used for these measurements.

**The first optimization target is the all-parameter doubling/interaction path: allocating matrix-product chains and copied derivative slices, repeated once per active parameter.** Selective Jacobian layouts already provide a substantial practical benefit. The unconditional documentation promise that a full forward-plus-Jacobian run costs less than twice a forward solve does not hold for these supported workloads.

## Method and measured scope

The retained [benchmark/profiling script](release_audit_2026_09_06/jacobian_profile.jl) derives a small scene from `test/test_parameters/JacobianTestFast.yaml`: one aerosol, five layers, Stokes I, three streams, embedded solar, a 6×6 diffuse operator, and a Lambertian surface. Absorption is removed to isolate RT propagation from line-by-line gas construction. The aerosol radius ceiling is 3 μm, lognormal median 0.15 μm and width 1.4. The spectral grid has either 64 or 512 points.

Full solves have fourteen native columns: surface pressure, seven aerosol parameters, five layer-resolved gas slots and surface albedo. The five gas slots are allocated even in this no-absorption case; they carry zero absorption tangents. A synthetic four-column `JacobianPlan` selects pressure, aerosol reference optical depth, aerosol profile location and surface albedo. This uses the existing `OCO_RRS_synth` vocabulary with a manually constructed one-band layout; it is **not** a benchmark of a complete three-band OCO retrieval.

The forward, full and selected solves use the same already-built model. Each exact path is warmed before three timed samples; reported times are medians, with compilation and model construction excluded. GC is requested before each sample and its subsequent cost remains included in elapsed time. CUDA is synchronized at measurement boundaries. Output is suppressed consistently; internal timing instrumentation remains enabled. Separate measurements cover model construction.

Environment: Julia 1.12.6; CUDA.jl 5.11.3; NVIDIA A100 PCIe 40 GB, driver 560.35.03; AMD EPYC 7H12; Float64; BLAS fixed to one thread. CPU runs use one or four Julia threads; CUDA runs use four host threads. This is a shared server with other active work, including the release test suite. Treat absolute times and especially thread-scaling differences as diagnostic measurements, not dedicated-machine performance guarantees. Sample values are retained in the logs.

## Warmed solve results: 64 spectral points

| Backend | Forward | Full, 14 columns | Selected, 4 columns | Full / forward | Full / selected |
|---|---:|---:|---:|---:|---:|
| CPU, 1 Julia thread | 0.0888 s | 1.568 s | 0.450 s | 17.65× | 3.48× |
| CPU, 4 Julia threads | 0.1273 s | 3.895 s | 1.268 s | 30.60× | 3.07× |
| CUDA A100 | 0.0844 s | 3.078 s | 0.852 s | 36.47× | 3.61× |

Forward and full-Jacobian radiances differ by at most `5.3e-18` on CPU and `3.3e-17` on CUDA. Selected radiances and selected derivative columns agree exactly with the corresponding full outputs on each backend. Thus the selection speedup is not obtained by changing the requested numerical result.

The single-thread CPU full solve allocates **1.550 GB** cumulatively versus **0.440 GB** for the selected solve and **0.093 GB** for forward. These are total host bytes allocated during a call, not peak resident memory. Median GC time is 102 ms for full and 15 ms for selected. Allocation reduction is therefore useful beyond GC time: constructing, copying and filling temporaries also consumes time.

CUDA full and selected solves allocate 263 MB and 69 MB on the host. Julia's `@timed.bytes` does **not** measure device memory traffic or peak VRAM; these numbers must not be presented as GPU memory usage.

## Larger spectral batch and a misleading forward baseline

| Backend, 512 spectral points | Forward | Full, 14 columns | Selected, 4 columns | Full / forward | Full / selected |
|---|---:|---:|---:|---:|---:|
| CPU, 1 Julia thread | 3.528 s | 11.417 s | 2.654 s | 3.24× | 4.30× |
| CUDA A100 | 3.120 s | 3.377 s | 0.955 s | 1.08× | 3.54× |

Selection remains exact. The CUDA full solve grows only about 10% when the wavelength count increases eightfold, consistent with substantial fixed submission overhead in the smaller case. CPU full allocations grow to 11.93 GB cumulatively; selected allocations are 3.43 GB. This batch demonstrates useful GPU throughput, while retaining a large benefit from selecting the requested columns.

However, the apparently favorable 1.08× CUDA full/forward ratio is partly a **forward-path cache defect**. A separate CPU forward profile attributes 74% of measured section time and 86% of allocations to `CacheInit` (3.14 s and 2.35 GiB in that sample). Sampling resolves this to `ZMomentTables` construction.

At `compEffectiveLayerProperties.jl:203`, the cache uses `length(aerosol_optics[iB][i].greek_coefs.β)` as the angular expansion length. Multiwavelength Greek arrays are matrices, so this multiplies the angular row count by the number of wavelengths. `ZMomentTables` then allocates six `(nμ, l_max, l_max)` arrays (`compute_Z_matrices.jl:53`), giving an unintended quadratic dependence on spectral count. Use the angular axis, `size(β,1)`, with equivalent handling for other Greek sources, and test matrix-valued coefficients. This is release finding R8, not merely a benchmark-tuning suggestion.

The [cache probe](release_audit_2026_09_06/cache_dimensions.jl) confirms `size(β) == (5,512)` versus `length(β) == 2560`: 1,887,436,800 bytes for the six inflated base tables instead of 7,200 bytes. Disabling the cache through its existing diagnostic flag reduces the CPU forward solve to **0.460 s**, allocating 416.6 MB, with **bit-identical R and T**. The flag is restored after the probe. This demonstrates the source of the baseline inflation without changing solver implementation.

The forward/Jacobian ratios above intentionally retain the unmodified candidate behavior. Correcting the cache will reduce the forward denominator and change those ratios. Do not use the defective baseline to substantiate the documentation's unconditional performance claim.

## Where the time and allocations go

The package's timer sections for a representative single-thread CPU full solve attribute 67.4% of measured section time to doubling, 22.3% to interaction, 9.0% to optical-property construction, and 1.2% to elemental kernels. Doubling and interaction account for approximately 64% and 28% of measured allocations. The two explicitly timed interaction inversions together take about 6.9 ms, or 0.4% of measured section time. These percentages describe the named sections, not an independently instrumented instruction-level partition.

The CUDA host timers similarly place roughly 70% in doubling and 21% in interaction. Nested host timers may include asynchronous submission or synchronization at different points; they are not a device-kernel timeline. The synchronized whole-call measurements above remain the performance comparison.

Julia CPU sampling finds BLAS matrix multiplication, array allocation/copying and GC. Four-thread sampling additionally contains substantial scheduler/condition-variable/lock waiting. Allocation sampling (`Profile.Allocs`, rate 0.005) traces the largest groups to the following production sites:

| Site | Mechanism and implication |
|---|---|
| `src/CoreRT/tools/cpu_batched.jl:77` | Every allocating dense `batched_mul` creates a new 3D result. This is the largest sampled allocation group. |
| `src/CoreRT/CoreKernel/doubling_lin.jl:274` | Each parameter evaluates a nested chain for the inverse tangent. `@views` avoids slice copies here, but does not make the `⊠` products allocation-free. |
| `src/CoreRT/CoreKernel/doubling_lin.jl:297` | Source-tangent updates build additional nested products per parameter and doubling step. |
| `src/CoreRT/CoreKernel/doubling_lin.jl:315` | R/T updates lack `@views`, so right-hand derivative slabs are copied in addition to allocating product results. |
| `src/CoreRT/CoreKernel/interaction_lin.jl:232` | Each layer interaction allocates several full 4D derivative work arrays, with further arrays for the reflected direction at line 276. |
| `src/CoreRT/CoreKernel/interaction_lin.jl:246` | Per-parameter interaction chains copy slices and allocate intermediate products; analogous source and mirrored-direction expressions repeat the pattern. |

The actual fused elemental path calls `doubling_allparams!` from `rt_kernel_lin.jl:112`; it does not propagate only three core-variable derivatives throughout the atmosphere. `doubling_allparams_helper!` loops over the active physical parameter count. Existing reusable doubling buffers and hoisted forward products are useful, but large allocating expression chains remain inside that loop.

On CPU, `cpu_batched.jl:78` launches a `Threads.@threads` loop for every dense batched product. At 6×6 matrices and 64 wavelengths, repeated thread scheduling can cost more than the arithmetic it distributes. Four threads were slower than one in this run; the shared-machine condition means this is evidence to investigate a serial threshold/coarser partition, not a universal recommendation to disable threading.

On CUDA, analogous per-parameter products submit many small batched GEMMs, broadcasts and slice copies. Moreover, `ext/gpu_batched_cuda.jl:222` explicitly implements `_as_cuarray3(view) = copy(view)` for the derivative-view overloads: adding `@views` alone does not remove device copies on this path. A workspace/in-place change must also support contiguous device views without materializing them.

A separate warmed **CUPTI device profile** confirms submission and transfer overhead for the 64-point full solve: 72,002 `cuLaunchKernel` calls, 35,215 calls of the dominant small GEMM kernel, 114,192 host-to-device copy calls, and 4,676 device-to-device copy calls. The instrumented trace lasts 6.03 s; CUDA APIs account for 2.65 s (44%) and recorded device activity for 1.15 s (19%). Profiling overhead increases the duration relative to the uninstrumented 3.08 s solve, so 19% is the trace's activity fraction, not an asserted production utilization statistic. Copy counts do not establish transfer byte volume or prove which data structures caused every transfer. Nevertheless, these counts give a concrete target: fewer small launches, copied views and repeated transfer setup. The raw summary is retained in `evidence/jacobian-cuda-device-profile.log`.

## Upstream model-construction cost

For the 64-point single-thread CPU scene, a forward build takes 4.76 ms, a full linearized build 93.83 ms, and a build with aerosol microphysics and H2O Jacobians disabled 7.88 ms. Allocations are respectively 6.33 MB, 88.84 MB and 7.82 MB. CUDA full construction is also about 120 ms here.

Upstream derivative construction is expensive relative to the forward builder, but still much smaller than the full RT solve in this scene. Start with RT propagation and reuse model/optics work when scenes permit. Disabling microphysics is only appropriate when those derivatives are not requested. The release audit's R1 normalization defect also means this flag currently needs a correctness fix before being treated as a universally physics-preserving performance switch.

`lin_model_from_parameters.jl:151` reserves `1 + N_var_gas` species slots, even when absorption is absent. Preserve legacy column semantics for callers, but allow explicit active layouts to avoid propagating unused/zero slots. The demonstrated 3.5× selection benefit includes avoiding both genuinely unrequested aerosol derivatives and these gas slots.

## Optimization order and acceptance criteria

1. **Use active layouts for retrievals with a small state vector.** Keep the selected/full radiance and column parity checks. Extend supported public selection vocabulary where real retrievals need it; do not silently remove legacy output columns.
2. **Replace allocating product chains with workspace-backed operations.** First repair the forward cache dimension error so comparisons use a valid baseline. For Jacobians, start at doubling lines 274–320 and interaction lines 232–304. Introduce views for derivative slabs and use the existing in-place batched machinery with reusable 3D scratch buffers. Respect aliasing: several formulas need old R/T/source values until all terms are evaluated. A blanket search-and-replace with in-place calls could corrupt tangents.
3. **Reuse interaction workspaces across layers and Fourier moments.** Size buffers for the active layout and keep lifetime/backend ownership explicit. Measure allocated bytes as well as runtime to show that the repeated 4D allocations were removed.
4. **Reduce per-parameter launch and scheduling overhead.** Evaluate coarser CPU work partitioning and a small-work serial path; for CUDA, consider batching over wavelength × parameter or fused derivative kernels. Preserve backend portability and benchmark real operator sizes before adding complexity.
5. **Optimize upstream construction after the propagation changes are measured.** Reuse invariant Mie/phase work and skip unrequested upstream derivatives through the existing plan boundary. Include gas-rich builds before extrapolating the no-absorption measurements to line-by-line retrieval workloads.

Every optimization should retain forward/full/selected parity and finite-difference checks for pressure, aerosol, gas and surface parameters, including the normalization and Cox–Munk fixes identified in the release report. Measure Float32/Float64, more streams and polarized operators, multiple layer counts, representative spectral grids and actual retrieval state sizes. The current audit diagnoses a real bottleneck and invalidates an unconditional performance promise; it does not establish an optimized implementation or a complete production performance envelope.

## Reproduction

From the candidate's `test/` directory, after instantiating its test environment:

```bash
AUDIT_BACKEND=cpu AUDIT_NSPEC=64 JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  julia --project=. ../docs/dev_notes/release_audit_2026_09_06/jacobian_profile.jl

AUDIT_BACKEND=cuda AUDIT_NSPEC=64 JULIA_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1 \
  julia --project=. ../docs/dev_notes/release_audit_2026_09_06/jacobian_profile.jl
```

Set `AUDIT_FORWARD_DIAGNOSTIC=true` to collect forward timer sections and a CPU sampling profile without rerunning the full benchmark matrix. Set `AUDIT_OUTPUT` to an existing writable directory for profiles. The script's default evidence directory is specific to this audit. Set `AUDIT_NSPEC=512` for the larger batch. Set `AUDIT_DEVICE_PROFILE=true` with the CUDA backend for a separate warmed CUPTI profile; this omits the benchmark loop. CPU allocation samples report sampled bytes only and must not be confused with the whole-call allocation totals.
