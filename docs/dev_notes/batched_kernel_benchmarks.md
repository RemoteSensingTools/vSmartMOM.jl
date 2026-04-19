# Custom KA Kernels vs cuBLAS: Batched Operation Benchmarks

**Date:** 2026-04-03  
**GPUs tested:** NVIDIA A100-PCIE-40GB, NVIDIA L40S  
**Benchmark scripts:** `test/benchmarks/batched_matmul_benchmark.jl`, `test/benchmarks/batched_inv_benchmark.jl`, `test/benchmarks/batched_fused_benchmark.jl`, `test/benchmarks/batched_fused_v2_benchmark.jl`, `test/benchmarks/interaction11_bundled_benchmark.jl`

## Motivation

vSmartMOM's RT solver performs thousands of batched matrix multiplications and inversions on 3D arrays of shape `(NquadN, NquadN, nSpec)`. These currently dispatch to cuBLAS (`gemm_strided_batched` for matmul, `getrf_strided_batched` + `getri_strided_batched` for inversion). For the small matrix sizes typical in RT (N=4-32), these operations are **launch-bound, not flop-bound** — cuBLAS kernel launch overhead dominates. Custom KernelAbstractions.jl kernels fuse operations and exploit shared memory to eliminate this overhead.

Typical matrix sizes in vSmartMOM:
- Stokes_I, Nquad=4-16: NquadN = 4-16
- Stokes_IQU, Nquad=4-16: NquadN = 12-48
- Stokes_IQUV, Nquad=4-16: NquadN = 16-64

---

## Part 1: Standalone Kernel Benchmarks

### Batched Matrix Multiplication

Two custom variants tested:

1. **KA global-memory**: One thread per output element `C[i,j,k]`, accumulates dot product.
2. **KA shared-memory (smem)**: One workgroup of N×N threads per batch element. Loads A and B tiles into shared memory. Requires N² ≤ 1024 (N ≤ 32).

### Batched Matrix Inversion

Three custom variants tested:

1. **Gauss-Jordan (GJ)**: N threads per workgroup, augmented matrix `[A|I]` in shared memory, partial pivoting. Shared memory: `N × 2N × sizeof(FT)`.
2. **LU sequential**: LU factorization + sequential forward/back substitution (thread 1 only). Stable but slow.
3. **LU parallel (LU-par)**: LU factorization, then each of N threads solves one column independently. Shared memory: `2 × N² × sizeof(FT) + N × sizeof(Int32)`.

### Results: A100

**Matmul — smem kernel speedup over cuBLAS:**

| N | Float32 | Float64 |
|---|---------|---------|
| 4 | 1.6-1.7x | 1.5-2.2x |
| 8 | 1.6-2.5x | 1.7-3.3x |
| 12 | 1.7-3.2x | 2.1-5.3x |
| 16 | 1.6-2.7x | 1.9-2.5x |
| 24 | 1.2-1.5x | 0.9-1.2x |
| 32 | 0.7-1.2x | 0.6-0.9x |

**Inversion — LU-par speedup over cuBLAS:**

| N | Float32 | Float64 |
|---|---------|---------|
| 4 | 2.1-2.7x | 1.9-2.4x |
| 8 | 1.7-2.7x | 1.5-2.6x |
| 12 | 1.6-2.8x | 1.6-2.6x |
| 16 | 0.8-2.5x | 0.6-2.0x |
| 24 | 1.3-1.5x | 1.3-4.2x |
| 32 | 0.3-0.6x | 0.7-1.4x |
| 48 | 0.6-0.9x | 1.4-1.7x |

### Results: L40S (Float32 only — no real FP64)

**Matmul smem:** Beats cuBLAS at every size up to N=32 (1.4–4.5x).  
**Inversion LU-par:** Beats cuBLAS at N ≤ 24 (1.5–3.0x).

---

## Part 2: Fused RT Microbenchmarks

### Key Insight

In the RT solver, inversions **never appear in isolation**. The dominant motifs are:

1. **Geometric progression** (doubling, `rt_helpers.jl:34-39`):
   `tt_gp = t⁺⁺ * inv(I - r⁻⁺ * r⁻⁺)`
2. **Interaction primitive**:
   `X = B * inv(I - A₁ * A₂)`

The scripts `test/benchmarks/batched_fused_benchmark.jl` and
`test/benchmarks/batched_fused_v2_benchmark.jl` benchmark these motifs in
isolation. They remain useful because they show that:

- launch overhead dominates for small `N`
- keeping intermediates in shared memory is a real win
- replacing "form inverse, then multiply" with "factor once, solve directly"
  lowers scratch usage and usually improves performance

However, these microbenchmarks are **not enough to choose production dispatch
for the full RT interaction path**. They do not include the full
`interaction_11` update, source terms, or all six output fields.

Practical takeaway from the microbenchmarks:

- fused custom kernels are clearly worthwhile to pursue for small matrices
- the main open question is **where they still win once the full
  `interaction_11` algebra is included**

---

## Part 3: Full `interaction_11` Benchmark

### What Was Added

New script:

- `test/benchmarks/interaction11_bundled_benchmark.jl`

This benchmark evaluates the actual forward `ScatteringInterface_11` update
from `src/CoreRT/CoreKernel/interaction.jl`, not just a simplified primitive.

It compares three implementations:

1. **Current**: the existing two-inverse expression tree
2. **Reduced algebra**: a one-inverse cuBLAS formulation
3. **Bundled**: a custom KA kernel using the reduced algebra, with one LU
   factorization and three solves per batch element

### Reduced Algebra

Using

- `(I - AB)⁻¹ A = A (I - BA)⁻¹`
- `(I - AB)⁻¹ = I + A (I - BA)⁻¹ B`

with `A = r⁻⁺`, `B = R⁺⁻`, the full `interaction_11` update can be rewritten to
use **one factorization of**

- `M = I - R⁺⁻ r⁻⁺`

and three solves:

- `X_T = M \ T⁺⁺`
- `X_R = M \ (R⁺⁻ t⁻⁻)`
- `x_J = M \ (J₀⁺ + R⁺⁻ j₀⁻)`

Then:

- `T⁺⁺_out = t⁺⁺ X_T`
- `R⁺⁻_out = r⁺⁻ + t⁺⁺ X_R`
- `J₀⁺_out = j₀⁺ + t⁺⁺ x_J`

and with `P = T⁻⁻ r⁻⁺`:

- `T⁻⁻_out = T⁻⁻ t⁻⁻ + P X_R`
- `R⁻⁺_out = R⁻⁺ + P X_T`
- `J₀⁻_out = J₀⁻ + T⁻⁻ j₀⁻ + P x_J`

This is the first optimization that should go into production even without a
custom kernel: it removes one full inverse/factorization from the current path.

### Benchmark Conditions

- GPU: NVIDIA A100-PCIE-40GB
- Batch sizes tested: `500`, `2000`, `5000`, `10000`, `50000`
- Float32 runs were done with:
  - `CUDA.math_mode!(CUDA.PEDANTIC_MATH)`
  - this disables possible TF32 tensor-core math in cuBLAS, so the Float32
    comparison is true FP32 vs true FP32

### Representative Results (A100, PEDANTIC_MATH for Float32)

#### Float32

| N | Batch | Current ms | Reduced ms | Bundled ms | Reduced speedup | Bundled speedup |
|---|-------|-----------:|-----------:|-----------:|----------------:|----------------:|
| 12 | 5000 | 1.7091 | 1.5944 | 0.4024 | 1.07x | **4.25x** |
| 24 | 5000 | 2.6911 | 2.1156 | 1.5831 | 1.27x | **1.70x** |
| 32 | 5000 | 2.9481 | 2.2876 | 4.0806 | 1.29x | 0.72x |
| 48 | 5000 | 6.2198 | 4.7606 | 10.5953 | 1.31x | 0.59x |
| 12 | 50000 | 8.4429 | 7.9708 | 1.8627 | 1.06x | **4.53x** |
| 24 | 50000 | 13.3489 | 10.6414 | 8.4275 | 1.25x | **1.58x** |
| 32 | 50000 | 20.0028 | 15.8034 | 36.7360 | 1.27x | 0.54x |
| 48 | 50000 | 58.5472 | 45.2393 | 100.5916 | 1.29x | 0.58x |

#### Float64

| N | Batch | Current ms | Reduced ms | Bundled ms | Reduced speedup | Bundled speedup |
|---|-------|-----------:|-----------:|-----------:|----------------:|----------------:|
| 12 | 5000 | 2.4986 | 2.3101 | 0.6103 | 1.08x | **4.09x** |
| 24 | 5000 | 4.8169 | 3.0392 | 1.8125 | 1.58x | **2.66x** |
| 32 | 5000 | 6.5137 | 4.1462 | 7.0349 | 1.57x | 0.93x |
| 48 | 5000 | 24.4122 | 14.2172 | 21.2337 | 1.72x | **1.15x** |
| 12 | 50000 | 11.3316 | 10.6271 | 2.7894 | 1.07x | **4.06x** |
| 24 | 50000 | 43.9747 | 27.1800 | 16.6502 | 1.62x | **2.64x** |
| 32 | 50000 | 59.7463 | 37.9413 | 64.2683 | 1.57x | 0.93x |
| 48 | 50000 | 242.4832 | 142.2510 | 192.8980 | 1.70x | **1.26x** |

### What These Results Mean

1. **The reduced algebra is worth integrating regardless of custom kernels.**
   It consistently beats the current two-inverse path, usually by:
   - `~1.1x-1.3x` in Float32
   - `~1.1x-1.7x` in Float64

2. **The bundled kernel is clearly valuable at small `N`.**
   On A100:
   - Float32: strong win through `N=24`
   - Float64: strong win through `N=24`

3. **The bundled kernel is not "always better".**
   Full `interaction_11` is more demanding than the earlier microbenchmarks.
   The bundled kernel loses for:
   - Float32: `N >= 32` in the tested range
   - Float64: `N=32` is roughly crossover / slightly worse

4. **There is still a high-batch Float64 opportunity at `N=48`.**
   The bundled kernel loses at `N=48, batch=500`, but wins at larger batches on
   A100 Float64. This is real, but should be treated as a backend-specific
   tuning opportunity, not the generic default.

### Conservative A100 Dispatch from the Full `interaction_11` Benchmark

| Dtype | Use bundled kernel | Use reduced cuBLAS |
|------|---------------------|--------------------|
| Float32 | `N <= 24` | `N >= 32` |
| Float64 | `N <= 24` | `N >= 32` |

Notes:

- `Float64, N=48` can win for large enough batches on A100, but that should be
  guarded by backend-specific tuning data before becoming default dispatch.
- These recommendations supersede the earlier over-optimistic claim that the
  fused interaction kernel wins "at nearly all RT-relevant sizes".

### Observed Numerical Agreement

For the random matrices tested in `interaction11_bundled_benchmark.jl`:

- Float32 `max abs diff`: typically `1e-7` to low `1e-6`
- Float64 `max abs diff`: typically `1e-15` to low `1e-14`

This is good, but it is **not yet sufficient for production acceptance**.

---

## Part 4: Caveats and Future Validation

### TF32 / TensorFloat on A100

The A100 may use TF32 tensor-core math for cuBLAS Float32 operations in default
math mode. That matters because:

- the custom KA kernel is true FP32 arithmetic
- the cuBLAS baseline may otherwise be using TF32 internally
- performance and numerical comparisons are not apples-to-apples in that case

For this reason, the full `interaction_11` Float32 benchmark was run with
`CUDA.PEDANTIC_MATH`.

Future benchmark tables should report both:

- **PEDANTIC_MATH** (true FP32 baseline)
- **DEFAULT_MATH** (production-like A100 baseline, TF32 may be used)

### Precision / Accuracy Tests Still Needed

Before production dispatch is added, the following validation is still needed:

1. **Conditioning sweep**
   - test matrices where `I - R⁺⁻r⁻⁺` ranges from well-conditioned to nearly
     singular
   - report both output error and solve residual `||MX-B||`

2. **Seeded reproducibility**
   - benchmark and accuracy runs should use fixed RNG seeds
   - the current benchmark is fine for trend-finding, but not yet a canonical
     regression dataset

3. **Kernel-level unit verification**
   - add unit tests that compare the bundled kernel and reduced-algebra path
     against the clean reference implementation in `interaction.jl`
   - verify all six `interaction_11` outputs across representative `N`,
     dtypes, and batch sizes
   - once a kernel fuses this much logic, it is easy for it to drift from the
     core code path where the math is still readable and directly auditable

4. **End-to-end RT accuracy**
   - compare radiances / reflectances after plugging the reduced path and the
     bundled path back into full forward RT
   - matrix-level agreement alone is not enough

5. **Full production RT benchmark run**
   - after integration, run the full forward RT production benchmarks, not just
     isolated matrix kernels
   - include realistic atmospheric cases, representative `nSpec`, and both CPU
     and GPU baselines
   - the production benchmark is needed to confirm that the matrix-kernel
     speedups survive the surrounding RT workflow and memory traffic

6. **Float32 TF32 sensitivity**
   - explicitly compare:
     - reduced-cuBLAS with `PEDANTIC_MATH`
     - reduced-cuBLAS with `DEFAULT_MATH`
     - bundled custom kernel

7. **Backend-specific revalidation**
   - CUDA A100 results should not be assumed to transfer to L40S, Metal, or ROCm

### Workspace Cleanup Before Integration

The current CUDA batched code (`ext/gpu_batched_cuda.jl`) still has allocation
waste in hot paths:

1. `batch_inv!(..., ws::RTWorkspace)` still allocates fresh `pivot`/`info`
2. `batch_solve!` allocates a temporary inverse buffer
3. `SubArray` batched multiply paths materialize copies

These should be fixed regardless of custom kernel adoption, since they slow
down both the baseline and the fallback paths.

### Backend Portability Notes

KernelAbstractions.jl provides **source portability**, not performance
portability.

- **CUDA**:
  - A100 needs explicit TF32 policy in benchmark and dispatch thinking
  - the full `interaction_11` benchmark suggests small-`N` custom kernels plus
    reduced-cuBLAS fallback

- **Metal**:
  - treat as Float32 / mixed-precision backend
  - vendor matmul should still be benchmarked first
  - custom reduced/bundled solve kernels are likely more valuable than on CUDA
    because public LU/inverse vendor paths are weaker

- **ROCm**:
  - expect similar small-`N` custom-kernel opportunities
  - do not reuse CUDA cutoffs without benchmarking

### Shared Memory Constraints

For the current bundled `interaction_11` kernel used in
`interaction11_bundled_benchmark.jl`, shared memory is:

- `2 × N² × sizeof(FT) + N × sizeof(FT) + N × sizeof(Int32)`

Representative usage:

| N | Float32 | Float64 |
|---|---------|---------|
| 24 | 4.7 KB | 9.3 KB |
| 32 | 8.2 KB | 16.4 KB |
| 48 | 18.4 KB | 36.6 KB |

So the current bundled kernel comfortably fits A100 shared-memory limits for
the tested sizes. The crossover is performance-driven, not capacity-driven.

---

## Recommendations (Updated)

1. **Integrate the reduced one-inverse algebra into `interaction_11` first.**
   This is the lowest-risk and most general improvement. It consistently beats
   the current two-inverse formulation.

2. **Add bundled-kernel dispatch only for the sizes where the full benchmark
   shows a clear win.**
   On A100, the conservative rule is:
   - Float32: bundled for `N <= 24`
   - Float64: bundled for `N <= 24`

3. **Treat larger-`N` custom dispatch as backend-specific tuning, not a global
   default.**
   The `Float64, N=48` A100 result is promising at large batch, but needs more
   validation before widening the cutoff.

4. **Do not make production decisions from the microbenchmarks alone.**
   `batched_fused_benchmark.jl` and `batched_fused_v2_benchmark.jl` are useful
   evidence, but the full `interaction_11` benchmark is the real dispatch guide.

5. **Run dedicated precision and accuracy tests before merging production
   dispatch changes.**
   Performance trends are now clear; correctness under realistic RT conditions
   still needs its own validation pass.

6. **Require both kernel unit tests and a full RT production benchmark before
   enabling default dispatch.**
   The bundled kernel is complex enough that it needs direct verification
   against the clean math implementation, and the final go/no-go decision
   should come from full RT timing, not matrix benchmarks alone.
