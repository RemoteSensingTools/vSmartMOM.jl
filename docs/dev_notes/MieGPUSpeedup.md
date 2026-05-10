# GPU-Accelerated Mie Scattering: Implementation & Testing Plan

## 1. Motivation

The Mie scattering module (`src/Scattering/`) computes aerosol optical properties
(Greek coefficients, SSA, extinction) by integrating over 1000--2500 size-distribution
quadrature points.  Each point requires a Mie series expansion up to nmax~300 terms.
On CPU this takes 1--4 seconds per aerosol per wavelength; for a retrieval with 20
aerosol types across 1000 spectral points, Mie computation dominates wall-clock time.

Each radius point is independent, making this a natural GPU target.  The challenge:
NVIDIA L40S and consumer GPUs have only 1/64 FP64 throughput, but the Dn logarithmic
derivative recursion requires more than Float32's 24 mantissa bits.

## 2. Current State (Prototype, March 2026)

Three new files exist on the `unified-vsmartmom` branch:

| File | Purpose | Status |
|------|---------|--------|
| `src/Scattering/gpu_precision.jl` | DoubleSingle, ComplexDS, Neumaier primitives | Working, tested |
| `src/Scattering/gpu_mie_kernels.jl` | KernelAbstractions @kernel definitions (Kernels 1--5) | Working on CPU backend |
| `src/Scattering/compute_NAI2_gpu.jl` | GPU NAI2 entry point + CPU Neumaier helpers | Working, validated vs CPU Float64 |
| `test/test_mie_gpu.jl` | Accuracy + performance test suite | All 62 tests pass |

**Validated accuracy** (CPU backend, nquad=2500, vs existing CPU Float64 reference):

| Scheme | SSA rel err | k_ext rel err | Greek coefs max abs err |
|--------|------------|--------------|------------------------|
| GPU-NativeFloat64 | 4.5e-16 | 2.0e-16 | ~1e-14 |
| GPU-DSEmulated | 4.5e-16 | 2.0e-16 | ~1e-14 |

Both schemes match at machine epsilon because the Float64 reference and the GPU
path use the same mathematical formulation; the only difference is the kernel
decomposition and Neumaier compensation.  On actual L40S hardware, the DS path
will show ~1e-7 differences due to the ~48-bit DS mantissa vs 53-bit Float64.

**Current CPU-backend speedup**: ~4--5x (from batched computation structure).
Expected GPU speedup: 50--100x for nquad >= 1000.

## 3. Precision Architecture

### Three Tiers

| Tier | Technique | Ops/elem | Where |
|------|-----------|----------|-------|
| **1 -- Emulated FP64** | DoubleSingle{Float32} pairs (~48 bits) | ~4x overhead | Dn downward recursion only |
| **2 -- Compensated Float32** | Neumaier accumulator (Float32 + compensation) | ~7 FLOP/add | S1/S2 sums, cross-section dots, Greek quadrature, bulk_f reduction |
| **3 -- Native Float32** | Plain Float32 | 1 FLOP | Riccati-Bessel recursion, pi/tau Legendre, an/bn from Dn, phase matrix elements |

### DoubleSingle Arithmetic

Represents a value as `hi + lo` where `|lo| <= ulp(hi)/2`:
- `TwoSum(a, b)`: error-free addition (6 FLOPs, Knuth)
- `TwoProd(a, b)`: error-free multiplication (2 FLOPs, uses hardware `fma`)
- `ComplexDS{T}`: 4 Float32 values per complex number (re_hi, re_lo, im_hi, im_lo)

Measured accuracy:
- DS mul/div: ~1e-14 relative error (near Float64)
- DS add: ~1e-11 worst-case (due to lo-term renormalization)
- Dn recursion (x=500, m=1.5-0.1i): 5.8e-14 max relative error vs Float64

### Neumaier Compensated Summation

Standard Kahan fails when |addend| > |running sum|.  Neumaier's variant handles
both cases.  Cost: 7 FLOPs per addition vs 1 for naive sum.

Measured improvement over naive Float32:
- Harmonic series (10K terms): 34x less error
- Alternating series (10K terms): 618x less error

### Precision Policy Dispatch

```julia
abstract type MiePrecisionPolicy end
struct NativeFloat64 <: MiePrecisionPolicy end   # A100, V100
struct DSEmulated    <: MiePrecisionPolicy end   # L40S, consumer GPUs
```

Only the Dn recursion kernel differs between policies.  Everything else uses
the model's FT parameter with Neumaier compensation where needed.

## 4. Kernel Decomposition

### Kernel 1: Mie Coefficients (an, bn)

- **Grid**: (nquad_radius,) -- one thread per radius point
- **Algorithm**: DS Dn downward recursion -> Float32 Riccati-Bessel forward -> an/bn
- **Precision**: Tier 1 (DS or native F64) for Dn, Tier 3 for psi/chi and an/bn
- **Memory**: Dn state is register-only; final an/bn stored in (nquad_radius, nmax_global) arrays
- **Warp divergence**: nmax varies from ~10 to ~300+ across radii.  See Section 5.

### Kernel 2+3 (fused): Amplitude Functions + Phase Matrix

- **Grid**: (n_mu, nquad_radius) -- one thread per (angle, radius) pair
- **Algorithm**: Neumaier-compensated S1/S2 sum over l, then f11/f33/f12/f34
- **Precision**: Tier 2 (Neumaier) for l-summation, Tier 3 for phase matrix products
- **Shared memory opportunity**: Load an[1:nmax_i], bn[1:nmax_i] per block column

### Kernel 4: Cross-Sections + Size-Distribution Reduction

- **Kernel 4a**: (nquad_radius,) -- per-radius C_sca, C_ext with Neumaier
- **Kernel 4b**: (n_mu,) -- bulk phase matrix reduction over size distribution
- **Kernel 4c**: scalar -- bulk_C_sca, bulk_C_ext weighted sums
- **Precision**: Tier 2 (Neumaier) -- the size-distribution weights span 10 orders of magnitude

### Kernel 5: Greek Coefficients

- **Grid**: (l_max,) -- one thread per Legendre order l
- **Algorithm**: Neumaier dot products w_mu' * (bulk_f .* P)
- **Shared memory**: Load w_mu and bulk_f vectors once per block

## 5. Warp Divergence Mitigation

nmax varies from ~10 (small particles) to ~300+ (large particles).  Threads in the
same warp would do vastly different amounts of work.

**Solution -- Sorted bucketing**:
1. Sort radius points by nmax into 4--8 buckets
2. Within each bucket, pad to the bucket's max nmax
3. Launch one kernel per bucket
4. Eliminates warp divergence within each bucket

This is straightforward to implement but not yet done in the prototype.  The
prototype uses a single kernel launch with per-thread nmax bounds, which works
correctly but wastes GPU cycles on threads that finish early.

## 6. Implementation Phases

### Phase 1: DS + Neumaier Primitives [DONE]

File: `src/Scattering/gpu_precision.jl`

- DoubleSingle{T} with TwoSum, TwoProd, add, mul, div, inv
- ComplexDS{T} with add, sub, mul, div, inv
- NeumaierAccum{T} with add, sum
- ComplexNeumaier{T} for S1/S2 complex sums
- MiePrecisionPolicy types (NativeFloat64, DSEmulated)

### Phase 2: KernelAbstractions Kernels [DONE -- CPU backend]

File: `src/Scattering/gpu_mie_kernels.jl`

- mie_coefficients_kernel_ds! (DSEmulated)
- mie_coefficients_kernel_f64! (NativeFloat64)
- amplitude_phase_kernel! (fused Kernel 2+3)
- cross_sections_kernel!
- size_reduction_kernel!
- greek_coefficients_kernel! (not yet used in pipeline)

### Phase 3: NAI2 GPU Entry Point [DONE -- CPU backend]

File: `src/Scattering/compute_NAI2_gpu.jl`

- compute_aerosol_optical_properties_gpu()
- neumaier_dot, neumaier_dot_3, neumaier_matvec helpers
- Greek coefficient computation (currently CPU-side; Kernel 5 ready but not wired)

### Phase 4: CUDA Backend Integration [TODO]

**What to do:**
1. Add CUDA kernel launch to `ext/vSmartMOMCUDAExt.jl`:
   ```julia
   function Scattering.compute_aerosol_optical_properties_gpu(
       model::MieModel{<:NAI2}, ::GPU; kwargs...)
       backend = CUDA.CUDABackend()
       # ... call existing GPU entry point with CUDABackend
   end
   ```
2. Test with `CUDA.CuArray` allocations instead of `KernelAbstractions.allocate`
3. Benchmark on A100 (NativeFloat64) vs L40S (DSEmulated)
4. Implement sorted bucketing for warp divergence mitigation
5. Move Greek coefficient computation to GPU (Kernel 5) and size-distribution
   reduction (Kernel 4b) to GPU -- currently done on CPU after copy-back

**Key decision: when to copy back to CPU.**
Current approach copies f11/f33/f12/f34 and C_sca/C_ext back to CPU for reduction.
If Greek coefficients also run on GPU, only the final 6 Greek coefficient vectors
(~2400 floats each) need to return.  This reduces PCIe transfer by ~100x.

### Phase 5: PCW GPU Path [TODO]

File to create: `src/Scattering/compute_PCW_gpu.jl`

**Kernel P1: Averaged Mie Products** (compute_avg_anbns!)
- Grid: (Nmax, Nmax) upper-triangular
- Per thread: Neumaier reduction over nquad_radius for one (n,m) pair
- Alternative: reshape as batched outer-product, use CUBLAS sgemm

**Kernel P2: Sl Terms** (compute_Sl)
- Grid: (l_max,) with one thread block per l
- The alternating-sign double sum is the most cancellation-prone computation.
  Use pairwise tree reduction in shared memory (O(log N) error growth)
  plus Neumaier for the inner accumulators.
- Branchless sign: `1 - 2 * ((l + n + m) & 1)`

### Phase 6: Linearized GPU Path [TODO -- see Section 8]

### Phase 7: Auto-dispatch [TODO]

Wire into the model construction pipeline:
```julia
function compute_aerosol_optical_properties(model::MieModel{<:NAI2}, arch::GPU, ...)
    policy = has_native_fp64(arch) ? NativeFloat64() : DSEmulated()
    compute_aerosol_optical_properties_gpu(model, devi(arch);
        precision_policy=policy)
end
```

Where `has_native_fp64(arch)` queries CUDA compute capability >= 6.0 with
FP64 throughput >= 1/4 of FP32 (true for V100/A100, false for L40S/RTX).

## 7. Testing Strategy

### Unit Tests (in test/test_mie_gpu.jl)

| Test | What it validates | Tolerance |
|------|------------------|-----------|
| TwoSum exactness | Error-free addition (a+b = hi+lo exactly) | Exact (bit-for-bit) |
| TwoProd exactness | Error-free multiplication | Exact |
| DS add/mul/div | Relative error vs Float64 | add < 1e-9, mul/div < 1e-12 |
| ComplexDS mul/div | Relative error vs Complex{Float64} | < 1e-12 |
| Neumaier harmonic | Compensated vs naive Float32 | < 1e-3 absolute |
| Neumaier alternating | Compensated vs naive Float32 | Better than naive |
| Neumaier wide range | 10-order dynamic range sum | < 1e-5 relative |
| Dn recursion (6 cases) | DS Dn vs Float64 Dn for x=10..500 | < 1e-7 relative |

### Integration Tests

| Test | What it validates | Tolerance |
|------|------------------|-----------|
| NAI2 GPU-F64 vs CPU | NativeFloat64 full pipeline | SSA/k: rtol 1e-6, Greek: atol 1e-6 |
| NAI2 GPU-DS vs CPU | DSEmulated full pipeline | SSA/k: rtol 1e-4, Greek: atol 1e-3 |
| Linearized GPU vs CPU | Jacobians for all 4 parameters | rtol 1e-4 per element |

### Cross-GPU Validation

Run identical test cases on:
- A100 (NativeFloat64) -- ground truth
- L40S (DSEmulated) -- verify DS precision
- CPU (existing code) -- regression baseline

Compare results pairwise:
- A100 vs CPU: should agree to Float64 epsilon
- L40S vs A100: should agree to ~1e-7 (DS precision)
- L40S vs CPU: should agree to ~1e-7

### Performance Benchmarks

Time the full pipeline at nquad_radius = {100, 500, 1000, 2500} for:
- CPU baseline (existing serial code)
- GPU NativeFloat64 (A100)
- GPU DSEmulated (L40S)

Report: wall-clock time, TFLOPS utilization, PCIe transfer overhead.

### End-to-End RT Validation

Run forward RT with GPU Mie -> verify reflectances match CPU to < 0.1%.
Use the existing Jacobian tests from `test/test_helpers.jl` with
`rel_errors()` and `fd_jacobian_R()` to validate through the full RT chain.

## 8. Linearization Impact

### Current Analytic Linearization

The existing CPU linearization (`compute_NAI2_lin.jl`, `mie_helper_functions_lin.jl`)
computes derivatives of all aerosol optical properties with respect to 4 parameters:

| Index | Parameter | Symbol | Domain |
|-------|-----------|--------|--------|
| 1 | Real refractive index | nr | Mie coefficients |
| 2 | Imaginary refractive index | ni | Mie coefficients |
| 3 | Size distribution mean | r_m | Size distribution weights |
| 4 | Size distribution width | sigma | Size distribution weights |

The derivative chain is:

```
nr, ni  -->  Dn_dot  -->  an_dot, bn_dot  -->  S1_dot, S2_dot  -->  f_dot  -->  Greek_dot
                                                                        |
r_m, sigma  -->  wx_dot  ------------------------------------------->  |
```

### GPU Linearization Strategy

The linearized Dn recursion uses the **same mathematical structure** as the forward
recursion, with additional derivative terms propagated alongside.  From
`mie_helper_functions_lin.jl`:

```julia
# Forward: Dn_prev = ratio - 1/(Dn_prev + ratio)
# Linearized (for each derivative index ctr):
#   Dn_dot[ctr] = ratio_dot[ctr] + denom_inv^2 * (Dn_dot_prev[ctr] + ratio_dot[ctr])
```

This means each thread needs **2x the register state** for the Dn recursion:
- Forward: 4 Float32 (ComplexDS re_hi, re_lo, im_hi, im_lo) for Dn_prev
- + 2 derivatives: 4 Float32 each = 8 additional Float32 for Dn_dot[1], Dn_dot[2]
- Total: 12 Float32 registers per thread for the DS variant

For the NativeFloat64 variant: 2 Float64 (Dn_prev) + 4 Float64 (Dn_dot[1], Dn_dot[2])
= 6 Float64 = 12 Float32-equivalent registers.

**Register pressure is manageable** -- A100 has 65536 registers per SM, and at
256 threads/block this is 256 registers/thread.  12 extra Float32 is fine.

### Implementation Plan for Linearized GPU Kernels

#### Kernel 1-Lin: Mie Coefficients + Derivatives

Extend `mie_coefficients_kernel_ds!` to also output:
- `an_dot[i, n, 1:2]` -- derivatives of an w.r.t. nr, ni
- `bn_dot[i, n, 1:2]` -- derivatives of bn w.r.t. nr, ni

The Dn derivative recursion must also use DS arithmetic on L40S, since it
has the same numerical structure (continued fraction) and cancellation issues.

Layout: `an_dot` is `(nquad_radius, nmax_global, 2)` -- the derivative index
is the last dimension for memory coalescing (threads access consecutive radii).

#### Kernel 2+3-Lin: Amplitude + Phase Matrix Derivatives

Extend `amplitude_phase_kernel!` to also compute:
```julia
# Forward:  S1 += prefac * (an * tau + bn * pi)
# Derivative: S1_dot[ctr] += prefac * (an_dot[ctr] * tau + bn_dot[ctr] * pi)
```

Phase matrix derivatives follow the same bilinear structure:
```julia
# Forward:  f11 = inv_x2 * (|s1|^2 + |s2|^2)
# Derivative: f11_dot[ctr] = inv_x2 * Re(s1_dot*conj(s1) + s1*conj(s1_dot)
#                                       + s2_dot*conj(s2) + s2*conj(s2_dot))
```

Output: `f11_dot[imu, ir, 1:2]`, etc.

**Register pressure doubles** for the S1/S2 Neumaier accumulators (2 complex
accumulators become 6), but this is still within budget.

#### Kernel 4-Lin: Cross-Section Derivatives

```julia
# Forward:  C_sca = (2pi/k^2) * sum((2n+1) * (|an|^2 + |bn|^2))
# Derivative: C_sca_dot[ctr] = (2pi/k^2) * sum((2n+1) * Re(
#   an_dot[ctr]*conj(an) + an*conj(an_dot[ctr]) + same for bn))
```

#### Size Distribution Derivatives (Kernel 4b-Lin)

Derivatives w.r.t. r_m and sigma affect only the weights, not the per-radius
Mie computations.  This is a separate reduction:
```julia
# bulk_f11_dot[ctr=3,4] = sum(f11[:,i] * wr_dot[ctr-2, i])
```

Where `wr_dot` are the derivatives of `4*pi*r^2*wx` w.r.t. distribution parameters.
These can be computed on CPU (cheap) and uploaded as a second weight vector.

#### Kernel 5-Lin: Greek Coefficient Derivatives

Same structure as forward kernel but with `bulk_f_dot` inputs:
```julia
# beta_dot[ctr, l] = (2l+1)/2 * w_mu' * (bulk_f11_dot[ctr,:] .* P[:,l])
```

### Linearization Summary

| Component | Forward cost | Linearized cost | Overhead |
|-----------|-------------|-----------------|----------|
| Dn recursion (DS) | 4 ComplexDS vals | + 4 ComplexDS vals (2 derivs) | 2x registers |
| Riccati-Bessel | 6 Float32 | unchanged (no nr/ni dependence) | 0 |
| an/bn computation | ~20 FLOPs/n | + ~40 FLOPs/n (quotient rule, 2 derivs) | ~2x FLOPs |
| S1/S2 accumulation | 2 ComplexNeumaier | + 4 ComplexNeumaier (2 derivs) | 3x accumulators |
| Phase matrix | 4 FLOPs/point | + 16 FLOPs/point (product rule, 2 derivs) | ~4x FLOPs |
| Size reduction | 1 Neumaier/row | + 4 Neumaier/row (4 derivs) | 5x |
| Greek coefficients | 6 Neumaier dots | + 24 Neumaier dots (4 derivs) | 5x |

**Expected overall overhead**: ~2.5--3x compared to forward-only GPU.
This is similar to the CPU linearization overhead ratio.

### ForwardDiff on GPU

The existing `phase_function_autodiff.jl` uses ForwardDiff with `Dual{Float64}` numbers.
On GPU, ForwardDiff with `Dual{Float32, Float32, 4}` could work in principle but has issues:

1. **Register explosion**: Each Dual carries 4 partials, so every Float32 becomes
   5 Float32 values.  For ComplexDS, this would be 20 Float32 per complex number.
2. **DS + Dual composition**: `DoubleSingle{Dual{Float32}}` is theoretically possible
   but has not been tested and would be very expensive (40 Float32 per complex DS Dual).
3. **NNlib batched_mul already supports Dual**: The CoreRT module handles this,
   but the Scattering module's inner loops have not been tested with Dual on GPU.

**Recommendation**: Use explicit analytical linearization for the GPU path (Phase 6).
ForwardDiff can still be used on CPU as a reference for validating the GPU Jacobians.

## 9. Memory Layout & Workspace

### MieWorkspace Struct

Add to `src/Scattering/types.jl`:

```julia
struct MieWorkspace{FT, AT <: AbstractArray}
    # Kernel 1 outputs
    an   :: AT  # (nquad_radius, nmax_global) Complex{FT}
    bn   :: AT  # (nquad_radius, nmax_global) Complex{FT}

    # Kernel 2+3 outputs
    f11  :: AT  # (n_mu, nquad_radius) FT
    f33  :: AT  # (n_mu, nquad_radius) FT
    f12  :: AT  # (n_mu, nquad_radius) FT
    f34  :: AT  # (n_mu, nquad_radius) FT

    # Kernel 4 outputs
    C_sca :: AT  # (nquad_radius,) FT
    C_ext :: AT  # (nquad_radius,) FT

    # Read-only inputs
    leg_pi  :: AT  # (n_mu, nmax_global) FT
    leg_tau :: AT  # (n_mu, nmax_global) FT
    x_params :: AT  # (nquad_radius,) FT
    nmax_per_r :: AbstractVector{Int}

    # Linearized outputs (optional, nothing if forward-only)
    an_dot :: Union{Nothing, AT}  # (nquad_radius, nmax_global, 2) Complex{FT}
    bn_dot :: Union{Nothing, AT}  # (nquad_radius, nmax_global, 2) Complex{FT}
    f11_dot :: Union{Nothing, AT}  # (n_mu, nquad_radius, 4) FT
    # ... etc
end
```

**Memory estimate** (nquad=2500, n_max=300, n_mu=599):

| Array | Shape | Bytes (FP32) |
|-------|-------|-------------|
| an, bn | (2500, 300) x2 | 12 MB |
| f11..f34 | (599, 2500) x4 | 24 MB |
| C_sca, C_ext | (2500,) x2 | 20 KB |
| leg_pi, leg_tau | (599, 300) x2 | 1.4 MB |
| **Forward total** | | **~38 MB** |
| an_dot, bn_dot | (2500, 300, 2) x2 | 24 MB |
| f_dot (x4 x4) | (599, 2500, 4) x4 | 96 MB |
| **Linearized total** | | **~158 MB** |

Both fit comfortably in GPU memory (A100: 40 GB, L40S: 48 GB).

### Workspace Reuse

For multi-wavelength computations, the workspace can be allocated once and
reused across wavelengths (nmax changes, but the maximum across wavelengths
can be pre-computed).  This avoids repeated GPU allocation overhead.

## 10. Integration with RTModel Pipeline

### Current Flow

```
YAML -> parameters_from_yaml() -> vSmartMOM_Parameters
     -> model_from_parameters() -> RTModel
        [internally: for each aerosol/wavelength:
           compute_aerosol_optical_properties(mie_model)]  <- THIS IS THE TARGET
     -> rt_run(model) -> (R, T)
```

### Modified Flow

```
model_from_parameters(params, arch=GPU()) ->
    for each aerosol/wavelength:
        if arch isa GPU
            compute_aerosol_optical_properties_gpu(mie_model, devi(arch);
                precision_policy=auto_detect_policy(arch))
        else
            compute_aerosol_optical_properties(mie_model)  # existing CPU path
        end
```

The GPU path returns the same `AerosolOptics` struct, so downstream RT code
is unchanged.  The linearized path returns `(AerosolOptics, linAerosolOptics)`
regardless of CPU/GPU backend.

## 11. Risk Assessment

| Risk | Severity | Mitigation |
|------|----------|------------|
| DS precision insufficient for extreme size parameters (x>1000) | Low | Validated up to x=500; DS provides 48 bits vs needed ~40 |
| Warp divergence kills GPU speedup | Medium | Sorted bucketing eliminates it; implementation straightforward |
| Register pressure too high for linearized DS kernel | Medium | A100 has 256 regs/thread at 256 threads/block; only need ~60 |
| Neumaier not sufficient for PCW alternating sums | Medium | Pairwise tree reduction as backup; PCW is Phase 5, not blocking |
| KernelAbstractions overhead on CUDA | Low | KA is thin wrapper; benchmark against raw CUDA.jl kernels |
| ForwardDiff Dual on GPU for validation | Low | Use CPU ForwardDiff as reference; don't need AD on GPU |

## 12. File Inventory

### New Files (to create)

| File | Phase | Description |
|------|-------|-------------|
| `src/Scattering/gpu_precision.jl` | 1 | DS, Neumaier, precision policies [DONE] |
| `src/Scattering/gpu_mie_kernels.jl` | 2 | @kernel Mie computations [DONE] |
| `src/Scattering/compute_NAI2_gpu.jl` | 3 | GPU NAI2 entry point [DONE] |
| `src/Scattering/compute_PCW_gpu.jl` | 5 | GPU PCW entry point |
| `src/Scattering/gpu_mie_kernels_lin.jl` | 6 | Linearized GPU kernels |
| `src/Scattering/compute_NAI2_gpu_lin.jl` | 6 | Linearized GPU NAI2 entry |
| `test/test_mie_gpu.jl` | Testing | Accuracy + performance tests [DONE] |
| `test/test_mie_gpu_lin.jl` | Testing | Linearized GPU tests |

### Existing Files to Modify (wiring only, no logic changes)

| File | Change |
|------|--------|
| `src/Scattering/Scattering.jl` | `include()` lines + exports [DONE] |
| `src/Scattering/types.jl` | Add `MieWorkspace` struct |
| `ext/vSmartMOMCUDAExt.jl` | Register GPU dispatch for Mie |
| `test/runtests.jl` | Add GPU Mie test set |
| `test/Project.toml` | `KernelAbstractions` dependency [DONE] |

### Files NOT Modified

All existing CPU code paths remain untouched:
- `src/Scattering/mie_helper_functions.jl`
- `src/Scattering/mie_helper_functions_lin.jl`
- `src/Scattering/compute_NAI2.jl`
- `src/Scattering/compute_NAI2_lin.jl`
- `src/Scattering/compute_PCW.jl`
- `src/Scattering/legendre_functions.jl`
- `src/Scattering/phase_function_autodiff.jl`

## 13. Recommended Next Steps

1. **Phase 4**: Wire CUDA backend, benchmark on A100 (1--2 days)
2. **Bucketing**: Implement sorted radius bucketing for warp divergence (1 day)
3. **Full GPU pipeline**: Move Greek coefficients + reduction to GPU (1 day)
4. **Phase 6**: Linearized GPU kernels for nr/ni derivatives (2--3 days)
5. **Phase 6b**: Size distribution derivatives (wx_dot) on GPU (1 day)
6. **Phase 5**: PCW GPU path with pairwise reduction (2--3 days)
7. **Phase 7**: Auto-dispatch based on GPU compute capability (0.5 day)
8. **Benchmarking**: Full comparison A100 vs L40S vs CPU across test suite (1 day)
9. **Phase W1--W5**: Wavelength-dependent optical properties (see Section 14)

## 14. Wavelength-Dependent Aerosol Optical Properties

### 14.1 Current Limitation

The existing pipeline computes Mie properties **once per spectral band** at the
band-center wavelength (`model_from_parameters.jl:165`):

```julia
mie_model = make_mie_model(...,
    (maximum(curr_band_λ) + minimum(curr_band_λ)) / 2,  # band center only
    ...)
```

This means:
- `AerosolOptics.k` and `AerosolOptics.ω̃` are **scalars** per band
- Greek coefficients are **constant** across all spectral points within a band
- Z matrices are computed once and **replicated** to `(NquadN, NquadN, nSpec)` via
  `expandOpticalProperties()` in `compEffectiveLayerProperties.jl`
- The refractive index (`nᵣ`, `nᵢ`) is **constant globally** (from a single `n_ref`
  value), with no wavelength dependence

For narrow bands (~1 nm) this is adequate.  For broader spectral ranges or
high-accuracy retrievals, this introduces systematic bias because:
- The refractive index of real aerosols varies with wavelength (e.g., ammonium
  sulfate: n varies from 1.53 at 0.3 um to 1.50 at 1.0 um)
- The size parameter `x = 2pi*r/lambda` changes across the band, shifting the
  Mie resonance structure
- SSA and extinction can vary by several percent across a wide band

### 14.2 Goal

Compute **per-spectral-point** aerosol optical properties:
- Greek coefficients `alpha(lambda), beta(lambda), ...`
- Single-scattering albedo `omega(lambda)`
- Extinction cross-section `k(lambda)`
- Z matrices `Z++(lambda), Z-+(lambda)`

with accuracy equivalent to running independent Mie computations at each
wavelength, but with overhead **< 2x** compared to the current single-band
computation when running on GPU.

### 14.3 Strategy: Three Options

#### Option A: Brute-Force GPU Batch (recommended first)

Compute full Mie at **every spectral point**, parallelized over the product
space `(nquad_radius x N_lambda)`.

```
GPU grid for Kernel 1:  (nquad_radius * N_lambda,)
                        = 2500 * 1000 = 2.5M threads

GPU grid for Kernel 2+3: (n_mu, nquad_radius * N_lambda)
```

**Advantages:**
- Exact: no interpolation error
- Simplest to implement: just add a wavelength loop dimension to existing kernels
- GPU fills easily: 2.5M threads for Kernel 1 is excellent occupancy

**Disadvantages:**
- Memory: storing an/bn for all (radius, wavelength) pairs simultaneously
  requires ~12 GB.  Must stream over wavelength chunks (see 14.8).
- Total FLOP count is N_lambda x single-wavelength, but GPU absorbs this
  via parallelism -- wall-clock time scales sub-linearly.

#### Option B: Anchor + Greek Coefficient Interpolation

1. Choose `N_anchor` wavelengths (5--10) spanning the spectral range
2. Run full Mie pipeline at each anchor wavelength
3. Interpolate the 6 Greek coefficient vectors to all target wavelengths
4. Interpolate SSA and k_ext (smooth in wavelength)
5. Compute Z matrices from interpolated Greek coefficients (cheap)

**Advantages:**
- Minimal Mie computation (N_anchor << N_lambda)
- Greek coefs are smooth functions of lambda when refractive index is smooth
- Z-matrix computation from Greek coefs is O(l_max^2 * Nquad^2) per wavelength,
  roughly 100x cheaper than full Mie

**Disadvantages:**
- Interpolation error depends on N_anchor and spectral range
- Need to validate: for what band widths and N_anchor is the error acceptable?
- Greek coefficients at high Legendre order can have fine structure

#### Option C: Anchor an/bn + Downstream Recomputation

1. Compute an/bn (Kernel 1 only) at N_anchor wavelengths for all radii
2. Interpolate an/bn to all target wavelengths (complex-valued, smooth)
3. Run Kernels 2--5 at each target wavelength from the interpolated an/bn

**Advantages:**
- More accurate than Option B (interpolates at the Mie coefficient level
  rather than the derived Greek coefficient level)
- an(x, m) varies more smoothly than the resulting Greek coefficients
  because no quadrature integration has been applied yet

**Disadvantages:**
- Still need to run Kernels 2--5 at every wavelength (though these are
  cheaper than Kernel 1)
- Interpolation of complex-valued an/bn requires care near Mie resonances

**Recommendation**: Start with Option A (brute-force GPU batch) as the baseline.
It produces exact results and validates the entire per-wavelength pipeline.
Then benchmark to see if Options B or C are needed for performance.

### 14.4 Refractive Index Model

The current `Aerosol` struct stores constant `nr`, `ni`.  For wavelength-dependent
properties we need a refractive index function `m(lambda)`.

**Proposed interface:**

```julia
abstract type AbstractRefractiveIndex end

# Current behavior: constant refractive index
struct ConstantRI{FT} <: AbstractRefractiveIndex
    nr::FT
    ni::FT
end

# Table-based: linear interpolation of tabulated data
struct TabulatedRI{FT} <: AbstractRefractiveIndex
    lambda::Vector{FT}   # wavelengths [um], sorted
    nr::Vector{FT}       # real part at each wavelength
    ni::Vector{FT}       # imaginary part at each wavelength
end

# Cauchy dispersion model: n = A + B/lambda^2 + C/lambda^4
struct CauchyRI{FT} <: AbstractRefractiveIndex
    A::FT
    B::FT
    C::FT  # often zero
    ni_model::AbstractRefractiveIndex  # imaginary part can be separate
end

# Evaluate at a wavelength
refractive_index(ri::ConstantRI, lambda) = (ri.nr, ri.ni)
refractive_index(ri::TabulatedRI, lambda) = (interp(ri.lambda, ri.nr, lambda),
                                              interp(ri.lambda, ri.ni, lambda))
```

**YAML config extension:**

```yaml
scattering:
  rt_aerosols:
    - type: "sulfate"
      refractive_index:
        type: "tabulated"
        file: "data/sulfate_ri.csv"   # columns: wavelength_um, nr, ni
      # OR:
      refractive_index:
        type: "cauchy"
        A: 1.52
        B: 0.0042
```

For backward compatibility, the current `nr`, `ni` scalar fields map to
`ConstantRI` automatically.

### 14.5 Type System Changes

The `AerosolOptics` struct already permits arrays in `omega`, `k`, `ft`:

```julia
struct AerosolOptics{FT}
    greek_coefs::GreekCoefs
    omega::Union{FT, AbstractArray{FT}}   # already supports vector
    k::Union{FT, AbstractArray{FT}}       # already supports vector
    ft::Union{FT, AbstractArray{FT}}      # already supports vector
end
```

**GreekCoefs extension** (two options):

Option 1 -- Make fields 2D:
```julia
mutable struct GreekCoefs{FT}
    alpha::AbstractArray{FT}  # (l_max,) or (l_max, N_lambda)
    beta::AbstractArray{FT}
    ...
end
```
Pro: minimal struct changes.  Con: all code accessing `beta[l]` must handle
both 1D and 2D cases.

Option 2 -- Vector of GreekCoefs:
```julia
struct SpectralAerosolOptics{FT}
    optics::Vector{AerosolOptics{FT}}  # one per spectral point
end
```
Pro: clean separation, no existing code changes.  Con: more allocations,
harder to batch on GPU.

**Recommendation**: Option 1 (2D fields) is more GPU-friendly because a single
contiguous `(l_max, N_lambda)` array can be transferred to the device in one copy.
Add a helper `is_spectral(gc::GreekCoefs) = ndims(gc.beta) == 2` for dispatch.

**AerosolState changes:**

```julia
struct AerosolState{FT, AO, AT}
    aerosol_optics::Vector{Vector{AO}}  # [iBand][iAer] -- AO can now be spectral
    tau_aer::Vector{AT}                 # [iBand] -> (nSpec, nAer, nLayers) or (nAer, nLayers)
end
```

The tau_aer computation in `model_from_parameters.jl:188-191` currently computes:
```julia
tau_aer[i_band][i_aer,:] = tau_ref * (k_band/k_ref) * profile_shape
```
With per-wavelength k(lambda), this becomes:
```julia
tau_aer[i_band][i_spec, i_aer, :] = tau_ref * (k[i_spec]/k_ref) * profile_shape
```

### 14.6 Z-Matrix Computation

Currently `compute_Z_moments()` returns `(Z++, Z-+)` with shape `(NquadN, NquadN)`,
which `expandOpticalProperties()` replicates to `(NquadN, NquadN, nSpec)`.

With per-wavelength Greek coefficients:
- Compute Z matrices per wavelength from per-wavelength Greek coefs
- `Z_aer` becomes `(NquadN, NquadN, nSpec)` natively -- no replication needed
- Cost: O(l_max^2 * Nquad^2) per wavelength per Fourier moment m

**Cost estimate** for l_max=200, Nquad=16 (scalar), N_lambda=1000:
- Per wavelength: ~200 * 16^2 = 51K FLOPs (trivial)
- Total: 51M FLOPs (< 1ms on CPU, negligible on GPU)

For Stokes_IQUV with Nquad=16: 200 * 64^2 = 820K FLOPs per wavelength,
total 820M FLOPs -- still < 100ms on CPU, < 1ms on GPU.

**GPU parallelization**: Z-matrix computation for all wavelengths can be batched
as a single kernel with grid `(N_lambda, l_max)` or similar.  The B matrix
construction and Pi matrix products are independent across wavelengths.

**Key insight**: The Legendre polynomials P, R, T depend only on the quadrature
angles mu, not on wavelength.  They can be computed **once** and reused for all
wavelengths.  Only the B matrices (which depend on Greek coefficients) change
per wavelength.

### 14.7 Linearization Implications

#### 14.7.1 Forward-Mode Chain Rule

With wavelength-dependent refractive index, the derivative chain becomes:

```
RI model params (A, B, ...)
  --> nr(lambda_j), ni(lambda_j)       via refractive_index()
    --> Dn(lambda_j)                   via downward recursion
      --> an(lambda_j), bn(lambda_j)   via Bohren-Huffman
        --> Greek(lambda_j)            via NAI2/PCW integration
          --> Z(lambda_j)              via compute_Z_moments
            --> R(lambda_j), T(lambda_j)   via RT solver
```

The existing 4-parameter Jacobian `[nr, ni, r_m, sigma]` changes:

| Parameter | Current | Per-wavelength |
|-----------|---------|----------------|
| nr | Scalar, same for all lambda | nr(lambda_j) from RI model, derivative via chain rule through model params |
| ni | Scalar, same for all lambda | ni(lambda_j), same as above |
| r_m | Scalar (distribution mean) | Unchanged -- distribution params are wavelength-independent |
| sigma | Scalar (distribution width) | Unchanged |

#### 14.7.2 If Using Tabulated RI

For tabulated nr(lambda), ni(lambda), the "parameters" are the table values
themselves.  If we want to retrieve nr at specific wavelengths, the Jacobian
is `d(Greek)/d(nr_j)` which is non-zero only at wavelength j (block-diagonal
in wavelength).  This is the simplest case.

#### 14.7.3 If Using a Dispersion Model

For Cauchy dispersion `nr(lambda) = A + B/lambda^2`, the derivatives are:
```
d(nr)/dA = 1                     (constant)
d(nr)/dB = 1/lambda^2            (wavelength-dependent)
```

The Jacobian `d(Greek(lambda_j))/dA = d(Greek)/d(nr) * d(nr)/dA` chains
through the existing `d(Greek)/d(nr)` derivative (already computed in
`compute_NAI2_lin.jl`) multiplied by the model derivative.

This is a **single matrix-vector product** per spectral point, adding
negligible overhead.

#### 14.7.4 GPU Linearization for Per-Wavelength

The linearized Dn recursion (`mie_helper_functions_lin.jl:54-70`) computes
`Dn_dot[1,n]` (d/d(nr)) and `Dn_dot[2,n]` (d/d(ni)) alongside the forward
values.  For per-wavelength computation:

- Each (radius, wavelength) thread computes its own Dn_dot values
- The derivative chain is **identical** to the single-wavelength case, just
  with different nr(lambda_j), ni(lambda_j) values
- Register overhead: same 2x as described in Section 8
- No new mathematical structure needed

The key difference is that the final Jacobian `d(R)/d(RI_params)` now sums
contributions from all wavelengths:
```
dR/dA = sum_j  dR/dGreek(lambda_j) * dGreek(lambda_j)/dnr(lambda_j) * dnr(lambda_j)/dA
```

This summation can be done after the per-wavelength GPU computation, on CPU,
since it involves only the compact RI model parameters.

### 14.8 Memory Management: Streaming Over Wavelengths

Full per-wavelength storage for nquad=2500, n_max=300, n_mu=599, N_lambda=1000:

| Array | Shape | FP32 Bytes |
|-------|-------|-----------|
| an, bn | (2500, 300, 1000) x 2 Complex | 12.0 GB |
| f11..f34 | (599, 2500, 1000) x 4 | 24.0 GB |
| C_sca, C_ext | (2500, 1000) x 2 | 20 MB |
| Greek coefs (6) | (599, 1000) x 6 | 14 MB |

**Problem**: an/bn alone require 12 GB, exceeding practical limits for a
single kernel launch alongside other working memory.

**Solution: Wavelength-chunked streaming**

Process wavelengths in chunks of `N_chunk` (e.g., 50--100):

```
for chunk in partition(1:N_lambda, N_chunk):
    Kernel 1: compute an/bn for all radii, chunk wavelengths
    Kernel 2+3: compute S1/S2, phase matrix for chunk
    Kernel 4: reduce over radii -> per-wavelength cross-sections, bulk phase matrix
    Kernel 5: Greek coefficients from bulk phase matrix
    Copy Greek coefs for chunk back to CPU
    # an/bn and f11..f34 can be discarded (or overwritten by next chunk)
end
```

Peak GPU memory per chunk (N_chunk=100):
| Array | Bytes |
|-------|-------|
| an, bn (2500, 300, 100) x 2 | 1.2 GB |
| f11..f34 (599, 2500, 100) x 4 | 2.4 GB |
| Other working memory | ~0.5 GB |
| **Total** | **~4.1 GB** |

This fits easily in both A100 (40 GB) and L40S (48 GB), leaving room for the
RT solver's own GPU allocations.

**Alternative: Fused Kernel Approach**

Instead of materializing the full f11..f34 arrays, fuse Kernels 1--4 so that
each thread block processes one (wavelength, angle_block) and reduces over
radii immediately in shared memory:

```
Per block: wavelength lambda_j, angles mu[start:end]
  For each radius i:
    Compute an/bn (Kernel 1, per-thread registers)
    Compute S1/S2 (Kernel 2, Neumaier in registers)
    Compute f11/f33/f12/f34 (Kernel 3, registers)
    Accumulate into shared-memory reducers (Kernel 4)
  Block-reduce -> bulk_f11[mu, lambda_j]
```

This eliminates the f11/f34 intermediate arrays entirely.  Peak memory drops
to just the output Greek coefficient arrays (~14 MB).  However, this fusion
requires rewriting the kernel decomposition and is a Phase W5 optimization.

### 14.9 Changes to the RT Pipeline

#### model_from_parameters.jl

Replace the band-center Mie call (line 165) with:

```julia
# All wavelengths in this band
curr_band_lambda = FT.(1e4 ./ params.spec_bands[i_band])

# Get per-wavelength refractive index
nr_vec = [refractive_index(aerosol.ri_model, lam)[1] for lam in curr_band_lambda]
ni_vec = [refractive_index(aerosol.ri_model, lam)[2] for lam in curr_band_lambda]

# GPU batch Mie at all wavelengths
aerosol_optics_raw = compute_aerosol_optical_properties_gpu_spectral(
    mie_model, curr_band_lambda, nr_vec, ni_vec, backend;
    precision_policy=policy)

# Truncate (now per-wavelength)
aerosol_optics[i_band][i_aer] = truncate_phase_spectral(truncation_type, aerosol_optics_raw)
```

#### compEffectiveLayerProperties.jl

Skip `expandOpticalProperties()` for aerosols with per-wavelength Z matrices:

```julia
# Current: Z++ is (NquadN, NquadN, 1), expanded to (NquadN, NquadN, nSpec)
# New: Z++ is already (NquadN, NquadN, nSpec) for spectral aerosols
if is_spectral(aerosol_optics)
    # Z matrices already have spectral dimension -- no expansion needed
    AerZ_pp, AerZ_mp = compute_Z_moments_spectral(pol_type, mu, aerosol_optics.greek_coefs, m)
else
    # Legacy: scalar Greek coefs, expand as before
    AerZ_pp, AerZ_mp = compute_Z_moments(pol_type, mu, aerosol_optics.greek_coefs, m)
end
```

### 14.10 Implementation Phases

#### Phase W1: Refractive Index Infrastructure (1--2 days)

- Define `AbstractRefractiveIndex`, `ConstantRI`, `TabulatedRI` types
- Add `ri_model` field to `Aerosol` struct (default: `ConstantRI(nr, ni)`)
- Implement `refractive_index(model, lambda)` dispatch
- YAML parsing for `refractive_index:` config block
- Test: verify backward compatibility with existing constant-RI configs

#### Phase W2: Batch GPU Mie Over Wavelengths (2--3 days)

- New function: `compute_aerosol_optical_properties_gpu_spectral()`
- Extend Kernel 1 to accept per-wavelength refractive indices
- Wavelength-chunked streaming loop (see 14.8)
- Output: `AerosolOptics` with vector-valued `k`, `omega`, 2D Greek coefs
- Test: compare against running independent single-wavelength GPU calls

#### Phase W3: Per-Wavelength Z Matrices (1--2 days)

- New function: `compute_Z_moments_spectral(pol_type, mu, greek_coefs_2d, m)`
- Batch B-matrix construction across wavelengths (reuse P, R, T)
- Modify `compEffectiveLayerProperties.jl` to handle spectral aerosol Z matrices
- Remove `expandOpticalProperties()` for spectral aerosols
- Test: end-to-end RT with per-wavelength aerosol, compare against wavelength-by-
  wavelength independent runs

#### Phase W4: Linearization Update (2--3 days)

- Extend `linAerosolOptics` to per-wavelength derivatives
- Compute `d(Greek(lambda_j))/d(nr(lambda_j))` at each spectral point (existing chain)
- Add chain rule for RI model parameters: `d(nr)/d(A)`, `d(nr)/d(B)`
- Validate: GPU linearized per-wavelength vs CPU finite-difference

#### Phase W5: Optimization (optional, 2--3 days)

- Anchor + interpolate (Options B/C from 14.3) if brute-force is too slow
- Fused kernel approach to eliminate intermediate f11/f34 arrays
- Benchmark: determine N_anchor needed for < 0.1% interpolation error
- Only pursue if Phase W2 benchmarks show brute-force is insufficient

### 14.11 Testing Strategy for Wavelength-Dependent Properties

| Test | Method | Acceptance Criterion |
|------|--------|---------------------|
| RI model consistency | `ConstantRI(1.3, 0.001)` vs current scalar nr/ni | Bit-identical results |
| Per-wavelength vs independent | Run N independent single-lambda Mie, compare to batch | < 1e-10 relative error |
| Spectral Z matrices | Compare `Z_spectral[:,:,j]` vs `compute_Z_moments(Greek[j])` | Bit-identical |
| End-to-end RT accuracy | Spectral aerosol RT vs band-center-only RT | Document differences (not a pass/fail -- this IS the improvement) |
| Interpolation accuracy (Phase W5) | Anchor=5,10,20 vs brute-force N_lambda=1000 | < 0.1% in reflectance |
| Linearized consistency | Analytic Jacobian vs finite-difference | rel_errors < 1e-4 per element |
| Performance | Wall-clock: per-wavelength GPU vs per-band CPU | < 2x overhead target |

### 14.12 Impact on Existing Code

**No changes to existing CPU Mie code.**  The per-wavelength pipeline is a
new code path alongside the existing one.  `ConstantRI` refractive index
plus a single-element wavelength vector reproduces current behavior exactly.

| File | Change Type | Description |
|------|-------------|-------------|
| `src/Scattering/types.jl` | Add types | `AbstractRefractiveIndex`, RI implementations, `is_spectral()` |
| `src/Scattering/compute_NAI2_gpu.jl` | Extend | Spectral batch entry point |
| `src/CoreRT/tools/model_from_parameters.jl` | Modify | Use per-wavelength Mie when RI model is spectral |
| `src/CoreRT/LayerOpticalProperties/` | Modify | Handle spectral aerosol Z matrices |
| `src/IO/Parameters.jl` | Extend | Parse `refractive_index:` YAML block |
| `src/Scattering/compute_Z_matrices.jl` | Extend | Batch Z computation over wavelengths |
