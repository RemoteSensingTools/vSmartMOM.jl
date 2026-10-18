# 7 · Architecture-Agnostic Code & GPU Acceleration

> **For:** anyone who cares why this package runs on CPU + CUDA + Metal from a single source tree, and why line-by-line RT on the GPU is the *design point* — not an afterthought. Power users, contributors adding kernels, performance-conscious retrieval developers.
>
> **Prev:** [6 · Linearization](06_linearization.md) · **Next:** [8 · Inelastic Extension](08_inelastic.md)

This page is the *differentiator*. The matrix-operator method is well-known
in the RT literature; many codes implement it. What sets vSmartMOM.jl apart
is *how* it implements it — one source tree, three GPU backends, ForwardDiff
flowing through the GPU path, polarization and architecture as types not
runtime branches, line-by-line on GPU at full bandwidth.

## All wavelengths in parallel: the design point

> Line-by-line RT on the GPU works because **the spectral axis IS the batch
> axis**. No correlated-k or k-distribution required. (You can still add
> correlated-k or SVD-based dimensionality reductions on top — they remain
> compatible — but they are not prerequisites.)

Every RT array in vSmartMOM has shape `(NquadN, NquadN, nSpec)` where
``\text{NquadN} = N_\mathrm{quad} \cdot n_\mathrm{stokes}``. The first two dimensions are the
discretized angular streams (Concepts/02); the third is the **spectral axis**
— the wavenumber grid. The elemental, doubling, and interaction kernels
operate on these arrays via `batched_mul` and `batch_inv!`, which broadcast
across that third dimension.

Concrete consequence: a typical OCO-2/3-style retrieval problem with
``N_\mathrm{quad} = 24``, ``n_\mathrm{stokes} = 3``, and ``n_\mathrm{Spec} = 5\,000`` line-by-line points
in the O₂ A-band runs as **one batched matmul** (shape `(72, 72, 5000)`) per
layer per Fourier moment. On a modern GPU, that single launch saturates
memory bandwidth. The CPU/MOM path that other RT codes take —
"loop over wavelengths, do a matmul for each" — would be 5000 sequential
launches.

The constant-`N_doubl` trick from [Concepts/03c](03c_mixing.md) is what
makes the batched call uniform across the spectral axis. Without it,
`N_doubl` would vary wavelength-to-wavelength (because `τ_λ` does), and the
batched call would have to be split into chunks.

```
   One batched call per layer per Fourier moment m:

   3D array shape       ──→     batched_mul ⊠     ──→    Broadcasts across
   (NquadN,                     (CUBLAS / KA / BLAS         nSpec ~10⁴–10⁶
    NquadN,                      depending on backend)      wavelengths
    nSpec)
```

## Architecture types and dispatch

The whole abstraction lives in [`src/Architectures.jl:1–189`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/src/Architectures.jl#L1-L189) (~190 lines). It
defines three architecture types and a four-method dispatch surface:

```julia
abstract type AbstractArchitecture end
struct CPU       <: AbstractArchitecture end
struct GPU       <: AbstractArchitecture end           # CUDA
struct MetalGPU  <: AbstractArchitecture end           # Apple Silicon

@inline devi(::CPU) = KernelAbstractions.CPU()
# devi(::GPU), devi(::MetalGPU) are injected by extensions when CUDA / Metal load.

@inline architecture(::Array) = CPU()
# architecture(::CuArray), architecture(::MtlArray) injected by extensions.

@inline array_type(::CPU) = Array
# array_type(::GPU), array_type(::MetalGPU) injected by extensions.

@inline synchronize_if_gpu() = sync_device()           # backend-aware sync
```

The four methods (`devi`, `architecture`, `array_type`, `synchronize_if_gpu`)
plus `default_architecture()` are all there is. CPU is wired up in this file;
GPU backends inject the rest at extension load.

## GPU backends via Julia 1.9 package extensions

vSmartMOM's CUDA backend code is not part of the main module. CUDA support
arrives through a package extension (Julia ≥ 1.9 feature) at
[`ext/vSmartMOMCUDAExt.jl:22–129`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/ext/vSmartMOMCUDAExt.jl#L22-L129). When the user does
`using vSmartMOM, CUDA` in the same session, Julia's extension loader runs
this code:

```julia
module vSmartMOMCUDAExt

using vSmartMOM
using CUDA

# These three lines are the entire backend dispatch:
Architectures.devi(::vSmartMOM.Architectures.GPU) = CUDA.CUDABackend(; always_inline=true)
Architectures.array_type(::vSmartMOM.Architectures.GPU) = CuArray
Architectures.architecture(::CuArray) = vSmartMOM.Architectures.GPU()

function __init__()
    if CUDA.functional()
        try
            CoreRT.CUBLAS_ref[] = CUDA.CUBLAS                # for batched matmul
            test_arr = CUDA.CuArray([1.0f0])
            CUDA.synchronize()
            Architectures._has_cuda[] = true
            Architectures._sync_gpu[] = CUDA.synchronize     # for synchronize_if_gpu
            CUDA.allowscalar(false)                          # enforce kernel discipline
        catch e
            @warn "vSmartMOM GPU initialization failed, falling back to CPU" exception=e
            Architectures._has_cuda[] = false
        end
    end
end

end
```

Three things to notice:

1. **CPU execution does not require a working CUDA device.** CUDA is
   currently a direct dependency in `Project.toml`, and is also required by
   CanopyOptics. Installing the package therefore resolves CUDA even on a
   CPU-only machine. Metal is an optional weak dependency.
2. **The dispatch is three method overrides.** The kernels themselves don't
   change.
3. **`CUDA.allowscalar(false)`** enforces that no kernel does `A[1,1,1] = …`
   on a `CuArray`. Scalar element access on GPU arrays is a performance
   cliff; the policy makes that cliff a hard error.

The Metal extension at [`ext/vSmartMOMMetalExt.jl:20–174`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/ext/vSmartMOMMetalExt.jl#L20-L174) follows the same
pattern, replacing `CuArray` with `MtlArray`.

## `KernelAbstractions.@kernel` is the abstraction layer

Where a CUDA-only code would have `@cuda_kernel function …`, vSmartMOM has
`@kernel function …` from KernelAbstractions.jl. The same source compiles
to CPU (multithreaded), CUDA, and Metal.

The D-matrix kernel from [`src/CoreRT/CoreKernel/doubling.jl:186–208`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/src/CoreRT/CoreKernel/doubling.jl#L186-L208) is a good
representative — small, complete, and identical across backends:

```julia
@kernel function apply_D!(n_stokes::Int, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻)
    iμ, jμ, n = @index(Global, NTuple)
    i = mod1(iμ, n_stokes)
    j = mod1(jμ, n_stokes)

    if (i > 2)
        r⁻⁺[iμ,jμ,n] = -r⁻⁺[iμ,jμ,n]
    end

    if ((i <= 2) & (j <= 2)) | ((i > 2) & (j > 2))
        r⁺⁻[iμ,jμ,n] =  r⁻⁺[iμ,jμ,n]
        t⁻⁻[iμ,jμ,n] =  t⁺⁺[iμ,jμ,n]
    else
        r⁺⁻[iμ,jμ,n] = -r⁻⁺[iμ,jμ,n]
        t⁻⁻[iμ,jμ,n] = -t⁺⁺[iμ,jμ,n]
    end
end
```

Same kernel runs on CPU + CUDA + Metal. The launcher
`apply_D_matrix!` (in [`src/CoreRT/CoreKernel/doubling.jl`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/src/CoreRT/CoreKernel/doubling.jl)) picks the backend
with `devi(architecture(r⁻⁺))` and instantiates the kernel for that backend. No backend-specific code path, no preprocessor, no
`#ifdef CUDA_AVAILABLE`.

The same pattern is used in **gas absorption** (line-shape kernels at
[`src/Absorption/compute_absorption_cross_section.jl:316–332`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/src/Absorption/compute_absorption_cross_section.jl#L316-L332)) and **Mie
scattering** (the optional GPU Mie path at
[`src/Scattering/gpu_mie_kernels.jl:30–158`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/src/Scattering/gpu_mie_kernels.jl#L30-L158)). Three different physical
modules, one kernel-abstraction discipline.

## Batched matrix algebra: `⊠` and `batch_inv!`

The two operators that show up everywhere in the kernel are `⊠ = CoreRT.batched_mul`
(package-owned, so vSmartMOM's backend methods do not change `NNlib.batched_mul`
for other packages) and `batch_inv!`. Each has three backend implementations sharing a single
user-visible interface.

**CPU** ([`src/CoreRT/tools/cpu_batched.jl:49–100`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/src/CoreRT/tools/cpu_batched.jl#L49-L100)) — threaded over the
spectral slice with one BLAS `gemm!` per slice:

```julia
function batched_mul(A::Array{T,3}, B::Array{T,3}) where {T<:LinearAlgebra.BLAS.BlasFloat}
    C = Array{T,3}(undef, size(A,1), size(B,2), size(A,3))
    Threads.@threads for k in 1:size(C,3)
        @views mul!(C[:,:,k], A[:,:,k], B[:,:,k])
    end
    return C
end

function batch_inv!(X::AbstractArray{FT,3}, A::AbstractArray{FT,3}) where {FT}
    Threads.@threads for i = 1:size(A, 3)
        @views X[:,:,i] = A[:,:,i]\I
    end
end
```

**CUDA** ([`ext/gpu_batched_cuda.jl:225–239`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/ext/gpu_batched_cuda.jl#L225-L239) for `batched_mul`, [`99–111`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/ext/gpu_batched_cuda.jl#L99-L111) for `batch_inv!`) — CUBLAS `gemm_strided_batched`
in a single launch covering all spectral points:

```julia
function vSmartMOM.CoreRT.batched_mul(A::CuArray{FT,3}, B::CuArray{FT,3}) where {FT}
    CUDA.CUBLAS.gemm_strided_batched('N', 'N', A, B)        # ONE launch, all wavelengths
end

function vSmartMOM.CoreRT.batch_inv!(X::CuArray{FT,3}, A::CuArray{FT,3}) where {FT}
    pivot, info = CUDA.CUBLAS.getrf_strided_batched!(A, true)
    CUDA.CUBLAS.getri_strided_batched!(A, X, pivot)
end
```

**Metal** ([`ext/vSmartMOMMetalExt.jl:102–141`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/ext/vSmartMOMMetalExt.jl#L102-L141) plus
`src/CoreRT/tools/ka_batched_kernels.jl:56–63, 128–206`) — portable
KernelAbstractions matmul + LU-with-pivoting in `@localmem` shared memory.
No vendor BLAS dependency:

```julia
@kernel function _batched_mul_kernel!(C, A, B, ::Val{K}) where {K}
    i, j, k = @index(Global, NTuple)
    s = zero(eltype(C))
    @inbounds for l in 1:K
        s += A[i, l, k] * B[l, j, k]
    end
    @inbounds C[i, j, k] = s
end

@kernel function _batched_inv_lu_par_kernel!(X, A, ::Val{N}) where {N}
    # full LU decomposition with partial pivoting in @localmem; ~80 lines.
    # See ka_batched_kernels.jl:128–206 for the body.
end
```

The user-visible operator `⊠` is the *same* across all three. Dispatch is at
the array type. Adding a fourth backend (e.g. ROCm/AMD) would mean writing
~50 lines of method overrides, not rewriting the kernels.

## ForwardDiff Duals through the GPU path

The CUDA batched algebra also supports ForwardDiff `Dual` arrays for callers
that need them. This is separate from the production analytic RT tangent
path, which propagates plain floating-point values and derivative arrays.
The overload is defined in
[`ext/gpu_batched_cuda.jl:255–291`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/ext/gpu_batched_cuda.jl#L255-L291):

```julia
function vSmartMOM.CoreRT.batched_mul(
        A::CuArray{ForwardDiff.Dual{T,V,N},3},
        B::CuArray{ForwardDiff.Dual{T,V,N},3}) where {T,V,N}

    Av = ForwardDiff.value.(A)                          # values on GPU
    Bv = ForwardDiff.value.(B)
    Cv = NNlib.batched_mul(Av, Bv)                      # batched matmul on values

    # ∂(A·B)/∂xᵢ = A·(∂B/∂xᵢ) + (∂A/∂xᵢ)·B  for each tag i
    dABdx = [
        NNlib.batched_mul(Av, ForwardDiff.partials.(B, i)) +
        NNlib.batched_mul(ForwardDiff.partials.(A, i), Bv)
        for i = 1:N
    ]
    dABdx = ForwardDiff.Partials.(tuple.(dABdx...))
    return eltype(A).(Cv, dABdx)
end
```

Both the value slab and each partial slab are batched-mul'd on GPU; nothing
goes back to host until the final result. The same overload is provided for
`batch_inv!` using the matrix identity ``\partial A^{-1}/\partial x = -A^{-1}(\partial A/\partial x)A^{-1}``.

Production aerosol microphysics derivatives use the hand-linearized Mie
`LinMode()` path and are handed to the analytic RT kernels as ordinary arrays.
The retired Mie `autodiff=true` wrapper is not enabled by these CUDA overloads.

## Allocation pattern

All the per-layer workspace arrays come through `array_type(arch)`:

```julia
# src/CoreRT/tools/rt_helper_functions.jl:114–184
_backend_zeros(FT, AT, dims...) = fill!(similar(AT(FT[]), dims...), zero(FT))
@inline default_matrix(FT, arr_type, dims, nSpec) =
    _backend_zeros(FT, arr_type, dims[1], dims[2], nSpec)
@inline default_J_matrix(FT, arr_type, dims, nSpec) =
    _backend_zeros(FT, arr_type, dims[1], 1, nSpec)

function make_added_layer(RS_type, FT, arr_type, dims, nSpec)
    AddedLayer(
        r⁻⁺ = default_matrix(FT, arr_type, dims, nSpec),
        t⁺⁺ = default_matrix(FT, arr_type, dims, nSpec),
        # … six more matrices …
    )
end
```

When `arr_type == Array`, you get a CPU `Array`. When `arr_type == CuArray`,
you get a `CuArray`. When `arr_type == MtlArray`, you get an `MtlArray`. The
rest of the code reads `r⁻⁺`, `t⁺⁺`, … without knowing or caring which
backend the storage lives on — the kernel methods specialize on the array
type at compile time.

## What is and isn't GPU-friendly today

Honest list of edges:

- **RT solver** (elemental → doubling → interaction) — fully GPU-friendly,
  including linearized variant. This is the hot path.
- **Gas absorption** (HITRAN line-by-line) — the standalone `Absorption`
  module contains a GPU Voigt kernel (`compute_absorption_cross_section.jl:316–332`)
  used by the direct σ(ν, T, p) API.  The RT production pipeline
  (`model_from_parameters` → `rt_run`) computes LBL absorption via
  [AtmosphericAbsorption.jl](https://github.com/RemoteSensingTools/AtmosphericAbsorption.jl),
  which has its own backend dispatch.
- **Mie scattering** — fully automatic via `make_mie_model(...; architecture=GPU())`.
  NAI-2 runs on GPU; PCW falls back to CPU with a one-time warning. The
  retired ForwardDiff Mie path (`autodiff=true`) raises an error; aerosol
  Jacobians use the hand-linearized `LinMode()` Mie path.
  Pass `architecture=CPU()` to force CPU Mie regardless of the RT backend.
- **Thermal emission** — not currently a parallel offline source path.
  Solar-only retrievals are unaffected.
- **Metal batched LU** — has a 32 KiB threadgroup-memory cap on the portable
  KA kernel. Very large stream/Stokes matrices (well beyond typical
  retrievals) may need to fall back to CUDA or CPU until a global-memory
  variant lands.

`config/quickstart.yaml` and the test suite cover the typical scientific
workloads.

## The Julia advantages

A short summary of *why* this is tractable in Julia and not (easily) in
Fortran/C++:

1. **Multiple dispatch + type-generic kernels.** `Stokes_I/IQ/IQU/IQUV`,
   `CPU/GPU/MetalGPU`, `Float32/Float64`, and `ForwardDiff.Dual` all
   specialize the *same kernel source* at compile time. No template
   metaprogramming, no preprocessor, no PIMPL idioms. The function
   `batched_mul(A, B)` resolves to a different machine-code instance for
   each `(eltype, array_type)` combination, with no source duplication.

2. **Package extensions (Julia 1.9+).** The CUDA and Metal backends load only
   when the user loads `CUDA` or `Metal`. Metal is a *weak* dependency; CUDA
   is currently a direct dependency (see above), but CPU runs don't need an
   NVIDIA driver or device. This is impossible in
   Fortran (single static dependency graph) and laborious in C++ (CMake
   toggles + ifdef hell).

3. **No two-language problem.** The kernel body in `doubling.jl` is the
   kernel that runs on the GPU. No Fortran shim, no Python wrapper, no
   pybind11 boundary. Profiling, debugging, and reading the source are
   the same activity at every layer.

4. **Composability via type promotion.** `ForwardDiff.Dual{T,V,N}` works
   inside `batched_mul` because the kernels are element-type generic.
   `CuArray{T}` works because Julia's `Array` interface is implemented by
   any `<:AbstractArray` subtype. `KernelAbstractions.@kernel` works
   because all device APIs implement the same `KernelAbstractions.Backend`
   protocol. None of this required touching the kernel code — it
   *composes*.

## What sets vSmartMOM apart, with file:line evidence

Recapping the eight differentiators from [Concepts/01](01_overview.md), now
with code anchors:

1. **Operator-level analytic linearization.** `src/CoreRT/CoreKernel/{elemental,doubling,interaction}_lin.jl`; chain-rule expansion fused into the elemental step via `get_elem_rt_fused!` / `get_elem_rt_SFI_fused!` (defined in `elemental_fused_lin.jl`, called from `elemental_lin.jl`).
2. **One `@kernel` source for CPU + CUDA + Metal.** [`src/Architectures.jl:33–187`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/src/Architectures.jl#L33-L187); injected backends in [`ext/vSmartMOMCUDAExt.jl:22–28`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/ext/vSmartMOMCUDAExt.jl#L22-L28) and [`ext/vSmartMOMMetalExt.jl:20–23`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/ext/vSmartMOMMetalExt.jl#L20-L23).
3. **Hybrid AD across the GPU boundary.** [`ext/gpu_batched_cuda.jl:255–291`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/ext/gpu_batched_cuda.jl#L255-L291).
4. **Polarization as a type, not a runtime branch.** [`src/Scattering/types.jl:142–215`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/src/Scattering/types.jl#L142-L215).
5. **Optical properties as algebra.** [`src/CoreRT/types.jl:1743–1795`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/src/CoreRT/types.jl#L1743-L1795).
6. **GPU backends as package extensions.** [`ext/vSmartMOMCUDAExt.jl`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/ext/vSmartMOMCUDAExt.jl) (entire file).
7. **Three `(τ, ϖ, Z)` core variables per layer for analytic Jacobians.** [`src/CoreRT/types_lin.jl:330–355`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/src/CoreRT/types_lin.jl#L330-L355).
8. **Exact finite-δ elemental, not the linear limit.** [`src/CoreRT/CoreKernel/elemental.jl:286–331`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/src/CoreRT/CoreKernel/elemental.jl#L286-L331). See [Concepts/04 § Elemental](04_mom_solver.md#Elemental-layer) for the side-by-side equations.

## Code anchors

| Concept | Source |
|---|---|
| Architecture types | [`src/Architectures.jl:33–187`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/src/Architectures.jl#L33-L187) |
| CUDA extension init | [`ext/vSmartMOMCUDAExt.jl:99–129`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/ext/vSmartMOMCUDAExt.jl#L99-L129) |
| Metal extension init | [`ext/vSmartMOMMetalExt.jl:142–174`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/ext/vSmartMOMMetalExt.jl#L142-L174) |
| `@kernel apply_D!` | [`src/CoreRT/CoreKernel/doubling.jl:186–208`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/src/CoreRT/CoreKernel/doubling.jl#L186-L208) |
| `@kernel line_shape_batch!` | [`src/Absorption/compute_absorption_cross_section.jl:316–332`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/src/Absorption/compute_absorption_cross_section.jl#L316-L332) |
| `@kernel mie_coefficients_kernel_ds!` | [`src/Scattering/gpu_mie_kernels.jl:30–158`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/src/Scattering/gpu_mie_kernels.jl#L30-L158) |
| `@kernel _batched_inv_lu_par_kernel!` (Metal) | [`src/CoreRT/tools/ka_batched_kernels.jl:128–206`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/src/CoreRT/tools/ka_batched_kernels.jl#L128-L206) |
| Batched matmul CPU | [`src/CoreRT/tools/cpu_batched.jl:49–100`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/src/CoreRT/tools/cpu_batched.jl#L49-L100) |
| Batched matmul CUDA (CUBLAS) | [`ext/gpu_batched_cuda.jl:225–239`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/ext/gpu_batched_cuda.jl#L225-L239) |
| Batched matmul Metal (KA) | [`ext/vSmartMOMMetalExt.jl:102–141`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/ext/vSmartMOMMetalExt.jl#L102-L141) |
| ForwardDiff.Dual through CUDA matmul | [`ext/gpu_batched_cuda.jl:255–291`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/ext/gpu_batched_cuda.jl#L255-L291) |
| Allocation by architecture | [`src/CoreRT/tools/rt_helper_functions.jl:114–184`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/src/CoreRT/tools/rt_helper_functions.jl#L114-L184) |

## Hands-on tutorials

Runnable examples with Plotly figures:

- [Hybrid AD across the GPU boundary](../tutorials/Tutorial_HybridAD.md)

## References

- [`KernelAbstractions.jl`](https://juliagpu.github.io/KernelAbstractions.jl), [`CUDA.jl`](https://juliagpu.org/cuda/), [`Metal.jl`](https://github.com/JuliaGPU/Metal.jl), [`NNlib.jl`](https://github.com/FluxML/NNlib.jl), [`ForwardDiff.jl`](https://github.com/JuliaDiff/ForwardDiff.jl) — see `Project.toml`.
- Bezanson et al. (2017), *Julia: A fresh approach to numerical computing*, SIAM Review **59**:65. (Multiple-dispatch motivation.)
- Crib sheet: `docs/dev_notes/theory_references.md` §J.
- User Guide: [Run on GPU](../gpu.md) for the runnable how-to.
