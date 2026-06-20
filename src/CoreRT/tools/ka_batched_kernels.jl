#=

Portable KernelAbstractions batched linear-algebra kernels.

These kernels are intentionally small-matrix focused: vSmartMOM's adding-
doubling core applies the same Nstream x Nstream matrix operation across many
spectral points. Vendor BLAS remains preferred when available (CUDA/CUBLAS),
but these routines give non-CUDA backends such as Metal a backend-native path.

=#

"""
    ka_batched_mul(A, B, backend)

Return the batched product `C[:, :, k] = A[:, :, k] * B[:, :, k]` using a
portable KernelAbstractions kernel on `backend`.
"""
function ka_batched_mul(A::AbstractArray{FT,3}, B::AbstractArray{FT,3}, backend) where {FT}
    @assert size(A, 3) == size(B, 3) "batch dim mismatch: $(size(A,3)) vs $(size(B,3))"
    @assert size(A, 2) == size(B, 1) "inner dim mismatch: $(size(A,2)) vs $(size(B,1))"
    C = similar(A, FT, (size(A, 1), size(B, 2), size(A, 3)))
    return ka_batched_mul!(C, A, B, backend)
end

"""
    ka_batched_mul!(C, A, B, backend)

Fill `C[:, :, k] = A[:, :, k] * B[:, :, k]` using a portable
KernelAbstractions kernel on `backend`.
"""
function ka_batched_mul!(C::AbstractArray{FT,3},
                         A::AbstractArray{FT,3},
                         B::AbstractArray{FT,3},
                         backend) where {FT}
    @assert size(C) == (size(A, 1), size(B, 2), size(A, 3))
    kernel! = _batched_mul_kernel!(backend)
    kernel!(C, A, B, Val(size(A, 2)); ndrange=size(C))
    KernelAbstractions.synchronize(backend)
    return C
end

"""
    _batched_mul_kernel!(C, A, B, Val(K))

KernelAbstractions device kernel for batched matrix multiplication. Each
workitem owns one `(i, j, k)` output element and accumulates
`C[i, j, k] = sum(A[i, l, k] * B[l, j, k], l=1:K)` in a scalar register.
The inputs are read-only and `K` is carried as a `Val` so the inner loop bound
is compile-time constant for the generated kernel.
"""
@kernel function _batched_mul_kernel!(C, @Const(A), @Const(B), ::Val{K}) where {K}
    i, j, k = @index(Global, NTuple)
    s = zero(eltype(C))
    @inbounds for l in 1:K
        s += A[i, l, k] * B[l, j, k]
    end
    @inbounds C[i, j, k] = s
end

"""
    ka_batch_inv_lu!(X, A, backend)

Fill `X[:, :, k] = inv(A[:, :, k])` using a portable KernelAbstractions LU
kernel with partial pivoting. One workgroup handles one matrix and uses
`N = size(A, 1)` workitems.

Set `max_localmem_bytes` for backends with fixed threadgroup/local-memory
budgets. The kernel uses two `N x N` local-memory arrays plus a length-`N`
`Int32` pivot array.
"""
function ka_batch_inv_lu!(X::AbstractArray{FT,3},
                          A::AbstractArray{FT,3},
                          backend;
                          max_localmem_bytes=nothing) where {FT}
    @assert size(A, 1) == size(A, 2) "batched inverse requires square matrices"
    @assert size(X) == size(A) "output size $(size(X)) must match input size $(size(A))"
    N = size(A, 1)
    check_ka_batch_inv_localmem(FT, N, max_localmem_bytes)
    batch = size(A, 3)
    kernel! = _batched_inv_lu_par_kernel!(backend, N)
    kernel!(X, A, Val(N); ndrange=(N * batch,))
    KernelAbstractions.synchronize(backend)
    return X
end

"""
    ka_batch_inv_localmem_bytes(FT, N)

Return the local-memory bytes required per workgroup by
`ka_batch_inv_lu!` for element type `FT` and matrix dimension `N`.
"""
ka_batch_inv_localmem_bytes(::Type{FT}, N::Integer) where {FT} =
    2 * N * N * sizeof(FT) + N * sizeof(Int32)

function check_ka_batch_inv_localmem(::Type{FT}, N::Integer, max_localmem_bytes) where {FT}
    max_localmem_bytes === nothing && return nothing

    required = ka_batch_inv_localmem_bytes(FT, N)
    required <= max_localmem_bytes && return nothing

    throw(ArgumentError(
        "portable batched inverse for N=$(N), eltype=$(FT) requires $(required) " *
        "bytes of local memory per workgroup, exceeding the backend limit of " *
        "$(max_localmem_bytes) bytes. Reduce the stream/Stokes dimension or use " *
        "a backend with a larger local-memory budget."
    ))
end

"""
    _batched_inv_lu_par_kernel!(X, A, Val(N))

KernelAbstractions device kernel for small batched matrix inversion using LU
factorization with partial pivoting. One workgroup handles one spectral
matrix `A[:, :, k]`, each local workitem owns one row/RHS column, and the
shared `LU`, `piv`, and `work` arrays live in kernel local memory. The kernel
writes `X[:, :, k] = inv(A[:, :, k])` without heap allocation.
"""
@kernel function _batched_inv_lu_par_kernel!(X, @Const(A), ::Val{N}) where {N}
    k = @index(Group, Linear)
    tid = @index(Local, Linear)

    LU   = @localmem eltype(X) (N, N)
    piv  = @localmem Int32 (N,)
    work = @localmem eltype(X) (N, N)

    @inbounds for col in 1:N
        LU[tid, col] = A[tid, col, k]
        work[tid, col] = zero(eltype(X))
    end
    @inbounds piv[tid] = Int32(tid)
    @synchronize()

    @inbounds for p in 1:N
        if tid == 1
            max_val = abs(LU[p, p])
            max_row = p
            for r in (p + 1):N
                v = abs(LU[r, p])
                if v > max_val
                    max_val = v
                    max_row = r
                end
            end
            if max_row != p
                for col in 1:N
                    tmp = LU[p, col]
                    LU[p, col] = LU[max_row, col]
                    LU[max_row, col] = tmp
                end
                tmp_p = piv[p]
                piv[p] = piv[max_row]
                piv[max_row] = tmp_p
            end
        end
        @synchronize()

        if tid > p
            LU[tid, p] /= LU[p, p]
        end
        @synchronize()

        if tid > p
            factor = LU[tid, p]
            for col in (p + 1):N
                LU[tid, col] -= factor * LU[p, col]
            end
        end
        @synchronize()
    end

    # Each workitem solves one column of the inverse.
    @inbounds for i in 1:N
        work[i, tid] = (piv[i] == Int32(tid)) ? one(eltype(X)) : zero(eltype(X))
    end

    @inbounds for i in 2:N
        s = zero(eltype(X))
        for j in 1:(i - 1)
            s += LU[i, j] * work[j, tid]
        end
        work[i, tid] -= s
    end

    @inbounds for i in N:-1:1
        s = zero(eltype(X))
        for j in (i + 1):N
            s += LU[i, j] * work[j, tid]
        end
        work[i, tid] = (work[i, tid] - s) / LU[i, i]
    end
    @synchronize()

    @inbounds for row in 1:N
        X[row, tid, k] = work[row, tid]
    end
end

#=
================================================================================
Fused geometric-progression solve for the matrix-operator-method (MOM)
adding-doubling core.

PHYSICS — why this term exists
------------------------------
The MOM combines two atmospheric layers by accounting for the *infinite* series
of internal reflections bouncing between them. Grant & Hunt's central result,
restated in Sanghavi et al. (2014, JQSRT 133:412–433, Eq. 1), is that this
series collapses into a single matrix inverse:

        E + X + X² + X³ + ⋯ = (E − X)⁻¹,

where `E` is the identity and `X` is the supermatrix for one pair of consecutive
reflections. Summing the series exactly (rather than truncating) is what lets
MOM handle thick, weakly-absorbing layers where the reflection series converges
slowly — its key advantage over many other RT methods.

When a homogeneous layer is *doubled* (combined with an identical copy), the
relevant pair of reflections is `X = R·R` (down-reflect, up-reflect), so the
factor is `(E − R·R)⁻¹`.

WHICH MOM EQUATIONS — the reused factor `T·(E − R·R)⁻¹`
------------------------------------------------------
The adding/doubling equations (Sanghavi et al. 2014, Eqs. 23–28) for the
transmission, reflection, and internal-source operators all share the *same*
pre-multiplied factor `T·(E − R·R)⁻¹`:

    Tᵈ = T·(E − R·R)⁻¹·T                         (transmission update, Eq. 23)
    Rᵈ = R + T·(E − R·R)⁻¹·R·T                   (reflection update,   Eq. 24)
    Jᵈ = J + T·(E − R·R)⁻¹·(J + R·J)             (source update,       Eq. 27/28)

so the adding-doubling kernel forms this factor **once per doubling step**
(the code's `tt_gp = t⁺⁺·(E − r⁻⁺·r⁻⁺)⁻¹`) and reuses it for the R, T, and J
updates. The matrices are small — `N = Nstream × Nstokes` (e.g. 18) — and the
same operation is applied across all `nSpec` spectral points (a large batch).

WHY THE FUSED KERNEL IS FASTER
------------------------------
Note `tt_gp` only needs `T·(inverse)`, never the bare inverse, so it is a linear
*solve*, not an inversion. The textbook batched path is three vendor calls —
`getrf` (LU) + `getri` (explicit inverse) + `gemm` (× T). For these tiny
matrices that path is dominated by per-call launch/sync overhead and global-
memory round-trips, materialises an explicit inverse we never need (≈ an extra
O(N³) stage), and forming the inverse is less numerically stable than a solve.

The fused kernel does the whole thing in **one launch**: one workgroup per
spectral matrix builds `M = E − R·R`, LU-factorises it with partial pivoting in
local (shared) memory, then **right-solves** `tt_gp·M = T` directly
(`tt_gp = T·U⁻¹·L⁻¹·P`). No explicit inverse, no separate GEMM, no inter-stage
global traffic. On the small RT matrices this is ~3–6× faster than the vendor
path at production batch sizes (`test/benchmarks/batched_fused_v2_benchmark.jl`).
The cost is a `3N²` local-memory tile; for large `N` that exceeds the device
budget and the caller falls back to the vendor inverse-then-multiply path.
================================================================================
=#

"""
    ka_fused_gp_solve_localmem_bytes(FT, N)

Local-memory bytes per workgroup required by [`ka_fused_gp_solve!`]: three
`N×N` tiles (`A`, the factorised `M = I − A·A`, and the in-place RHS/solution
`Bt`) plus a length-`N` `Int32` pivot vector.
"""
ka_fused_gp_solve_localmem_bytes(::Type{FT}, N::Integer) where {FT} =
    3 * N * N * sizeof(FT) + N * sizeof(Int32)

"""
    ka_fused_gp_solve!(X, A, B, backend; max_localmem_bytes=nothing)

Compute the adding-doubling geometric-progression factor

    X[:, :, k] = B[:, :, k] * (I − A[:, :, k] * A[:, :, k])⁻¹    (= t⁺⁺·(E − r⁻⁺·r⁻⁺)⁻¹)

in a *single* fused KernelAbstractions launch. One workgroup handles one
spectral matrix and `N = size(A, 1)` workitems cooperate in local memory: build
`M = I − A·A`, LU-factorise it with partial pivoting, then **right-solve**
`X·M = B` directly (`X = B·U⁻¹·L⁻¹·P`). The explicit inverse is never formed and
there is no separate batched GEMM — for the small Nstream×Nstokes matrices of the
RT core this is several× faster than the vendor `getrf + getri + gemm` path
(see `test/benchmarks/batched_fused_v2_benchmark.jl`).

Pass `max_localmem_bytes` to assert the `3N²`-tile footprint fits the backend's
local-memory budget; callers should fall back to the inverse-then-multiply path
when it does not (large `N`).
"""
function ka_fused_gp_solve!(X::AbstractArray{FT,3},
                            A::AbstractArray{FT,3},
                            B::AbstractArray{FT,3},
                            backend;
                            max_localmem_bytes=nothing) where {FT}
    @assert size(A, 1) == size(A, 2) "fused GP solve requires square A"
    @assert size(X) == size(A) == size(B) "X, A, B must share dimensions"
    N = size(A, 1)
    if max_localmem_bytes !== nothing
        required = ka_fused_gp_solve_localmem_bytes(FT, N)
        required <= max_localmem_bytes || throw(ArgumentError(
            "fused GP solve for N=$(N), eltype=$(FT) needs $(required) bytes of " *
            "local memory per workgroup, exceeding the backend limit of " *
            "$(max_localmem_bytes) bytes; fall back to the inverse path."))
    end
    batch = size(A, 3)
    kernel! = _fused_gp_solve_kernel!(backend, N)
    kernel!(X, A, B, Val(N); ndrange=(N * batch,))
    KernelAbstractions.synchronize(backend)
    return X
end

"""
    _fused_gp_solve_kernel!(X, A, B, Val(N))

Device kernel for [`ka_fused_gp_solve!`]. One workgroup per spectral matrix,
one workitem per row. Local tiles `At`/`M`/`Bt` hold the inputs, the LU factors,
and the in-place solution; `piv` holds the pivot permutation. After the
cooperative LU, each workitem solves its own row independently (`X = B·U⁻¹·L⁻¹`),
then scatters through the column pivot.
"""
@kernel function _fused_gp_solve_kernel!(X, @Const(A), @Const(B), ::Val{N}) where {N}
    k   = @index(Group, Linear)
    tid = @index(Local, Linear)

    At  = @localmem eltype(X) (N, N)
    M   = @localmem eltype(X) (N, N)
    piv = @localmem Int32 (N,)
    Bt  = @localmem eltype(X) (N, N)

    # Load A and the RHS B into local memory.
    @inbounds for j in 1:N
        At[tid, j] = A[tid, j, k]
        Bt[tid, j] = B[tid, j, k]
    end
    @synchronize()

    # M = I − A·A (row `tid`).
    @inbounds for j in 1:N
        s = zero(eltype(X))
        for l in 1:N
            s += At[tid, l] * At[l, j]
        end
        M[tid, j] = ((tid == j) ? one(eltype(X)) : zero(eltype(X))) - s
    end
    @inbounds piv[tid] = Int32(tid)
    @synchronize()

    # LU factorisation of M with partial pivoting (workitems cooperate).
    @inbounds for p in 1:N
        if tid == 1
            max_val = abs(M[p, p])
            max_row = p
            for r in (p + 1):N
                v = abs(M[r, p])
                if v > max_val
                    max_val = v
                    max_row = r
                end
            end
            if max_row != p
                for col in 1:N
                    tmp = M[p, col]
                    M[p, col] = M[max_row, col]
                    M[max_row, col] = tmp
                end
                tmp_p = piv[p]
                piv[p] = piv[max_row]
                piv[max_row] = tmp_p
            end
        end
        @synchronize()
        if tid > p
            M[tid, p] /= M[p, p]
        end
        @synchronize()
        if tid > p
            factor = M[tid, p]
            for col in (p + 1):N
                M[tid, col] -= factor * M[p, col]
            end
        end
        @synchronize()
    end

    # Right-side direct solve X = B·U⁻¹·L⁻¹·P, in place on row `tid` of Bt.
    # (a) forward solve Z·U = B
    @inbounds for j in 1:N
        s = zero(eltype(X))
        for l in 1:(j - 1)
            s += Bt[tid, l] * M[l, j]
        end
        Bt[tid, j] = (Bt[tid, j] - s) / M[j, j]
    end
    # (b) backward solve W·L = Z (unit lower diagonal)
    @inbounds for j in (N - 1):-1:1
        s = zero(eltype(X))
        for l in (j + 1):N
            s += Bt[tid, l] * M[l, j]
        end
        Bt[tid, j] -= s
    end
    # (c) apply the column permutation
    @inbounds for j in 1:N
        X[tid, piv[j], k] = Bt[tid, j]
    end
end

"""
    _gp_fused_localmem_limit(backend)

Per-workgroup local-memory budget (bytes) usable by [`ka_fused_gp_solve!`] on
`backend`. Conservative portable default (Metal-class threadgroup memory); the
CUDA extension raises it to the 48 KiB static shared-memory limit.
"""
_gp_fused_localmem_limit(::KernelAbstractions.Backend) = 32 * 1024

"""
    _FUSED_GP_ENABLED

Runtime kill-switch for the fused geometric-progression solve. Default `true`.
Set `_FUSED_GP_ENABLED[] = false` to force the portable inverse-then-multiply
fallback (used for A/B validation and as a safety override).
"""
const _FUSED_GP_ENABLED = Ref(true)

"""
    _use_fused_gp(X) -> Bool

Whether the fused geometric-progression solve ([`ka_fused_gp_solve!`]) should be
used for output array `X`. True on GPU backends when [`_FUSED_GP_ENABLED`] is set
and the `3N²` local-memory tile fits [`_gp_fused_localmem_limit`]; false on the
CPU (the threaded BLAS inverse path is already efficient there) and when the tile
is too large for the device (the caller falls back to inverse-then-multiply).
"""
@inline function _use_fused_gp(X::AbstractArray{FT,3}) where {FT}
    _FUSED_GP_ENABLED[] || return false
    backend = KernelAbstractions.get_backend(X)
    backend isa KernelAbstractions.CPU && return false
    return ka_fused_gp_solve_localmem_bytes(FT, size(X, 1)) <= _gp_fused_localmem_limit(backend)
end
