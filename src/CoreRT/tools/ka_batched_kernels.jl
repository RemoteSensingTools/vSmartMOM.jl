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
    # Host-side sync only where it is semantically needed: KA's CPU backend
    # runs kernels as async tasks, so consumers need the barrier. On GPU the
    # pipeline is single-stream — stream ordering already sequences every
    # consumer — and the eager sync both wastes a full-stream round trip per
    # call AND invalidates CUDA stream capture (graph replay). Phase-0 fix.
    backend isa KernelAbstractions.CPU && KernelAbstractions.synchronize(backend)
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
    # Host-side sync only where it is semantically needed: KA's CPU backend
    # runs kernels as async tasks, so consumers need the barrier. On GPU the
    # pipeline is single-stream — stream ordering already sequences every
    # consumer — and the eager sync both wastes a full-stream round trip per
    # call AND invalidates CUDA stream capture (graph replay). Phase-0 fix.
    backend isa KernelAbstractions.CPU && KernelAbstractions.synchronize(backend)
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
    # Host-side sync only where it is semantically needed: KA's CPU backend
    # runs kernels as async tasks, so consumers need the barrier. On GPU the
    # pipeline is single-stream — stream ordering already sequences every
    # consumer — and the eager sync both wastes a full-stream round trip per
    # call AND invalidates CUDA stream capture (graph replay). Phase-0 fix.
    backend isa KernelAbstractions.CPU && KernelAbstractions.synchronize(backend)
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

# ------------------------------------------------------------------------------
# Two-matrix variant for the MOM adding (interaction) step.
#
# The composite/added-layer combine (Sanghavi et al. 2014, Eqs. 23–28) needs the
# same `B·(I − A₁·A₂)⁻¹` factor, but with two DISTINCT reflection matrices, e.g.
#   T01_inv = T⁻⁻ · (E − r⁻⁺·R⁺⁻)⁻¹      and   T21_inv = t⁺⁺ · (E − R⁺⁻·r⁻⁺)⁻¹,
# each reused across the layer's R/T/J updates. As in the doubling case the bare
# inverse is only an intermediate, so it is again a fused LU + right-solve. With
# A₁≠A₂ there is no single tile to cache, so this kernel needs only `2N²` local
# memory (vs `3N²` for the geometric-progression variant) and fits larger `N`.
# ------------------------------------------------------------------------------

"""
    ka_fused_solve_localmem_bytes(FT, N)

Local-memory bytes per workgroup required by [`ka_fused_solve!`]: two `N×N` tiles
(the factorised `M = I − A₁·A₂` and the in-place RHS/solution `Bt`) plus a
length-`N` `Int32` pivot vector.
"""
ka_fused_solve_localmem_bytes(::Type{FT}, N::Integer) where {FT} =
    2 * N * N * sizeof(FT) + N * sizeof(Int32)

"""
    ka_fused_solve!(X, A1, A2, B, backend; max_localmem_bytes=nothing)

Compute `X[:, :, k] = B[:, :, k] * (I − A1[:, :, k] * A2[:, :, k])⁻¹` for two
distinct batched matrices `A1`, `A2` in a single fused KernelAbstractions launch
(build `M = I − A1·A2`, LU-factorise with partial pivoting in local memory, then
right-solve `X·M = B` directly — no explicit inverse, no separate GEMM). This is
the adding-step (`interaction!`) analogue of [`ka_fused_gp_solve!`]; see the
physics header above and `test/benchmarks/batched_fused_v2_benchmark.jl`.
"""
function ka_fused_solve!(X::AbstractArray{FT,3},
                         A1::AbstractArray{FT,3},
                         A2::AbstractArray{FT,3},
                         B::AbstractArray{FT,3},
                         backend;
                         max_localmem_bytes=nothing) where {FT}
    @assert size(A1, 1) == size(A1, 2) "fused solve requires square A1·A2"
    @assert size(X) == size(A1) == size(A2) == size(B) "X, A1, A2, B must share dimensions"
    N = size(A1, 1)
    if max_localmem_bytes !== nothing
        required = ka_fused_solve_localmem_bytes(FT, N)
        required <= max_localmem_bytes || throw(ArgumentError(
            "fused solve for N=$(N), eltype=$(FT) needs $(required) bytes of local " *
            "memory per workgroup, exceeding the backend limit of $(max_localmem_bytes) " *
            "bytes; fall back to the inverse path."))
    end
    batch = size(A1, 3)
    kernel! = _fused_solve_kernel!(backend, N)
    kernel!(X, A1, A2, B, Val(N); ndrange=(N * batch,))
    # Host-side sync only where it is semantically needed: KA's CPU backend
    # runs kernels as async tasks, so consumers need the barrier. On GPU the
    # pipeline is single-stream — stream ordering already sequences every
    # consumer — and the eager sync both wastes a full-stream round trip per
    # call AND invalidates CUDA stream capture (graph replay). Phase-0 fix.
    backend isa KernelAbstractions.CPU && KernelAbstractions.synchronize(backend)
    return X
end

"""
    _fused_solve_kernel!(X, A1, A2, B, Val(N))

Device kernel for [`ka_fused_solve!`]. Identical to [`_fused_gp_solve_kernel!`]
except `M = I − A1·A2` is built from two distinct inputs (no cached `A` tile).
"""
@kernel function _fused_solve_kernel!(X, @Const(A1), @Const(A2), @Const(B), ::Val{N}) where {N}
    k   = @index(Group, Linear)
    tid = @index(Local, Linear)

    M   = @localmem eltype(X) (N, N)
    piv = @localmem Int32 (N,)
    Bt  = @localmem eltype(X) (N, N)

    # M = I − A1·A2 (row `tid`), and load the RHS B.
    @inbounds for j in 1:N
        s = zero(eltype(X))
        for l in 1:N
            s += A1[tid, l, k] * A2[l, j, k]
        end
        M[tid, j] = ((tid == j) ? one(eltype(X)) : zero(eltype(X))) - s
        Bt[tid, j] = B[tid, j, k]
    end
    @inbounds piv[tid] = Int32(tid)
    @synchronize()

    # LU factorisation of M with partial pivoting.
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

    # Right-side direct solve X = B·U⁻¹·L⁻¹·P (in place on row `tid` of Bt).
    @inbounds for j in 1:N
        s = zero(eltype(X))
        for l in 1:(j - 1)
            s += Bt[tid, l] * M[l, j]
        end
        Bt[tid, j] = (Bt[tid, j] - s) / M[j, j]
    end
    @inbounds for j in (N - 1):-1:1
        s = zero(eltype(X))
        for l in (j + 1):N
            s += Bt[tid, l] * M[l, j]
        end
        Bt[tid, j] -= s
    end
    @inbounds for j in 1:N
        X[tid, piv[j], k] = Bt[tid, j]
    end
end

"""
    _use_fused_solve(X) -> Bool

Whether the two-matrix fused solve ([`ka_fused_solve!`]) should be used for output
array `X`: same policy as [`_use_fused_gp`] but with the smaller `2N²` tile.
"""
@inline function _use_fused_solve(X::AbstractArray{FT,3}) where {FT}
    _FUSED_GP_ENABLED[] || return false
    backend = KernelAbstractions.get_backend(X)
    backend isa KernelAbstractions.CPU && return false
    return ka_fused_solve_localmem_bytes(FT, size(X, 1)) <= _gp_fused_localmem_limit(backend)
end

# ==============================================================================
# Phase-1 fused doubling-step kernels (perf/fused-adding-kernels).
#
# The Phase-0 in-place doubling updates still cost ~12 (rt) and ~18 (source)
# device events per call — each batched product is a cuBLAS pointer-array
# gemm (3 pageable pointer-table HtoD uploads + 1 kernel for the tiny
# Nstream×Nstokes tiles) plus separate broadcast kernels, and the cuBLAS
# nodes make CUDA-graph replay impossible (single-shot graphs). These two
# kernels collapse each update into ONE backend-agnostic KA launch:
#
#   rt:      R ← R + (A·R)·T ;  T ← A·T ;  expk ← expk²
#   source:  j₁± = j₀±·e ;  j₀⁻ += A·(j₁⁻ + R·j₀⁺) ;  j₀⁺ = j₁⁺ + A·(j₀⁺ + R·j₁⁻)
#
# with A = tt_gp, R = r⁻⁺, T = t⁺⁺. The rt kernel uses the associativity
# (A·R)·T[:,j] = A·(R·T[:,j]) so one workitem owns one output COLUMN via two
# mat-vecs against shared tiles — same 3N² local-memory budget as the fused
# GP solve, no fourth tile. The source kernel is row-parallel (vectors).
#
# EQUIVALENCE CONTRACT: same math, different reduction order than cuBLAS —
# tolerance-equivalent (NOT bitwise) to the `_bmm!` path on every backend;
# validated in test/test_fused_doubling.jl. Kill-switches restore the
# `_bmm!` path exactly.
# ==============================================================================

"Local-memory bytes for [`ka_fused_doubling_rt!`]: three N×N tiles."
ka_fused_doubling_rt_localmem_bytes(::Type{FT}, N::Integer) where {FT} =
    3 * N * N * sizeof(FT)

"Local-memory bytes for [`ka_fused_doubling_source!`]: two N×N tiles + six N-vectors."
ka_fused_doubling_source_localmem_bytes(::Type{FT}, N::Integer) where {FT} =
    (2 * N * N + 6 * N) * sizeof(FT)

"Runtime kill-switch for the fused doubling-step kernels. Default `true`."
const _FUSED_DOUBLING_ENABLED = Ref(true)

"""
    _use_fused_doubling_rt(X) / _use_fused_doubling_source(X) -> Bool

Gate for the fused doubling kernels: GPU backends only (the CPU threaded
BLAS path is already efficient), kill-switch on, and the local-memory tile
fits the device budget. `X` is any of the update's matrix arguments.
"""
@inline function _use_fused_doubling_rt(X::AbstractArray{FT,3}) where {FT}
    _FUSED_DOUBLING_ENABLED[] || return false
    backend = KernelAbstractions.get_backend(X)
    backend isa KernelAbstractions.CPU && return false
    return ka_fused_doubling_rt_localmem_bytes(FT, size(X, 1)) <= _gp_fused_localmem_limit(backend)
end
@inline function _use_fused_doubling_source(X::AbstractArray{FT,3}) where {FT}
    _FUSED_DOUBLING_ENABLED[] || return false
    backend = KernelAbstractions.get_backend(X)
    backend isa KernelAbstractions.CPU && return false
    return ka_fused_doubling_source_localmem_bytes(FT, size(X, 1)) <= _gp_fused_localmem_limit(backend)
end

"""
    ka_fused_doubling_rt!(r⁻⁺, t⁺⁺, tt_gp, expk, backend)

One-launch doubling R/T update:

    r⁻⁺[:,:,k] ← r⁻⁺[:,:,k] + (tt_gp·r⁻⁺)·t⁺⁺   (per spectral point k)
    t⁺⁺[:,:,k] ← tt_gp[:,:,k]·t⁺⁺[:,:,k]
    expk[k]    ← expk[k]²

One workgroup per k, one workitem per output column; shared tiles carry
the (pre-update) A/R/T so global writes never race the reads.
"""
function ka_fused_doubling_rt!(r⁻⁺::AbstractArray{FT,3}, t⁺⁺::AbstractArray{FT,3},
                               tt_gp::AbstractArray{FT,3}, expk::AbstractArray{FT,1},
                               backend) where {FT}
    N = size(r⁻⁺, 1)
    @assert size(r⁻⁺) == size(t⁺⁺) == size(tt_gp)
    @assert size(r⁻⁺, 2) == N
    batch = size(r⁻⁺, 3)
    kernel! = _fused_doubling_rt_kernel!(backend, N)
    kernel!(r⁻⁺, t⁺⁺, tt_gp, expk, Val(N); ndrange=(N * batch,))
    backend isa KernelAbstractions.CPU && KernelAbstractions.synchronize(backend)
    return nothing
end

@kernel function _fused_doubling_rt_kernel!(r⁻⁺, t⁺⁺, @Const(tt_gp), expk, ::Val{N}) where {N}
    k   = @index(Group, Linear)
    tid = @index(Local, Linear)

    Ash = @localmem eltype(r⁻⁺) (N, N)
    Rsh = @localmem eltype(r⁻⁺) (N, N)
    Tsh = @localmem eltype(r⁻⁺) (N, N)
    u   = @private eltype(r⁻⁺) (N,)

    @inbounds for j in 1:N
        Ash[tid, j] = tt_gp[tid, j, k]
        Rsh[tid, j] = r⁻⁺[tid, j, k]
        Tsh[tid, j] = t⁺⁺[tid, j, k]
    end
    @synchronize()

    # Doubling R/T adding equations for two IDENTICAL sub-slabs
    # (Sanghavi et al. 2014, Eqs. 23–24, collapsed; see doubling.jl header):
    #
    #     r⁻⁺ ← r⁻⁺ + tt_gp · r⁻⁺ · t⁺⁺        with tt_gp = t⁺⁺·(E − r⁻⁺r⁻⁺)⁻¹
    #     t⁺⁺ ← tt_gp · t⁺⁺
    #
    # This workitem owns output COLUMN `tid`. Associativity turns the
    # triple product into two mat-vecs:  (tt_gp·r⁻⁺)·t⁺⁺[:,tid] =
    # tt_gp·(r⁻⁺·t⁺⁺[:,tid]), so no N×N intermediate is materialised.
    @inbounds begin
        # u = r⁻⁺ · t⁺⁺[:,tid]                  (first mat-vec)
        for l in 1:N
            acc = zero(eltype(r⁻⁺))
            for m in 1:N
                acc += Rsh[l, m] * Tsh[m, tid]
            end
            u[l] = acc
        end
        # r⁻⁺[:,tid] ← r⁻⁺[:,tid] + tt_gp·u    (second mat-vec, accumulate)
        # t⁺⁺[:,tid] ← tt_gp · t⁺⁺[:,tid]      (independent mat-vec, fused loop)
        for l in 1:N
            acc_r = zero(eltype(r⁻⁺))
            acc_t = zero(eltype(r⁻⁺))
            for m in 1:N
                acc_r += Ash[l, m] * u[m]
                acc_t += Ash[l, m] * Tsh[m, tid]
            end
            r⁻⁺[l, tid, k] = Rsh[l, tid] + acc_r
            t⁺⁺[l, tid, k] = acc_t
        end
        # expk ← expk² — the attenuation factor tracks the doubled optical
        # thickness (exp(-2δτ/μ₀) = exp(-δτ/μ₀)²); one workitem per k.
        if tid == 1
            expk[k] = expk[k]^2
        end
    end
end

"""
    ka_fused_doubling_source!(j₀⁺, j₀⁻, j₁⁺, j₁⁻, r⁻⁺, tt_gp, expk, backend)

One-launch doubling source update (per spectral point k, e = expk[k]):

    j₁± = j₀±·e
    j₀⁻ ← j₀⁻ + tt_gp·(j₁⁻ + r⁻⁺·j₀⁺)
    j₀⁺ ← j₁⁺ + tt_gp·(j₀⁺ + r⁻⁺·j₁⁻)     (pre-update j₀⁺ on the RHS)

j₁± are also written to their global workspace buffers (same contract as
the `_bmm!` path). One workgroup per k, one workitem per vector row.
"""
function ka_fused_doubling_source!(j₀⁺::AbstractArray{FT,3}, j₀⁻::AbstractArray{FT,3},
                                   j₁⁺::AbstractArray{FT,3}, j₁⁻::AbstractArray{FT,3},
                                   r⁻⁺::AbstractArray{FT,3}, tt_gp::AbstractArray{FT,3},
                                   expk::AbstractArray{FT,1}, backend) where {FT}
    N = size(r⁻⁺, 1)
    batch = size(r⁻⁺, 3)
    kernel! = _fused_doubling_source_kernel!(backend, N)
    kernel!(j₀⁺, j₀⁻, j₁⁺, j₁⁻, r⁻⁺, tt_gp, expk, Val(N); ndrange=(N * batch,))
    backend isa KernelAbstractions.CPU && KernelAbstractions.synchronize(backend)
    return nothing
end

@kernel function _fused_doubling_source_kernel!(j₀⁺, j₀⁻, j₁⁺, j₁⁻,
                                                @Const(r⁻⁺), @Const(tt_gp),
                                                @Const(expk), ::Val{N}) where {N}
    k   = @index(Group, Linear)
    tid = @index(Local, Linear)

    Ash = @localmem eltype(j₀⁺) (N, N)
    Rsh = @localmem eltype(j₀⁺) (N, N)
    xp  = @localmem eltype(j₀⁺) (N,)   # pre-update j₀⁺
    xm  = @localmem eltype(j₀⁺) (N,)   # pre-update j₀⁻
    w1  = @localmem eltype(j₀⁺) (N,)   # j₁⁻ + R·j₀⁺
    w2  = @localmem eltype(j₀⁺) (N,)   # j₀⁺ + R·j₁⁻

    @inbounds begin
        for j in 1:N
            Ash[tid, j] = tt_gp[tid, j, k]
            Rsh[tid, j] = r⁻⁺[tid, j, k]
        end
        xp[tid] = j₀⁺[tid, 1, k]
        xm[tid] = j₀⁻[tid, 1, k]
    end
    @synchronize()

    # Doubling source-cascade equations (Sanghavi et al. 2014, Eqs. 27–28,
    # restated for identical sub-slabs as Eqs. 8 of Sanghavi & Frankenberg
    # 2023; see doubling_source_update!'s docstring):
    #
    #     j₁± = j₀± · expk            (beam-attenuated copies, expk = e^{-δτ/μ₀})
    #     j₀⁻ ← j₀⁻ + tt_gp·(j₁⁻ + r⁻⁺·j₀⁺)
    #     j₀⁺ ← j₁⁺ + tt_gp·(j₀⁺ + r⁻⁺·j₁⁻)    (pre-update j₀⁺ on the RHS)
    #
    # Stage 1 — workitem `tid` builds row `tid` of both inner vectors:
    #     w1 = j₁⁻ + r⁻⁺·j₀⁺
    #     w2 = j₀⁺ + r⁻⁺·j₁⁻
    # (NOTE: the scalar e is re-read from expk[k] in each stage — scalar
    # locals do not survive KA's @synchronize loop-splitting.)
    @inbounds begin
        e = expk[k]
        rj0p = zero(eltype(j₀⁺))   # (r⁻⁺·j₀⁺)[tid]
        rj1m = zero(eltype(j₀⁺))   # (r⁻⁺·j₁⁻)[tid], with j₁⁻[m] = xm[m]·e
        for m in 1:N
            rj0p += Rsh[tid, m] * xp[m]
            rj1m += Rsh[tid, m] * (xm[m] * e)
        end
        w1[tid] = xm[tid] * e + rj0p     # j₁⁻[tid] + (r⁻⁺·j₀⁺)[tid]
        w2[tid] = xp[tid] + rj1m         # j₀⁺[tid] + (r⁻⁺·j₁⁻)[tid]
    end
    @synchronize()

    # Stage 2 — apply the geometric-progression factor and assemble:
    #     j₀⁻[tid] = j₀⁻[tid] + (tt_gp·w1)[tid]
    #     j₀⁺[tid] = j₁⁺[tid] + (tt_gp·w2)[tid]
    @inbounds begin
        e = expk[k]
        tw1 = zero(eltype(j₀⁺))    # (tt_gp·w1)[tid]
        tw2 = zero(eltype(j₀⁺))    # (tt_gp·w2)[tid]
        for m in 1:N
            tw1 += Ash[tid, m] * w1[m]
            tw2 += Ash[tid, m] * w2[m]
        end
        j₁⁺[tid, 1, k] = xp[tid] * e     # persisted workspace, same contract
        j₁⁻[tid, 1, k] = xm[tid] * e     # as the _bmm! path
        j₀⁻[tid, 1, k] = xm[tid] + tw1
        j₀⁺[tid, 1, k] = xp[tid] * e + tw2
    end
end
