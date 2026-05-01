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
