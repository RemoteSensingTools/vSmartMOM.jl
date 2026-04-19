#!/usr/bin/env julia
#
# Benchmark: Custom KernelAbstractions.jl batched matrix inversion vs cuBLAS getrf+getri
# Tests three approaches:
#   1. cuBLAS getrf + getri (two kernel launches)
#   2. Custom Gauss-Jordan elimination (single kernel, N threads)
#   3. Custom LU decomposition + triangular solves (single kernel, N threads)
#
# Usage:  julia --project=test test/benchmarks/batched_inv_benchmark.jl
#

using CUDA
using KernelAbstractions
using Printf
using LinearAlgebra

@assert CUDA.functional() "CUDA must be available to run this benchmark"

# ============================================================================
# Gauss-Jordan elimination with partial pivoting (N threads per workgroup)
# Shared memory: [A | I] augmented matrix, N × 2N
# ============================================================================
@kernel function batched_inv_gj_kernel!(X, A, ::Val{N}) where {N}
    k = @index(Group, Linear)
    row = @index(Local, Linear)

    aug = @localmem eltype(X) (N, 2 * N)

    # Load [A | I]
    @inbounds for col in 1:N
        aug[row, col] = A[row, col, k]
        aug[row, N + col] = (row == col) ? one(eltype(X)) : zero(eltype(X))
    end
    @synchronize()

    @inbounds for p in 1:N
        # Thread 1 does pivoting
        if row == 1
            max_val = abs(aug[p, p])
            max_row = p
            for r in (p + 1):N
                v = abs(aug[r, p])
                if v > max_val
                    max_val = v
                    max_row = r
                end
            end
            if max_row != p
                for col in 1:(2 * N)
                    tmp = aug[p, col]
                    aug[p, col] = aug[max_row, col]
                    aug[max_row, col] = tmp
                end
            end
        end
        @synchronize()

        # Scale pivot row
        if row == p
            pivot_val = aug[p, p]
            for col in 1:(2 * N)
                aug[p, col] /= pivot_val
            end
        end
        @synchronize()

        # Eliminate from all other rows
        if row != p
            factor = aug[row, p]
            for col in 1:(2 * N)
                aug[row, col] -= factor * aug[p, col]
            end
        end
        @synchronize()
    end

    @inbounds for col in 1:N
        X[row, col, k] = aug[row, N + col]
    end
end

function ka_batch_inv_gj!(X::CuArray{FT,3}, A::CuArray{FT,3}) where {FT}
    N = size(A, 1)
    batch = size(A, 3)
    backend = CUDABackend()
    batched_inv_gj_kernel!(backend, N)(X, A, Val(N); ndrange=(N * batch,))
    KernelAbstractions.synchronize(backend)
    return X
end

# ============================================================================
# LU decomposition + forward/back substitution (N threads per workgroup)
# Single fused kernel: PA = LU, then solve L(UX) = PI column by column.
# Shared memory: LU matrix (N×N) + permutation (N) + temp column (N)
# Much less shared memory than GJ: N×N + 2N vs N×2N
# ============================================================================
@kernel function batched_inv_lu_kernel!(X, A, ::Val{N}) where {N}
    k = @index(Group, Linear)
    row = @index(Local, Linear)

    # Shared memory for LU (in-place over A), permutation vector, and temp column
    LU  = @localmem eltype(X) (N, N)
    piv = @localmem Int32 (N,)
    col_buf = @localmem eltype(X) (N,)  # for forward/back substitution

    # Load A into shared memory
    @inbounds for col in 1:N
        LU[row, col] = A[row, col, k]
    end
    @inbounds piv[row] = Int32(row)
    @synchronize()

    # ---- LU factorization with partial pivoting ----
    @inbounds for p in 1:N
        # Thread 1 finds pivot and swaps rows
        if row == 1
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
                # Swap rows in LU
                for col in 1:N
                    tmp = LU[p, col]
                    LU[p, col] = LU[max_row, col]
                    LU[max_row, col] = tmp
                end
                # Swap permutation
                tmp_p = piv[p]
                piv[p] = piv[max_row]
                piv[max_row] = tmp_p
            end
        end
        @synchronize()

        # Compute L multipliers for rows below pivot
        if row > p
            LU[row, p] /= LU[p, p]
        end
        @synchronize()

        # Update trailing submatrix: A[i,j] -= L[i,p] * U[p,j] for i,j > p
        if row > p
            factor = LU[row, p]
            for col in (p + 1):N
                LU[row, col] -= factor * LU[p, col]
            end
        end
        @synchronize()
    end

    # ---- Solve for each column of the inverse ----
    # For column `c` of X: solve LU * x = P * e_c
    # Step 1: Forward substitution: L * y = P * e_c
    # Step 2: Back substitution: U * x = y
    @inbounds for c in 1:N

        # Load permuted identity column into col_buf
        # (P * e_c)[row] = 1 if piv[row] == c, else 0
        col_buf[row] = (piv[row] == Int32(c)) ? one(eltype(X)) : zero(eltype(X))
        @synchronize()

        # Forward substitution: y[i] = b[i] - sum(L[i,j]*y[j], j=1..i-1)
        # Must be done sequentially by row (thread 1 drives, others read)
        for i in 2:N
            if row == 1
                s = zero(eltype(X))
                for j in 1:(i - 1)
                    s += LU[i, j] * col_buf[j]
                end
                col_buf[i] -= s
            end
            @synchronize()
        end

        # Back substitution: x[i] = (y[i] - sum(U[i,j]*x[j], j=i+1..N)) / U[i,i]
        for i in N:-1:1
            if row == 1
                s = zero(eltype(X))
                for j in (i + 1):N
                    s += LU[i, j] * col_buf[j]
                end
                col_buf[i] = (col_buf[i] - s) / LU[i, i]
            end
            @synchronize()
        end

        # Write result column
        X[row, c, k] = col_buf[row]
        @synchronize()
    end
end

function ka_batch_inv_lu!(X::CuArray{FT,3}, A::CuArray{FT,3}) where {FT}
    N = size(A, 1)
    batch = size(A, 3)
    backend = CUDABackend()
    batched_inv_lu_kernel!(backend, N)(X, A, Val(N); ndrange=(N * batch,))
    KernelAbstractions.synchronize(backend)
    return X
end

# ============================================================================
# LU with parallel triangular solves (N threads solve N columns simultaneously)
# Each thread handles one column of the identity for forward/back substitution.
# Shared memory: LU (N×N) + piv (N) + work matrix (N×N)
# ============================================================================
@kernel function batched_inv_lu_par_kernel!(X, A, ::Val{N}) where {N}
    k = @index(Group, Linear)
    tid = @index(Local, Linear)  # tid ∈ 1:N, each thread owns one RHS column

    LU   = @localmem eltype(X) (N, N)
    piv  = @localmem Int32 (N,)
    work = @localmem eltype(X) (N, N)  # work[:,c] = solution column c

    # Load A into shared memory (each thread loads one row)
    @inbounds for col in 1:N
        LU[tid, col] = A[tid, col, k]
    end
    @inbounds piv[tid] = Int32(tid)
    @synchronize()

    # ---- LU factorization with partial pivoting (same as above) ----
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

    # ---- Parallel triangular solves: each thread solves for one RHS column ----
    # Thread `tid` solves: L * y = P * e_{tid}, then U * x = y
    # `c` is the column index this thread is responsible for
    c = tid

    # Initialize work[:,c] = P * e_c
    @inbounds for i in 1:N
        work[i, c] = (piv[i] == Int32(c)) ? one(eltype(X)) : zero(eltype(X))
    end
    @synchronize()

    # Forward substitution (sequential over rows, but each thread does its own column)
    @inbounds for i in 2:N
        s = zero(eltype(X))
        for j in 1:(i - 1)
            s += LU[i, j] * work[j, c]
        end
        work[i, c] -= s
        # No sync needed: each thread only reads/writes its own column
    end

    # Back substitution
    @inbounds for i in N:-1:1
        s = zero(eltype(X))
        for j in (i + 1):N
            s += LU[i, j] * work[j, c]
        end
        work[i, c] = (work[i, c] - s) / LU[i, i]
    end

    # Write result: thread `tid` writes column `c` of X
    @inbounds for i in 1:N
        X[i, c, k] = work[i, c]
    end
end

function ka_batch_inv_lu_par!(X::CuArray{FT,3}, A::CuArray{FT,3}) where {FT}
    N = size(A, 1)
    batch = size(A, 3)
    backend = CUDABackend()
    batched_inv_lu_par_kernel!(backend, N)(X, A, Val(N); ndrange=(N * batch,))
    KernelAbstractions.synchronize(backend)
    return X
end

# ============================================================================
# cuBLAS reference (mirrors ext/gpu_batched_cuda.jl)
# ============================================================================
function cublas_batch_inv!(X::CuArray{FT,3}, A::CuArray{FT,3}) where {FT}
    pivot, info = CUDA.CUBLAS.getrf_strided_batched!(A, true)
    CUDA.synchronize()
    CUDA.CUBLAS.getri_strided_batched!(A, X, pivot)
    CUDA.synchronize()
    return X
end

# ============================================================================
# Benchmarking utilities
# ============================================================================
function make_invertible_batch(FT, N, batch)
    A_cpu = Array{FT}(undef, N, N, batch)
    for k in 1:batch
        A_cpu[:, :, k] = FT.(I(N)) + FT(0.3) * randn(FT, N, N)
    end
    return CuArray(A_cpu)
end

function time_kernel(f!, X, A_orig; nwarmup=3, nruns=20)
    for _ in 1:nwarmup
        A = copy(A_orig)
        f!(X, A)
    end
    CUDA.synchronize()

    times = Float64[]
    for _ in 1:nruns
        A = copy(A_orig)
        CUDA.synchronize()
        t = CUDA.@elapsed f!(X, A)
        push!(times, t)
    end
    return sort(times)[nruns ÷ 2]
end

function check_identity_residual(X, A_orig, N, batch)
    A_h = Array(A_orig)
    X_h = Array(X)
    max_err = zero(eltype(X_h))
    n_check = min(batch, 50)
    for k in 1:n_check
        res = A_h[:, :, k] * X_h[:, :, k] - I(N)
        max_err = max(max_err, maximum(abs.(res)))
    end
    return max_err
end

# Estimate shared memory in bytes for each kernel
smem_gj(N, FT) = N * 2N * sizeof(FT)                         # augmented [A|I]
smem_lu(N, FT) = N * N * sizeof(FT) + N * sizeof(Int32) + N * sizeof(FT)  # LU + piv + col_buf
smem_lu_par(N, FT) = 2 * N * N * sizeof(FT) + N * sizeof(Int32)           # LU + work + piv

const SMEM_LIMIT = 48 * 1024  # 48 KB shared memory limit (sm_80+)

# ============================================================================
# Main benchmark
# ============================================================================
function run_benchmark()
    matrix_sizes = [4, 8, 12, 16, 24, 32, 48, 64]
    batch_sizes  = [500, 2000, 5000]
    float_types  = [Float32, Float64]

    println("=" ^ 130)
    println("Batched Matrix Inversion Benchmark: cuBLAS vs Gauss-Jordan vs LU vs LU-parallel")
    println("GPU: ", CUDA.name(CUDA.device()))
    println("=" ^ 130)

    for FT in float_types
        println("\n", "-" ^ 130)
        @printf("%-7s | %4s | %5s | %10s | %10s | %10s | %10s | %7s | %7s | %7s | %9s | %9s | %9s\n",
                "Type", "N", "Batch", "cuBLAS ms", "GJ ms", "LU ms", "LU-par ms",
                "GJ spd", "LU spd", "LUP spd", "GJ resid", "LU resid", "LUP resid")
        println("-" ^ 130)

        for N in matrix_sizes
            for batch in batch_sizes
                A_orig = make_invertible_batch(FT, N, batch)
                X_cublas  = similar(A_orig)
                X_gj      = similar(A_orig)
                X_lu      = similar(A_orig)
                X_lu_par  = similar(A_orig)

                # cuBLAS reference
                t_cublas = time_kernel(cublas_batch_inv!, X_cublas, A_orig)
                A_tmp = copy(A_orig)
                cublas_batch_inv!(X_cublas, A_tmp)

                # Gauss-Jordan
                can_gj = smem_gj(N, FT) <= SMEM_LIMIT
                if can_gj
                    t_gj = time_kernel(ka_batch_inv_gj!, X_gj, A_orig)
                    A_tmp = copy(A_orig)
                    ka_batch_inv_gj!(X_gj, A_tmp)
                    res_gj = check_identity_residual(X_gj, A_orig, N, batch)
                else
                    t_gj = NaN; res_gj = NaN
                end

                # LU sequential solve
                can_lu = smem_lu(N, FT) <= SMEM_LIMIT
                if can_lu
                    t_lu = time_kernel(ka_batch_inv_lu!, X_lu, A_orig)
                    A_tmp = copy(A_orig)
                    ka_batch_inv_lu!(X_lu, A_tmp)
                    res_lu = check_identity_residual(X_lu, A_orig, N, batch)
                else
                    t_lu = NaN; res_lu = NaN
                end

                # LU parallel solve
                can_lu_par = smem_lu_par(N, FT) <= SMEM_LIMIT
                if can_lu_par
                    t_lu_par = time_kernel(ka_batch_inv_lu_par!, X_lu_par, A_orig)
                    A_tmp = copy(A_orig)
                    ka_batch_inv_lu_par!(X_lu_par, A_tmp)
                    res_lu_par = check_identity_residual(X_lu_par, A_orig, N, batch)
                else
                    t_lu_par = NaN; res_lu_par = NaN
                end

                spd_gj     = t_cublas / t_gj
                spd_lu     = t_cublas / t_lu
                spd_lu_par = t_cublas / t_lu_par

                @printf("%-7s | %4d | %5d | %10.4f | %10.4f | %10.4f | %10.4f | %6.2fx | %6.2fx | %6.2fx | %9.2e | %9.2e | %9.2e\n",
                        string(FT), N, batch,
                        t_cublas * 1000, t_gj * 1000, t_lu * 1000, t_lu_par * 1000,
                        spd_gj, spd_lu, spd_lu_par,
                        res_gj, res_lu, res_lu_par)

                CUDA.unsafe_free!(A_orig)
                CUDA.unsafe_free!(X_cublas)
                CUDA.unsafe_free!(X_gj)
                CUDA.unsafe_free!(X_lu)
                CUDA.unsafe_free!(X_lu_par)
            end
        end
    end

    println("\n", "=" ^ 130)
    println("Speedup > 1.0 means custom kernel is FASTER than cuBLAS")
    println("'resid' = max |A·A⁻¹ - I| over sampled batch elements (lower is better)")
    println("NaN = kernel skipped (shared memory exceeds $(SMEM_LIMIT ÷ 1024) KB limit)")
    println("=" ^ 130)
end

run_benchmark()
