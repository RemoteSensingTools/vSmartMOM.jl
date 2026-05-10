#!/usr/bin/env julia
#
# Benchmark: Fused direct-solve kernels (v2) vs fused inv+mul (v1) vs cuBLAS
#
# v1: LU factorize M, solve for all columns of M⁻¹, then multiply B * M⁻¹
# v2: LU factorize M, then right-solve X * M = B directly (never form M⁻¹)
#
# v2 saves one N×N shared memory tile and one O(N³) stage.
# Shared memory: v1 = 4×N²×sizeof(FT), v2 = 2×N²×sizeof(FT) + N×sizeof(Int32)
#
# Usage:  julia --project=test test/benchmarks/batched_fused_v2_benchmark.jl
#

using CUDA
using KernelAbstractions
using Printf
using LinearAlgebra

@assert CUDA.functional() "CUDA must be available to run this benchmark"

# ============================================================================
# v2: Fused direct-solve kernel for X = B * inv(I - A₁ * A₂)
#
# After LU factorizing M = PLU, we want X = B * M⁻¹ = B * U⁻¹ * L⁻¹ * P.
# Each thread handles one ROW of X independently:
#   Step a: Forward solve  Z[tid,:] * U = B[tid,:]  (process j=1..N)
#   Step b: Backward solve W[tid,:] * L = Z[tid,:]  (process j=N..1, unit diag)
#   Step c: Permute columns: X[tid, piv[j]] = W[tid, j]
#
# Shared memory: M (N×N) for LU + Bt (N×N) for B/solve workspace + piv (N)
# Total: 2 × N² × sizeof(FT) + N × sizeof(Int32)
# ============================================================================
@kernel function fused_direct_solve_kernel!(X, @Const(A1), @Const(A2), @Const(B), ::Val{N}) where {N}
    k = @index(Group, Linear)
    tid = @index(Local, Linear)

    M   = @localmem eltype(X) (N, N)
    piv = @localmem Int32 (N,)
    Bt  = @localmem eltype(X) (N, N)   # B tile, then in-place solve workspace

    # Step 1: Compute M = I - A₁*A₂, load B
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

    # Step 2: LU factorization with partial pivoting (all threads cooperate)
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

    # Step 3: Right-side direct solve — each thread handles its own row independently
    # X = B * M⁻¹ = B * U⁻¹ * L⁻¹ * P
    # We solve in-place on Bt[tid, :]

    # (a) Forward solve: Z * U = B  →  Z[tid,j] = (B[tid,j] - Σ_{l=1}^{j-1} Z[tid,l]*U[l,j]) / U[j,j]
    @inbounds for j in 1:N
        s = zero(eltype(X))
        for l in 1:(j - 1)
            s += Bt[tid, l] * M[l, j]   # M[l,j] = U[l,j] for l ≤ j
        end
        Bt[tid, j] = (Bt[tid, j] - s) / M[j, j]
    end
    # No sync needed — each thread only reads/writes its own row of Bt,
    # and M (the LU factors) is read-only after factorization

    # (b) Backward solve: W * L = Z  →  W[tid,j] = Z[tid,j] - Σ_{l=j+1}^{N} W[tid,l]*L[l,j]
    # L has unit diagonal, so no division needed
    @inbounds for j in (N - 1):-1:1
        s = zero(eltype(X))
        for l in (j + 1):N
            s += Bt[tid, l] * M[l, j]   # M[l,j] = L[l,j] for l > j
        end
        Bt[tid, j] -= s
    end

    # (c) Permute columns: X[tid, piv[j]] = W[tid, j]
    # piv maps original row → pivoted row, so we write to piv[j]-th column
    @inbounds for j in 1:N
        X[tid, piv[j], k] = Bt[tid, j]
    end
end

function ka_fused_direct_solve!(X::CuArray{FT,3}, A1::CuArray{FT,3}, A2::CuArray{FT,3}, B::CuArray{FT,3}) where {FT}
    N = size(A1, 1)
    batch = size(A1, 3)
    backend = CUDABackend()
    fused_direct_solve_kernel!(backend, N)(X, A1, A2, B, Val(N); ndrange=(N * batch,))
    KernelAbstractions.synchronize(backend)
    return X
end

# ============================================================================
# v2 GP variant: X = B * inv(I - A * A)  (A₁ = A₂ = A, load A once)
# ============================================================================
@kernel function fused_direct_solve_gp_kernel!(X, @Const(A), @Const(B), ::Val{N}) where {N}
    k = @index(Group, Linear)
    tid = @index(Local, Linear)

    At  = @localmem eltype(X) (N, N)   # A tile (reused as M after matmul)
    M   = @localmem eltype(X) (N, N)
    piv = @localmem Int32 (N,)
    Bt  = @localmem eltype(X) (N, N)

    # Load A and B
    @inbounds for j in 1:N
        At[tid, j] = A[tid, j, k]
        Bt[tid, j] = B[tid, j, k]
    end
    @synchronize()

    # M = I - A*A (using shared A tile)
    @inbounds for j in 1:N
        s = zero(eltype(X))
        for l in 1:N
            s += At[tid, l] * At[l, j]
        end
        M[tid, j] = ((tid == j) ? one(eltype(X)) : zero(eltype(X))) - s
    end
    @inbounds piv[tid] = Int32(tid)
    @synchronize()

    # LU factorization
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

    # Direct solve: X = B * U⁻¹ * L⁻¹ * P (in-place on Bt)
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

function ka_fused_direct_solve_gp!(X::CuArray{FT,3}, A::CuArray{FT,3}, B::CuArray{FT,3}) where {FT}
    N = size(A, 1)
    batch = size(A, 3)
    backend = CUDABackend()
    fused_direct_solve_gp_kernel!(backend, N)(X, A, B, Val(N); ndrange=(N * batch,))
    KernelAbstractions.synchronize(backend)
    return X
end

# ============================================================================
# v1: Fused inv+mul kernel (from batched_fused_benchmark.jl) for comparison
# ============================================================================
@kernel function fused_inv_mul_v1_kernel!(X, @Const(A1), @Const(A2), @Const(B), ::Val{N}) where {N}
    k = @index(Group, Linear)
    tid = @index(Local, Linear)

    M    = @localmem eltype(X) (N, N)
    piv  = @localmem Int32 (N,)
    work = @localmem eltype(X) (N, N)
    Bt   = @localmem eltype(X) (N, N)

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

    # LU factorization
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

    # Parallel column solve for M⁻¹
    c = tid
    @inbounds for i in 1:N
        work[i, c] = (piv[i] == Int32(c)) ? one(eltype(X)) : zero(eltype(X))
    end
    @synchronize()

    @inbounds for i in 2:N
        s = zero(eltype(X))
        for j in 1:(i - 1)
            s += M[i, j] * work[j, c]
        end
        work[i, c] -= s
    end

    @inbounds for i in N:-1:1
        s = zero(eltype(X))
        for j in (i + 1):N
            s += M[i, j] * work[j, c]
        end
        work[i, c] = (work[i, c] - s) / M[i, i]
    end
    @synchronize()

    # X = B * M⁻¹
    @inbounds for j in 1:N
        s = zero(eltype(X))
        for l in 1:N
            s += Bt[tid, l] * work[l, j]
        end
        X[tid, j, k] = s
    end
end

function ka_fused_v1!(X::CuArray{FT,3}, A1::CuArray{FT,3}, A2::CuArray{FT,3}, B::CuArray{FT,3}) where {FT}
    N = size(A1, 1)
    batch = size(A1, 3)
    backend = CUDABackend()
    fused_inv_mul_v1_kernel!(backend, N)(X, A1, A2, B, Val(N); ndrange=(N * batch,))
    KernelAbstractions.synchronize(backend)
    return X
end

# ============================================================================
# cuBLAS reference
# ============================================================================
function cublas_3step!(X::CuArray{FT,3}, A1::CuArray{FT,3}, A2::CuArray{FT,3}, B::CuArray{FT,3},
                       I_static::CuArray{FT,3}, temp::CuArray{FT,3}) where {FT}
    CUDA.CUBLAS.gemm_strided_batched!('N', 'N', one(FT), A1, A2, zero(FT), temp)
    temp .= I_static .- temp
    CUDA.synchronize()
    pivot, info = CUDA.CUBLAS.getrf_strided_batched!(temp, true)
    CUDA.synchronize()
    CUDA.CUBLAS.getri_strided_batched!(temp, X, pivot)
    CUDA.synchronize()
    CUDA.CUBLAS.gemm_strided_batched!('N', 'N', one(FT), B, X, zero(FT), temp)
    CUDA.synchronize()
    X .= temp
    return X
end

# ============================================================================
# Benchmarking
# ============================================================================
function make_rt_matrices(FT, N, batch)
    r = CuArray(FT(0.3) * randn(FT, N, N, batch))
    t = CuArray(FT.(I(N)) .+ FT(0.1) * randn(FT, N, N, batch))
    I_static = CuArray(repeat(FT.(I(N)), 1, 1, batch))
    return r, t, I_static
end

function time_kernel(f!, X, args...; nwarmup=3, nruns=20)
    for _ in 1:nwarmup
        f!(X, args...)
    end
    CUDA.synchronize()
    times = Float64[]
    for _ in 1:nruns
        CUDA.synchronize()
        t = CUDA.@elapsed f!(X, args...)
        push!(times, t)
    end
    return sort(times)[nruns ÷ 2]
end

function check_residual(X, A1, A2, B, N, batch)
    A1_h = Array(A1); A2_h = Array(A2); B_h = Array(B); X_h = Array(X)
    max_err = zero(eltype(X_h))
    n_check = min(batch, 30)
    for k in 1:n_check
        M = I(N) - A1_h[:,:,k] * A2_h[:,:,k]
        X_ref = B_h[:,:,k] * inv(M)
        max_err = max(max_err, maximum(abs.(X_h[:,:,k] .- X_ref)))
    end
    return max_err
end

# Shared memory estimates
smem_v1(N, FT) = 4 * N * N * sizeof(FT) + N * sizeof(Int32)      # M + work + Bt + ...
smem_v2(N, FT) = 2 * N * N * sizeof(FT) + N * sizeof(Int32)      # M + Bt + piv
smem_v2_gp(N, FT) = 3 * N * N * sizeof(FT) + N * sizeof(Int32)   # At + M + Bt + piv
const SMEM_LIMIT = 48 * 1024

# ============================================================================
# Main benchmark
# ============================================================================
function run_benchmark()
    matrix_sizes = [4, 8, 12, 16, 24, 32, 48, 64]
    batch_sizes  = [500, 2000, 5000, 10000, 50000]
    float_types  = [Float32, Float64]

    println("=" ^ 130)
    println("Fused Direct-Solve (v2) vs Fused Inv+Mul (v1) vs cuBLAS")
    println("GPU: ", CUDA.name(CUDA.device()))
    println("=" ^ 130)

    for FT in float_types
        # --- Interaction pattern: X = B * inv(I - A₁ * A₂) ---
        println("\n", "-" ^ 130)
        println("--- Interaction: X = T⁻⁻ * inv(I - r⁻⁺ * R⁺⁻)  [$FT] ---")
        @printf("%-4s | %5s | %10s | %10s | %10s | %8s | %8s | %9s | %9s | %6s | %6s\n",
                "N", "Batch", "cuBLAS ms", "v1 ms", "v2 ms", "v1 spd", "v2 spd",
                "v1 resid", "v2 resid", "v1 KB", "v2 KB")
        println("-" ^ 130)

        for N in matrix_sizes
            can_v1 = smem_v1(N, FT) <= SMEM_LIMIT
            can_v2 = smem_v2(N, FT) <= SMEM_LIMIT

            for batch in batch_sizes
                r, t, I_s = make_rt_matrices(FT, N, batch)
                R = CuArray(FT(0.3) * randn(FT, N, N, batch))
                T_comp = CuArray(FT.(I(N)) .+ FT(0.1) * randn(FT, N, N, batch))
                X_cublas = similar(r)
                X_v1     = similar(r)
                X_v2     = similar(r)
                temp     = similar(r)

                t_cublas = time_kernel(cublas_3step!, X_cublas, r, R, T_comp, I_s, temp)

                if can_v1
                    t_v1 = time_kernel(ka_fused_v1!, X_v1, r, R, T_comp)
                    res_v1 = check_residual(X_v1, r, R, T_comp, N, batch)
                else
                    t_v1 = NaN; res_v1 = NaN
                end

                if can_v2
                    t_v2 = time_kernel(ka_fused_direct_solve!, X_v2, r, R, T_comp)
                    res_v2 = check_residual(X_v2, r, R, T_comp, N, batch)
                else
                    t_v2 = NaN; res_v2 = NaN
                end

                @printf("%4d | %5d | %10.4f | %10.4f | %10.4f | %7.2fx | %7.2fx | %9.2e | %9.2e | %5.1f | %5.1f\n",
                        N, batch,
                        t_cublas * 1000, t_v1 * 1000, t_v2 * 1000,
                        t_cublas / t_v1, t_cublas / t_v2,
                        res_v1, res_v2,
                        smem_v1(N, FT) / 1024, smem_v2(N, FT) / 1024)

                CUDA.unsafe_free!(r); CUDA.unsafe_free!(t); CUDA.unsafe_free!(I_s)
                CUDA.unsafe_free!(R); CUDA.unsafe_free!(T_comp)
                CUDA.unsafe_free!(X_cublas); CUDA.unsafe_free!(X_v1); CUDA.unsafe_free!(X_v2)
                CUDA.unsafe_free!(temp)
            end
        end

        # --- GP pattern: X = B * inv(I - A * A) ---
        println()
        println("--- Geometric Progression: X = t⁺⁺ * inv(I - r⁻⁺ * r⁻⁺)  [$FT] ---")
        @printf("%-4s | %5s | %10s | %10s | %10s | %8s | %8s | %9s | %9s | %6s | %6s\n",
                "N", "Batch", "cuBLAS ms", "v1 ms", "v2 ms", "v1 spd", "v2 spd",
                "v1 resid", "v2 resid", "v1 KB", "v2 KB")
        println("-" ^ 130)

        for N in matrix_sizes
            can_v1 = smem_v1(N, FT) <= SMEM_LIMIT
            can_v2_gp = smem_v2_gp(N, FT) <= SMEM_LIMIT

            for batch in batch_sizes
                r, t, I_s = make_rt_matrices(FT, N, batch)
                X_cublas = similar(r)
                X_v1     = similar(r)
                X_v2     = similar(r)
                temp     = similar(r)

                t_cublas = time_kernel(cublas_3step!, X_cublas, r, r, t, I_s, temp)

                if can_v1
                    t_v1 = time_kernel(ka_fused_v1!, X_v1, r, r, t)
                    res_v1 = check_residual(X_v1, r, r, t, N, batch)
                else
                    t_v1 = NaN; res_v1 = NaN
                end

                if can_v2_gp
                    t_v2 = time_kernel(ka_fused_direct_solve_gp!, X_v2, r, t)
                    res_v2 = check_residual(X_v2, r, r, t, N, batch)
                else
                    t_v2 = NaN; res_v2 = NaN
                end

                @printf("%4d | %5d | %10.4f | %10.4f | %10.4f | %7.2fx | %7.2fx | %9.2e | %9.2e | %5.1f | %5.1f\n",
                        N, batch,
                        t_cublas * 1000, t_v1 * 1000, t_v2 * 1000,
                        t_cublas / t_v1, t_cublas / t_v2,
                        res_v1, res_v2,
                        smem_v1(N, FT) / 1024, smem_v2_gp(N, FT) / 1024)

                CUDA.unsafe_free!(r); CUDA.unsafe_free!(t); CUDA.unsafe_free!(I_s)
                CUDA.unsafe_free!(X_cublas); CUDA.unsafe_free!(X_v1); CUDA.unsafe_free!(X_v2)
                CUDA.unsafe_free!(temp)
            end
        end
    end

    println("\n", "=" ^ 130)
    println("v1 = fused inv+mul (4 tiles smem)   v2 = fused direct solve (2 tiles smem)")
    println("Speedup is vs cuBLAS multi-step. NaN = skipped (shared memory limit).")
    println("=" ^ 130)
end

run_benchmark()
