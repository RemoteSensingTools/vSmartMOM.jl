#!/usr/bin/env julia
#
# Benchmark: Fused batched kernels for RT solver patterns
#
# Tests the two dominant inversion patterns in vSmartMOM's RT solver:
#   Pattern 1 (geometric progression): X = B * inv(I - A*A)       [doubling]
#   Pattern 2 (interaction):           X = B * inv(I - A₁*A₂)     [adding]
#
# Compares:
#   (a) cuBLAS 3-step:  matmul → sub → getrf+getri → matmul
#   (b) Custom KA 3-step: KA matmul → sub → KA LU-par inv → KA matmul
#   (c) Fused single kernel: everything in one launch
#
# Usage:  julia --project=test test/benchmarks/batched_fused_benchmark.jl
#

using CUDA
using KernelAbstractions
using Printf
using LinearAlgebra

@assert CUDA.functional() "CUDA must be available to run this benchmark"

# ============================================================================
# Fused kernel: X = B * inv(I - A₁ * A₂)
# One workgroup of N threads per batch element.
# Shared memory layout: M (N×N) + work (N×N) + piv (N) + B_tile (N×N)
#
# Steps (all in shared memory):
#   1. Compute M = I - A₁*A₂
#   2. LU factorize M with partial pivoting
#   3. Parallel triangular solve: each thread solves one column of M⁻¹
#   4. Compute X = B * M⁻¹ (each thread computes one row of the result)
# ============================================================================
@kernel function fused_inv_mul_kernel!(X, @Const(A1), @Const(A2), @Const(B), ::Val{N}) where {N}
    k = @index(Group, Linear)
    tid = @index(Local, Linear)  # tid ∈ 1:N

    M    = @localmem eltype(X) (N, N)
    piv  = @localmem Int32 (N,)
    work = @localmem eltype(X) (N, N)   # columns of M⁻¹
    Bt   = @localmem eltype(X) (N, N)   # B tile

    # Step 1: M[tid, :] = I[tid, :] - (A1[tid, :, k] * A2[:, :, k])[tid, :]
    # Also load B into shared memory
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

    # Step 2: LU factorization with partial pivoting
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

    # Step 3: Parallel triangular solve — each thread solves column `tid` of M⁻¹
    c = tid

    # Initialize: work[:,c] = P * e_c
    @inbounds for i in 1:N
        work[i, c] = (piv[i] == Int32(c)) ? one(eltype(X)) : zero(eltype(X))
    end

    # Forward substitution (L * y = Pb, each thread its own column)
    @inbounds for i in 2:N
        s = zero(eltype(X))
        for j in 1:(i - 1)
            s += M[i, j] * work[j, c]
        end
        work[i, c] -= s
    end

    # Back substitution (U * x = y)
    @inbounds for i in N:-1:1
        s = zero(eltype(X))
        for j in (i + 1):N
            s += M[i, j] * work[j, c]
        end
        work[i, c] = (work[i, c] - s) / M[i, i]
    end
    @synchronize()

    # Step 4: X = B * M⁻¹ — each thread computes row `tid` of the output
    @inbounds for j in 1:N
        s = zero(eltype(X))
        for l in 1:N
            s += Bt[tid, l] * work[l, j]
        end
        X[tid, j, k] = s
    end
end

function ka_fused_inv_mul!(X::CuArray{FT,3}, A1::CuArray{FT,3}, A2::CuArray{FT,3}, B::CuArray{FT,3}) where {FT}
    N = size(A1, 1)
    batch = size(A1, 3)
    backend = CUDABackend()
    fused_inv_mul_kernel!(backend, N)(X, A1, A2, B, Val(N); ndrange=(N * batch,))
    KernelAbstractions.synchronize(backend)
    return X
end

# ============================================================================
# Fused kernel for geometric progression: X = B * inv(I - A * A)
# Same as above but A₁ = A₂ = A, saves one global memory read
# ============================================================================
@kernel function fused_gp_kernel!(X, @Const(A), @Const(B), ::Val{N}) where {N}
    k = @index(Group, Linear)
    tid = @index(Local, Linear)

    M    = @localmem eltype(X) (N, N)
    piv  = @localmem Int32 (N,)
    work = @localmem eltype(X) (N, N)
    At   = @localmem eltype(X) (N, N)   # A tile (read once, used for A*A and later)
    Bt   = @localmem eltype(X) (N, N)

    # Load A and B into shared memory
    @inbounds for j in 1:N
        At[tid, j] = A[tid, j, k]
        Bt[tid, j] = B[tid, j, k]
    end
    @synchronize()

    # Step 1: M = I - A*A (using shared A tile)
    @inbounds for j in 1:N
        s = zero(eltype(X))
        for l in 1:N
            s += At[tid, l] * At[l, j]
        end
        M[tid, j] = ((tid == j) ? one(eltype(X)) : zero(eltype(X))) - s
    end
    @inbounds piv[tid] = Int32(tid)
    @synchronize()

    # Step 2: LU factorization with partial pivoting
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

    # Step 3: Parallel triangular solve
    c = tid
    @inbounds for i in 1:N
        work[i, c] = (piv[i] == Int32(c)) ? one(eltype(X)) : zero(eltype(X))
    end

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

    # Step 4: X = B * M⁻¹
    @inbounds for j in 1:N
        s = zero(eltype(X))
        for l in 1:N
            s += Bt[tid, l] * work[l, j]
        end
        X[tid, j, k] = s
    end
end

function ka_fused_gp!(X::CuArray{FT,3}, A::CuArray{FT,3}, B::CuArray{FT,3}) where {FT}
    N = size(A, 1)
    batch = size(A, 3)
    backend = CUDABackend()
    fused_gp_kernel!(backend, N)(X, A, B, Val(N); ndrange=(N * batch,))
    KernelAbstractions.synchronize(backend)
    return X
end

# ============================================================================
# Reference implementations: cuBLAS 3-step
# ============================================================================

"""cuBLAS 3-step: X = B * inv(I - A₁ * A₂)"""
function cublas_3step_inv_mul!(X::CuArray{FT,3}, A1::CuArray{FT,3}, A2::CuArray{FT,3}, B::CuArray{FT,3},
                                I_static::CuArray{FT,3}, temp::CuArray{FT,3}) where {FT}
    # Step 1: temp = A1 * A2
    CUDA.CUBLAS.gemm_strided_batched!('N', 'N', one(FT), A1, A2, zero(FT), temp)
    # Step 2: temp = I - temp
    temp .= I_static .- temp
    CUDA.synchronize()
    # Step 3: X = inv(temp)
    pivot, info = CUDA.CUBLAS.getrf_strided_batched!(temp, true)
    CUDA.synchronize()
    CUDA.CUBLAS.getri_strided_batched!(temp, X, pivot)
    CUDA.synchronize()
    # Step 4: X = B * inv(M)
    CUDA.CUBLAS.gemm_strided_batched!('N', 'N', one(FT), B, X, zero(FT), temp)
    CUDA.synchronize()
    X .= temp
    return X
end

"""cuBLAS 3-step for geometric progression: X = B * inv(I - A * A)"""
function cublas_3step_gp!(X::CuArray{FT,3}, A::CuArray{FT,3}, B::CuArray{FT,3},
                           I_static::CuArray{FT,3}, temp::CuArray{FT,3}) where {FT}
    cublas_3step_inv_mul!(X, A, A, B, I_static, temp)
end

# ============================================================================
# Benchmarking utilities
# ============================================================================
function make_rt_matrices(FT, N, batch)
    # Create matrices representative of RT solver:
    # r (reflectance) should have spectral radius < 1 for I-r*r to be invertible
    # t (transmittance) is close to identity
    r = CuArray(FT(0.3) * randn(FT, N, N, batch))  # small reflectance
    t = CuArray(FT.(I(N)) .+ FT(0.1) * randn(FT, N, N, batch))  # near-identity transmission
    I_static = CuArray(repeat(FT.(I(N)), 1, 1, batch))
    return r, t, I_static
end

function time_fused(f!, X, args...; nwarmup=3, nruns=20)
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

function time_cublas_3step(f!, X, args...; nwarmup=3, nruns=20)
    # cuBLAS getrf modifies input, so we need copies
    # args contains (A1, A2, B, I_static, temp) or (A, B, I_static, temp)
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

function check_residual(X, A1, A2, B, I_static, N, batch)
    # Check ||X - B * inv(I - A1*A2)||_max using CPU computation
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

# Shared memory estimate (bytes) for fused kernels: 4 × N² × sizeof(FT) + N × sizeof(Int32)
smem_fused(N, FT) = 4 * N * N * sizeof(FT) + N * sizeof(Int32)
const SMEM_LIMIT = 48 * 1024

# ============================================================================
# Main benchmark
# ============================================================================
function run_benchmark()
    matrix_sizes = [4, 8, 12, 16, 24, 32, 48]
    batch_sizes  = [500, 2000, 5000, 10000, 50000]
    float_types  = [Float32, Float64]

    println("=" ^ 120)
    println("Fused Batched Kernel Benchmark: B * inv(I - A₁*A₂)")
    println("Pattern: geometric progression (doubling) and interaction inv+mul")
    println("GPU: ", CUDA.name(CUDA.device()))
    println("=" ^ 120)

    for FT in float_types
        println("\n", "-" ^ 120)
        println("--- Geometric Progression: X = t⁺⁺ * inv(I - r⁻⁺ * r⁻⁺) ---")
        @printf("%-7s | %4s | %5s | %10s | %10s | %8s | %9s | %9s\n",
                "Type", "N", "Batch", "cuBLAS ms", "fused ms", "speedup", "fused res", "smem KB")
        println("-" ^ 120)

        for N in matrix_sizes
            can_fuse = smem_fused(N, FT) <= SMEM_LIMIT
            for batch in batch_sizes
                r, t, I_s = make_rt_matrices(FT, N, batch)
                X_cublas = similar(r)
                X_fused  = similar(r)
                temp     = similar(r)

                # cuBLAS 3-step (uses copies internally since getrf is in-place)
                # We time the full pattern: matmul + sub + inv + matmul
                t_cublas = time_cublas_3step(cublas_3step_gp!, X_cublas, r, t, I_s, temp)

                if can_fuse
                    t_fused = time_fused(ka_fused_gp!, X_fused, r, t)
                    res_fused = check_residual(X_fused, r, r, t, I_s, N, batch)
                    spd = t_cublas / t_fused
                else
                    t_fused = NaN; res_fused = NaN; spd = NaN
                end

                smem_kb = smem_fused(N, FT) / 1024

                @printf("%-7s | %4d | %5d | %10.4f | %10.4f | %7.2fx | %9.2e | %8.1f\n",
                        string(FT), N, batch,
                        t_cublas * 1000, t_fused * 1000, spd, res_fused, smem_kb)

                CUDA.unsafe_free!(r); CUDA.unsafe_free!(t); CUDA.unsafe_free!(I_s)
                CUDA.unsafe_free!(X_cublas); CUDA.unsafe_free!(X_fused); CUDA.unsafe_free!(temp)
            end
        end

        println()
        println("--- Interaction: X = T⁻⁻ * inv(I - r⁻⁺ * R⁺⁻) ---")
        @printf("%-7s | %4s | %5s | %10s | %10s | %8s | %9s | %9s\n",
                "Type", "N", "Batch", "cuBLAS ms", "fused ms", "speedup", "fused res", "smem KB")
        println("-" ^ 120)

        for N in matrix_sizes
            can_fuse = smem_fused(N, FT) <= SMEM_LIMIT
            for batch in batch_sizes
                # For interaction: A1 = r⁻⁺, A2 = R⁺⁻ (different matrices)
                r, t, I_s = make_rt_matrices(FT, N, batch)
                R = CuArray(FT(0.3) * randn(FT, N, N, batch))  # second reflectance (composite)
                T_comp = CuArray(FT.(I(N)) .+ FT(0.1) * randn(FT, N, N, batch))  # composite transmission

                X_cublas = similar(r)
                X_fused  = similar(r)
                temp     = similar(r)

                # cuBLAS 3-step: T_comp * inv(I - r * R)
                t_cublas = time_cublas_3step(cublas_3step_inv_mul!, X_cublas, r, R, T_comp, I_s, temp)

                if can_fuse
                    t_fused = time_fused(ka_fused_inv_mul!, X_fused, r, R, T_comp)
                    res_fused = check_residual(X_fused, r, R, T_comp, I_s, N, batch)
                    spd = t_cublas / t_fused
                else
                    t_fused = NaN; res_fused = NaN; spd = NaN
                end

                smem_kb = smem_fused(N, FT) / 1024

                @printf("%-7s | %4d | %5d | %10.4f | %10.4f | %7.2fx | %9.2e | %8.1f\n",
                        string(FT), N, batch,
                        t_cublas * 1000, t_fused * 1000, spd, res_fused, smem_kb)

                CUDA.unsafe_free!(r); CUDA.unsafe_free!(t); CUDA.unsafe_free!(I_s)
                CUDA.unsafe_free!(R); CUDA.unsafe_free!(T_comp)
                CUDA.unsafe_free!(X_cublas); CUDA.unsafe_free!(X_fused); CUDA.unsafe_free!(temp)
            end
        end
    end

    println("\n", "=" ^ 120)
    println("Speedup > 1.0 means fused kernel is FASTER than cuBLAS 3-step")
    println("NaN = kernel skipped (shared memory exceeds $(SMEM_LIMIT ÷ 1024) KB)")
    println("=" ^ 120)
end

run_benchmark()
