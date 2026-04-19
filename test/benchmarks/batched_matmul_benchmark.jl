#!/usr/bin/env julia
#
# Benchmark: Custom KernelAbstractions.jl batched matmul vs cuBLAS gemm_strided_batched
#
# Usage:  julia --project=test test/benchmarks/batched_matmul_benchmark.jl
#

using CUDA
using KernelAbstractions
using Printf

@assert CUDA.functional() "CUDA must be available to run this benchmark"

# ============================================================================
# Custom KA kernel: one thread per output element C[i,j,k]
# ============================================================================
@kernel function batched_matmul_kernel!(C, @Const(A), @Const(B), ::Val{N}) where {N}
    i, j, k = @index(Global, NTuple)

    # Accumulate dot product A[i,:,k] * B[:,j,k]
    acc = zero(eltype(C))
    @inbounds for l in 1:N
        acc += A[i, l, k] * B[l, j, k]
    end
    @inbounds C[i, j, k] = acc
end

"""
    ka_batched_mul!(C, A, B)

Custom KernelAbstractions batched matrix multiply: C[:,:,k] = A[:,:,k] * B[:,:,k].
"""
function ka_batched_mul!(C::CuArray{FT,3}, A::CuArray{FT,3}, B::CuArray{FT,3}) where {FT}
    N = size(A, 1)
    batch = size(A, 3)
    backend = CUDABackend()
    # Launch N×N×batch threads; workgroup size tuned for small matrices
    wg = min(N * N, 256)
    batched_matmul_kernel!(backend, wg)(C, A, B, Val(N); ndrange=(N, N, batch))
    KernelAbstractions.synchronize(backend)
    return C
end

# ============================================================================
# Shared-memory variant for small matrices: load A,B tiles into shared mem
# ============================================================================
@kernel function batched_matmul_smem_kernel!(C, @Const(A), @Const(B), ::Val{N}) where {N}
    # k = batch index (one workgroup per batch element)
    k = @index(Group, Linear)
    # (i, j) = local thread index within the N×N workgroup
    li = @index(Local, Linear)
    i = ((li - 1) % N) + 1
    j = ((li - 1) ÷ N) + 1

    # Allocate shared memory for one A and B slice
    A_tile = @localmem eltype(C) (N, N)
    B_tile = @localmem eltype(C) (N, N)

    # Load tiles
    @inbounds A_tile[i, j] = A[i, j, k]
    @inbounds B_tile[i, j] = B[i, j, k]
    @synchronize()

    # Compute C[i,j,k]
    acc = zero(eltype(C))
    @inbounds for l in 1:N
        acc += A_tile[i, l] * B_tile[l, j]
    end
    @inbounds C[i, j, k] = acc
end

"""
    ka_batched_mul_smem!(C, A, B)

Shared-memory KA batched matmul. One workgroup (N×N threads) per batch element.
Best for small N where N² ≤ 1024.
"""
function ka_batched_mul_smem!(C::CuArray{FT,3}, A::CuArray{FT,3}, B::CuArray{FT,3}) where {FT}
    N = size(A, 1)
    batch = size(A, 3)
    @assert N * N <= 1024 "Shared-memory variant requires N² ≤ 1024 (got N=$N)"
    backend = CUDABackend()
    # N² threads per workgroup, batch workgroups
    batched_matmul_smem_kernel!(backend, N * N)(C, A, B, Val(N); ndrange=(N * N * batch,))
    KernelAbstractions.synchronize(backend)
    return C
end

# ============================================================================
# cuBLAS reference
# ============================================================================
function cublas_batched_mul!(C::CuArray{FT,3}, A::CuArray{FT,3}, B::CuArray{FT,3}) where {FT}
    CUDA.CUBLAS.gemm_strided_batched!('N', 'N', one(FT), A, B, zero(FT), C)
    CUDA.synchronize()
    return C
end

# ============================================================================
# Benchmarking utilities
# ============================================================================
function time_kernel(f!, C, A, B; nwarmup=3, nruns=20)
    # Warmup
    for _ in 1:nwarmup
        f!(C, A, B)
    end
    CUDA.synchronize()

    # Timed runs
    times = Float64[]
    for _ in 1:nruns
        CUDA.synchronize()
        t = CUDA.@elapsed f!(C, A, B)
        push!(times, t)
    end
    return sort(times)[nruns ÷ 2]  # median
end

function check_correctness(C_ref, C_test)
    return maximum(abs.(Array(C_ref) .- Array(C_test)))
end

# ============================================================================
# Main benchmark
# ============================================================================
function run_benchmark()
    matrix_sizes = [4, 8, 12, 16, 24, 32, 48, 64]
    batch_sizes  = [500, 2000, 5000]
    float_types  = [Float32, Float64]

    println("=" ^ 105)
    println("Batched Matrix Multiplication Benchmark: Custom KA kernels vs cuBLAS")
    println("GPU: ", CUDA.name(CUDA.device()))
    println("=" ^ 105)

    for FT in float_types
        println("\n", "-" ^ 105)
        @printf("%-10s | %-6s | %-6s | %12s | %12s | %12s | %8s | %8s | %10s | %10s\n",
                "Type", "N", "Batch", "cuBLAS (ms)", "KA (ms)", "KA-smem (ms)",
                "KA spd", "smem spd", "KA err", "smem err")
        println("-" ^ 105)

        for N in matrix_sizes
            for batch in batch_sizes
                # Allocate random matrices
                A = CUDA.rand(FT, N, N, batch)
                B = CUDA.rand(FT, N, N, batch)
                C_cublas = similar(A)
                C_ka     = similar(A)
                C_smem   = similar(A)

                # cuBLAS reference
                t_cublas = time_kernel(cublas_batched_mul!, C_cublas, A, B)

                # KA global-memory kernel
                t_ka = time_kernel(ka_batched_mul!, C_ka, A, B)
                err_ka = check_correctness(C_cublas, C_ka)

                # KA shared-memory kernel (only if N² ≤ 1024)
                if N * N <= 1024
                    t_smem = time_kernel(ka_batched_mul_smem!, C_smem, A, B)
                    err_smem = check_correctness(C_cublas, C_smem)
                else
                    t_smem = NaN
                    err_smem = NaN
                end

                spd_ka   = t_cublas / t_ka
                spd_smem = N * N <= 1024 ? t_cublas / t_smem : NaN

                @printf("%-10s | %6d | %6d | %12.4f | %12.4f | %12.4f | %7.2fx | %7.2fx | %10.2e | %10.2e\n",
                        string(FT), N, batch,
                        t_cublas * 1000, t_ka * 1000, t_smem * 1000,
                        spd_ka, spd_smem, err_ka, err_smem)

                # Free GPU memory
                CUDA.unsafe_free!(A)
                CUDA.unsafe_free!(B)
                CUDA.unsafe_free!(C_cublas)
                CUDA.unsafe_free!(C_ka)
                CUDA.unsafe_free!(C_smem)
            end
        end
    end

    println("\n", "=" ^ 105)
    println("Speedup > 1.0 means custom kernel is FASTER than cuBLAS")
    println("=" ^ 105)
end

run_benchmark()
