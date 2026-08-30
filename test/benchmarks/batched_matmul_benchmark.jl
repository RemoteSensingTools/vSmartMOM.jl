#!/usr/bin/env julia
#
# Benchmark: Custom KernelAbstractions.jl batched matmul vs cuBLAS gemm_strided_batched
#
# Usage:  julia --project=test test/benchmarks/batched_matmul_benchmark.jl
#
# Production-shape Float32 sweep:
#   VSMARTMOM_BMM_MATRIX_SIZES=16,20,36,48,64 \
#   VSMARTMOM_BMM_BATCH_SIZES=4001,6001 \
#   VSMARTMOM_BMM_FLOAT_TYPES=Float32 \
#   VSMARTMOM_BMM_MATH_MODE=default \
#     julia --project=test test/benchmarks/batched_matmul_benchmark.jl
#
# `VSMARTMOM_BMM_MATH_MODE` accepts default, pedantic, or tf32.
#

using CUDA
using KernelAbstractions
using Printf
using StaticArrays

@assert CUDA.functional() "CUDA must be available to run this benchmark"

function configure_math_mode!()
    mode = lowercase(get(ENV, "VSMARTMOM_BMM_MATH_MODE", "default"))
    if mode == "default"
        CUDA.math_mode!(CUDA.DEFAULT_MATH)
    elseif mode == "pedantic"
        CUDA.math_mode!(CUDA.PEDANTIC_MATH)
    elseif mode == "tf32"
        CUDA.math_mode!(CUDA.FAST_MATH; precision=:TensorFloat32)
    else
        throw(ArgumentError(
            "VSMARTMOM_BMM_MATH_MODE must be default, pedantic, or tf32 (got $mode)"
        ))
    end
    return mode
end

function env_int_list(name, default)
    value = get(ENV, name, "")
    isempty(value) && return default
    return parse.(Int, split(value, ','))
end

function env_float_types(default)
    value = get(ENV, "VSMARTMOM_BMM_FLOAT_TYPES", "")
    isempty(value) && return default
    names = split(value, ',')
    all(name -> name in ("Float32", "Float64"), names) ||
        throw(ArgumentError("VSMARTMOM_BMM_FLOAT_TYPES accepts Float32 and Float64"))
    return map(name -> name == "Float32" ? Float32 : Float64, names)
end

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

# ==============================================================================
# Shared-memory variant with a fixed-size workgroup
# ==============================================================================
# Unlike `batched_matmul_smem_kernel!`, this kernel does not require one thread
# per output element.  A fixed-size workgroup cooperatively loads both matrices,
# then each thread computes one or more complete output dot products.  This keeps
# the simple, auditable GEMM equation while allowing N^2 > 1024 (notably N=36,
# 48, and 64 in the polarized RT configurations).
@kernel function batched_matmul_smem_tiled_kernel!(
    C, @Const(A), @Const(B), ::Val{N}, ::Val{WG}
) where {N,WG}
    k = @index(Group, Linear)
    tid = @index(Local, Linear)

    A_tile = @localmem eltype(C) (N, N)
    B_tile = @localmem eltype(C) (N, N)

    # Coalesced, strip-mined load: consecutive workitems read consecutive
    # column-major entries from each spectral matrix.
    idx = tid
    while idx <= N * N
        i = ((idx - 1) % N) + 1
        j = ((idx - 1) ÷ N) + 1
        @inbounds A_tile[i, j] = A[i, j, k]
        @inbounds B_tile[i, j] = B[i, j, k]
        idx += WG
    end
    @synchronize()

    # Strip-mine C over the same workgroup.  Each dot product stays in one
    # register accumulator, so no inter-thread reduction or second barrier is
    # required.
    out = tid
    while out <= N * N
        i = ((out - 1) % N) + 1
        j = ((out - 1) ÷ N) + 1
        acc = zero(eltype(C))
        @inbounds for l in 1:N
            acc += A_tile[i, l] * B_tile[l, j]
        end
        @inbounds C[i, j, k] = acc
        out += WG
    end
end

"""
    ka_batched_mul_smem_tiled!(C, A, B)

Shared-memory batched matrix multiply using at most 256 workitems per spectral
matrix.  Requires `2N²*sizeof(FT)` bytes of shared memory, but does not require
`N² <= 1024` threads.
"""
function ka_batched_mul_smem_tiled!(
    C::CuArray{FT,3}, A::CuArray{FT,3}, B::CuArray{FT,3}
) where {FT}
    N = size(A, 1)
    batch = size(A, 3)
    backend = CUDABackend()
    wg = min(256, max(32, cld(N * N, 32) * 32))
    batched_matmul_smem_tiled_kernel!(backend, wg)(
        C, A, B, Val(N), Val(wg); ndrange=wg * batch
    )
    KernelAbstractions.synchronize(backend)
    return C
end

# ==============================================================================
# Shared-memory + register column micro-tile variant
# ==============================================================================
# Each workitem owns one matrix row and NC adjacent output columns.  The A value
# is therefore loaded once per inner-loop iteration and reused across NC FMAs.
# The accumulator has compile-time length, allowing the GPU compiler to keep it
# in registers rather than shared memory.
@kernel function batched_matmul_smem_coltile_kernel!(
    C, @Const(A), @Const(B), ::Val{N}, ::Val{WG}, ::Val{NC}
) where {N,WG,NC}
    k = @index(Group, Linear)
    tid = @index(Local, Linear)

    A_tile = @localmem eltype(C) (N, N)
    B_tile = @localmem eltype(C) (N, N)

    idx = tid
    while idx <= N * N
        i = ((idx - 1) % N) + 1
        j = ((idx - 1) ÷ N) + 1
        @inbounds A_tile[i, j] = A[i, j, k]
        @inbounds B_tile[i, j] = B[i, j, k]
        idx += WG
    end
    @synchronize()

    ncol_tiles = cld(N, NC)
    nworkers = N * ncol_tiles
    worker = tid
    while worker <= nworkers
        i = ((worker - 1) % N) + 1
        j0 = ((worker - 1) ÷ N) * NC + 1
        acc = MVector{NC,eltype(C)}(undef)
        @inbounds for q in 1:NC
            acc[q] = zero(eltype(C))
        end

        @inbounds for l in 1:N
            a = A_tile[i, l]
            for q in 1:NC
                j = j0 + q - 1
                j <= N && (acc[q] += a * B_tile[l, j])
            end
        end

        @inbounds for q in 1:NC
            j = j0 + q - 1
            j <= N && (C[i, j, k] = acc[q])
        end
        worker += WG
    end
end

"""
    ka_batched_mul_smem_coltile!(C, A, B, Val(NC), Val(MAXWG))

Shared-memory GEMM with `NC` adjacent columns accumulated per workitem.  The
benchmark sweeps `NC` and maximum-workgroup values because register pressure
and occupancy cross over differently on A100 and L40S.
"""
function ka_batched_mul_smem_coltile!(
    C::CuArray{FT,3}, A::CuArray{FT,3}, B::CuArray{FT,3},
    ::Val{NC}, ::Val{MAXWG}
) where {FT,NC,MAXWG}
    N = size(A, 1)
    batch = size(A, 3)
    backend = CUDABackend()
    nworkers = N * cld(N, NC)
    wg = min(MAXWG, max(32, cld(nworkers, 32) * 32))
    batched_matmul_smem_coltile_kernel!(backend, wg)(
        C, A, B, Val(N), Val(wg), Val(NC); ndrange=wg * batch
    )
    KernelAbstractions.synchronize(backend)
    return C
end

# ==============================================================================
# Shared-memory output-column block variant
# ==============================================================================
# One spectral matrix is split across several workgroups.  Each workgroup loads
# all of A but only NBC columns of B, reducing shared memory and increasing the
# number of resident/independent blocks at N=48-64.  The extra A traffic is a
# deliberate trade for occupancy; the benchmark sweeps NBC to measure it.
@kernel function batched_matmul_smem_colblock_kernel!(
    C, @Const(A), @Const(B), ::Val{N}, ::Val{WG}, ::Val{NBC}
) where {N,WG,NBC}
    group = @index(Group, Linear)
    tid = @index(Local, Linear)
    ncol_blocks = cld(N, NBC)
    k = ((group - 1) ÷ ncol_blocks) + 1
    col_block = ((group - 1) % ncol_blocks) + 1
    j0 = (col_block - 1) * NBC + 1

    A_tile = @localmem eltype(C) (N, N)
    B_tile = @localmem eltype(C) (N, NBC)

    idx = tid
    while idx <= N * N
        i = ((idx - 1) % N) + 1
        j = ((idx - 1) ÷ N) + 1
        @inbounds A_tile[i, j] = A[i, j, k]
        idx += WG
    end

    idx = tid
    while idx <= N * NBC
        l = ((idx - 1) % N) + 1
        jlocal = ((idx - 1) ÷ N) + 1
        j = j0 + jlocal - 1
        @inbounds B_tile[l, jlocal] = j <= N ? B[l, j, k] : zero(eltype(C))
        idx += WG
    end
    @synchronize()

    out = tid
    while out <= N * NBC
        i = ((out - 1) % N) + 1
        jlocal = ((out - 1) ÷ N) + 1
        j = j0 + jlocal - 1
        if j <= N
            acc = zero(eltype(C))
            @inbounds for l in 1:N
                acc += A_tile[i, l] * B_tile[l, jlocal]
            end
            @inbounds C[i, j, k] = acc
        end
        out += WG
    end
end

function ka_batched_mul_smem_colblock!(
    C::CuArray{FT,3}, A::CuArray{FT,3}, B::CuArray{FT,3}, ::Val{NBC}
) where {FT,NBC}
    N = size(A, 1)
    batch = size(A, 3)
    backend = CUDABackend()
    wg = min(256, max(32, cld(N * NBC, 32) * 32))
    ncol_blocks = cld(N, NBC)
    batched_matmul_smem_colblock_kernel!(backend, wg)(
        C, A, B, Val(N), Val(wg), Val(NBC);
        ndrange=wg * batch * ncol_blocks
    )
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
function time_kernel(f!, args...; nwarmup=3, nruns=20)
    # Warmup
    for _ in 1:nwarmup
        f!(args...)
    end
    CUDA.synchronize()

    # Timed runs
    times = Float64[]
    for _ in 1:nruns
        CUDA.synchronize()
        t = CUDA.@elapsed f!(args...)
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
    matrix_sizes = env_int_list(
        "VSMARTMOM_BMM_MATRIX_SIZES", [4, 8, 12, 16, 24, 32, 48, 64]
    )
    batch_sizes = env_int_list("VSMARTMOM_BMM_BATCH_SIZES", [500, 2000, 5000])
    float_types = env_float_types([Float32, Float64])
    math_mode = configure_math_mode!()
    CUDA.seed!(0x5eed)

    println("=" ^ 211)
    println("Batched Matrix Multiplication Benchmark: Custom KA kernels vs cuBLAS")
    println("GPU: ", CUDA.name(CUDA.device()))
    println("cuBLAS math mode: ", math_mode)
    println("CUDA RNG seed: 0x5eed")
    println("=" ^ 211)

    for FT in float_types
        println("\n", "-" ^ 211)
        @printf("%-10s | %-6s | %-6s | %12s | %12s | %12s | %15s | %13s | %13s | %6s | %6s | %6s | %8s | %8s | %10s | %9s | %9s | %10s | %10s | %10s | %9s | %9s\n",
                "Type", "N", "Batch", "cuBLAS (ms)", "KA (ms)", "KA-smem (ms)", "KA-tiled (ms)", "KA-col (ms)",
                "KA-block (ms)", "cols", "wg", "block", "KA spd", "smem spd", "tiled spd", "col spd", "block spd",
                "KA err", "smem err", "tiled err", "col err", "block err")
        println("-" ^ 211)

        for N in matrix_sizes
            for batch in batch_sizes
                # Allocate random matrices
                A = CUDA.rand(FT, N, N, batch)
                B = CUDA.rand(FT, N, N, batch)
                C_cublas = similar(A)
                C_ka     = similar(A)
                C_smem   = similar(A)
                C_tiled  = similar(A)
                C_col    = similar(A)
                C_block  = similar(A)

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

                # Fixed-workgroup shared-memory kernel.  Stay within the
                # portable (non-opt-in) per-block shared-memory limit.
                smem_bytes = 2 * N * N * sizeof(FT)
                smem_limit = CUDA.attribute(
                    CUDA.device(), CUDA.DEVICE_ATTRIBUTE_MAX_SHARED_MEMORY_PER_BLOCK
                )
                if smem_bytes <= smem_limit
                    t_tiled = time_kernel(ka_batched_mul_smem_tiled!, C_tiled, A, B)
                    err_tiled = check_correctness(C_cublas, C_tiled)
                else
                    t_tiled = NaN
                    err_tiled = NaN
                end


                # Register micro-tile sweep for the sizes where N² workitems no
                # longer fit in one block.  Report the fastest width, but retain
                # every candidate's compilation and correctness check here.
                t_col = Inf
                err_col = NaN
                best_cols = 0
                best_col_wg = 0
                if N * N > 1024 && smem_bytes <= smem_limit
                    for NC in (2, 4, 8, 16)
                        for MAXWG in (128, 256, 512)
                            t_candidate = time_kernel(
                                ka_batched_mul_smem_coltile!, C_col, A, B,
                                Val(NC), Val(MAXWG)
                            )
                            err_candidate = check_correctness(C_cublas, C_col)
                            if t_candidate < t_col
                                t_col = t_candidate
                                err_col = err_candidate
                                best_cols = NC
                                best_col_wg = MAXWG
                            end
                        end
                    end
                end
                isfinite(t_col) || (t_col = NaN)

                # Split the output columns over multiple workgroups.  This is
                # primarily an N=48-64 occupancy experiment.
                t_block = Inf
                err_block = NaN
                best_block_cols = 0
                if N * N > 1024
                    for NBC in (8, 16, 32)
                        block_smem_bytes = (N * N + N * NBC) * sizeof(FT)
                        block_smem_bytes <= smem_limit || continue
                        t_candidate = time_kernel(
                            ka_batched_mul_smem_colblock!, C_block, A, B, Val(NBC)
                        )
                        err_candidate = check_correctness(C_cublas, C_block)
                        if t_candidate < t_block
                            t_block = t_candidate
                            err_block = err_candidate
                            best_block_cols = NBC
                        end
                    end
                end
                isfinite(t_block) || (t_block = NaN)

                spd_ka   = t_cublas / t_ka
                spd_smem = N * N <= 1024 ? t_cublas / t_smem : NaN
                spd_tiled = smem_bytes <= smem_limit ? t_cublas / t_tiled : NaN
                spd_col = best_cols > 0 ? t_cublas / t_col : NaN
                spd_block = best_block_cols > 0 ? t_cublas / t_block : NaN

                @printf("%-10s | %6d | %6d | %12.4f | %12.4f | %12.4f | %15.4f | %13.4f | %13.4f | %6d | %6d | %6d | %7.2fx | %7.2fx | %9.2fx | %8.2fx | %8.2fx | %10.2e | %10.2e | %10.2e | %9.2e | %9.2e\n",
                        string(FT), N, batch,
                        t_cublas * 1000, t_ka * 1000, t_smem * 1000, t_tiled * 1000,
                        t_col * 1000, t_block * 1000, best_cols, best_col_wg, best_block_cols,
                        spd_ka, spd_smem, spd_tiled, spd_col, spd_block,
                        err_ka, err_smem, err_tiled, err_col, err_block)

                # Free GPU memory
                CUDA.unsafe_free!(A)
                CUDA.unsafe_free!(B)
                CUDA.unsafe_free!(C_cublas)
                CUDA.unsafe_free!(C_ka)
                CUDA.unsafe_free!(C_smem)
                CUDA.unsafe_free!(C_tiled)
                CUDA.unsafe_free!(C_col)
                CUDA.unsafe_free!(C_block)
            end
        end
    end

    println("\n", "=" ^ 211)
    println("Speedup > 1.0 means custom kernel is FASTER than cuBLAS")
    println("=" ^ 211)
end

run_benchmark()
