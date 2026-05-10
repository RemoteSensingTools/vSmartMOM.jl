#!/usr/bin/env julia
#=
Toy benchmark mimicking the Raman (inelastic) doubling/interaction memory
and compute patterns in vSmartMOM.jl.

Goal: Identify bottlenecks among:
  1. GPU memory allocation of 4D tensors
  2. Batched matrix multiply (⊠) throughput on 3D slices
  3. Batched matrix inversion throughput
  4. CPU↔GPU transfer cost for slice-at-a-time staging
  5. In-loop temporary allocation overhead

Run:  julia raman_memory_benchmark.jl [--small|--medium|--large|--production]
=#

using CUDA
using NNlib
using LinearAlgebra
using Printf
using Statistics

# ── Dimension presets ──────────────────────────────────────────────────────────

const PRESETS = Dict(
    "small"      => (NquadN=15, nSpec=68,   nRaman=20,  ndoubl=8, nLayers=12, nFourier=3),
    "medium"     => (NquadN=15, nSpec=500,  nRaman=80,  ndoubl=8, nLayers=12, nFourier=3),
    "large"      => (NquadN=15, nSpec=2000, nRaman=172, ndoubl=8, nLayers=12, nFourier=3),
    "production" => (NquadN=15, nSpec=5424, nRaman=172, ndoubl=8, nLayers=12, nFourier=3),
)

function parse_preset()
    for a in ARGS
        key = replace(a, "--" => "")
        haskey(PRESETS, key) && return PRESETS[key]
    end
    return PRESETS["medium"]  # default
end

# ── Helpers ────────────────────────────────────────────────────────────────────

# Batched multiply shorthand (same as vSmartMOM uses)
const ⊠ = NNlib.batched_mul

"""Batched inversion via LU on GPU — mimics batch_inv! in vSmartMOM"""
function batch_inv_gpu!(out::CuArray{T,3}, A::CuArray{T,3}) where T
    # Simple per-slice inversion using CUBLAS batched getrf/getri would be ideal.
    # For this benchmark, approximate with batched_mul of pre-inverted slices
    # since we're measuring memory/transfer, not algebra accuracy.
    # In practice vSmartMOM calls CUBLAS batched LU.
    n, _, nspec = size(A)
    # Just do A \ I via broadcast for timing purposes
    I_mat = CUDA.CuArray(Matrix{T}(I, n, n))
    for k in 1:nspec
        @views out[:,:,k] .= A[:,:,k]  # placeholder — real code does LU
    end
    CUDA.synchronize()
end

"""Simple batched inversion approximation for benchmarking"""
function batch_inv_simple!(out, A)
    out .= A  # placeholder for timing envelope
    CUDA.synchronize()
end

function gpu_memory_gb()
    free, total = CUDA.Mem.info()
    used = (total - free) / 1e9
    return (used=used, free=free/1e9, total=total/1e9)
end

"""Copy a 3D slice from CPU 4D array to GPU 3D buffer (contiguous copy)"""
function h2d_slice!(gpu_buf::CuArray{T,3}, cpu_4d::Array{T,4}, Δn::Int) where T
    n1, n2, n3, _ = size(cpu_4d)
    # Create a contiguous view that CUDA can bulk-copy
    offset = (Δn - 1) * n1 * n2 * n3
    src = unsafe_wrap(Array, pointer(cpu_4d, offset + 1), (n1, n2, n3))
    copyto!(gpu_buf, src)
end

"""Copy GPU 3D buffer back to a 3D slice in CPU 4D array"""
function d2h_slice!(cpu_4d::Array{T,4}, gpu_buf::CuArray{T,3}, Δn::Int) where T
    n1, n2, n3, _ = size(cpu_4d)
    offset = (Δn - 1) * n1 * n2 * n3
    dest = unsafe_wrap(Array, pointer(cpu_4d, offset + 1), (n1, n2, n3))
    copyto!(dest, gpu_buf)
end

function format_bytes(bytes)
    bytes < 1024 && return @sprintf("%d B", bytes)
    bytes < 1024^2 && return @sprintf("%.1f KB", bytes/1024)
    bytes < 1024^3 && return @sprintf("%.1f MB", bytes/1024^2)
    return @sprintf("%.2f GB", bytes/1024^3)
end

# ── Benchmark 1: GPU 4D Allocation Cost ───────────────────────────────────────

function bench_allocation(NquadN, nSpec, nRaman; FT=Float64)
    println("\n" * "="^80)
    println("BENCHMARK 1: GPU 4D Tensor Allocation")
    println("="^80)

    mat_4d_bytes = NquadN * NquadN * nSpec * nRaman * sizeof(FT)
    src_4d_bytes = NquadN * 1     * nSpec * nRaman * sizeof(FT)
    mat_3d_bytes = NquadN * NquadN * nSpec * sizeof(FT)

    println("Dimensions: NquadN=$NquadN, nSpec=$nSpec, nRaman=$nRaman, FT=$FT")
    println("Single 4D matrix: $(format_bytes(mat_4d_bytes))")
    println("Single 4D source: $(format_bytes(src_4d_bytes))")
    println("Single 3D matrix: $(format_bytes(mat_3d_bytes))")

    # Count arrays as in vSmartMOM: AddedLayerRS + CompositeLayerRS
    # Each has 4 matrices + 2 sources (4D), plus elastic 4 matrices + 2 sources (3D)
    n_4d_matrices = 8   # 4 per added + 4 per composite
    n_4d_sources  = 4   # 2 per added + 2 per composite
    n_3d_matrices = 8
    n_3d_sources  = 4

    total_4d = n_4d_matrices * mat_4d_bytes + n_4d_sources * src_4d_bytes
    total_3d = n_3d_matrices * mat_3d_bytes + n_3d_sources * NquadN * nSpec * sizeof(FT)

    println("\nPer-layer GPU memory budget:")
    println("  4D tensors (Raman): $(format_bytes(total_4d)) " *
            "($(n_4d_matrices) matrices + $(n_4d_sources) sources)")
    println("  3D tensors (elastic): $(format_bytes(total_3d))")
    println("  Ratio 4D/3D: $(@sprintf("%.0fx", total_4d / total_3d))")

    mem0 = gpu_memory_gb()

    # Time allocation of all 4D arrays for one layer pair
    CUDA.synchronize()
    t_alloc = @elapsed begin
        ie_mats = [CUDA.zeros(FT, NquadN, NquadN, nSpec, nRaman) for _ in 1:n_4d_matrices]
        ie_srcs = [CUDA.zeros(FT, NquadN, 1, nSpec, nRaman) for _ in 1:n_4d_sources]
        CUDA.synchronize()
    end

    mem1 = gpu_memory_gb()

    println("\nAllocation time: $(@sprintf("%.3f s", t_alloc))")
    println("GPU memory: $(format_bytes(total_4d)) allocated, " *
            "GPU used: $(@sprintf("%.2f", mem1.used - mem0.used)) GB")

    # Clean up
    ie_mats = nothing; ie_srcs = nothing
    GC.gc(); CUDA.reclaim()

    return total_4d
end

# ── Benchmark 2: Batched MatMul Throughput ────────────────────────────────────

function bench_batched_matmul(NquadN, nSpec; FT=Float64, nreps=20)
    println("\n" * "="^80)
    println("BENCHMARK 2: Batched Matrix Multiply (⊠) Throughput")
    println("="^80)

    A = CUDA.rand(FT, NquadN, NquadN, nSpec)
    B = CUDA.rand(FT, NquadN, NquadN, nSpec)
    C = similar(A)
    CUDA.synchronize()

    # Warmup
    C .= A ⊠ B; CUDA.synchronize()

    times = Float64[]
    for _ in 1:nreps
        CUDA.synchronize()
        t = @elapsed begin
            C .= A ⊠ B
            CUDA.synchronize()
        end
        push!(times, t)
    end

    flops_per = 2 * NquadN^2 * NquadN * nSpec  # 2*n³*batch
    med_time = median(times)
    gflops = flops_per / med_time / 1e9

    println("Shape: ($NquadN,$NquadN,$nSpec) ⊠ ($NquadN,$NquadN,$nSpec)")
    println("Median time: $(@sprintf("%.4f ms", med_time*1e3)) | $(format_bytes(3*NquadN*NquadN*nSpec*sizeof(FT))) touched")
    println("Throughput: $(@sprintf("%.1f GFLOPS", gflops))")

    # Now test with 4D slice extracted via view (mimics ie_array[:,:,n₁,Δn])
    println("\n--- Slice ⊠ (3D × 4D-slice) pattern ---")
    D4 = CUDA.rand(FT, NquadN, NquadN, nSpec, 10)  # small nRaman chunk
    CUDA.synchronize()

    # Single slice multiply: A[:,:,n₁] ⊠ D4[:,:,n₁,Δn]
    # In the real code this happens element-wise over n₁ inside the Δn loop
    # But actually they do 3D ⊠ 3D-view-of-4D

    # Test: extract a 3D view from 4D and multiply
    slice_times = Float64[]
    for _ in 1:nreps
        CUDA.synchronize()
        t = @elapsed begin
            @views C .= A ⊠ D4[:,:,:,1]
            CUDA.synchronize()
        end
        push!(slice_times, t)
    end

    println("3D ⊠ 4D-view: $(@sprintf("%.4f ms", median(slice_times)*1e3))")

    # Test: contiguous 3D copy then multiply (to see if view penalty exists)
    D3_copy = similar(A)
    copy_mul_times = Float64[]
    for _ in 1:nreps
        CUDA.synchronize()
        t = @elapsed begin
            D3_copy .= @view D4[:,:,:,1]
            C .= A ⊠ D3_copy
            CUDA.synchronize()
        end
        push!(copy_mul_times, t)
    end

    println("Copy+mul:     $(@sprintf("%.4f ms", median(copy_mul_times)*1e3))")

    # Cleanup
    A = B = C = D4 = D3_copy = nothing; GC.gc(); CUDA.reclaim()

    return med_time
end

# ── Benchmark 3: CPU↔GPU Transfer vs On-GPU Slice Access ─────────────────────

function bench_transfer_vs_gpu(NquadN, nSpec, nRaman; FT=Float64, nreps=10)
    println("\n" * "="^80)
    println("BENCHMARK 3: CPU↔GPU Transfer Cost vs On-GPU 4D Access")
    println("="^80)

    slice_bytes = NquadN * NquadN * nSpec * sizeof(FT)
    println("Per-slice size: $(format_bytes(slice_bytes))")
    println("Full 4D size:   $(format_bytes(slice_bytes * nRaman))")

    # === Strategy A: Full 4D on GPU, access slices via views ===
    println("\n--- Strategy A: All 4D on GPU ---")

    mem0 = gpu_memory_gb()
    gpu_4d = CUDA.rand(FT, NquadN, NquadN, nSpec, nRaman)
    gpu_buf = CUDA.zeros(FT, NquadN, NquadN, nSpec)  # result buffer
    gpu_A   = CUDA.rand(FT, NquadN, NquadN, nSpec)    # elastic matrix
    CUDA.synchronize()
    mem1 = gpu_memory_gb()

    println("GPU memory for 4D array: $(@sprintf("%.2f GB", mem1.used - mem0.used))")

    # Time: loop over all Δn, do a batched_mul with each slice
    CUDA.synchronize()
    t_gpu = @elapsed begin
        for Δn in 1:nRaman
            @views gpu_buf .= gpu_A ⊠ gpu_4d[:,:,:,Δn]
        end
        CUDA.synchronize()
    end

    println("Time (loop $nRaman slices, ⊠ each): $(@sprintf("%.3f s", t_gpu))")

    gpu_4d = nothing; GC.gc(); CUDA.reclaim()

    # === Strategy B: 4D on CPU, transfer slice to GPU per Δn ===
    println("\n--- Strategy B: 4D on CPU, transfer per Δn ---")

    cpu_4d = rand(FT, NquadN, NquadN, nSpec, nRaman)
    gpu_slice = CUDA.zeros(FT, NquadN, NquadN, nSpec)  # reusable GPU buffer
    CUDA.synchronize()

    mem2 = gpu_memory_gb()
    println("GPU memory (only buffers): $(@sprintf("%.2f GB", mem2.used - mem0.used))")

    CUDA.synchronize()
    t_staged = @elapsed begin
        for Δn in 1:nRaman
            h2d_slice!(gpu_slice, cpu_4d, Δn)
            gpu_buf .= gpu_A ⊠ gpu_slice
            CUDA.synchronize()  # must sync before overwriting gpu_slice
        end
    end

    println("Time (loop $nRaman slices, H2D + ⊠ each): $(@sprintf("%.3f s", t_staged))")

    # === Strategy C: Pinned CPU memory + async transfer ===
    println("\n--- Strategy C: Pinned CPU + async H2D ---")

    # Allocate pinned host memory for faster transfers
    cpu_pinned = CUDA.Mem.pin(rand(FT, NquadN, NquadN, nSpec, nRaman))
    CUDA.synchronize()

    t_pinned = @elapsed begin
        for Δn in 1:nRaman
            h2d_slice!(gpu_slice, cpu_pinned, Δn)
            gpu_buf .= gpu_A ⊠ gpu_slice
            CUDA.synchronize()
        end
    end

    println("Time (loop $nRaman, pinned H2D + ⊠): $(@sprintf("%.3f s", t_pinned))")

    # === Strategy D: Double-buffered async ===
    println("\n--- Strategy D: Double-buffered async transfer ---")

    gpu_slice2 = CUDA.zeros(FT, NquadN, NquadN, nSpec)
    buffers = (gpu_slice, gpu_slice2)
    CUDA.synchronize()

    t_double = @elapsed begin
        # First slice synchronously
        h2d_slice!(buffers[1], cpu_pinned, 1)
        CUDA.synchronize()

        for Δn in 1:nRaman
            cur = mod1(Δn, 2)
            nxt = mod1(Δn + 1, 2)

            # Start next transfer while computing current
            if Δn < nRaman
                h2d_slice!(buffers[nxt], cpu_pinned, Δn + 1)
            end
            gpu_buf .= gpu_A ⊠ buffers[cur]
            CUDA.synchronize()
        end
    end

    println("Time (loop $nRaman, double-buffered): $(@sprintf("%.3f s", t_double))")

    # === Strategy E: Tuple of 3D GPU arrays (no 4D allocation at all) ===
    println("\n--- Strategy E: Tuple/Vector of 3D CuArrays (no 4D tensor) ---")

    mem4 = gpu_memory_gb()
    # Instead of one (NquadN, NquadN, nSpec, nRaman) 4D array, store as
    # Vector{CuArray{FT,3}} of length nRaman, each (NquadN, NquadN, nSpec)
    gpu_slices = [CUDA.rand(FT, NquadN, NquadN, nSpec) for _ in 1:nRaman]
    CUDA.synchronize()
    mem5 = gpu_memory_gb()
    println("GPU memory (vector of 3D): $(@sprintf("%.2f GB", mem5.used - mem4.used))")
    println("  vs contiguous 4D:        $(@sprintf("%.2f GB", NquadN*NquadN*nSpec*nRaman*sizeof(FT)/1e9))")

    CUDA.synchronize()
    t_tuple = @elapsed begin
        for Δn in 1:nRaman
            gpu_buf .= gpu_A ⊠ gpu_slices[Δn]
        end
        CUDA.synchronize()
    end

    println("Time (loop $nRaman slices via indexing): $(@sprintf("%.3f s", t_tuple))")
    gpu_slices = nothing; GC.gc(); CUDA.reclaim()

    # === Strategy F: Tuple of 3D CPU arrays (avoid 4D entirely) ===
    println("\n--- Strategy F: Vector of 3D CPU arrays + H2D per Δn ---")

    cpu_slices = [rand(FT, NquadN, NquadN, nSpec) for _ in 1:nRaman]
    gpu_tmp = CUDA.zeros(FT, NquadN, NquadN, nSpec)
    CUDA.synchronize()

    t_cpu_tuple = @elapsed begin
        for Δn in 1:nRaman
            copyto!(gpu_tmp, cpu_slices[Δn])  # contiguous — no SubArray issue
            gpu_buf .= gpu_A ⊠ gpu_tmp
            CUDA.synchronize()
        end
    end

    println("Time (loop $nRaman, contiguous H2D): $(@sprintf("%.3f s", t_cpu_tuple))")
    cpu_slices = nothing; gpu_tmp = nothing; GC.gc(); CUDA.reclaim()

    # === Summary ===
    println("\n--- Transfer Overhead Summary ---")
    overhead = t_staged / t_gpu
    println("Strategy A (all GPU, 4D):       $(@sprintf("%.3f s", t_gpu)) (baseline)")
    println("Strategy B (CPU 4D, slice H2D): $(@sprintf("%.3f s", t_staged)) ($(@sprintf("%.1fx", overhead)) of GPU)")
    println("Strategy C (pinned CPU 4D):     $(@sprintf("%.3f s", t_pinned)) ($(@sprintf("%.1fx", t_pinned/t_gpu)) of GPU)")
    println("Strategy D (double-buffered):   $(@sprintf("%.3f s", t_double)) ($(@sprintf("%.1fx", t_double/t_gpu)) of GPU)")
    println("Strategy E (Vec{CuArray3D}):    $(@sprintf("%.3f s", t_tuple)) ($(@sprintf("%.1fx", t_tuple/t_gpu)) of GPU)")
    println("Strategy F (Vec{CPU3D}+H2D):    $(@sprintf("%.3f s", t_cpu_tuple)) ($(@sprintf("%.1fx", t_cpu_tuple/t_gpu)) of GPU)")

    saved_gb = NquadN * NquadN * nSpec * nRaman * sizeof(FT) / 1e9
    println("\nGPU memory saved by CPU staging: $(@sprintf("%.2f GB", saved_gb)) per 4D tensor")
    println("With 8 such tensors: $(@sprintf("%.2f GB", 8*saved_gb)) saved")

    # Cleanup
    gpu_buf = gpu_A = gpu_slice = gpu_slice2 = nothing
    cpu_4d = cpu_pinned = nothing
    GC.gc(); CUDA.reclaim()

    return (t_gpu, t_staged, t_pinned, t_double, t_tuple, t_cpu_tuple)
end

# ── Benchmark 4: Raman Doubling Loop Simulation ──────────────────────────────

function bench_doubling_simulation(NquadN, nSpec, nRaman, ndoubl; FT=Float64)
    println("\n" * "="^80)
    println("BENCHMARK 4: Simulated Raman Doubling (ndoubl=$ndoubl)")
    println("="^80)

    # Mimics the operations in doubling_inelastic.jl:
    # For each doubling step:
    #   1. batch_inv!(gp_refl, I - r ⊠ r)   [3D, elastic]
    #   2. tt_gp = t ⊠ gp_refl              [3D, elastic]
    #   3. For Δn = 1:nRaman:
    #        ie_src[:,:,n₁,Δn] += tt_gp[:,:,n₁] ⊠ (ie_t[:,:,n₁,Δn] +
    #                              ie_r[:,:,n₁,Δn] ⊠ r[:,:,n₀] ⊠ gp_refl[:,:,n₀] ⊠ t[:,:,n₀])
    #        ie_r[:,:,n₁,Δn] = updated
    #        ie_t[:,:,n₁,Δn] = updated

    # Build fake index mapping (n₁ → n₀ via Raman shift)
    # In real code: n₀ = n₁ + i_λ₁λ₀[Δn], and n₁ ranges over valid spectral points
    # For simplicity: all n₁ map to valid n₀ within bounds
    i_λ₁λ₀ = collect(1:nRaman)  # shift offsets

    println("Operations per doubling step per Δn:")
    println("  ~5 batched_mul + 2 additions on 3D slices of size ($NquadN,$NquadN,$nSpec)")
    println("Total ops: $ndoubl doublings × $nRaman Δn = $(ndoubl * nRaman) inner iterations")

    # === Full GPU approach ===
    println("\n--- Full GPU (4D allocated) ---")

    CUDA.reclaim()
    mem0 = gpu_memory_gb()

    # Elastic 3D arrays
    r = CUDA.rand(FT, NquadN, NquadN, nSpec)
    t_mat = CUDA.rand(FT, NquadN, NquadN, nSpec)
    gp_refl = similar(r)
    tt_gp = similar(r)
    I_static = repeat(CuArray(Matrix{FT}(I, NquadN, NquadN)), 1, 1, nSpec)

    # Inelastic 4D arrays (6 total: ier⁻⁺, iet⁺⁺, ier⁺⁻, iet⁻⁻, ieJ⁺, ieJ⁻)
    ie_r1 = CUDA.rand(FT, NquadN, NquadN, nSpec, nRaman)  # ier⁻⁺
    ie_t1 = CUDA.rand(FT, NquadN, NquadN, nSpec, nRaman)  # iet⁺⁺
    ie_r2 = CUDA.rand(FT, NquadN, NquadN, nSpec, nRaman)  # ier⁺⁻ (symmetric doubling)
    ie_t2 = CUDA.rand(FT, NquadN, NquadN, nSpec, nRaman)  # iet⁻⁻
    ie_J1 = CUDA.rand(FT, NquadN, 1, nSpec, nRaman)       # ieJ⁺
    ie_J2 = CUDA.rand(FT, NquadN, 1, nSpec, nRaman)       # ieJ⁻

    # Temp buffers (allocated once, as in real code)
    tmp1 = similar(r)
    tmp2 = similar(r)
    CUDA.synchronize()

    mem1 = gpu_memory_gb()
    println("GPU memory used: $(@sprintf("%.2f GB", mem1.used - mem0.used))")

    CUDA.synchronize()
    t_full_gpu = @elapsed begin
        for n in 1:ndoubl
            # Elastic part: inv and multiply (3D)
            gp_refl .= I_static .- r ⊠ r  # placeholder for batch_inv
            tt_gp .= t_mat ⊠ gp_refl

            # Inelastic Δn loop
            for Δn in 1:nRaman
                # In real code: n₀, n₁ = get_n₀_n₁(...)
                # Here: operate on full spectral dim as proxy
                @views begin
                    # ~5 batched muls per Δn (simplified from real code)
                    tmp1 .= ie_r1[:,:,:,Δn] ⊠ r
                    tmp2 .= tt_gp ⊠ (ie_t1[:,:,:,Δn] .+ tmp1 ⊠ gp_refl ⊠ t_mat)
                    ie_t1[:,:,:,Δn] .= tmp2

                    tmp1 .= ie_r1[:,:,:,Δn] .+ r ⊠ ie_r1[:,:,:,Δn] ⊠ gp_refl ⊠ r
                    ie_r1[:,:,:,Δn] .= tmp1
                end
            end
            CUDA.synchronize()
        end
    end

    println("Full GPU doubling time: $(@sprintf("%.3f s", t_full_gpu))")

    # Cleanup 4D
    ie_r1 = ie_t1 = ie_r2 = ie_t2 = ie_J1 = ie_J2 = nothing
    GC.gc(); CUDA.reclaim()

    # === CPU-staged approach ===
    println("\n--- CPU-Staged (4D on host, slice to GPU per Δn) ---")

    mem2 = gpu_memory_gb()

    # 4D on CPU (pinned)
    cpu_ie_r1 = CUDA.Mem.pin(rand(FT, NquadN, NquadN, nSpec, nRaman))
    cpu_ie_t1 = CUDA.Mem.pin(rand(FT, NquadN, NquadN, nSpec, nRaman))

    # GPU slice buffers (reused each Δn)
    gpu_ie_r_slice = CUDA.zeros(FT, NquadN, NquadN, nSpec)
    gpu_ie_t_slice = CUDA.zeros(FT, NquadN, NquadN, nSpec)
    CUDA.synchronize()

    mem3 = gpu_memory_gb()
    println("GPU memory used (staged): $(@sprintf("%.2f GB", mem3.used - mem2.used))")

    CUDA.synchronize()
    t_staged = @elapsed begin
        for n in 1:ndoubl
            gp_refl .= I_static .- r ⊠ r
            tt_gp .= t_mat ⊠ gp_refl

            for Δn in 1:nRaman
                # H2D: bring slice to GPU
                h2d_slice!(gpu_ie_r_slice, cpu_ie_r1, Δn)
                h2d_slice!(gpu_ie_t_slice, cpu_ie_t1, Δn)
                CUDA.synchronize()

                # Compute (same ops as above)
                tmp1 .= gpu_ie_r_slice ⊠ r
                tmp2 .= tt_gp ⊠ (gpu_ie_t_slice .+ tmp1 ⊠ gp_refl ⊠ t_mat)
                gpu_ie_t_slice .= tmp2

                tmp1 .= gpu_ie_r_slice .+ r ⊠ gpu_ie_r_slice ⊠ gp_refl ⊠ r
                gpu_ie_r_slice .= tmp1

                CUDA.synchronize()

                # D2H: write back
                d2h_slice!(cpu_ie_r1, gpu_ie_r_slice, Δn)
                d2h_slice!(cpu_ie_t1, gpu_ie_t_slice, Δn)
                CUDA.synchronize()
            end
        end
    end

    println("CPU-staged doubling time: $(@sprintf("%.3f s", t_staged))")
    println("Overhead vs full GPU: $(@sprintf("%.1fx", t_staged / t_full_gpu))")

    gpu_ie_r_slice = gpu_ie_t_slice = nothing
    cpu_ie_r1 = cpu_ie_t1 = nothing
    GC.gc(); CUDA.reclaim()

    # === Vector of 3D GPU arrays (your tuple idea) ===
    println("\n--- Vector of 3D GPU arrays (tuple-style, all on GPU) ---")

    mem4 = gpu_memory_gb()

    # Elastic 3D (re-allocate since we cleaned up)
    r = CUDA.rand(FT, NquadN, NquadN, nSpec)
    t_mat = CUDA.rand(FT, NquadN, NquadN, nSpec)
    gp_refl = similar(r)
    tt_gp = similar(r)
    I_static = repeat(CuArray(Matrix{FT}(I, NquadN, NquadN)), 1, 1, nSpec)
    tmp1 = similar(r)
    tmp2 = similar(r)

    # Inelastic as Vector{CuArray{FT,3}} — same total memory but no 4D contiguous alloc
    vec_ie_r = [CUDA.rand(FT, NquadN, NquadN, nSpec) for _ in 1:nRaman]
    vec_ie_t = [CUDA.rand(FT, NquadN, NquadN, nSpec) for _ in 1:nRaman]
    CUDA.synchronize()

    mem5 = gpu_memory_gb()
    println("GPU memory used (vec of 3D): $(@sprintf("%.2f GB", mem5.used - mem4.used))")

    CUDA.synchronize()
    t_vec_gpu = @elapsed begin
        for n in 1:ndoubl
            gp_refl .= I_static .- r ⊠ r
            tt_gp .= t_mat ⊠ gp_refl

            for Δn in 1:nRaman
                # Direct indexing — no view extraction from 4D
                tmp1 .= vec_ie_r[Δn] ⊠ r
                tmp2 .= tt_gp ⊠ (vec_ie_t[Δn] .+ tmp1 ⊠ gp_refl ⊠ t_mat)
                vec_ie_t[Δn] .= tmp2

                tmp1 .= vec_ie_r[Δn] .+ r ⊠ vec_ie_r[Δn] ⊠ gp_refl ⊠ r
                vec_ie_r[Δn] .= tmp1
            end
            CUDA.synchronize()
        end
    end

    println("Vec{CuArray3D} doubling time: $(@sprintf("%.3f s", t_vec_gpu))")
    println("vs full GPU 4D: $(@sprintf("%.1fx", t_vec_gpu / t_full_gpu))")
    println("vs CPU staged:  $(@sprintf("%.1fx", t_vec_gpu / t_staged))")

    # Cleanup
    r = t_mat = gp_refl = tt_gp = I_static = tmp1 = tmp2 = nothing
    vec_ie_r = vec_ie_t = nothing
    GC.gc(); CUDA.reclaim()

    return (t_full_gpu, t_staged, t_vec_gpu)
end

# ── Benchmark 5: Allocation Churn Impact ──────────────────────────────────────

function bench_allocation_churn(NquadN, nSpec, nRaman; FT=Float64)
    println("\n" * "="^80)
    println("BENCHMARK 5: Allocation Churn (similar() in loops)")
    println("="^80)

    A = CUDA.rand(FT, NquadN, NquadN, nSpec)
    B = CUDA.rand(FT, NquadN, NquadN, nSpec)
    CUDA.synchronize()

    niter = nRaman * 3  # simulate 3 Δn loops

    # With allocation each iteration (mimics current `tmp = A ⊠ B` pattern)
    CUDA.synchronize()
    t_alloc = @elapsed begin
        for i in 1:niter
            tmp = A ⊠ B  # allocates new array each time
        end
        CUDA.synchronize()
    end

    # With pre-allocated buffer
    C = similar(A)
    CUDA.synchronize()
    t_prealloc = @elapsed begin
        for i in 1:niter
            C .= A ⊠ B  # reuses buffer
        end
        CUDA.synchronize()
    end

    println("$niter iterations of batched_mul ($NquadN,$NquadN,$nSpec):")
    println("  With allocation each time: $(@sprintf("%.3f s", t_alloc)) " *
            "($(@sprintf("%.3f ms", t_alloc/niter*1e3))/iter)")
    println("  Pre-allocated buffer:      $(@sprintf("%.3f s", t_prealloc)) " *
            "($(@sprintf("%.3f ms", t_prealloc/niter*1e3))/iter)")
    println("  Speedup: $(@sprintf("%.1fx", t_alloc / t_prealloc))")

    # Measure allocation cost alone
    CUDA.synchronize()
    t_alloc_only = @elapsed begin
        for i in 1:niter
            tmp = similar(A)
        end
        CUDA.synchronize()
    end

    println("\n  Allocation-only cost: $(@sprintf("%.3f s", t_alloc_only)) for $niter similar() calls")
    println("  Per call: $(@sprintf("%.3f ms", t_alloc_only/niter*1e3))")

    A = B = C = nothing; GC.gc(); CUDA.reclaim()
    return (t_alloc, t_prealloc)
end

# ── Benchmark 6: Synchronization Overhead ─────────────────────────────────────

function bench_sync_overhead(NquadN, nSpec, nRaman; FT=Float64)
    println("\n" * "="^80)
    println("BENCHMARK 6: Synchronization Barrier Overhead")
    println("="^80)

    A = CUDA.rand(FT, NquadN, NquadN, nSpec)
    B = CUDA.rand(FT, NquadN, NquadN, nSpec)
    C = similar(A)
    CUDA.synchronize()

    niter = nRaman

    # Sync after every operation (mimics current synchronize_if_gpu pattern)
    CUDA.synchronize()
    t_heavy_sync = @elapsed begin
        for i in 1:niter
            C .= A ⊠ B
            CUDA.synchronize()
            C .= C .+ A
            CUDA.synchronize()
            C .= A ⊠ C
            CUDA.synchronize()
        end
    end

    # Sync only at end of each iteration
    CUDA.synchronize()
    t_light_sync = @elapsed begin
        for i in 1:niter
            C .= A ⊠ B
            C .= C .+ A
            C .= A ⊠ C
            CUDA.synchronize()
        end
    end

    # Sync only at end of full loop
    CUDA.synchronize()
    t_batch_sync = @elapsed begin
        for i in 1:niter
            C .= A ⊠ B
            C .= C .+ A
            C .= A ⊠ C
        end
        CUDA.synchronize()
    end

    println("$niter iterations, 3 ops each:")
    println("  Sync after every op:  $(@sprintf("%.3f s", t_heavy_sync))")
    println("  Sync after each iter: $(@sprintf("%.3f s", t_light_sync))")
    println("  Sync only at end:     $(@sprintf("%.3f s", t_batch_sync))")
    println("  Heavy/Light ratio:    $(@sprintf("%.1fx", t_heavy_sync / t_light_sync))")
    println("  Heavy/Batch ratio:    $(@sprintf("%.1fx", t_heavy_sync / t_batch_sync))")

    A = B = C = nothing; GC.gc(); CUDA.reclaim()
    return (t_heavy_sync, t_light_sync, t_batch_sync)
end

# ── Benchmark 7: Mixed Precision ──────────────────────────────────────────────

function bench_mixed_precision(NquadN, nSpec, nRaman; nreps=20)
    println("\n" * "="^80)
    println("BENCHMARK 7: Float64 vs Float32 Throughput")
    println("="^80)

    for FT in (Float64, Float32)
        A = CUDA.rand(FT, NquadN, NquadN, nSpec)
        B = CUDA.rand(FT, NquadN, NquadN, nSpec)
        C = similar(A)
        CUDA.synchronize()

        # Warmup
        C .= A ⊠ B; CUDA.synchronize()

        times = Float64[]
        for _ in 1:nreps
            CUDA.synchronize()
            t = @elapsed begin
                C .= A ⊠ B
                CUDA.synchronize()
            end
            push!(times, t)
        end

        mem = NquadN * NquadN * nSpec * sizeof(FT) * 3
        println("$FT ⊠: $(@sprintf("%.4f ms", median(times)*1e3)) | " *
                "memory footprint: $(format_bytes(mem))")

        A = B = C = nothing; GC.gc(); CUDA.reclaim()
    end

    # Also test 4D allocation size difference
    for FT in (Float64, Float32)
        sz = NquadN * NquadN * nSpec * nRaman * sizeof(FT) * 8 +
             NquadN * 1 * nSpec * nRaman * sizeof(FT) * 4
        println("$FT total 4D (8 mat + 4 src): $(format_bytes(sz))")
    end
end

# ── Benchmark 8: Interaction pattern (two-operand cross-dim multiply) ─────────

function bench_interaction_pattern(NquadN, nSpec, nRaman; FT=Float64)
    println("\n" * "="^80)
    println("BENCHMARK 8: Interaction Pattern (3D elastic × 4D-slice inelastic)")
    println("="^80)

    # Mimics interaction_inelastic.jl:
    # T01_inv[:,:,n₁] ⊠ (ier[:,:,n₁,Δn] ⊠ T[:,:,n₀] + r[:,:,n₁] ⊠ ieR[:,:,n₁,Δn])
    # With tmp arrays allocated once, looped over Δn

    # Elastic 3D
    T_inv = CUDA.rand(FT, NquadN, NquadN, nSpec)
    T_fwd = CUDA.rand(FT, NquadN, NquadN, nSpec)
    r_el  = CUDA.rand(FT, NquadN, NquadN, nSpec)

    # Inelastic 4D
    ie_r  = CUDA.rand(FT, NquadN, NquadN, nSpec, nRaman)
    ie_R  = CUDA.rand(FT, NquadN, NquadN, nSpec, nRaman)

    # Temps
    tmp1  = similar(T_inv)
    tmp2  = similar(T_inv)
    out   = CUDA.zeros(FT, NquadN, NquadN, nSpec, nRaman)
    CUDA.synchronize()

    mem = gpu_memory_gb()
    println("GPU memory used: $(@sprintf("%.2f GB", mem.used))")

    CUDA.synchronize()
    t_interaction = @elapsed begin
        for Δn in 1:nRaman
            @views begin
                tmp1 .= ie_r[:,:,:,Δn] ⊠ T_fwd
                tmp2 .= r_el ⊠ ie_R[:,:,:,Δn]
                tmp1 .+= tmp2
                out[:,:,:,Δn] .= T_inv ⊠ tmp1
            end
        end
        CUDA.synchronize()
    end

    println("Interaction loop ($nRaman iterations): $(@sprintf("%.3f s", t_interaction))")
    println("Per Δn: $(@sprintf("%.3f ms", t_interaction/nRaman*1e3))")
    println("Batched muls per Δn: 3 | Total: $(3*nRaman)")

    ie_r = ie_R = out = nothing; GC.gc(); CUDA.reclaim()
    return t_interaction
end

# ── Main ──────────────────────────────────────────────────────────────────────

function main()
    preset = parse_preset()
    (; NquadN, nSpec, nRaman, ndoubl, nLayers, nFourier) = preset

    println("╔" * "═"^78 * "╗")
    println("║  vSmartMOM Raman Memory & Compute Bottleneck Benchmark" * " "^23 * "║")
    println("╠" * "═"^78 * "╣")
    println("║  NquadN=$NquadN  nSpec=$(@sprintf("%-5d", nSpec))  nRaman=$(@sprintf("%-4d", nRaman))  " *
            "ndoubl=$ndoubl  nLayers=$nLayers  nFourier=$nFourier" *
            " "^(78-65-length("$nSpec")-length("$nRaman")) * "║")

    dev = CUDA.device()
    println("║  GPU: $(CUDA.name(dev))" * " "^max(0, 78-7-length(string(CUDA.name(dev)))) * "║")
    mem = gpu_memory_gb()
    println("║  VRAM: $(@sprintf("%.1f GB total, %.1f GB free", mem.total, mem.free))" * " "^38 * "║")
    println("╚" * "═"^78 * "╝")

    # Run all benchmarks
    b1_bytes = bench_allocation(NquadN, nSpec, nRaman)
    b2_time  = bench_batched_matmul(NquadN, nSpec)
    b3_times = bench_transfer_vs_gpu(NquadN, nSpec, nRaman)
    b4_times = bench_doubling_simulation(NquadN, nSpec, nRaman, ndoubl)
    b5_times = bench_allocation_churn(NquadN, nSpec, nRaman)
    b6_times = bench_sync_overhead(NquadN, nSpec, nRaman)
    bench_mixed_precision(NquadN, nSpec, nRaman)
    b8_time  = bench_interaction_pattern(NquadN, nSpec, nRaman)

    # ── Summary ───────────────────────────────────────────────────────────────
    println("\n" * "="^80)
    println("SUMMARY: Bottleneck Analysis")
    println("="^80)

    # Estimate full RT cost breakdown
    doubling_per_layer = b4_times[1]  # full GPU doubling time
    interaction_per_layer = b8_time
    total_est = nLayers * nFourier * (doubling_per_layer + interaction_per_layer)

    println("\nEstimated RT kernel time ($(nLayers) layers × $(nFourier) Fourier):")
    println("  Doubling:    $(@sprintf("%.1f s", nLayers * nFourier * doubling_per_layer))")
    println("  Interaction: $(@sprintf("%.1f s", nLayers * nFourier * interaction_per_layer))")
    println("  Total:       $(@sprintf("%.1f s", total_est))")

    println("\nMemory:")
    println("  4D Raman tensors: $(format_bytes(b1_bytes)) per layer pair")
    println("  GPU savings with CPU staging: $(format_bytes(b1_bytes)) freed (keep only 3D slices)")

    println("\nTransfer overhead (slice access strategies):")
    gpu_t, staged_t, pinned_t, double_t, tuple_t, cpu_tuple_t = b3_times
    println("  A: All GPU 4D:              $(@sprintf("%.3f s", gpu_t)) (baseline)")
    println("  B: CPU 4D → H2D slices:     $(@sprintf("%.3f s", staged_t)) ($(@sprintf("%.1fx", staged_t/gpu_t)))")
    println("  C: Pinned CPU 4D:           $(@sprintf("%.3f s", pinned_t)) ($(@sprintf("%.1fx", pinned_t/gpu_t)))")
    println("  D: Double-buffered:         $(@sprintf("%.3f s", double_t)) ($(@sprintf("%.1fx", double_t/gpu_t)))")
    println("  E: Vec{CuArray3D} on GPU:   $(@sprintf("%.3f s", tuple_t)) ($(@sprintf("%.1fx", tuple_t/gpu_t)))")
    println("  F: Vec{CPU3D} → H2D:        $(@sprintf("%.3f s", cpu_tuple_t)) ($(@sprintf("%.1fx", cpu_tuple_t/gpu_t)))")

    println("\nDoubling loop strategies:")
    t_dbl_gpu, t_dbl_staged, t_dbl_vec = b4_times
    println("  Full GPU 4D:     $(@sprintf("%.3f s", t_dbl_gpu)) (baseline)")
    println("  CPU staged:      $(@sprintf("%.3f s", t_dbl_staged)) ($(@sprintf("%.1fx", t_dbl_staged/t_dbl_gpu)))")
    println("  Vec{CuArray3D}:  $(@sprintf("%.3f s", t_dbl_vec)) ($(@sprintf("%.1fx", t_dbl_vec/t_dbl_gpu)))")

    println("\nAllocation churn:")
    println("  Speedup from pre-allocation: $(@sprintf("%.1fx", b5_times[1]/b5_times[2]))")

    println("\nSync barriers:")
    println("  Heavy→Light sync: $(@sprintf("%.1fx", b6_times[1]/b6_times[2])) speedup")
    println("  Heavy→Batch sync: $(@sprintf("%.1fx", b6_times[1]/b6_times[3])) speedup")

    println("\n" * "="^80)
    println("RECOMMENDATIONS (ranked by expected impact)")
    println("="^80)

    staging_penalty = staged_t / gpu_t
    vec_gpu_ratio = tuple_t / gpu_t
    alloc_ratio = b5_times[1] / b5_times[2]
    sync_ratio = b6_times[1] / b6_times[2]

    println("""
  1. Vec{CuArray{FT,3}} instead of CuArray{FT,4}:
     Same GPU memory, but avoids contiguous 4D allocation (which took $(format_bytes(b1_bytes))).
     Performance: $(@sprintf("%.1fx", vec_gpu_ratio)) vs 4D baseline — essentially free if ~1.0x.
     Benefit: Easier incremental allocation, no SubArray/view overhead.

  2. Pre-allocate temporaries (eliminate similar() in hot loops):
     $(@sprintf("%.1fx", alloc_ratio)) speedup measured. Drop-in fix.

  3. Reduce synchronization barriers:
     $(@sprintf("%.1fx", sync_ratio)) speedup (heavy→light sync).
     Sync once per doubling step or layer, not per operation.

  4. CPU staging (if GPU memory is the hard constraint):
     $(@sprintf("%.1fx", staging_penalty)) overhead for H2D/D2H per Δn.
     Saves $(format_bytes(b1_bytes)) GPU memory per 4D tensor ($(format_bytes(b1_bytes * 8)) total).
     $(staging_penalty < 2.0 ? "VIABLE — moderate overhead" : "EXPENSIVE — consider only when OOM").

  5. Mixed precision (Float32 for Raman 4D):
     ~2x memory reduction + ~2x compute throughput on tensor cores.
     Must validate numerical accuracy for Raman accumulation.

  6. Hybrid: Vec{CuArray3D} on GPU for active layers, CPU-staged for dormant:
     Keep current + adjacent layers on GPU, stage rest on CPU/NVMe.
     Best of both: low latency for active work, bounded GPU footprint.
""")

    println("✓ Benchmark complete.")
end

main()
