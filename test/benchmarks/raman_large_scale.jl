#!/usr/bin/env julia
#=
Large-scale Raman memory benchmark:
  NquadN=30, nSpec=30000, nRaman=100, ndoubl=8

Tests strategies for handling 4D tensors that may exceed GPU memory.
=#

using CUDA, NNlib, LinearAlgebra, Printf, Statistics

const ⊠ = NNlib.batched_mul

function gpu_mem()
    free, total = CUDA.Mem.info()
    return (used=(total-free)/1e9, free=free/1e9, total=total/1e9)
end

fb(bytes) = bytes < 1024^2 ? @sprintf("%.1f KB",bytes/1024) :
            bytes < 1024^3 ? @sprintf("%.1f MB",bytes/1024^2) :
                             @sprintf("%.2f GB",bytes/1024^3)

function h2d!(gpu, cpu4d, Δn)
    n1,n2,n3,_ = size(cpu4d)
    src = unsafe_wrap(Array, pointer(cpu4d, (Δn-1)*n1*n2*n3+1), (n1,n2,n3))
    copyto!(gpu, src)
end

function d2h!(cpu4d, gpu, Δn)
    n1,n2,n3,_ = size(cpu4d)
    dst = unsafe_wrap(Array, pointer(cpu4d, (Δn-1)*n1*n2*n3+1), (n1,n2,n3))
    copyto!(dst, gpu)
end

# ── Parameters ──
NquadN = 30; nSpec = 30000; nRaman = 100; ndoubl = 8; FT = Float64

slice_bytes = NquadN * NquadN * nSpec * sizeof(FT)
total_4d    = slice_bytes * nRaman

println("╔══════════════════════════════════════════════════════════════════╗")
println("║  Large-Scale Raman Memory Benchmark                            ║")
println("╠══════════════════════════════════════════════════════════════════╣")
@printf("║  NquadN=%d  nSpec=%d  nRaman=%d  ndoubl=%d\n", NquadN, nSpec, nRaman, ndoubl)
println("║  GPU: $(CUDA.name(CUDA.device()))")
m = gpu_mem()
@printf("║  VRAM: %.1f GB total, %.1f GB free\n", m.total, m.free)
println("║  Single 3D slice: $(fb(slice_bytes))")
println("║  Full nRaman×3D:  $(fb(total_4d))")
println("║  18 tensors peak: $(fb(total_4d*18))")
println("╚══════════════════════════════════════════════════════════════════╝")

# Shared elastic 3D arrays (always on GPU)
r     = CUDA.rand(FT, NquadN, NquadN, nSpec)
t_mat = CUDA.rand(FT, NquadN, NquadN, nSpec)
gp    = similar(r); tt_gp = similar(r)
I_s   = repeat(CuArray(Matrix{FT}(I, NquadN, NquadN)), 1, 1, nSpec)
tmp1  = similar(r); tmp2 = similar(r)
CUDA.synchronize()

# ═══════════════════════════════════════════════════════════════════════
# Helper: run doubling loop with given access pattern
# ═══════════════════════════════════════════════════════════════════════
function doubling_loop!(get_r, get_t, set_r!, set_t!, nRaman, ndoubl, r, t_mat, gp, tt_gp, I_s, tmp1, tmp2)
    for n in 1:ndoubl
        gp .= I_s .- r ⊠ r
        tt_gp .= t_mat ⊠ gp
        for Δn in 1:nRaman
            gr = get_r(Δn); gt = get_t(Δn)
            tmp1 .= gr ⊠ r
            tmp2 .= tt_gp ⊠ (gt .+ tmp1 ⊠ gp ⊠ t_mat)
            set_t!(Δn, tmp2)
            tmp1 .= gr .+ r ⊠ gr ⊠ gp ⊠ r
            set_r!(Δn, tmp1)
        end
        CUDA.synchronize()
    end
end

# ═══════════════════════════════════════════════════════════════════════
# Strategy A: Full 4D CuArray
# ═══════════════════════════════════════════════════════════════════════
can_fit = total_4d * 2 < gpu_mem().free * 1e9 * 0.85
if can_fit
    println("\n--- Strategy A: Full 4D CuArray on GPU ---")
    GC.gc(); CUDA.reclaim(); m0 = gpu_mem()
    ie_r4 = CUDA.rand(FT, NquadN, NquadN, nSpec, nRaman)
    ie_t4 = CUDA.rand(FT, NquadN, NquadN, nSpec, nRaman)
    CUDA.synchronize(); m1 = gpu_mem()
    @printf("  GPU used: %.2f GB\n", m1.used - m0.used)

    CUDA.synchronize()
    t_A = @elapsed doubling_loop!(
        Δn -> @view(ie_r4[:,:,:,Δn]),
        Δn -> @view(ie_t4[:,:,:,Δn]),
        (Δn,v) -> (ie_r4[:,:,:,Δn] .= v),
        (Δn,v) -> (ie_t4[:,:,:,Δn] .= v),
        nRaman, ndoubl, r, t_mat, gp, tt_gp, I_s, tmp1, tmp2)
    @printf("  Doubling time: %.3f s\n", t_A)
    ie_r4 = ie_t4 = nothing; GC.gc(); CUDA.reclaim()
else
    println("\n--- Strategy A: SKIPPED (need $(fb(Int(total_4d*2))), free $(fb(Int(gpu_mem().free*1e9)))) ---")
    t_A = NaN
end

# ═══════════════════════════════════════════════════════════════════════
# Strategy B: Vector{CuArray{FT,3}} on GPU
# ═══════════════════════════════════════════════════════════════════════
println("\n--- Strategy B: Vector{CuArray{FT,3}} (all on GPU) ---")
GC.gc(); CUDA.reclaim(); m0 = gpu_mem()
vr = [CUDA.rand(FT, NquadN, NquadN, nSpec) for _ in 1:nRaman]
vt = [CUDA.rand(FT, NquadN, NquadN, nSpec) for _ in 1:nRaman]
CUDA.synchronize(); m1 = gpu_mem()
@printf("  GPU used: %.2f GB\n", m1.used - m0.used)

CUDA.synchronize()
t_B = @elapsed doubling_loop!(
    Δn -> vr[Δn], Δn -> vt[Δn],
    (Δn,v) -> (vr[Δn] .= v), (Δn,v) -> (vt[Δn] .= v),
    nRaman, ndoubl, r, t_mat, gp, tt_gp, I_s, tmp1, tmp2)
@printf("  Doubling time: %.3f s\n", t_B)
vr = vt = nothing; GC.gc(); CUDA.reclaim()

# ═══════════════════════════════════════════════════════════════════════
# Strategy C: CPU pinned 4D + H2D/D2H per Δn
# ═══════════════════════════════════════════════════════════════════════
println("\n--- Strategy C: CPU pinned 4D + H2D/D2H per Δn ---")
GC.gc(); CUDA.reclaim(); m0 = gpu_mem()
cpu_r4 = CUDA.Mem.pin(rand(FT, NquadN, NquadN, nSpec, nRaman))
cpu_t4 = CUDA.Mem.pin(rand(FT, NquadN, NquadN, nSpec, nRaman))
gb_r = CUDA.zeros(FT, NquadN, NquadN, nSpec)
gb_t = CUDA.zeros(FT, NquadN, NquadN, nSpec)
CUDA.synchronize(); m1 = gpu_mem()
@printf("  GPU used: %.2f GB  |  CPU pinned: %s\n", m1.used - m0.used, fb(2*total_4d))

CUDA.synchronize()
t_C = @elapsed begin
    for n in 1:ndoubl
        gp .= I_s .- r ⊠ r
        tt_gp .= t_mat ⊠ gp
        for Δn in 1:nRaman
            h2d!(gb_r, cpu_r4, Δn); h2d!(gb_t, cpu_t4, Δn); CUDA.synchronize()
            tmp1 .= gb_r ⊠ r
            tmp2 .= tt_gp ⊠ (gb_t .+ tmp1 ⊠ gp ⊠ t_mat)
            gb_t .= tmp2
            tmp1 .= gb_r .+ r ⊠ gb_r ⊠ gp ⊠ r
            gb_r .= tmp1
            CUDA.synchronize()
            d2h!(cpu_r4, gb_r, Δn); d2h!(cpu_t4, gb_t, Δn); CUDA.synchronize()
        end
    end
end
@printf("  Doubling time: %.3f s\n", t_C)
cpu_r4 = cpu_t4 = gb_r = gb_t = nothing; GC.gc(); CUDA.reclaim()

# ═══════════════════════════════════════════════════════════════════════
# Strategy D: CPU pinned Vec{Array3D} + contiguous H2D/D2H
# ═══════════════════════════════════════════════════════════════════════
println("\n--- Strategy D: CPU pinned Vec{Array3D} + contiguous H2D/D2H ---")
GC.gc(); CUDA.reclaim(); m0 = gpu_mem()
cpvr = [CUDA.Mem.pin(rand(FT, NquadN, NquadN, nSpec)) for _ in 1:nRaman]
cpvt = [CUDA.Mem.pin(rand(FT, NquadN, NquadN, nSpec)) for _ in 1:nRaman]
gb_r = CUDA.zeros(FT, NquadN, NquadN, nSpec)
gb_t = CUDA.zeros(FT, NquadN, NquadN, nSpec)
CUDA.synchronize(); m1 = gpu_mem()
@printf("  GPU used: %.2f GB\n", m1.used - m0.used)

CUDA.synchronize()
t_D = @elapsed begin
    for n in 1:ndoubl
        gp .= I_s .- r ⊠ r
        tt_gp .= t_mat ⊠ gp
        for Δn in 1:nRaman
            copyto!(gb_r, cpvr[Δn]); copyto!(gb_t, cpvt[Δn]); CUDA.synchronize()
            tmp1 .= gb_r ⊠ r
            tmp2 .= tt_gp ⊠ (gb_t .+ tmp1 ⊠ gp ⊠ t_mat)
            gb_t .= tmp2
            tmp1 .= gb_r .+ r ⊠ gb_r ⊠ gp ⊠ r
            gb_r .= tmp1
            CUDA.synchronize()
            copyto!(cpvr[Δn], gb_r); copyto!(cpvt[Δn], gb_t); CUDA.synchronize()
        end
    end
end
@printf("  Doubling time: %.3f s\n", t_D)
cpvr = cpvt = gb_r = gb_t = nothing; GC.gc(); CUDA.reclaim()

# ═══════════════════════════════════════════════════════════════════════
# Strategy E: Tiled — tile_size Δn on GPU, rest on CPU
# ═══════════════════════════════════════════════════════════════════════
for tile_size in [4, 16, 32]
    @printf("\n--- Strategy E (tile=%d): %d on GPU, rest on CPU ---\n", tile_size, tile_size)
    GC.gc(); CUDA.reclaim(); m0 = gpu_mem()

    cpu_r = [CUDA.Mem.pin(rand(FT, NquadN, NquadN, nSpec)) for _ in 1:nRaman]
    cpu_t = [CUDA.Mem.pin(rand(FT, NquadN, NquadN, nSpec)) for _ in 1:nRaman]
    gpu_rt = [CUDA.zeros(FT, NquadN, NquadN, nSpec) for _ in 1:tile_size]
    gpu_tt = [CUDA.zeros(FT, NquadN, NquadN, nSpec) for _ in 1:tile_size]
    CUDA.synchronize(); m1 = gpu_mem()
    @printf("  GPU used: %.2f GB  (vs %.2f GB all-on-GPU)\n", m1.used-m0.used, 2*nRaman*slice_bytes/1e9)

    CUDA.synchronize()
    t_E = @elapsed begin
        for n in 1:ndoubl
            gp .= I_s .- r ⊠ r
            tt_gp .= t_mat ⊠ gp

            for ts in 1:tile_size:nRaman
                te = min(ts + tile_size - 1, nRaman)
                # H2D tile
                for (li, Δn) in enumerate(ts:te)
                    copyto!(gpu_rt[li], cpu_r[Δn])
                    copyto!(gpu_tt[li], cpu_t[Δn])
                end
                CUDA.synchronize()
                # Compute tile
                for (li, _) in enumerate(ts:te)
                    tmp1 .= gpu_rt[li] ⊠ r
                    tmp2 .= tt_gp ⊠ (gpu_tt[li] .+ tmp1 ⊠ gp ⊠ t_mat)
                    gpu_tt[li] .= tmp2
                    tmp1 .= gpu_rt[li] .+ r ⊠ gpu_rt[li] ⊠ gp ⊠ r
                    gpu_rt[li] .= tmp1
                end
                CUDA.synchronize()
                # D2H tile
                for (li, Δn) in enumerate(ts:te)
                    copyto!(cpu_r[Δn], gpu_rt[li])
                    copyto!(cpu_t[Δn], gpu_tt[li])
                end
                CUDA.synchronize()
            end
        end
    end
    @printf("  Doubling time: %.3f s\n", t_E)
    cpu_r = cpu_t = gpu_rt = gpu_tt = nothing; GC.gc(); CUDA.reclaim()
end

# ═══════════════════════════════════════════════════════════════════════
# Float32 comparison
# ═══════════════════════════════════════════════════════════════════════
println("\n--- Float32 comparison: Vec{CuArray{Float32,3}} on GPU ---")
GC.gc(); CUDA.reclaim(); m0 = gpu_mem()
r32 = CUDA.rand(Float32, NquadN, NquadN, nSpec)
t32 = CUDA.rand(Float32, NquadN, NquadN, nSpec)
gp32 = similar(r32); tt32 = similar(r32)
I32 = repeat(CuArray(Matrix{Float32}(I, NquadN, NquadN)), 1, 1, nSpec)
tmp132 = similar(r32); tmp232 = similar(r32)
vr32 = [CUDA.rand(Float32, NquadN, NquadN, nSpec) for _ in 1:nRaman]
vt32 = [CUDA.rand(Float32, NquadN, NquadN, nSpec) for _ in 1:nRaman]
CUDA.synchronize(); m1 = gpu_mem()
@printf("  GPU used: %.2f GB (vs %.2f GB Float64)\n", m1.used-m0.used, 2*nRaman*slice_bytes/1e9)

CUDA.synchronize()
t_F32 = @elapsed doubling_loop!(
    Δn -> vr32[Δn], Δn -> vt32[Δn],
    (Δn,v) -> (vr32[Δn] .= v), (Δn,v) -> (vt32[Δn] .= v),
    nRaman, ndoubl, r32, t32, gp32, tt32, I32, tmp132, tmp232)
@printf("  Doubling time: %.3f s\n", t_F32)

# ═══════════════════════════════════════════════════════════════════════
# Summary table
# ═══════════════════════════════════════════════════════════════════════
println("\n" * "="^66)
println("SUMMARY (NquadN=$NquadN, nSpec=$nSpec, nRaman=$nRaman)")
println("="^66)
println("  Single 3D slice: $(fb(slice_bytes))")
println("  Full nRaman×3D:  $(fb(total_4d)) (×2 for r+t = $(fb(2*total_4d)))")
println("  18 tensors peak: $(fb(total_4d*18))")
println()
@printf("  %-40s  %8s  %10s\n", "Strategy", "Time (s)", "GPU (GB)")
println("  " * "-"^62)
if !isnan(t_A)
    @printf("  %-40s  %8.3f  %10.2f\n", "A: Full 4D CuArray", t_A, 2*total_4d/1e9)
end
@printf("  %-40s  %8.3f  %10.2f\n", "B: Vec{CuArray3D} all GPU", t_B, 2*nRaman*slice_bytes/1e9)
@printf("  %-40s  %8.3f  %10.2f\n", "C: CPU pinned 4D + H2D/D2H", t_C, 2*slice_bytes/1e9)
@printf("  %-40s  %8.3f  %10.2f\n", "D: CPU pinned Vec3D + H2D/D2H", t_D, 2*slice_bytes/1e9)
@printf("  %-40s  %8.3f  %10.2f\n", "E: Tiled (tile=16)", 0.0, 2*16*slice_bytes/1e9)
@printf("  %-40s  %8.3f  %10.2f\n", "B-F32: Vec{CuArray3D} Float32", t_F32, nRaman*NquadN^2*nSpec*4*2/1e9)
println()
if !isnan(t_A)
    @printf("  Vec3D vs 4D speedup: %.1fx\n", t_A / t_B)
end
@printf("  CPU-staged overhead vs Vec3D: %.1fx\n", t_D / t_B)
@printf("  Float32 vs Float64 speedup:  %.1fx\n", t_B / t_F32)
