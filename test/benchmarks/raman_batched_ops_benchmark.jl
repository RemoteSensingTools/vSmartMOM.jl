#=
Micro-benchmark: Batched GPU operations for Raman optimization
================================================================
Tests different approaches for the inner-loop matrix operations:
1. Broadcasting with ⊠ (current approach - creates temporaries)
2. NNlib.batched_mul! in-place on pre-allocated buffers
3. 4D array views vs explicit 3D copies
4. Loop of 2D mul! vs batched 3D operations

This informs the best strategy for optimizing doubling_inelastic.jl.

Typical dimensions from the Raman code:
- NquadN = 9 (l_trunc=5, Stokes_IQU)
- nSpec ≈ 5500 (spectral points per microwindow)
- n_Raman ≈ 963 (Raman transitions)
- The inner loop operates on 3D sub-arrays of shape (9, 9, ~nL) and (9, 1, ~nL)
  where nL ≤ nSpec varies per Raman transition.

Usage:
  julia --project=/home/sanghavi/code/github/vSmartMOM.jl \
        /home/sanghavi/code/github/vSmartMOM.jl/test/benchmarks/raman_batched_ops_benchmark.jl
=#

using CUDA
device!(1)
using NNlib
using LinearAlgebra
using Dates

println("="^70)
println("Batched GPU Operations Micro-Benchmark for Raman Optimization")
println("="^70)
println("GPU: ", CUDA.name(CUDA.device()))
println("Free GPU memory: ", round(CUDA.available_memory() / 1e9, digits=2), " GB")
println()

# ── Realistic dimensions ──
const NQ = 9       # NquadN (3 quad points × 3 Stokes)
const NSPEC = 5500  # spectral points per microwindow
const NRAMAN = 963  # Raman transitions
const NL = 5000     # typical length of n₁ range (≤ NSPEC)

println("Dimensions: NquadN=$NQ, nSpec=$NSPEC, n_Raman=$NRAMAN, nL=$NL")
println()

# ── Helper: Unicode batched_mul ──
⊠(A, B) = NNlib.batched_mul(A, B)

# =================================================================
# Test 1: Single expression with temporaries (current approach)
# =================================================================
println("-"^70)
println("Test 1: Single ⊠ expression (current approach - allocating)")
println("-"^70)

function test_allocating_expr(A_mat, B_mat, C_mat, D_vec, E_vec, F_vec; niter=100)
    for _ in 1:niter
        # Mimics: tmp = A ⊠ (D + B ⊠ E + C ⊠ F)
        #  This creates ~5 temporary GPU arrays per iteration
        result = A_mat ⊠ (D_vec .+ B_mat ⊠ E_vec .+ C_mat ⊠ F_vec)
    end
    CUDA.synchronize()
end

A_mat = CUDA.rand(Float64, NQ, NQ, NL)
B_mat = CUDA.rand(Float64, NQ, NQ, NL)
C_mat = CUDA.rand(Float64, NQ, NQ, NL)
D_vec = CUDA.rand(Float64, NQ, 1, NL)
E_vec = CUDA.rand(Float64, NQ, 1, NL)
F_vec = CUDA.rand(Float64, NQ, 1, NL)

# Warmup
test_allocating_expr(A_mat, B_mat, C_mat, D_vec, E_vec, F_vec; niter=2)
CUDA.reclaim()

t1 = @elapsed test_allocating_expr(A_mat, B_mat, C_mat, D_vec, E_vec, F_vec; niter=100)
println("  100 iterations: $(round(t1*1000, digits=1)) ms")
alloc1 = @allocated test_allocating_expr(A_mat, B_mat, C_mat, D_vec, E_vec, F_vec; niter=100)
println("  Allocations: $(round(alloc1/1e6, digits=1)) MB")

# =================================================================
# Test 2: In-place batched_mul! with pre-allocated buffers
# =================================================================
println()
println("-"^70)
println("Test 2: In-place batched_mul! (pre-allocated buffers)")
println("-"^70)

function test_inplace_expr(A_mat, B_mat, C_mat, D_vec, E_vec, F_vec,
                           ws_vec1, ws_vec2, ws_result; niter=100)
    for _ in 1:niter
        # Same computation: result = A ⊠ (D + B ⊠ E + C ⊠ F)
        # Step 1: ws_vec1 = D (copy)
        ws_vec1 .= D_vec
        # Step 2: ws_vec2 = B ⊠ E
        NNlib.batched_mul!(ws_vec2, B_mat, E_vec)
        ws_vec1 .+= ws_vec2
        # Step 3: ws_vec2 = C ⊠ F
        NNlib.batched_mul!(ws_vec2, C_mat, F_vec)
        ws_vec1 .+= ws_vec2
        # Step 4: result = A ⊠ ws_vec1
        NNlib.batched_mul!(ws_result, A_mat, ws_vec1)
    end
    CUDA.synchronize()
end

ws_vec1   = CUDA.zeros(Float64, NQ, 1, NL)
ws_vec2   = CUDA.zeros(Float64, NQ, 1, NL)
ws_result = CUDA.zeros(Float64, NQ, 1, NL)

# Warmup
test_inplace_expr(A_mat, B_mat, C_mat, D_vec, E_vec, F_vec, ws_vec1, ws_vec2, ws_result; niter=2)
CUDA.reclaim()

t2 = @elapsed test_inplace_expr(A_mat, B_mat, C_mat, D_vec, E_vec, F_vec, ws_vec1, ws_vec2, ws_result; niter=100)
println("  100 iterations: $(round(t2*1000, digits=1)) ms")
alloc2 = @allocated test_inplace_expr(A_mat, B_mat, C_mat, D_vec, E_vec, F_vec, ws_vec1, ws_vec2, ws_result; niter=100)
println("  Allocations: $(round(alloc2/1e6, digits=1)) MB")
println("  Speedup vs allocating: $(round(t1/t2, digits=2))x")

# =================================================================
# Test 3: 4D array views — batched_mul! on views of 4D CuArrays
# =================================================================
println()
println("-"^70)
println("Test 3: batched_mul! on views of 4D CuArrays")
println("-"^70)

function test_4d_views(A4d, B4d, C_vec4d, D_vec4d, ws_vec1, ws_vec2, ws_result;
                       niter=100, Δn_range=1:10)
    for _ in 1:niter
        for Δn in Δn_range
            # Mimic: result = A[:,:,n₁,Δn] ⊠ (C[:,:,n₁,Δn] + B[:,:,n₁,Δn])
            nL = size(A4d, 3)
            @views NNlib.batched_mul!(ws_vec1[:,:,1:nL], A4d[:,:,:,Δn], C_vec4d[:,:,:,Δn])
            @views ws_vec1[:,:,1:nL] .+= D_vec4d[:,:,:,Δn]
            @views NNlib.batched_mul!(ws_result[:,:,1:nL], A4d[:,:,:,Δn], ws_vec1[:,:,1:nL])
        end
    end
    CUDA.synchronize()
end

# Create 4D arrays mimicking inelastic R/T arrays
n_test_raman = 10  # test with small n_Raman to keep memory reasonable
A4d     = CUDA.rand(Float64, NQ, NQ, NL, n_test_raman)
B4d     = CUDA.rand(Float64, NQ, NQ, NL, n_test_raman)
C_vec4d = CUDA.rand(Float64, NQ, 1, NL, n_test_raman)
D_vec4d = CUDA.rand(Float64, NQ, 1, NL, n_test_raman)

# Warmup
test_4d_views(A4d, B4d, C_vec4d, D_vec4d, ws_vec1, ws_vec2, ws_result;
              niter=1, Δn_range=1:n_test_raman)
CUDA.reclaim()

t3 = @elapsed test_4d_views(A4d, B4d, C_vec4d, D_vec4d, ws_vec1, ws_vec2, ws_result;
                             niter=100, Δn_range=1:n_test_raman)
alloc3 = @allocated test_4d_views(A4d, B4d, C_vec4d, D_vec4d, ws_vec1, ws_vec2, ws_result;
                                   niter=100, Δn_range=1:n_test_raman)
println("  100 iterations × $(n_test_raman) Δn: $(round(t3*1000, digits=1)) ms")
println("  Allocations: $(round(alloc3/1e6, digits=1)) MB")

# Free 4D arrays
A4d = nothing; B4d = nothing; C_vec4d = nothing; D_vec4d = nothing
CUDA.reclaim()

# =================================================================
# Test 4: Loop of 2D mul! vs single 3D batched_mul!
# =================================================================
println()
println("-"^70)
println("Test 4: Loop of 2D mul! vs single 3D batched_mul!")
println("-"^70)

function test_loop_2d_mul(A3d, B3d, C3d; niter=100)
    nL = size(A3d, 3)
    for _ in 1:niter
        for k in 1:nL
            @views mul!(C3d[:,:,k], A3d[:,:,k], B3d[:,:,k])
        end
    end
    CUDA.synchronize()
end

function test_batched_3d(A3d, B3d, C3d; niter=100)
    for _ in 1:niter
        NNlib.batched_mul!(C3d, A3d, B3d)
    end
    CUDA.synchronize()
end

# mat × mat case
A3d_mm = CUDA.rand(Float64, NQ, NQ, NL)
B3d_mm = CUDA.rand(Float64, NQ, NQ, NL)
C3d_mm = CUDA.zeros(Float64, NQ, NQ, NL)

# Warmup
test_loop_2d_mul(A3d_mm, B3d_mm, C3d_mm; niter=1)
test_batched_3d(A3d_mm, B3d_mm, C3d_mm; niter=1)

t4a = @elapsed test_loop_2d_mul(A3d_mm, B3d_mm, C3d_mm; niter=100)
alloc4a = @allocated test_loop_2d_mul(A3d_mm, B3d_mm, C3d_mm; niter=100)
println("  2D mul! loop (mat×mat, 100 iter): $(round(t4a*1000, digits=1)) ms, alloc=$(round(alloc4a/1e6, digits=1)) MB")

t4b = @elapsed test_batched_3d(A3d_mm, B3d_mm, C3d_mm; niter=100)
alloc4b = @allocated test_batched_3d(A3d_mm, B3d_mm, C3d_mm; niter=100)
println("  batched_mul!  (mat×mat, 100 iter): $(round(t4b*1000, digits=1)) ms, alloc=$(round(alloc4b/1e6, digits=1)) MB")
println("  Batched speedup: $(round(t4a/t4b, digits=2))x")

# mat × vec case
A3d_mv = CUDA.rand(Float64, NQ, NQ, NL)
B3d_mv = CUDA.rand(Float64, NQ, 1, NL)
C3d_mv = CUDA.zeros(Float64, NQ, 1, NL)

test_loop_2d_mul(A3d_mv, B3d_mv, C3d_mv; niter=1)
test_batched_3d(A3d_mv, B3d_mv, C3d_mv; niter=1)

t4c = @elapsed test_loop_2d_mul(A3d_mv, B3d_mv, C3d_mv; niter=100)
t4d = @elapsed test_batched_3d(A3d_mv, B3d_mv, C3d_mv; niter=100)
println("  2D mul! loop (mat×vec, 100 iter): $(round(t4c*1000, digits=1)) ms")
println("  batched_mul!  (mat×vec, 100 iter): $(round(t4d*1000, digits=1)) ms")
println("  Batched speedup: $(round(t4c/t4d, digits=2))x")

# =================================================================
# Test 5: batch_inv! performance (as used in doubling)
# =================================================================
println()
println("-"^70)
println("Test 5: batch_inv! (from gpu_batched.jl)")
println("-"^70)

using CUDA.CUBLAS

function batch_inv_test!(X::CuArray{FT,3}, A::CuArray{FT,3}) where {FT}
    pivot, info = CUBLAS.getrf_strided_batched!(A, true); CUDA.synchronize()
    CUBLAS.getri_strided_batched!(A, X, pivot); CUDA.synchronize()
end

A_inv = CUDA.rand(Float64, NQ, NQ, NSPEC) .+ CuArray(Float64.(I(NQ)))  # ensure invertible
X_inv = similar(A_inv)

# Warmup
batch_inv_test!(X_inv, copy(A_inv))

t5 = @elapsed begin
    for _ in 1:10
        batch_inv_test!(X_inv, copy(A_inv))
    end
end
println("  batch_inv! ($(NQ)×$(NQ)×$(NSPEC), 10 iter): $(round(t5*1000, digits=1)) ms")

# =================================================================
# Test 6: Full doubling-like loop simulation
# =================================================================
println()
println("-"^70)
println("Test 6: Full doubling-like loop (963 Raman transitions)")
println("  Compares: allocating ⊠ expressions vs in-place batched_mul!")
println("-"^70)

# Simulate the first Raman SFI loop from doubling_inelastic.jl
# For each Δn: result = A[:,:,n₁] ⊠ (D[:,:,n₁,Δn] + B[:,:,n₁] ⊠ E[:,:,n₁,Δn])
# where n₁ is a range of length ~nL

function simulate_raman_loop_allocating(A3d, B3d, D4d, E4d, nRaman; nL=NL)
    for Δn = 1:nRaman
        @views result = A3d[:,:,1:nL] ⊠ (D4d[:,:,1:nL,Δn] .+ B3d[:,:,1:nL] ⊠ E4d[:,:,1:nL,Δn])
        @views D4d[:,:,1:nL,Δn] .= result
    end
    CUDA.synchronize()
end

function simulate_raman_loop_inplace(A3d, B3d, D4d, E4d, ws1, ws2, ws_out, nRaman; nL=NL)
    for Δn = 1:nRaman
        @views ws1[:,:,1:nL] .= D4d[:,:,1:nL,Δn]
        @views NNlib.batched_mul!(ws2[:,:,1:nL], B3d[:,:,1:nL], E4d[:,:,1:nL,Δn])
        @views ws1[:,:,1:nL] .+= ws2[:,:,1:nL]
        @views NNlib.batched_mul!(ws_out[:,:,1:nL], A3d[:,:,1:nL], ws1[:,:,1:nL])
        @views D4d[:,:,1:nL,Δn] .= ws_out[:,:,1:nL]
    end
    CUDA.synchronize()
end

# Use smaller n_Raman for memory reasons (scale results)
n_test = min(NRAMAN, 100)  # test with 100 transitions
nL_test = min(NL, 3000)    # smaller spectral range for memory

A3d_sim = CUDA.rand(Float64, NQ, NQ, NSPEC)
B3d_sim = CUDA.rand(Float64, NQ, NQ, NSPEC)
D4d_sim = CUDA.rand(Float64, NQ, 1, NSPEC, n_test)
E4d_sim = CUDA.rand(Float64, NQ, 1, NSPEC, n_test)
ws1_sim = CUDA.zeros(Float64, NQ, 1, NSPEC)
ws2_sim = CUDA.zeros(Float64, NQ, 1, NSPEC)
ws_out_sim = CUDA.zeros(Float64, NQ, 1, NSPEC)

# Warmup
D4d_bak = copy(D4d_sim)
simulate_raman_loop_allocating(A3d_sim, B3d_sim, D4d_sim, E4d_sim, 2; nL=nL_test)
D4d_sim .= D4d_bak
simulate_raman_loop_inplace(A3d_sim, B3d_sim, D4d_sim, E4d_sim, ws1_sim, ws2_sim, ws_out_sim, 2; nL=nL_test)
CUDA.reclaim()

# Correctness check
D4d_sim .= D4d_bak
D4d_check = copy(D4d_sim)

simulate_raman_loop_allocating(A3d_sim, B3d_sim, D4d_sim, E4d_sim, n_test; nL=nL_test)
result_alloc = Array(D4d_sim)

D4d_sim .= D4d_check
simulate_raman_loop_inplace(A3d_sim, B3d_sim, D4d_sim, E4d_sim, ws1_sim, ws2_sim, ws_out_sim, n_test; nL=nL_test)
result_inplace = Array(D4d_sim)

match = all(result_alloc .== result_inplace)
if !match
    maxdiff = maximum(abs.(result_alloc .- result_inplace))
    println("  WARNING: Results differ! Max abs diff = $maxdiff")
else
    println("  Correctness: BIT-EXACT MATCH between allocating and in-place")
end

# Timed runs
D4d_sim .= D4d_bak
CUDA.reclaim()
t6a = @elapsed simulate_raman_loop_allocating(A3d_sim, B3d_sim, D4d_sim, E4d_sim, n_test; nL=nL_test)
alloc6a = @allocated begin
    D4d_sim .= D4d_bak
    simulate_raman_loop_allocating(A3d_sim, B3d_sim, D4d_sim, E4d_sim, n_test; nL=nL_test)
end

D4d_sim .= D4d_bak
CUDA.reclaim()
t6b = @elapsed simulate_raman_loop_inplace(A3d_sim, B3d_sim, D4d_sim, E4d_sim, ws1_sim, ws2_sim, ws_out_sim, n_test; nL=nL_test)
alloc6b = @allocated begin
    D4d_sim .= D4d_bak
    simulate_raman_loop_inplace(A3d_sim, B3d_sim, D4d_sim, E4d_sim, ws1_sim, ws2_sim, ws_out_sim, n_test; nL=nL_test)
end

println("  Allocating ($(n_test) Δn): $(round(t6a*1000, digits=1)) ms, alloc=$(round(alloc6a/1e6, digits=1)) MB")
println("  In-place   ($(n_test) Δn): $(round(t6b*1000, digits=1)) ms, alloc=$(round(alloc6b/1e6, digits=1)) MB")
println("  Speedup: $(round(t6a/t6b, digits=2))x")
println("  Alloc reduction: $(round(alloc6a/max(alloc6b,1), digits=0))x")
println("  Extrapolated for $(NRAMAN) Δn: alloc ~$(round(t6a/n_test*NRAMAN*1000, digits=0)) ms vs ~$(round(t6b/n_test*NRAMAN*1000, digits=0)) ms")

# =================================================================
# Summary
# =================================================================
println()
println("="^70)
println("SUMMARY")
println("="^70)
println("  Approach                     | Time    | Allocs   | Notes")
println("  -----------------------------|---------|----------|------")
println("  Broadcasting ⊠ (allocating)  | $(lpad(round(t1*1000, digits=1), 6)) ms | $(lpad(round(alloc1/1e6, digits=1), 7)) MB | Current code")
println("  batched_mul! (in-place)      | $(lpad(round(t2*1000, digits=1), 6)) ms | $(lpad(round(alloc2/1e6, digits=1), 7)) MB | Proposed")
println("  2D mul! loop (mat×mat)       | $(lpad(round(t4a*1000, digits=1), 6)) ms | $(lpad(round(alloc4a/1e6, digits=1), 7)) MB | Alternative")
println("  batched_mul! (mat×mat)       | $(lpad(round(t4b*1000, digits=1), 6)) ms | $(lpad(round(alloc4b/1e6, digits=1), 7)) MB | Alternative")
println("  Raman loop alloc ($(n_test) Δn)   | $(lpad(round(t6a*1000, digits=1), 6)) ms | $(lpad(round(alloc6a/1e6, digits=1), 7)) MB | Realistic sim")
println("  Raman loop inplace ($(n_test) Δn) | $(lpad(round(t6b*1000, digits=1), 6)) ms | $(lpad(round(alloc6b/1e6, digits=1), 7)) MB | Realistic sim")
println()
println("Key findings for optimization strategy:")
println("  - batched_mul! speedup over ⊠: $(round(t1/t2, digits=2))x")
println("  - 3D batched vs 2D loop (mat×mat): $(round(t4a/t4b, digits=2))x")
println("  - 3D batched vs 2D loop (mat×vec): $(round(t4c/t4d, digits=2))x")
println("  - Raman loop speedup: $(round(t6a/t6b, digits=2))x")
println("="^70)
