# Phase-2 fused interaction kernels vs the Phase-0 _bmm! path (_11 case).
# Contract: tolerance-equivalent (reduction order differs from cuBLAS).
using Test
using vSmartMOM
using vSmartMOM.CoreRT: ka_fused_interaction_down!, ka_fused_interaction_up!,
                        interaction_helper!, ScatteringInterface_11,
                        AddedLayer, CompositeLayer,
                        _FUSED_INTERACTION_ENABLED, _FUSED_GP_ENABLED
using KernelAbstractions
using LinearAlgebra
using Random
Random.seed!(23)

const N, S = 24, 129

function mk_layers(FT, to_dev)
    mk() = to_dev(FT(0.05) .* rand(FT, N, N, S))
    mkv() = to_dev(rand(FT, N, 1, S))
    a = AddedLayer(
        r⁻⁺ = mk(), t⁺⁺ = mk() .+ FT(0.5), r⁺⁻ = mk(), t⁻⁻ = mk() .+ FT(0.5),
        j₀⁺ = mkv(), j₀⁻ = mkv(),
        temp1 = to_dev(zeros(FT, N, N, S)), temp2 = to_dev(zeros(FT, N, N, S)),
        temp1_ptr = nothing, temp2_ptr = nothing,
        dbl_gp_refl = to_dev(zeros(FT, N, N, S)),
        dbl_j₁⁺ = to_dev(zeros(FT, N, 1, S)), dbl_j₁⁻ = to_dev(zeros(FT, N, 1, S)),
        dbl_v1 = to_dev(zeros(FT, N, 1, S)), dbl_v2 = to_dev(zeros(FT, N, 1, S)))
    c = CompositeLayer(
        R⁻⁺ = mk(), R⁺⁻ = mk(), T⁺⁺ = mk() .+ FT(0.5), T⁻⁻ = mk() .+ FT(0.5),
        J₀⁺ = mkv(), J₀⁻ = mkv())
    return a, c
end
copy_comp(c) = CompositeLayer(R⁻⁺ = copy(c.R⁻⁺), R⁺⁻ = copy(c.R⁺⁻),
                              T⁺⁺ = copy(c.T⁺⁺), T⁻⁻ = copy(c.T⁻⁻),
                              J₀⁺ = copy(c.J₀⁺), J₀⁻ = copy(c.J₀⁻))

function run_case(FT, to_dev, I_static; rtol)
    a, c = mk_layers(FT, to_dev)
    c_ref = copy_comp(c)
    # reference: kill-switch forces the Phase-0 _bmm! interaction path
    _FUSED_INTERACTION_ENABLED[] = false
    interaction_helper!(ScatteringInterface_11(), true, c_ref, a, I_static)
    _FUSED_INTERACTION_ENABLED[] = true
    interaction_helper!(ScatteringInterface_11(), true, c, a, I_static)
    for f in (:R⁻⁺, :R⁺⁻, :T⁺⁺, :T⁻⁻, :J₀⁺, :J₀⁻)
        @test isapprox(Array(getproperty(c, f)), Array(getproperty(c_ref, f)); rtol = rtol)
    end
end

@testset "fused interaction _11 vs _bmm! (CPU — gate forces fallback, sanity)" begin
    # On CPU the gate is off for BOTH runs; this asserts determinism of the
    # reference path (guards the test harness itself).
    I_static = repeat(Matrix{Float64}(I, N, N), 1, 1, 1)
    a, c = mk_layers(Float64, identity)
    c2 = copy_comp(c)
    interaction_helper!(ScatteringInterface_11(), true, c, a, I_static)
    interaction_helper!(ScatteringInterface_11(), true, c2, a, I_static)
    @test c.R⁻⁺ == c2.R⁻⁺ && c.J₀⁺ == c2.J₀⁺
end

using CUDA
if CUDA.functional()
    CUDA.device!(parse(Int, get(ENV, "BENCH_DEV", "1")))
    @testset "fused interaction _11 vs _bmm! (CUDA)" begin
        I64 = CuArray(repeat(Matrix{Float64}(I, N, N), 1, 1, 1))
        I32 = CuArray(repeat(Matrix{Float32}(I, N, N), 1, 1, 1))
        run_case(Float64, CuArray, I64; rtol = 1e-12)
        run_case(Float32, CuArray, I32; rtol = 3e-5)
    end
else
    @info "CUDA unavailable — GPU arm skipped"
end
println("fused interaction tests: done")
