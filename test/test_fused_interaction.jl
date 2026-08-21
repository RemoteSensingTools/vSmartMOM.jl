# Phase-2 fused interaction kernels vs the Phase-0 _bmm! path (_11 case).
# Contract: tolerance-equivalent (reduction order differs from cuBLAS).
using Test
using vSmartMOM
using vSmartMOM.CoreRT: ka_fused_interaction_down!, ka_fused_interaction_up!,
                        interaction_helper!, _interaction_11_fused!,
                        _use_fused_interaction, _fused_interaction_default,
                        ScatteringInterface_11, AddedLayer, CompositeLayer,
                        _FUSED_INTERACTION_MODE, _FUSED_GP_ENABLED
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
    c_dir = copy_comp(c)
    # Mode policy: default :auto defers to the backend — the CUDA ext opts
    # OUT (cuBLAS wins at production sizes), so :auto must gate false here.
    # :on/:off must flip the gate both ways on this device — the guard that
    # was missing when the kernels were left unwired: the gate result must
    # actually track the mode, not silently stay on one path.
    default = _FUSED_INTERACTION_MODE[]
    @test default === :auto
    @test _fused_interaction_default(KernelAbstractions.get_backend(a.r⁻⁺)) == false  # CUDA ext override
    @test !_use_fused_interaction(a.r⁻⁺)                 # :auto → CUDA backend policy → off
    _FUSED_INTERACTION_MODE[] = :off
    @test !_use_fused_interaction(a.r⁻⁺)
    interaction_helper!(ScatteringInterface_11(), true, c_ref, a, I_static)  # _bmm! ladder
    _FUSED_INTERACTION_MODE[] = :on
    @test _use_fused_interaction(a.r⁻⁺)
    interaction_helper!(ScatteringInterface_11(), true, c, a, I_static)      # fused branch
    # the dedicated fused function, called directly, must agree with the
    # gated dispatch (same inputs → same launches)
    _interaction_11_fused!(c_dir, a, I_static)
    _FUSED_INTERACTION_MODE[] = default
    for f in (:R⁻⁺, :R⁺⁻, :T⁺⁺, :T⁻⁻, :J₀⁺, :J₀⁻)
        @test isapprox(Array(getproperty(c, f)), Array(getproperty(c_ref, f)); rtol = rtol)
        @test Array(getproperty(c_dir, f)) == Array(getproperty(c, f))
    end
end

@testset "fused kernels on the KA-CPU backend vs plain matrix algebra" begin
    # Direct kernel calls on KernelAbstractions.CPU() — GPU-less CI compiles
    # and validates the retained kernels (the gated production path never
    # reaches them on CPU, so without this they'd be dead code in CI).
    # Reference = the adding equations written as plain matrix algebra:
    #   DOWN: J₀⁻ += G₁(r⁻⁺J₀⁺ + j₀⁻);  R⁻⁺ += G₁r⁻⁺T⁺⁺;  T⁻⁻ = G₁t⁻⁻
    #   UP:   J₀⁺ = j₀⁺ + G₂(J₀⁺ + R⁺⁻j₀⁻);  T⁺⁺ = G₂T⁺⁺;  R⁺⁻ = r⁺⁻ + G₂R⁺⁻t⁻⁻
    n, s = 6, 4
    cpu = KernelAbstractions.CPU()
    mk()  = 0.1 .* rand(n, n, s)
    mkv() = rand(n, 1, s)
    G₁, G₂ = mk(), mk()
    r⁻⁺, r⁺⁻, t⁻⁻, T⁺⁺, R⁻⁺, R⁺⁻ = mk(), mk(), mk(), mk(), mk(), mk()
    J₀⁺, J₀⁻, j₀⁺, j₀⁻ = mkv(), mkv(), mkv(), mkv()

    # DOWN half
    J₀⁻d, R⁻⁺d, T⁻⁻d = copy(J₀⁻), copy(R⁻⁺), copy(t⁻⁻ .* 0)
    T⁻⁻d .= mk()   # arbitrary pre-state; kernel overwrites it
    ka_fused_interaction_down!(J₀⁻d, R⁻⁺d, T⁻⁻d, G₁, r⁻⁺, J₀⁺, j₀⁻, T⁺⁺, t⁻⁻, cpu)
    for k in 1:s
        w = r⁻⁺[:, :, k] * J₀⁺[:, 1, k] .+ j₀⁻[:, 1, k]
        @test J₀⁻d[:, 1, k] ≈ J₀⁻[:, 1, k] .+ G₁[:, :, k] * w            rtol = 1e-12
        @test R⁻⁺d[:, :, k] ≈ R⁻⁺[:, :, k] .+ G₁[:, :, k] * r⁻⁺[:, :, k] * T⁺⁺[:, :, k] rtol = 1e-12
        @test T⁻⁻d[:, :, k] ≈ G₁[:, :, k] * t⁻⁻[:, :, k]                 rtol = 1e-12
    end

    # UP half (kernel overwrites J₀⁺/T⁺⁺/R⁺⁻ — feed copies, keep originals as pre-state)
    J₀⁺u, T⁺⁺u, R⁺⁻u = copy(J₀⁺), copy(T⁺⁺), copy(R⁺⁻)
    ka_fused_interaction_up!(J₀⁺u, T⁺⁺u, R⁺⁻u, G₂, j₀⁺, j₀⁻, r⁺⁻, t⁻⁻, cpu)
    for k in 1:s
        w = J₀⁺[:, 1, k] .+ R⁺⁻[:, :, k] * j₀⁻[:, 1, k]
        @test J₀⁺u[:, 1, k] ≈ j₀⁺[:, 1, k] .+ G₂[:, :, k] * w            rtol = 1e-12
        @test T⁺⁺u[:, :, k] ≈ G₂[:, :, k] * T⁺⁺[:, :, k]                 rtol = 1e-12
        @test R⁺⁻u[:, :, k] ≈ r⁺⁻[:, :, k] .+ G₂[:, :, k] * R⁺⁻[:, :, k] * t⁻⁻[:, :, k] rtol = 1e-12
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
