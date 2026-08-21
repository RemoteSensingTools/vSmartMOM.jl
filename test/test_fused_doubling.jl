# Phase-1 fused doubling kernels vs the Phase-0 _bmm! path.
# Contract: tolerance-equivalent (NOT bitwise — different reduction order
# than cuBLAS/BLAS), on KA-CPU (explicit call; the production gate is
# GPU-only) and on CUDA when available, F64 and F32.
using Test
using vSmartMOM
using vSmartMOM.CoreRT: doubling_source_update!, doubling_rt_update!,
                        ka_fused_doubling_rt!, ka_fused_doubling_source!,
                        _FUSED_DOUBLING_ENABLED
using KernelAbstractions
using Random
Random.seed!(11)

function bmm_reference!(r, t, tt, jp, jm, j1p, j1m, expk)
    v1, v2 = similar(jp), similar(jp)
    m1, m2 = similar(r), similar(r)
    was = _FUSED_DOUBLING_ENABLED[]
    _FUSED_DOUBLING_ENABLED[] = false
    try
        doubling_source_update!(jp, jm, j1p, j1m, r, tt, expk, v1, v2)
        doubling_rt_update!(r, t, tt, expk, m1, m2)
    finally
        _FUSED_DOUBLING_ENABLED[] = was
    end
end

function fused!(r, t, tt, jp, jm, j1p, j1m, expk, backend)
    ka_fused_doubling_source!(jp, jm, j1p, j1m, r, tt, expk, backend)
    ka_fused_doubling_rt!(r, t, tt, expk, backend)
end

function run_case(FT, to_dev, backend; rtol)
    N, S = 24, 129
    mk() = to_dev(FT(0.05) .* rand(FT, N, N, S))
    mkv() = to_dev(rand(FT, N, 1, S))
    r, t, tt = mk(), mk(), mk()
    jp, jm, j1p, j1m = mkv(), mkv(), mkv(), mkv()
    expk = to_dev(rand(FT, S))

    Rr, Rt, Rtt = copy(r), copy(t), copy(tt)
    Rjp, Rjm, Rj1p, Rj1m = copy(jp), copy(jm), copy(j1p), copy(j1m)
    Rexpk = copy(expk)
    bmm_reference!(Rr, Rt, Rtt, Rjp, Rjm, Rj1p, Rj1m, Rexpk)

    fused!(r, t, tt, jp, jm, j1p, j1m, expk, backend)
    KernelAbstractions.synchronize(backend)

    for (got, want, name) in ((r, Rr, "r⁻⁺"), (t, Rt, "t⁺⁺"),
                              (jp, Rjp, "j₀⁺"), (jm, Rjm, "j₀⁻"),
                              (j1p, Rj1p, "j₁⁺"), (j1m, Rj1m, "j₁⁻"),
                              (expk, Rexpk, "expk"))
        @test isapprox(Array(got), Array(want); rtol = rtol)
    end
end

@testset "fused doubling kernels vs _bmm! (KA-CPU)" begin
    run_case(Float64, identity, KernelAbstractions.CPU(); rtol = 1e-13)
    run_case(Float32, identity, KernelAbstractions.CPU(); rtol = 2e-5)
end

using CUDA
if CUDA.functional()
    CUDA.device!(parse(Int, get(ENV, "BENCH_DEV", "1")))
    @testset "fused doubling kernels vs _bmm! (CUDA)" begin
        run_case(Float64, CuArray, CUDA.CUDABackend(); rtol = 1e-13)
        run_case(Float32, CuArray, CUDA.CUDABackend(); rtol = 2e-5)
    end
    @testset "wire-in: gate uses the fused path on GPU" begin
        using vSmartMOM.CoreRT: _use_fused_doubling_rt, _use_fused_doubling_source
        x = CUDA.rand(Float32, 24, 24, 8)
        @test _use_fused_doubling_rt(x)
        @test _use_fused_doubling_source(x)
        _FUSED_DOUBLING_ENABLED[] = false
        @test !_use_fused_doubling_rt(x)
        _FUSED_DOUBLING_ENABLED[] = true
    end
else
    @info "CUDA unavailable — GPU arm skipped"
end
println("fused doubling kernel tests: done")
