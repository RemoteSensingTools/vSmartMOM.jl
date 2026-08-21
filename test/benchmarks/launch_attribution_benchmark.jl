# Kernel-launch attribution for the vSmartMOM adding-doubling hot lines.
# Production shapes at nstr=8/IQU: supermatrices (N,N,nspec), N = 24-ish.
# Counts kernels + host time per expression via CUDA.@profile.
# Run from the repo root:  julia --project=. test/benchmarks/launch_attribution_benchmark.jl
using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."); io = devnull)
using CUDA, vSmartMOM, Printf
const bmm = vSmartMOM.CoreRT.batched_mul   # NNlib ⊠

CUDA.device!(parse(Int, get(ENV, "BENCH_DEV", "0")))
N, S = 24, 801
FT = Float32
A()  = CUDA.rand(FT, N, N, S)
V()  = CUDA.rand(FT, N, 1, S)

r, t, tt, tmp = A(), A(), A(), A()
j0p, j0m, j1p, j1m = V(), V(), V(), V()
expk = CUDA.rand(FT, S)
Istat = CuArray(repeat(reshape(Matrix{FT}(vSmartMOM.CoreRT.LinearAlgebra.I, N, N), N, N, 1), 1, 1, S))

function prof(name, f)
    f(); CUDA.synchronize()                     # warm/compile
    p = CUDA.@profile begin
        for _ in 1:10; f(); end
        CUDA.synchronize()
    end
    nk = length(p.device.id) / 10
    tdev = sum(p.device.stop .- p.device.start) / 10 * 1e3   # ms per iter
    thost = sum(p.host.stop .- p.host.start) / 10 * 1e3
    @printf("%-46s : %5.1f kernels/iter | GPU %7.3f ms | host(API) %7.3f ms\n",
            name, nk, tdev, thost)
end

println("shapes: ($N,$N,$S) matrices, ($N,1,$S) sources, Float32\n")

prof("pure dotted chain  a .= a .+ b .* c",       () -> (tmp .= tmp .+ r .* t))
prof("single ⊠           tmp .= r ⊠ t",           () -> (tmp .= bmm(r, t)))
prof("rt_update line     r .= r .+ tt ⊠ r ⊠ t",   () -> (r .= r .+ bmm(bmm(tt, r), t)))
prof("source line        j₀⁻.=j₀⁻.+tt⊠(j₁⁻.+r⊠j₀⁺)",
     () -> (j0m .= j0m .+ bmm(tt, (j1m .+ bmm(r, j0p)))))
prof("expk view scaling  j₁[:,1,:] .= j₀[:,1,:].*k'",
     () -> (@views j1p[:, 1, :] .= j0p[:, 1, :] .* expk'))
prof("geometric progression (fused-or-fallback)",
     () -> vSmartMOM.CoreRT.compute_geometric_progression!(tmp, tt, r, t, Istat, A(), nothing, nothing))
prof("Bool reduce        any(isnan, r)",          () -> any(isnan, r))
prof("full doubling_source_update! (allocating)",
     () -> vSmartMOM.CoreRT.doubling_source_update!(j0p, j0m, j1p, j1m, r, tt, expk))
prof("full doubling_rt_update!     (allocating)",
     () -> vSmartMOM.CoreRT.doubling_rt_update!(r, t, tt, expk))
# Phase-0 in-place variants (preallocated temps)
let sv1 = V(), sv2 = V(), sm1 = A(), sm2 = A()
    prof("doubling_source_update! (in-place, Phase 0)",
         () -> vSmartMOM.CoreRT.doubling_source_update!(j0p, j0m, j1p, j1m, r, tt, expk, sv1, sv2))
    prof("doubling_rt_update!     (in-place, Phase 0)",
         () -> vSmartMOM.CoreRT.doubling_rt_update!(r, t, tt, expk, sm1, sm2))
end
