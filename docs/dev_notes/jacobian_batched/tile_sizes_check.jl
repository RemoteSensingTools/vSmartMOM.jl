# Measure the small-operator tile crossover, including external-solar IQU N=15.
using vSmartMOM, CUDA, KernelAbstractions, Statistics, Random
const C=vSmartMOM.CoreRT
CUDA.allowscalar(false)
function measure_tiles(f)
    f();CUDA.synchronize()
    samples=map(1:3) do _
        GC.gc();CUDA.synchronize()
        (@elapsed begin for _ in 1:20;f();end;CUDA.synchronize();end)/20
    end
    median(samples)
end
function main()
    rng=Xoshiro(1518)
    for FT in (Float64,Float32), n in (6,9,12,15,18)
        ns,np=parse(Int,get(ENV,"AUDIT_NSPEC","10000")),3
        A,B=map(_->CuArray(randn(rng,FT,n,n,ns)),1:2)
        dA,dB=map(_->CuArray(randn(rng,FT,n,n,ns,np)),1:2)
        out=similar(dA);backend=KernelAbstractions.get_backend(out)
        plain()=C._jac_product_kernel!(backend)(out,A,dA,B,dB,Val(n);ndrange=(n,n,ns*np))
        tiled()=C._jac_product_tiled!(backend,n*n)(out,A,dA,B,dB,Val(n);ndrange=n*n*ns*np)
        plain();reference=Array(out)
        tiled();result=Array(out)
        @assert isapprox(result,reference;rtol=200eps(FT),atol=200eps(FT))
        pt,tt=measure_tiles(plain),measure_tiles(tiled)
        println("TILE_CROSSOVER FT=$FT n=$n ns=$ns directions=$np plain=$pt tiled=$tt speedup=$(pt/tt) max_error=$(maximum(abs,result-reference))")
        flush(stdout)
    end
end
main()
