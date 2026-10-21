using vSmartMOM, CUDA, NNlib, Statistics
const C = vSmartMOM.CoreRT
# Standalone historical probe: compare the per-element kernel with alternatives.
C._TILED_JACOBIANS_ENABLED[] = false
function vendor!(out,A,dA,B,dB)
    for p in axes(out,4)
        NNlib.batched_mul!(@view(out[:,:,:,p]), @view(dA[:,:,:,p]), B)
        NNlib.batched_mul!(@view(out[:,:,:,p]), A, @view(dB[:,:,:,p]), 1.0, 1.0)
    end
end
function bench(f)
    f(); CUDA.synchronize()
    median([CUDA.@elapsed(f()) for _ in 1:5])
end
CUDA.allowscalar(false)
for n in (6,18), ns in (64,512,10000), cols in (1,n)
    A=CUDA.rand(Float64,n,n,ns); B=CUDA.rand(Float64,n,cols,ns)
    dA=CUDA.rand(Float64,n,n,ns,14); dB=CUDA.rand(Float64,n,cols,ns,14)
    out=similar(dB); ref=similar(dB)
    fast=bench(()->C._jprod!(out,A,dA,B,dB))
    vendor=bench(()->vendor!(ref,A,dA,B,dB))
    @assert isapprox(out,ref;rtol=1e-12,atol=1e-12)
    println("PRODUCT n=$n ns=$ns cols=$cols portable=$fast vendor=$vendor ratio=$(fast/vendor)");flush(stdout)
end
