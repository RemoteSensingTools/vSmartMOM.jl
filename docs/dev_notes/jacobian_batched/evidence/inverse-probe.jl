using vSmartMOM, CUDA, LinearAlgebra, Statistics
const C=vSmartMOM.CoreRT
CUDA.allowscalar(false)
function measure(f)
    f();CUDA.synchronize()
    median([CUDA.@elapsed(f()) for _ in 1:7])
end
for n in (6,18,32)
    ns=10000
    A=CUDA.rand(Float64,n,n,ns).*0.001
    E=CuArray(repeat(Matrix{Float64}(I,n,n),1,1,ns))
    M=E.-A
    temp=similar(M); out=similar(M); ref=similar(M)
    vendor=measure(()->begin copyto!(temp,M);C.batch_inv!(out,temp) end)
    ka=measure(()->C.ka_batch_inv_lu!(ref,M,CUDA.CUDABackend()))
    @assert isapprox(out,ref;rtol=1e-12,atol=1e-12)
    fused=measure(()->C.ka_fused_solve!(ref,A,A,E,CUDA.CUDABackend()))
    println("INV n=$n vendor=$vendor ka=$ka fused_product_inverse=$fused");flush(stdout)
end
