# Compare blocked product rules with one cuBLAS pointer batch across (λ,p).
# Experimental only: a production path would need workspace-owned pointer plans.
using vSmartMOM, CUDA, KernelAbstractions, LinearAlgebra, Statistics, Random
const C= vSmartMOM.CoreRT
CUDA.allowscalar(false)
BLAS.set_num_threads(1)

@kernel function spectral_batch_pointers!(ptrs, base, matrix_bytes, ns, shared)
    b=@index(Global,Linear)
    s=shared ? mod1(b,ns) : b
    @inbounds ptrs[b]=reinterpret(eltype(ptrs),base+UInt((s-1)*matrix_bytes))
end
function batch_pointers(A,ns,np)
    ptrs=CuArray{CuPtr{eltype(A)}}(undef,ns*np)
    spectral_batch_pointers!(KernelAbstractions.get_backend(ptrs))(ptrs,UInt(pointer(A)),
        sizeof(eltype(A))*size(A,1)*size(A,2),ns,ndims(A)==3;ndrange=length(ptrs))
    return ptrs
end
function vendor_plan(out,A,dA,B,dB)
    ns,np=size(out,3),size(out,4)
    # Retain owning arrays with their pointer vectors; no global pointer cache.
    (;out,A,dA,B,dB,pointers=map(x->batch_pointers(x,ns,np),(out,A,dA,B,dB)))
end
for (FT,gemm) in ((Float64,:cublasDgemmBatched),(Float32,:cublasSgemmBatched))
    @eval function vendor_product!(out::CuArray{$FT,4},plan)
        (;A,dA,B,dB,pointers)=plan
        Cp,Ap,dAp,Bp,dBp=pointers
        n=size(out,1); batches=size(out,3)*size(out,4)
        CUDA.CUBLAS.$gemm(CUDA.CUBLAS.handle(),'N','N',n,n,n,
            one($FT),dAp,n,Bp,n,zero($FT),Cp,n,batches)
        CUDA.CUBLAS.$gemm(CUDA.CUBLAS.handle(),'N','N',n,n,n,
            one($FT),Ap,n,dBp,n,one($FT),Cp,n,batches)
        return out
    end
end
function measured(f)
    f();CUDA.synchronize()
    samples=map(1:3) do _
        GC.gc();CUDA.synchronize()
        (@elapsed begin for _ in 1:20;f();end;CUDA.synchronize();end)/20
    end
    (;median=median(samples),samples)
end
function main()
    rng=Xoshiro(608)
    for FT in (Float64,Float32), (n,ns,np) in ((18,10000,3),(33,512,7),(57,512,7),(64,512,7))
        A,B=map(_->CuArray(randn(rng,FT,n,n,ns)),1:2)
        dA,dB=map(_->CuArray(randn(rng,FT,n,n,ns,np)),1:2)
        output=similar(dA)
        plan=vendor_plan(output,A,dA,B,dB)
        blocked()=C._jprod!(output,A,dA,B,dB)
        vendor()=vendor_product!(output,plan)
        blocked();reference=Array(output)
        vendor();result=Array(output)
        @assert isapprox(result,reference;rtol=200eps(FT),atol=200eps(FT))
        bt=measured(blocked);vt=measured(vendor)
        println("VENDOR_PRODUCT FT=$FT n=$n ns=$ns directions=$np blocked=$bt vendor=$vt speedup=$(bt.median/vt.median) max_error=$(maximum(abs,result-reference))")
        flush(stdout)
    end
end
main()
