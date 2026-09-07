using vSmartMOM, CUDA, KernelAbstractions, Statistics
const C = vSmartMOM.CoreRT
# Standalone historical probe: compare the per-element kernel with alternatives.
C._TILED_JACOBIANS_ENABLED[] = false
@kernel function tiled!(out,@Const(A),@Const(dA),@Const(B),@Const(dB),::Val{N}) where N
    batch=@index(Group,Linear)
    tid=@index(Local,Linear)
    s=mod1(batch,size(out,3)); p=(batch-1)÷size(out,3)+1
    i=mod1(tid,N); j=(tid-1)÷N+1
    a=@localmem eltype(out) (N,N)
    da=@localmem eltype(out) (N,N)
    b=@localmem eltype(out) (N,N)
    db=@localmem eltype(out) (N,N)
    @inbounds begin
        a[i,j]=A[i,j,s]; da[i,j]=dA[i,j,s,p]
        b[i,j]=B[i,j,s]; db[i,j]=dB[i,j,s,p]
    end
    @synchronize
    x=zero(eltype(out)); y=zero(eltype(out))
    @inbounds for k in 1:N
        x+=da[i,k]*b[k,j]; y+=a[i,k]*db[k,j]
    end
    @inbounds out[i,j,s,p]=x+y
end
@kernel function flat!(out,@Const(A),@Const(dA),@Const(B),@Const(dB),::Val{N}) where N
    idx=@index(Global,Linear)
    i=mod1(idx,N); j=mod((idx-1)÷N,N)+1; batch=(idx-1)÷(N*N)+1
    s=mod1(batch,size(out,3)); p=(batch-1)÷size(out,3)+1
    x=zero(eltype(out)); y=zero(eltype(out))
    @inbounds for k in 1:N
        x+=dA[i,k,s,p]*B[k,j,s]; y+=A[i,k,s]*dB[k,j,s,p]
    end
    @inbounds out[i,j,s,p]=x+y
end
function bench(f)
    f();CUDA.synchronize()
    median([CUDA.@elapsed(f()) for _ in 1:5])
end
CUDA.allowscalar(false)
for n in (6,18,32), ns in (512,10000)
    A=CUDA.rand(Float64,n,n,ns);B=CUDA.rand(Float64,n,n,ns)
    dA=CUDA.rand(Float64,n,n,ns,14);dB=CUDA.rand(Float64,n,n,ns,14)
    out=similar(dB);ref=similar(dB)
    old=bench(()->C._jprod!(ref,A,dA,B,dB))
    tile=bench(()->tiled!(get_backend(out),n*n)(out,A,dA,B,dB,Val(n);ndrange=n*n*ns*14))
    @assert isapprox(out,ref;rtol=1e-12,atol=1e-12)
    linear=bench(()->flat!(get_backend(out),256)(out,A,dA,B,dB,Val(n);ndrange=length(out)))
    @assert isapprox(out,ref;rtol=1e-12,atol=1e-12)
    println("TILE n=$n ns=$ns old=$old tile=$tile flat=$linear");flush(stdout)
end
