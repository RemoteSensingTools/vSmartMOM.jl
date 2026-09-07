# cuBLAS pointer batches for medium RT operators. Each plan folds (λ,p) into
# one batch, while 3D forward operands repeat their λ pointers across p.
# Pointer vectors are generated on the device and belong to one solve's fixed
# work arrays. There is no global cache or host list of per-column matrices.
struct CUDAJacobianProducts{P}
    pointers::Dict{NTuple{5,UInt},Tuple{Any,P}}
end

function CoreRT.make_jacobian_products(A::CuArray{FT,3,CUDA.DeviceMemory}, np) where {FT}
    CoreRT._VENDOR_JACOBIANS_ENABLED[] && FT <: Union{Float32,Float64} &&
        32 < size(A,1) <= 64 && size(A,3) >= 512 || return nothing
    P = typeof(CuArray{CuPtr{FT}}(undef,0))
    CUDAJacobianProducts(Dict{NTuple{5,UInt},Tuple{Any,P}}())
end

@kernel function _jacobian_batch_pointers!(out, base, matrix_bytes, ns, shared)
    b = @index(Global,Linear)
    s = shared ? mod1(b,ns) : b
    # Passing the address as UInt prevents CUDA adaptation from changing a
    # host CuPtr into an LLVMPtr of a different type than this pointer vector.
    @inbounds out[b] = reinterpret(eltype(out),base+UInt((s-1)*matrix_bytes))
end

function _jacobian_pointers!(cache::CUDAJacobianProducts,A,ns,np)
    # Active-prefix views share the same address and strides between layers.
    # Key by that geometry, not view-object identity, so the cache stays
    # bounded when a fresh equivalent view is constructed at each layer.
    bytes = sizeof(eltype(A))*stride(A,3)
    key = (UInt(pointer(A)),UInt(bytes),UInt(ns),UInt(np),UInt(ndims(A)))
    entry = get!(cache.pointers,key) do
        ptrs = CuArray{CuPtr{eltype(A)}}(undef,ns*np)
        _jacobian_batch_pointers!(KernelAbstractions.get_backend(ptrs))(
            ptrs,key[1],bytes,ns,ndims(A)==3;ndrange=length(ptrs))
        (A,ptrs) # Retain the backing allocation for as long as its plan lives.
    end
    return entry[2]
end

# Source-vector products keep the existing fused kernel. Medium square
# products were benchmarked at N=33/57/64; N≤32 retains its faster KA tiles.
@inline _vendor_batch_layout(A,ns,np) =
    stride(A,1)==1 && size(A,3)==ns && (ndims(A)==3 ||
        (ndims(A)==4 && size(A,4)>=np && stride(A,4)==ns*stride(A,3)))

@inline function _vendor_jacobian_square(out,A,B)
    n,ns,np = size(out,1),size(out,3),size(out,4)
    return !isempty(out) && 32<n<=64 && ns>=512 &&
        n==size(out,2)==size(A,1)==size(A,2)==size(B,1)==size(B,2) &&
        stride(out,2)==stride(A,2)==stride(B,2)==n &&
        _vendor_batch_layout(out,ns,np) && _vendor_batch_layout(A,ns,np) &&
        _vendor_batch_layout(B,ns,np)
end

for (FT,gemm) in ((Float64,:cublasDgemmBatched),(Float32,:cublasSgemmBatched))
    @eval function _vendor_jacobian_mul!(cache,out::AbstractArray{$FT},A,B,β)
        n,ns,np = size(out,1),size(out,3),size(out,4)
        Ap = _jacobian_pointers!(cache,A,ns,np)
        Bp = _jacobian_pointers!(cache,B,ns,np)
        Cp = _jacobian_pointers!(cache,out,ns,np)
        CUDA.CUBLAS.$gemm(CUDA.CUBLAS.handle(),'N','N',n,n,n,
            one($FT),Ap,n,Bp,n,$FT(β),Cp,n,ns*np)
        return out
    end
end

function CoreRT._jmul_with_plan!(cache::CUDAJacobianProducts,out,A,B,β)
    _vendor_jacobian_square(out,A,B) || return CoreRT._jmul!(out,A,B,β)
    return _vendor_jacobian_mul!(cache,out,A,B,β)
end

function CoreRT._jprod_with_plan!(cache::CUDAJacobianProducts,out,A,dA,B,dB)
    (_vendor_jacobian_square(out,dA,B) && _vendor_jacobian_square(out,A,dB)) ||
        return CoreRT._jprod!(out,A,dA,B,dB)
    # S2014 (C.6): d(AB)=dA B+A dB. β=0 then β=1 avoids an intermediate
    # tensor, and both calls stay on the same CUDA stream. No inverse changes.
    _vendor_jacobian_mul!(cache,out,dA,B,zero(eltype(out)))
    return _vendor_jacobian_mul!(cache,out,A,dB,one(eltype(out)))
end

# A copied workspace owns new GPU allocations. Rebuild its pointer metadata
# lazily instead of copying device addresses that refer to the original arrays.
function Base.deepcopy_internal(cache::CUDAJacobianProducts, stack::IdDict)
    haskey(stack,cache) && return stack[cache]
    copy = CUDAJacobianProducts(empty(cache.pointers))
    stack[cache] = copy
    return copy
end
