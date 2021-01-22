# batch Matrix inversion for GPU
function batch_solve!(X::CuArray{FT,3}, A::CuArray{FT,3}, B::CuArray{FT,3}) where {FT}
    temp = similar(A)
    pivot, info   = CUBLAS.getrf_strided_batched!(A, true);
    getri_strided_batched(A, temp, pivot) # inv stored in aux1
    X = temp ⊠ B  
    # synchronize()
end

# batch Matrix inversion for CPU
function batch_solve!(X::AbstractArray{FT,3}, A::AbstractArray{FT,3}, B::AbstractArray{FT,3}) where {FT}
    for i = 1:size(A, 3)
        @views ldiv!(X[:,:,i], qr!(A[:,:,i]), B[:,:,i])
    end
end

# CUDA has no strided batched getri, but we can at least avoid constructing costly views (copied this over from gertf)
function getri_strided_batched(A::AbstractArray{Float32,3}, C::AbstractArray{Float32,3}, pivotArray::CuMatrix{Cint})
    m, n = size(A, 1), size(A, 2)
    if m != n
        throw(DimensionMismatch("All matrices must be square!"))
    end
    n = size(A, 1)
    ldc = max(1, stride(C, 2))
    lda = max(1, stride(A, 2))
    info = CUDA.zeros(Cint, size(A, 3))
    Cptrs = CUBLAS.unsafe_strided_batch(C)
    Aptrs = CUBLAS.unsafe_strided_batch(A)
    CUBLAS.cublasSgetriBatched(CUBLAS.handle(), n, Aptrs, lda, pivotArray, Cptrs, ldc, info, size(A, 3))
    return nothing
end

# CUDA has no strided batched getri, but we can at least avoid constructing costly views (copied this over from gertf)
function getri_strided_batched(A::AbstractArray{Float64,3}, C::AbstractArray{Float64,3}, pivotArray::CuMatrix{Cint})
    m, n = size(A, 1), size(A, 2)
    if m != n
        throw(DimensionMismatch("All matrices must be square!"))
    end
    n = size(A, 1)
    ldc = max(1, stride(C, 2))
    lda = max(1, stride(A, 2))
    info = CUDA.zeros(Cint, size(A, 3))
    # Cptrs = CUBLAS.unsafe_strided_batch(C)
    Cptrs = unsafe_strided_batch_fake(C, size(A, 3))
    Aptrs = CUBLAS.unsafe_strided_batch(A)
    CUBLAS.cublasDgetriBatched(CUBLAS.handle(), n, Aptrs, lda, pivotArray, Cptrs, ldc, info, size(A, 3))
    return nothing
end

@inline function unsafe_strided_batch_fake(strided::CuArray{T}, batchsize::Int) where {T}
   # batchsize = last(size(strided))
   stride = prod(size(strided)[1:end - 1])
   ptrs = [pointer(strided, stride + 1) for i in 1:batchsize]
   return CuArray(ptrs)
end