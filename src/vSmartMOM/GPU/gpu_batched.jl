# Batched Linear Algebra code

function batch_solve!(X::CuArray{FT,3}, A::CuArray{FT,3}, B::CuArray{FT,3}) where {FT}
    temp = similar(A)
    pivot, info   = CUBLAS.getrf_strided_batched!(A, true);
    synchronize()
    getri_strided_batched!(A, temp, pivot); # inv stored in aux1
    NNlib.batched_mul!(X, temp, B)
    synchronize()
end

# batch Matrix inversion for CPU
function batch_solve!(X::AbstractArray{FT,3}, A::AbstractArray{FT,3}, B::AbstractArray{FT,3}) where {FT}
    Threads.@threads for i = 1:size(A, 3)
        @views ldiv!(X[:,:,i], qr!(A[:,:,i]), B[:,:,i])
    end
end

function batch_inv!(X::CuArray{FT,3}, A::CuArray{FT,3}) where {FT}
    pivot, info   = CUBLAS.getrf_strided_batched!(A, true);
    synchronize()
    getri_strided_batched!(A, X, pivot); # inv stored in aux1
    synchronize()
    return nothing
end

# batch Matrix inversion for CPU
function batch_inv!(X::AbstractArray{FT,3}, A::AbstractArray{FT,3}) where {FT}
    Threads.@threads for i = 1:size(A, 3)
        @views X[:,:,i] = A[:,:,i]\I;
    end
end

# batched matrix multiply (overwrite NNlib definition)
function batched_mul(A::CuArray{FT,3}, B::CuArray{FT,3}) where {FT}
    CUBLAS.gemm_strided_batched('N', 'N', A, B)
end
