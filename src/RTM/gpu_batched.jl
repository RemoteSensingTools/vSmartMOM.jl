# batch Matrix inversion for GPU
# function batch_solve!(X::CuArray{FT,3}, A::CuArray{FT,3}, B::CuArray{FT,3}) where {FT}
#    temp = similar(A)
#    pivot, info   = CUBLAS.getrf_strided_batched!(A, true);
#    getri_strided_batched(A, temp, pivot) # inv stored in aux1
#    X = temp ⊠ B  
#    # synchronize()
# end

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
    for i = 1:size(A, 3)
        @views ldiv!(X[:,:,i], qr!(A[:,:,i]), B[:,:,i])
    end
end

