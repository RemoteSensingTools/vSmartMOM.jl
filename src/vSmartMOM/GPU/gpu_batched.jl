# Batched Linear Algebra code

"Given 3D CuArrays A and B, fill in X[:,:,k] = A[:,:,k] \\ B[:,:,k]" 
function batch_solve!(X::CuArray{FT,3}, A::CuArray{FT,3}, B::CuArray{FT,3}) where {FT}

    # Temporary factorization matrix
    temp = similar(A)

    # LU-factorize A
    pivot, info   = CUBLAS.getrf_strided_batched!(A, true);
    synchronize()

    # Invert LU factorization of A
    getri_strided_batched!(A, temp, pivot);

    # X = inv(A) * B
    NNlib.batched_mul!(X, temp, B)
    synchronize()

end

"Given 3D Julia Arrays A and B, fill in X[:,:,k] = A[:,:,k] \\ B[:,:,k]" 
function batch_solve!(X::AbstractArray{FT,3}, A::AbstractArray{FT,3}, B::AbstractArray{FT,3}) where {FT}
    Threads.@threads for i = 1:size(A, 3)
        @views ldiv!(X[:,:,i], qr!(A[:,:,i]), B[:,:,i])
    end
end

"Given 3D CuArray A, fill in X[:,:,k] = A[:,:,k] \\ I" 
function batch_inv!(X::CuArray{FT,3}, A::CuArray{FT,3}) where {FT}

    # LU-factorize A
    pivot, info   = CUBLAS.getrf_strided_batched!(A, true);
    synchronize()

    # Invert LU factorization of A
    getri_strided_batched!(A, X, pivot); 
    synchronize()

    return nothing
end

"Given 3D Julia Array A, fill in X[:,:,k] = A[:,:,k] \\ I" 
function batch_inv!(X::AbstractArray{FT,3}, A::AbstractArray{FT,3}) where {FT}
    Threads.@threads for i = 1:size(A, 3)
        @views X[:,:,i] = A[:,:,i]\I;
    end
end

"Batched matrix multiply (overwrite NNlib definition)"
function batched_mul(A::CuArray{FT,3}, B::CuArray{FT,3}) where {FT}
    CUBLAS.gemm_strided_batched('N', 'N', A, B)
end
