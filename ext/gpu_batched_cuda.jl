# GPU-specific batched operations for vSmartMOM
# This file defines CuArray-specific methods when CUDA is available

using TimerOutputs
using NNlib
using ForwardDiff

# Import the ⊠ operator (batched multiply) from CoreRT if needed
import NNlib: batched_mul

# Note: synchronize() in CoreRT already uses Architectures.synchronize_if_gpu()
# which calls CUDA.synchronize() when CUDA is available (via _sync_gpu Ref)

"Given 3D CuArrays A and B, fill in X[:,:,k] = A[:,:,k] \\ B[:,:,k]" 
function vSmartMOM.CoreRT.batch_solve!(X::CuArray{FT,3}, A::CuArray{FT,3}, B::CuArray{FT,3}) where {FT}
    # Temporary factorization matrix
    temp = similar(A)

    # LU-factorize A
    pivot, info = CUDA.CUBLAS.getrf_strided_batched!(A, true)
    vSmartMOM.CoreRT.synchronize()

    # Invert LU factorization of A
    CUDA.CUBLAS.getri_strided_batched!(A, temp, pivot)

    # X = inv(A) * B
    NNlib.batched_mul!(X, temp, B)
    vSmartMOM.CoreRT.synchronize()
end

"Given 3D CuArray A, fill in X[:,:,k] = A[:,:,k] \\ I" 
function vSmartMOM.CoreRT.batch_inv!(X::CuArray{FT,3}, A::CuArray{FT,3}) where {FT}
    # LU-factorize A
    @timeit "getrf_strided" pivot, info = CUDA.CUBLAS.getrf_strided_batched!(A, true)
    vSmartMOM.CoreRT.synchronize()
    # Invert LU factorization of A
    @timeit "getri_strided" CUDA.CUBLAS.getri_strided_batched!(A, X, pivot)
    vSmartMOM.CoreRT.synchronize()
end

"Given 3D CuArray A with pointers, fill in X[:,:,k] = A[:,:,k] \\ I (Float32 version)" 
function vSmartMOM.CoreRT.batch_inv!(X::CuArray{FT,3}, A::CuArray{FT,3}, Xptrs, Aptrs) where {FT<:Float32}
    n = size(A,1)
    lda = max(1,stride(A,2))
    batchSize = length(Aptrs)
    
    @timeit "info" info = CuArray{Cint}(undef, batchSize)
    @timeit "pivot" pivot = CUDA.zeros(Cint, (n, batchSize))

    CUDA.CUBLAS.cublasSgetrfBatched(CUDA.CUBLAS.handle(), n, Aptrs, lda, pivot, info, batchSize)
    vSmartMOM.CoreRT.synchronize()
    
    # Invert LU factorization of A
    CUDA.CUBLAS.cublasSgetriBatched(CUDA.CUBLAS.handle(), n, Aptrs, lda, pivot, Xptrs, lda, info, batchSize)
    vSmartMOM.CoreRT.synchronize()
end

"Given 3D CuArray A with pointers, fill in X[:,:,k] = A[:,:,k] \\ I (Float64 version)" 
function vSmartMOM.CoreRT.batch_inv!(X::CuArray{FT,3}, A::CuArray{FT,3}, Xptrs, Aptrs) where {FT<:Float64}
    n = size(A,1)
    lda = max(1,stride(A,2))
    batchSize = length(Aptrs)
    
    @timeit "info" info = CuArray{Cint}(undef, batchSize)
    @timeit "pivot" pivot = CUDA.zeros(Cint, (n, batchSize))

    CUDA.CUBLAS.cublasDgetrfBatched(CUDA.CUBLAS.handle(), n, Aptrs, lda, pivot, info, batchSize)
    vSmartMOM.CoreRT.synchronize()
    
    # Invert LU factorization of A
    CUDA.CUBLAS.cublasDgetriBatched(CUDA.CUBLAS.handle(), n, Aptrs, lda, pivot, Xptrs, lda, info, batchSize)
    vSmartMOM.CoreRT.synchronize()
end

"Batched matrix multiply using CUBLAS (overwrite NNlib definition)"
function vSmartMOM.CoreRT.batched_mul(A::CuArray{FT,3}, B::CuArray{FT,3}) where {FT}
    CUDA.CUBLAS.gemm_strided_batched('N', 'N', A, B)
end

"Define batched matrix multiply for GPU and Duals"
function vSmartMOM.CoreRT.batched_mul(A::CuArray{ForwardDiff.Dual{T,V,N},3}, B::CuArray{ForwardDiff.Dual{T,V,N},3}) where {T,V,N}
    # Extract values:
    Av = ForwardDiff.value.(A)
    Bv = ForwardDiff.value.(B)
    
    # Use strided batch for A*B (using ⊠ operator which calls gemm_strided_batched):
    Cv = NNlib.batched_mul(Av, Bv)
    
    # Compute derivatives ∂(AB)/∂x = A * ∂B/∂x + ∂A/∂x * B;
    dABdx = [NNlib.batched_mul(Av, ForwardDiff.partials.(B,i)) + 
             NNlib.batched_mul(ForwardDiff.partials.(A,i), Bv) for i=1:N]
    dABdx = ForwardDiff.Partials.(tuple.(dABdx...))
    
    return eltype(A).(Cv, dABdx)
end

"Overload of batch_inv! for Dual numbers on GPU"
function vSmartMOM.CoreRT.batch_inv!(X::CuArray{ForwardDiff.Dual{T,V,N},3}, A::CuArray{ForwardDiff.Dual{T,V,N},3}) where {T,V,N}
    Atemp = ForwardDiff.value.(A)
    invA  = similar(Atemp)
    
    # Set invA=A⁻¹
    vSmartMOM.CoreRT.batch_inv!(invA, Atemp)

    # Get partial derivatives
    K = [ForwardDiff.partials.(A,i) for i=1:N]

    # Compute derivatives ∂A⁻¹/∂x = -A⁻¹ * ∂A/∂x * A⁻¹; using NNlib batched matrix multiply
    @timeit "InvDerivs" dAdx = [
        -NNlib.batched_mul(NNlib.batched_mul(invA, K[i]), invA) for i=1:N
    ]
    
    # Pack into tuples again
    dAdx = ForwardDiff.Partials.(tuple.(dAdx...))
    X .= eltype(X).(invA, dAdx)
end
