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
    vSmartMOM.Architectures.synchronize_if_gpu()

    # Invert LU factorization of A
    CUDA.CUBLAS.getri_strided_batched!(A, temp, pivot)

    # X = inv(A) * B
    NNlib.batched_mul!(X, temp, B)
    vSmartMOM.Architectures.synchronize_if_gpu()
end

"Given 3D CuArray A, fill in X[:,:,k] = A[:,:,k] \\ I" 
function vSmartMOM.CoreRT.batch_inv!(X::CuArray{FT,3}, A::CuArray{FT,3}) where {FT}
    # LU-factorize A
    @timeit "getrf_strided" pivot, info = CUDA.CUBLAS.getrf_strided_batched!(A, true)
    vSmartMOM.Architectures.synchronize_if_gpu()
    # Invert LU factorization of A
    @timeit "getri_strided" CUDA.CUBLAS.getri_strided_batched!(A, X, pivot)
    vSmartMOM.Architectures.synchronize_if_gpu()
end

"""
    batch_inv!(X, A, ws::RTWorkspace)

Workspace-aware batch inversion using pre-allocated pivot/info arrays.
Eliminates repeated GPU memory allocations in tight loops.
"""
function vSmartMOM.CoreRT.batch_inv!(X::CuArray{FT,3}, A::CuArray{FT,3}, 
                                      ws::vSmartMOM.CoreRT.RTWorkspace) where {FT}
    # LU-factorize A using pre-allocated pivot/info
    pivot, info = CUDA.CUBLAS.getrf_strided_batched!(A, true)
    vSmartMOM.Architectures.synchronize_if_gpu()
    # Invert LU factorization of A using pre-allocated output
    CUDA.CUBLAS.getri_strided_batched!(A, X, pivot)
    vSmartMOM.Architectures.synchronize_if_gpu()
end

"Given 3D CuArray A with pointers, fill in X[:,:,k] = A[:,:,k] \\ I (Float32 version)" 
function vSmartMOM.CoreRT.batch_inv!(X::CuArray{FT,3}, A::CuArray{FT,3}, Xptrs, Aptrs) where {FT<:Float32}
    n = size(A,1)
    lda = max(1,stride(A,2))
    batchSize = length(Aptrs)
    
    @timeit "info" info = CuArray{Cint}(undef, batchSize)
    @timeit "pivot" pivot = CUDA.zeros(Cint, (n, batchSize))

    CUDA.CUBLAS.cublasSgetrfBatched(CUDA.CUBLAS.handle(), n, Aptrs, lda, pivot, info, batchSize)
    vSmartMOM.Architectures.synchronize_if_gpu()
    
    # Invert LU factorization of A
    CUDA.CUBLAS.cublasSgetriBatched(CUDA.CUBLAS.handle(), n, Aptrs, lda, pivot, Xptrs, lda, info, batchSize)
    vSmartMOM.Architectures.synchronize_if_gpu()
end

"Given 3D CuArray A with pointers, fill in X[:,:,k] = A[:,:,k] \\ I (Float64 version)" 
function vSmartMOM.CoreRT.batch_inv!(X::CuArray{FT,3}, A::CuArray{FT,3}, Xptrs, Aptrs) where {FT<:Float64}
    n = size(A,1)
    lda = max(1,stride(A,2))
    batchSize = length(Aptrs)
    
    @timeit "info" info = CuArray{Cint}(undef, batchSize)
    @timeit "pivot" pivot = CUDA.zeros(Cint, (n, batchSize))

    CUDA.CUBLAS.cublasDgetrfBatched(CUDA.CUBLAS.handle(), n, Aptrs, lda, pivot, info, batchSize)
    vSmartMOM.Architectures.synchronize_if_gpu()
    
    # Invert LU factorization of A
    CUDA.CUBLAS.cublasDgetriBatched(CUDA.CUBLAS.handle(), n, Aptrs, lda, pivot, Xptrs, lda, info, batchSize)
    vSmartMOM.Architectures.synchronize_if_gpu()
end

"""
    make_gpu_rt_workspace(FT, NquadN, nSpec)

Create an RTWorkspace with CuArray-backed temporaries for GPU computation.
"""
function vSmartMOM.CoreRT.make_gpu_rt_workspace(FT::Type, NquadN::Int, nSpec::Int)
    dims3 = (NquadN, NquadN, nSpec)
    dims_J = (NquadN, 1, nSpec)
    
    vSmartMOM.CoreRT.RTWorkspace(
        CuArray(zeros(FT, dims3)),      # gp_refl
        CuArray(zeros(FT, dims3)),      # tt_gp
        CuArray(zeros(FT, dims_J)),     # J₁⁺
        CuArray(zeros(FT, dims_J)),     # J₁⁻
        CUDA.zeros(Cint, NquadN, nSpec), # pivot
        CUDA.zeros(Cint, nSpec),         # info
        CuArray(zeros(FT, dims3)),      # tmp_inv
        CuArray(zeros(FT, dims3)),      # T_inv
        CuArray(zeros(FT, dims3)),      # tmp3d_a
        CuArray(zeros(FT, dims3)),      # tmp3d_b
    )
end

"Batched matrix multiply using CUBLAS (overwrite NNlib definition)"
function vSmartMOM.CoreRT.batched_mul(A::CuArray{FT,3}, B::CuArray{FT,3}) where {FT}
    CUDA.CUBLAS.gemm_strided_batched('N', 'N', A, B)
end

"Batched multiply for 3D CuArray views (materialize to contiguous CuArray first)."
@inline _as_cuarray3(A::SubArray{FT,3,<:CuArray}) where {FT} = copy(A)

function vSmartMOM.CoreRT.batched_mul(A::SubArray{FT,3,<:CuArray}, B::CuArray{FT,3}) where {FT}
    CUDA.CUBLAS.gemm_strided_batched('N', 'N', _as_cuarray3(A), B)
end

function vSmartMOM.CoreRT.batched_mul(A::CuArray{FT,3}, B::SubArray{FT,3,<:CuArray}) where {FT}
    CUDA.CUBLAS.gemm_strided_batched('N', 'N', A, _as_cuarray3(B))
end

function vSmartMOM.CoreRT.batched_mul(A::SubArray{FT,3,<:CuArray}, B::SubArray{FT,3,<:CuArray}) where {FT}
    CUDA.CUBLAS.gemm_strided_batched('N', 'N', _as_cuarray3(A), _as_cuarray3(B))
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
