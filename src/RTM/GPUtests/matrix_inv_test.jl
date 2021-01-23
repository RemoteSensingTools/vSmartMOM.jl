using Revise
using CUDA
using Test
using BenchmarkTools
using LinearAlgebra
using NNlib

include("CUDA_getri.jl")

# Needs some warning if memory is getting too large !
FT = Float32
n = 32
nSpec = 10000

# Create CPU Matrices:
A_ = (randn(FT, n, n, nSpec));
B_ = (randn(FT, n, n, nSpec));
X_ = (randn(FT, n, n, nSpec));

orig_A = deepcopy(A_)
orig_B = deepcopy(B_)
orig_X = deepcopy(X_)

# And move to GPU as CuArray
A = CuArray(A_);
B = CuArray(B_);
X = CuArray(X_);

# batch Matrix inversion for CPU
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
        # @views X[:,:,i] = inv(A[:,:,i])
        @views ldiv!(X[:,:,i], qr!(A[:,:,i]), B[:,:,i])
    end
end  


@time batch_solve!(X, A, B, );
@time batch_solve!(X_, A_, B_, );

# @test X_ ≈ Array(X) rtol = 1e-2

@time batch_solve!(X, A, B, );
@time batch_solve!(X_, A_, B_, );

