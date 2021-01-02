using CUDA
using Test
# using TensorOperations
using KernelAbstractions
using KernelAbstractions.Extras
using StaticArrays
using BenchmarkTools
using TensorOperations
using LinearAlgebra
using NNlib
CUDA.allowscalar(false)

FT = Float32
n = 32
nSpec = 30000

# Create CPU Matrices:
r⁻⁺_ = (randn(FT, n, n, nSpec));
t⁺⁺_ = (randn(FT, n, n, nSpec));
r⁺⁻_ = (randn(FT, n, n, nSpec));
t⁻⁻_ = (randn(FT, n, n, nSpec));
R⁻⁺_ = (randn(FT, n, n, nSpec));
T⁺⁺_ = (randn(FT, n, n, nSpec));
R⁺⁻_ = (randn(FT, n, n, nSpec));
T⁻⁻_ = (randn(FT, n, n, nSpec));

# And move to GPU as CuArray
r⁻⁺ = CuArray(r⁻⁺_);
t⁺⁺ = CuArray(t⁺⁺_);
r⁺⁻ = CuArray(r⁺⁻_);
t⁻⁻ = CuArray(t⁻⁻_);
R⁻⁺ = CuArray(R⁻⁺_);
T⁺⁺ = CuArray(T⁺⁺_);
R⁺⁻ = CuArray(R⁺⁻_);
T⁻⁻ = CuArray(T⁻⁻_);

# This is still a bit stupid but noy sure how to make this better...
I_static_ = zeros(FT, n, n, nSpec);
for i = 1:nSpec
    for j = 1:n
        I_static_[j,j,i] = 1
    end
end
I_static = CuArray(I_static_);

function rt_interaction!(R⁻⁺::AbstractArray{FT,3}, 
                         T⁺⁺::AbstractArray{FT,3}, 
                         R⁺⁻::AbstractArray{FT,3}, 
                         T⁻⁻::AbstractArray{FT,3}, 
                         r⁻⁺::AbstractArray{FT,3}, 
                         t⁺⁺::AbstractArray{FT,3}, 
                         r⁺⁻::AbstractArray{FT,3}, 
                         t⁻⁻::AbstractArray{FT,3}, 
                         I_static::AbstractArray{FT,3}) where {FT}
    aux1 = similar(R⁻⁺)
    aux2 = similar(R⁻⁺)
    aux3 = similar(R⁻⁺)
    # Compute M1 = (I - R⁺⁻ * r⁻⁺) \ T⁺⁺
    aux1 = I_static - R⁺⁻ ⊠ r⁻⁺
    batch_solve!(aux2, aux1, T⁺⁺)
  
    # Compute t_R⁻⁺ = R⁻⁺ + T⁻⁻ * r⁻⁺ * M1
    aux1 = r⁻⁺ ⊠ aux2  # r⁻⁺ * M1
    aux3 = T⁻⁻ ⊠ aux1  # 
    R⁻⁺  = R⁻⁺ + aux3
    
    # T⁺⁺ = t⁺⁺ * M1
    T⁺⁺ = t⁺⁺ ⊠ aux2
    
    # Repeating for mirror-reflected directions
    # Compute M1 = (I - r⁻⁺ * R⁺⁻) \ t⁻⁻
    aux1 = I_static - r⁻⁺ ⊠ R⁺⁻
    batch_solve!(aux2, aux1, t⁻⁻)

    # t_R⁺⁻ = r⁺⁻ + t⁺⁺ * R⁺⁻ * M1
    aux3 = R⁺⁻ ⊠ aux2
    aux1 = t⁺⁺ ⊠ aux3
    R⁺⁻  = r⁺⁻ + aux1
    T⁻⁻  = T⁺⁺ ⊠ aux2
    # synchronize()
end

# batch Matrix inversion for CPU
function batch_solve!(X::CuArray{FT,3}, A::CuArray{FT,3}, B::CuArray{FT,3}) where {FT}
    temp = similar(A)
    pivot, info   = CUBLAS.getrf_strided_batched!(A, true);
    getri_strided_batched(A, temp, pivot) # inv stored in aux1
    X = temp ⊠ B  # inv(I - R⁺⁻ * r⁻⁺) * T⁺⁺
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

@show n, nSpec
function run_rt_interaction(R⁻⁺, T⁺⁺, R⁺⁻, T⁻⁻, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻,  I_static)
    rt_interaction!(R⁻⁺, T⁺⁺, R⁺⁻, T⁻⁻, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻,  I_static)
    synchronize()
end

# Testing time:
println("GPU time:")
@btime run_rt_interaction(R⁻⁺, T⁺⁺, R⁺⁻, T⁻⁻, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻,  I_static)
println("CPU time (uses multi-threading in LAPACK/BLAS!):")
@btime run_rt_interaction(R⁻⁺_, T⁺⁺_, R⁺⁻_, T⁻⁻_, r⁻⁺_, t⁺⁺_, r⁺⁻_, t⁻⁻_, I_static_)

@testset "GPU-CPU consistency" begin
    for (matGPU, matCPU) in ((R⁻⁺, R⁻⁺_),
                             (T⁺⁺, T⁺⁺_),
                             (R⁺⁻, R⁺⁻_), 
                             (T⁻⁻, T⁻⁻_),
                             (r⁻⁺, r⁻⁺_),
                             (t⁺⁺, t⁺⁺_), 
                             (r⁺⁻, r⁺⁻_), 
                             (t⁻⁻, t⁻⁻_))
        @test Array(matGPU) ≈ matCPU    
    end
end