using CUDA
# using TensorOperations
using KernelAbstractions
using KernelAbstractions.Extras
using StaticArrays
using BenchmarkTools
using TensorOperations
using LinearAlgebra
CUDA.allowscalar(false)

FT = Float32
n = 64
nSpec = 2^14


r⁻⁺_ = (randn(FT, n, n, nSpec));
t⁺⁺_ = (randn(FT, n, n, nSpec));
r⁺⁻_ = (randn(FT, n, n, nSpec));
t⁻⁻_ = (randn(FT, n, n, nSpec));
R⁻⁺_ = (randn(FT, n, n, nSpec));
T⁺⁺_ = (randn(FT, n, n, nSpec));
R⁺⁻_ = (randn(FT, n, n, nSpec));
T⁻⁻_ = (randn(FT, n, n, nSpec));
M1_ = (r⁺⁻_);
M2_ = (r⁺⁻_);
M3_ = (r⁺⁻_);

r⁻⁺ = CuArray(r⁻⁺_);
t⁺⁺ = CuArray(t⁺⁺_);
r⁺⁻ = CuArray(r⁺⁻_);
t⁻⁻ = CuArray(t⁻⁻_);
R⁻⁺ = CuArray(R⁻⁺_);
T⁺⁺ = CuArray(T⁺⁺_);
R⁺⁻ = CuArray(R⁺⁻_);
T⁻⁻ = CuArray(T⁻⁻_);
M1 = similar(r⁺⁻);
M2 = similar(r⁺⁻);
M3 = similar(r⁺⁻);

I_static_ = zeros(FT, n, n, nSpec);
for i=1:nSpec
    for j=1:n
        I_static_[j,j,i] = 1
    end
end
I_static = CuArray(I_static_);

function rt_interaction_GPU!(R⁻⁺, T⁺⁺, R⁺⁻, T⁻⁻, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻, aux1, aux2, aux3,I_static)

    # Compute M1 = (I - R⁺⁻ * r⁻⁺) \ T⁺⁺
    aux2 = CUBLAS.gemm_strided_batched('N','N',R⁺⁻, r⁻⁺) #mat_multiply!(aux1, R⁺⁻, r⁻⁺,n)                 # aux1  =  R⁺⁻ * r⁻⁺
    aux2 = I_static - aux2
    
    pivot, info      = CUBLAS.getrf_strided_batched!(aux2, true)
    getri_strided_batched(aux2, aux1, pivot) # inv stored in aux1
    aux2 = CUBLAS.gemm_strided_batched('N','N',aux1,T⁺⁺)  # inv(I - R⁺⁻ * r⁻⁺) * T⁺⁺
  
    # Compute t_R⁻⁺ = R⁻⁺ + T⁻⁻ * r⁻⁺ * M1
    aux1 = CUBLAS.gemm_strided_batched('N','N',r⁻⁺, aux2)  # r⁻⁺ * M1
    aux3 = CUBLAS.gemm_strided_batched('N','N',T⁻⁻, aux1)  # 
    R⁻⁺ = R⁻⁺ + aux3
    
    # t_T⁺⁺ = t⁺⁺ * M1
    T⁺⁺ = CUBLAS.gemm_strided_batched('N','N',t⁺⁺, aux2)
    
    # Repeating for mirror-reflected directions
    # Compute M1 = (I - r⁻⁺ * R⁺⁻) \ t⁻⁻
    aux2 = CUBLAS.gemm_strided_batched('N','N',r⁻⁺, R⁺⁻)
    aux2 = I_static -aux2 
    pivot, info      = CUBLAS.getrf_strided_batched!(aux2, true)
    getri_strided_batched(aux2, aux1, pivot) # inv stored in aux1
    aux2 = CUBLAS.gemm_strided_batched('N','N',aux1,t⁻⁻)

    # t_R⁺⁻ = r⁺⁻ + t⁺⁺ * R⁺⁻ * M1
    aux3 = CUBLAS.gemm_strided_batched('N','N',R⁺⁻, aux2)
    aux1 = CUBLAS.gemm_strided_batched('N','N',t⁺⁺, aux3)
    R⁺⁻ =  r⁺⁻ + aux1
    T⁻⁻ = CUBLAS.gemm_strided_batched('N','N',T⁺⁺, aux2)
end

# CUDA has no strided batched getrf, but we can at least avoid constructing costly views
function getri_strided_batched(A::AbstractArray{Float32, 3},C::AbstractArray{Float32, 3},pivotArray::CuMatrix{Cint})
    m,n = size(A,1), size(A,2)
    if m != n
        throw(DimensionMismatch("All matrices must be square!"))
    end
    n = size(A,1)
    ldc = max(1,stride(C,2))
    lda = max(1,stride(A,2))
    info = CUDA.zeros(Cint,size(A,3))
    Cptrs = CUBLAS.unsafe_strided_batch(C)
    Aptrs = CUBLAS.unsafe_strided_batch(A)
    CUBLAS.cublasSgetriBatched(CUBLAS.handle(), n, Aptrs, lda, pivotArray, Cptrs, ldc, info, size(A,3))
    return nothing
end

@show n, nSpec
function rt_gpu(R⁻⁺, T⁺⁺, R⁺⁻, T⁻⁻, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻, M1, M2, M3,I_static)
rt_interaction_GPU!(R⁻⁺, T⁺⁺, R⁺⁻, T⁻⁻, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻, M1, M2, M3,I_static)
synchronize()
end
@time rt_gpu(R⁻⁺, T⁺⁺, R⁺⁻, T⁻⁻, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻, M1, M2, M3,I_static)


function rt_interaction!(R⁻⁺, T⁺⁺, R⁺⁻, T⁻⁻, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻,aux1,aux2,aux3,I_static)
    # M1 = (I - R⁺⁻ * r⁻⁺) \ T⁺⁺;aux1 = similar(R⁺⁻)

        mul!(aux1, R⁺⁻, r⁻⁺)        # R⁺⁻ * r⁻⁺
        aux1[:] = I_static - aux1   # (I - R⁺⁻ * r⁻⁺)
        ldiv!(aux2, qr!(aux1), T⁺⁺) # M1 = (I - R⁺⁻ * r⁻⁺) \ T⁺⁺

        # t_R⁻⁺ = R⁻⁺ + T⁻⁻ * r⁻⁺ * M1
        #@show size(aux1), size(r⁻⁺), size(aux2)
        aux1[:] = r⁻⁺ * aux2
        aux3[:] = T⁻⁻ * aux1
        #mul!(aux1, r⁻⁺, aux2)   # r⁻⁺ * M1
        #mul!(aux3, T⁻⁻, aux1)   # T⁻⁻ * r⁻⁺ * M1
        R⁻⁺[:] = R⁻⁺ + aux3     # t_R⁻⁺ = R⁻⁺ + T⁻⁻ * r⁻⁺ * M1
        
        # t_T⁺⁺ = t⁺⁺ * M1
        mul!(T⁺⁺, t⁺⁺, aux2)

        # Repeating for mirror-reflected directions

        # M1 = (I - r⁻⁺ * R⁺⁻) \ t⁻⁻
        mul!(aux1, r⁻⁺, R⁺⁻)        # r⁻⁺ * R⁺⁻
        aux1[:] = I_static - aux1   # (I - r⁻⁺ * R⁺⁻)
        ldiv!(aux2, qr!(aux1), t⁻⁻) # M1 = (I - r⁻⁺ * R⁺⁻) \ t⁻⁻

        # t_R⁺⁻ = r⁺⁻ + t⁺⁺ * R⁺⁻ * M1
        aux3[:] = R⁺⁻ * aux2
        aux1[:] = t⁺⁺ * aux3
        #mul!(aux3, R⁺⁻, aux2)   # R⁺⁻ * M1
        #mul!(aux1, t⁺⁺, aux3)   # t⁺⁺ * R⁺⁻ * M1
        R⁺⁻[:] = r⁺⁻ + aux1     # t_R⁺⁻ = r⁺⁻ + t⁺⁺ * R⁺⁻ * M1

        # t_T⁻⁻ = T⁺⁺ * M1
        mul!(T⁻⁻, T⁺⁺, aux2)
                 
        return nothing 
end

function runAllCPU!(R⁻⁺, T⁺⁺, R⁺⁻, T⁻⁻, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻,aux1,aux2,aux3, I_static)
    for i=1:size(R⁻⁺,3)
        rt_interaction!(view(R⁻⁺,:,:,i), view(T⁺⁺,:,:,i), view(R⁺⁻,:,:,i), view(T⁻⁻,:,:,i), view(r⁻⁺,:,:,i), view(t⁺⁺,:,:,i), view(r⁺⁻,:,:,i), view(t⁻⁻,:,:,i),view(aux1,:,:,i),view(aux2,:,:,i),view(aux3,:,:,i), view(I_static,:,:,i))
    end
    return nothing
end

@time runAllCPU!(R⁻⁺_, T⁺⁺_, R⁺⁻_, T⁻⁻_, r⁻⁺_, t⁺⁺_, r⁺⁻_, t⁻⁻_,M1_,M2_,M3_, I_static_)