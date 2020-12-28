using CUDA
# using TensorOperations
using Einsum
using KernelAbstractions
using KernelAbstractions.Extras
using StaticArrays
using BenchmarkTools

FT = Float32
n = 200;
nSpec = 1000;

C = CuArray(randn(FT, n, n, nSpec));
B = CuArray(randn(FT, n, n, nSpec));
A = CuArray(randn(FT, n, n, nSpec));
X = CuArray(randn(FT, n, n, nSpec));
Y = CuArray(randn(FT, n, n, nSpec));

r⁻⁺ = CuArray(randn(FT, n, n, nSpec));
t⁺⁺ = CuArray(randn(FT, n, n, nSpec));
r⁺⁻ = CuArray(randn(FT, n, n, nSpec));
t⁻⁻ = CuArray(randn(FT, n, n, nSpec));
R⁻⁺ = CuArray(randn(FT, n, n, nSpec));
T⁺⁺ = CuArray(randn(FT, n, n, nSpec));
R⁺⁻ = CuArray(randn(FT, n, n, nSpec));
T⁻⁻ = CuArray(randn(FT, n, n, nSpec));
M1 = similar(r⁺⁻);
M2 = similar(r⁺⁻);
M3 = similar(r⁺⁻);

@kernel function rt_interaction!(R⁻⁺, T⁺⁺, R⁺⁻, T⁻⁻, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻, M1, M2, M3, n)
    N = @index(Global)
    
    # Views:
    R⁻⁺_ = view(R⁻⁺, :, :, N)
    T⁺⁺_ = view(T⁺⁺, :, :, N)
    R⁺⁻_ = view(R⁺⁻, :, :, N)
    T⁻⁻_ = view(T⁻⁻, :, :, N) 
    r⁻⁺_ = view(r⁻⁺, :, :, N)
    t⁺⁺_ = view(t⁺⁺, :, :, N)
    r⁺⁻_ = view(r⁺⁻, :, :, N)
    t⁻⁻_ = view(t⁻⁻, :, :, N)
    M1_  = view(M1, :, :, N)
    M2_  = view(M2, :, :, N)
    M3_  = view(M3, :, :, N)

    # M1 = (I - R⁺⁻ * r⁻⁺) \ T⁺⁺
    mat_multiply!(M1_, R⁺⁻_, r⁻⁺_,  n)
    @inbounds for i = 1:n
        M1_[i,i] += 1;
    end
    mat_inv_mul!(M1_, T⁺⁺_, M3_, M2_, n)

    # t_R⁻⁺ = R⁻⁺ + T⁻⁻ * r⁻⁺ * M1
    mat_multiply!(M2_, r⁻⁺_, M1_,  n)
    mat_multiply!(M3_, T⁻⁻_, M2_,  n)
    @inbounds for i = 1:n
        for j = 1:n
            R⁻⁺_[j,i] = R⁻⁺_[j,i] + M3_[j,i];
        end
    end

end

@kernel function mat_mul!(out, in1, in2, n)
    N = @index(Global)
    mat_multiply!(view(out, :, :, N), view(in1, :, :, N), view(in2, :, :, N),  n)
end

@kernel function mat_LU!(A, n)
    N = @index(Global)
    mat_LU_!(view(A, :, :, N), n)
    # n,  m = size(A[:,:,1]);
    
end

@kernel function mat_inv_multiply!(A, B, X, Y, n)
    N = @index(Global)
    mat_inv_mul!(view(A, :, :, N), view(B, :, :, N), view(X, :, :, N), view(Y, :, :, N), n)    
end

function mat_multiply!(out, in1, in2, n)
    @inbounds begin
        @unroll for i = 1:n # for each column but the last
            for j = 1:n
                out[i,j] = in1[i,1] * in2[1,j];
                for k = 2:n
                    out[i, j] += in1[i, k] * in2[k, j];
                end
            end
        end
    end
end

function mat_LU_!(A, n)
    @inbounds begin
        @unroll for col = 1:n - 1 # for each column but the last
            piv = A[col, col]; # pivot
          # col1 = col + 1;
            ipiv = 1 / piv
            @unroll for r = col + 1:n # col of L
                A[r, col] *= ipiv
            end
          # update submatrix A[col1:n,col1:n]
            @unroll for c = col + 1:n
                nacolc = -A[col, c]
                @unroll for r = col + 1:n
                    A[r, c] += A[r, col] * nacolc;
                end
            end
        end
    end
end

# Compute X = A⁻¹ * B
function mat_inv_mul!(A, B, X, Y, n)
    FT = eltype(A)
    # Naive LU decomposition:
    fill!(X, 0)
    fill!(Y, 0)
    @inbounds begin
        @unroll for i = 1:n 
            for j = 1:(i - 1)
                alpha = A[i,j]
                for k = 1:(j - 1)
                    alpha = alpha - A[i, k] * A[k, j];
                end
                A[i, j] = alpha / A[j, j];
            end
            for j = i:n
                alpha = A[i, j];
                for k = 1:(i - 1)
                    alpha = alpha - A[i, k] * A[k, j];
                end
                A[i, j] = alpha;
            end
        end
        #  A = L+U-I
        
        for col = 1:n
            # find solution of Ly = b   
            for i = 1:n
                alpha = FT(0);
                for k = 1:i
                    alpha = alpha + A[i, k] * Y[k,col];
                end
                Y[i, col] = B[i,col] - alpha;
            end
            # find solution of Ux = y
            for i = n:(-1):1
                alpha = FT(0);
                for k = (i + 1):1:n
                    alpha = alpha + A[i, k] * X[k,col];
                end
                X[i,col] = (Y[i,col] - alpha) / A[i, i];
            end
        end
    end
end

function testLU!(A, n)
    event = ker2!(A, n, ndrange=nSpec)
    wait(CUDADevice(), event)
    return nothing
end

ker!  = mat_mul!(CUDADevice())
function testGPU_mult!(A, B, C, n)
    event = ker!(A, B, C, n, ndrange=nSpec)
    wait(CUDADevice(), event)
    return nothing
end
@time testGPU_mult!(A, B, C, n)
@time testGPU_mult!(A, B, C, n)

ker2! = mat_inv_multiply!(CUDADevice())
function testGPU_inv!(A, B, X, Y, n)
    event = ker2!(A, B, X, Y, n, ndrange=nSpec)
    wait(CUDADevice(), event)
    return nothing
end
@time testGPU_inv!(A, B, X, Y, n)
@time testGPU_inv!(A, B, X, Y, n)


ker!  = rt_interaction!(CUDADevice())
function test_interaction_GPU!(R⁻⁺, T⁺⁺, R⁺⁻, T⁻⁻, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻, M1, M2, M3, n, nSpec)
    event = ker!(R⁻⁺, T⁺⁺, R⁺⁻, T⁻⁻, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻, M1, M2, M3, n, ndrange=nSpec)
    wait(CUDADevice(), event)
    return nothing
end
# ker!  = mat_test!(CUDADevice())
# ker2! = mat_LU!(CUDADevice())

# @time testGPU!(A, B, C)

# Bechmark against 
# A = [rand(FT, n, n) for i in 1:nSpec]
# move to device
# d_A = CuArray{FT,2}[]
# for i in 1:length(A)
#    push!(d_A, CuArray(A[i]))
# nd

# function invGPU(d_A)
#    info, d_C = CUBLAS.matinv_batched(d_A)
#    return collect(d_C)
# end

# Test method:
# A = randn(Float32, n, n);
# B = randn(Float32, n, n);
# X = randn(Float32, n, n);
# Y = randn(Float32, n, n);

# full_ldiv = A \ B

# mat_inv_mul!(A, B, X, Y, n)
