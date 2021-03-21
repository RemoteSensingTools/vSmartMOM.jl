using RadiativeTransfer.vSmartMOM
using RadiativeTransfer
using RadiativeTransfer.Scattering
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

if has_cuda_gpu()
    CUDA.allowscalar(false)
end

# Needs some warning if memory is getting too large !
FT = Float64
n = 16
nSpec = 100

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
# I_static_ = Diagonal{FT}(ones(n))
# I_static = CuArray(repeat(I_static, 1, 1))
I_static_ = Diagonal{FT}(ones(n))
I_static  = Diagonal(CuArray(I_static_));



function rt_interaction!(R⁻⁺::AbstractArray{FT,3}, 
                         T⁺⁺::AbstractArray{FT,3}, 
                         R⁺⁻::AbstractArray{FT,3}, 
                         T⁻⁻::AbstractArray{FT,3}, 
                         r⁻⁺::AbstractArray{FT,3}, 
                         t⁺⁺::AbstractArray{FT,3}, 
                         r⁺⁻::AbstractArray{FT,3}, 
                         t⁻⁻::AbstractArray{FT,3},
                         I_static::AbstractArray) where {FT}
    # Used to store `(I - R⁺⁻ * r⁻⁺)⁻¹ * T⁺⁺`
    tmp_inv = similar(t⁺⁺)

    # Compute and store `(I - R⁺⁻ * r⁻⁺)⁻¹ * T⁺⁺`
    vSmartMOM.batch_solve!(tmp_inv, I_static .- R⁺⁻ ⊠ r⁻⁺, T⁺⁺)
    # synchronize()
    R⁻⁺ += (T⁻⁻ ⊠ r⁻⁺ ⊠ tmp_inv)
    # synchronize()
    T⁺⁺[:] = t⁺⁺ ⊠ tmp_inv

    # Repeating for mirror-reflected directions

    # Compute and store `(I - r⁻⁺ * R⁺⁻)⁻¹ * t⁻⁻`
    # synchronize()
    vSmartMOM.batch_solve!(tmp_inv, I_static .- r⁻⁺ ⊠ R⁺⁻, t⁻⁻)
    synchronize()
    R⁺⁻ .= r⁺⁻ + t⁺⁺ ⊠ R⁺⁻ ⊠ tmp_inv
    T⁻⁻[:] = tmp_inv
    # synchronize()
    T⁻⁻[:] = T⁺⁺ ⊠ tmp_inv

    return nothing
end


@show n, nSpec
function run_rt_interaction(R⁻⁺, T⁺⁺, R⁺⁻, T⁻⁻, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻,  I_static)
    rt_interaction!(R⁻⁺, T⁺⁺, R⁺⁻, T⁻⁻, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻,  I_static)
    synchronize()
end

# Testing RT Interaction time:
println("RT Interaction GPU time:")
@time run_rt_interaction(R⁻⁺, T⁺⁺, R⁺⁻, T⁻⁻, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻,  I_static)
println("RT Interaction CPU time (uses multi-threading in LAPACK/BLAS!):")
@time run_rt_interaction(R⁻⁺_, T⁺⁺_, R⁺⁻_, T⁻⁻_, r⁻⁺_, t⁺⁺_, r⁺⁻_, t⁻⁻_, I_static_)

@testset "GPU-CPU consistency" begin
    for (matGPU, matCPU) in ((R⁻⁺, R⁻⁺_),
                             (T⁺⁺, T⁺⁺_),
                             (R⁺⁻, R⁺⁻_), 
                             (T⁻⁻, T⁻⁻_),
                             (r⁻⁺, r⁻⁺_),
                             (t⁺⁺, t⁺⁺_), 
                             (r⁺⁻, r⁺⁻_), 
                             (t⁻⁻, t⁻⁻_))
        @show Array(matGPU) ≈ matCPU
        @test Array(matGPU) ≈ matCPU    
    end
end
