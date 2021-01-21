using RadiativeTransfer.RTM
using RadiativeTransfer
using RadiativeTransfer.PhaseFunction
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
FT = Float32
n = 32
nSpec = 20000

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
I_static_ = Diagonal{FT}(ones(n))
I_static  = Diagonal(CuArray(I_static_));

# This should come from poltype later:
D_IQUV  = FT[1, 1, -1, -1]

# D = Diagonal(repeat(pol_type.D, size(qp_μ)[1]))
# Needs to be using Architecture framework in rtm main code and pol_type.D as above:
# Pol numbers:
nP = 4
D_ = (zeros(FT, n, n, 1));
_arr = (repeat(D_IQUV, n ÷ nP));
for i = 1:n
    D_[i,i,1] = _arr[i];
end
D = CuArray(D_)


function rt_interaction!(R⁻⁺::AbstractArray{FT,3}, 
                         T⁺⁺::AbstractArray{FT,3}, 
                         R⁺⁻::AbstractArray{FT,3}, 
                         T⁻⁻::AbstractArray{FT,3}, 
                         r⁻⁺::AbstractArray{FT,3}, 
                         t⁺⁺::AbstractArray{FT,3}, 
                         r⁺⁻::AbstractArray{FT,3}, 
                         t⁻⁻::AbstractArray{FT,3},
                         I_static::AbstractArray) where {FT}
    aux1 = similar(R⁻⁺)
    aux2 = similar(R⁻⁺)
    aux3 = similar(R⁻⁺)
    # Compute M1 = (I - R⁺⁻ * r⁻⁺) \ T⁺⁺
    aux1 = I_static .- R⁺⁻ ⊠ r⁻⁺
    batch_solve!(aux2, aux1, T⁺⁺)
  
    # Compute t_R⁻⁺ = R⁻⁺ + T⁻⁻ * r⁻⁺ * M1
    aux1 = r⁻⁺ ⊠ aux2  # r⁻⁺ * M1
    aux3 = T⁻⁻ ⊠ aux1  # 
    R⁻⁺  = R⁻⁺ + aux3
    
    # T⁺⁺ = t⁺⁺ * M1
    T⁺⁺ = t⁺⁺ ⊠ aux2
    
    # Repeating for mirror-reflected directions
    # Compute M1 = (I - r⁻⁺ * R⁺⁻) \ t⁻⁻
    aux1 = I_static .- r⁻⁺ ⊠ R⁺⁻
    batch_solve!(aux2, aux1, t⁻⁻)

    # t_R⁺⁻ = r⁺⁻ + t⁺⁺ * R⁺⁻ * M1
    aux3 = R⁺⁻ ⊠ aux2
    aux1 = t⁺⁺ ⊠ aux3
    R⁺⁻  = r⁺⁻ + aux1
    T⁻⁻  = T⁺⁺ ⊠ aux2
    
    # synchronize()
end

function rt_doubling!(ndoubl::Int, 
                      r⁻⁺::AbstractArray{FT,3}, 
                      t⁺⁺::AbstractArray{FT,3}, 
                      r⁺⁻::AbstractArray{FT,3}, 
                      t⁻⁻::AbstractArray{FT,3},
                      D::AbstractArray{FT,3},
                      I_static::AbstractArray) where {FT}
    # # ToDo: Important output doubling applied to elemental layer, using same variables r⁺⁻, r⁻⁺, t⁻⁻, t⁺⁺ (can be renamed to t⁺⁺, etc)
    # Need to check with paper nomenclature. This is basically eqs. 23-28 in vSmartMOM but using simplifications in eq. 29-32)
    
    aux1 = similar(t⁺⁺)
    aux2 = similar(t⁺⁺)
    aux3 = similar(t⁺⁺)

    for n = 1:ndoubl

        # M1 = (I - r⁻⁺ * r⁻⁺) \ t⁺⁺
        aux1 = I_static .- r⁻⁺ ⊠ r⁻⁺     # (I - r⁻⁺ * r⁻⁺)      
        batch_solve!(aux2, aux1, t⁺⁺)   # M1 = (I - r⁻⁺ * r⁻⁺) \ t⁺⁺

        # r⁻⁺[:] = r⁻⁺ + t⁺⁺ * r⁻⁺ * M1
        aux1 = r⁻⁺ ⊠ aux2               # r⁻⁺ * M1
        aux3 = t⁺⁺ ⊠ aux1               # t⁺⁺ * r⁻⁺ * M1
        r⁻⁺  = r⁻⁺ + aux3               # r⁻⁺[:] = r⁻⁺ + t⁺⁺ * r⁻⁺ * M1

        # t⁺⁺[:] = t⁺⁺ * M1 
        aux1 = t⁺⁺ ⊠ aux2           # t⁺⁺ * M1 
        t⁺⁺  = aux1                   # t⁺⁺[:] = t⁺⁺ * M1 
    end
    # After doubling, revert D(DR)->R, where D = Diagonal{1,1,-1,-1}
    r⁻⁺ = D ⊠ r⁻⁺ ⊠ D
    # Using r⁺⁻ = Dr⁻⁺D
    r⁺⁻ = D ⊠ r⁻⁺ ⊠ D
    # Using t⁻⁻ = Dt⁺⁺D
    t⁻⁻ = D ⊠ t⁺⁺ ⊠ D
    return nothing 
end

function rt_elemental!(pol_type, dτ, ϖ, Z⁺⁺, Z⁻⁺, m, 
                       ndoubl::Int, qp_μ, wt_μ, 
                       r⁻⁺::AbstractArray{FT,3}, 
                       t⁺⁺::AbstractArray{FT,3}, 
                       r⁺⁻::AbstractArray{FT,3}, 
                       t⁻⁻::AbstractArray{FT,3}, 
                       D::AbstractArray{FT,3}, 
                       I_static::AbstractArray) where {FT}


    aux1 = similar(Z⁻⁺)

    qp_μ4 = repeat(qp_μ, 3)
    wt_μ4 = repeat(wt_μ, 3)
    # println(typeof(qp_μ4))
    #  reduce(vcat, (fill.(qp_μ,[pol_type.n])))

    # wt_μ4 = reduce(vcat, (fill.(wt_μ,[pol_type.n])))
    Nquadn = length(qp_μ4)

    wct = m==0 ? 0.50 * ϖ * wt_μ4  : 0.25 * ϖ * wt_μ4

    d_qp = r⁺⁻ isa CuArray{Float32,3} ? CuArray(Diagonal(1 ./ Array(qp_μ4))) : Diagonal(1 ./ Array(qp_μ4))
    # d_wct = r⁺⁻ isa CuArray{Float32,3} ? CuArray(Diagonal(1 ./ Array(qp_μ4))) : Diagonal(1 ./ Array(qp_μ4))
    d_wct = r⁺⁻ isa CuArray{Float32,3} ? CuArray(Diagonal(Array(wct))) : Diagonal(Array(wct))

    aux1 = d_qp ⊠ Z⁻⁺
    r⁻⁺ = aux1 ⊠ (d_wct * dτ)

    t⁺⁺ = I_static .- (d_qp ⊠ ((I_static .- Z⁺⁺ ⊠ d_wct) * dτ))

    if ndoubl<1
        r⁺⁻ = D ⊠ r⁻⁺ ⊠ D
        t⁻⁻ = D ⊠ t⁺⁺ ⊠ D
    else
        r⁻⁺ = D ⊠ r⁻⁺
    end
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

# CUDA has no strided batched getri, but we can at least avoid constructing costly views (copied this over from gertf)
function getri_strided_batched(A::AbstractArray{Float64,3}, C::AbstractArray{Float64,3}, pivotArray::CuMatrix{Cint})
    m, n = size(A, 1), size(A, 2)
    if m != n
        throw(DimensionMismatch("All matrices must be square!"))
    end
    n = size(A, 1)
    ldc = max(1, stride(C, 2))
    lda = max(1, stride(A, 2))
    info = CUDA.zeros(Cint, size(A, 3))
    # Cptrs = CUBLAS.unsafe_strided_batch(C)
    Cptrs = unsafe_strided_batch_fake(C, size(A, 3))
    Aptrs = CUBLAS.unsafe_strided_batch(A)
    CUBLAS.cublasDgetriBatched(CUBLAS.handle(), n, Aptrs, lda, pivotArray, Cptrs, ldc, info, size(A, 3))
    return nothing
end

@show n, nSpec
function run_rt_interaction(R⁻⁺, T⁺⁺, R⁺⁻, T⁻⁻, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻,  I_static)
    rt_interaction!(R⁻⁺, T⁺⁺, R⁺⁻, T⁻⁻, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻,  I_static)
    synchronize()
end

function run_rt_doubling(nDoubling, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻, D, I_static)
    rt_doubling!(nDoubling, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻, D, I_static)
    synchronize()
end

function run_rt_elemental(pol_type, dτ, ϖ, Z⁺⁺, Z⁻⁺, m, 
                          ndoubl, qp_μ, wt_μ, 
                          r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻, D, I_static)

    rt_elemental!(pol_type, dτ, ϖ, Z⁺⁺, Z⁻⁺, m, 
                  ndoubl, qp_μ, wt_μ, 
                  r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻, D, I_static)
    synchronize()
end

# Testing RT Interaction time:
println("RT Interaction GPU time:")
@btime run_rt_interaction(R⁻⁺, T⁺⁺, R⁺⁻, T⁻⁻, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻,  I_static)
println("RT Interaction CPU time (uses multi-threading in LAPACK/BLAS!):")
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

nDoubling = 4
@show nDoubling
# Testing Doubling:
println("RT Doubling GPU time:")
@btime run_rt_doubling(nDoubling, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻, D, I_static)
println("RT Doubling CPU time:")
@btime run_rt_doubling(nDoubling, r⁻⁺_, t⁺⁺_, r⁺⁻_, t⁻⁻_, D_, I_static_)


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


pol_type = Stokes_IQU{FT}()


m = 0
ndoubl = 0
weight = m == 0 ? 0.5 : 1.0
dτ = 2.458321474082143e-7
ϖ = 1.0
Z⁺⁺ = randn(57,57,3)
Z⁻⁺ = randn(57,57,3)
vza = [60., 45., 30., 15., 0., 15., 30., 45., 60.]
vaz = [180., 180., 180., 180., 0., 0., 0., 0., 0.]
sza = 60.
Ltrunc = 14
Nquad, qp_μ, wt_μ = rt_set_streams(RTM.RadauQuad(), Ltrunc, sza, vza);
I_static_ = Diagonal{FT}(ones(57))
I_static  = Diagonal(CuArray(I_static_));

n = 57
D_ = (zeros(FT, n, n, 1));
_arr = (repeat(D_IQUV, n ÷ nP + 1));
for i = 1:n
    D_[i,i,1] = _arr[i];
end
D = CuArray(D_)

println("RT Elemental GPU time:")
@btime run_rt_elemental(pol_type, dτ, ϖ, CuArray(Z⁺⁺), CuArray(Z⁻⁺), m, 
                        ndoubl, CuArray(qp_μ), CuArray(wt_μ), 
                        r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻, D, I_static)

println("RT Elemental CPU time:")
@btime run_rt_elemental(pol_type, dτ, ϖ, Z⁺⁺, Z⁻⁺, m, 
                        ndoubl, qp_μ, wt_μ, 
                        r⁻⁺_, t⁺⁺_, r⁺⁻_, t⁻⁻_, D_, I_static_)

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


# @inline function unsafe_strided_batch_fake(strided::CuArray{T}, batchsize::Int) where {T}
#    # batchsize = last(size(strided))
#    stride = prod(size(strided)[1:end - 1])
#    ptrs = [pointer(strided, stride + 1) for i in 1:batchsize]
#    return CuArray(ptrs)
# end