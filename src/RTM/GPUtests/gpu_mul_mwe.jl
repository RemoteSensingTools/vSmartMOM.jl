using CUDA
# using TensorOperations
using Einsum
using KernelAbstractions
using KernelAbstractions.Extras
using StaticArrays
using BenchmarkTools

FT = Float32
n = 100;
nSpec = 10000;

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
pivot = CuArray(zeros(Int, n, nSpec));



ker!  = rt_interaction!(CUDADevice(),5)
function test_interaction_GPU!(R⁻⁺, T⁺⁺, R⁺⁻, T⁻⁻, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻, M1, M2, M3,pivot, n, nSpec)
    event = ker!(R⁻⁺, T⁺⁺, R⁺⁻, T⁻⁻, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻, M1, M2, M3,pivot, n, ndrange=nSpec)
    wait(CUDADevice(), event)
    synchronize()
    return nothing
end

nSpec = 1000

r⁻⁺ = (randn(FT, n, n, nSpec));
t⁺⁺ = (randn(FT, n, n, nSpec));
r⁺⁻ = (randn(FT, n, n, nSpec));
t⁻⁻ = (randn(FT, n, n, nSpec));
R⁻⁺ = (randn(FT, n, n, nSpec));
T⁺⁺ = (randn(FT, n, n, nSpec));
R⁺⁻ = (randn(FT, n, n, nSpec));
T⁻⁻ = (randn(FT, n, n, nSpec));
M1 = similar(r⁺⁻);
M2 = similar(r⁺⁻);
M3 = similar(r⁺⁻);
pivot = (zeros(Int, n, nSpec));
ker!  = rt_interaction!(CPU(),4)
function test_interaction_GPU!(R⁻⁺, T⁺⁺, R⁺⁻, T⁻⁻, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻, M1, M2, M3,pivot, n, nSpec)
    event = ker!(R⁻⁺, T⁺⁺, R⁺⁻, T⁻⁻, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻, M1, M2, M3,pivot, n, ndrange=nSpec)
    wait(CPU(), event)
    synchronize()
    return nothing
end


ker!  = rt_interaction_rup!(CPU(),4)
function test_interaction_rup!(R⁻⁺, T⁺⁺, R⁺⁻, T⁻⁻, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻, n, nSpec)
    event = ker!(R⁻⁺, T⁺⁺, R⁺⁻, T⁻⁻, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻,  ndrange=nSpec)
    wait(CPU(), event)
    synchronize()
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
A = randn(Float32, n, n);
B = randn(Float32, n, n);
X = randn(Float32, n, n);
Y = randn(Float32, n, n);

# l, u, p = lu(A)
# LU_decomposition!(A,pivot,Val(n))

full_ldiv = A \ B

@time LU_decomposition!(A, pivot, Val(n))
@time LU_solve!(A, B, X, Y, pivot, Val(n)) 

