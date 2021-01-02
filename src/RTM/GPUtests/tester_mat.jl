using CUDA
# using TensorOperations
using Einsum
using KernelAbstractions
using KernelAbstractions.Extras
using StaticArrays
using BenchmarkTools

FT = Float32
n = 400;
nSpec = 1000;

array_type = CuArray
include("src/RTM/matrix_linear_algebra.jl")

C = (randn(FT, n, n, nSpec));
B = (randn(FT, n, n, nSpec));
A = (randn(FT, n, n, nSpec));
X = (randn(FT, n, n, nSpec));
Y = (randn(FT, n, n, nSpec));
pivot = (zeros(Int, n, nSpec));

C_ = array_type(C)
B_ = array_type(B)
A_ = array_type(A)
X_ = array_type(X)
Y_ = array_type(Y)
pivot_ = array_type(pivot)


tester = (zeros(FT, n, n, nSpec));
tester_ = array_type(tester)
for i = 1:nSpec
    tester[:,:,i]  = A[:,:,i] \ B[:,:,i]
    tester_[:,:,i] = A_[:,:,i] \ B_[:,:,i]
end

@kernel function mat_inv_multiply!(A, B, X, Y, pivot, n)
    N = @index(Global)
    LU_decomposition!(view(A,:,:,N), view(pivot,:,N), n)
    #LU_solve!(view(A,:,:,N), view(B,:,:,N),view(X,:,:,N), view(Y,:,:,N), view(pivot,:,N), n)
    #@synchronize    
end

device = CUDADevice()
ker2! = mat_inv_multiply!(device,1024)
function testGPU_inv!(A, B, X, Y, pivot, n)
    event = ker2!(A, B, X, Y, pivot, n, ndrange=nSpec)
    wait(device, event)
    synchronize()
    return nothing
end
@time testGPU_inv!(A_, B_, X_, Y_, pivot_, n)

for i = 1:nSpec
    T = Array(X_[:,:,i])
    diff = (T .- tester[:,:,i])

    @show i, tester[:,:,i] ≈ Array(X_[:,:,i]), maximum(diff), maximum(diff ./ tester[:,:,i] * 100)
end
    


device = CPU()
ker2! = mat_inv_multiply!(device)
function testGPU_inv!(A, B, X, Y, pivot, n)
    event = ker2!(A, B, X, Y, pivot, n, ndrange=nSpec)
    wait(device, event)
    return nothing
end
@time testGPU_inv!(A, B, X, Y, pivot, n)