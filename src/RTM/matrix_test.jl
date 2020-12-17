using LinearAlgebra
using CUDA
using StaticArrays
using BenchmarkTools
using NNlib

n = 32
nSpec = 10000
A = randn(Float32, n, n, nSpec);
B = randn(Float32, n, n, nSpec);
AC = CuArray(A);
AB = CuArray(B);

@btime batched_mul($A, $B);
@btime batched_mul($AC, $AB);

n = 30
nSpec = 10000
elty = Float32
# generate matrices
A = [rand(elty, n, n) for i in 1:nSpec]
# move to device
d_A = CuArray{elty,2}[]
for i in 1:length(A)
    push!(d_A, CuArray(A[i]))
end

function invGPU(A)
    info, d_C = CUBLAS.matinv_batched(d_A)
    return collect(d_C)
end

function invCPU!(A, C)
    for Cs in 1:length(A)
        C[Cs]   = inv(A[Cs])
    end
    # return C
end

@btime invCPU(A)
@btime invGPU(d_A)