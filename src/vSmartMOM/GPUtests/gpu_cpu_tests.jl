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
BLAS.set_num_threads(1)

if has_cuda_gpu()
    CUDA.allowscalar(false)
end

# Needs some warning if memory is getting too large !
FT = Float64
n = 32
nSpec  = 20000
ndoubl = 5

dims = (n,n)
pol_type  = Stokes_IQU{FT}()
I_static  = Diagonal{FT}(ones(n));
I_static_ = Diagonal(CuArray(Diagonal{FT}(ones(n))));

SFI = true
expk  = exp.(-randn(nSpec));
expk_ = CuArray(expk);

added_layer_CPU     = vSmartMOM.make_added_layer_rand(FT, Array, dims, nSpec);
#composite_layer_CPU = vSmartMOM.make_composite_layer(FT, Array, dims, nSpec)
added_layer_GPU     = vSmartMOM.make_added_layer_rand(FT, CuArray, dims, nSpec); 
#composite_layer_GPU = vSmartMOM.make_composite_layer(FT, CuArray, dims, nSpec)

println("CPU runs:")
@btime vSmartMOM.doubling!(pol_type, SFI, expk , ndoubl, added_layer_CPU, I_static , RadiativeTransfer.Architectures.CPU())
println("GPU runs:")
@btime vSmartMOM.doubling!(pol_type, SFI, expk_, ndoubl, added_layer_GPU, I_static_, RadiativeTransfer.Architectures.GPU())
