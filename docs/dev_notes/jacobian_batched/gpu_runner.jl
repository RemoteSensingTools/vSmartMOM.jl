# Audit harness: execute the shipped GPU tests with absolute include paths.
# The candidate's runner repeats local/gpu in its include paths.
using vSmartMOM, vSmartMOM.Architectures, vSmartMOM.Scattering
using vSmartMOM.CoreRT, vSmartMOM.InelasticScattering
using Test, Printf, Logging, LinearAlgebra, Statistics, Distributions, CUDA
CUDA.functional() || error("Audit requires a functional CUDA GPU")
println("AUDIT GPU: ", CUDA.name(CUDA.device()))
BLAS.set_num_threads(1)
cd(joinpath(pkgdir(vSmartMOM), "test"))
@testset "GPU tests with corrected include paths" begin
    for filename in ["test_mie_gpu.jl", "test_raman_fused_kernels.jl",
                     "test_multisensor_heights_gpu.jl", "test_forward_raman_gpu.jl",
                     "test_jacobians_GPU.jl"]
        @testset "$filename" begin
            println("AUDIT_GPU_BEGIN ", filename); flush(stdout)
            include(joinpath(pwd(), "local", "gpu", filename))
        end
    end
end
