# GPU / Metal test suite — NOT part of the CI unit suite (test/runtests.jl).
#
# Run locally with:
#     julia --project=test test/local/gpu/runtests.jl
# (or as part of the full local suite: test/local/runtests.jl)
#
# - test_mie_gpu.jl runs on the KernelAbstractions CPU backend (no GPU required)
#   and validates the DoubleSingle/Neumaier precision layer + NAI2 pipeline.
# - test_forward_raman_gpu.jl and test_jacobians_GPU.jl need a functional CUDA
#   device; they are guarded and skip cleanly when CUDA is unavailable.
#
# We run from the test/ root (two levels up) so the GPU test files' relative
# paths resolve: CI configs at test_parameters/, local-data configs (e.g.
# O2Parameters_GPU) at local/test_parameters/.

using vSmartMOM
using vSmartMOM.Architectures
using vSmartMOM.Scattering
using vSmartMOM.CoreRT
using vSmartMOM.InelasticScattering
using Test
using Printf
using Logging
using LinearAlgebra
using Statistics
using Distributions

const _GPU_ORIG_CWD = pwd()
cd(normpath(joinpath(@__DIR__, "..", "..")))
try
    @testset "GPU / Metal" begin
        # GPU Mie kernels via the KA CPU backend — no CUDA required.
        @testset "Mie GPU kernels" begin include("local/gpu/test_mie_gpu.jl") end

        CUDA_AVAILABLE = try
            using CUDA
            CUDA.functional()
        catch
            false
        end

        if CUDA_AVAILABLE
            @testset "Raman GPU"     begin include("local/gpu/test_forward_raman_gpu.jl") end
            @testset "Jacobians GPU" begin include("local/gpu/test_jacobians_GPU.jl")     end
        else
            @info "CUDA not available — skipping CUDA-only GPU tests (Raman GPU, Jacobians GPU)."
        end
    end
finally
    cd(_GPU_ORIG_CWD)
end
