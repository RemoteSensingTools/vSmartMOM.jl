#!/usr/bin/env julia
# ==========================================================================
# Jacobian tests on GPU (CUDA)
# ==========================================================================
# Runs the same linearized RT as test_jacobians_unit.jl but with
# architecture = GPU(). Checks that R and dR are finite and optionally
# that results match CPU within tolerance (GPU can have small float differences).
#
# Requires CUDA-capable device and CUDA.jl. Skips with @test_skip if no GPU.
# ==========================================================================

using Test
using vSmartMOM, vSmartMOM.CoreRT
using Distributions, Statistics

const YAML_FAST = "test_parameters/JacobianTestFast.yaml"

# Load CUDA at top level if available; then check device
CUDA_LOADED = false
try
    using CUDA
    global CUDA_LOADED = true
catch
end
can_use_gpu() = CUDA_LOADED ? CUDA.functional() : false

function run_lin_rt_gpu(params)
    model, lin_model = model_from_parameters(LinMode(), params)
    NAer = length(params.scattering_params.rt_aerosols)
    NGas = size(lin_model.τ̇_abs[1], 1)
    NSurf = 1
    R, T, dR, dT = rt_run(model, lin_model, NAer, NGas, NSurf)
    return R, dR, NAer, NGas, NSurf
end

# =====================================================================
@testset "Jacobian GPU" begin
# =====================================================================

    has_gpu = can_use_gpu()
    if !has_gpu
        @test_skip "CUDA not available or not functional; skipping GPU Jacobian tests"
        println("  Skipping GPU tests (no CUDA or device).")
    else

    # Load params and switch to GPU
    params_cpu = parameters_from_yaml(YAML_FAST)
    params_gpu = parameters_from_yaml(YAML_FAST)
    params_gpu.architecture = vSmartMOM.GPU()

    println("Running linearized RT on GPU...")
    try
        R_gpu, dR_gpu, NAer, NGas, NSurf = run_lin_rt_gpu(params_gpu)
        Nparams = CoreRT.n_total(CoreRT.ParameterLayout(
            n_aerosols=NAer, n_gases=NGas, n_surface=NSurf))
        R_gpu_a = Array(R_gpu)
        dR_gpu_a = Array(dR_gpu)

        @testset "GPU outputs finite" begin
            @test all(isfinite.(R_gpu_a))
            @test all(isfinite.(dR_gpu_a))
        end

        @testset "GPU vs CPU consistency (same config)" begin
            R_cpu, dR_cpu, _, _, _ = run_lin_rt_gpu(params_cpu)
            R_cpu_a = Array(R_cpu)
            dR_cpu_a = Array(dR_cpu)
            rtol = 1e-4
            @test all(isapprox.(R_gpu_a, R_cpu_a; rtol=rtol, nans=true))
            @test all(isapprox.(dR_gpu_a, dR_cpu_a; rtol=rtol, nans=true))
            println("  GPU vs CPU: R and dR match within rtol=$rtol")
        end
        println("  ✓ GPU Jacobian tests passed.")
    catch e
        # Functional hardware must expose compilation, scalar-indexing, and
        # numerical failures to the runner. Only absent hardware is skipped.
        rethrow()
    end
    end  # has_gpu
end

println("\nGPU Jacobian tests complete.")
