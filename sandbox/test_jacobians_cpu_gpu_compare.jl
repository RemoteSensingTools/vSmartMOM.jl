#!/usr/bin/env julia
# ==========================================================================
# CPU vs GPU Jacobian comparison test
# ==========================================================================
# Runs linearized RT twice (CPU and GPU) on the same configuration and checks
# that reflectance Jacobians are numerically close.
# ==========================================================================

using Test
using vSmartMOM, vSmartMOM.CoreRT
using Statistics

const YAML_FAST = "test/test_parameters/JacobianTestFast.yaml"
const RTOL_R = 1e-4
const RTOL_dR = 1e-4
const ATOL_R = 1e-9
const ATOL_dR = 1e-9

CUDA_LOADED = false
try
    using CUDA
    global CUDA_LOADED = true
catch
end

can_use_gpu() = CUDA_LOADED ? CUDA.functional() : false

function run_lin_rt(params)
    model, lin_model = model_from_parameters(LinMode(), params)
    NAer = length(params.scattering_params.rt_aerosols)
    NGas = size(lin_model.τ̇_abs[1], 1)
    NSurf = 1
    R, _, dR, _ = rt_run(model, lin_model, NAer, NGas, NSurf)
    return R, dR
end

function rel_metrics(reference, candidate; floor=1e-12)
    diff = candidate .- reference
    abs_diff = abs.(diff)
    abs_max = maximum(abs_diff)
    abs_mean = mean(abs_diff)

    mask = abs.(reference) .> floor
    if any(mask)
        rel = abs_diff[mask] ./ abs.(reference[mask])
        return (abs_max=abs_max, abs_mean=abs_mean, rel_max=maximum(rel), rel_mean=mean(rel))
    end

    return (abs_max=abs_max, abs_mean=abs_mean, rel_max=0.0, rel_mean=0.0)
end

@testset "Jacobian CPU vs GPU comparison" begin
    if !can_use_gpu()
        @test_skip "CUDA not available or not functional; skipping CPU/GPU Jacobian comparison"
        println("Skipping comparison test (no functional CUDA).")
    else
        params_cpu = parameters_from_yaml(YAML_FAST)
        params_cpu.architecture = vSmartMOM.CPU()
        params_gpu = parameters_from_yaml(YAML_FAST)
        params_gpu.architecture = vSmartMOM.GPU()

        println("Running linearized RT on CPU...")
        R_cpu, dR_cpu = run_lin_rt(params_cpu)
        R_cpu_a = Array(R_cpu)
        dR_cpu_a = Array(dR_cpu)

        println("Running linearized RT on GPU...")
        R_gpu, dR_gpu = run_lin_rt(params_gpu)
        R_gpu_a = Array(R_gpu)
        dR_gpu_a = Array(dR_gpu)

        @test size(R_cpu_a) == size(R_gpu_a)
        @test size(dR_cpu_a) == size(dR_gpu_a)
        @test all(isfinite.(R_cpu_a))
        @test all(isfinite.(R_gpu_a))
        @test all(isfinite.(dR_cpu_a))
        @test all(isfinite.(dR_gpu_a))

        r_metrics = rel_metrics(R_cpu_a, R_gpu_a)
        dr_metrics = rel_metrics(dR_cpu_a, dR_gpu_a)
        println("R  CPU vs GPU: max_rel=$(r_metrics.rel_max), mean_rel=$(r_metrics.rel_mean), max_abs=$(r_metrics.abs_max)")
        println("dR CPU vs GPU: max_rel=$(dr_metrics.rel_max), mean_rel=$(dr_metrics.rel_mean), max_abs=$(dr_metrics.abs_max)")

        @test all(isapprox.(R_gpu_a, R_cpu_a; rtol=RTOL_R, atol=ATOL_R, nans=true))
        @test all(isapprox.(dR_gpu_a, dR_cpu_a; rtol=RTOL_dR, atol=ATOL_dR, nans=true))
    end
end

println("\nCPU/GPU Jacobian comparison complete.")
