#!/usr/bin/env julia
# ==========================================================================
# Hybrid AD Tests: CPU/GPU Linearized RT + ForwardDiff Mie Consistency
# ==========================================================================

using Test
using vSmartMOM, vSmartMOM.CoreRT, vSmartMOM.Scattering
using Distributions, Statistics, LinearAlgebra

const YAML_FAST = "test_parameters/JacobianTestFast.yaml"

CUDA_OK = false
try
    using CUDA
    global CUDA_OK = CUDA.functional()
catch
end

# ---------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------
function run_lin(params)
    model, lin_model = model_from_parameters(LinMode(), params)
    NAer  = length(params.scattering_params.rt_aerosols)
    NGas  = size(lin_model.τ̇_abs[1], 1)
    NSurf = 1
    R, T, dR, dT = rt_run(model, lin_model, NAer, NGas, NSurf)
    return R, dR, model, lin_model, NAer, NGas, NSurf
end

function run_fwd(params)
    R, _, _, _, _, _, _ = run_lin(params)
    return R
end

function rel_errors(analytic, fd; threshold=1e-12)
    mask = abs.(fd) .> threshold
    if any(mask)
        errs = abs.(analytic[mask] .- fd[mask]) ./ abs.(fd[mask])
        return (max=maximum(errs), mean=mean(errs))
    else
        return (max=maximum(abs.(analytic .- fd)), mean=mean(abs.(analytic .- fd)))
    end
end

# =====================================================================
@testset "Hybrid AD" begin
# =====================================================================

# ------------------------------------------------------------------
@testset "CPU Linearized RT" begin
    params = parameters_from_yaml(YAML_FAST)
    params.architecture = vSmartMOM.Architectures.CPU()
    R, dR, model, lin_model, NAer, NGas, NSurf = run_lin(params)

    @test all(isfinite.(R))
    @test all(isfinite.(dR))

    layout = CoreRT.ParameterLayout(aerosol_params=7, n_aerosols=NAer,
                                     n_gases=NGas, n_surface=NSurf)

    # Surface albedo: finite-difference check
    idx_alb = CoreRT.surface_index(layout, 1)
    eps_alb = 1e-4
    params_p = parameters_from_yaml(YAML_FAST)
    params_p.architecture = vSmartMOM.Architectures.CPU()
    brdf = params_p.brdf[1]
    params_p.brdf[1] =
        LambertianSurfaceScalar{Float64}(brdf.albedo + eps_alb)
    R_p = run_fwd(params_p)
    dR_fd_alb = (R_p .- R) ./ eps_alb

    errs = rel_errors(dR[idx_alb, :, :, :], dR_fd_alb)
    @test errs.max < 0.05
    println("  CPU albedo Jacobian: max rel err = $(round(errs.max, sigdigits=3))")
end

# ------------------------------------------------------------------
if CUDA_OK
    @testset "GPU Linearized RT" begin
        params_gpu = parameters_from_yaml(YAML_FAST)
        params_gpu.architecture = vSmartMOM.Architectures.GPU()
        R_gpu, dR_gpu, _, _, NAer, NGas, NSurf = run_lin(params_gpu)

        R_gpu_a  = Array(R_gpu)
        dR_gpu_a = Array(dR_gpu)

        @test all(isfinite.(R_gpu_a))
        @test all(isfinite.(dR_gpu_a))

        # Compare with CPU
        params_cpu = parameters_from_yaml(YAML_FAST)
        params_cpu.architecture = vSmartMOM.Architectures.CPU()
        R_cpu, dR_cpu, _, _, _, _, _ = run_lin(params_cpu)

        @test all(isapprox.(R_gpu_a, R_cpu; rtol=1e-4, nans=true))
        @test all(isapprox.(dR_gpu_a, dR_cpu; rtol=1e-3, nans=true))
        println("  GPU vs CPU: R and dR match.")
    end

    @testset "GPU Forward noRS" begin
        params_gpu = parameters_from_yaml(YAML_FAST)
        params_gpu.architecture = vSmartMOM.Architectures.GPU()
        R_gpu = run_fwd(params_gpu)

        params_cpu = parameters_from_yaml(YAML_FAST)
        params_cpu.architecture = vSmartMOM.Architectures.CPU()
        R_cpu = run_fwd(params_cpu)

        R_gpu_a = Array(R_gpu)
        @test all(isfinite.(R_gpu_a))
        @test all(isapprox.(R_gpu_a, R_cpu; rtol=1e-4, nans=true))
        println("  GPU Forward: matches CPU.")
    end
else
    @test_skip "CUDA not available - GPU tests skipped"
end

# ------------------------------------------------------------------
@testset "ForwardDiff Mie Consistency" begin
    params = parameters_from_yaml(YAML_FAST)
    rt_aer = params.scattering_params.rt_aerosols[1]
    truncation_type = Scattering.δBGE{Float64}(params.l_trunc, params.Δ_angle)
    mie_model = make_mie_model(params.scattering_params.decomp_type,
                                rt_aer.aerosol, params.scattering_params.λ_ref,
                                params.polarization_type,
                                truncation_type,
                                params.scattering_params.r_max,
                                params.scattering_params.nquad_radius)

    aer_fwd, lin_aer = compute_aerosol_optical_properties(LinMode(), mie_model, Float64)

    aer_ad = compute_aerosol_optical_properties(mie_model; autodiff=true)

    @test isapprox(aer_fwd.ω̃, aer_ad.ω̃; rtol=1e-6)
    @test isapprox(aer_fwd.k, aer_ad.k; rtol=1e-6)

    @test isapprox(aer_fwd.greek_coefs.α, aer_ad.greek_coefs.α; rtol=1e-5)
    @test isapprox(aer_fwd.greek_coefs.β, aer_ad.greek_coefs.β; rtol=1e-5)

    println("  ForwardDiff vs analytic Mie: ω̃ and k match.")
end

# =====================================================================
end # Hybrid AD
println("\nHybrid AD tests complete.")
