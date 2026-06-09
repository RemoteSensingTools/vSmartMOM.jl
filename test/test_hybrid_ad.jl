#!/usr/bin/env julia
# ==========================================================================
# Hybrid AD Tests: CPU/GPU Linearized RT + ForwardDiff Mie Consistency
# ==========================================================================

using Test
using vSmartMOM, vSmartMOM.CoreRT, vSmartMOM.Scattering
using Distributions, Statistics, LinearAlgebra, ForwardDiff

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

    errs = rel_errors(dR[:, :, :, idx_alb], dR_fd_alb)
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
    # Bug B fix test: verify that compute_mie_ab! no longer throws MethodError
    # when called with ForwardDiff Dual types.
    # We test directly at the compute_mie_ab! level to avoid a pre-existing
    # DimensionMismatch bug in phase_function_autodiff.jl's JacobianResult
    # preallocator (it uses gauleg/r_max directly while compute_NAI2 uses
    # gauleg_log with r_min from quantile, giving a different n_max).

    params = parameters_from_yaml(YAML_FAST)
    rt_aer = params.scattering_params.rt_aerosols[1]
    (; size_distribution, nᵣ, nᵢ) = rt_aer.aerosol

    # Test that Dual numbers flow through compute_mie_ab! without MethodError.
    # Use a fixed size_param and pre-compute the integer buffer sizes in Float64
    # (avoiding round(Int, Dual) issues) by fixing nᵣ/nᵢ for sizing only.
    fixed_size_param = 100.0   # representative
    n_max_fixed = Scattering.get_n_max(fixed_size_param)
    nmx_fixed = round(Int, max(n_max_fixed, fixed_size_param * max(nᵣ, nᵢ)) + 51)

    function mie_C_sca(x)
        # x = [nᵣ, nᵢ]; sizes are fixed integers from outer scope
        m_ref = x[1] - x[2] * im
        an = zeros(Complex{eltype(x)}, n_max_fixed)
        bn = zeros(Complex{eltype(x)}, n_max_fixed)
        Dn = zeros(Complex{eltype(x)}, nmx_fixed)
        Scattering.compute_mie_ab!(fixed_size_param, m_ref, an, bn, Dn)
        k = 2π / 0.77
        n_ = eltype(x).(2 .* (1:n_max_fixed) .+ 1)
        return [2π / k^2 * sum(n_ .* (abs2.(an) .+ abs2.(bn)));  # C_sca
                2π / k^2 * sum(n_ .* real.(an .+ bn))]           # C_ext
    end

    x0 = [nᵣ, nᵢ]
    # Should not throw MethodError after Bug B fix:
    J = ForwardDiff.jacobian(mie_C_sca, x0)
    @test size(J) == (2, 2)
    @test all(isfinite.(J))
    println("  ForwardDiff Mie ab: Jacobian computed successfully, shape=$(size(J))")

    # Compare ForwardDiff ∂C_sca/∂nᵣ against central finite differences
    eps_nr = 1e-6
    Cp = mie_C_sca([nᵣ + eps_nr, nᵢ])
    Cm = mie_C_sca([nᵣ - eps_nr, nᵢ])
    dCsca_fd = (Cp[1] - Cm[1]) / (2 * eps_nr)
    dCsca_ad = J[1, 1]
    errs_ad = abs(dCsca_ad - dCsca_fd) / (abs(dCsca_fd) + 1e-30)
    @test errs_ad < 1e-4
    println("  ∂C_sca/∂nᵣ: AD=$(round(dCsca_ad, sigdigits=5)), FD=$(round(dCsca_fd, sigdigits=5)), rel_err=$(round(errs_ad, sigdigits=3))")

    # NOTE: compute_aerosol_optical_properties(model; autodiff=true) has a
    # pre-existing DimensionMismatch bug in phase_function_autodiff.jl:
    # JacobianResult is preallocated using get_n_max(2π*r_max/λ), but
    # compute_NAI2 uses gauleg_log with r_min=max(quantile(...,1e-8),1e-6*r_max),
    # so maximum(x_size_param) < 2π*r_max/λ when run in Float64 but the Dual
    # path may give a different n_max. Fix requires editing
    # src/Scattering/phase_function_autodiff.jl (outside the scope of Bug B).
end

# =====================================================================
end # Hybrid AD
println("\nHybrid AD tests complete.")
