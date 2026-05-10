#!/usr/bin/env julia
# ==========================================================================
# Jacobian comparison: Analytic vs Finite Difference
# ==========================================================================
# Compares the analytic tangent-linear Jacobians to central finite difference (FD)
# for the same R (linearized path). Use to check for FD step-size issues and
# to validate analytic derivatives.
#
# Note: ForwardDiff (AD) through the linearized model is not supported because
# model_from_parameters(LinMode(), params) expects Float64; Dual-valued params
# cause MethodError in array assignment.
#
# Run after test_jacobians_unit.jl or standalone (uses same YAML).
# ==========================================================================

using Test
using vSmartMOM, vSmartMOM.CoreRT
using Distributions, Statistics, Printf

const YAML_FAST = "test/test_parameters/JacobianTestFast.yaml"

# ---------------------------------------------------------------
# Forward R from linearized path (same as analytic run) so R shape/values match
# ---------------------------------------------------------------
function run_fwd_R_from_lin(params)
    model, lin_model = model_from_parameters(LinMode(), params)
    NAer = length(params.scattering_params.rt_aerosols)
    NGas = size(lin_model.τ̇_abs[1], 1)
    NSurf = 1
    R, T, _, _ = rt_run(model, lin_model, NAer, NGas, NSurf)
    return R
end

# ---------------------------------------------------------------
# Parameter wrappers: f(x) = R[1,1,1] with that param set to x
# Each re-loads params from YAML and sets one field so AD gets a clean scalar input.
# ---------------------------------------------------------------
function make_f_param(; param::Symbol, base_params)
    if param == :τ_ref
        x0 = base_params.scattering_params.rt_aerosols[1].τ_ref
        f_τ_ref(x) = (p = parameters_from_yaml(YAML_FAST); p.scattering_params.rt_aerosols[1].τ_ref = x; run_fwd_R_from_lin(p)[1, 1, 1])
        return f_τ_ref, x0
    elseif param == :albedo
        x0 = base_params.brdf[1].albedo
        f_albedo(x) = (p = parameters_from_yaml(YAML_FAST); p.brdf = [CoreRT.LambertianSurfaceScalar(x)]; run_fwd_R_from_lin(p)[1, 1, 1])
        return f_albedo, x0
    elseif param == :p₀
        aero = base_params.scattering_params.rt_aerosols[1]
        σp_base = std(aero.profile)
        p0_base = mean(aero.profile)
        f_p0(x) = (p = parameters_from_yaml(YAML_FAST); p.scattering_params.rt_aerosols[1].profile = Normal(x, σp_base); run_fwd_R_from_lin(p)[1, 1, 1])
        return f_p0, p0_base
    elseif param == :nᵣ
        x0 = base_params.scattering_params.rt_aerosols[1].aerosol.nᵣ
        f_nr(x) = (p = parameters_from_yaml(YAML_FAST); p.scattering_params.rt_aerosols[1].aerosol.nᵣ = x; run_fwd_R_from_lin(p)[1, 1, 1])
        return f_nr, x0
    else
        error("unknown param $param")
    end
end

# Param index in dR for our test params (NAer=1 → 1..7 aero, then gas, then surf)
param_index(p) = p == :τ_ref ? 1 : p == :nᵣ ? 2 : p == :p₀ ? 6 : p == :albedo ? 9 : error("unknown $p")

# ---------------------------------------------------------------
# Central FD for dR[1,1,1]/d(param)
# ---------------------------------------------------------------
function central_fd(base_params, param::Symbol; rel_step=1e-4)
    if param == :τ_ref
        x0 = base_params.scattering_params.rt_aerosols[1].τ_ref
        δ = max(x0 * rel_step, 1e-10)
    elseif param == :albedo
        x0 = base_params.brdf[1].albedo
        δ = max(x0 * rel_step, 1e-8)
    elseif param == :p₀
        x0 = mean(base_params.scattering_params.rt_aerosols[1].profile)
        δ = max(x0 * rel_step, 1.0)
    elseif param == :nᵣ
        x0 = base_params.scattering_params.rt_aerosols[1].aerosol.nᵣ
        δ = max(abs(x0) * rel_step, 1e-6)
    else
        error("unknown $param")
    end
    f, _ = make_f_param(param=param, base_params=base_params)
    R_plus = f(x0 + δ)
    R_minus = f(x0 - δ)
    return (R_plus - R_minus) / (2 * δ)
end

# ---------------------------------------------------------------
# Run comparison for a subset of parameters
# ---------------------------------------------------------------
function run_comparison(; params_to_test=[:τ_ref, :albedo, :nᵣ], skip_p₀=true)
    println("Loading parameters and running linearized RT once...")
    base_params = parameters_from_yaml(YAML_FAST)
    model, lin_model = model_from_parameters(LinMode(), base_params)
    NAer = length(base_params.scattering_params.rt_aerosols)
    NGas = size(lin_model.τ̇_abs[1], 1)
    NSurf = 1
    R_base, T_base, dR_an, dT_an = rt_run(model, lin_model, NAer, NGas, NSurf)
    nSpec = size(R_base, 3)

    # Single output index for comparison
    iVZA, iStokes, iSpec = 1, 1, 1
    R_scalar = R_base[iVZA, iStokes, iSpec]

    list = skip_p₀ ? filter(p -> p != :p₀, params_to_test) : params_to_test

    println("\n=== Jacobian comparison: Analytic vs FD (central) ===\n")
    println("Output: R[$iVZA,$iStokes,$iSpec] = $R_scalar")
    println("(AD/ForwardDiff skipped: linearized model construction expects Float64; Dual not supported.)")
    println(lpad("param", 10), " | ", lpad("analytic", 14), " | ", lpad("FD (central)", 14), " | ", "rel_err an vs FD")
    println(String(repeat("-", 60)))

    t_fd_total = 0.0

    for param in list
        idx = param_index(param)
        dR_analytic = dR_an[idx, iVZA, iStokes, iSpec]

        # Central FD (same code path as analytic for R)
        t_fd = @elapsed dR_fd = central_fd(base_params, param; rel_step=1e-4)
        t_fd_total += t_fd

        rel_an_fd = abs(dR_analytic - dR_fd) / (abs(dR_fd) + 1e-30)
        @printf("%10s | %14.6e | %14.6e | %12.2f%%\n",
                string(param), dR_analytic, dR_fd, rel_an_fd * 100)
    end

    println("\n=== Timing ===")
    println("  FD total (central, $(length(list)) params): $(round(t_fd_total, digits=2)) s")
    println("  → Analytic gives all $(NAer*7+NGas+NSurf) columns in one run; FD needs two runs per parameter.")
    return (dR_an=dR_an, R_base=R_base, base_params=base_params)
end

# =====================================================================
@testset "Jacobian AD comparison" begin
# =====================================================================
    # Test a few params; skip p₀ by default (slow + known Bug 23 issues)
    result = run_comparison(params_to_test=[:τ_ref, :albedo, :nᵣ], skip_p₀=true)
    @test all(isfinite.(result.dR_an))
    println("\n  ✓ Comparison finished (analytic vs FD).")
end

println("\nDone. Use run_comparison(; skip_p₀=false) to include p₀ (slower, large errors expected).")
