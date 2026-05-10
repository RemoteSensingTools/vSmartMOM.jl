# =================================================================
# Shared test utilities for vSmartMOM Jacobian and forward model tests
# Include this file in test scripts: include("test_helpers.jl")
# =================================================================

using vSmartMOM, vSmartMOM.CoreRT
using Test
using Statistics: mean as stats_mean
using Distributions: Normal, LogNormal

"""
    run_lin_rt(params)

Run the linearized RT model and return reflectance, Jacobians, and model objects.
Returns: (R, dR, NAer, NGas, NSurf, model, lin_model)
"""
function run_lin_rt(params)
    model, lin_model = model_from_parameters(LinMode(), params)
    NAer = length(params.scattering_params.rt_aerosols)
    NGas = size(lin_model.τ̇_abs[1], 1)
    NSurf = 1
    R, T, dR, dT = rt_run(model, lin_model, NAer, NGas, NSurf)
    return R, dR, NAer, NGas, NSurf, model, lin_model
end

"""
    run_fwd_only(params)

Run the forward RT model (using the linearized code path for consistency).
Returns: R (reflectance only)
"""
function run_fwd_only(params)
    R, _, _, _, _, _, _ = run_lin_rt(params)
    return R
end

"""
    run_fwd_noRS(params; i_band=1)

Run the forward RT model without Raman scattering (pure elastic).
Returns: (R, model)
"""
function run_fwd_noRS(params; i_band=1)
    model = model_from_parameters(params)
    R = rt_run(model; i_band=i_band)
    return R, model
end

"""
    rel_errors(analytic, fd; threshold=1e-12)

Compute max and mean relative error between analytic and finite-difference Jacobians,
ignoring entries where |fd| < threshold to avoid division by near-zero.
Returns: (max_err, mean_err)
"""
function rel_errors(analytic, fd; threshold=1e-12)
    mask = abs.(fd) .> threshold
    if !any(mask)
        return (NaN, NaN)
    end
    rel = abs.((analytic[mask] .- fd[mask]) ./ fd[mask])
    return (maximum(rel), stats_mean(rel))
end

"""
    fd_jacobian_R(params, R_base, perturber!, get_delta)

Generic finite-difference Jacobian for TOA reflectance.
- `R_base`: baseline reflectance from `run_fwd_only(params)`
- `perturber!(tmp_params)`: modifies tmp_params in-place
- `get_delta(params, tmp_params)`: returns the scalar perturbation Δ

Returns: (K_FD, Δ)
"""
function fd_jacobian_R(params, R_base, perturber!, get_delta)
    tmp_params = deepcopy(params)
    perturber!(tmp_params)
    R_pert = run_fwd_only(tmp_params)
    Δ = get_delta(params, tmp_params)
    K_FD = (R_pert .- R_base) ./ Δ
    return K_FD, Δ
end
