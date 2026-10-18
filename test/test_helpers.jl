# =================================================================
# Shared support for vSmartMOM Jacobian and forward-model tests.
# Support file, not a test suite: runtests.jl does not include it directly.
#
# Load it with the include guard below, which keeps one copy of the module
# when several test files are included into the same process:
#
#     isdefined(@__MODULE__, :VSmartMOMTestSupport) || include("test_helpers.jl")
#     using .VSmartMOMTestSupport
#
# Contracts:
#   run_linearized_rt  linearized solver (`LinMode()`); radiances and Jacobians
#   linearized_R       TOA radiance from the linearized solver only
#   run_forward_rt     ordinary forward solver (`model_from_parameters(params)`)
#   forward_R          TOA radiance from the ordinary forward solver only
#   rel_errors         relative errors where |fd| > threshold, absolute check
#                      elsewhere; NaN on non-finite input, Inf on a near-zero
#                      mismatch
# =================================================================

module VSmartMOMTestSupport

using vSmartMOM, vSmartMOM.CoreRT
using Statistics: mean, median

export run_linearized_rt, linearized_R, run_forward_rt, forward_R,
       rel_errors, fd_jacobian_R

"""
    run_linearized_rt(params)

Build `model_from_parameters(LinMode(), params)` and run the linearized solver
`rt_run(model, lin_model, NAer, NGas, NSurf)` for band 1.

Returns the named tuple `(R, dR, NAer, NGas, NSurf, model, lin_model, T, dT)`.
Positional destructuring of its first seven fields matches the historical
`run_lin_rt` order.
"""
function run_linearized_rt(params)
    model, lin_model = model_from_parameters(LinMode(), params)
    NAer = isnothing(params.scattering_params) ? 0 :
           length(params.scattering_params.rt_aerosols)
    NGas = size(lin_model.τ̇_abs[1], 1)
    NSurf = 1
    R, T, dR, dT = rt_run(model, lin_model, NAer, NGas, NSurf)
    return (R=R, dR=dR, NAer=NAer, NGas=NGas, NSurf=NSurf,
            model=model, lin_model=lin_model, T=T, dT=dT)
end

"""
    linearized_R(params)

TOA radiance from the *linearized* solver (`run_linearized_rt(params).R`).
Use it for finite differences around a linearized base state; it is not an
independent ordinary-forward check (use [`forward_R`](@ref) for that).
"""
linearized_R(params) = run_linearized_rt(params).R

"""
    run_forward_rt(params; i_band=1)

Build the ordinary forward model (`model_from_parameters(params)`) and run
`rt_run(model; i_band)` (pure elastic, no Raman).

Returns `(R, T, result, model)`: `R` and `T` are the TOA-upwelling and
BOA-downwelling radiance arrays `[nVZA × nStokes × nSpec]`, and `result` is
the full `ObserverRTResult`.
"""
function run_forward_rt(params; i_band::Integer=1)
    model = model_from_parameters(params)
    result = rt_run(model; i_band=i_band)
    R, T = result
    return (R=R, T=T, result=result, model=model)
end

"""
    forward_R(params; i_band=1)

TOA radiance array from the ordinary forward solver (`run_forward_rt(params).R`).
"""
forward_R(params; i_band::Integer=1) = run_forward_rt(params; i_band=i_band).R

"""
    rel_errors(analytic, fd; threshold=1e-12, atol=threshold)

Relative errors `|analytic - fd| / |fd|` over the entries with
`|fd| > threshold`. Entries with `|fd| ≤ threshold` are checked on absolute
error instead: each must satisfy `|analytic - fd| ≤ atol`.

Returns `(max, mean, median, n, max_abs_masked)`: relative-error statistics,
the number `n` of entries compared relatively, and the largest absolute error
among the near-zero entries (`0.0` if there are none).
- Non-finite input in either array: every statistic is `NaN` (`n == 0`).
- A near-zero entry with `|analytic - fd| > atol`: `max`, `mean` and `median`
  are `Inf`, so a large analytic value against `fd ≈ 0` cannot pass.
- All entries near zero and within `atol` (a legitimately all-zero
  Jacobian): `max`, `mean` and `median` are `0.0` with `n == 0`.
Any `err.max < tol` assertion therefore fails on NaN, Inf, or a masked
mismatch; matching NaNs can never pass an accuracy check. Positional
destructuring `max_e, mean_e = rel_errors(...)` takes `(max, mean)`.
"""
function rel_errors(analytic, fd; threshold=1e-12, atol=threshold)
    size(analytic) == size(fd) ||
        throw(DimensionMismatch("analytic $(size(analytic)) vs fd $(size(fd))"))
    (all(isfinite, analytic) && all(isfinite, fd)) ||
        return (max=NaN, mean=NaN, median=NaN, n=0, max_abs_masked=NaN)
    mask = abs.(fd) .> threshold
    max_abs_masked = maximum(abs.(analytic[.!mask] .- fd[.!mask]); init=0.0)
    n = count(mask)
    max_abs_masked > atol &&
        return (max=Inf, mean=Inf, median=Inf, n=n, max_abs_masked=max_abs_masked)
    n == 0 && return (max=0.0, mean=0.0, median=0.0, n=0, max_abs_masked=max_abs_masked)
    errs = abs.(analytic[mask] .- fd[mask]) ./ abs.(fd[mask])
    return (max=maximum(errs), mean=mean(errs), median=median(errs), n=n,
            max_abs_masked=max_abs_masked)
end

"""
    fd_jacobian_R(params, R_base, perturber!, get_delta)

One-sided finite-difference Jacobian of the TOA radiance.
- `R_base`: base radiance from `linearized_R(params)`
- `perturber!(tmp_params)`: modifies `tmp_params` in place
- `get_delta(params, tmp_params)`: returns the scalar perturbation Δ

The perturbed radiance also comes from `linearized_R`, so base and perturbed
states use the same solver. Returns `(K_FD, Δ)`.
"""
function fd_jacobian_R(params, R_base, perturber!, get_delta)
    tmp_params = deepcopy(params)
    perturber!(tmp_params)
    R_pert = linearized_R(tmp_params)
    Δ = get_delta(params, tmp_params)
    K_FD = (R_pert .- R_base) ./ Δ
    return K_FD, Δ
end

end # module VSmartMOMTestSupport
