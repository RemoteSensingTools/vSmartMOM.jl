module OptimalEstimation

using LinearAlgebra

export OESettings,
       ForwardEvaluation,
       IterationRecord,
       OEResult,
       solve_optimal_estimation,
       band_reduced_chi_squared,
       posterior_products

"""Configurable Connor/ACOS-style optimal-estimation controls."""
Base.@kwdef struct OESettings
    convergence_threshold::Float64 = 2.0
    maximum_iterations::Int = 7
    maximum_divergences::Int = 2
    maximum_band_chi_squared::Float64 = 1.4
    initial_gamma::Float64 = 10.0
    poor_reduction_ratio::Float64 = 0.25
    good_reduction_ratio::Float64 = 0.75
    rejection_ratio::Float64 = 1.0e-4
    minimum_gamma::Float64 = 1.0e-8
end

"""
One forward-model evaluation in measurement space.

`jacobian` follows the OE convention `(measurement, active_parameter)`.
`band_ranges` are one-based inclusive ranges into `measurement`. The optional
timing tuple is carried into the iteration record without imposing a particular
forward-model decomposition.
"""
struct ForwardEvaluation
    measurement::Vector{Float64}
    jacobian::Matrix{Float64}
    band_ranges::Vector{UnitRange{Int}}
    timing::NamedTuple
end

function ForwardEvaluation(measurement::AbstractVector,
                           jacobian::AbstractMatrix,
                           band_ranges::AbstractVector{<:UnitRange};
                           timing=(forward_seconds=NaN,
                                   jacobian_seconds=NaN,
                                   instrument_seconds=NaN))
    y = Float64.(measurement)
    K = Float64.(jacobian)
    ranges = UnitRange{Int}[Int(first(r)):Int(last(r)) for r in band_ranges]
    size(K, 1) == length(y) || throw(DimensionMismatch(
        "Jacobian has $(size(K, 1)) measurement rows but F(x) has $(length(y))"))
    isempty(ranges) && throw(ArgumentError("at least one spectral band is required"))
    reduce(vcat, collect.(ranges)) == collect(eachindex(y)) || throw(ArgumentError(
        "band ranges must partition the measurement vector in order"))
    return ForwardEvaluation(y, K, ranges, timing)
end

"""Diagnostics for one evaluated LM trial state."""
struct IterationRecord
    trial::Int
    iteration::Int
    state::Vector{Float64}
    proposed_step::Vector{Float64}
    gamma::Float64
    reduction_ratio::Float64
    accepted::Bool
    divergent::Bool
    measurement_cost::Float64
    prior_cost::Float64
    total_cost::Float64
    forecast_cost::Float64
    d_sigma_sq::Float64
    d_sigma_sq_scaled::Float64
    effective_parameter_count::Int
    band_chi_squared::Vector{Float64}
    evaluation_seconds::Float64
    linear_algebra_seconds::Float64
    model_timing::NamedTuple
end

"""Complete terminal retrieval result and the requested error-analysis matrices."""
struct OEResult
    final_state::Vector{Float64}
    final_measurement::Vector{Float64}
    final_jacobian::Matrix{Float64}
    gain_matrix::Matrix{Float64}
    posterior_covariance::Matrix{Float64}
    averaging_kernel::Matrix{Float64}
    final_band_chi_squared::Vector{Float64}
    records::Vector{IterationRecord}
    outcome::Int
    converged::Bool
    fit_quality_ok::Bool
    divergence_count::Int
    final_evaluation_seconds::Float64
end

function _validate_settings(settings::OESettings)
    settings.convergence_threshold > 0 || throw(ArgumentError(
        "convergence threshold must be positive"))
    settings.maximum_iterations >= 1 || throw(ArgumentError(
        "maximum_iterations must be at least one"))
    settings.maximum_divergences >= 1 || throw(ArgumentError(
        "maximum_divergences must be at least one"))
    settings.initial_gamma >= 0 || throw(ArgumentError(
        "initial_gamma must be non-negative"))
    0 <= settings.rejection_ratio < settings.poor_reduction_ratio <
        settings.good_reduction_ratio || throw(ArgumentError(
        "reduction-ratio thresholds must be ordered"))
    return nothing
end

"""Per-band measurement chi-square divided by the number of band samples."""
function band_reduced_chi_squared(residual::AbstractVector,
                                  variance::AbstractVector,
                                  band_ranges::AbstractVector{<:UnitRange})
    length(residual) == length(variance) || throw(DimensionMismatch(
        "residual and variance lengths differ"))
    return [sum(abs2, @view(residual[r]) ./ sqrt.(@view(variance[r]))) / length(r)
            for r in band_ranges]
end

function _costs(x, xa, residual, inv_variance, Sa_inverse)
    measurement_cost = dot(residual, inv_variance .* residual)
    prior_delta = x - xa
    prior_cost = dot(prior_delta, Sa_inverse * prior_delta)
    return measurement_cost, prior_cost, measurement_cost + prior_cost
end

function _inversion_step(x, xa, evaluation, residual, inv_variance,
                         sigma_ap, Sa_scaled_inverse, gamma)
    K = evaluation.jacobian
    weighted_K = K .* reshape(sqrt.(inv_variance), :, 1)
    KtSeK = weighted_K' * weighted_K
    scaled_KtSeK = KtSeK .* (sigma_ap * sigma_ap')
    lhs = Symmetric((1 + gamma) .* Sa_scaled_inverse .+ scaled_KtSeK)
    rhs = sigma_ap .* (K' * (inv_variance .* (-residual))) .+
          Sa_scaled_inverse * ((xa - x) ./ sigma_ap)

    dx_scaled = cholesky(lhs) \ rhs
    dx = sigma_ap .* dx_scaled
    d_sigma_sq = dot(dx_scaled, rhs)
    sensitive = vec(maximum(abs.(KtSeK); dims=2)) .> 1.0e-20
    effective_count = count(sensitive)
    effective_count > 0 || throw(ArgumentError(
        "the measurement Jacobian has no active parameter sensitivity"))
    d_sigma_sq_scaled = d_sigma_sq / effective_count

    residual_forecast = residual + K * dx
    x_forecast = x + dx
    return (; dx, d_sigma_sq, d_sigma_sq_scaled, effective_count,
            KtSeK, residual_forecast, x_forecast)
end

function _forecast_cost(step, xa, inv_variance, Sa_inverse)
    measurement_cost, prior_cost, total_cost = _costs(
        step.x_forecast, xa, step.residual_forecast, inv_variance, Sa_inverse)
    return (; measurement_cost, prior_cost, total_cost)
end

function _reduction_ratio(actual_cost, base_cost, forecast_cost)
    actual_change = actual_cost - base_cost
    predicted_change = forecast_cost - base_cost
    actual_change >= 0 && return 0.0
    abs(predicted_change) < -actual_change * 1.0e-20 && return 1.0
    return actual_change / predicted_change
end

function _updated_gamma(gamma, ratio, iteration, settings)
    iteration <= 1 && return gamma
    if ratio < settings.poor_reduction_ratio
        return gamma > settings.minimum_gamma ? 10gamma : one(gamma)
    elseif ratio > settings.good_reduction_ratio
        return gamma / 2
    end
    return gamma
end

function _evaluate_timed(evaluate, x)
    value = nothing
    elapsed = @elapsed value = evaluate(x)
    value isa ForwardEvaluation || throw(ArgumentError(
        "forward callback must return ForwardEvaluation"))
    return value, elapsed
end

"""
    posterior_products(K, Se_diagonal, Sa)

Return the undamped OE posterior covariance, gain matrix, and averaging kernel
at a terminal state. Dense `Se` is never formed.
"""
function posterior_products(K::AbstractMatrix,
                            Se_diagonal::AbstractVector,
                            Sa::AbstractMatrix)
    nmeasurement, nstate = size(K)
    length(Se_diagonal) == nmeasurement || throw(DimensionMismatch(
        "Se diagonal and Jacobian measurement dimensions differ"))
    size(Sa) == (nstate, nstate) || throw(DimensionMismatch(
        "Sa and Jacobian state dimensions differ"))
    inv_variance = 1.0 ./ Float64.(Se_diagonal)
    Sa_inverse = inv(Symmetric(Float64.(Sa)))
    weighted_K = Float64.(K) .* reshape(sqrt.(inv_variance), :, 1)
    information = Symmetric(Sa_inverse + weighted_K' * weighted_K)
    factor = cholesky(information)
    posterior = Matrix(factor \ Matrix{Float64}(I, nstate, nstate))
    posterior = (posterior + posterior') / 2
    gain = posterior * (Float64.(K)' .* reshape(inv_variance, 1, :))
    averaging_kernel = gain * Float64.(K)
    return (; posterior, gain, averaging_kernel)
end

"""
    solve_optimal_estimation(evaluate, measurement, Se_diagonal, xa, Sa;
                             initial_state=xa, settings=OESettings())

Run a Connor/ACOS-style damped optimal-estimation retrieval. The callback is
evaluated only at trial states and must return the forward measurement and
analytic Jacobian after the complete instrument operator. Measurement noise
is frozen through the supplied diagonal `Se_diagonal`.

Outcome codes follow the public OCO convention: 1 converged/all bands pass,
2 converged/at least one band fails, 3 maximum iterations, and 4 maximum
divergent steps. Rejected trial states remain in `records` with
`accepted=false`.
"""
function solve_optimal_estimation(evaluate,
                                  measurement::AbstractVector,
                                  Se_diagonal::AbstractVector,
                                  xa::AbstractVector,
                                  Sa::AbstractMatrix;
                                  initial_state::AbstractVector=xa,
                                  settings::OESettings=OESettings(),
                                  record_callback::Function=record -> nothing)
    _validate_settings(settings)
    y = Float64.(measurement)
    variance = Float64.(Se_diagonal)
    xa64 = Float64.(xa)
    x = Float64.(initial_state)
    Sa64 = Float64.(Sa)
    nstate = length(xa64)
    length(x) == nstate || throw(DimensionMismatch(
        "initial state and prior state lengths differ"))
    size(Sa64) == (nstate, nstate) || throw(DimensionMismatch(
        "Sa must be square with one row per active parameter"))
    length(variance) == length(y) || throw(DimensionMismatch(
        "measurement and Se diagonal lengths differ"))
    all(isfinite, y) || throw(ArgumentError("measurement contains non-finite values"))
    all(x -> isfinite(x) && x > 0, variance) || throw(ArgumentError(
        "Se diagonal must contain finite positive variances"))
    isposdef(Symmetric(Sa64)) || throw(ArgumentError(
        "active-state prior covariance must be positive definite"))

    sigma_ap = sqrt.(diag(Sa64))
    all(>(0), sigma_ap) || throw(ArgumentError(
        "active-state prior variances must be positive"))
    Sa_inverse = inv(Symmetric(Sa64))
    Sa_scaled = Sa64 ./ (sigma_ap * sigma_ap')
    Sa_scaled_inverse = inv(Symmetric(Sa_scaled))
    inv_variance = 1.0 ./ variance

    records = IterationRecord[]
    gamma = Float64(settings.initial_gamma)
    divergence_count = 0
    accepted_iteration = 0
    trial = 0
    base = nothing
    pending = nothing
    terminal_evaluation = nothing
    terminal_evaluation_seconds = NaN
    terminal_state = copy(x)
    outcome = 3

    while true
        trial += 1
        evaluation, evaluation_seconds = _evaluate_timed(evaluate, x)
        length(evaluation.measurement) == length(y) || throw(DimensionMismatch(
            "forward callback measurement length changed during retrieval"))
        size(evaluation.jacobian, 2) == nstate || throw(DimensionMismatch(
            "forward callback returned $(size(evaluation.jacobian, 2)) Jacobian columns; " *
            "expected $nstate"))
        residual = evaluation.measurement - y
        measurement_cost, prior_cost, total_cost = _costs(
            x, xa64, residual, inv_variance, Sa_inverse)
        chi_squared = band_reduced_chi_squared(
            residual, variance, evaluation.band_ranges)

        gamma_for_step = gamma
        step = nothing
        linear_seconds = @elapsed step = _inversion_step(
            x, xa64, evaluation, residual, inv_variance,
            sigma_ap, Sa_scaled_inverse, gamma_for_step)
        forecast = _forecast_cost(step, xa64, inv_variance, Sa_inverse)

        ratio = NaN
        accepted = true
        divergent = false
        if base !== nothing
            ratio = _reduction_ratio(total_cost, base.total_cost,
                                     pending.forecast.total_cost)
            gamma = _updated_gamma(gamma, ratio, accepted_iteration + 1, settings)
            divergent = ratio <= settings.rejection_ratio
            accepted = !divergent
        end

        record = IterationRecord(
            trial, accepted_iteration + 1, copy(x), copy(step.dx),
            gamma_for_step, ratio, accepted, divergent,
            measurement_cost, prior_cost, total_cost, forecast.total_cost,
            step.d_sigma_sq, step.d_sigma_sq_scaled, step.effective_count,
            chi_squared, evaluation_seconds, linear_seconds, evaluation.timing)
        push!(records, record)
        record_callback(record)

        if divergent
            divergence_count += 1
            if divergence_count > settings.maximum_divergences
                outcome = 4
                terminal_state = copy(base.x)
                terminal_evaluation = base.evaluation
                terminal_evaluation_seconds = base.evaluation_seconds
                break
            end
            # Recompute the proposal from the last accepted state with the
            # newly increased damping, without paying another forward call.
            retry = _inversion_step(
                base.x, xa64, base.evaluation, base.residual, inv_variance,
                sigma_ap, Sa_scaled_inverse, gamma)
            retry_forecast = _forecast_cost(retry, xa64, inv_variance, Sa_inverse)
            pending = (; step=retry, forecast=retry_forecast)
            x = base.x + retry.dx
            continue
        end

        accepted_iteration += 1
        base = (; x=copy(x), evaluation, residual, total_cost,
                evaluation_seconds)
        pending = (; step, forecast)

        if step.d_sigma_sq_scaled < settings.convergence_threshold
            terminal_state = x + step.dx
            terminal_evaluation, terminal_evaluation_seconds =
                _evaluate_timed(evaluate, terminal_state)
            outcome = 1 # updated from terminal band chi-square below
            break
        elseif accepted_iteration >= settings.maximum_iterations
            terminal_state = copy(x)
            terminal_evaluation = evaluation
            terminal_evaluation_seconds = evaluation_seconds
            outcome = 3
            break
        end

        x = x + step.dx
    end

    final_residual = terminal_evaluation.measurement - y
    final_chi_squared = band_reduced_chi_squared(
        final_residual, variance, terminal_evaluation.band_ranges)
    converged = outcome == 1
    fit_quality_ok = all(<(settings.maximum_band_chi_squared), final_chi_squared)
    if converged && !fit_quality_ok
        outcome = 2
    end
    converged = outcome in (1, 2)

    products = posterior_products(
        terminal_evaluation.jacobian, variance, Sa64)
    return OEResult(
        terminal_state,
        terminal_evaluation.measurement,
        terminal_evaluation.jacobian,
        products.gain,
        products.posterior,
        products.averaging_kernel,
        final_chi_squared,
        records,
        outcome,
        converged,
        fit_quality_ok,
        divergence_count,
        terminal_evaluation_seconds,
    )
end

end # module OptimalEstimation
