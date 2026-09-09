#!/usr/bin/env julia

"""
Validate the analytical surface-pressure column for a saved Round-5 migration
benchmark using a central finite difference of the complete instrument-space
forward model.

This is intentionally a real three-band retrieval-state check, not a synthetic
kernel fixture. `ROUND5_FD_STEP_HPA` defaults to 2 hPa. The smaller 0.5 hPa
step used by isolated Float32 pressure tests is not large enough for stable
subtraction of the complete instrument-space radiance in this scene; its
0.02 hPa real-spectroscopy probe is Float64 and is still less transferable.
"""

include(joinpath(@__DIR__, "benchmark_round5_migration_case.jl"))

using LinearAlgebra
using NCDatasets
using Printf

function fd_relative_l2(delta, reference)
    denominator = norm(reference)
    return iszero(denominator) ? norm(delta) : norm(delta) / denominator
end

function main_pressure_fd()
    baseline_path = required_path("ROUND5_FD_BASELINE")
    noise_path = required_path("ROUND5_BENCH_NOISE")
    truth_table = required_path("ROUND5_TRUTH_TABLE")
    prior_off = required_path("ROUND5_PRIOR_OFF")
    prior_on = required_path("ROUND5_PRIOR_ON")
    coefficient_path = required_path("ROUND5_STOKES_COEFFICIENTS")
    component_path = required_path("ROUND5_SCENE_COMPONENTS")
    output_path = abspath(get(ENV, "ROUND5_FD_OUTPUT", ""))
    isempty(output_path) && error("ROUND5_FD_OUTPUT is required")
    step_hpa = parse(Float64, get(ENV, "ROUND5_FD_STEP_HPA", "2.0"))
    step_hpa > 0 || error("ROUND5_FD_STEP_HPA must be positive")

    baseline = load(baseline_path)
    implementation = String(baseline["implementation"])
    assert_implementation(
        joinpath(INVERSION_ROOT, "VSmartMOMForward.jl"), implementation)
    state_index = Int(baseline["state_index"])
    truth = only(filter(
        case -> case.state_index == state_index,
        read_truth_cases(truth_table)))
    map = resolve_round5_sif_map([truth])
    prior_path = map.sif_on ? prior_on : prior_off
    prior = validate_round5_prior(
        load_retrieval_prior(truth.surface; path=prior_path), map;
        path=prior_path)
    prior.xa == baseline["active_state"] || error(
        "finite-difference state differs from saved benchmark state")

    evaluator = OCOForwardEvaluator(
        architecture=:GPU,
        float_type=Float32,
        nstreams=9,
        coefficient_path=coefficient_path,
        component_path=component_path)
    set_fixed_upper_co2_ppm!(evaluator, truth.fixed_upper_co2_ppm)
    round5_evaluator = Round5ForwardEvaluator(evaluator, map)

    plus_state = copy(prior.xa)
    minus_state = copy(prior.xa)
    plus_state[1] += step_hpa
    minus_state[1] -= step_hpa
    println("Computing +$step_hpa hPa instrument-space forward model...")
    plus = round5_evaluator(plus_state)
    println("Computing -$step_hpa hPa instrument-space forward model...")
    minus = round5_evaluator(minus_state)

    finite_difference = (plus.measurement .- minus.measurement) ./ (2step_hpa)
    analytical = Float64.(baseline["jacobian"][:, 1])
    delta = analytical .- finite_difference
    noise = NCDataset(noise_path, "r") do dataset
        Float64.(dataset["noise_std_corrected"][:])
    end
    length(noise) == length(delta) || error("noise vector length differs")

    relative_error = fd_relative_l2(delta, finite_difference)
    noise_weighted_relative_error = fd_relative_l2(
        delta ./ noise, finite_difference ./ noise)
    max_absolute_error = maximum(abs, delta)
    gate = 0.005
    passed = relative_error <= gate
    acceptance = passed ? "PASS" : "FAIL"

    starts = Int.(baseline["band_starts"])
    stops = Int.(baseline["band_stops"])
    band_relative_errors = [fd_relative_l2(
        @view(delta[first_index:last_index]),
        @view(finite_difference[first_index:last_index]))
        for (first_index, last_index) in zip(starts, stops)]

    mkpath(dirname(output_path))
    jldsave(output_path;
        implementation,
        state_index,
        step_hpa,
        analytical,
        finite_difference,
        delta,
        relative_error,
        noise_weighted_relative_error,
        max_absolute_error,
        band_relative_errors,
        gate,
        passed)
    @printf("implementation=%s state=%03d step_hpa=%.6g\n",
            implementation, state_index, step_hpa)
    @printf("pressure relative L2 error=%.9e\n", relative_error)
    @printf("pressure noise-weighted relative L2 error=%.9e\n",
            noise_weighted_relative_error)
    @printf("pressure max absolute error=%.9e\n", max_absolute_error)
    @printf("per-band relative L2=(%.9e, %.9e, %.9e)\n",
            band_relative_errors...)
    println("acceptance=$acceptance gate=$gate")
    println("Saved $output_path")
    require_pass = lowercase(get(ENV, "ROUND5_FD_REQUIRE_PASS", "0")) in
        ("1", "true", "yes", "on")
    require_pass && !passed && error(
        "surface-pressure finite-difference gate failed")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main_pressure_fd()
