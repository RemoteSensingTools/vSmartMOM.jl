#!/usr/bin/env julia

"""
Evaluate one deterministic, randomly selected Round-5 aerosol scene for a
legacy or optimized checkout and save the complete measurement-space forward
model and analytical Jacobian.

This benchmark deliberately evaluates the surface-specific a-priori state,
which is the first forward/Jacobian call made by every retrieval. The truth
case supplies the fixed upper-atmosphere CO2 value, surface class, and fixed
SIF mode. Selection is reproducible from `ROUND5_BENCH_SEED` and is restricted
to aerosol truth scenes so the expensive aerosol path is exercised.

Required environment variables:

- `ROUND5_INVERSION_ROOT`: workflow inversion directory for this checkout.
- `ROUND5_TRUTH_TABLE`: corrected bottom-layer truth-state table.
- `ROUND5_PRIOR_OFF` / `ROUND5_PRIOR_ON`: generated Round-5 priors.
- `ROUND5_STOKES_COEFFICIENTS`: representative OCO Mueller coefficients.
- `ROUND5_SCENE_COMPONENTS`: aerosol component catalog.
- `ROUND5_BENCH_OUTPUT`: destination JLD2 file.
- `ROUND5_BENCH_IMPLEMENTATION`: `legacy`, `candidate_reference`, or
  `optimized`.
- `ROUND5_BENCH_COMMIT`: exact checkout commit under test.

The benchmark performs one untimed warm-up evaluation and then
`ROUND5_BENCH_REPEATS` timed evaluations (default one). Julia compilation and
input loading are outside the reported timing boundary.
"""

using JLD2
using Random
using SHA
using Statistics

function required_path(name::AbstractString; directory::Bool=false)
    value = get(ENV, name, "")
    isempty(value) && error("$name is required")
    valid = directory ? isdir(value) : isfile(value)
    valid || error("$name does not identify a $(directory ? "directory" : "file"): $value")
    return realpath(value)
end

function file_sha256(path::AbstractString)
    return open(path, "r") do stream
        bytes2hex(sha256(stream))
    end
end

const INVERSION_ROOT = required_path("ROUND5_INVERSION_ROOT"; directory=true)

include(joinpath(INVERSION_ROOT, "OptimalEstimation.jl"))
include(joinpath(INVERSION_ROOT, "RetrievalCases.jl"))
include(joinpath(INVERSION_ROOT, "RetrievalState.jl"))
include(joinpath(INVERSION_ROOT, "VSmartMOMForward.jl"))
include(joinpath(INVERSION_ROOT, "Round4SIFTruthConvention.jl"))
include(joinpath(INVERSION_ROOT, "Round5FixedSIF.jl"))
include(joinpath(INVERSION_ROOT, "Round5RetrievalCampaign.jl"))

using .RetrievalCases
using .RetrievalState
using .VSmartMOMForward
using .Round5FixedSIF
using .Round5RetrievalCampaign

function assert_implementation(source_path::AbstractString,
                               implementation::AbstractString)
    source = read(source_path, String)
    optimized_copy = occursin(
        "copy_parameters(evaluator.base_parameters; share_luts=true)", source)
    optimized_basis = occursin("jacobian_basis=:local", source)
    optimized_adding = occursin("jacobian_adding=:source", source)
    legacy_copy = occursin("deepcopy(evaluator.base_parameters)", source)
    if implementation == "optimized"
        all((optimized_copy, optimized_basis, optimized_adding)) || error(
            "optimized benchmark did not select all three optimized adapter paths")
        legacy_copy && error("optimized adapter still contains the legacy deepcopy")
    elseif implementation in ("legacy", "candidate_reference")
        legacy_copy || error("legacy benchmark did not select the deepcopy adapter")
        any((optimized_copy, optimized_basis, optimized_adding)) && error(
            "legacy adapter unexpectedly contains optimized paths")
    else
        error("ROUND5_BENCH_IMPLEMENTATION must be legacy, " *
              "candidate_reference, or optimized")
    end
    return nothing
end

function timing_value(timing::NamedTuple, name::Symbol)
    return haskey(timing, name) ? Float64(getproperty(timing, name)) : NaN
end

function main()
    truth_table = required_path("ROUND5_TRUTH_TABLE")
    prior_off = required_path("ROUND5_PRIOR_OFF")
    prior_on = required_path("ROUND5_PRIOR_ON")
    coefficient_path = required_path("ROUND5_STOKES_COEFFICIENTS")
    component_path = required_path("ROUND5_SCENE_COMPONENTS")
    output_path = abspath(get(ENV, "ROUND5_BENCH_OUTPUT", ""))
    isempty(output_path) && error("ROUND5_BENCH_OUTPUT is required")
    implementation = lowercase(get(ENV, "ROUND5_BENCH_IMPLEMENTATION", ""))
    checkpoint = lowercase(get(ENV, "ROUND5_BENCH_COMMIT", ""))
    occursin(r"^[0-9a-f]{40}$", checkpoint) || error(
        "ROUND5_BENCH_COMMIT must be a full lowercase commit SHA")
    assert_implementation(
        joinpath(INVERSION_ROOT, "VSmartMOMForward.jl"), implementation)

    seed = parse(UInt64, get(ENV, "ROUND5_BENCH_SEED", "20260909"))
    repeats = parse(Int, get(ENV, "ROUND5_BENCH_REPEATS", "1"))
    repeats >= 1 || error("ROUND5_BENCH_REPEATS must be positive")

    cases = filter(read_truth_cases(truth_table)) do truth
        truth.campaign == :bottom_layer_XCO2 && truth.aerosol_case != :none
    end
    isempty(cases) && error("truth table contains no bottom-layer aerosol cases")
    rng = Xoshiro(seed)
    truth = rand(rng, cases)
    map = resolve_round5_sif_map([truth])
    prior_path = map.sif_on ? prior_on : prior_off
    prior = validate_round5_prior(
        load_retrieval_prior(truth.surface; path=prior_path), map;
        path=prior_path)

    evaluator = OCOForwardEvaluator(
        architecture=:GPU,
        float_type=Float32,
        nstreams=9,
        coefficient_path=coefficient_path,
        component_path=component_path)
    set_fixed_upper_co2_ppm!(evaluator, truth.fixed_upper_co2_ppm)
    round5_evaluator = Round5ForwardEvaluator(evaluator, map)

    println("Round-5 migration benchmark")
    println("implementation=$implementation checkpoint=$checkpoint seed=$seed")
    println("state=$(truth.state_index) surface=$(truth.surface) " *
            "aerosol=$(truth.aerosol_case) sif=$(truth.sif_case) " *
            "bottom_co2_ppm=$(truth.bottom_co2_ppm)")
    println("Warm-up evaluation (excluded from timed result)...")
    warmup = nothing
    warmup_seconds = @elapsed warmup = round5_evaluator(prior.xa)

    wall_seconds = Vector{Float64}(undef, repeats)
    model_build_seconds = Vector{Float64}(undef, repeats)
    rt_linearized_seconds = Vector{Float64}(undef, repeats)
    instrument_seconds = Vector{Float64}(undef, repeats)
    evaluation = nothing
    repeat_forward_spread = 0.0
    repeat_jacobian_spread = 0.0
    for repetition in 1:repeats
        GC.gc(true)
        wall_seconds[repetition] = @elapsed begin
            evaluation = round5_evaluator(prior.xa)
        end
        model_build_seconds[repetition] = timing_value(
            evaluation.timing, :model_build_seconds)
        rt_linearized_seconds[repetition] = timing_value(
            evaluation.timing, :rt_linearized_seconds)
        instrument_seconds[repetition] = timing_value(
            evaluation.timing, :instrument_seconds)
        repeat_forward_spread = max(
            repeat_forward_spread,
            maximum(abs.(evaluation.measurement .- warmup.measurement)))
        repeat_jacobian_spread = max(
            repeat_jacobian_spread,
            maximum(abs.(evaluation.jacobian .- warmup.jacobian)))
        println("timed_repetition=$repetition wall_seconds=$(wall_seconds[repetition])")
    end

    mkpath(dirname(output_path))
    jldsave(output_path;
        implementation,
        checkpoint,
        seed,
        state_index=truth.state_index,
        surface=String(truth.surface),
        aerosol_case=String(truth.aerosol_case),
        sif_case=String(truth.sif_case),
        bottom_co2_ppm=truth.bottom_co2_ppm,
        column_xco2_ppm=truth.xco2_ppm,
        fixed_upper_co2_ppm=truth.fixed_upper_co2_ppm,
        active_state=prior.xa,
        active_parameter_names=prior.parameter_names,
        active_to_full=prior.active_to_full,
        measurement=evaluation.measurement,
        jacobian=evaluation.jacobian,
        band_starts=first.(evaluation.band_ranges),
        band_stops=last.(evaluation.band_ranges),
        warmup_seconds,
        wall_seconds,
        model_build_seconds,
        rt_linearized_seconds,
        instrument_seconds,
        repeat_forward_spread,
        repeat_jacobian_spread,
        truth_table_sha256=file_sha256(truth_table),
        prior_sha256=file_sha256(prior_path),
        coefficient_sha256=file_sha256(coefficient_path),
        component_sha256=file_sha256(component_path),
        config_sha256=file_sha256(VSmartMOMForward.RRSXCO2Common.CONFIG))
    println("Saved $output_path")
    println("median_timed_wall_seconds=$(median(wall_seconds))")
    println("repeat_forward_spread=$repeat_forward_spread")
    println("repeat_jacobian_spread=$repeat_jacobian_spread")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
