#!/usr/bin/env julia

"""
Run a state/perturbation subset of one retrieval class.

Production examples (device 1 is the default):

    RETRIEVAL_CLASS=corrected julia --project=. \
        RRS_XCO2/inversion/run_retrievals.jl
    RETRIEVAL_CLASS=uncorrected julia --project=. \
        RRS_XCO2/inversion/run_retrievals.jl

Environment controls:

- `RETRIEVAL_CLASS=corrected|uncorrected|paired` (required). `paired` runs
  adjacent corrected/uncorrected experiments on
  the same host, recreating their identical deterministic normalized draw.
- `FIRST_STATE`, `LAST_STATE` (default 1,64; only no-SIF states are selected)
- `AEROSOL_CASE_FILTER=all|none|aerosol` (default `all`)
- `FIRST_PERTURBATION`, `LAST_PERTURBATION` (default 1,10)
- `RETRIEVAL_ARCH=GPU|CPU` (default GPU)
- `RETRIEVAL_FLOAT_TYPE=Float32|Float64` (default Float32)
- `CUDA_DEVICE` (default 1)
- `FORCE=1` to replace complete products
- `FAIL_FAST=0` to continue after a failed retrieval (default is fail fast)
"""

using Dates
using NCDatasets
using Printf
using Sockets: gethostname

include(joinpath(@__DIR__, "OptimalEstimation.jl"))
include(joinpath(@__DIR__, "RetrievalCases.jl"))
include(joinpath(@__DIR__, "RetrievalState.jl"))
include(joinpath(@__DIR__, "VSmartMOMForward.jl"))
include(joinpath(@__DIR__, "RetrievalOutput.jl"))
using .OptimalEstimation
using .RetrievalCases
using .RetrievalState
using .VSmartMOMForward
using .RetrievalOutput

env_flag(name, default="0") = lowercase(get(ENV, name, default)) in
    ("1", "true", "yes", "on")

function output_complete(path)
    isfile(path) || return false
    try
        return NCDataset(path) do dataset
            get(dataset.attrib, "retrieval_complete", 0) == 1
        end
    catch
        return false
    end
end

function selected_class()
    value = get(ENV, "RETRIEVAL_CLASS", "")
    value in ("corrected", "uncorrected", "paired") || error(
        "RETRIEVAL_CLASS must be corrected, uncorrected, or paired")
    return Symbol(value)
end

function selected_float_type()
    value = get(ENV, "RETRIEVAL_FLOAT_TYPE", "Float32")
    value == "Float32" && return Float32
    value == "Float64" && return Float64
    error("RETRIEVAL_FLOAT_TYPE must be Float32 or Float64")
end

function main()
    measurement_class = selected_class()
    first_state = parse(Int, get(ENV, "FIRST_STATE", "1"))
    last_state = parse(Int, get(ENV, "LAST_STATE", "64"))
    first_perturbation = parse(Int, get(ENV, "FIRST_PERTURBATION", "1"))
    last_perturbation = parse(Int, get(ENV, "LAST_PERTURBATION", "10"))
    1 <= first_perturbation <= last_perturbation <= 10 || error(
        "perturbation limits must lie in 1:10")
    architecture_name = uppercase(get(ENV, "RETRIEVAL_ARCH", "GPU"))
    architecture_name in ("CPU", "GPU") || error(
        "RETRIEVAL_ARCH must be CPU or GPU")
    architecture = Symbol(architecture_name)
    float_type = selected_float_type()
    force = env_flag("FORCE")
    fail_fast = env_flag("FAIL_FAST", "1")
    aerosol_case_filter = lowercase(get(ENV, "AEROSOL_CASE_FILTER", "all"))
    aerosol_case_filter in ("all", "none", "aerosol") || error(
        "AEROSOL_CASE_FILTER must be all, none, or aerosol")

    truth_cases = filter(read_no_sif_truth_cases()) do truth
        truth_has_aerosol = truth.aerosol_case != :none
        aerosol_case_filter == "all" ||
            (aerosol_case_filter == "aerosol") == truth_has_aerosol
    end
    all_experiments = build_experiments(truth_cases)
    experiments = filter(all_experiments) do experiment
        (measurement_class == :paired ||
         experiment.measurement_class == measurement_class) &&
        first_state <= experiment.truth.state_index <= last_state &&
        first_perturbation <= experiment.noise_index <= last_perturbation
    end
    isempty(experiments) && error("the requested subset contains no experiments")
    settings = OESettings()

    println("="^78)
    println("RRS-XCO2 retrieval suite")
    cuda_device = get(ENV, "CUDA_DEVICE", "1")
    println("host=$(gethostname()) class=$measurement_class architecture=$architecture " *
            "float_type=$float_type CUDA_DEVICE=$cuda_device")
    println("experiments=$(length(experiments)) state_range=$first_state:$last_state " *
            "perturbations=$first_perturbation:$last_perturbation " *
            "aerosol_case_filter=$aerosol_case_filter nstreams=9")
    println("started_utc=$(now(UTC))")
    println("="^78)

    evaluator = OCOForwardEvaluator(;
        architecture, float_type, nstreams=9)
    failures = 0
    completed = 0
    for (sequence, experiment) in enumerate(experiments)
        output_path = retrieval_output_path(experiment)
        if output_complete(output_path) && !force
            println("[$sequence/$(length(experiments))] skip complete $output_path")
            continue
        end
        truth = experiment.truth
        # Layers 1:4 have zero prior variance and are not part of the active
        # state. The synthetic truth profiles are vertically uniform, so fix
        # these layers to this scene's truth value rather than an unconditional
        # 400 ppm. This keeps every 380/400/420/440 ppm scene representable.
        set_fixed_upper_co2_ppm!(evaluator, truth.xco2_ppm)
        @printf("[%d/%d] state=%03d perturbation=%02d class=%s surface=%s aerosol=%s XCO2=%.0f\n",
                sequence, length(experiments), truth.state_index,
                experiment.noise_index, String(experiment.measurement_class),
                String(truth.surface), String(truth.aerosol_case), truth.xco2_ppm)
        @printf("  fixed_upper_co2_layers=1:4 fixed_upper_co2_ppm=%.1f\n",
                truth.xco2_ppm)
        try
            prior = load_retrieval_prior(truth.surface)
            realization = load_measurement_realization(experiment)
            callback = record -> @printf(
                "  trial=%d iteration=%d accepted=%d gamma=%.4g d_sigma_scaled=%.5g chi2=(%.4g,%.4g,%.4g) time=%.3fs\n",
                record.trial, record.iteration, record.accepted,
                record.gamma, record.d_sigma_sq_scaled,
                record.band_chi_squared..., record.evaluation_seconds)
            result = solve_optimal_estimation(
                evaluator,
                realization.perturbed,
                realization.variance,
                prior.xa,
                prior.Sa;
                settings,
                record_callback=callback)
            xco2_diagnostics = (
                a_priori_ppm=column_averaged_co2_ppm(evaluator, prior.xa),
                trial_ppm=[column_averaged_co2_ppm(evaluator, record.state)
                           for record in result.records],
                final_ppm=column_averaged_co2_ppm(evaluator, result.final_state),
            )
            provenance = Dict{String,Any}()
            VSmartMOMForward.RRSXCO2Common.write_absco_provenance!(provenance)
            VSmartMOMForward.RRSXCO2Common.write_fourier_convergence_provenance!(
                provenance)
            write_retrieval_result(
                experiment, realization, result, prior.xa, prior.Sa,
                prior.parameter_names; output_path, settings,
                solar_spectrum_path=VSmartMOMForward.RRSXCO2Common.SOLAR_OUT,
                xco2_diagnostics,
                provenance,
                overwrite=true)
            completed += 1
            @printf("  outcome=%d converged=%d fit_ok=%d XCO2=%.6f ppm final_chi2=(%.4g,%.4g,%.4g) output=%s\n",
                    result.outcome, result.converged, result.fit_quality_ok,
                    xco2_diagnostics.final_ppm,
                    result.final_band_chi_squared..., output_path)
        catch exception
            failures += 1
            showerror(stderr, exception, catch_backtrace())
            println(stderr)
            error_path = replace(output_path, r"\.nc$" => ".error.log")
            open(error_path, "w") do io
                println(io, "failed_utc=$(now(UTC))")
                println(io, "host=$(gethostname())")
                showerror(io, exception, catch_backtrace())
                println(io)
            end
            fail_fast && rethrow()
        end
    end
    println("finished_utc=$(now(UTC)) completed=$completed failures=$failures")
    failures == 0 || error("$failures retrievals failed")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
