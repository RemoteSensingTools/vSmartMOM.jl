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
- `FIRST_STATE`, `LAST_STATE` (default 1,64)
- `SIF_CASE_FILTER=off|on|all` (default `off`). `on` selects every nonzero-SIF
  truth case; `off` preserves the original no-SIF campaign behavior.
- `AEROSOL_CASE_FILTER=all|none|aerosol` (default `all`)
- `RETRIEVAL_TRUTH_TABLE`: alternate campaign truth-state table.
- `RETRIEVAL_MEASUREMENT_DIR`: alternate synthetic OCO-radiance directory.
- `RETRIEVAL_NOISE_DIR`: alternate frozen-noise directory.
- `RETRIEVAL_PRIOR_PATH`: alternate generated a-priori NetCDF.
- `RETRIEVAL_OUTPUT_ROOT`: alternate root containing `corrected/`,
  `uncorrected/`, and the retrieval manifest. Campaigns must use distinct
  roots; this prevents products with reused state indices from colliding.
- `RETRIEVAL_EXTERNAL_SIF_OWNERSHIP_MARKER`: optional path to a durable
  handoff marker. It defaults to
  `<RETRIEVAL_OUTPUT_ROOT>/.control/sif_owned_externally`. If the marker exists,
  any selection containing SIF-on work is refused before model allocation;
  no-SIF selections remain unaffected.
- `RETRIEVAL_WRITE_MANIFEST=0` suppresses manifest writing for a secondary
  worker sharing an output root (default `1`).
- `FIRST_PERTURBATION`, `LAST_PERTURBATION` (default 1,11). Indices 1:10
  inject paired noise; index 11 is the exact unperturbed measurement. When it
  is selected, index 11 is computed first for each state, followed by 1:10.
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

function output_complete(path, truth=nothing)
    isfile(path) || return false
    try
        return NCDataset(path) do dataset
            get(dataset.attrib, "retrieval_complete", 0) == 1 || return false
            isnothing(truth) || validated_sif_provenance(
                dataset.attrib, truth.sif_case; source=path)
            return true
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
    last_perturbation = parse(Int, get(ENV, "LAST_PERTURBATION", "11"))
    1 <= first_perturbation <= last_perturbation <= UNPERTURBED_INDEX || error(
        "perturbation limits must lie in 1:$UNPERTURBED_INDEX")
    architecture_name = uppercase(get(ENV, "RETRIEVAL_ARCH", "GPU"))
    architecture_name in ("CPU", "GPU") || error(
        "RETRIEVAL_ARCH must be CPU or GPU")
    architecture = Symbol(architecture_name)
    float_type = selected_float_type()
    force = env_flag("FORCE")
    fail_fast = env_flag("FAIL_FAST", "1")
    sif_case_filter = lowercase(get(ENV, "SIF_CASE_FILTER", "off"))
    sif_case_filter in ("off", "on", "all") || error(
        "SIF_CASE_FILTER must be off, on, or all")
    aerosol_case_filter = lowercase(get(ENV, "AEROSOL_CASE_FILTER", "all"))
    aerosol_case_filter in ("all", "none", "aerosol") || error(
        "AEROSOL_CASE_FILTER must be all, none, or aerosol")

    truth_table = abspath(get(
        ENV, "RETRIEVAL_TRUTH_TABLE", default_truth_table()))
    measurement_directory = abspath(get(
        ENV, "RETRIEVAL_MEASUREMENT_DIR", default_measurement_directory()))
    noise_directory = abspath(get(
        ENV, "RETRIEVAL_NOISE_DIR", default_noise_directory()))
    prior_path = abspath(get(
        ENV, "RETRIEVAL_PRIOR_PATH", RetrievalState.DEFAULT_PRIOR_PATH))
    output_root = abspath(get(ENV, "RETRIEVAL_OUTPUT_ROOT", @__DIR__))

    parsed_truth_cases = read_truth_cases(truth_table)
    if only(unique(case.co2_profile_mode for case in parsed_truth_cases)) ==
            :bottom_layer
        # Reused state indices make accidental fallback to the full-column
        # inputs especially dangerous: those files exist and would otherwise
        # pass a simple existence test. Require every campaign-local path.
        for variable in ("RETRIEVAL_MEASUREMENT_DIR", "RETRIEVAL_NOISE_DIR",
                         "RETRIEVAL_PRIOR_PATH", "RETRIEVAL_OUTPUT_ROOT")
            haskey(ENV, variable) || error(
                "$variable is required for a bottom-layer CO2 campaign")
        end
    end
    truth_cases = filter(
            select_sif_truth_cases(parsed_truth_cases, sif_case_filter)) do truth
        truth_has_aerosol = truth.aerosol_case != :none
        aerosol_case_filter == "all" ||
            (aerosol_case_filter == "aerosol") == truth_has_aerosol
    end
    all_experiments = build_experiments(
        truth_cases; measurement_directory, noise_directory)
    experiments = filter(all_experiments) do experiment
        (measurement_class == :paired ||
         experiment.measurement_class == measurement_class) &&
        first_state <= experiment.truth.state_index <= last_state &&
        first_perturbation <= experiment.noise_index <= last_perturbation
    end
    isempty(experiments) && error("the requested subset contains no experiments")
    enforce_sif_ownership(output_root, experiments)
    settings = OESettings()

    println("="^78)
    println("RRS-XCO2 retrieval suite")
    cuda_device = get(ENV, "CUDA_DEVICE", "1")
    println("host=$(gethostname()) class=$measurement_class architecture=$architecture " *
            "float_type=$float_type CUDA_DEVICE=$cuda_device")
    println("experiments=$(length(experiments)) state_range=$first_state:$last_state " *
            "perturbations=$first_perturbation:$last_perturbation " *
            "sif_case_filter=$sif_case_filter " *
            "aerosol_case_filter=$aerosol_case_filter nstreams=9")
    println("truth_table=$truth_table")
    println("measurement_directory=$measurement_directory")
    println("noise_directory=$noise_directory")
    println("prior_path=$prior_path")
    println("output_root=$output_root")
    println("started_utc=$(now(UTC))")
    println("="^78)

    if env_flag("RETRIEVAL_WRITE_MANIFEST", "1")
        write_experiment_manifest(
            all_experiments;
            output_path=joinpath(output_root, "retrieval_manifest.dat"),
            inversion_root=output_root)
    end

    evaluator = OCOForwardEvaluator(;
        architecture, float_type, nstreams=9)
    failures = 0
    completed = 0
    for (sequence, experiment) in enumerate(experiments)
        output_path = retrieval_output_path(
            experiment; inversion_root=output_root)
        if output_complete(output_path, experiment.truth) && !force
            println("[$sequence/$(length(experiments))] skip complete $output_path")
            continue
        end
        truth = experiment.truth
        # Layers 1:4 have zero prior variance and are not part of the active
        # state. Uniform-column cases use their scene VMR; bottom-layer cases
        # use the separately tabulated 400 ppm background. Never infer this
        # fixed value from column XCO2, which is not a layer VMR.
        set_fixed_upper_co2_ppm!(evaluator, truth.fixed_upper_co2_ppm)
        @printf("[%d/%d] state=%03d perturbation=%02d class=%s surface=%s aerosol=%s sif=%s XCO2=%.6f\n",
                sequence, length(experiments), truth.state_index,
                experiment.noise_index, String(experiment.measurement_class),
                String(truth.surface), String(truth.aerosol_case),
                String(truth.sif_case), truth.xco2_ppm)
        @printf("  fixed_upper_co2_layers=1:4 fixed_upper_co2_ppm=%.1f\n",
                truth.fixed_upper_co2_ppm)
        try
            prior = load_retrieval_prior(truth.surface; path=prior_path)
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
            provenance["retrieval_campaign"] = String(truth.campaign)
            provenance["source_truth_table"] = truth_table
            provenance["source_apriori"] = prior_path
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
