#!/usr/bin/env julia

"""
Run the bottom-layer-XCO2 round-4 retrieval model.

Round 4 knows SIF exactly at 759 nm while retaining the already validated
30-column `OCO_RRS_synth` forward/Jacobian implementation internally.  The
boundary wrapper expands the 28-coordinate SIF-off or 29-coordinate SIF-on
state and applies the exact Jacobian chain rule.  Truth spectra are read only;
this program never regenerates them.

This runner accepts the ordinary `run_retrievals.jl` selection controls plus:

- `ROUND4_SIF_MODE=auto|off|on` (default `auto`). `auto` is allowed only for
  a homogeneous selected truth subset. An explicit value must agree with it.
- `ROUND4_KNOWN_WAVELENGTH_NM=759` (optional assertion; every other value is
  rejected).
- `RETRIEVAL_PRIOR_PATH` is required and must identify the matching generated
  28- or 29-coordinate round-4 prior.
- `RETRIEVAL_OUTPUT_ROOT` is required and must contain both `round4` and
  `sif759` in its path, preventing collision with earlier retrieval rounds.
- `RETRIEVAL_STOKES_COEFFICIENT_PATH`, `RETRIEVAL_SCENE_COMPONENTS_PATH`,
  and `RRS_XCO2_SIF_TEMPLATE_PATH` are required explicit files. This avoids
  hidden data dependencies in a fresh source-only Gattaca checkout.

The SIF-off model fixes both `(SIF760,mSIF)` to exact zero. The SIF-on model
fixes the full corrected-v2 truth-template radiance at 759 nm and retrieves
only `mSIF`; `SIF760` at the core's 760-nm reference is derived from it.
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
include(joinpath(@__DIR__, "Round4KnownSIF.jl"))
include(joinpath(@__DIR__, "Round4SIFTruthConvention.jl"))
include(joinpath(@__DIR__, "Round4RetrievalCampaign.jl"))

using .OptimalEstimation
using .RetrievalCases
using .RetrievalState
using .VSmartMOMForward
using .RetrievalOutput
using .Round4KnownSIF
using .Round4RetrievalCampaign

env_flag(name, default="0") = lowercase(get(ENV, name, default)) in
    ("1", "true", "yes", "on")

function required_input_file(name::AbstractString)
    haskey(ENV, name) && !isempty(ENV[name]) || error(
        "$name is required by the isolated round-4 runner")
    isfile(ENV[name]) || error("$name does not name a file: $(ENV[name])")
    return realpath(ENV[name])
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

function output_complete(path, truth, map, prior_path, campaign_identity)
    isfile(path) || return false
    try
        return NCDataset(path) do dataset
            get(dataset.attrib, "retrieval_complete", 0) == 1 || return false
            validated_sif_provenance(
                dataset.attrib, truth.sif_case; source=path)
            validate_round4_output_provenance(
                dataset.attrib, map; prior_path, campaign_identity,
                source=path)
            haskey(dataset, "active_core_parameter_index") || return false
            Int.(dataset["active_core_parameter_index"][:]) ==
                active_core_indices(map) || return false
            return true
        end
    catch
        return false
    end
end

function VSmartMOMForward.column_averaged_co2_ppm(
        evaluator::Round4ForwardEvaluator, state::AbstractVector)
    return VSmartMOMForward.column_averaged_co2_ppm(
        evaluator.core_evaluator,
        expand_round4_state(evaluator.map, state))
end

function VSmartMOMForward.set_fixed_upper_co2_ppm!(
        evaluator::Round4ForwardEvaluator, ppm::Real)
    VSmartMOMForward.set_fixed_upper_co2_ppm!(
        evaluator.core_evaluator, ppm)
    return evaluator
end

function main()
    measurement_class = selected_class()
    first_state = parse(Int, get(ENV, "FIRST_STATE", "1"))
    last_state = parse(Int, get(ENV, "LAST_STATE", "80"))
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

    for variable in ("RETRIEVAL_TRUTH_TABLE", "RETRIEVAL_MEASUREMENT_DIR",
                     "RETRIEVAL_NOISE_DIR", "RETRIEVAL_PRIOR_PATH",
                     "RETRIEVAL_OUTPUT_ROOT",
                     "RETRIEVAL_STOKES_COEFFICIENT_PATH",
                     "RETRIEVAL_SCENE_COMPONENTS_PATH",
                     "RRS_XCO2_SIF_TEMPLATE_PATH")
        haskey(ENV, variable) || error(
            "$variable is required by the isolated round-4 runner")
    end
    truth_table = abspath(ENV["RETRIEVAL_TRUTH_TABLE"])
    measurement_directory = abspath(ENV["RETRIEVAL_MEASUREMENT_DIR"])
    noise_directory = abspath(ENV["RETRIEVAL_NOISE_DIR"])
    prior_path = abspath(ENV["RETRIEVAL_PRIOR_PATH"])
    stokes_coefficient_path = required_input_file(
        "RETRIEVAL_STOKES_COEFFICIENT_PATH")
    scene_components_path = required_input_file(
        "RETRIEVAL_SCENE_COMPONENTS_PATH")
    sif_template_path = required_input_file("RRS_XCO2_SIF_TEMPLATE_PATH")
    output_root = validate_round4_output_root(ENV["RETRIEVAL_OUTPUT_ROOT"])
    campaign_identity = round4_campaign_identity()

    parsed_truth_cases = read_truth_cases(truth_table)
    all(case -> case.campaign == :bottom_layer_XCO2,
        parsed_truth_cases) || error(
        "round-4 runner requires the bottom-layer-XCO2 truth campaign")
    selected_truth_cases = filter(
            select_sif_truth_cases(parsed_truth_cases, sif_case_filter)) do truth
        truth_has_aerosol = truth.aerosol_case != :none
        aerosol_case_filter == "all" ||
            (aerosol_case_filter == "aerosol") == truth_has_aerosol
    end
    isempty(selected_truth_cases) && error(
        "the requested SIF/aerosol subset contains no truth scenes")
    known_wavelength = parse(Float64, get(
        ENV, "ROUND4_KNOWN_WAVELENGTH_NM", "759"))
    map = resolve_round4_sif_map(
        selected_truth_cases;
        requested=get(ENV, "ROUND4_SIF_MODE", "auto"),
        known_wavelength_nm=known_wavelength)

    # This barrier validates the existing truth -> OCO measurement -> frozen
    # noise chain. It performs no truth computation.
    require_sif_release_barrier(
        truth_table, selected_truth_cases,
        measurement_directory, noise_directory)
    all_experiments = build_experiments(
        selected_truth_cases; measurement_directory, noise_directory)
    experiments = filter(all_experiments) do experiment
        (measurement_class == :paired ||
         experiment.measurement_class == measurement_class) &&
        first_state <= experiment.truth.state_index <= last_state &&
        first_perturbation <= experiment.noise_index <= last_perturbation
    end
    isempty(experiments) && error("the requested subset contains no experiments")
    enforce_sif_ownership(output_root, experiments)

    # Validate every surface prior before allocating the expensive RT model.
    prior_cache = Dict{Symbol,RetrievalPrior}()
    for surface in unique(experiment.truth.surface for experiment in experiments)
        prior = load_retrieval_prior(surface; path=prior_path)
        prior_cache[surface] = validate_round4_prior(prior, map; path=prior_path)
    end

    settings = OESettings()
    println("="^78)
    println("RRS-XCO2 round-4 known-SIF759 retrieval suite")
    cuda_device = get(ENV, "CUDA_DEVICE", "1")
    sif_mode_name = map.sif_on ? "on" : "off"
    println("host=$(gethostname()) class=$measurement_class architecture=$architecture " *
            "float_type=$float_type CUDA_DEVICE=$cuda_device")
    println("state_model=$ROUND4_STATE_MODEL sif_mode=$sif_mode_name " *
            "active_state_count=$(active_state_count(map)) core_state_count=$CORE_STATE_COUNT")
    @printf("known_SIF wavelength=%.1f nm Lnu=%.17g mW m-2 sr-1 (cm-1)-1\n",
            ROUND4_KNOWN_WAVELENGTH_NM, map.Lν759)
    println("experiments=$(length(experiments)) state_range=$first_state:$last_state " *
            "perturbations=$first_perturbation:$last_perturbation " *
            "aerosol_case_filter=$aerosol_case_filter nstreams=9")
    println("truth_table=$truth_table")
    println("measurement_directory=$measurement_directory")
    println("noise_directory=$noise_directory")
    println("prior_path=$prior_path")
    println("stokes_coefficient_path=$stokes_coefficient_path")
    println("scene_components_path=$scene_components_path")
    println("sif_template_path=$sif_template_path")
    println("output_root=$output_root")
    checkpoint_sha = campaign_identity["round4_code_checkpoint_sha"]
    identity_sha = campaign_identity["round4_campaign_identity_sha256"]
    println("code_checkpoint=$checkpoint_sha")
    println("campaign_identity_sha256=$identity_sha")
    println("started_utc=$(now(UTC))")
    println("="^78)

    if env_flag("RETRIEVAL_WRITE_MANIFEST", "1")
        write_experiment_manifest(
            all_experiments;
            output_path=joinpath(output_root, "retrieval_manifest.dat"),
            inversion_root=output_root)
    end

    core_evaluator = OCOForwardEvaluator(;
        architecture, float_type, nstreams=9,
        coefficient_path=stokes_coefficient_path,
        component_path=scene_components_path)
    evaluator = Round4ForwardEvaluator(core_evaluator, map)
    active_to_core = active_core_indices(map)
    fixed_provenance = round4_output_provenance(
        map; prior_path, campaign_identity)

    failures = 0
    completed = 0
    for (sequence, experiment) in enumerate(experiments)
        output_path = retrieval_output_path(
            experiment; inversion_root=output_root)
        if output_complete(
                output_path, experiment.truth, map, prior_path,
                campaign_identity) && !force
            println("[$sequence/$(length(experiments))] skip complete $output_path")
            continue
        end

        truth = experiment.truth
        set_fixed_upper_co2_ppm!(evaluator, truth.fixed_upper_co2_ppm)
        @printf("[%d/%d] state=%03d perturbation=%02d class=%s surface=%s aerosol=%s sif=%s XCO2=%.6f\n",
                sequence, length(experiments), truth.state_index,
                experiment.noise_index, String(experiment.measurement_class),
                String(truth.surface), String(truth.aerosol_case),
                String(truth.sif_case), truth.xco2_ppm)
        @printf("  fixed_upper_co2_layers=1:4 fixed_upper_co2_ppm=%.1f\n",
                truth.fixed_upper_co2_ppm)
        try
            prior = prior_cache[truth.surface]
            realization = load_measurement_realization(experiment)
            validate_round4_realization(
                map, truth, realization;
                source="truth state $(truth.state_index) realization")
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
                final_ppm=column_averaged_co2_ppm(
                    evaluator, result.final_state),
            )
            provenance = copy(fixed_provenance)
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
                jacobian_flavor=ROUND4_JACOBIAN_FLAVOR,
                active_to_core,
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
            mkpath(dirname(error_path))
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
