#!/usr/bin/env julia

#=
Round 7 changes ONLY the observation: y_uncorrected - H[Cabannes+RRS-Rayleigh].
Use --project=/home/sanghavi/code/github/uni_vSmartMOM_round5_jacobian (or its
identical deployment checkout). All numerical modules, including Round5FixedSIF
and Round5RetrievalCampaign, are included unchanged from that active project.

Required: ROUND7_OBSERVATION_DIR (80 OCO2round7_NNN.nc + provenance.json + SHA256SUMS),
RETRIEVAL_TRUTH_TABLE, RETRIEVAL_PRIOR_PATH (actual Round-5 fixed-SIF prior),
RETRIEVAL_STOKES_COEFFICIENT_PATH, RETRIEVAL_SCENE_COMPONENTS_PATH,
RRS_XCO2_SIF_TEMPLATE_PATH, SOLAR_OUT, RRS_XCO2_CONFIG, RETRIEVAL_OUTPUT_ROOT.
ROUND7_SOURCE_PRIOR can relocate the source_prior_path recorded in the prior.
Output root must end in retrievals_nosif or retrievals_sif in a round7 tree.
SIF_CASE_FILTER=off|on; AEROSOL_CASE_FILTER=all|none|aerosol; FIRST/LAST_STATE
and FIRST/LAST_PERTURBATION select work. Only RETRIEVAL_CLASS=corrected is
accepted: corrected/ is a compatibility name, never an ideal-correction claim.

Identity recipe (see round7_identity_values): base_tree_sha256 is SHA256 of
`git ls-tree -r --full-tree HEAD` stdout. Code files are ROUND7_CODE_FILES in
their declared order. Input files are, in order: truth table, Round-5 prior,
its source prior, Stokes coefficients, scene components, SIF template, solar
spectrum, RT config, observation SHA256SUMS. Paths are not hashed. Export
ROUND7_CODE_CHECKPOINT_SHA, ROUND7_CODESET_SHA256, ROUND7_INPUT_SET_SHA256,
ROUND7_CAMPAIGN_IDENTITY_SHA256, ROUND7_OBSERVATION_MANIFEST_SHA256 (alias:
OBSERVATION_MANIFEST_SHA256). ROUND7_PRINT_IDENTITY=1 prints these exports
without writing anything, for deployment tooling; it does not certify inputs.

ROUND7_PREFLIGHT_ONLY=1 verifies all 80 checksums, every selected observation
(all 11 stored members), priors, maps, and existing outputs without loading
CUDA/the forward evaluator. ROUND7_SMOKE=1 allows at most four solves and
requires unchanged strict convergence AND fit-quality gates, including on
resume. A failed fit can be scientifically legitimate; its output is retained.
FORCE is deliberately unsupported. No old outputs are silently overwritten.
=#

using Dates
using NCDatasets
using Printf
using SHA
using Sockets: gethostname

const BASE_PROJECT_DIRECTORY = dirname(Base.active_project())
const BASE_INVERSION_DIRECTORY = joinpath(BASE_PROJECT_DIRECTORY,
    "sandbox", "workflows", "RRS_XCO2", "inversion")
for name in ("OptimalEstimation.jl", "RetrievalCases.jl", "RetrievalState.jl",
             "RetrievalOutput.jl", "Round4SIFTruthConvention.jl",
             "Round5FixedSIF.jl", "Round5RetrievalCampaign.jl")
    include(joinpath(BASE_INVERSION_DIRECTORY, name))
end
include(joinpath(@__DIR__, "Round7RetrievalCampaign.jl"))

using .OptimalEstimation
using .RetrievalCases
using .RetrievalState
using .RetrievalOutput
using .Round5FixedSIF
import .Round5RetrievalCampaign as R5
using .Round7RetrievalCampaign

env_flag(name, default="0") = lowercase(get(ENV, name, default)) in
    ("1", "true", "yes", "on")

# Even importing the GPU forward module is avoided in CPU-only preflight.
if !env_flag("ROUND7_PREFLIGHT_ONLY") && !env_flag("ROUND7_PRINT_IDENTITY")
    include(joinpath(BASE_INVERSION_DIRECTORY, "VSmartMOMForward.jl"))
    using .VSmartMOMForward

    function VSmartMOMForward.column_averaged_co2_ppm(
            evaluator::Round5ForwardEvaluator, state::AbstractVector)
        return VSmartMOMForward.column_averaged_co2_ppm(
            evaluator.core_evaluator, expand_round5_state(evaluator.map, state))
    end

    function VSmartMOMForward.set_fixed_upper_co2_ppm!(
            evaluator::Round5ForwardEvaluator, ppm::Real)
        VSmartMOMForward.set_fixed_upper_co2_ppm!(evaluator.core_evaluator, ppm)
        return evaluator
    end
end

function required_input_file(name)
    path = get(ENV, name, "")
    isfile(path) || error("$name must name an existing file: $path")
    return realpath(path)
end

function round7_inputs()
    truth = required_input_file("RETRIEVAL_TRUTH_TABLE")
    prior = required_input_file("RETRIEVAL_PRIOR_PATH")
    source_prior = NCDataset(prior) do ds
        path = get(ENV, "ROUND7_SOURCE_PRIOR", get(ds.attrib, "source_prior_path", ""))
        isfile(path) || error("missing Round-5 source prior; set ROUND7_SOURCE_PRIOR")
        realpath(path)
    end
    stokes = required_input_file("RETRIEVAL_STOKES_COEFFICIENT_PATH")
    components = required_input_file("RETRIEVAL_SCENE_COMPONENTS_PATH")
    sif = required_input_file("RRS_XCO2_SIF_TEMPLATE_PATH")
    solar = required_input_file("SOLAR_OUT")
    config = required_input_file("RRS_XCO2_CONFIG")
    observations = get(ENV, "ROUND7_OBSERVATION_DIR", "")
    isdir(observations) || error("ROUND7_OBSERVATION_DIR must name an existing directory")
    observations = realpath(observations)
    manifest = joinpath(observations, "SHA256SUMS")
    return (; truth, prior, source_prior, stokes, components, sif, solar, config,
            observations, manifest,
            ordered=[truth, prior, source_prior, stokes, components, sif, solar, config, manifest])
end

function checked_round7_identity(inputs)
    round7_runtime_provenance() # Check the loaded package, not just the checkout.
    checkpoint = strip(read(`git -C $BASE_PROJECT_DIRECTORY rev-parse HEAD`, String))
    checkpoint == ROUND7_BASE_CHECKPOINT || error("active project is not the pinned accelerated Round-5 checkpoint")
    isempty(read(`git -C $BASE_PROJECT_DIRECTORY status --porcelain --untracked-files=all`, String)) ||
        error("pinned Round-5 checkout is dirty; refusing changed base modules")
    tree = bytes2hex(sha256(read(`git -C $BASE_PROJECT_DIRECTORY ls-tree -r --full-tree HEAD`)))
    return round7_identity_values(; checkpoint, base_tree_sha256=tree,
        code_paths=[joinpath(@__DIR__, name) for name in ROUND7_CODE_FILES],
        input_paths=inputs.ordered, observation_manifest=inputs.manifest)
end

function round7_runtime_provenance()
    loaded_module = Round4SIFTruthConvention.RRSXCO2Common.vSmartMOM
    loaded_source = realpath(pathof(loaded_module))
    expected_source = realpath(joinpath(BASE_PROJECT_DIRECTORY, "src", "vSmartMOM.jl"))
    loaded_source == expected_source || error(
        "loaded vSmartMOM source $loaded_source is not the pinned project source $expected_source")
    return Dict{String,Any}(
        "round7_base_project_path" => realpath(BASE_PROJECT_DIRECTORY),
        "round7_vsmartmom_source_path" => loaded_source,
        "round7_vsmartmom_source_sha256" => file_sha256(loaded_source),
        "round7_project_sha256" => file_sha256(joinpath(BASE_PROJECT_DIRECTORY, "Project.toml")),
        "round7_julia_manifest_sha256" => file_sha256(joinpath(BASE_PROJECT_DIRECTORY, "Manifest.toml")))
end

function prepare_round7()
    env_flag("FORCE") && error("FORCE is unsupported: Round-7 outputs are immutable")
    get(ENV, "RETRIEVAL_CLASS", "corrected") == "corrected" ||
        error("Round 7 runs only imperfectly corrected observations; no uncorrected rerun")
    get(ENV, "RETRIEVAL_FLOAT_TYPE", "Float32") == "Float32" ||
        error("Round 7 retains the Round-5 Float32 forward model")
    get(ENV, "RETRIEVAL_NSTREAMS", "9") == "9" || error("Round 7 requires nstreams=9")
    architecture = Symbol(uppercase(get(ENV, "RETRIEVAL_ARCH", "GPU")))
    architecture in (:CPU, :GPU) || error("RETRIEVAL_ARCH must be CPU or GPU")
    first_state = parse(Int, get(ENV, "FIRST_STATE", "1"))
    last_state = parse(Int, get(ENV, "LAST_STATE", "80"))
    first_perturbation = parse(Int, get(ENV, "FIRST_PERTURBATION", "1"))
    last_perturbation = parse(Int, get(ENV, "LAST_PERTURBATION", "11"))
    1 <= first_state <= last_state <= 80 || error("state range must lie in 1:80")
    1 <= first_perturbation <= last_perturbation <= 11 || error("perturbation range must lie in 1:11")
    sif_mode = lowercase(get(ENV, "SIF_CASE_FILTER", "off"))
    sif_mode in ("on", "off") || error("SIF_CASE_FILTER must be on or off")
    aerosol_mode = lowercase(get(ENV, "AEROSOL_CASE_FILTER", "all"))
    aerosol_mode in ("all", "none", "aerosol") || error("invalid AEROSOL_CASE_FILTER")
    inputs = round7_inputs()
    expected_identity = checked_round7_identity(inputs)
    if env_flag("ROUND7_PRINT_IDENTITY")
        for key in sort(collect(keys(expected_identity)))
            println("export $key=$(expected_identity[key])")
        end
        return nothing
    end
    identity = round7_campaign_identity(; expected=expected_identity)
    manifest = validate_round7_manifest(inputs.observations;
        expected_sha256=identity["round7_observation_manifest_sha256"])
    truths = read_truth_cases(inputs.truth)
    length(truths) == 80 && sort([t.state_index for t in truths]) == collect(1:80) &&
        all(t -> t.campaign == :bottom_layer_XCO2 &&
            t.sif_case in (:off, :angular_integral760_0p5), truths) ||
        error("Round 7 requires the validated 80-scene bottom-layer truth table")
    selected = filter(select_sif_truth_cases(truths, sif_mode)) do truth
        first_state <= truth.state_index <= last_state &&
            (aerosol_mode == "all" || (aerosol_mode == "aerosol") == (truth.aerosol_case != :none))
    end
    isempty(selected) && error("selected Round-7 subset is empty")
    map = R5.resolve_round5_sif_map(selected; requested=sif_mode)
    haskey(ENV, "RETRIEVAL_OUTPUT_ROOT") || error("RETRIEVAL_OUTPUT_ROOT is required")
    output_root = validate_round7_output_root(ENV["RETRIEVAL_OUTPUT_ROOT"]; sif_on=map.sif_on)
    # Build from all truths to keep IDs stable across partition/filter choices.
    selected_indices = Set(t.state_index for t in selected)
    experiments = filter(round7_experiments(truths, inputs.observations)) do e
        e.truth.state_index in selected_indices &&
            first_perturbation <= e.noise_index <= last_perturbation
    end
    enforce_sif_ownership(output_root, experiments)
    smoke = env_flag("ROUND7_SMOKE")
    smoke && length(experiments) > 4 && error("ROUND7_SMOKE requires an explicit subset of at most four solves")
    priors = Dict(surface => validate_round7_prior(
        load_retrieval_prior(surface; path=inputs.prior), map;
        path=inputs.prior, source_prior_path=inputs.source_prior)
        for surface in unique(t.surface for t in selected))
    truth_sha256 = file_sha256(inputs.truth)
    observations = Dict(t.state_index => load_round7_observation(manifest, t, map;
        truth_table_sha256=truth_sha256, prior_sha256=file_sha256(inputs.prior),
        analyzer_sha256=file_sha256(inputs.stokes)) for t in selected)
    # Preserve the original stored UInt64 seed, including values above 2^53.
    # The base constructor's canonical seeds are not authoritative for this release.
    experiments = [RetrievalExperiment(e.retrieval_index, e.pair_index, e.truth,
        e.noise_index, :corrected, observations[e.truth.state_index].seeds[e.noise_index],
        e.measurement_path, e.noise_path) for e in experiments]
    provenance = Dict(t.state_index => round7_output_provenance(map;
        prior_path=inputs.prior, observation=observations[t.state_index],
        campaign_identity=identity, truth_table_sha256=truth_sha256) for t in selected)
    for record in values(provenance)
        record["round7_architecture"] = String(architecture)
        merge!(record, round7_runtime_provenance())
    end
    # Inspect every selected existing output before any expensive solve/write.
    complete = Dict{String,Bool}()
    for e in experiments
        path = retrieval_output_path(e; inversion_root=output_root)
        complete[path] = validate_round7_existing_output(path, e,
            round7_realization(observations[e.truth.state_index], e.noise_index),
            priors[e.truth.surface], map;
            provenance=provenance[e.truth.state_index], smoke)
    end
    println("Round 7 preflight passed: $(length(selected)) scenes, $(length(experiments)) solves; " *
        "80 observation hashes; fixed Round-5 SIF/tight UTLS; Float32; nstreams=9")
    return (; inputs, identity, manifest, experiments, priors, observations, provenance,
            map, output_root, architecture, smoke, complete)
end

function run_round7(c)
    all(values(c.complete)) && return println("All selected Round-7 outputs already validated; no evaluator needed")
    core = OCOForwardEvaluator(; architecture=c.architecture, float_type=Float32,
        nstreams=9, coefficient_path=c.inputs.stokes, component_path=c.inputs.components)
    evaluator = Round5ForwardEvaluator(core, c.map)
    settings = OESettings()
    completed = 0
    failures = 0
    for (sequence, experiment) in enumerate(c.experiments)
        truth = experiment.truth
        output_path = retrieval_output_path(experiment; inversion_root=c.output_root)
        realization = round7_realization(c.observations[truth.state_index], experiment.noise_index)
        prior = c.priors[truth.surface]
        provenance = copy(c.provenance[truth.state_index])
        # Recheck on resume, including smoke gates; never reinterpret mismatch as missing.
        if validate_round7_existing_output(output_path, experiment, realization, prior, c.map;
                                          provenance, smoke=c.smoke)
            println("[$sequence/$(length(c.experiments))] skip validated $output_path")
            continue
        end
        println("[$sequence/$(length(c.experiments))] state=$(truth.state_index) " *
            "perturbation=$(experiment.noise_index) imperfectly_corrected")
        set_fixed_upper_co2_ppm!(evaluator, truth.fixed_upper_co2_ppm)
        try
            callback = r -> @printf("  trial=%d iteration=%d accepted=%d d_sigma=%.5g chi2=(%.4g,%.4g,%.4g)\n",
                r.trial, r.iteration, r.accepted, r.d_sigma_sq_scaled, r.band_chi_squared...)
            result = solve_optimal_estimation(evaluator, realization.perturbed,
                realization.variance, prior.xa, prior.Sa; settings, record_callback=callback)
            diagnostics = (
                a_priori_ppm=column_averaged_co2_ppm(evaluator, prior.xa),
                trial_ppm=[column_averaged_co2_ppm(evaluator, r.state) for r in result.records],
                final_ppm=column_averaged_co2_ppm(evaluator, result.final_state))
            VSmartMOMForward.RRSXCO2Common.write_absco_provenance!(provenance)
            VSmartMOMForward.RRSXCO2Common.write_fourier_convergence_provenance!(provenance)
            provenance["retrieval_campaign"] = String(truth.campaign)
            provenance["source_truth_table"] = c.inputs.truth
            provenance["source_apriori"] = c.inputs.prior
            # Private staging directory + non-overwriting link publishes a complete
            # result atomically. A competing worker can never clobber an old file.
            mkpath(dirname(output_path))
            mktempdir(dirname(output_path); prefix=".round7-") do staging
                temporary = joinpath(staging, "result.nc")
                write_retrieval_result(experiment, realization, result, prior.xa, prior.Sa,
                    prior.parameter_names; output_path=temporary, settings,
                    solar_spectrum_path=c.inputs.solar, jacobian_flavor=R5.ROUND5_JACOBIAN_FLAVOR,
                    active_to_core=active_core_indices(c.map), xco2_diagnostics=diagnostics,
                    provenance, overwrite=false)
                validate_round7_existing_output(temporary, experiment, realization, prior, c.map;
                                                 provenance=c.provenance[truth.state_index])
                hardlink(temporary, output_path)
            end
            completed += 1
            @printf("  converged=%d fit_ok=%d XCO2=%.6f chi2=(%.4g,%.4g,%.4g) output=%s\n",
                result.converged, result.fit_quality_ok, diagnostics.final_ppm,
                result.final_band_chi_squared..., output_path)
            c.smoke && round7_smoke_gate(result.converged, result.fit_quality_ok; source=output_path)
        catch exception
            failures += 1
            showerror(stderr, exception, catch_backtrace())
            println(stderr)
            # No error-log overwrite, no output deletion, no relaxed smoke threshold.
            (c.smoke || env_flag("FAIL_FAST", "1")) && rethrow()
        end
    end
    println("finished_utc=$(now(UTC)) host=$(gethostname()) completed=$completed failures=$failures")
    failures == 0 || error("$failures Round-7 retrievals failed")
end

function main()
    campaign = prepare_round7()
    isnothing(campaign) && return
    env_flag("ROUND7_PREFLIGHT_ONLY") && return campaign
    return run_round7(campaign)
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
