#!/usr/bin/env julia

"""
Fail-closed preflight for the local round-4, known-SIF-at-759-nm, SIF-off
bottom-layer retrieval campaign.

The corresponding shell launcher owns the scientific schedule.  This file
validates the immutable inputs, the reduced 28-coordinate prior, any products
that would be resumed, and the common campaign identity before a GPU model is
allocated.  It never runs a retrieval.
"""

using LinearAlgebra
using NCDatasets
using Printf
using SHA

include(joinpath(@__DIR__, "RetrievalCases.jl"))
include(joinpath(@__DIR__, "RetrievalState.jl"))
using .RetrievalCases
using .RetrievalState

module RadianceValidation
include(joinpath(@__DIR__, "instrument", "validate_oco_radiances.jl"))
end

module NoiseValidation
include(joinpath(@__DIR__, "instrument", "validate_noise_covariances.jl"))
end

const ROUND4_CAMPAIGN_ID = "bottom_layer_round4_known_sif759_nosif_v1"
const ROUND4_IDENTITY_SCHEMA = 1
const ROUND4_STATE_MODEL = "round4_known_sif759"
const ROUND4_CO2_MODEL = "acos_mapped_tapered_vertical_correlation"
const ROUND4_JACOBIAN_FLAVOR =
    "OCO_RRS_synth_round4_known_sif759_boundary_chain"
const ROUND4_KNOWN_WAVELENGTH_NM = 759.0
const ROUND4_ACTIVE_TO_FULL = vcat(1, collect(6:32))
const ROUND4_ACTIVE_TO_CORE = collect(1:28)
const ROUND4_GLOBAL_STATE_SPEC =
    "1-5,11-15,21-25,31-35,41-45,51-55,61-65,71-75"
const ROUND4_SCHEDULE = join((
    "curry1:none:1-5,21-25,41-45,61-65",
    "wurst0:aerosol:11-15,51-55",
    "wurst1:aerosol:31-35,71-75",
), ";")

function parse_state_spec(value::AbstractString)
    states = Int[]
    for raw_token in split(value, ',')
        token = strip(raw_token)
        isempty(token) && error("state specification contains an empty token")
        bounds = split(token, '-'; limit=2)
        first_state = parse(Int, first(bounds))
        last_state = length(bounds) == 1 ? first_state : parse(Int, last(bounds))
        first_state <= last_state || error(
            "descending state range is not allowed: $token")
        append!(states, first_state:last_state)
    end
    isempty(states) && error("state specification selected no states")
    length(states) == length(unique(states)) || error(
        "state specification contains duplicate indices")
    return states
end

function required_path(name::AbstractString)
    haskey(ENV, name) || error("$name is required")
    return abspath(ENV[name])
end

file_sha256(path::AbstractString) = open(path, "r") do stream
    bytes2hex(sha256(stream))
end

function required_sha256(name::AbstractString)
    haskey(ENV, name) || error("$name is required")
    value = ENV[name]
    occursin(r"^[0-9a-f]{64}$", value) || error(
        "$name must be a lowercase 64-character SHA-256")
    return value
end

function required_git_sha(name::AbstractString)
    haskey(ENV, name) || error("$name is required")
    value = ENV[name]
    occursin(r"^[0-9a-f]{40}$", value) || error(
        "$name must be a lowercase 40-character Git SHA")
    return value
end

function path_is_within(path::AbstractString, root::AbstractString)
    relative = relpath(abspath(path), abspath(root))
    return relative == "." ||
        (relative != ".." && !startswith(relative, ".." * Base.Filesystem.path_separator))
end

function validate_path_isolation(campaign_root, prior_path, source_prior_path,
                                 output_root)
    active_output = joinpath(campaign_root, "retrievals")
    round3_output = joinpath(
        campaign_root,
        "retrievals_acos_mapped_tapered_vertical_correlation_nosif")
    round4_root = joinpath(campaign_root, "round4_known_sif759")
    prior_path != source_prior_path || error(
        "round 4 may not use its round-3 source prior directly")
    path_is_within(prior_path, round4_root) || error(
        "round-4 prior must be stored inside $round4_root")
    path_is_within(output_root, round4_root) || error(
        "round-4 output must be stored inside $round4_root")
    lowercase(basename(output_root)) == "retrievals_nosif" || error(
        "round-4 no-SIF output basename must be retrievals_nosif")
    for old_root in (active_output, round3_output)
        (path_is_within(output_root, old_root) ||
         path_is_within(old_root, output_root)) && error(
            "round-4 output overlaps an earlier retrieval namespace: $old_root")
    end
    return nothing
end

function _all_values(dataset, variable)
    value = dataset[variable]
    indices = ntuple(_ -> Colon(), ndims(value))
    return Array(value[indices...])
end

function _numeric_attribute(dataset, name)
    haskey(dataset.attrib, name) || error("prior is missing attribute $name")
    return Float64(dataset.attrib[name])
end

"""Validate the no-SIF prior and exact invariance of every non-SIF block."""
function validate_round4_nosif_prior(prior_path;
                                     expected_sha256,
                                     source_prior_path,
                                     source_prior_sha256)
    isfile(prior_path) || error("missing round-4 no-SIF prior: $prior_path")
    isfile(source_prior_path) || error(
        "missing approved round-3 source prior: $source_prior_path")
    file_sha256(prior_path) == expected_sha256 || error(
        "round-4 prior SHA mismatch")
    file_sha256(source_prior_path) == source_prior_sha256 || error(
        "round-3 source-prior SHA mismatch")

    NCDataset(prior_path) do selected
        get(selected.attrib, "apriori_complete", 0) == 1 || error(
            "round-4 prior is not marked apriori_complete=1")
        String(get(selected.attrib, "retrieval_state_model", "")) ==
            ROUND4_STATE_MODEL || error("wrong round-4 retrieval_state_model")
        String(get(selected.attrib, "round4_sif_case", "")) == "off" ||
            error("round-4 local launcher requires the SIF-off prior")
        _numeric_attribute(selected, "known_sif_wavelength_nm") ==
            ROUND4_KNOWN_WAVELENGTH_NM || error(
                "round-4 prior does not fix SIF at 759 nm")
        iszero(_numeric_attribute(
            selected, "known_sif_Lnu_mW_m-2_sr-1_per_cm-1")) || error(
                "SIF-off prior has nonzero known Lnu(759)")
        iszero(_numeric_attribute(
            selected, "known_sif_Llambda_mW_m-2_sr-1_nm-1")) || error(
                "SIF-off prior has nonzero known Llambda(759)")
        Int(get(selected.attrib, "active_state_count", -1)) == 28 || error(
            "SIF-off round-4 prior must advertise 28 active coordinates")
        String(get(selected.attrib, "co2_covariance_model", "")) ==
            ROUND4_CO2_MODEL || error("wrong round-4 CO2 covariance model")
        String(get(selected.attrib, "source_prior_sha256", "")) ==
            source_prior_sha256 || error(
                "round-4 prior does not name the approved source-prior SHA")
        abspath(String(get(selected.attrib, "source_prior_path", ""))) ==
            abspath(source_prior_path) || error(
                "round-4 prior records a different source-prior path")

        active = Int.(_all_values(selected, "active_parameter_index"))
        active == ROUND4_ACTIVE_TO_FULL || error(
            "round-4 no-SIF active/full mapping is not [1;6:32]")
        mask = Int.(_all_values(selected, "active_mask"))
        expected_mask = zeros(Int, length(mask))
        expected_mask[ROUND4_ACTIVE_TO_FULL] .= 1
        mask == expected_mask || error("round-4 active mask is inconsistent")

        xa = Float64.(_all_values(selected, "xa"))
        Sa = Float64.(_all_values(selected, "Sa"))
        Sa_active = Float64.(_all_values(selected, "Sa_active"))
        size(xa, 1) == 34 || error("round-4 full state must have 34 coordinates")
        size(Sa)[1:2] == (34, 34) || error(
            "round-4 full covariance must be 34x34 per surface")
        size(Sa_active)[1:2] == (28, 28) || error(
            "round-4 active covariance must be 28x28 per surface")
        all(iszero, @view xa[33:34, :]) || error(
            "SIF-off prior means must fix SIF760 and mSIF to zero")
        all(iszero, @view Sa[33:34, :, :]) || error(
            "SIF-off prior covariance rows 33:34 must be zero")
        all(iszero, @view Sa[:, 33:34, :]) || error(
            "SIF-off prior covariance columns 33:34 must be zero")

        NCDataset(source_prior_path) do source
            source_xa = Float64.(_all_values(source, "xa"))
            source_Sa = Float64.(_all_values(source, "Sa"))
            xa[1:32, :] == source_xa[1:32, :] || error(
                "round 4 changed a non-SIF prior mean")
            Sa[1:32, 1:32, :] == source_Sa[1:32, 1:32, :] || error(
                "round 4 changed a non-SIF covariance")
            for attribute in ("surface_order", "parameter_names",
                              "parameter_units", "state_order",
                              "aerosol_coordinate_transform")
                get(selected.attrib, attribute, nothing) ==
                    get(source.attrib, attribute, nothing) || error(
                        "round 4 changed non-SIF metadata $attribute")
            end
        end

        for surface in axes(Sa_active, 3)
            expected = Sa[ROUND4_ACTIVE_TO_FULL,
                          ROUND4_ACTIVE_TO_FULL, surface]
            Sa_active[:, :, surface] == expected || error(
                "Sa_active is inconsistent for surface $surface")
            isposdef(Symmetric(Sa_active[:, :, surface])) || error(
                "round-4 active covariance is not positive definite for surface $surface")
        end
    end
    return (; prior_path=abspath(prior_path), prior_sha256=expected_sha256,
            source_prior_path=abspath(source_prior_path),
            source_prior_sha256)
end

function truth_scene_path(truth, truth_root::AbstractString)
    directory = truth.aerosol_case == :none ? truth_root :
        joinpath(truth_root, "aerosol_chunked")
    return joinpath(directory,
                    @sprintf("hiressim_%03d.nc", truth.state_index))
end

function validate_selected_truth(states, cases, expected_scene_class)
    expected_scene_class in ("none", "aerosol") || error(
        "ROUND4_EXPECTED_SCENE_CLASS must be none or aerosol")
    by_index = Dict(case.state_index => case for case in cases)
    return map(states) do state
        haskey(by_index, state) || error("truth table does not contain state $state")
        truth = by_index[state]
        truth.sif_case == :off || error("state $state is not SIF-off")
        has_aerosol = truth.aerosol_case != :none
        has_aerosol == (expected_scene_class == "aerosol") || error(
            "state $state does not belong to the scheduled $expected_scene_class class")
        truth
    end
end

function round4_input_set_sha256(cases, truth_table, oco_root,
                                 noise_root, coefficients_path)
    records = [
        "truth_table $(file_sha256(truth_table))",
        "snr_coefficients $(file_sha256(coefficients_path))",
    ]
    selected = sort!(filter(case -> case.sif_case == :off, collect(cases));
                     by=case -> case.state_index)
    length(selected) == 40 || error(
        "bottom-layer campaign must contain exactly 40 SIF-off states")
    for truth in selected
        index = truth.state_index
        # The retrieval reads only the state table and these instrument-space
        # measurement/noise products. validate_file below independently
        # follows and validates each product's campaign-aware source_truth_scene
        # provenance (clear root versus aerosol_chunked); duplicating every
        # high-resolution truth file in this digest would be costly and would
        # pin data the retrieval itself never opens.
        for (kind, path) in (
                ("measurement", joinpath(
                    oco_root, @sprintf("OCO2sims_%03d.nc", index))),
                ("noise", joinpath(
                    noise_root, @sprintf("OCO2noise_%03d.nc", index))))
            isfile(path) || error("missing round-4 $kind input for state $index: $path")
            push!(records, "$kind $index $(file_sha256(path))")
        end
    end
    return bytes2hex(sha256(codeunits(join(records, '\n') * "\n")))
end

function round4_codeset_paths(repo_root)
    paths = String[]
    for relative in ("Project.toml", "Manifest.toml",
                     "config/oco_grass_3aerosol.yaml")
        path = joinpath(repo_root, relative)
        isfile(path) && push!(paths, path)
    end
    for relative_root in ("src", "ext", "RRS_XCO2/scripts")
        root = joinpath(repo_root, relative_root)
        isdir(root) || continue
        for (directory, _, files) in walkdir(root), file in files
            endswith(file, ".jl") && push!(paths, joinpath(directory, file))
        end
    end
    for relative in (
            "RRS_XCO2/inversion/OptimalEstimation.jl",
            "RRS_XCO2/inversion/RetrievalCases.jl",
            "RRS_XCO2/inversion/RetrievalState.jl",
            "RRS_XCO2/inversion/RetrievalOutput.jl",
            "RRS_XCO2/inversion/VSmartMOMForward.jl",
            "RRS_XCO2/inversion/Round4KnownSIF.jl",
            "RRS_XCO2/inversion/Round4SIFTruthConvention.jl",
            "RRS_XCO2/inversion/Round4RetrievalCampaign.jl",
            "RRS_XCO2/inversion/instrument/SyntheticOCO2.jl",
            "RRS_XCO2/inversion/run_round4_known_sif_retrievals.jl",
            "RRS_XCO2/inversion/preflight_round4_known_sif759_no_sif_retrievals.jl",
            "RRS_XCO2/inversion/run_round4_known_sif759_no_sif_partition.sh")
        path = joinpath(repo_root, relative)
        isfile(path) || error("round-4 code-set file is missing: $path")
        push!(paths, path)
    end
    return sort!(unique(paths))
end

function round4_codeset_sha256(repo_root)
    records = String[]
    for path in round4_codeset_paths(repo_root)
        push!(records, "$(relpath(path, repo_root)) $(file_sha256(path))")
    end
    return bytes2hex(sha256(codeunits(join(records, '\n') * "\n")))
end

function validate_round3_complete(campaign_root)
    output_root = joinpath(
        campaign_root,
        "retrievals_acos_mapped_tapered_vertical_correlation_nosif")
    for claim_name in (".state_claims", ".partition_claims")
        claim_root = joinpath(output_root, claim_name)
        if isdir(claim_root) && any(entry -> isdir(joinpath(claim_root, entry)),
                                   readdir(claim_root))
            error("round 3 still has active claims in $claim_root")
        end
    end
    states = parse_state_spec(ROUND4_GLOBAL_STATE_SPEC)
    for measurement_class in ("corrected", "uncorrected"), state in states,
            perturbation in vcat(11, collect(1:10))
        path = joinpath(output_root, measurement_class,
            @sprintf("retrieval_state%03d_perturbation%02d.nc",
                     state, perturbation))
        isfile(path) || error("round 3 is incomplete; missing $path")
        NCDataset(path) do dataset
            Int(get(dataset.attrib, "retrieval_complete", 0)) == 1 || error(
                "round-3 product is incomplete: $path")
        end
    end
    return true
end

function validate_output_campaign_identity(attributes;
                                           code_checkpoint,
                                           codeset_sha256,
                                           input_set_sha256,
                                           campaign_identity_sha256)
    expected = (
        "round4_campaign_identity_status" => "complete",
        "round4_code_checkpoint_sha" => code_checkpoint,
        "round4_codeset_sha256" => codeset_sha256,
        "round4_input_set_sha256" => input_set_sha256,
        "round4_campaign_identity_sha256" => campaign_identity_sha256,
    )
    for (name, value) in expected
        String(get(attributes, name, "")) == value || error(
            "existing product has a different $name")
    end
    return true
end

function validate_existing_outputs(output_root, selected_truth, prior_path;
                                   code_checkpoint,
                                   codeset_sha256,
                                   input_set_sha256,
                                   campaign_identity_sha256)
    checked = 0
    isdir(output_root) || return checked
    prior_cache = Dict{Symbol,RetrievalPrior}()
    for truth in selected_truth, measurement_class in ("corrected", "uncorrected"),
            perturbation in vcat(11, collect(1:10))
        path = joinpath(output_root, measurement_class,
            @sprintf("retrieval_state%03d_perturbation%02d.nc",
                     truth.state_index, perturbation))
        isfile(path) || continue
        prior = get!(prior_cache, truth.surface) do
            load_retrieval_prior(truth.surface; path=prior_path)
        end
        NCDataset(path) do dataset
            Int(get(dataset.attrib, "retrieval_complete", 0)) == 1 || error(
                "partial round-4 output requires inspection: $path")
            Int(get(dataset.attrib, "truth_state_index", -1)) ==
                truth.state_index || error("truth-state mismatch in $path")
            Int(get(dataset.attrib, "perturbation_index", -1)) ==
                perturbation || error("perturbation mismatch in $path")
            String(get(dataset.attrib, "measurement_class", "")) ==
                measurement_class || error("measurement-class mismatch in $path")
            String(get(dataset.attrib, "retrieval_state_model", "")) ==
                ROUND4_STATE_MODEL || error("wrong state model in $path")
            String(get(dataset.attrib, "round4_sif_case", "")) == "off" ||
                error("non-off SIF product found in local no-SIF output: $path")
            Float64(get(dataset.attrib,
                        "round4_known_sif_wavelength_nm", NaN)) == 759.0 ||
                error("wrong known-SIF wavelength in $path")
            String(get(dataset.attrib, "jacobian_flavor", "")) ==
                ROUND4_JACOBIAN_FLAVOR || error("wrong Jacobian flavor in $path")
            Int(get(dataset.attrib, "state_dimension", -1)) == 28 || error(
                "wrong active state dimension in $path")
            abspath(String(get(dataset.attrib, "source_apriori", ""))) ==
                abspath(prior_path) || error("wrong prior provenance in $path")
            String(get(dataset.attrib, "round4_prior_sha256", "")) ==
                file_sha256(prior_path) || error("wrong prior hash in $path")
            validate_output_campaign_identity(dataset.attrib;
                code_checkpoint, codeset_sha256, input_set_sha256,
                campaign_identity_sha256)
            Int.(_all_values(dataset, "active_core_parameter_index")) ==
                ROUND4_ACTIVE_TO_CORE || error(
                    "wrong active/core mapping in $path")
            Float64.(_all_values(dataset, "a_priori_state")) == prior.xa ||
                error("embedded prior mean mismatch in $path")
            Float64.(_all_values(dataset, "a_priori_covariance")) == prior.Sa ||
                error("embedded prior covariance mismatch in $path")
        end
        checked += 1
    end
    return checked
end

function campaign_identity_text(; campaign_root, truth_table, oco_root,
                                noise_root, output_root, prior_identity,
                                input_set_sha256, code_checkpoint,
                                codeset_sha256)
    fields = (
        "identity_schema" => string(ROUND4_IDENTITY_SCHEMA),
        "campaign_id" => ROUND4_CAMPAIGN_ID,
        "global_schedule" => ROUND4_SCHEDULE,
        "sif_case" => "off",
        "known_sif_wavelength_nm" => "759",
        "known_sif_Lnu" => "0",
        "known_sif_Llambda" => "0",
        "active_state_count" => "28",
        "active_to_full" => join(ROUND4_ACTIVE_TO_FULL, ','),
        "active_to_core" => join(ROUND4_ACTIVE_TO_CORE, ','),
        "co2_covariance_model" => ROUND4_CO2_MODEL,
        "prior_sha256" => prior_identity.prior_sha256,
        "prior_path" => prior_identity.prior_path,
        "source_prior_sha256" => prior_identity.source_prior_sha256,
        "source_prior_path" => prior_identity.source_prior_path,
        "input_set_sha256" => input_set_sha256,
        "code_checkpoint" => code_checkpoint,
        "codeset_sha256" => codeset_sha256,
        "truth_table" => abspath(truth_table),
        "measurement_directory" => abspath(oco_root),
        "noise_directory" => abspath(noise_root),
        "campaign_root" => abspath(campaign_root),
        "output_root" => abspath(output_root),
    )
    return join(("$key=$value" for (key, value) in fields), '\n') * "\n"
end

"""Atomically create or byte-for-byte verify the common campaign identity."""
function initialize_campaign_identity!(output_root, expected;
                                       wait_seconds::Real=30)
    control = joinpath(output_root, ".control")
    identity = joinpath(control, "campaign_identity.dat")
    lock = joinpath(control, ".identity_lock")
    mkpath(control)
    if isfile(identity)
        read(identity, String) == expected || error(
            "existing round-4 campaign identity differs")
        return identity
    end
    isempty(filter(path -> endswith(path, ".nc"),
                   collect(Iterators.flatten((
                       isdir(joinpath(output_root, class)) ?
                           readdir(joinpath(output_root, class); join=true) : String[]
                       for class in ("corrected", "uncorrected")))))) || error(
        "round-4 products exist without a campaign identity")
    acquired = try
        mkdir(lock)
        true
    catch exception
        exception isa Base.IOError || exception isa SystemError || rethrow()
        false
    end
    if !acquired
        deadline = time() + wait_seconds
        while time() < deadline && !isfile(identity)
            sleep(0.1)
        end
        isfile(identity) || error(
            "identity lock exists without a published identity: $lock")
        read(identity, String) == expected || error(
            "concurrent round-4 campaign identity differs")
        return identity
    end
    temporary = identity * ".tmp.$(getpid())"
    try
        open(temporary, "w") do io
            write(io, expected)
        end
        mv(temporary, identity)
    finally
        isfile(temporary) && rm(temporary)
        isdir(lock) && rm(lock)
    end
    return identity
end

function git_checkpoint(repo_root)
    return strip(read(`git -C $repo_root rev-parse HEAD`, String))
end

function main()
    states = parse_state_spec(get(ENV, "ROUND4_EXPECTED_STATES", ""))
    expected_scene_class = get(ENV, "ROUND4_EXPECTED_SCENE_CLASS", "")
    campaign_root = required_path("BOTTOM_RETRIEVAL_CAMPAIGN_ROOT")
    prior_path = required_path("RETRIEVAL_PRIOR_PATH")
    source_prior_path = required_path("ROUND4_SOURCE_PRIOR_PATH")
    output_root = required_path("RETRIEVAL_OUTPUT_ROOT")
    repo_root = required_path("ROUND4_REPO_ROOT")
    validate_path_isolation(
        campaign_root, prior_path, source_prior_path, output_root)

    truth_root = joinpath(campaign_root, "truth")
    truth_table = joinpath(truth_root, "true_states.dat")
    oco_root = joinpath(truth_root, "OCO_radiances")
    noise_root = joinpath(oco_root, "noise_covariances")
    coefficients_path = joinpath(
        @__DIR__, "instrument", "representative_snr_coefficients.nc")
    cases = read_truth_cases(truth_table)
    selected_truth = validate_selected_truth(
        states, cases, expected_scene_class)

    actual_prior_sha = file_sha256(prior_path)
    actual_source_sha = file_sha256(source_prior_path)
    actual_input_sha = round4_input_set_sha256(
        cases, truth_table, oco_root, noise_root, coefficients_path)
    actual_checkpoint = git_checkpoint(repo_root)
    actual_codeset_sha = round4_codeset_sha256(repo_root)

    # Candidate hashes are useful only after the scientific shape and exact
    # non-SIF inheritance of the proposed prior have passed validation.
    candidate_prior_identity = validate_round4_nosif_prior(
        prior_path; expected_sha256=actual_prior_sha,
        source_prior_path, source_prior_sha256=actual_source_sha)

    if get(ENV, "ROUND4_NOSIF_PRINT_CANDIDATE_HASHES", "0") == "1"
        println("ROUND4_PRIOR_SHA256=$actual_prior_sha")
        println("ROUND4_SOURCE_PRIOR_SHA256=$actual_source_sha")
        println("ROUND4_INPUT_SET_SHA256=$actual_input_sha")
        println("ROUND4_CODE_CHECKPOINT_SHA=$actual_checkpoint")
        println("ROUND4_CODESET_SHA256=$actual_codeset_sha")
        return
    end

    expected_prior_sha = required_sha256("ROUND4_PRIOR_SHA256")
    expected_source_sha = required_sha256("ROUND4_SOURCE_PRIOR_SHA256")
    expected_input_sha = required_sha256("ROUND4_INPUT_SET_SHA256")
    expected_codeset_sha = required_sha256("ROUND4_CODESET_SHA256")
    expected_checkpoint = required_git_sha("ROUND4_CODE_CHECKPOINT_SHA")
    actual_checkpoint == expected_checkpoint || error(
        "current Git checkpoint differs from ROUND4_CODE_CHECKPOINT_SHA")
    actual_codeset_sha == expected_codeset_sha || error(
        "round-4 code-set digest differs from ROUND4_CODESET_SHA256")
    actual_input_sha == expected_input_sha || error(
        "round-4 truth/measurement/noise input digest differs")

    candidate_prior_identity.prior_sha256 == expected_prior_sha || error(
        "round-4 prior SHA differs from ROUND4_PRIOR_SHA256")
    candidate_prior_identity.source_prior_sha256 == expected_source_sha ||
        error("source-prior SHA differs from ROUND4_SOURCE_PRIOR_SHA256")
    prior_identity = candidate_prior_identity

    coefficients = NoiseValidation.read_representative_snr_coefficients(
        coefficients_path)
    for truth in selected_truth
        state = truth.state_index
        radiance_path = joinpath(oco_root, @sprintf("OCO2sims_%03d.nc", state))
        noise_path = joinpath(
            noise_root, @sprintf("OCO2noise_%03d.nc", state))
        RadianceValidation.validate_file(radiance_path, state)
        NoiseValidation.validate_file(
            noise_path, radiance_path, coefficients, state)
    end

    get(ENV, "ROUND4_REQUIRE_ROUND3_COMPLETE", "1") in ("0", "1") ||
        error("ROUND4_REQUIRE_ROUND3_COMPLETE must be 0 or 1")
    get(ENV, "ROUND4_REQUIRE_ROUND3_COMPLETE", "1") == "1" &&
        validate_round3_complete(campaign_root)
    identity_text = campaign_identity_text(
        ; campaign_root, truth_table, oco_root, noise_root, output_root,
        prior_identity, input_set_sha256=actual_input_sha,
        code_checkpoint=actual_checkpoint, codeset_sha256=actual_codeset_sha)
    identity_sha256 = bytes2hex(sha256(codeunits(identity_text)))
    existing = validate_existing_outputs(
        output_root, selected_truth, prior_path;
        code_checkpoint=actual_checkpoint,
        codeset_sha256=actual_codeset_sha,
        input_set_sha256=actual_input_sha,
        campaign_identity_sha256=identity_sha256)
    initialize = get(ENV, "ROUND4_NOSIF_INITIALIZE_IDENTITY", "0")
    initialize in ("0", "1") || error(
        "ROUND4_NOSIF_INITIALIZE_IDENTITY must be 0 or 1")
    identity = joinpath(output_root, ".control", "campaign_identity.dat")
    if initialize == "1"
        identity = initialize_campaign_identity!(output_root, identity_text)
    elseif isfile(identity)
        read(identity, String) == identity_text || error(
            "existing round-4 campaign identity differs")
    end
    println("round-4 no-SIF preflight passed: states=$(join(states, ',')) " *
            "class=$expected_scene_class resumable_outputs=$existing " *
            "prior_sha256=$actual_prior_sha input_set_sha256=$actual_input_sha " *
            "codeset_sha256=$actual_codeset_sha identity=$identity")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
