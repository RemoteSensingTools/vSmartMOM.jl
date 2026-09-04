#!/usr/bin/env julia

"""
Fail-closed preflight for an isolated bottom-layer, no-SIF retrieval partition
using the opt-in tapered mapped-ACOS CO2 prior.

Required environment variables:

- `EXPECTED_STATES`: comma-separated state indices/ranges.
- `BOTTOM_RETRIEVAL_CAMPAIGN_ROOT`: bottom-layer truth campaign root.
- `RETRIEVAL_PRIOR_PATH`: generated tapered-prior NetCDF.
- `CURRENT_PRIOR_PATH`: approved untapered reference-prior NetCDF.
- `TAPERED_PRIOR_SHA256` and `CURRENT_PRIOR_SHA256`: caller-approved hashes
  that pin both priors.
- `RETRIEVAL_OUTPUT_ROOT`: isolated output root for this retrieval round.

This program validates only inputs and any already-present resumable outputs;
it does not run a forward model or write a retrieval product. With
`TAPERED_NOSIF_INITIALIZE_IDENTITY=1`, it atomically creates or verifies the
small campaign-identity receipt below the output root before workers start.
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

const EXPECTED_CO2_MODEL =
    "acos_mapped_tapered_vertical_correlation"
const EXPECTED_CO2_BASE_MODEL = "acos_mapped"
const EXPECTED_TAPER_TYPE = "nonstationary_ar1_schur"
const CAMPAIGN_ID =
    "bottom_layer_nosif_acos_mapped_tapered_vertical_correlation_v1"
const CAMPAIGN_IDENTITY_SCHEMA = 1
const EXPECTED_TAPER_RETENTION = [
    0.98 + (index - 1) * (0.55 - 0.98) / 10 for index in 1:11
]
const EXPECTED_ADJACENT_CORRELATION = [
    0.9599841058533383,
    0.9199699332712215,
    0.8724550511563143,
    0.8248973950782327,
    0.7881857693079238,
    0.7454300983633975,
    0.7008696412340173,
    0.6547195361487076,
    0.6024021854175606,
    0.5617510691168982,
    0.5013049539297996,
]

function parse_state_spec(value::AbstractString)
    states = Int[]
    for raw_token in split(value, ',')
        token = strip(raw_token)
        isempty(token) && error("EXPECTED_STATES contains an empty token")
        bounds = split(token, '-'; limit=2)
        first_state = parse(Int, first(bounds))
        last_state = length(bounds) == 1 ? first_state :
            parse(Int, last(bounds))
        first_state <= last_state || error(
            "descending state range is not allowed: $token")
        append!(states, first_state:last_state)
    end
    isempty(states) && error("EXPECTED_STATES selected no states")
    length(states) == length(unique(states)) || error(
        "EXPECTED_STATES contains duplicate state indices")
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

function validate_path_isolation(campaign_root, prior_path, output_root)
    active_prior = abspath(joinpath(
        campaign_root, "retrieval_setup", "apriori_states.nc"))
    active_output = abspath(joinpath(campaign_root, "retrievals"))
    prior_path != active_prior || error(
        "tapered retrievals may not use the active campaign prior $active_prior")
    output_root != active_output || error(
        "tapered retrievals may not use the active campaign output $active_output")
    startswith(output_root, active_output * Base.Filesystem.path_separator) &&
        error("tapered output must not be nested inside the active output root")
    return nothing
end

function validate_tapered_prior(prior_path;
                                expected_sha256,
                                reference_prior_path=joinpath(
                                    dirname(prior_path), "apriori_states.nc"),
                                reference_sha256)
    isfile(prior_path) || error(
        "missing tapered prior $prior_path; generate it with " *
        "CO2_COVARIANCE_MODEL=$EXPECTED_CO2_MODEL and a distinct APRIORI_OUTPUT")
    isfile(reference_prior_path) || error(
        "missing current approved reference prior $reference_prior_path")
    occursin(r"^[0-9a-f]{64}$", expected_sha256) || error(
        "tapered-prior SHA must be a lowercase 64-character SHA-256")
    occursin(r"^[0-9a-f]{64}$", reference_sha256) || error(
        "reference-prior SHA must be a lowercase 64-character SHA-256")
    actual_sha256 = file_sha256(prior_path)
    actual_sha256 == expected_sha256 || error(
        "tapered-prior SHA mismatch for $prior_path")
    actual_reference_sha256 = file_sha256(reference_prior_path)
    actual_reference_sha256 == reference_sha256 || error(
        "current reference-prior SHA mismatch for $reference_prior_path")
    NCDataset(prior_path) do dataset
        get(dataset.attrib, "apriori_complete", 0) == 1 || error(
            "prior is not marked apriori_complete=1: $prior_path")
        String(get(dataset.attrib, "co2_covariance_model", "")) ==
            EXPECTED_CO2_MODEL || error(
                "wrong CO2 covariance model in $prior_path")
        String(get(dataset.attrib, "co2_covariance_base_model", "")) ==
            EXPECTED_CO2_BASE_MODEL || error(
                "wrong CO2 covariance base model in $prior_path")
        String(get(dataset.attrib, "co2_correlation_taper_type", "")) ==
            EXPECTED_TAPER_TYPE || error(
                "wrong CO2 correlation taper in $prior_path")
        adjacent = Float64.(dataset[
            "co2_selected_adjacent_correlation"][:])
        isapprox(adjacent, EXPECTED_ADJACENT_CORRELATION;
                 atol=2e-14, rtol=0) || error(
            "stored tapered adjacent-layer correlations differ from the approved values")
        active_indices = Int.(dataset["active_parameter_index"][:])
        active_indices == vcat(1, collect(6:34)) || error(
            "tapered prior has the wrong 30-parameter active-state mapping")
        size(dataset["Sa_active"]) == (30, 30, 4) || error(
            "tapered prior does not contain four 30x30 active covariances")
        Float64(get(dataset.attrib, "aerosol_ln_aod_sigma", NaN)) == 0.75 ||
            error("tapered prior must retain sigma(ln AOD)=0.75")
        Float64(get(dataset.attrib,
                    "sif_wavelength_slope_prior_mw_m2_sr_nm2", NaN)) == 0.0 ||
            error("tapered prior must retain the zero-centered SIF slope")
        Float64(get(dataset.attrib,
                    "sif_wavelength_slope_sigma_mw_m2_sr_nm2", NaN)) ==
            0.002625 || error(
            "tapered prior must retain the approved SIF-slope sigma=0.002625")
        for order in ("p1", "p2"), band in
                ("o2a", "weak_co2", "strong_co2")
            name = "surface_$(order)_sigma_$(band)"
            Float64(get(dataset.attrib, name, NaN)) == 0.002 || error(
                "tapered prior must retain $name=0.002")
        end
        for surface_index in 1:4
            covariance = Float64.(dataset[
                "Sa_active"][:, :, surface_index])
            isposdef(Symmetric(covariance)) || error(
                "surface $surface_index tapered active covariance is not positive definite")
        end

        # The scientific change approved for this round is only the CO2
        # interlayer correlation.  Compare against the current campaign prior
        # directly so no aerosol, height, surface, SIF, state-mean, marginal
        # variance, or active-coordinate change can slip in through a builder
        # default.
        NCDataset(reference_prior_path) do reference
            selected_xa = Float64.(dataset["xa"][:, :])
            reference_xa = Float64.(reference["xa"][:, :])
            selected_xa == reference_xa || error(
                "tapered and current priors have different state means")
            selected_active = Int.(dataset["active_parameter_index"][:])
            reference_active = Int.(reference["active_parameter_index"][:])
            selected_active == reference_active || error(
                "tapered and current priors have different active coordinates")
            for variable in (
                    "active_mask", "prior_sigma", "co2_layer_center_height",
                    "co2_layer_center_pressure", "sif_wavelength_state",
                    "sif_wavelength_sigma")
                selected_values = Array(dataset[variable][:])
                reference_values = Array(reference[variable][:])
                isequal(selected_values, reference_values) || error(
                    "tapered and current priors differ in $variable")
            end
            for attribute in (
                    "surface_order", "parameter_names", "parameter_units",
                    "state_order", "aerosol_coordinate_transform")
                get(dataset.attrib, attribute, nothing) ==
                    get(reference.attrib, attribute, nothing) || error(
                    "tapered and current priors differ in $attribute")
            end

            active_co2_positions = findall(
                index -> index in 2:17, selected_active)
            length(active_co2_positions) == 12 || error(
                "expected exactly 12 active CO2 layers")
            taper = Matrix{Float64}(I, 12, 12)
            for row in 1:11
                retained = 1.0
                for column in row + 1:12
                    retained *= EXPECTED_TAPER_RETENTION[column - 1]
                    taper[row, column] = taper[column, row] = retained
                end
            end
            changed_co2_offdiagonal = false
            for surface_index in 1:4
                selected = Float64.(dataset["Sa"][:, :, surface_index])
                current = Float64.(reference["Sa"][:, :, surface_index])
                diag(selected) == diag(current) || error(
                    "tapered prior changed a marginal variance for surface $surface_index")
                difference = selected - current
                difference[2:17, 2:17] .= 0
                all(iszero, difference) || error(
                    "tapered prior changed a non-CO2 covariance for surface $surface_index")
                fixed_co2 = setdiff(collect(2:17), selected_active)
                selected[fixed_co2, 2:17] == current[fixed_co2, 2:17] ||
                    error("tapered prior changed fixed-upper CO2 covariance for surface $surface_index")
                selected[2:17, fixed_co2] == current[2:17, fixed_co2] ||
                    error("tapered prior changed fixed-upper CO2 covariance for surface $surface_index")

                selected_active_covariance = Float64.(dataset[
                    "Sa_active"][:, :, surface_index])
                reference_active_covariance = Float64.(reference[
                    "Sa_active"][:, :, surface_index])
                selected_active_covariance ==
                    selected[selected_active, selected_active] || error(
                    "tapered Sa_active is inconsistent for surface $surface_index")
                reference_active_covariance ==
                    current[reference_active, reference_active] || error(
                    "reference Sa_active is inconsistent for surface $surface_index")
                selected_co2 = selected_active_covariance[
                    active_co2_positions, active_co2_positions]
                reference_co2 = reference_active_covariance[
                    active_co2_positions, active_co2_positions]
                expected_co2 = reference_co2 .* taper
                for index in axes(expected_co2, 1)
                    expected_co2[index, index] = reference_co2[index, index]
                end
                isapprox(selected_co2, expected_co2;
                         atol=0.0, rtol=32eps(Float64)) || error(
                    "tapered prior has an incorrect non-adjacent CO2 covariance for surface $surface_index")
                changed_co2_offdiagonal |= selected_co2 != reference_co2
            end
            changed_co2_offdiagonal || error(
                "tapered prior contains no CO2 off-diagonal change")
        end
    end
    return (; prior_sha256=actual_sha256,
            reference_prior_sha256=actual_reference_sha256,
            prior_path=abspath(prior_path),
            reference_prior_path=abspath(reference_prior_path))
end

function existing_retrieval_products(output_root)
    isdir(output_root) || return String[]
    products = String[]
    for measurement_class in ("corrected", "uncorrected")
        directory = joinpath(output_root, measurement_class)
        isdir(directory) || continue
        append!(products, filter(
            path -> occursin(
                r"^retrieval_state\d{3}_perturbation\d{2}\.nc$",
                basename(path)),
            joinpath.(directory, readdir(directory))))
    end
    return products
end

function campaign_identity_text(campaign_root, truth_table, oco_root,
                                noise_root, prior_identity,
                                input_set_sha256,
                                snr_coefficients_sha256)
    fields = (
        "identity_schema" => string(CAMPAIGN_IDENTITY_SCHEMA),
        "campaign_id" => CAMPAIGN_ID,
        "sif_filter" => "off",
        "co2_covariance_model" => EXPECTED_CO2_MODEL,
        "prior_sha256" => prior_identity.prior_sha256,
        "prior_path" => prior_identity.prior_path,
        "reference_prior_sha256" => prior_identity.reference_prior_sha256,
        "reference_prior_path" => prior_identity.reference_prior_path,
        "truth_table_sha256" => file_sha256(truth_table),
        "no_sif_input_set_sha256" => input_set_sha256,
        "snr_coefficients_sha256" => snr_coefficients_sha256,
        "truth_table" => abspath(truth_table),
        "measurement_directory" => abspath(oco_root),
        "noise_directory" => abspath(noise_root),
        "campaign_root" => abspath(campaign_root),
    )
    return join(("$key=$value" for (key, value) in fields), '\n') * "\n"
end

"""Atomically initialize or verify the immutable local campaign identity."""
function initialize_campaign_identity!(output_root, expected;
                                       wait_seconds::Real=30)
    control = joinpath(output_root, ".control")
    identity = joinpath(control, "campaign_identity.dat")
    lock = joinpath(control, ".identity_lock")
    mkpath(control)
    if isfile(identity)
        read(identity, String) == expected || error(
            "existing tapered no-SIF campaign identity does not match current inputs")
        return identity
    end
    isempty(existing_retrieval_products(output_root)) || error(
        "tapered no-SIF output contains retrievals but no campaign identity")

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
            "campaign identity lock exists without a published identity: $lock")
        read(identity, String) == expected || error(
            "concurrent tapered no-SIF campaign identity differs")
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

function validate_existing_outputs(output_root, selected_truth, prior_path)
    isdir(output_root) || return 0
    checked = 0
    for truth in selected_truth, measurement_class in (:corrected, :uncorrected),
            perturbation in vcat(11, collect(1:10))
        path = joinpath(
            output_root, String(measurement_class),
            @sprintf("retrieval_state%03d_perturbation%02d.nc",
                     truth.state_index, perturbation))
        isfile(path) || continue
        prior = load_retrieval_prior(truth.surface; path=prior_path)
        NCDataset(path) do dataset
            get(dataset.attrib, "retrieval_complete", 0) == 1 || error(
                "partial retrieval output requires inspection before resume: $path")
            abspath(String(get(dataset.attrib, "source_apriori", ""))) ==
                prior_path || error(
                    "existing output was made with a different prior: $path")
            Float64.(dataset["a_priori_state"][:]) == prior.xa || error(
                "existing output embeds a different prior state: $path")
            Float64.(dataset["a_priori_covariance"][:, :]) == prior.Sa || error(
                "existing output embeds a different prior covariance: $path")
        end
        checked += 1
    end
    return checked
end

function truth_scene_path(truth, truth_root::AbstractString)
    directory = truth.aerosol_case == :none ? truth_root :
        joinpath(truth_root, "aerosol_chunked")
    return joinpath(
        directory, @sprintf("hiressim_%03d.nc", truth.state_index))
end

function no_sif_input_set_sha256(cases, truth_root, oco_root, noise_root)
    records = String[]
    selected = sort!(filter(case -> case.sif_case == :off, collect(cases));
                     by=case -> case.state_index)
    length(selected) == 40 || error(
        "bottom-layer campaign must contain exactly 40 no-SIF states")
    for truth in selected
        index = truth.state_index
        for (kind, path) in (
                ("truth", truth_scene_path(truth, truth_root)),
                ("measurement", joinpath(
                    oco_root, @sprintf("OCO2sims_%03d.nc", index))),
                ("noise", joinpath(
                    noise_root, @sprintf("OCO2noise_%03d.nc", index))))
            isfile(path) || error(
                "missing no-SIF $kind input for state $index: $path")
            push!(records, "$kind $index $(file_sha256(path))")
        end
    end
    return bytes2hex(sha256(codeunits(join(records, '\n') * "\n")))
end

function main()
    states = parse_state_spec(get(ENV, "EXPECTED_STATES", ""))
    campaign_root = required_path("BOTTOM_RETRIEVAL_CAMPAIGN_ROOT")
    prior_path = required_path("RETRIEVAL_PRIOR_PATH")
    reference_prior_path = required_path("CURRENT_PRIOR_PATH")
    output_root = required_path("RETRIEVAL_OUTPUT_ROOT")
    expected_prior_sha256 = required_sha256("TAPERED_PRIOR_SHA256")
    expected_reference_sha256 = required_sha256("CURRENT_PRIOR_SHA256")
    validate_path_isolation(campaign_root, prior_path, output_root)
    prior_identity = validate_tapered_prior(
        prior_path; expected_sha256=expected_prior_sha256,
        reference_prior_path,
        reference_sha256=expected_reference_sha256)

    truth_root = joinpath(campaign_root, "truth")
    truth_table = joinpath(truth_root, "true_states.dat")
    oco_root = joinpath(truth_root, "OCO_radiances")
    noise_root = joinpath(oco_root, "noise_covariances")
    truth_cases = read_truth_cases(truth_table)
    truth_by_index = Dict(truth.state_index => truth for truth in truth_cases)
    selected_truth = map(states) do state
        haskey(truth_by_index, state) || error(
            "state $state is absent from $truth_table")
        truth = truth_by_index[state]
        truth.sif_case == :off || error(
            "state $state is SIF-on; this launcher is permanently no-SIF")
        truth
    end

    coefficients_path = joinpath(
        @__DIR__, "instrument", "representative_snr_coefficients.nc")
    coefficients = NoiseValidation.read_representative_snr_coefficients(
        coefficients_path)
    for truth in selected_truth
        state = truth.state_index
        truth_path = truth_scene_path(truth, truth_root)
        radiance_path = joinpath(
            oco_root, @sprintf("OCO2sims_%03d.nc", state))
        noise_path = joinpath(
            noise_root, @sprintf("OCO2noise_%03d.nc", state))
        isfile(truth_path) || error("missing truth scene $truth_path")
        isfile(radiance_path) || error("missing OCO radiance $radiance_path")
        isfile(noise_path) || error("missing frozen-noise product $noise_path")
        RadianceValidation.validate_file(radiance_path, state)
        NoiseValidation.validate_file(
            noise_path, radiance_path, coefficients, state)
    end

    existing = validate_existing_outputs(
        output_root, selected_truth, prior_path)
    input_set_sha256 = no_sif_input_set_sha256(
        truth_cases, truth_root, oco_root, noise_root)
    identity_text = campaign_identity_text(
        campaign_root, truth_table, oco_root, noise_root, prior_identity,
        input_set_sha256, file_sha256(coefficients_path))
    initialize_identity = get(
        ENV, "TAPERED_NOSIF_INITIALIZE_IDENTITY", "0")
    initialize_identity in ("0", "1") || error(
        "TAPERED_NOSIF_INITIALIZE_IDENTITY must be 0 or 1")
    identity_was_initialized = initialize_identity == "1"
    identity = identity_was_initialized ?
        initialize_campaign_identity!(output_root, identity_text) :
        joinpath(output_root, ".control", "campaign_identity.dat")
    if initialize_identity == "0" && isfile(identity)
        read(identity, String) == identity_text || error(
            "existing tapered no-SIF campaign identity does not match current inputs")
    end
    println("tapered no-SIF preflight passed: states=$(join(states, ',')) " *
            "inputs=$(length(states)) resumable_outputs=$existing " *
            "prior=$prior_path prior_sha256=$(prior_identity.prior_sha256) " *
            "input_set_sha256=$input_set_sha256 " *
            "output=$output_root identity=$identity " *
            "identity_initialized=$identity_was_initialized")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
