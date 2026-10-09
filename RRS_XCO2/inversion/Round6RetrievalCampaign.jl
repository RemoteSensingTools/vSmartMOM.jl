module Round6RetrievalCampaign

using LinearAlgebra
using NCDatasets
using SHA

using ..Round4SIFTruthConvention
using ..Round6FixedSIF

export ROUND6_STATE_MODEL,
       ROUND6_JACOBIAN_FLAVOR,
       ROUND6_KNOWN_WAVELENGTH_NM,
       expected_active_to_full,
       resolve_round6_sif_map,
       validate_round6_prior,
       validate_round6_realization,
       validate_round6_output_root,
       round6_campaign_identity,
       round6_output_provenance,
       validate_round6_output_provenance,
       file_sha256

const ROUND6_STATE_MODEL = "round6_fixed_sif"
const ROUND6_JACOBIAN_FLAVOR =
    "OCO_RRS_synth_round6_fixed_sif_boundary"
const ROUND6_KNOWN_WAVELENGTH_NM = 759.0
const ROUND6_ACTIVE_TO_FULL = vcat(1, collect(6:32))
const EXPECTED_PRIOR_MODEL =
    "round6_fixed_sif759_msif_standard_utls_acos_mapped_tapered_vertical_correlation"
const EXPECTED_CO2_MODEL = "acos_mapped_tapered_vertical_correlation"

file_sha256(path::AbstractString) = open(path, "r") do stream
    bytes2hex(sha256(stream))
end

expected_active_to_full(::Round6SIFMap) = copy(ROUND6_ACTIVE_TO_FULL)

function _scene_sif_on(scene)
    hasproperty(scene, :sif_case) || throw(ArgumentError(
        "round-6 scene is missing sif_case"))
    return Symbol(getproperty(scene, :sif_case)) != :off
end

"""Resolve one homogeneous truth subset to its fixed round-6 SIF map."""
function resolve_round6_sif_map(
        scenes;
        requested::AbstractString="auto",
        known_wavelength_nm::Real=ROUND6_KNOWN_WAVELENGTH_NM)
    Float64(known_wavelength_nm) == ROUND6_KNOWN_WAVELENGTH_NM ||
        throw(ArgumentError(
            "round 6 fixes SIF at exactly 759 nm"))
    isempty(scenes) && throw(ArgumentError(
        "cannot select a round-6 SIF map from an empty scene set"))
    modes = unique(_scene_sif_on(scene) for scene in scenes)
    length(modes) == 1 || throw(ArgumentError(
        "one round-6 process cannot mix SIF-on and SIF-off truth scenes"))
    truth_on = only(modes)
    normalized = lowercase(strip(requested))
    normalized in ("auto", "on", "off") || throw(ArgumentError(
        "ROUND6_SIF_MODE must be auto, on, or off"))
    requested_on = normalized == "auto" ? truth_on : normalized == "on"
    requested_on == truth_on || throw(ArgumentError(
        "ROUND6_SIF_MODE=$normalized disagrees with selected truth"))
    if truth_on
        convention = validate_round4_sif_truth_convention()
        return Round6SIFMap(true;
            Lnu759=convention.diagnostic.Lnu,
            mSIF=convention.mSIF)
    end
    return Round6SIFMap(false)
end

function _require_attribute(attributes, name::AbstractString, source)
    haskey(attributes, name) || error(
        "$source is missing required round-6 attribute '$name'")
    return attributes[name]
end

function _require_equal(attributes, name, expected, source)
    actual = _require_attribute(attributes, name, source)
    matches = expected isa Real && actual isa Real ?
        isapprox(Float64(actual), Float64(expected); atol=2e-14, rtol=2e-12) :
        string(actual) == string(expected)
    matches || error(
        "$source has $name=$(repr(actual)); expected $(repr(expected))")
    return actual
end

"""Fail closed unless `prior` is the requested 28-coordinate round-6 prior."""
function validate_round6_prior(prior, map::Round6SIFMap;
                               path::AbstractString)
    isfile(path) || throw(ArgumentError("missing round-6 prior: $path"))
    prior.active_to_full == expected_active_to_full(map) || error(
        "round-6 prior has an incorrect active/full map")
    length(prior.xa) == ROUND6_STATE_COUNT || error(
        "round-6 prior must contain 28 active coordinates")
    size(prior.Sa) == (ROUND6_STATE_COUNT, ROUND6_STATE_COUNT) || error(
        "round-6 active covariance must be 28 by 28")
    isposdef(Symmetric(prior.Sa)) || error(
        "round-6 active covariance is not positive definite")
    any(name -> name in ("SIF760", "mSIF"), prior.parameter_names) &&
        error("round-6 prior must exclude both fixed SIF coordinates")

    source = "round-6 prior $path"
    NCDataset(path, "r") do dataset
        attributes = dataset.attrib
        _require_equal(attributes, "retrieval_state_model",
                       ROUND6_STATE_MODEL, source)
        _require_equal(attributes, "round6_prior_model",
                       EXPECTED_PRIOR_MODEL, source)
        _require_equal(attributes, "round6_prior_definition_version", 1,
                       source)
        _require_equal(attributes, "round6_sif_case",
                       map.sif_on ? "on" : "off", source)
        _require_equal(attributes, "known_sif_wavelength_nm",
                       ROUND6_KNOWN_WAVELENGTH_NM, source)
        _require_equal(attributes,
                       "known_sif_Lnu_mW_m-2_sr-1_per_cm-1",
                       map.Lnu759, source)
        _require_equal(attributes,
                       "fixed_mSIF_mW_m-2_sr-1_per_cm-2",
                       map.mSIF, source)
        _require_equal(attributes, "active_state_count",
                       ROUND6_STATE_COUNT, source)
        _require_equal(attributes, "co2_covariance_model",
                       EXPECTED_CO2_MODEL, source)
        _require_equal(attributes, "round6_active_to_full",
                       join(expected_active_to_full(map), " "), source)
        _require_equal(attributes, "round6_active_to_core",
                       join(active_core_indices(map), " "), source)
        _require_equal(attributes, "stratospheric_aerosol_sigma_scale",
                       1.0, source)
        Int.(dataset["active_core_parameter_index"][:]) ==
            active_core_indices(map) || error(
                "$source has an incorrect active/core mapping")
        Float64.(dataset["aerosol_ln_aod_sigma_by_species"][:]) ≈
            [0.75, 0.75, 0.75] || error(
                "$source has incorrect aerosol AOD sigmas")
        Float64.(dataset["aerosol_ln_z0_sigma_by_species"][:]) ≈
            [0.10, 0.10, 0.10] || error(
                "$source has incorrect aerosol-height sigmas")
        source_hash = String(_require_attribute(
            attributes, "source_prior_sha256", source))
        occursin(r"^[0-9a-f]{64}$", source_hash) || error(
            "$source has a malformed source_prior_sha256")
        source_path = get(ENV, "ROUND6_SOURCE_PRIOR",
                          String(_require_attribute(attributes, "source_prior_path", source)))
        file_sha256(source_path) == source_hash || error("round-6 source prior hash differs")
        NCDataset(source_path, "r") do original
            dataset["xa"][1:32, :] == original["xa"][1:32, :] ||
                error("round 6 changed a non-SIF prior mean")
            dataset["Sa"][1:32, 1:32, :] == original["Sa"][1:32, 1:32, :] ||
                error("round 6 changed a non-SIF prior covariance")
        end
        all(iszero, dataset["Sa"][33:34, :, :]) &&
            all(iszero, dataset["Sa"][:, 33:34, :]) || error("round-6 SIF covariance is not zero")
    end
    return prior
end

"""Validate truth/noise ownership against the fixed SIF selection."""
function validate_round6_realization(map::Round6SIFMap, truth, realization;
                                     source::AbstractString="measurement")
    _scene_sif_on(truth) == map.sif_on || error(
        "$source SIF mode disagrees with the round-6 state map")
    if map.sif_on
        validate_round4_sif_provenance(
            realization.provenance; enabled=true, source)
    else
        iszero(map.Lnu759) && iszero(map.mSIF) || error(
            "SIF-off round-6 map has nonzero fixed SIF")
    end
    return realization
end

function validate_round6_output_root(path::AbstractString)
    root = abspath(path)
    normalized = lowercase(replace(root, '-' => '_'))
    occursin("round6", normalized) && occursin("fixed_sif", normalized) ||
        throw(ArgumentError(
            "RETRIEVAL_OUTPUT_ROOT must be an isolated round6/fixed-SIF namespace"))
    return root
end

function round6_campaign_identity(environment=ENV; required::Bool=true)
    specifications = (
        ("ROUND6_CODE_CHECKPOINT_SHA", "round6_code_checkpoint_sha", 40),
        ("ROUND6_CODESET_SHA256", "round6_codeset_sha256", 64),
        ("ROUND6_INPUT_SET_SHA256", "round6_input_set_sha256", 64),
        ("ROUND6_CAMPAIGN_IDENTITY_SHA256",
         "round6_campaign_identity_sha256", 64),
    )
    present = [haskey(environment, source)
               for (source, _, _) in specifications]
    if !any(present)
        required && error("round-6 production identity variables are required")
        return Dict{String,Any}(
            "round6_campaign_identity_status" => "not_provided_nonproduction")
    end
    all(present) || error("round-6 campaign identity is incomplete")
    result = Dict{String,Any}(
        "round6_campaign_identity_status" => "complete")
    for (source, destination, expected_length) in specifications
        value = String(environment[source])
        occursin(Regex("^[0-9a-f]{$expected_length}\$"), value) || error(
            "$source must contain $expected_length lowercase hex characters")
        result[destination] = value
    end
    return result
end

function round6_output_provenance(
        map::Round6SIFMap;
        prior_path::AbstractString,
        campaign_identity::AbstractDict=Dict{String,Any}())
    provenance = Dict{String,Any}(
        "retrieval_state_model" => ROUND6_STATE_MODEL,
        "round6_sif_case" => map.sif_on ? "on" : "off",
        "round6_known_sif_wavelength_nm" => ROUND6_KNOWN_WAVELENGTH_NM,
        "round6_known_sif_Lnu_mW_m-2_sr-1_per_cm-1" => map.Lnu759,
        "round6_fixed_mSIF_mW_m-2_sr-1_per_cm-2" => map.mSIF,
        "round6_core_SIF760_mW_m-2_sr-1_per_cm-1" =>
            expand_round6_state(map, zeros(ROUND6_STATE_COUNT))[
                CORE_SIF760_INDEX],
        "round6_sif_parameter_status" =>
            "SIF759 and mSIF fixed; both excluded from active solve",
        "round6_core_state_dimension" => CORE_STATE_COUNT,
        "round6_active_state_dimension" => ROUND6_STATE_COUNT,
        "round6_active_to_core_parameter_index" =>
            join(active_core_indices(map), " "),
        "round6_active_to_full_parameter_index" =>
            join(expected_active_to_full(map), " "),
        "round6_stratospheric_aerosol_sigma_scale" => 1.0,
        "round6_non_sif_prior_status" => "unchanged from round3/round4 source prior",
        "round6_co2_prior_status" =>
            "unchanged from round4 source prior",
        "round6_prior_sha256" => file_sha256(prior_path),
    )
    for (name, value) in campaign_identity
        key = String(name)
        haskey(provenance, key) && error(
            "campaign identity would replace round-6 provenance key $key")
        provenance[key] = value
    end
    return provenance
end

function validate_round6_output_provenance(
        attributes, map::Round6SIFMap;
        prior_path::AbstractString,
        campaign_identity::AbstractDict=Dict{String,Any}(),
        source::AbstractString="output")
    expected = round6_output_provenance(
        map; prior_path, campaign_identity)
    for (name, value) in expected
        _require_equal(attributes, name, value, source)
    end
    _require_equal(attributes, "jacobian_flavor",
                   ROUND6_JACOBIAN_FLAVOR, source)
    _require_equal(attributes, "state_dimension",
                   ROUND6_STATE_COUNT, source)
    return true
end

end # module Round6RetrievalCampaign
