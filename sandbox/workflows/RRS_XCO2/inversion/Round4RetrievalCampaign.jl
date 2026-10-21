module Round4RetrievalCampaign

using LinearAlgebra
using NCDatasets
using SHA

using ..Round4KnownSIF
using ..Round4SIFTruthConvention

export ROUND4_STATE_MODEL,
       ROUND4_JACOBIAN_FLAVOR,
       ROUND4_KNOWN_WAVELENGTH_NM,
       expected_active_to_full,
       resolve_round4_sif_map,
       validate_round4_prior,
       validate_round4_realization,
       validate_round4_output_root,
       round4_campaign_identity,
       round4_output_provenance,
       validate_round4_output_provenance,
       file_sha256

const ROUND4_STATE_MODEL = "round4_known_sif759"
const ROUND4_JACOBIAN_FLAVOR =
    "OCO_RRS_synth_round4_known_sif759_boundary_chain"
const ROUND4_KNOWN_WAVELENGTH_NM = 759.0

# The prior files index the original 34-coordinate state, whereas
# Round4KnownSIF indexes the already reduced 30-column OCO_RRS_synth basis.
# Full coordinates 2:5 (CO2 layers 1:4) are fixed in every retrieval.
const BASE_ACTIVE_TO_FULL = vcat(1, collect(6:32))
const MSIF_FULL_INDEX = 34

file_sha256(path::AbstractString) = open(path, "r") do stream
    bytes2hex(sha256(stream))
end

expected_active_to_full(map::Round4SIFMap) = map.sif_on ?
    vcat(BASE_ACTIVE_TO_FULL, MSIF_FULL_INDEX) : copy(BASE_ACTIVE_TO_FULL)

function _scene_sif_on(scene)
    hasproperty(scene, :sif_case) || throw(ArgumentError(
        "round-4 scene is missing sif_case"))
    return Symbol(getproperty(scene, :sif_case)) != :off
end

"""
    resolve_round4_sif_map(scenes; requested="auto", known_wavelength_nm=759)

Select one homogeneous round-4 SIF state model.  `auto` is accepted only when
all selected truth scenes agree.  An explicit `on` or `off` must still match
every scene; it cannot relabel truth data.
"""
function resolve_round4_sif_map(
        scenes;
        requested::AbstractString="auto",
        known_wavelength_nm::Real=ROUND4_KNOWN_WAVELENGTH_NM)
    wavelength = Float64(known_wavelength_nm)
    wavelength == ROUND4_KNOWN_WAVELENGTH_NM || throw(ArgumentError(
        "round 4 fixes SIF at exactly 759 nm; received $wavelength nm"))
    isempty(scenes) && throw(ArgumentError(
        "cannot select a round-4 SIF map from an empty scene set"))

    scene_modes = unique(_scene_sif_on(scene) for scene in scenes)
    length(scene_modes) == 1 || throw(ArgumentError(
        "one round-4 process cannot mix SIF-on and SIF-off truth scenes"))
    scene_on = only(scene_modes)
    normalized = lowercase(strip(requested))
    normalized in ("auto", "on", "off") || throw(ArgumentError(
        "ROUND4_SIF_MODE must be auto, on, or off"))
    requested_on = normalized == "auto" ? scene_on : normalized == "on"
    requested_on == scene_on || throw(ArgumentError(
        "ROUND4_SIF_MODE=$normalized disagrees with the selected truth scenes"))

    if scene_on
        convention = validate_round4_sif_truth_convention()
        return Round4SIFMap(true; Lν759=convention.diagnostic.Lnu)
    end
    return Round4SIFMap(false)
end

function _require_attribute(attributes, name::AbstractString, source)
    haskey(attributes, name) || error(
        "$source is missing required round-4 attribute '$name'")
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

"""Fail closed unless a generated prior is the exact round-4 SIF variant."""
function validate_round4_prior(prior, map::Round4SIFMap;
                               path::AbstractString)
    isfile(path) || throw(ArgumentError("missing round-4 prior: $path"))
    expected_full = expected_active_to_full(map)
    prior.active_to_full == expected_full || error(
        "round-4 prior active/full map is $(prior.active_to_full); expected " *
        string(expected_full))
    length(prior.xa) == active_state_count(map) || error(
        "round-4 prior state length disagrees with the selected SIF mode")
    size(prior.Sa) == (active_state_count(map), active_state_count(map)) ||
        error("round-4 prior covariance shape disagrees with the selected SIF mode")
    isposdef(Symmetric(prior.Sa)) || error(
        "round-4 active prior covariance is not positive definite")
    if map.sif_on
        last(prior.parameter_names) == "mSIF" || error(
            "SIF-on round-4 prior must retain mSIF as its final coordinate")
    else
        any(name -> name in ("SIF760", "mSIF"), prior.parameter_names) &&
            error("SIF-off round-4 prior must exclude both SIF coordinates")
    end

    source = "round-4 prior $path"
    NCDataset(path) do dataset
        attributes = dataset.attrib
        _require_equal(attributes, "retrieval_state_model",
                       ROUND4_STATE_MODEL, source)
        _require_equal(attributes, "round4_prior_model",
                       "round4_known_sif759_acos_mapped_tapered_vertical_correlation",
                       source)
        _require_equal(attributes, "round4_prior_definition_version", 1,
                       source)
        _require_equal(attributes, "round4_sif_case",
                       map.sif_on ? "on" : "off", source)
        _require_equal(attributes, "known_sif_wavelength_nm",
                       ROUND4_KNOWN_WAVELENGTH_NM, source)
        _require_equal(attributes,
                       "known_sif_Lnu_mW_m-2_sr-1_per_cm-1",
                       map.Lν759, source)
        expected_Llambda = map.sif_on ?
            lnu_to_llambda(map.Lν759, ROUND4_KNOWN_WAVELENGTH_NM) : 0.0
        _require_equal(attributes,
                       "known_sif_Llambda_mW_m-2_sr-1_nm-1",
                       expected_Llambda, source)
        _require_equal(attributes, "active_state_count",
                       active_state_count(map), source)
        _require_equal(attributes, "co2_covariance_model",
                       "acos_mapped_tapered_vertical_correlation", source)
        haskey(dataset, "active_core_parameter_index") || error(
            "$source is missing active_core_parameter_index")
        Int.(dataset["active_core_parameter_index"][:]) ==
            active_core_indices(map) || error(
            "$source has an incorrect active/core mapping")
        _require_equal(attributes, "round4_active_to_full",
                       join(expected_full, " "), source)
        _require_equal(attributes, "round4_active_to_core",
                       join(active_core_indices(map), " "), source)
        source_hash = String(_require_attribute(
            attributes, "source_prior_sha256", source))
        occursin(r"^[0-9a-f]{64}$", source_hash) || error(
            "$source has a malformed source_prior_sha256")
        !isempty(strip(String(_require_attribute(
            attributes, "source_prior_path", source)))) || error(
            "$source has an empty source_prior_path")
    end
    return prior
end

"""Validate that measurement/noise provenance matches the selected SIF map."""
function validate_round4_realization(map::Round4SIFMap, truth, realization;
                                     source::AbstractString="measurement")
    _scene_sif_on(truth) == map.sif_on || error(
        "$source SIF mode disagrees with the round-4 state map")
    if map.sif_on
        validate_round4_sif_provenance(
            realization.provenance; enabled=true, source)
    else
        iszero(map.Lν759) || error(
            "SIF-off round-4 map has a nonzero fixed SIF value")
    end
    return realization
end

"""Require a visibly isolated output namespace for the new retrieval model."""
function validate_round4_output_root(path::AbstractString)
    root = abspath(path)
    normalized = lowercase(replace(root, '-' => '_'))
    occursin("round4", normalized) && occursin("sif759", normalized) ||
        throw(ArgumentError(
            "RETRIEVAL_OUTPUT_ROOT must be an isolated round4/known-SIF759 " *
            "namespace; received $root"))
    return root
end

"""Read and strictly validate the portable production-campaign identity."""
function round4_campaign_identity(environment=ENV; required::Bool=true)
    specifications = (
        ("ROUND4_CODE_CHECKPOINT_SHA", "round4_code_checkpoint_sha", 40),
        ("ROUND4_CODESET_SHA256", "round4_codeset_sha256", 64),
        ("ROUND4_INPUT_SET_SHA256", "round4_input_set_sha256", 64),
        ("ROUND4_CAMPAIGN_IDENTITY_SHA256",
         "round4_campaign_identity_sha256", 64),
    )
    present = [haskey(environment, source) for (source, _, _) in specifications]
    if !any(present)
        required && error(
            "round-4 production identity variables are required")
        return Dict{String,Any}(
            "round4_campaign_identity_status" =>
                "not_provided_nonproduction")
    end
    all(present) || error(
        "round-4 campaign identity is incomplete; provide all four " *
        "ROUND4_* SHA variables")
    result = Dict{String,Any}(
        "round4_campaign_identity_status" => "complete")
    for (source, destination, length_expected) in specifications
        value = String(environment[source])
        occursin(Regex("^[0-9a-f]{$length_expected}\$"), value) || error(
            "$source must contain exactly $length_expected lowercase hex characters")
        result[destination] = value
    end
    return result
end

function round4_output_provenance(
        map::Round4SIFMap;
        prior_path::AbstractString,
        campaign_identity::AbstractDict=Dict{String,Any}(),
        convention_loader::Function=validate_round4_sif_truth_convention)
    # SIF-off retrievals must remain independent of the private corrected-SIF
    # template. Only the SIF-on model needs to evaluate that truth resource.
    convention = map.sif_on ? convention_loader() : nothing
    known_Llambda = map.sif_on ? convention.diagnostic.Llambda : 0.0
    provenance = Dict{String,Any}(
        "retrieval_state_model" => ROUND4_STATE_MODEL,
        "round4_known_sif_enabled" => 1,
        "round4_sif_case" => map.sif_on ? "on" : "off",
        "round4_known_sif_wavelength_nm" => ROUND4_KNOWN_WAVELENGTH_NM,
        "round4_known_sif_wavenumber_cm-1" => NU_759_CM1,
        "round4_core_sif_reference_wavenumber_cm-1" => NU_760_CM1,
        "round4_delta_nu_760_minus_759_cm-1" =>
            DELTA_NU_760_MINUS_759_CM1,
        "round4_known_sif_Lnu_mW_m-2_sr-1_per_cm-1" => map.Lν759,
        "round4_known_sif_Llambda_mW_m-2_sr-1_nm-1" => known_Llambda,
        "round4_known_sif_source" => map.sif_on ?
            "corrected-v2 full SIF truth template evaluated at 759 nm" :
            "exact zero for SIF-off truth",
        "round4_core_sif_reference_wavelength_nm" => 760.0,
        "round4_core_state_dimension" => CORE_STATE_COUNT,
        "round4_active_state_dimension" => active_state_count(map),
        "round4_active_to_core_parameter_index" =>
            join(active_core_indices(map), " "),
        "round4_active_to_full_parameter_index" =>
            join(expected_active_to_full(map), " "),
        "round4_sif_slope_status" => map.sif_on ?
            "active; prior-constrained" : "fixed to exact zero; zero variance",
        "round4_state_expansion" => map.sif_on ?
            "SIF760=Lnu759+mSIF*(nu760-nu759); core mSIF=active mSIF" :
            "SIF760=0; mSIF=0",
        "round4_jacobian_chain_rule" => map.sif_on ?
            "K_mSIF_round4=K_SIF760*(nu760-nu759)+K_mSIF_core" :
            "retain core columns 1:28; omit both fixed SIF columns",
        "round4_prior_sha256" => file_sha256(prior_path),
    )
    if !map.sif_on
        provenance["round4_canonical_mSIF_fixed_mean"] = 0.0
        provenance["round4_canonical_mSIF_fixed_variance"] = 0.0
        provenance["round4_canonical_mSIF_fixed_units"] =
            "mW m-2 sr-1 (cm-1)-2"
        provenance["round4_canonical_mSIF_excluded_from_active_solve"] = 1
    end
    for (name, value) in campaign_identity
        key = String(name)
        haskey(provenance, key) && error(
            "campaign identity would replace round-4 provenance key $key")
        provenance[key] = value
    end
    return provenance
end

"""Validate the round-4 identity used when deciding whether output is resumable."""
function validate_round4_output_provenance(attributes, map::Round4SIFMap;
                                           prior_path::AbstractString,
                                           campaign_identity::AbstractDict=
                                               Dict{String,Any}(),
                                           source::AbstractString="output")
    expected = round4_output_provenance(
        map; prior_path, campaign_identity)
    for (name, value) in expected
        _require_equal(attributes, name, value, source)
    end
    _require_equal(attributes, "jacobian_flavor",
                   ROUND4_JACOBIAN_FLAVOR, source)
    _require_equal(attributes, "state_dimension",
                   active_state_count(map), source)
    return true
end

end # module Round4RetrievalCampaign
