module Round7AltGrdRetrievalCampaign

#=
Round 7_altgrd: the completed Round 7 retrieval (fixed SIF, tight UTLS,
imperfectly LUT-corrected O2 observations), repeated with ONE change: the
vertical grid is held fixed in geometric altitude when the retrieved surface
pressure changes (`FixedAltitudeGrid`), instead of fixing the pressure
interfaces above the surface (`FixedPressureGrid`, Round 7).

Everything that defines the Round 7 inputs is reused unchanged: the sealed
observations, stored noise and covariance, priors, LUT-derived corrections,
settings and validators all come from `Round7RetrievalCampaign.jl`, which is a
byte-identical copy of the deployed Round 7 module. This module only adds the
new state-model name, its own code identity, an isolated output namespace,
and the extra provenance describing the grid convention.
=#

using SHA
using ..Round7RetrievalCampaign
import ..Round7RetrievalCampaign as R7

export ROUND7_ALTGRD_STATE_MODEL, ROUND7_ALTGRD_CODE_FILES,
       ROUND7_ALTGRD_REFERENCE_STATE_MODEL, ROUND7_ALTGRD_VERTICAL_GRID,
       round7_altgrd_identity_values, round7_altgrd_campaign_identity,
       validate_round7_altgrd_output_root, round7_altgrd_output_provenance,
       parse_perturbation_selection

const ROUND7_ALTGRD_STATE_MODEL = "round7_altgrd_fixed_sif_imperfect_correction"
const ROUND7_ALTGRD_REFERENCE_STATE_MODEL = R7.ROUND7_STATE_MODEL
const ROUND7_ALTGRD_VERTICAL_GRID = "FixedAltitudeGrid"
# The Round 7 module is part of the code identity: it defines the inputs.
const ROUND7_ALTGRD_CODE_FILES = ("Round7AltGrdRetrievalCampaign.jl",
    "Round7RetrievalCampaign.jl", "run_round7_altgrd_retrievals.jl",
    "test_round7_altgrd_retrieval_campaign.jl")
const _IDENTITY_NAMES = ("ROUND7_ALTGRD_CODE_CHECKPOINT_SHA",
    "ROUND7_ALTGRD_CODESET_SHA256", "ROUND7_ALTGRD_INPUT_SET_SHA256",
    "ROUND7_ALTGRD_CAMPAIGN_IDENTITY_SHA256",
    "ROUND7_ALTGRD_OBSERVATION_MANIFEST_SHA256")

_digest_lines(values) = bytes2hex(sha256(join(values, "\n") * "\n"))

function _hash(value, name; length=64)
    value isa AbstractString && occursin(Regex("^[0-9a-f]{$length}\$"), value) ||
        error("$name must be $length lowercase hexadecimal characters")
    return String(value)
end

"""
Portable deployment identity, using the Round 7 recipe with the Round
7_altgrd state model. `checkpoint` is the clean commit of this branch that
the pilot runs on; `code_paths` are `ROUND7_ALTGRD_CODE_FILES` in order;
`input_paths` are the Round 7 inputs in the Round 7 order followed by the six
spectroscopy tables (whose contents Round 7 did not hash).
"""
function round7_altgrd_identity_values(; checkpoint, base_tree_sha256, code_paths,
                                        input_paths, observation_manifest)
    _hash(checkpoint, "checkpoint"; length=40)
    _hash(base_tree_sha256, "base_tree_sha256")
    code = _digest_lines(vcat(base_tree_sha256, file_sha256.(code_paths)))
    inputs = _digest_lines(file_sha256.(input_paths))
    observations = file_sha256(observation_manifest)
    campaign = _digest_lines([ROUND7_ALTGRD_STATE_MODEL, checkpoint, code,
                              inputs, observations])
    return Dict(
        "ROUND7_ALTGRD_CODE_CHECKPOINT_SHA" => checkpoint,
        "ROUND7_ALTGRD_CODESET_SHA256" => code,
        "ROUND7_ALTGRD_INPUT_SET_SHA256" => inputs,
        "ROUND7_ALTGRD_CAMPAIGN_IDENTITY_SHA256" => campaign,
        "ROUND7_ALTGRD_OBSERVATION_MANIFEST_SHA256" => observations)
end

"""Read the exported identity and, when given, require the recomputed one."""
function round7_altgrd_campaign_identity(environment=ENV; expected=nothing)
    result = Dict{String,Any}("round7_altgrd_campaign_identity_status" => "complete")
    for name in _IDENTITY_NAMES
        haskey(environment, name) || error("$name is required")
        value = _hash(environment[name], name;
                      length=endswith(name, "CHECKPOINT_SHA") ? 40 : 64)
        isnothing(expected) || value == expected[name] ||
            error("$name differs from the actual code/input identity")
        result[lowercase(name)] = value
    end
    return result
end

"""
The output must be a new `round7_altgrd` tree ending in `retrievals_sif` or
`retrievals_nosif`; the completed Round 7 tree can never be selected.
"""
function validate_round7_altgrd_output_root(path; sif_on)
    root = R7.validate_round7_output_root(path; sif_on)
    occursin("round7_altgrd", lowercase(root)) || throw(ArgumentError(
        "Round 7_altgrd output must lie in an isolated round7_altgrd namespace"))
    occursin(ROUND7_ALTGRD_REFERENCE_STATE_MODEL * "/", root) && throw(ArgumentError(
        "refusing to write into the completed Round 7 output tree"))
    return root
end

"""
Round 7 provenance, relabelled with the Round 7_altgrd state model and
identity, plus an explicit record of the vertical-grid convention. Only
per-result facts belong here: resume validation compares every entry, so
launch-specific choices (member selection, GPU) are recorded at write time.
"""
function round7_altgrd_output_provenance(map; prior_path, observation,
                                         campaign_identity, truth_table_sha256)
    result = round7_output_provenance(map; prior_path, observation,
        campaign_identity, truth_table_sha256)
    result["retrieval_state_model"] = ROUND7_ALTGRD_STATE_MODEL
    merge!(result, Dict{String,Any}(
        "round7_altgrd_reference_state_model" => ROUND7_ALTGRD_REFERENCE_STATE_MODEL,
        "round7_altgrd_vertical_grid_convention" => ROUND7_ALTGRD_VERTICAL_GRID,
        "round7_altgrd_reference_vertical_grid_convention" => "FixedPressureGrid",
        "round7_altgrd_difference_from_reference" =>
            "vertical grid convention for surface pressure only",
        "round7_altgrd_psurf_profile_rule" =>
            "p_half(psurf) = p_half_ref * psurf / psurf_ref; layer T, q, VMR and geometric thickness fixed",
        "round7_altgrd_aerosol_profile_rule" =>
            "altitude-form LogNormal aerosol fixed in geometric altitude above the surface; zero surface-pressure tangent",
        "round7_altgrd_xco2_weights" =>
            "dry-air layer columns all scale with psurf; normalized weights independent of psurf"))
    return result
end

"""Parse a comma/space separated list of perturbation members, e.g. `1,11`."""
function parse_perturbation_selection(text::AbstractString)
    members = sort!(unique!(parse.(Int, split(text, r"[,\s]+"; keepempty=false))))
    isempty(members) && throw(ArgumentError("perturbation selection is empty"))
    all(in(1:11), members) || throw(ArgumentError("perturbations must lie in 1:11"))
    return members
end

end # module
