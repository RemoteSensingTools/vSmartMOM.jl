module GattacaRound4SIFReadiness

"""
Fail-closed Gattaca2 gate for the round-4 SIF-on retrieval campaign.

The existing corrected-v2 truth publication remains read-only. This module
reuses its validated publication barrier, then independently pins the new
29-coordinate known-SIF759 prior, code/input hashes, and a private output
namespace that cannot collide with round 3.
"""

using LinearAlgebra
using NCDatasets
using Printf
using SHA

include(joinpath(@__DIR__, "GattacaTaperedSIFReadiness.jl"))
include(joinpath(@__DIR__, "Round4SIFTruthConvention.jl"))

using .GattacaTaperedSIFReadiness
using .Round4SIFTruthConvention

export CAMPAIGN_ID,
       EXPECTED_SIF_STATES,
       Round4ReadinessPaths,
       legacy_input_repo_root,
       required_output_root,
       required_prior_path,
       required_summary_path,
       file_sha256,
       validate_checkout_separation,
       validate_output_isolation,
       validate_prior_identity,
       validate_release_barrier,
       codeset_sha256,
       prepare_readiness!,
       validate_state_outputs,
       main

const CAMPAIGN_ID =
    "bottom_layer_round4_known_sif759_sif_on_" *
    "acos_mapped_tapered_vertical_correlation_v1"
const PRIOR_FILENAME =
    "apriori_states_round4_known_sif759_on_" *
    "acos_mapped_tapered_vertical_correlation.nc"
const SUMMARY_FILENAME = replace(PRIOR_FILENAME, r"\.nc$" => ".dat")
const SOURCE_PRIOR_FILENAME =
    "source_apriori_states_acos_mapped_tapered_vertical_correlation.nc"
const EXPECTED_SIF_STATES = copy(
    GattacaTaperedSIFReadiness.EXPECTED_SIF_STATES)
const EXPECTED_ACTIVE_TO_FULL = vcat(1, collect(6:32), 34)
const EXPECTED_ACTIVE_TO_CORE = vcat(collect(1:28), 30)
const EXPECTED_CO2_MODEL = "acos_mapped_tapered_vertical_correlation"
const EXPECTED_STATE_MODEL = "round4_known_sif759"
const EXPECTED_PRIOR_MODEL =
    "round4_known_sif759_acos_mapped_tapered_vertical_correlation"
const EXPECTED_JACOBIAN_FLAVOR =
    "OCO_RRS_synth_round4_known_sif759_boundary_chain"
const EXPECTED_SLOPE_MEAN = 0.0
const EXPECTED_SLOPE_SIGMA = 0.002625
const EXPECTED_NATIVE_SLOPE_MEAN = -7.34275712494799e-7
const EXPECTED_NATIVE_SLOPE_SIGMA = 8.780708772523119e-6
const IDENTITY_SCHEMA = 1

struct Round4ReadinessPaths
    repo_root::String
    private_root::String
    full_truth_root::String
    restart_root::String
    bottom_campaign_root::String
    bottom_truth_table::String
    measurement_directory::String
    noise_directory::String
    source_prior_path::String
    prior_path::String
    summary_path::String
    output_root::String
    stokes_coefficient_path::String
    scene_components_path::String
    sif_template_path::String
end

resolved_target(path::AbstractString) =
    GattacaTaperedSIFReadiness.resolved_target(path)
contains_path(parent::AbstractString, child::AbstractString) =
    GattacaTaperedSIFReadiness.contains_path(parent, child)

function _required_environment_path(name::AbstractString)
    value = get(ENV, name, "")
    isempty(value) && error(
        "$name must be set explicitly for the isolated round-4 checkout")
    return value
end

required_result_root(private_root::AbstractString) = resolved_target(joinpath(
    private_root, "results", CAMPAIGN_ID))
required_output_root(private_root::AbstractString) = resolved_target(joinpath(
    private_root, "results", CAMPAIGN_ID, "retrievals"))
required_prior_path(private_root::AbstractString) = resolved_target(joinpath(
    private_root, "results", CAMPAIGN_ID, "retrieval_setup", PRIOR_FILENAME))
required_summary_path(private_root::AbstractString) = resolved_target(joinpath(
    private_root, "results", CAMPAIGN_ID, "retrieval_setup", SUMMARY_FILENAME))
required_source_prior_path(private_root::AbstractString) = resolved_target(joinpath(
    private_root, "results", CAMPAIGN_ID,
    "retrieval_setup", SOURCE_PRIOR_FILENAME))

function Round4ReadinessPaths(;
        repo_root=_required_environment_path("RRS_REPO"),
        private_root=get(ENV, "RRS_PRIVATE_ROOT",
                         joinpath(homedir(), "RRS_XCO2_private")),
        full_truth_root=_required_environment_path("FULL_COLUMN_TRUTH_ROOT"),
        restart_root=get(ENV, "SIF_RESTART_ROOT",
                         joinpath(private_root, "results", ".sif_v2_restart")),
        bottom_campaign_root=
            _required_environment_path("BOTTOM_XCO2_CAMPAIGN_ROOT"),
        bottom_truth_table=get(
            ENV, "RETRIEVAL_TRUTH_TABLE",
            joinpath(bottom_campaign_root, "truth", "true_states.dat")),
        measurement_directory=get(
            ENV, "RETRIEVAL_MEASUREMENT_DIR",
            joinpath(bottom_campaign_root, "truth", "OCO_radiances")),
        noise_directory=get(
            ENV, "RETRIEVAL_NOISE_DIR",
            joinpath(measurement_directory, "noise_covariances")),
        source_prior_path=get(
            ENV, "ROUND4_SOURCE_PRIOR_PATH",
            required_source_prior_path(private_root)),
        prior_path=get(
            ENV, "RETRIEVAL_PRIOR_PATH", required_prior_path(private_root)),
        summary_path=get(
            ENV, "ROUND4_PRIOR_SUMMARY_PATH",
            required_summary_path(private_root)),
        output_root=get(
            ENV, "RETRIEVAL_OUTPUT_ROOT", required_output_root(private_root)),
        stokes_coefficient_path=
            _required_environment_path("RETRIEVAL_STOKES_COEFFICIENT_PATH"),
        scene_components_path=
            _required_environment_path("RETRIEVAL_SCENE_COMPONENTS_PATH"),
        sif_template_path=
            _required_environment_path("RRS_XCO2_SIF_TEMPLATE_PATH"))
    return Round4ReadinessPaths(
        resolved_target(repo_root), resolved_target(private_root),
        resolved_target(full_truth_root), resolved_target(restart_root),
        resolved_target(bottom_campaign_root), resolved_target(bottom_truth_table),
        resolved_target(measurement_directory), resolved_target(noise_directory),
        resolved_target(source_prior_path), resolved_target(prior_path),
        resolved_target(summary_path), resolved_target(output_root),
        resolved_target(stokes_coefficient_path),
        resolved_target(scene_components_path),
        resolved_target(sif_template_path))
end

function _git_toplevel(path::AbstractString, description::AbstractString)
    isdir(path) || error("missing $description: $path")
    command = `git -C $path rev-parse --show-toplevel`
    top = try
        readchomp(command)
    catch exception
        error("$description is not inside a Git checkout: $path ($exception)")
    end
    return realpath(top)
end

"Return the canonical data-bearing round-3 checkout inferred from its campaign root."
function legacy_input_repo_root(paths::Round4ReadinessPaths)
    return _git_toplevel(paths.bottom_campaign_root,
                         "bottom-layer campaign root")
end

"""
Require two distinct checkouts: round-4 code in `RRS_REPO`, and immutable
round-3/full-column inputs in the checkout named by the two explicit input
roots.  Exact canonical layouts prevent a symlink or an environment fallback
from silently redirecting either job while round 3 is still active.
"""
function validate_checkout_separation(paths::Round4ReadinessPaths)
    code_repo = _git_toplevel(paths.repo_root, "round-4 code root")
    code_repo == paths.repo_root || error(
        "RRS_REPO must name the canonical top level of the round-4 checkout; " *
        "received $(paths.repo_root), top level is $code_repo")

    legacy_repo = legacy_input_repo_root(paths)
    (code_repo != legacy_repo &&
     !contains_path(code_repo, legacy_repo) &&
     !contains_path(legacy_repo, code_repo)) || error(
        "round 4 must use a separate, non-nested checkout from the active " *
        "round-3 input checkout; code=$code_repo legacy=$legacy_repo")
    (!contains_path(code_repo, paths.private_root) &&
     !contains_path(paths.private_root, code_repo) &&
     !contains_path(legacy_repo, paths.private_root) &&
     !contains_path(paths.private_root, legacy_repo)) || error(
        "RRS_PRIVATE_ROOT must be disjoint from both Git checkouts")

    expected_bottom = resolved_target(joinpath(
        legacy_repo, "RRS_XCO2", "bottom_layer_XCO2_retrievals"))
    expected_full = resolved_target(joinpath(
        legacy_repo, "RRS_XCO2", "truth_map"))
    paths.bottom_campaign_root == expected_bottom || error(
        "BOTTOM_XCO2_CAMPAIGN_ROOT must be $expected_bottom; received " *
        paths.bottom_campaign_root)
    paths.full_truth_root == expected_full || error(
        "FULL_COLUMN_TRUTH_ROOT must be $expected_full; received " *
        paths.full_truth_root)
    isdir(paths.full_truth_root) || error(
        "missing full-column truth root: $(paths.full_truth_root)")

    expected_table = resolved_target(joinpath(
        expected_bottom, "truth", "true_states.dat"))
    expected_measurements = resolved_target(joinpath(
        expected_bottom, "truth", "OCO_radiances"))
    expected_noise = resolved_target(joinpath(
        expected_measurements, "noise_covariances"))
    paths.bottom_truth_table == expected_table || error(
        "RETRIEVAL_TRUTH_TABLE must remain the canonical legacy campaign table")
    paths.measurement_directory == expected_measurements || error(
        "RETRIEVAL_MEASUREMENT_DIR must remain the canonical legacy measurement directory")
    paths.noise_directory == expected_noise || error(
        "RETRIEVAL_NOISE_DIR must remain the canonical legacy noise directory")

    expected_stokes = resolved_target(joinpath(
        legacy_repo, "RRS_XCO2", "inversion", "instrument",
        "representative_stokes_coefficients.nc"))
    expected_components = resolved_target(joinpath(
        expected_bottom, "truth", "scene_components.dat"))
    expected_sif_template = resolved_target(joinpath(
        legacy_repo, "src", "SIF_emission", "sif-spectra.csv"))
    paths.stokes_coefficient_path == expected_stokes || error(
        "RETRIEVAL_STOKES_COEFFICIENT_PATH must be $expected_stokes")
    paths.scene_components_path == expected_components || error(
        "RETRIEVAL_SCENE_COMPONENTS_PATH must be $expected_components")
    paths.sif_template_path == expected_sif_template || error(
        "RRS_XCO2_SIF_TEMPLATE_PATH must be $expected_sif_template")
    for (description, path) in (
            ("representative Stokes coefficients", expected_stokes),
            ("bottom-layer scene components", expected_components),
            ("SIF spectral template", expected_sif_template))
        isfile(path) || error("missing $description: $path")
    end

    for (description, path) in (
            ("round-4 output", paths.output_root),
            ("round-4 prior", paths.prior_path),
            ("round-4 prior summary", paths.summary_path),
            ("byte-exact round-4 source-prior copy", paths.source_prior_path))
        contains_path(code_repo, path) && error(
            "$description must not be stored in the round-4 Git checkout: $path")
        contains_path(legacy_repo, path) && error(
            "$description must not be stored in the data-bearing round-3 Git checkout: $path")
    end
    contains_path(paths.full_truth_root, paths.output_root) && error(
        "round-4 output must not be written below the full-column truth root")
    contains_path(paths.bottom_campaign_root, paths.output_root) && error(
        "round-4 output must not be written below the bottom-layer campaign root")
    return (; code_repo, legacy_repo, expected_bottom, expected_full)
end

function _round4_truth_convention(paths::Round4ReadinessPaths)
    state = Round4SIFTruthConvention.RRSXCO2Common.campaign_sif_state(
        sif_template_path=paths.sif_template_path)
    return validate_round4_sif_truth_convention(
        round4_sif_truth_convention(; state))
end

file_sha256(path::AbstractString) = open(path, "r") do stream
    bytes2hex(sha256(stream))
end

function _require_hash(value::AbstractString, name::AbstractString,
                       length_expected::Int=64)
    occursin(Regex("^[0-9a-f]{$length_expected}\$"), value) || error(
        "$name must contain exactly $length_expected lowercase hexadecimal characters")
    return value
end

function _require_file(path, description)
    isfile(path) || error("missing $description: $path")
    return path
end

function _require_attribute(attributes, name, path)
    haskey(attributes, name) || error("$path is missing attribute $name")
    return attributes[name]
end

function _require_exact_attribute(attributes, name, expected, path)
    actual = _require_attribute(attributes, name, path)
    matches = expected isa Real && actual isa Real ?
        isapprox(Float64(actual), Float64(expected);
                 atol=4eps(max(abs(Float64(expected)), 1.0)),
                 rtol=4eps(Float64)) : string(actual) == string(expected)
    matches || error("$path has $name=$(repr(actual)); expected $(repr(expected))")
    return actual
end

"Require the one immutable private namespace assigned to round 4."
function validate_output_isolation(paths::Round4ReadinessPaths)
    expected = required_output_root(paths.private_root)
    paths.output_root == expected || error(
        "RETRIEVAL_OUTPUT_ROOT must be $expected; received $(paths.output_root)")
    contains_path(paths.private_root, paths.output_root) || error(
        "round-4 output must remain below RRS_PRIVATE_ROOT")
    contains_path(paths.repo_root, paths.output_root) && error(
        "round-4 output must not be inside the Git checkout")
    legacy_repo = legacy_input_repo_root(paths)
    contains_path(legacy_repo, paths.output_root) && error(
        "round-4 output must not be inside the data-bearing round-3 checkout")
    old = GattacaTaperedSIFReadiness.required_output_root(paths.private_root)
    paths.output_root != old || error("round 4 may not reuse the round-3 output")
    return expected
end

function _validate_prior_isolation(paths::Round4ReadinessPaths)
    paths.prior_path == required_prior_path(paths.private_root) || error(
        "RETRIEVAL_PRIOR_PATH must be the dedicated round-4 SIF-on prior path")
    paths.summary_path == required_summary_path(paths.private_root) || error(
        "ROUND4_PRIOR_SUMMARY_PATH must be the dedicated round-4 summary path")
    paths.source_prior_path == required_source_prior_path(paths.private_root) ||
        error("ROUND4_SOURCE_PRIOR_PATH must be the dedicated round-4 source-prior copy")
    for path in (paths.prior_path, paths.summary_path, paths.source_prior_path)
        contains_path(paths.private_root, path) || error(
            "private prior input is outside RRS_PRIVATE_ROOT: $path")
        contains_path(paths.repo_root, path) && error(
            "private prior input is inside the Git checkout: $path")
        contains_path(legacy_input_repo_root(paths), path) && error(
            "private prior input is inside the data-bearing round-3 checkout: $path")
    end
    return nothing
end

function _all_values(dataset, name)
    haskey(dataset, name) || error("NetCDF file is missing variable $name")
    return Array(nomissing(Array(dataset[name]), NaN))
end

"Validate the exact 29-D production prior and its approved round-3 source."
function validate_prior_identity(paths::Round4ReadinessPaths;
                                 prior_sha256::AbstractString,
                                 summary_sha256::AbstractString,
                                 source_prior_sha256::AbstractString)
    _validate_prior_isolation(paths)
    _require_hash(prior_sha256, "ROUND4_PRIOR_SHA256")
    _require_hash(summary_sha256, "ROUND4_PRIOR_SUMMARY_SHA256")
    _require_hash(source_prior_sha256, "ROUND4_SOURCE_PRIOR_SHA256")
    _require_file(paths.prior_path, "round-4 SIF-on prior")
    _require_file(paths.summary_path, "round-4 SIF-on prior summary")
    _require_file(paths.source_prior_path,
                  "byte-exact round-4 source-prior copy")
    basename(paths.prior_path) == PRIOR_FILENAME || error(
        "round-4 prior must use filename $PRIOR_FILENAME")
    basename(paths.summary_path) == SUMMARY_FILENAME || error(
        "round-4 summary must use filename $SUMMARY_FILENAME")
    file_sha256(paths.prior_path) == prior_sha256 || error(
        "round-4 prior SHA-256 mismatch")
    file_sha256(paths.summary_path) == summary_sha256 || error(
        "round-4 prior-summary SHA-256 mismatch")
    file_sha256(paths.source_prior_path) == source_prior_sha256 || error(
        "round-4 source-prior-copy SHA-256 mismatch")

    convention = _round4_truth_convention(paths)
    known_Lnu = Float64(convention.diagnostic.Lnu)
    known_Llambda = Float64(convention.diagnostic.Llambda)
    delta_nu = 1.0e7 / 760.0 - 1.0e7 / 759.0

    source_xa, source_Sa = NCDataset(paths.source_prior_path, "r") do source
        _require_exact_attribute(source.attrib, "co2_covariance_model",
                                 EXPECTED_CO2_MODEL, paths.source_prior_path)
        (Float64.(_all_values(source, "xa")),
         Float64.(_all_values(source, "Sa")))
    end

    result = NCDataset(paths.prior_path, "r") do prior
        attributes = prior.attrib
        for (name, expected) in (
                "apriori_complete" => 1,
                "retrieval_state_model" => EXPECTED_STATE_MODEL,
                "round4_prior_model" => EXPECTED_PRIOR_MODEL,
                "round4_prior_definition_version" => 1,
                "round4_sif_case" => "on",
                "known_sif_wavelength_nm" => 759.0,
                "known_sif_Lnu_mW_m-2_sr-1_per_cm-1" => known_Lnu,
                "known_sif_Llambda_mW_m-2_sr-1_nm-1" => known_Llambda,
                "core_sif_reference_wavelength_nm" => 760.0,
                "active_state_count" => 29,
                "co2_covariance_model" => EXPECTED_CO2_MODEL,
                "source_prior_sha256" => source_prior_sha256,
                "round4_slope_prior_mean_mW_m-2_sr-1_nm-2" =>
                    EXPECTED_SLOPE_MEAN,
                "round4_slope_prior_sigma_mW_m-2_sr-1_nm-2" =>
                    EXPECTED_SLOPE_SIGMA,
                "round4_native_mSIF_prior_mean_mW_m-2_sr-1_per_cm-2" =>
                    EXPECTED_NATIVE_SLOPE_MEAN,
                "round4_native_mSIF_prior_sigma_mW_m-2_sr-1_per_cm-2" =>
                    EXPECTED_NATIVE_SLOPE_SIGMA)
            _require_exact_attribute(attributes, name, expected,
                                     paths.prior_path)
        end
        active = Int.(_all_values(prior, "active_parameter_index"))
        active_core = Int.(_all_values(prior, "active_core_parameter_index"))
        active == EXPECTED_ACTIVE_TO_FULL || error(
            "round-4 prior has the wrong active/full mapping")
        active_core == EXPECTED_ACTIVE_TO_CORE || error(
            "round-4 prior has the wrong active/core mapping")
        xa = Float64.(_all_values(prior, "xa"))
        Sa = Float64.(_all_values(prior, "Sa"))
        Sa_active = Float64.(_all_values(prior, "Sa_active"))
        size(xa) == (34, 4) || error("round-4 xa must have shape 34x4")
        size(Sa) == (34, 34, 4) || error("round-4 Sa must have shape 34x34x4")
        size(Sa_active) == (29, 29, 4) || error(
            "round-4 Sa_active must have shape 29x29x4")
        xa[1:32, :] == source_xa[1:32, :] || error(
            "round-4 prior changed a non-SIF mean")
        Sa[1:32, 1:32, :] == source_Sa[1:32, 1:32, :] || error(
            "round-4 prior changed a non-SIF covariance")
        all(iszero, Sa[33:34, 1:32, :]) || error(
            "round-4 prior couples SIF to a non-SIF coordinate")
        all(iszero, Sa[1:32, 33:34, :]) || error(
            "round-4 prior couples non-SIF to SIF")

        native_sigma = Float64(_require_attribute(attributes,
            "round4_native_mSIF_prior_sigma_mW_m-2_sr-1_per_cm-2",
            paths.prior_path))
        expected_sif_covariance = native_sigma^2 .* [
            delta_nu^2 delta_nu
            delta_nu   1.0
        ]
        for surface in 1:4
            isapprox(xa[33, surface] - delta_nu * xa[34, surface],
                     known_Lnu; atol=3e-18, rtol=4eps(Float64)) || error(
                "surface $surface prior mean violates exact Lnu759")
            isapprox(Sa[33:34, 33:34, surface], expected_sif_covariance;
                     atol=3e-22, rtol=16eps(Float64)) || error(
                "surface $surface has the wrong one-dimensional SIF covariance")
            Sa_active[:, :, surface] == Sa[active, active, surface] || error(
                "surface $surface Sa_active does not equal its full-state block")
            isposdef(Symmetric(Sa_active[:, :, surface])) || error(
                "surface $surface active covariance is not positive definite")
        end
        return (; xa, Sa, Sa_active, active, active_core)
    end

    summary = read(paths.summary_path, String)
    occursin("Exact SIF anchor: lambda=759.0 nm", summary) || error(
        "round-4 prior summary does not identify the exact 759-nm anchor")
    occursin("Source prior SHA-256: $source_prior_sha256", summary) || error(
        "round-4 prior summary names the wrong source prior")
    occursin("Gaussian conditioning", summary) || error(
        "round-4 prior summary omits the production slope-semantics warning")
    return merge(result, (;
        sha256=prior_sha256,
        summary_sha256,
        source_prior_sha256,
        path=paths.prior_path,
        summary_path=paths.summary_path,
        source_path=paths.source_prior_path,
        known_Lnu,
        known_Llambda,
    ))
end

function _legacy_paths(paths::Round4ReadinessPaths)
    return GattacaTaperedSIFReadiness.ReadinessPaths(
        legacy_input_repo_root(paths), paths.private_root, paths.full_truth_root,
        paths.restart_root, paths.bottom_campaign_root,
        paths.bottom_truth_table, paths.measurement_directory,
        paths.noise_directory, paths.prior_path, paths.output_root)
end

"Reuse the complete, already tested corrected-v2 publication barrier."
function validate_release_barrier(
        paths::Round4ReadinessPaths;
        full_validator=GattacaTaperedSIFReadiness.validate_published_sif_release,
        bottom_validator=
            GattacaTaperedSIFReadiness.validate_published_bottom_sif_release,
        input_validator=GattacaTaperedSIFReadiness.validate_bottom_sif_inputs)
    legacy = _legacy_paths(paths)
    full = full_validator(legacy)
    bottom = bottom_validator(legacy, full)
    inputs = input_validator(legacy, full, bottom)
    inputs.states == EXPECTED_SIF_STATES || error(
        "corrected-v2 barrier did not validate the canonical 40 SIF-on states")
    return (; full, bottom, inputs)
end

function _codeset_paths(repo_root)
    paths = String[]
    for relative in ("Project.toml", "Manifest.toml",
                     "RRS_XCO2/config/oco_grass_3aerosol.yaml")
        path = joinpath(repo_root, relative)
        _require_file(path, "round-4 code-set input")
        push!(paths, path)
    end
    for relative_root in ("src", "ext", "RRS_XCO2/scripts")
        root = joinpath(repo_root, relative_root)
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
            "RRS_XCO2/inversion/GattacaTaperedSIFReadiness.jl",
            "RRS_XCO2/inversion/GattacaRound4SIFReadiness.jl",
            "RRS_XCO2/inversion/instrument/SyntheticOCO2.jl",
            "RRS_XCO2/inversion/run_round4_known_sif_retrievals.jl",
            "RRS_XCO2/inversion/retrieval_setup/build_round4_known_sif_apriori.jl",
            "RRS_XCO2/inversion/gattaca_round4_known_sif759_retrievals.sbatch",
            "RRS_XCO2/inversion/submit_gattaca_round4_known_sif759_retrievals.sh")
        path = joinpath(repo_root, relative)
        _require_file(path, "round-4 code-set input")
        push!(paths, path)
    end
    return sort!(unique(paths))
end

function codeset_sha256(repo_root::AbstractString)
    records = ["$(relpath(path, repo_root)) $(file_sha256(path))"
               for path in _codeset_paths(repo_root)]
    return bytes2hex(sha256(codeunits(join(records, '\n') * "\n")))
end

function _input_set_sha256(paths, prior, release)
    records = (
        "campaign $CAMPAIGN_ID",
        "prior $(prior.sha256)",
        "prior_summary $(prior.summary_sha256)",
        "source_prior $(prior.source_prior_sha256)",
        "known_Lnu759 $(repr(prior.known_Lnu))",
        "full_release $(release.full.validation_receipt_sha256)",
        "bottom_release $(release.bottom.receipt_sha256)",
        "bottom_release_inputs $(release.bottom.input_set_sha256)",
        "validated_sif_inputs $(release.inputs.input_set_sha256)",
        "bottom_truth_table $(release.inputs.table_sha256)",
        "representative_stokes_coefficients " *
            file_sha256(paths.stokes_coefficient_path),
        "bottom_scene_components " * file_sha256(paths.scene_components_path),
        "sif_spectral_template " * file_sha256(paths.sif_template_path),
    )
    return bytes2hex(sha256(codeunits(join(records, '\n') * "\n")))
end

function _identity_base(paths, prior, release, checkpoint, codeset, input_set)
    legacy_repo = legacy_input_repo_root(paths)
    fields = (
        "identity_schema" => string(IDENTITY_SCHEMA),
        "campaign_id" => CAMPAIGN_ID,
        "source_checkpoint_sha" => checkpoint,
        "codeset_sha256" => codeset,
        "input_set_sha256" => input_set,
        "prior_sha256" => prior.sha256,
        "prior_summary_sha256" => prior.summary_sha256,
        "source_prior_sha256" => prior.source_prior_sha256,
        "known_sif_wavelength_nm" => "759",
        "known_sif_Lnu" => repr(prior.known_Lnu),
        "active_state_count" => "29",
        "active_to_full" => join(EXPECTED_ACTIVE_TO_FULL, ','),
        "active_to_core" => join(EXPECTED_ACTIVE_TO_CORE, ','),
        "sif_release_validation_sha256" =>
            release.full.validation_receipt_sha256,
        "bottom_sif_release_receipt_sha256" =>
            release.bottom.receipt_sha256,
        "bottom_sif_input_set_sha256" => release.inputs.input_set_sha256,
        "measurement_classes" => "corrected,uncorrected",
        "perturbation_order" => "11,1,2,3,4,5,6,7,8,9,10",
        "round4_code_repo_root" => paths.repo_root,
        "legacy_input_repo_root" => legacy_repo,
        "bottom_campaign_root" => paths.bottom_campaign_root,
        "full_column_truth_root" => paths.full_truth_root,
        "representative_stokes_coefficients" =>
            paths.stokes_coefficient_path,
        "representative_stokes_coefficients_sha256" =>
            file_sha256(paths.stokes_coefficient_path),
        "bottom_scene_components" => paths.scene_components_path,
        "bottom_scene_components_sha256" =>
            file_sha256(paths.scene_components_path),
        "sif_spectral_template" => paths.sif_template_path,
        "sif_spectral_template_sha256" =>
            file_sha256(paths.sif_template_path),
        "output_root" => paths.output_root,
    )
    return join(("$name=$value" for (name, value) in fields), '\n') * "\n"
end

function _existing_products(output_root)
    isdir(output_root) || return String[]
    products = String[]
    for class in ("corrected", "uncorrected")
        directory = joinpath(output_root, class)
        isdir(directory) || continue
        append!(products, filter(path -> occursin(
            r"^retrieval_state\d{3}_perturbation\d{2}\.nc$", basename(path)),
            joinpath.(directory, readdir(directory))))
    end
    return products
end

function _initialize_identity!(paths, text; wait_seconds::Real=60)
    control = joinpath(paths.output_root, ".control")
    identity = joinpath(control, "campaign_identity.dat")
    lock = joinpath(control, ".identity_lock")
    mkpath(control)
    if isfile(identity)
        read(identity, String) == text || error(
            "existing round-4 campaign identity differs")
        return identity
    end
    isempty(_existing_products(paths.output_root)) || error(
        "round-4 output contains products but no campaign identity")
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
            "round-4 identity lock exists without a published identity")
        read(identity, String) == text || error(
            "concurrent round-4 campaign identity differs")
        return identity
    end
    temporary = identity * ".tmp.$(getpid())"
    try
        write(temporary, text)
        mv(temporary, identity)
    finally
        isfile(temporary) && rm(temporary)
        isdir(lock) && rm(lock)
    end
    return identity
end

function _write_identity_environment!(paths, checkpoint, codeset,
                                      input_set, identity_sha256)
    destination = joinpath(dirname(paths.output_root),
                           "retrieval_setup", "campaign_identity.env")
    mkpath(dirname(destination))
    text = join((
        "ROUND4_CODE_CHECKPOINT_SHA=$checkpoint",
        "ROUND4_CODESET_SHA256=$codeset",
        "ROUND4_INPUT_SET_SHA256=$input_set",
        "ROUND4_CAMPAIGN_IDENTITY_SHA256=$identity_sha256",
    ), '\n') * "\n"
    if isfile(destination)
        read(destination, String) == text || error(
            "existing round-4 identity environment differs")
        return destination
    end
    temporary = destination * ".tmp.$(getpid())"
    write(temporary, text)
    mv(temporary, destination)
    return destination
end

function prepare_readiness!(paths::Round4ReadinessPaths;
                            prior_sha256::AbstractString,
                            summary_sha256::AbstractString,
                            source_prior_sha256::AbstractString,
                            checkpoint::AbstractString,
                            barrier_keywords...)
    validate_checkout_separation(paths)
    validate_output_isolation(paths)
    _require_hash(checkpoint, "ROUND4_CODE_CHECKPOINT_SHA", 40)
    prior = validate_prior_identity(paths;
        prior_sha256, summary_sha256, source_prior_sha256)
    release = validate_release_barrier(paths; barrier_keywords...)
    codeset = codeset_sha256(paths.repo_root)
    input_set = _input_set_sha256(paths, prior, release)
    base = _identity_base(
        paths, prior, release, checkpoint, codeset, input_set)
    identity_sha256 = bytes2hex(sha256(codeunits(base)))
    identity_text = base * "campaign_identity_sha256=$identity_sha256\n"
    identity = _initialize_identity!(paths, identity_text)
    environment = _write_identity_environment!(
        paths, checkpoint, codeset, input_set, identity_sha256)
    return (; prior, release, codeset, input_set, identity_sha256,
            identity, environment)
end

function _read_environment(path)
    values = Dict{String,String}()
    for line in eachline(_require_file(path, "campaign identity environment"))
        fields = split(strip(line), '='; limit=2)
        length(fields) == 2 || error("malformed identity environment: $path")
        values[fields[1]] = fields[2]
    end
    return values
end

"Validate one completed state after the smoke or production run."
function validate_state_outputs(paths::Round4ReadinessPaths, state::Integer;
                                perturbations=vcat(11, collect(1:10)))
    validate_checkout_separation(paths)
    validate_output_isolation(paths)
    state in EXPECTED_SIF_STATES || error("state $state is not a SIF-on state")
    all(index -> index in 1:11, perturbations) || error(
        "perturbations must lie in 1:11")
    identity_environment = _read_environment(joinpath(
        dirname(paths.output_root), "retrieval_setup", "campaign_identity.env"))
    prior_sha = file_sha256(paths.prior_path)
    prior_record = NCDataset(paths.prior_path, "r") do prior
        (;
            surfaces=split(String(prior.attrib["surface_order"])),
            xa=Float64.(_all_values(prior, "xa")),
            Sa_active=Float64.(_all_values(prior, "Sa_active")),
        )
    end
    checked = 0
    for class in ("corrected", "uncorrected"), perturbation in perturbations
        path = joinpath(paths.output_root, class,
            @sprintf("retrieval_state%03d_perturbation%02d.nc",
                     state, perturbation))
        _require_file(path, "round-4 retrieval output")
        NCDataset(path, "r") do dataset
            attributes = dataset.attrib
            for (name, expected) in (
                    "retrieval_complete" => 1,
                    "truth_state_index" => state,
                    "perturbation_index" => perturbation,
                    "measurement_class" => class,
                    "retrieval_state_model" => EXPECTED_STATE_MODEL,
                    "round4_sif_case" => "on",
                    "round4_known_sif_wavelength_nm" => 759.0,
                    "round4_known_sif_Lnu_mW_m-2_sr-1_per_cm-1" =>
                        _round4_truth_convention(paths).diagnostic.Lnu,
                    "state_dimension" => 29,
                    "jacobian_flavor" => EXPECTED_JACOBIAN_FLAVOR,
                    "round4_prior_sha256" => prior_sha,
                    "round4_code_checkpoint_sha" => identity_environment[
                        "ROUND4_CODE_CHECKPOINT_SHA"],
                    "round4_codeset_sha256" => identity_environment[
                        "ROUND4_CODESET_SHA256"],
                    "round4_input_set_sha256" => identity_environment[
                        "ROUND4_INPUT_SET_SHA256"],
                    "round4_campaign_identity_sha256" => identity_environment[
                        "ROUND4_CAMPAIGN_IDENTITY_SHA256"])
                _require_exact_attribute(attributes, name, expected, path)
            end
            Int.(_all_values(dataset, "active_core_parameter_index")) ==
                EXPECTED_ACTIVE_TO_CORE || error(
                "$path has the wrong active/core mapping")
            surface_name = String(_require_attribute(
                attributes, "surface", path))
            surface_index = findfirst(==(surface_name), prior_record.surfaces)
            isnothing(surface_index) && error(
                "$path names an unknown prior surface $surface_name")
            Float64.(_all_values(dataset, "a_priori_state")) ==
                prior_record.xa[EXPECTED_ACTIVE_TO_FULL, surface_index] ||
                error("$path embeds the wrong round-4 prior mean")
            Float64.(_all_values(dataset, "a_priori_covariance")) ==
                prior_record.Sa_active[:, :, surface_index] || error(
                "$path embeds the wrong round-4 prior covariance")
            if perturbation == 11
                all(iszero, Float64.(_all_values(
                    dataset, "injected_measurement_noise"))) || error(
                    "$path noiseless perturbation contains injected noise")
            end
        end
        checked += 1
    end
    return checked
end

function main()
    action = get(ENV, "GATTACA_ROUND4_ACTION", "prepare")
    paths = Round4ReadinessPaths()
    if action == "prepare"
        readiness = prepare_readiness!(paths;
            prior_sha256=get(ENV, "ROUND4_PRIOR_SHA256", ""),
            summary_sha256=get(ENV, "ROUND4_PRIOR_SUMMARY_SHA256", ""),
            source_prior_sha256=get(ENV, "ROUND4_SOURCE_PRIOR_SHA256", ""),
            checkpoint=get(ENV, "ROUND4_CODE_CHECKPOINT_SHA", ""))
        println("Gattaca round-4 SIF-on readiness: PASSED")
        println("campaign_id=$CAMPAIGN_ID")
        println("known_Lnu759=$(readiness.prior.known_Lnu)")
        println("active_state_count=29")
        println("codeset_sha256=$(readiness.codeset)")
        println("input_set_sha256=$(readiness.input_set)")
        println("campaign_identity_sha256=$(readiness.identity_sha256)")
        println("identity_environment=$(readiness.environment)")
        return
    end
    action == "validate-output" || error(
        "GATTACA_ROUND4_ACTION must be prepare or validate-output")
    state = parse(Int, get(ENV, "EXPECTED_STATE", ""))
    phase = get(ENV, "GATTACA_ROUND4_PHASE", "")
    perturbations = phase == "smoke" ? [11] :
        phase == "production" ? vcat(11, collect(1:10)) :
        error("GATTACA_ROUND4_PHASE must be smoke or production")
    count = validate_state_outputs(paths, state; perturbations)
    println("Gattaca round-4 output validation: PASSED ($count files)")
end

end # module GattacaRound4SIFReadiness

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    GattacaRound4SIFReadiness.main()
end
