module GattacaTaperedSIFReadiness

"""
Fail-closed readiness checks for the isolated Gattaca bottom-layer SIF-on
retrieval campaign that uses the vertically tapered ACOS CO2 covariance.

This module deliberately does not generate truth, instrument products, noise,
or priors.  It proves that those separately published inputs form one coherent
version-2 SIF chain before a GPU retrieval may start.  In particular, the
existence of a NetCDF file is never treated as readiness: its metadata and
content provenance must trace back to the published direct-SIF truth release.
"""

using NCDatasets
using LinearAlgebra
using Printf
using SHA

include(joinpath(@__DIR__, "RetrievalCases.jl"))
using .RetrievalCases

export CAMPAIGN_ID,
       EXPECTED_PRIOR_MODEL,
       EXPECTED_PRIOR_BASE_MODEL,
       EXPECTED_SIF_STATES,
       ReadinessPaths,
       file_sha256,
       required_output_root,
       required_prior_path,
       validate_output_isolation,
       validate_prior_isolation,
       validate_prior_identity,
       validate_published_sif_release,
       validate_published_bottom_sif_release,
       validate_bottom_sif_inputs,
       validate_readiness,
       initialize_campaign_identity!,
       main

const CAMPAIGN_ID =
    "bottom_layer_sif_acos_mapped_tapered_vertical_correlation_v1"
const EXPECTED_PRIOR_MODEL = "acos_mapped_tapered_vertical_correlation"
const EXPECTED_PRIOR_BASE_MODEL = "acos_mapped"
const EXPECTED_SIF_STATES = reduce(vcat, [
    collect(first_state:(first_state + 4))
    for first_state in (6, 16, 26, 36, 46, 56, 66, 76)
])
const EXPECTED_MEASUREMENT_COUNT = 2742
const IDENTITY_SCHEMA = 1
const EXPECTED_ACTIVE_PARAMETER_INDEX = vcat(1, collect(6:34))
const FULL_CO2_PARAMETER_INDEX = 2:17
const EXPECTED_AEROSOL_LN_AOD_SIGMA = 0.75
const EXPECTED_SURFACE_P1_SIGMA = 0.002
const EXPECTED_SURFACE_P2_SIGMA = 0.002
const EXPECTED_SIF_REFERENCE_RADIANCE = 0.1
const EXPECTED_SIF_SLOPE_MEAN = 0.0
const EXPECTED_SIF_SLOPE_SIGMA = 0.002625
const EXPECTED_TAPER_RETENTION = [
    0.98 + (index - 1) * (0.55 - 0.98) / 10 for index in 1:11
]

const BOTTOM_PROVENANCE_ATTRIBUTES = (
    "campaign",
    "source_state_table",
    "column_xco2_ppm",
    "background_co2_ppm",
    "bottom_co2_ppm",
    "bottom_co2_index",
    "bottom_co2_layer_index",
    "co2_profile_order",
    "co2_profile_definition",
    "co2_profile_ppm",
    "bottom_layer_dry_column_fraction",
    "state_table_sha256",
    "full_column_source_state",
    "full_column_source_scene",
    "full_column_source_sha256",
    "full_column_state_table_sha256",
    "bottom_co2_truth_mode",
    "co2_truth_reuse_source",
    "producer_script",
    "producer_script_sha256",
    "bottom_layer_truth_complete",
)

struct ReadinessPaths
    repo_root::String
    private_root::String
    full_truth_root::String
    restart_root::String
    bottom_campaign_root::String
    bottom_truth_table::String
    measurement_directory::String
    noise_directory::String
    prior_path::String
    output_root::String
end

"Resolve symlinks through the longest existing prefix of a prospective path."
function resolved_target(path::AbstractString)
    target = normpath(abspath(path))
    suffix = String[]
    cursor = target
    while !ispath(cursor)
        parent = dirname(cursor)
        parent == cursor && error("cannot resolve path ancestry for $target")
        push!(suffix, basename(cursor))
        cursor = parent
    end
    resolved = realpath(cursor)
    for component in reverse(suffix)
        resolved = joinpath(resolved, component)
    end
    return normpath(resolved)
end

function contains_path(parent::AbstractString, child::AbstractString)
    relative = relpath(child, parent)
    return relative == "." ||
           !(relative == ".." ||
             startswith(relative, "..$(Base.Filesystem.path_separator)"))
end

"Match a canonical receipt destination without binding it to one host's root."
function _has_path_tail(path::AbstractString, tail::AbstractString)
    path_parts = splitpath(normpath(path))
    tail_parts = splitpath(normpath(tail))
    length(path_parts) >= length(tail_parts) || return false
    return path_parts[end-length(tail_parts)+1:end] == tail_parts
end

required_output_root(private_root::AbstractString) = resolved_target(joinpath(
    private_root, "results", CAMPAIGN_ID, "retrievals"))
required_prior_path(private_root::AbstractString) = resolved_target(joinpath(
    private_root, "results", CAMPAIGN_ID, "retrieval_setup",
    "apriori_states_$(EXPECTED_PRIOR_MODEL).nc"))

function ReadinessPaths(;
        repo_root=get(ENV, "RRS_REPO", joinpath(homedir(), "code", "uni_vSmartMOM")),
        private_root=get(ENV, "RRS_PRIVATE_ROOT", joinpath(homedir(), "RRS_XCO2_private")),
        full_truth_root=get(ENV, "FULL_COLUMN_TRUTH_ROOT",
                            joinpath(repo_root, "RRS_XCO2", "truth_map")),
        restart_root=get(ENV, "SIF_RESTART_ROOT",
                         joinpath(private_root, "results", ".sif_v2_restart")),
        bottom_campaign_root=get(
            ENV, "BOTTOM_XCO2_CAMPAIGN_ROOT",
            joinpath(repo_root, "RRS_XCO2", "bottom_layer_XCO2_retrievals")),
        bottom_truth_table=get(
            ENV, "RETRIEVAL_TRUTH_TABLE",
            joinpath(bottom_campaign_root, "truth", "true_states.dat")),
        measurement_directory=get(
            ENV, "RETRIEVAL_MEASUREMENT_DIR",
            joinpath(bottom_campaign_root, "truth", "OCO_radiances")),
        noise_directory=get(
            ENV, "RETRIEVAL_NOISE_DIR",
            joinpath(measurement_directory, "noise_covariances")),
        prior_path=get(
            ENV, "RETRIEVAL_PRIOR_PATH",
            required_prior_path(private_root)),
        output_root=get(
            ENV, "RETRIEVAL_OUTPUT_ROOT",
            required_output_root(private_root)))
    return ReadinessPaths(
        resolved_target(repo_root), resolved_target(private_root),
        resolved_target(full_truth_root), resolved_target(restart_root),
        resolved_target(bottom_campaign_root), resolved_target(bottom_truth_table),
        resolved_target(measurement_directory), resolved_target(noise_directory),
        resolved_target(prior_path), resolved_target(output_root))
end

file_sha256(path::AbstractString) = open(path, "r") do stream
    bytes2hex(sha256(stream))
end

function _require_file(path, description)
    isfile(path) || error("missing $description: $path")
    return path
end

function _require_attr(attributes, key, path)
    haskey(attributes, key) || error("$path is missing required attribute $key")
    return attributes[key]
end

function _metadata(path::AbstractString)
    _require_file(path, "receipt")
    result = Dict{String,String}()
    for raw in eachline(path)
        line = strip(raw)
        isempty(line) && continue
        startswith(line, '#') && (line = strip(line[2:end]))
        fields = split(line; limit=2)
        length(fields) == 2 || continue
        result[fields[1]] = fields[2]
    end
    return result
end

function _required_metadata(metadata, key, path)
    haskey(metadata, key) || error("$path is missing receipt field $key")
    return metadata[key]
end

function _release_rows(path::AbstractString)
    rows = Dict{Int,NamedTuple}()
    for raw in eachline(path)
        line = strip(raw)
        (isempty(line) || startswith(line, '#')) && continue
        fields = split(line)
        length(fields) == 5 || error("malformed SIF release row in $path: $line")
        index = try
            parse(Int, fields[1])
        catch
            error("invalid SIF release state index in $path: $(fields[1])")
        end
        haskey(rows, index) && error("duplicate state $index in $path")
        rows[index] = (
            class=fields[2],
            staged_sha256=fields[3],
            archived_sha256=fields[4],
            destination=fields[5],
        )
    end
    return rows
end

function _bottom_release_rows(path::AbstractString)
    rows = Dict{Int,NamedTuple}()
    for raw in eachline(path)
        line = strip(raw)
        (isempty(line) || startswith(line, '#')) && continue
        fields = split(line)
        length(fields) == 4 || error(
            "malformed bottom-layer SIF release row in $path: $line")
        index = try
            parse(Int, fields[1])
        catch
            error("invalid bottom-layer release state index in $path: $(fields[1])")
        end
        haskey(rows, index) && error("duplicate state $index in $path")
        hashes = fields[2:4]
        all(hash -> occursin(r"^[0-9a-f]{64}$", hash), hashes) || error(
            "invalid SHA-256 in bottom-layer release row $index")
        rows[index] = (
            truth_sha256=hashes[1],
            measurement_sha256=hashes[2],
            noise_sha256=hashes[3],
        )
    end
    return rows
end

"Require the one immutable private result namespace assigned to this campaign."
function validate_output_isolation(paths::ReadinessPaths)
    expected = required_output_root(paths.private_root)
    paths.output_root == expected || error(
        "RETRIEVAL_OUTPUT_ROOT must be the isolated tapered-SIF path $expected; " *
        "received $(paths.output_root)")
    contains_path(paths.private_root, paths.output_root) || error(
        "tapered-SIF output must remain below RRS_PRIVATE_ROOT")
    contains_path(paths.repo_root, paths.output_root) && error(
        "refusing to write retrieval products inside the Git checkout")
    return expected
end

"Require the generated tapered prior to share this campaign's private root."
function validate_prior_isolation(paths::ReadinessPaths)
    expected = required_prior_path(paths.private_root)
    paths.prior_path == expected || error(
        "RETRIEVAL_PRIOR_PATH must be the isolated tapered-SIF path $expected; " *
        "received $(paths.prior_path)")
    contains_path(paths.private_root, paths.prior_path) || error(
        "tapered-SIF prior must remain below RRS_PRIVATE_ROOT")
    contains_path(paths.repo_root, paths.prior_path) && error(
        "refusing to use a generated tapered prior inside the Git checkout")
    return expected
end

function _require_exact_numeric_attr(attributes, key, expected, path)
    value = Float64(_require_attr(attributes, key, path))
    value == Float64(expected) || error(
        "$path has $key=$value; expected exactly $expected")
    return value
end

function _validate_approved_nonco2_attributes(dataset, path)
    attributes = dataset.attrib
    _require_exact_numeric_attr(
        attributes, "aerosol_ln_aod_sigma",
        EXPECTED_AEROSOL_LN_AOD_SIGMA, path)
    for band in ("o2a", "weak_co2", "strong_co2")
        _require_exact_numeric_attr(
            attributes, "surface_p1_sigma_$band",
            EXPECTED_SURFACE_P1_SIGMA, path)
        _require_exact_numeric_attr(
            attributes, "surface_p2_sigma_$band",
            EXPECTED_SURFACE_P2_SIGMA, path)
    end
    _require_exact_numeric_attr(
        attributes, "sif_slope_reference_radiance_mw_m2_sr_nm",
        EXPECTED_SIF_REFERENCE_RADIANCE, path)
    _require_exact_numeric_attr(
        attributes, "sif_wavelength_slope_prior_mw_m2_sr_nm2",
        EXPECTED_SIF_SLOPE_MEAN, path)
    _require_exact_numeric_attr(
        attributes, "sif_wavelength_slope_sigma_mw_m2_sr_nm2",
        EXPECTED_SIF_SLOPE_SIGMA, path)
    return nothing
end

function _prior_snapshot(path)
    NCDataset(path) do dataset
        _validate_approved_nonco2_attributes(dataset, path)
        active = Int.(dataset["active_parameter_index"][:])
        active == EXPECTED_ACTIVE_PARAMETER_INDEX || error(
            "$path has the wrong active-state mapping; expected full-state " *
            "indices 1 and 6:34")
        xa = Float64.(dataset["xa"][:, :])
        Sa = Float64.(dataset["Sa"][:, :, :])
        Sa_active = Float64.(dataset["Sa_active"][:, :, :])
        size(xa) == (34, 4) || error("$path xa must have shape 34x4")
        size(Sa) == (34, 34, 4) || error("$path Sa must have shape 34x34x4")
        size(Sa_active) == (30, 30, 4) || error(
            "$path Sa_active must have shape 30x30x4")
        for surface in axes(Sa_active, 3)
            Sa_active[:, :, surface] == Sa[active, active, surface] || error(
                "$path Sa_active does not exactly equal its active full-state block")
            isposdef(Symmetric(Sa_active[:, :, surface])) || error(
                "$path has a non-positive-definite active covariance")
        end
        required_arrays = Dict(
            "active_mask" => Array(dataset["active_mask"][:]),
            "prior_sigma" => Float64.(dataset["prior_sigma"][:, :]),
            "co2_layer_center_height" =>
                Float64.(dataset["co2_layer_center_height"][:]),
            "co2_layer_center_pressure" =>
                Float64.(dataset["co2_layer_center_pressure"][:]),
            "sif_wavelength_state" =>
                Float64.(dataset["sif_wavelength_state"][:, :]),
            "sif_wavelength_sigma" =>
                Float64.(dataset["sif_wavelength_sigma"][:, :]),
        )
        required_attributes = Dict(
            key => String(_require_attr(dataset.attrib, key, path))
            for key in ("surface_order", "parameter_names", "parameter_units",
                        "state_order", "aerosol_coordinate_transform"))
        return (; active, xa, Sa, Sa_active, required_arrays,
                required_attributes)
    end
end

function _validate_against_current_prior(target, reference,
                                         target_path, reference_path)
    target.active == reference.active || error(
        "tapered and current priors have different active-state mappings")
    target.xa == reference.xa || error(
        "tapered prior xa differs from the approved current campaign prior")
    isequal(target.required_arrays, reference.required_arrays) || error(
        "tapered prior changes an approved non-CO2/layout array")
    target.required_attributes == reference.required_attributes || error(
        "tapered prior changes state layout metadata")

    full_nonco2 = trues(34, 34)
    full_nonco2[FULL_CO2_PARAMETER_INDEX, FULL_CO2_PARAMETER_INDEX] .= false
    active_co2 = findall(index -> index in FULL_CO2_PARAMETER_INDEX,
                         target.active)
    length(active_co2) == 12 || error(
        "expected exactly 12 active CO2 layers in the retrieval prior")
    active_nonco2 = trues(30, 30)
    active_nonco2[active_co2, active_co2] .= false
    taper = Matrix{Float64}(I, length(active_co2), length(active_co2))
    for row in 1:length(active_co2)-1
        retained = 1.0
        for column in row+1:length(active_co2)
            retained *= EXPECTED_TAPER_RETENTION[column - 1]
            taper[row, column] = taper[column, row] = retained
        end
    end
    changed_co2_offdiagonal = false
    for surface in 1:4
        target.Sa[:, :, surface][full_nonco2] ==
            reference.Sa[:, :, surface][full_nonco2] || error(
            "tapered prior changes Sa outside the CO2 block for surface $surface")
        target.Sa_active[:, :, surface][active_nonco2] ==
            reference.Sa_active[:, :, surface][active_nonco2] || error(
            "tapered prior changes active Sa outside the CO2 block for surface $surface")
        diag(target.Sa[FULL_CO2_PARAMETER_INDEX,
                       FULL_CO2_PARAMETER_INDEX, surface]) ==
            diag(reference.Sa[FULL_CO2_PARAMETER_INDEX,
                              FULL_CO2_PARAMETER_INDEX, surface]) || error(
            "tapered prior changes a CO2 marginal variance for surface $surface")
        target_co2 = target.Sa[FULL_CO2_PARAMETER_INDEX,
                               FULL_CO2_PARAMETER_INDEX, surface]
        reference_co2 = reference.Sa[FULL_CO2_PARAMETER_INDEX,
                                     FULL_CO2_PARAMETER_INDEX, surface]
        target_active_co2 = target.Sa_active[
            active_co2, active_co2, surface]
        reference_active_co2 = reference.Sa_active[
            active_co2, active_co2, surface]
        expected_active_co2 = reference_active_co2 .* taper
        # The builder explicitly restores the reference diagonal after its
        # D*R*D reconstruction, so marginal variances remain bit-identical.
        for index in axes(expected_active_co2, 1)
            expected_active_co2[index, index] =
                reference_active_co2[index, index]
        end
        isapprox(target_active_co2, expected_active_co2;
                 atol=0.0, rtol=32eps(Float64)) || error(
            "tapered prior CO2 block does not equal the complete approved " *
            "12-layer covariance times the nonstationary AR(1) taper for " *
            "surface $surface")

        fixed_co2 = setdiff(collect(FULL_CO2_PARAMETER_INDEX), target.active)
        target.Sa[fixed_co2, FULL_CO2_PARAMETER_INDEX, surface] ==
            reference.Sa[fixed_co2, FULL_CO2_PARAMETER_INDEX, surface] || error(
            "tapered prior changes fixed-upper CO2 covariance for surface $surface")
        target.Sa[FULL_CO2_PARAMETER_INDEX, fixed_co2, surface] ==
            reference.Sa[FULL_CO2_PARAMETER_INDEX, fixed_co2, surface] || error(
            "tapered prior changes fixed-upper CO2 covariance for surface $surface")
        changed_co2_offdiagonal |= target_co2 != reference_co2
    end
    changed_co2_offdiagonal || error(
        "tapered prior is byte-identical to the current CO2 covariance; " *
        "no correlation taper was applied")
    return nothing
end

"""
Require the generated prior's filename/model/hash and unchanged non-CO2 state.

The approved current bottom-layer prior is itself hash-pinned.  The new prior
must have identical `xa`, layout, SIF arrays, and every covariance element
outside the CO2 block.  Inside that block only off-diagonal entries may change:
all CO2 marginal variances must remain exactly equal to the current prior.
"""
function validate_prior_identity(path::AbstractString,
                                 expected_sha256::AbstractString;
                                 reference_path::AbstractString,
                                 reference_sha256::AbstractString)
    _require_file(path, "tapered-CO2 prior")
    _require_file(reference_path, "approved current bottom-layer prior")
    occursin(r"^[0-9a-f]{64}$", expected_sha256) || error(
        "TAPERED_PRIOR_SHA256 must be a lowercase 64-character SHA-256")
    occursin(r"^[0-9a-f]{64}$", reference_sha256) || error(
        "CURRENT_PRIOR_SHA256 must be a lowercase 64-character SHA-256")
    expected_name = "apriori_states_$(EXPECTED_PRIOR_MODEL).nc"
    basename(path) == expected_name || error(
        "the tapered prior must use the explicit filename $expected_name; " *
        "received $(basename(path))")
    actual_sha256 = file_sha256(path)
    actual_sha256 == expected_sha256 || error(
        "tapered prior SHA-256 mismatch: expected $expected_sha256, " *
        "received $actual_sha256 for $path")
    actual_reference_sha256 = file_sha256(reference_path)
    actual_reference_sha256 == reference_sha256 || error(
        "approved current prior SHA-256 mismatch: expected $reference_sha256, " *
        "received $actual_reference_sha256 for $reference_path")

    NCDataset(path) do dataset
        Int(_require_attr(dataset.attrib, "apriori_complete", path)) == 1 ||
            error("prior is not marked apriori_complete=1: $path")
        model = String(_require_attr(
            dataset.attrib, "co2_covariance_model", path))
        model == EXPECTED_PRIOR_MODEL || error(
            "prior uses co2_covariance_model=$model; expected " *
            EXPECTED_PRIOR_MODEL)
        base_model = String(_require_attr(
            dataset.attrib, "co2_covariance_base_model", path))
        base_model == EXPECTED_PRIOR_BASE_MODEL || error(
            "tapered prior uses base model $base_model; expected " *
            EXPECTED_PRIOR_BASE_MODEL)
        String(_require_attr(
            dataset.attrib, "co2_correlation_taper_type", path)) ==
            "nonstationary_ar1_schur" || error(
            "tapered prior uses the wrong correlation-taper construction")
        retention = Float64.(
            dataset["co2_correlation_taper_adjacent_retention"][:])
        retention == EXPECTED_TAPER_RETENTION || error(
            "tapered prior uses unexpected adjacent retention factors")
        base = Float64.(dataset["co2_base_adjacent_correlation"][:])
        selected = Float64.(
            dataset["co2_selected_adjacent_correlation"][:])
        isapprox(selected, base .* retention; atol=2e-15, rtol=2e-14) || error(
            "selected adjacent CO2 correlations do not equal base .* taper")
    end
    target = _prior_snapshot(path)
    reference = _prior_snapshot(reference_path)
    _validate_against_current_prior(
        target, reference, path, reference_path)
    return (; model=EXPECTED_PRIOR_MODEL, sha256=actual_sha256,
            path=resolved_target(path),
            reference_path=resolved_target(reference_path),
            reference_sha256=actual_reference_sha256)
end

"""
Validate the completed direct-SIF publication and return its immutable hashes.

The canonical completion receipt is linked cryptographically to the detailed
validation receipt.  Every published SIF-on scene is then compared with the
validated staged hash.  This prevents a downstream bottom-layer product from
being blessed by a stale receipt after its full-column source was replaced.
"""
function validate_published_sif_release(paths::ReadinessPaths)
    publication_marker = joinpath(
        paths.full_truth_root, ".sif_v2_publication_in_progress")
    !ispath(publication_marker) || error(
        "corrected SIF truth publication is incomplete: $publication_marker")
    complete_receipt = _require_file(
        joinpath(paths.full_truth_root, "sif_v2_release_complete.dat"),
        "corrected SIF publication receipt")
    validation_receipt = _require_file(
        joinpath(paths.restart_root, "sif_v2_release_validation.dat"),
        "corrected SIF validation receipt")

    complete = _metadata(complete_receipt)
    validation = _metadata(validation_receipt)
    validation_sha256 = file_sha256(validation_receipt)
    _required_metadata(
        complete, "validation_receipt_sha256", complete_receipt) ==
        validation_sha256 || error(
        "publication receipt does not match the current validation receipt")
    release_id = _required_metadata(complete, "release_id", complete_receipt)
    release_id == validation_sha256[1:16] || error(
        "publication release_id is inconsistent with its validation receipt")
    parse(Int, _required_metadata(
        validation, "release_schema", validation_receipt)) == 1 || error(
        "unsupported corrected-SIF release schema")
    parse(Int, _required_metadata(
        validation, "sif_definition_version", validation_receipt)) == 2 ||
        error("corrected-SIF release does not use definition version 2")
    _required_metadata(validation, "sif_case", validation_receipt) ==
        "angular_integral760_0p5" || error(
        "corrected-SIF release uses the wrong SIF case")

    full_table = _require_file(
        joinpath(paths.full_truth_root, "true_states.dat"),
        "published full-column truth table")
    table_sha256 = file_sha256(full_table)
    table_sha256 == _required_metadata(
        validation, "corrected_state_table_sha256", validation_receipt) ||
        error("published full-column truth table differs from the validated table")
    full_cases = read_truth_cases(full_table)
    sif_cases = select_sif_truth_cases(full_cases, :on)
    length(sif_cases) == 32 || error(
        "published full-column table does not contain 32 SIF-on states")

    rows = _release_rows(validation_receipt)
    expected_indices = sort!(getfield.(sif_cases, :state_index))
    sort!(collect(keys(rows))) == expected_indices || error(
        "corrected-SIF validation receipt does not enumerate the 32 published states")
    source_hashes = Dict{Int,String}()
    by_index = Dict(case.state_index => case for case in sif_cases)
    for index in expected_indices
        case = by_index[index]
        row = rows[index]
        expected_class = case.aerosol_case == :none ? "clear" : "aerosol"
        row.class == expected_class || error(
            "release class mismatch for full-column state $index")
        relative_destination = joinpath(
            case.aerosol_case == :none ? "" : "aerosol_chunked",
            @sprintf("hiressim_%03d.nc", index))
        expected_destination = joinpath(
            paths.full_truth_root, relative_destination)
        # The validation receipt is deliberately transferable from the
        # publication host to Gattaca. Its absolute root therefore may differ,
        # but the canonical truth-map-relative destination may not.
        receipt_tail = joinpath("truth_map", relative_destination)
        _has_path_tail(row.destination, receipt_tail) || error(
            "release destination mismatch for full-column state $index")
        _require_file(expected_destination, "published corrected-SIF scene")
        actual = file_sha256(expected_destination)
        actual == row.staged_sha256 || error(
            "published corrected-SIF state $index differs from its validated hash")
        source_hashes[index] = actual
    end
    return (; release_id, validation_receipt_sha256=validation_sha256,
            complete_receipt_sha256=file_sha256(complete_receipt),
            table_sha256, source_hashes)
end

function _truth_scene(paths::ReadinessPaths, truth)
    root = truth.aerosol_case == :none ?
        joinpath(paths.bottom_campaign_root, "truth") :
        joinpath(paths.bottom_campaign_root, "truth", "aerosol_chunked")
    return joinpath(root, @sprintf("hiressim_%03d.nc", truth.state_index))
end

function _no_sif_triplet_set_sha256(paths::ReadinessPaths)
    cases = sort!(select_sif_truth_cases(
        read_truth_cases(paths.bottom_truth_table), :off);
        by=truth -> truth.state_index)
    length(cases) == 40 || error(
        "bottom-layer table does not contain exactly 40 no-SIF states")
    lines = String[]
    for truth in cases
        index = truth.state_index
        products = (
            truth=_require_file(
                _truth_scene(paths, truth), "released no-SIF truth"),
            measurement=_require_file(joinpath(
                paths.measurement_directory,
                @sprintf("OCO2sims_%03d.nc", index)),
                "released no-SIF measurement"),
            noise=_require_file(joinpath(
                paths.noise_directory,
                @sprintf("OCO2noise_%03d.nc", index)),
                "released no-SIF noise covariance"),
        )
        for kind in (:truth, :measurement, :noise)
            push!(lines, "$(String(kind)) $index " *
                         file_sha256(getproperty(products, kind)))
        end
    end
    return bytes2hex(sha256(codeunits(join(lines, '\n') * "\n")))
end

"""
Require the dedicated bottom-layer SIF-v2 release receipt and verify all files.

This is a second release boundary after the full-column publication: it covers
assembly of bottom-layer truth, synthetic OCO processing, and frozen-noise
generation.  Retrievals are forbidden while its interruption marker exists or
when even one current file differs from the receipt.
"""
function validate_published_bottom_sif_release(paths::ReadinessPaths,
                                               full_release)
    truth_root = joinpath(paths.bottom_campaign_root, "truth")
    marker = joinpath(truth_root, ".bottom_layer_sif_v2_publication_in_progress")
    !ispath(marker) || error(
        "bottom-layer SIF-v2 publication is incomplete: $marker")
    receipt = _require_file(
        joinpath(truth_root, "bottom_layer_sif_v2_release_complete.dat"),
        "bottom-layer SIF-v2 release receipt")
    metadata = _metadata(receipt)
    parse(Int, _required_metadata(metadata, "release_schema", receipt)) == 1 ||
        error("unsupported bottom-layer SIF release schema")
    parse(Int, _required_metadata(
        metadata, "sif_definition_version", receipt)) == 2 || error(
        "bottom-layer SIF release does not use definition version 2")
    _required_metadata(
        metadata, "full_column_release_receipt_sha256", receipt) ==
        full_release.complete_receipt_sha256 || error(
        "bottom-layer release is not bound to the current full-column release")
    _required_metadata(
        metadata, "full_column_state_table_sha256", receipt) ==
        full_release.table_sha256 || error(
        "bottom-layer release is not bound to the current full-column table")

    bottom_table_sha256 = file_sha256(paths.bottom_truth_table)
    _required_metadata(metadata, "bottom_state_table_sha256", receipt) ==
        bottom_table_sha256 || error(
        "bottom-layer release is not bound to the current bottom state table")
    input_set_sha256 = _required_metadata(metadata, "input_set_sha256", receipt)
    occursin(r"^[0-9a-f]{64}$", input_set_sha256) || error(
        "bottom-layer release has an invalid input_set_sha256")

    archive_root = resolved_target(get(
        ENV, "BOTTOM_SIF_V1_ARCHIVE_ROOT",
        joinpath(paths.private_root, "archive",
                 "sif_wavelength_integral_0p5_20260904",
                 "bottom_layer_XCO2_retrievals")))
    contains_path(paths.private_root, archive_root) || error(
        "bottom-layer legacy archive must remain under RRS_PRIVATE_ROOT")
    contains_path(paths.repo_root, archive_root) && error(
        "bottom-layer legacy archive must not be stored in Git")
    legacy_table_relative = _required_metadata(
        metadata, "legacy_bottom_state_table_archive_relative", receipt)
    legacy_table_relative == joinpath("truth", "true_states.dat") || error(
        "bottom-layer release names an unexpected legacy state-table path")
    archive_manifest_relative = _required_metadata(
        metadata, "legacy_archive_manifest_relative", receipt)
    archive_manifest_relative == "bottom_layer_legacy_manifest.dat" || error(
        "bottom-layer release names an unexpected legacy archive manifest")
    legacy_table = _require_file(joinpath(
        archive_root, legacy_table_relative), "archived legacy state table")
    archive_manifest = _require_file(joinpath(
        archive_root, archive_manifest_relative), "legacy archive manifest")
    file_sha256(legacy_table) == _required_metadata(
        metadata, "legacy_bottom_state_table_sha256", receipt) || error(
        "archived legacy state table differs from the bottom release receipt")
    file_sha256(archive_manifest) == _required_metadata(
        metadata, "legacy_archive_manifest_sha256", receipt) || error(
        "legacy archive manifest differs from the bottom release receipt")
    _required_metadata(
        metadata, "no_sif_byte_preservation_policy", receipt) ==
        "preserve_truth_measurement_noise_bytes_and_legacy_state_table_hash_attributes" ||
        error("bottom release does not guarantee no-SIF byte preservation")
    no_sif_triplet_set_sha256 = _no_sif_triplet_set_sha256(paths)
    _required_metadata(metadata, "no_sif_triplet_set_sha256", receipt) ==
        no_sif_triplet_set_sha256 || error(
        "current no-SIF inputs differ from the preserved release set")

    rows = _bottom_release_rows(receipt)
    sort!(collect(keys(rows))) == EXPECTED_SIF_STATES || error(
        "bottom-layer SIF release does not enumerate the canonical 40 states")
    for truth in select_sif_truth_cases(
            read_truth_cases(paths.bottom_truth_table), :on)
        index = truth.state_index
        row = rows[index]
        truth_path = _require_file(
            _truth_scene(paths, truth), "released bottom-layer SIF truth")
        measurement_path = _require_file(joinpath(
            paths.measurement_directory, @sprintf("OCO2sims_%03d.nc", index)),
            "released bottom-layer SIF measurement")
        noise_path = _require_file(joinpath(
            paths.noise_directory, @sprintf("OCO2noise_%03d.nc", index)),
            "released bottom-layer SIF noise covariance")
        for (path, expected, description) in (
                (truth_path, row.truth_sha256, "truth"),
                (measurement_path, row.measurement_sha256, "measurement"),
                (noise_path, row.noise_sha256, "noise"))
            file_sha256(path) == expected || error(
                "bottom-layer SIF state $index $description differs from " *
                "its release receipt")
        end
    end
    return (; receipt_sha256=file_sha256(receipt), input_set_sha256,
            bottom_table_sha256, no_sif_triplet_set_sha256,
            archive_manifest_sha256=file_sha256(archive_manifest),
            rows, path=receipt)
end

function _read_attributes(path::AbstractString)
    NCDataset(path) do dataset
        Dict{String,Any}(String(key) => dataset.attrib[key]
                         for key in keys(dataset.attrib))
    end
end

function _require_bottom_metadata(attributes, truth, path, table_sha256,
                                  full_table_sha256, released_source_hashes)
    Int(_require_attr(attributes, "state_index", path)) == truth.state_index ||
        error("state metadata mismatch in $path")
    String(_require_attr(attributes, "campaign", path)) ==
        "bottom_layer_XCO2" || error("wrong campaign in $path")
    Int(_require_attr(attributes, "bottom_layer_truth_complete", path)) == 1 ||
        error("bottom-layer truth is not marked complete in $path")
    String(_require_attr(attributes, "sif_case", path)) ==
        String(truth.sif_case) || error("SIF case mismatch in $path")
    String(_require_attr(attributes, "state_table_sha256", path)) ==
        table_sha256 || error("bottom-layer state-table hash mismatch in $path")
    String(_require_attr(
        attributes, "full_column_state_table_sha256", path)) ==
        full_table_sha256 || error(
        "full-column state-table hash mismatch in $path")
    source_index = Int(_require_attr(
        attributes, "full_column_source_state", path))
    haskey(released_source_hashes, source_index) || error(
        "$path traces to full-column state $source_index, which is absent " *
        "from the corrected-SIF release")
    String(_require_attr(attributes, "full_column_source_sha256", path)) ==
        released_source_hashes[source_index] || error(
        "$path does not trace to the published corrected-SIF source hash")
    validated_sif_provenance(attributes, truth.sif_case; source=path)
    return source_index
end

"""Validate all 40 bottom-layer version-2 SIF measurement/noise pairs."""
function validate_bottom_sif_inputs(paths::ReadinessPaths, release,
                                    bottom_release)
    _require_file(paths.bottom_truth_table, "bottom-layer truth table")
    cases = read_truth_cases(paths.bottom_truth_table)
    selected = select_sif_truth_cases(cases, :on)
    indices = sort!(getfield.(selected, :state_index))
    indices == EXPECTED_SIF_STATES || error(
        "bottom-layer truth table does not contain the canonical 40 SIF-on states")
    table_sha256 = file_sha256(paths.bottom_truth_table)
    input_records = String[
        "bottom_truth_table $table_sha256 $(paths.bottom_truth_table)",
        "full_truth_table $(release.table_sha256) " *
        joinpath(paths.full_truth_root, "true_states.dat"),
    ]
    instrument_processor = _require_file(joinpath(
        paths.repo_root, "RRS_XCO2", "inversion", "instrument",
        "process_truth_map.jl"), "instrument processor source")
    noise_processor = _require_file(joinpath(
        paths.repo_root, "RRS_XCO2", "inversion", "instrument",
        "generate_noise_covariances.jl"), "noise processor source")
    instrument_processor_sha256 = file_sha256(instrument_processor)
    noise_processor_sha256 = file_sha256(noise_processor)
    append!(input_records, (
        "instrument_processor $instrument_processor_sha256 $instrument_processor",
        "noise_processor $noise_processor_sha256 $noise_processor",
    ))

    for truth in selected
        index = truth.state_index
        truth_path = _require_file(
            _truth_scene(paths, truth), "bottom-layer SIF truth scene")
        measurement_path = _require_file(joinpath(
            paths.measurement_directory, @sprintf("OCO2sims_%03d.nc", index)),
            "bottom-layer SIF OCO measurement")
        noise_path = _require_file(joinpath(
            paths.noise_directory, @sprintf("OCO2noise_%03d.nc", index)),
            "bottom-layer SIF noise covariance")

        truth_attributes = _read_attributes(truth_path)
        measurement_attributes = _read_attributes(measurement_path)
        noise_attributes = _read_attributes(noise_path)
        _require_bottom_metadata(
            truth_attributes, truth, truth_path, table_sha256,
            release.table_sha256, release.source_hashes)
        _require_bottom_metadata(
            measurement_attributes, truth, measurement_path, table_sha256,
            release.table_sha256, release.source_hashes)
        _require_bottom_metadata(
            noise_attributes, truth, noise_path, table_sha256,
            release.table_sha256, release.source_hashes)
        Int(_require_attr(
            measurement_attributes, "instrument_processing_complete",
            measurement_path)) == 1 || error(
            "instrument processing is incomplete in $measurement_path")
        Int(_require_attr(
            noise_attributes, "noise_covariance_complete", noise_path)) == 1 ||
            error("noise covariance is incomplete in $noise_path")

        source_truth = String(_require_attr(
            measurement_attributes, "source_truth_scene", measurement_path))
        source_truth_tail = joinpath(
            "bottom_layer_XCO2_retrievals", "truth",
            truth.aerosol_case == :none ? "" : "aerosol_chunked",
            basename(truth_path))
        _has_path_tail(source_truth, source_truth_tail) || error(
            "measurement $measurement_path does not name its canonical " *
            "bottom-layer source scene")
        String(_require_attr(
            measurement_attributes, "source_truth_sha256",
            measurement_path)) == bottom_release.rows[index].truth_sha256 ||
            error("measurement $measurement_path is not bound to the " *
                  "released truth bytes")
        String(_require_attr(
            measurement_attributes, "instrument_processor_script_sha256",
            measurement_path)) == instrument_processor_sha256 || error(
            "measurement $measurement_path was made by a different " *
            "instrument processor")
        _has_path_tail(String(_require_attr(
            measurement_attributes, "instrument_processor_script",
            measurement_path)),
            joinpath("RRS_XCO2", "inversion", "instrument",
                     "process_truth_map.jl")) || error(
            "measurement $measurement_path names the wrong instrument processor")

        measurement_source = String(_require_attr(
            noise_attributes, "source_synthetic_measurement", noise_path))
        _has_path_tail(measurement_source, joinpath(
            "bottom_layer_XCO2_retrievals", "truth", "OCO_radiances",
            basename(measurement_path))) || error(
            "noise $noise_path does not name its canonical measurement")
        String(_require_attr(
            noise_attributes, "source_synthetic_measurement_sha256",
            noise_path)) == bottom_release.rows[index].measurement_sha256 ||
            error("noise $noise_path is not bound to the released " *
                  "measurement bytes")
        String(_require_attr(
            noise_attributes, "noise_processor_script_sha256",
            noise_path)) == noise_processor_sha256 || error(
            "noise $noise_path was made by a different covariance processor")
        _has_path_tail(String(_require_attr(
            noise_attributes, "noise_processor_script", noise_path)),
            joinpath("RRS_XCO2", "inversion", "instrument",
                     "generate_noise_covariances.jl")) || error(
            "noise $noise_path names the wrong covariance processor")
        for key in BOTTOM_PROVENANCE_ATTRIBUTES
            haskey(truth_attributes, key) || error(
                "$truth_path lacks bottom-layer provenance attribute $key")
            for (attributes, path) in (
                    (measurement_attributes, measurement_path),
                    (noise_attributes, noise_path))
                _require_attr(attributes, key, path) == truth_attributes[key] ||
                    error("$key differs between $path and $truth_path")
            end
        end

        # Exercise the same loader that the production retrieval uses.  This
        # checks both class vectors, frozen Se, spectral partition, and exact
        # version-2 SIF provenance rather than only examining file attributes.
        experiments = build_experiments(
            [truth]; measurement_directory=paths.measurement_directory,
            noise_directory=paths.noise_directory)
        for experiment in experiments
            experiment.noise_index == UNPERTURBED_INDEX || continue
            realization = load_measurement_realization(experiment)
            length(realization.noiseless) == EXPECTED_MEASUREMENT_COUNT || error(
                "state $index has $(length(realization.noiseless)) samples; " *
                "expected $EXPECTED_MEASUREMENT_COUNT")
        end

        push!(input_records,
              @sprintf("state %03d truth %s measurement %s noise %s",
                       index, bottom_release.rows[index].truth_sha256,
                       bottom_release.rows[index].measurement_sha256,
                       bottom_release.rows[index].noise_sha256))
    end
    digest = bytes2hex(sha256(codeunits(join(input_records, '\n'))))
    return (; states=indices, table_sha256, input_set_sha256=digest,
            records=input_records)
end

function validate_readiness(paths::ReadinessPaths;
                            prior_sha256::AbstractString,
                            reference_prior_path::AbstractString,
                            reference_prior_sha256::AbstractString)
    validate_output_isolation(paths)
    validate_prior_isolation(paths)
    prior = validate_prior_identity(
        paths.prior_path, prior_sha256;
        reference_path=reference_prior_path,
        reference_sha256=reference_prior_sha256)
    release = validate_published_sif_release(paths)
    bottom_release = validate_published_bottom_sif_release(paths, release)
    inputs = validate_bottom_sif_inputs(paths, release, bottom_release)
    return (; prior, release, bottom_release, inputs)
end

function _identity_text(readiness, source_sha256)
    fields = (
        "identity_schema" => string(IDENTITY_SCHEMA),
        "campaign_id" => CAMPAIGN_ID,
        "source_checkpoint_sha" => source_sha256,
        "prior_model" => readiness.prior.model,
        "prior_sha256" => readiness.prior.sha256,
        "prior_path" => readiness.prior.path,
        "reference_prior_sha256" => readiness.prior.reference_sha256,
        "reference_prior_path" => readiness.prior.reference_path,
        "sif_release_id" => readiness.release.release_id,
        "sif_release_validation_sha256" =>
            readiness.release.validation_receipt_sha256,
        "bottom_sif_release_receipt_sha256" =>
            readiness.bottom_release.receipt_sha256,
        "bottom_sif_release_input_set_sha256" =>
            readiness.bottom_release.input_set_sha256,
        "bottom_sif_no_sif_triplet_set_sha256" =>
            readiness.bottom_release.no_sif_triplet_set_sha256,
        "bottom_sif_legacy_archive_manifest_sha256" =>
            readiness.bottom_release.archive_manifest_sha256,
        "bottom_sif_input_set_sha256" => readiness.inputs.input_set_sha256,
        "bottom_truth_table_sha256" => readiness.inputs.table_sha256,
        "sif_filter" => "on",
        "measurement_classes" => "corrected,uncorrected",
    )
    return join(("$key=$value" for (key, value) in fields), '\n') * "\n"
end

function _existing_products(output_root)
    isdir(output_root) || return String[]
    products = String[]
    for class in ("corrected", "uncorrected")
        directory = joinpath(output_root, class)
        isdir(directory) || continue
        append!(products, filter(
            path -> occursin(r"^retrieval_state\d{3}_perturbation\d{2}\.nc$",
                             basename(path)),
            joinpath.(directory, readdir(directory))))
    end
    return products
end

"""
Initialize or verify the immutable campaign identity beside the isolated output.

Concurrent Slurm tasks serialize initialization with a directory lock.  A
missing identity alongside existing retrieval products is always an error; no
task may adopt unlabelled output from an earlier campaign.
"""
function initialize_campaign_identity!(paths::ReadinessPaths, readiness,
                                       source_sha256::AbstractString;
                                       wait_seconds::Real=30)
    occursin(r"^[0-9a-f]{40}$", source_sha256) || error(
        "source checkpoint must be a lowercase 40-character Git SHA")
    validate_output_isolation(paths)
    expected = _identity_text(readiness, source_sha256)
    control = joinpath(paths.output_root, ".control")
    identity = joinpath(control, "campaign_identity.dat")
    lock = joinpath(control, ".identity_lock")
    mkpath(control)

    if isfile(identity)
        read(identity, String) == expected || error(
            "existing campaign identity does not match current source/prior/inputs: " *
            identity)
        return identity
    end
    isempty(_existing_products(paths.output_root)) || error(
        "isolated output contains retrievals but has no campaign identity")

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
            "concurrent campaign identity differs from current source/prior/inputs")
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

function main()
    paths = ReadinessPaths()
    prior_sha256 = get(ENV, "TAPERED_PRIOR_SHA256", "")
    reference_prior_sha256 = get(ENV, "CURRENT_PRIOR_SHA256", "")
    reference_prior_path = get(
        ENV, "CURRENT_PRIOR_PATH",
        joinpath(paths.bottom_campaign_root, "retrieval_setup",
                 "apriori_states.nc"))
    source_sha256 = get(ENV, "CAMPAIGN_CHECKPOINT_SHA", "")
    readiness = validate_readiness(
        paths; prior_sha256, reference_prior_path,
        reference_prior_sha256)
    identity = initialize_campaign_identity!(
        paths, readiness, source_sha256)
    println("Gattaca tapered-SIF retrieval readiness: PASSED")
    println("campaign_id=$(CAMPAIGN_ID)")
    println("prior_model=$(readiness.prior.model)")
    println("prior_sha256=$(readiness.prior.sha256)")
    println("reference_prior_sha256=$(readiness.prior.reference_sha256)")
    println("sif_release_id=$(readiness.release.release_id)")
    println("bottom_sif_release_receipt_sha256=" *
            readiness.bottom_release.receipt_sha256)
    println("bottom_sif_no_sif_triplet_set_sha256=" *
            readiness.bottom_release.no_sif_triplet_set_sha256)
    println("bottom_sif_input_set_sha256=$(readiness.inputs.input_set_sha256)")
    println("states=$(join(readiness.inputs.states, ','))")
    println("output_root=$(paths.output_root)")
    println("campaign_identity=$identity")
end

end # module GattacaTaperedSIFReadiness

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    GattacaTaperedSIFReadiness.main()
end
