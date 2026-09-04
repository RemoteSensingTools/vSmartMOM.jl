#!/usr/bin/env julia

"""
Stage, validate, and transactionally publish the corrected bottom-layer SIF
products without recomputing radiative transfer.

The default action is `validate`; it never modifies the canonical campaign.
The complete sequence is:

1. `BOTTOM_SIF_V2_ACTION=stage` assembles 40 corrected SIF-on bottom-layer
   truth scenes from the *published* full-column SIF-v2 A band plus the
   existing bottom-layer no-SIF CO2 bands.  It then invokes the existing
   synthetic-OCO and noise producers inside an isolated staging tree.
2. the default `validate` action verifies all 40 truth/measurement/noise
   triplets, exact O2 and CO2 component reuse, strict SIF-v2 provenance, and
   byte preservation of all 40 no-SIF products;
3. `BOTTOM_SIF_V2_ACTION=publish` additionally requires
   `CONFIRM_BOTTOM_SIF_V2_PUBLICATION=publish_40_bottom_sif_scenes`.  It first
   archives the legacy table and every active legacy SIF truth, measurement,
   noise, and retrieval product.  Verified same-directory sidecars are then
   promoted, with the authoritative state table last; and
4. a canonical release receipt binds the full-column release, corrected
   bottom-layer table, inputs, and all 40 published product triplets.

No-SIF truth, synthetic radiance, noise, or retrieval files are ever replaced.
Their SHA-256 hashes are captured before staging/publication and rechecked
afterward.  The publisher refuses active retrieval claims and requires the
external-SIF ownership marker, so no local SIF writer can race the migration.

Because POSIX cannot atomically replace 124 independent files, an uncatchable
process death can interrupt the rename sequence.  The persistent
`.bottom_layer_sif_v2_publication_in_progress` marker is therefore a mandatory
downstream barrier.  `BOTTOM_SIF_V2_ACTION=restore` plus
`CONFIRM_BOTTOM_SIF_V2_RESTORE=restore_archived_bottom_sif_v1` restores the
fully checksummed legacy archive after such an interruption.  Never remove the
marker by hand.

Environment paths:

* `BOTTOM_XCO2_CAMPAIGN_ROOT` (default
  `RRS_XCO2/bottom_layer_XCO2_retrievals`)
* `FULL_COLUMN_TRUTH_ROOT` (default `RRS_XCO2/truth_map`)
* `BOTTOM_SIF_V2_RESTART_ROOT` (default `<campaign>/.sif_v2_restart`)
* `BOTTOM_SIF_V1_ARCHIVE_ROOT` (default
  `RRS_XCO2/obsolete/sif_wavelength_integral_0p5_20260904/`
  `bottom_layer_XCO2_retrievals`)
* `FULL_COLUMN_SIF_VALIDATION_RECEIPT` (default
  `<full truth>/.sif_v2_restart/sif_v2_release_validation.dat`)
* `OCO_STOKES_COEFFICIENTS` and `OCO_SNR_COEFFICIENTS` select the exact
  instrument coefficient files and are included in the release checksum.
* `BREAK_STALE_BOTTOM_SIF_V2_LOCK=1` is permitted only when no publication
  marker exists.  Set it after verifying that a publisher killed before
  marker creation is no longer active; otherwise the lock remains fail-closed.
"""

module BottomLayerSIFV2Migration

using Dates
using DelimitedFiles
using NCDatasets
using Printf
using SHA
using UUIDs

include(joinpath(@__DIR__, "common.jl"))
using .RRSXCO2Common
include(joinpath(@__DIR__, "bottom_layer_truth_common.jl"))
using .BottomLayerTruthCommon
include(joinpath(@__DIR__, "..", "inversion", "RetrievalCases.jl"))
using .RetrievalCases: validated_sif_provenance, matching_sif_provenance

export MigrationPaths, stage_release, validate_release, publish_release,
       restore_legacy, main

const CONFIRM_PUBLICATION = "publish_40_bottom_sif_scenes"
const CONFIRM_RESTORE = "restore_archived_bottom_sif_v1"
const RELEASE_SCHEMA = 1
const FULL_RELEASE_RECEIPT_NAME = "sif_v2_release_complete.dat"
const FULL_RELEASE_MARKER_NAME = ".sif_v2_publication_in_progress"
const BOTTOM_RELEASE_RECEIPT_NAME = "bottom_layer_sif_v2_release_complete.dat"
const BOTTOM_RELEASE_MARKER_NAME =
    ".bottom_layer_sif_v2_publication_in_progress"
const EXTERNAL_OWNERSHIP_RELATIVE = joinpath(
    "retrievals", ".control", "sif_owned_externally")
const SIF_TRUTH_VARIABLES = (
    "radiance_rayleigh_o2a",
    "radiance_cabannes_o2a",
    "radiance_rrs_o2a",
)
const CO2_TRUTH_VARIABLES = (
    "radiance_rayleigh_weak_co2",
    "radiance_rayleigh_strong_co2",
)
const CO2_SOURCE_ATTRIBUTES = (
    "co2_absco_completed",
    "co2_absco_regeneration_complete",
    "weak_h2o_absco_lut",
    "weak_co2_absco_lut",
    "strong_h2o_absco_lut",
    "strong_co2_absco_lut",
    "h2o_line_absorption_by_band",
    "strong_co2_convolution_support_sigma",
    "strong_co2_short_shoulder_merged",
    "strong_co2_short_shoulder_points",
    "strong_co2_shoulder_merged_utc",
    "strong_co2_shoulder_source",
)
const ALL_TRUTH_VARIABLES = (SIF_TRUTH_VARIABLES..., CO2_TRUTH_VARIABLES...)
const STATIC_TABLES = ("true_states.dat", "control_reuse_map.dat",
                       "scene_components.dat")
const SIF_STATE_RANGES = (
    6:10, 16:20, 26:30, 36:40,
    46:50, 56:60, 66:70, 76:80,
)
const EXPECTED_SIF_STATES = sort!(vcat(collect.(SIF_STATE_RANGES)...))
const EXPECTED_NO_SIF_STATES = sort!(setdiff(collect(1:80),
                                             EXPECTED_SIF_STATES))
const EXPECTED_FULL_COLUMN_SIF_STATES = sort!(vcat(
    [collect(offset .+ (5:8)) for offset in (0, 8, 16, 24, 32, 40, 48, 56)]...))
const LEGACY_SIF_CASE = "total_0p5"
const INSTRUMENT_ROOT = normpath(joinpath(
    @__DIR__, "..", "inversion", "instrument"))
const INSTRUMENT_PROCESSOR = joinpath(INSTRUMENT_ROOT, "process_truth_map.jl")
const NOISE_PROCESSOR = joinpath(
    INSTRUMENT_ROOT, "generate_noise_covariances.jl")
const STOKES_COEFFICIENTS = joinpath(
    INSTRUMENT_ROOT, "representative_stokes_coefficients.nc")
const SNR_COEFFICIENTS = joinpath(
    INSTRUMENT_ROOT, "representative_snr_coefficients.nc")

struct MigrationPaths
    campaign_root::String
    full_truth_root::String
    restart_root::String
    archive_root::String
    full_validation_receipt::String
    stokes_coefficients::String
    snr_coefficients::String
end

"Resolve symlinks in the longest existing prefix of a possibly-new path."
function resolved_target(path)
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

function contains_path(parent, child)
    relative = relpath(child, parent)
    return relative == "." ||
           !(relative == ".." || startswith(
               relative, "..$(Base.Filesystem.path_separator)"))
end

function MigrationPaths(; campaign_root=nothing, full_truth_root=nothing,
                          restart_root=nothing, archive_root=nothing,
                          full_validation_receipt=nothing,
                          stokes_coefficients=nothing,
                          snr_coefficients=nothing)
    rrs_root = normpath(joinpath(@__DIR__, ".."))
    campaign = resolved_target(something(campaign_root, get(
        ENV, "BOTTOM_XCO2_CAMPAIGN_ROOT",
        joinpath(rrs_root, "bottom_layer_XCO2_retrievals"))))
    full = resolved_target(something(full_truth_root, get(
        ENV, "FULL_COLUMN_TRUTH_ROOT", joinpath(rrs_root, "truth_map"))))
    restart = resolved_target(something(restart_root, get(
        ENV, "BOTTOM_SIF_V2_RESTART_ROOT",
        joinpath(campaign, ".sif_v2_restart"))))
    archive = resolved_target(something(archive_root, get(
        ENV, "BOTTOM_SIF_V1_ARCHIVE_ROOT",
        joinpath(rrs_root, "obsolete", "sif_wavelength_integral_0p5_20260904",
                 "bottom_layer_XCO2_retrievals"))))
    validation = resolved_target(something(full_validation_receipt, get(
        ENV, "FULL_COLUMN_SIF_VALIDATION_RECEIPT",
        joinpath(full, ".sif_v2_restart", "sif_v2_release_validation.dat"))))
    stokes = resolved_target(something(stokes_coefficients, get(
        ENV, "OCO_STOKES_COEFFICIENTS", STOKES_COEFFICIENTS)))
    snr = resolved_target(something(snr_coefficients, get(
        ENV, "OCO_SNR_COEFFICIENTS", SNR_COEFFICIENTS)))

    (contains_path(restart, full) || contains_path(full, restart)) && error(
        "bottom restart root and full-column truth must not overlap")
    bottom_truth = resolved_target(joinpath(campaign, "truth"))
    (contains_path(restart, bottom_truth) ||
     contains_path(bottom_truth, restart)) && error(
        "bottom restart root and canonical bottom truth must not overlap")
    (contains_path(archive, campaign) || contains_path(campaign, archive)) &&
        error("archive and canonical bottom campaign must not overlap")
    (contains_path(archive, restart) || contains_path(restart, archive)) &&
        error("archive and bottom restart root must not overlap")
    return MigrationPaths(campaign, full, restart, archive, validation,
                          stokes, snr)
end

truth_root(paths::MigrationPaths) = joinpath(paths.campaign_root, "truth")
oco_root(paths::MigrationPaths) = joinpath(truth_root(paths), "OCO_radiances")
noise_root(paths::MigrationPaths) = joinpath(oco_root(paths), "noise_covariances")
retrieval_root(paths::MigrationPaths) = joinpath(paths.campaign_root, "retrievals")
stage_campaign_root(paths::MigrationPaths) = joinpath(paths.restart_root,
                                                       "campaign")
stage_truth_root(paths::MigrationPaths) = joinpath(stage_campaign_root(paths),
                                                   "truth")
stage_oco_root(paths::MigrationPaths) = joinpath(stage_truth_root(paths),
                                                "OCO_radiances")
stage_noise_root(paths::MigrationPaths) = joinpath(stage_oco_root(paths),
                                                  "noise_covariances")
canonical_receipt(paths::MigrationPaths) =
    joinpath(truth_root(paths), BOTTOM_RELEASE_RECEIPT_NAME)
publication_marker(paths::MigrationPaths) =
    joinpath(truth_root(paths), BOTTOM_RELEASE_MARKER_NAME)
publication_lock(paths::MigrationPaths) =
    joinpath(paths.restart_root, ".bottom_layer_sif_v2_publication_lock")
stage_receipt(paths::MigrationPaths) =
    joinpath(paths.restart_root, "bottom_layer_sif_v2_staging_validation.dat")
archive_manifest(paths::MigrationPaths) =
    joinpath(paths.archive_root, "bottom_layer_legacy_manifest.dat")

file_sha256(path::AbstractString) = open(path, "r") do io
    bytes2hex(sha256(io))
end

function write_atomic(writer::Function, path::AbstractString)
    mkpath(dirname(path))
    temporary = path * ".tmp.$(getpid()).$(uuid4())"
    try
        open(writer, temporary, "w")
        Base.Filesystem.rename(temporary, path)
    finally
        isfile(temporary) && rm(temporary)
    end
    return path
end

function copy_verified(source::AbstractString, destination::AbstractString;
                       expected_hash=file_sha256(source))
    isfile(source) || error("missing source file: $source")
    file_sha256(source) == expected_hash || error(
        "source checksum changed before copy: $source")
    if isfile(destination)
        file_sha256(destination) == expected_hash || error(
            "existing destination differs from its source: $destination")
        return destination
    end
    mkpath(dirname(destination))
    temporary = destination * ".tmp.$(getpid()).$(uuid4())"
    try
        cp(source, temporary)
        file_sha256(temporary) == expected_hash || error(
            "verified copy failed: $source -> $temporary")
        Base.Filesystem.rename(temporary, destination)
    finally
        isfile(temporary) && rm(temporary)
    end
    return destination
end

function scene_root(root::AbstractString, aerosol_index::Integer)
    return aerosol_index == 2 ? joinpath(root, "aerosol_chunked") : root
end

truth_scene(root, state) = joinpath(
    scene_root(root, state.aerosol_index),
    @sprintf("hiressim_%03d.nc", state.index))
oco_scene(root, state_index::Integer) =
    joinpath(root, @sprintf("OCO2sims_%03d.nc", state_index))
noise_scene(root, state_index::Integer) =
    joinpath(root, @sprintf("OCO2noise_%03d.nc", state_index))

function full_source_scene(paths::MigrationPaths, state)
    old = BottomLayerTruthCommon.old_control_index(state)
    root = state.aerosol_index == 2 ?
        joinpath(paths.full_truth_root, "aerosol_chunked") :
        paths.full_truth_root
    return joinpath(root, @sprintf("hiressim_%03d.nc", old))
end

function required_attribute(attributes, key, path)
    haskey(attributes, key) || error("$path lacks required attribute $key")
    return attributes[key]
end

function finite_array(dataset, name, path)
    haskey(dataset, name) || error("$path lacks $name")
    raw = dataset[name][:, :]
    any(ismissing, raw) && error("$path contains missing $name values")
    values = Array(raw)
    size(values, 1) == 3 || error("$path:$name must have three Stokes rows")
    all(isfinite, values) || error("$path:$name contains non-finite values")
    maximum(abs, values) < 1e30 || error("$path:$name contains fill values")
    return values
end

function same_bits(left::AbstractArray, right::AbstractArray)
    size(left) == size(right) || return false
    eltype(left) == eltype(right) || return false
    isbitstype(eltype(left)) || return left == right
    return reinterpret(UInt8, vec(left)) == reinterpret(UInt8, vec(right))
end

function arrays_match(left_path, right_path, variables)
    return NCDataset(left_path) do left
        NCDataset(right_path) do right
            all(variable -> same_bits(
                    finite_array(left, variable, left_path),
                    finite_array(right, variable, right_path)), variables)
        end
    end
end

function read_keyed_receipt(path)
    isfile(path) || error("missing release receipt: $path")
    values = Dict{String,String}()
    rows = String[]
    for line in eachline(path)
        fields = split(strip(line))
        isempty(fields) && continue
        if first(fields) == "#"
            length(fields) == 3 && (values[fields[2]] = fields[3])
        elseif length(fields) == 2
            values[fields[1]] = fields[2]
        else
            push!(rows, line)
        end
    end
    return values, rows
end

function validate_full_column_release(paths::MigrationPaths, states)
    marker = joinpath(paths.full_truth_root, FULL_RELEASE_MARKER_NAME)
    ispath(marker) && error(
        "full-column SIF publication is still in progress: $marker")
    receipt = joinpath(paths.full_truth_root, FULL_RELEASE_RECEIPT_NAME)
    receipt_values, _ = read_keyed_receipt(receipt)
    validation_hash = get(receipt_values, "validation_receipt_sha256", "")
    length(validation_hash) == 64 || error(
        "full-column release receipt lacks its validation-receipt checksum")
    file_sha256(paths.full_validation_receipt) == validation_hash || error(
        "full-column validation receipt does not match its canonical release receipt")

    validation_values, validation_rows = read_keyed_receipt(
        paths.full_validation_receipt)
    get(validation_values, "release_schema", "") == "1" || error(
        "unsupported full-column release schema")
    get(validation_values, "sif_definition_version", "") ==
        string(RRSXCO2Common.SIF_DEFINITION_VERSION) || error(
        "full-column release is not SIF definition version 2")
    full_table = joinpath(paths.full_truth_root, "true_states.dat")
    get(validation_values, "corrected_state_table_sha256", "") ==
        file_sha256(full_table) || error(
        "published full-column state table differs from its release validation")

    released_hashes = Dict{Int,String}()
    for line in validation_rows
        fields = split(strip(line))
        length(fields) == 5 || error(
            "malformed full-column validation row: $line")
        index = parse(Int, fields[1])
        haskey(released_hashes, index) && error(
            "duplicate full-column release row for state $index")
        released_hashes[index] = fields[3]
    end
    sort!(collect(keys(released_hashes))) == EXPECTED_FULL_COLUMN_SIF_STATES ||
        error("full-column validation receipt has the wrong 32-state SIF set")

    source_states = unique(BottomLayerTruthCommon.old_control_index.(states))
    length(source_states) == 8 || error(
        "bottom-layer SIF assembly should use exactly eight full-column controls")
    state_by_old = Dict(BottomLayerTruthCommon.old_control_index(state) => state
                        for state in states)
    source_hashes = Dict{Int,String}()
    for old_index in sort(source_states)
        state = state_by_old[old_index]
        source = full_source_scene(paths, state)
        isfile(source) || error("missing published full-column source: $source")
        digest = file_sha256(source)
        get(released_hashes, old_index, "") == digest || error(
            "full-column source $old_index differs from its release receipt")
        NCDataset(source) do dataset
            Int(required_attribute(dataset.attrib, "simulation_complete", source)) == 1 ||
                error("full-column source is incomplete: $source")
            Int(required_attribute(dataset.attrib, "o2_truth_regenerated", source)) == 1 ||
                error("full-column source is not direct regenerated truth: $source")
            Int(required_attribute(dataset.attrib, "o2_truth_reused", source)) == 0 ||
                error("full-column source still claims reused O2 truth: $source")
            validated_sif_provenance(
                dataset.attrib, Symbol(state.sif_case); source)
            for variable in ALL_TRUTH_VARIABLES
                finite_array(dataset, variable, source)
            end
        end
        source_hashes[old_index] = digest
    end
    wavelength_hashes = Dict{Int,String}()
    for aerosol_index in 1:2
        full_wavelength = joinpath(
            scene_root(paths.full_truth_root, aerosol_index),
            "sim_wavelength.nc")
        bottom_wavelength = joinpath(
            scene_root(truth_root(paths), aerosol_index),
            "sim_wavelength.nc")
        isfile(full_wavelength) || error(
            "missing full-column wavelength grid: $full_wavelength")
        isfile(bottom_wavelength) || error(
            "missing bottom-layer wavelength grid: $bottom_wavelength")
        digest = file_sha256(full_wavelength)
        file_sha256(bottom_wavelength) == digest || error(
            "bottom/full wavelength grids differ for aerosol index $aerosol_index")
        wavelength_hashes[aerosol_index] = digest
    end
    return (; receipt, receipt_hash=file_sha256(receipt), full_table,
            full_table_hash=file_sha256(full_table), released_hashes,
            source_hashes, wavelength_hashes)
end

function write_stage_tables!(paths::MigrationPaths)
    resolved_target(BottomLayerTruthCommon.FULL_COLUMN_TRUTH_ROOT) ==
        paths.full_truth_root || error(
        "FULL_COLUMN_TRUTH_ROOT must be set before Julia starts so the " *
        "bottom-layer table writer and migration paths use the same source")
    root = stage_truth_root(paths)
    mkpath(root)
    writers = Dict(
        "true_states.dat" => BottomLayerTruthCommon._write_state_table,
        "control_reuse_map.dat" => BottomLayerTruthCommon._write_reuse_map,
        "scene_components.dat" => BottomLayerTruthCommon._write_component_catalog,
    )
    for name in STATIC_TABLES
        destination = joinpath(root, name)
        temporary = destination * ".tmp.$(getpid()).$(uuid4())"
        writers[name](temporary)
        if isfile(destination)
            read(destination) == read(temporary) || error(
                "existing staged table differs: $destination")
            rm(temporary)
        else
            Base.Filesystem.rename(temporary, destination)
        end
    end
    return root
end

function compare_state_tables(paths::MigrationPaths)
    legacy_path = joinpath(truth_root(paths), "true_states.dat")
    corrected_path = joinpath(stage_truth_root(paths), "true_states.dat")
    isfile(legacy_path) || error("missing legacy bottom state table: $legacy_path")
    isfile(corrected_path) || error("missing corrected staged state table: $corrected_path")
    old = readdlm(legacy_path; comments=true)
    new = readdlm(corrected_path; comments=true)
    size(old) == (80, 31) && size(new) == (80, 31) || error(
        "bottom state tables must both be 80x31")
    fixed_columns = vcat(collect(1:6), collect(8:19), collect(23:31))
    for index in 1:80
        all(old[index, column] == new[index, column] for column in fixed_columns) ||
            error("non-SIF bottom state-table fields changed in state $index")
        if index in EXPECTED_NO_SIF_STATES
            all(old[index, column] == new[index, column] for column in 20:22) ||
                error("no-SIF fields changed in bottom state $index")
            String(old[index, 7]) == "off" && String(new[index, 7]) == "off" ||
                error("no-SIF label changed in bottom state $index")
        else
            String(old[index, 7]) == LEGACY_SIF_CASE || error(
                "bottom state $index is not the expected legacy SIF case")
            String(new[index, 7]) == RRSXCO2Common.SIF_CASE_ON || error(
                "bottom state $index lacks the corrected SIF-v2 label")
        end
    end
    states = BottomLayerTruthCommon.read_bottom_states(corrected_path)
    getfield.(filter(state -> state.sif_index == 2, states), :index) ==
        EXPECTED_SIF_STATES || error("corrected table has the wrong SIF state set")
    return states
end

function copy_no_sif_dependencies!(paths::MigrationPaths, states)
    for aerosol_index in 1:2
        canonical = scene_root(truth_root(paths), aerosol_index)
        staged = scene_root(stage_truth_root(paths), aerosol_index)
        copy_verified(joinpath(canonical, "sim_wavelength.nc"),
                      joinpath(staged, "sim_wavelength.nc"))
    end
    for state in filter(state -> state.sif_index == 1, states)
        copy_verified(truth_scene(truth_root(paths), state),
                      truth_scene(stage_truth_root(paths), state))
        copy_verified(oco_scene(oco_root(paths), state.index),
                      oco_scene(stage_oco_root(paths), state.index))
        copy_verified(noise_scene(noise_root(paths), state.index),
                      noise_scene(stage_noise_root(paths), state.index))
    end
end

function relative_to_rrs(path)
    return relpath(abspath(path), normpath(joinpath(@__DIR__, "..")))
end

function assemble_truth_scene!(paths::MigrationPaths, state, corrected_table_hash;
                               force=false)
    state.sif_index == 2 || error("assemble_truth_scene! accepts only SIF-on states")
    destination = truth_scene(stage_truth_root(paths), state)
    if isfile(destination) && !force
        return destination
    end
    source = full_source_scene(paths, state)
    partner_index = BottomLayerTruthCommon.paired_no_sif_index(state)
    partner_states = BottomLayerTruthCommon.read_bottom_states(
        joinpath(stage_truth_root(paths), "true_states.dat"))
    partner = partner_states[partner_index]
    partner.sif_index == 1 || error("state $(state.index) has no no-SIF partner")
    co2_source = truth_scene(truth_root(paths), partner)
    isfile(source) && isfile(co2_source) || error(
        "missing source inputs for bottom SIF state $(state.index)")

    mkpath(dirname(destination))
    temporary = destination * ".tmp.$(getpid()).$(uuid4())"
    cp(source, temporary; force=true)
    try
        NCDataset(co2_source) do co2
            NCDataset(temporary, "a") do output
                for variable in CO2_TRUTH_VARIABLES
                    output[variable][:, :] = co2[variable][:, :]
                end
                attributes = output.attrib
                # The copied full-column container carries CO2 ancillary
                # metadata for a different profile. Replace it from the exact
                # paired bottom-layer CO2 source, or remove it if the source
                # does not define that optional field.
                for key in CO2_SOURCE_ATTRIBUTES
                    if haskey(co2.attrib, key)
                        attributes[key] = co2.attrib[key]
                    elseif haskey(attributes, key)
                        delete!(attributes, key)
                    end
                end
                attributes["campaign"] = "bottom_layer_XCO2"
                attributes["state_index"] = Int32(state.index)
                attributes["surface"] = state.surface
                attributes["aerosol_case"] = state.aerosol_case
                attributes["sif_case"] = state.sif_case
                RRSXCO2Common.write_sif_provenance!(attributes, true)
                attributes["xco2_ppm"] = state.xco2_ppm
                attributes["column_xco2_ppm"] = state.xco2_ppm
                attributes["background_co2_ppm"] = state.background_co2_ppm
                attributes["bottom_co2_ppm"] = state.bottom_co2_ppm
                attributes["bottom_co2_index"] = Int32(state.bottom_co2_index)
                attributes["bottom_co2_layer_index"] =
                    Int32(state.bottom_layer_index)
                attributes["co2_profile_order"] = "TOA-to-BOA"
                attributes["co2_profile_definition"] =
                    "layers 1:15 fixed at background_co2_ppm; layer 16 set to bottom_co2_ppm"
                profile = fill(state.background_co2_ppm, 16)
                profile[state.bottom_layer_index] = state.bottom_co2_ppm
                attributes["co2_profile_ppm"] = profile
                attributes["bottom_layer_dry_column_fraction"] =
                    BottomLayerTruthCommon.BOTTOM_DRY_COLUMN_FRACTION
                attributes["source_state_table"] = state.aerosol_index == 2 ?
                    "../true_states.dat" : "true_states.dat"
                attributes["state_table_sha256"] = corrected_table_hash
                attributes["full_column_source_state"] = Int32(
                    BottomLayerTruthCommon.old_control_index(state))
                attributes["full_column_source_scene"] = relative_to_rrs(source)
                attributes["full_column_source_sha256"] = file_sha256(source)
                attributes["full_column_state_table_sha256"] = file_sha256(
                    joinpath(paths.full_truth_root, "true_states.dat"))
                attributes["o2_truth_reused"] = Int32(1)
                attributes["o2_truth_reuse_source"] = relative_to_rrs(source)
                attributes["bottom_co2_truth_mode"] =
                    "published SIF-v2 O2 reused; weak/strong CO2 reused from paired no-SIF bottom-layer scene"
                attributes["co2_truth_reuse_source"] = relative_to_rrs(co2_source)
                attributes["co2_truth_reuse_sha256"] = file_sha256(co2_source)
                attributes["producer_script"] = relative_to_rrs(@__FILE__)
                attributes["producer_script_sha256"] = file_sha256(@__FILE__)
                attributes["bottom_layer_truth_complete"] = Int32(1)
                attributes["simulation_complete"] = Int32(1)
                haskey(attributes, "chunked_simulation_complete") &&
                    (attributes["chunked_simulation_complete"] = Int32(1))
                completed = string(now(UTC))
                attributes["created"] = completed
                attributes["completed"] = completed
            end
        end
        isfile(destination) && rm(destination)
        Base.Filesystem.rename(temporary, destination)
    finally
        isfile(temporary) && rm(temporary)
    end
    return destination
end

function julia_command(script)
    project = dirname(something(Base.active_project(),
        joinpath(normpath(joinpath(@__DIR__, "..", "..")), "Project.toml")))
    return `$(Base.julia_cmd()) --project=$project --startup-file=no $script`
end

function run_stage_process(script, environment)
    command = julia_command(script)
    run(addenv(command, environment...))
end

function bind_output_source!(path, key, canonical_source, hash_key,
                             staged_source, processor_key,
                             processor_script)
    NCDataset(path, "a") do dataset
        dataset.attrib[key] = abspath(canonical_source)
        dataset.attrib[hash_key] = file_sha256(staged_source)
        dataset.attrib[processor_key] = abspath(processor_script)
        dataset.attrib["$(processor_key)_sha256"] =
            file_sha256(processor_script)
    end
end

function build_instrument_and_noise!(paths::MigrationPaths, states)
    truth_dirs = join((stage_truth_root(paths),
                       joinpath(stage_truth_root(paths), "aerosol_chunked")), ':')
    for range in SIF_STATE_RANGES
        run_stage_process(INSTRUMENT_PROCESSOR, (
            "TRUTH_INPUT_DIRS" => truth_dirs,
            "SYNTHETIC_OCO_OUT" => stage_oco_root(paths),
            "SUPPLEMENTAL_SHOULDER_DIR" => "none",
            "OCO_STOKES_COEFFICIENTS" => paths.stokes_coefficients,
            "FIRST_STATE" => string(first(range)),
            "LAST_STATE" => string(last(range)),
            "BANDS" => "o2a,weak_co2,strong_co2",
            # SIF stage products are deliberately regenerated every time.
            # Provenance equality alone cannot prove that an existing
            # convolution was made from the current staged truth bytes.
            "FORCE" => "1",
        ))
    end
    for state in filter(state -> state.sif_index == 2, states)
        bind_output_source!(
            oco_scene(stage_oco_root(paths), state.index),
            "source_truth_scene",
            truth_scene(truth_root(paths), state),
            "source_truth_sha256",
            truth_scene(stage_truth_root(paths), state),
            "instrument_processor_script",
            INSTRUMENT_PROCESSOR)
    end

    # Regenerate every SIF noise product from the just-regenerated measurement.
    # The subsequent all-state pass skips products and only emits the complete
    # 80-state manifest, preserving copied no-SIF bytes.
    for range in SIF_STATE_RANGES
        run_stage_process(NOISE_PROCESSOR, (
            "SYNTHETIC_OCO_DIR" => stage_oco_root(paths),
            "OCO_NOISE_OUT" => stage_noise_root(paths),
            "OCO_SNR_COEFFICIENTS" => paths.snr_coefficients,
            "FIRST_STATE" => string(first(range)),
            "LAST_STATE" => string(last(range)),
            "FORCE" => "1",
        ))
    end
    for state in filter(state -> state.sif_index == 2, states)
        bind_output_source!(
            noise_scene(stage_noise_root(paths), state.index),
            "source_synthetic_measurement",
            oco_scene(oco_root(paths), state.index),
            "source_synthetic_measurement_sha256",
            oco_scene(stage_oco_root(paths), state.index),
            "noise_processor_script",
            NOISE_PROCESSOR)
    end
    run_stage_process(NOISE_PROCESSOR, (
        "SYNTHETIC_OCO_DIR" => stage_oco_root(paths),
        "OCO_NOISE_OUT" => stage_noise_root(paths),
        "OCO_SNR_COEFFICIENTS" => paths.snr_coefficients,
        "FIRST_STATE" => "1",
        "LAST_STATE" => "80",
        "FORCE" => "0",
    ))
end

function no_sif_hashes(paths::MigrationPaths, states)
    result = Dict{Tuple{Symbol,Int},String}()
    for state in filter(state -> state.sif_index == 1, states)
        result[(:truth, state.index)] = file_sha256(
            truth_scene(truth_root(paths), state))
        result[(:measurement, state.index)] = file_sha256(
            oco_scene(oco_root(paths), state.index))
        result[(:noise, state.index)] = file_sha256(
            noise_scene(noise_root(paths), state.index))
    end
    return result
end

function validate_no_sif_unchanged(paths::MigrationPaths, states, hashes)
    for state in filter(state -> state.sif_index == 1, states)
        for (kind, path) in (
                (:truth, truth_scene(truth_root(paths), state)),
                (:measurement, oco_scene(oco_root(paths), state.index)),
                (:noise, noise_scene(noise_root(paths), state.index)))
            file_sha256(path) == hashes[(kind, state.index)] || error(
                "no-SIF $kind state $(state.index) changed during migration")
        end
    end
end

function validate_truth_triplet(paths::MigrationPaths, state, source_info;
                                root=stage_truth_root(paths))
    path = truth_scene(root, state)
    source = full_source_scene(paths, state)
    partner_states = BottomLayerTruthCommon.read_bottom_states(
        joinpath(stage_truth_root(paths), "true_states.dat"))
    partner = partner_states[BottomLayerTruthCommon.paired_no_sif_index(state)]
    co2_source = truth_scene(truth_root(paths), partner)
    isfile(path) || error("missing staged bottom SIF truth: $path")
    NCDataset(path) do dataset
        Int(required_attribute(dataset.attrib, "simulation_complete", path)) == 1 ||
            error("incomplete bottom SIF truth: $path")
        Int(required_attribute(dataset.attrib, "bottom_layer_truth_complete", path)) == 1 ||
            error("bottom-layer completion marker is missing: $path")
        Int(required_attribute(dataset.attrib, "state_index", path)) == state.index ||
            error("bottom SIF truth state metadata mismatch: $path")
        String(required_attribute(dataset.attrib, "sif_case", path)) ==
            state.sif_case || error("bottom SIF truth case mismatch: $path")
        validated_sif_provenance(dataset.attrib, Symbol(state.sif_case); source=path)
        String(required_attribute(dataset.attrib, "state_table_sha256", path)) ==
            file_sha256(joinpath(stage_truth_root(paths), "true_states.dat")) ||
            error("bottom SIF truth was not assembled from the corrected table: $path")
        old = BottomLayerTruthCommon.old_control_index(state)
        String(required_attribute(dataset.attrib, "full_column_source_sha256", path)) ==
            source_info.source_hashes[old] || error(
            "bottom SIF truth has wrong full-column source hash: $path")
        String(required_attribute(dataset.attrib, "co2_truth_reuse_sha256", path)) ==
            file_sha256(co2_source) || error(
            "bottom SIF truth has wrong no-SIF CO2 source hash: $path")
        String(required_attribute(
            dataset.attrib, "producer_script_sha256", path)) ==
            file_sha256(@__FILE__) || error(
            "bottom SIF truth was assembled by a different migration script: $path")
        for variable in ALL_TRUTH_VARIABLES
            finite_array(dataset, variable, path)
        end
    end
    arrays_match(path, source, SIF_TRUTH_VARIABLES) || error(
        "bottom state $(state.index) does not exactly reuse published SIF-v2 O2")
    arrays_match(path, co2_source, CO2_TRUTH_VARIABLES) || error(
        "bottom state $(state.index) does not exactly preserve paired no-SIF CO2")
    NCDataset(path) do output
        NCDataset(co2_source) do co2
            for key in CO2_SOURCE_ATTRIBUTES
                haskey(output.attrib, key) == haskey(co2.attrib, key) || error(
                    "$path does not preserve CO2 metadata presence for $key")
                haskey(co2.attrib, key) || continue
                output.attrib[key] == co2.attrib[key] || error(
                    "$path does not preserve paired no-SIF CO2 metadata $key")
            end
        end
    end
    return path
end

function validate_measurement(paths::MigrationPaths, state;
                              root=stage_oco_root(paths))
    path = oco_scene(root, state.index)
    isfile(path) || error("missing staged bottom SIF measurement: $path")
    NCDataset(path) do dataset
        Int(required_attribute(dataset.attrib,
                               "instrument_processing_complete", path)) == 1 ||
            error("incomplete synthetic OCO product: $path")
        Int(required_attribute(dataset.attrib, "state_index", path)) == state.index ||
            error("measurement state metadata mismatch: $path")
        String(required_attribute(dataset.attrib, "source_truth_scene", path)) ==
            abspath(truth_scene(truth_root(paths), state)) || error(
            "measurement does not name its canonical bottom truth source: $path")
        String(required_attribute(dataset.attrib, "source_truth_sha256", path)) ==
            file_sha256(truth_scene(stage_truth_root(paths), state)) || error(
            "measurement is not bound to the current staged truth bytes: $path")
        String(required_attribute(
            dataset.attrib, "instrument_processor_script_sha256", path)) ==
            file_sha256(INSTRUMENT_PROCESSOR) || error(
            "measurement used an unexpected instrument processor: $path")
        coefficient_path = String(required_attribute(
            dataset.attrib, "representative_stokes_coefficients", path))
        isfile(coefficient_path) || error(
            "measurement coefficient source is unavailable: $coefficient_path")
        file_sha256(coefficient_path) ==
            file_sha256(paths.stokes_coefficients) || error(
            "measurement used unexpected Stokes coefficients: $path")
        provenance = validated_sif_provenance(
            dataset.attrib, Symbol(state.sif_case); source=path)
        truth_path = truth_scene(stage_truth_root(paths), state)
        truth_provenance = NCDataset(truth_path) do truth
            validated_sif_provenance(
                truth.attrib, Symbol(state.sif_case); source=truth_path)
        end
        matching_sif_provenance(
            truth_provenance, provenance;
            first_source=truth_path, second_source=path)
        for band in ("o2a", "weak_co2", "strong_co2")
            rayleigh = vec(finite_array_1d(dataset,
                "I_OCO_rayleigh_$(band)", path))
            corrected = vec(finite_array_1d(dataset,
                "I_OCO_corrected_$(band)", path))
            same_bits(rayleigh, corrected) || error(
                "$path has a non-Rayleigh corrected $band measurement")
            uncorrected = vec(finite_array_1d(dataset,
                "I_OCO_uncorrected_$(band)", path))
            if band == "o2a"
                cabannes = vec(finite_array_1d(dataset,
                    "I_OCO_cabannes_o2a", path))
                rrs = vec(finite_array_1d(dataset, "I_OCO_rrs_o2a", path))
                uncorrected == cabannes .+ rrs || error(
                    "$path has an invalid uncorrected O2 measurement")
            else
                same_bits(uncorrected, rayleigh) || error(
                    "$path has unequal corrected/uncorrected $band")
            end
        end
    end
    return path
end

function finite_array_1d(dataset, name, path)
    haskey(dataset, name) || error("$path lacks $name")
    raw = dataset[name][:]
    any(ismissing, raw) && error("$path contains missing $name values")
    values = Array(raw)
    all(isfinite, values) || error("$path:$name contains non-finite values")
    maximum(abs, values) < 1e30 || error("$path:$name contains fill values")
    return values
end

function validate_noise(paths::MigrationPaths, state;
                        root=stage_noise_root(paths),
                        measurement_root=stage_oco_root(paths))
    path = noise_scene(root, state.index)
    measurement_path = oco_scene(measurement_root, state.index)
    isfile(path) || error("missing staged bottom SIF noise: $path")
    measurement_provenance = NCDataset(measurement_path) do dataset
        validated_sif_provenance(
            dataset.attrib, Symbol(state.sif_case); source=measurement_path)
    end
    NCDataset(path) do dataset
        Int(required_attribute(dataset.attrib,
                               "noise_covariance_complete", path)) == 1 ||
            error("incomplete noise product: $path")
        Int(required_attribute(dataset.attrib, "state_index", path)) == state.index ||
            error("noise state metadata mismatch: $path")
        String(required_attribute(
            dataset.attrib, "source_synthetic_measurement", path)) ==
            abspath(oco_scene(oco_root(paths), state.index)) || error(
            "noise does not name its canonical synthetic source: $path")
        String(required_attribute(
            dataset.attrib, "source_synthetic_measurement_sha256", path)) ==
            file_sha256(measurement_path) || error(
            "noise is not bound to the current staged measurement bytes: $path")
        String(required_attribute(
            dataset.attrib, "noise_processor_script_sha256", path)) ==
            file_sha256(NOISE_PROCESSOR) || error(
            "noise used an unexpected covariance processor: $path")
        coefficient_path = String(required_attribute(
            dataset.attrib, "representative_snr_coefficients", path))
        isfile(coefficient_path) || error(
            "noise coefficient source is unavailable: $coefficient_path")
        file_sha256(coefficient_path) ==
            file_sha256(paths.snr_coefficients) || error(
            "noise used unexpected SNR coefficients: $path")
        provenance = validated_sif_provenance(
            dataset.attrib, Symbol(state.sif_case); source=path)
        matching_sif_provenance(
            measurement_provenance, provenance;
            first_source=measurement_path, second_source=path)
        for class in ("corrected", "uncorrected")
            measurement = finite_array_1d(dataset,
                "measurement_$(class)", path)
            sigma = finite_array_1d(dataset, "noise_std_$(class)", path)
            variance = finite_array_1d(dataset,
                "Se_diagonal_$(class)", path)
            all(variance .> 0) || error("$path has non-positive $class noise")
            variance == sigma .^ 2 || error(
                "$path has inconsistent $class variance")
            NCDataset(measurement_path) do synthetic
                expected = vcat((finite_array_1d(synthetic,
                    "I_OCO_$(class)_$(band)", measurement_path)
                    for band in ("o2a", "weak_co2", "strong_co2"))...)
                same_bits(measurement, expected) || error(
                    "$path stores the wrong $class measurement")
            end
        end
    end
    return path
end

function no_sif_triplet_set_hash(states, no_sif)
    lines = String[]
    for state in sort!(filter(state -> state.sif_index == 1, states);
                       by=state -> state.index)
        for kind in (:truth, :measurement, :noise)
            push!(lines, "$(String(kind)) $(state.index) " *
                         no_sif[(kind, state.index)])
        end
    end
    return bytes2hex(sha256(join(lines, '\n') * "\n"))
end

function input_set_hash(paths, states, source_info, no_sif,
                        legacy_table_hash)
    lines = String[
        "full_release_receipt $(source_info.receipt_hash)",
        "full_state_table $(source_info.full_table_hash)",
        "bottom_state_table $(file_sha256(joinpath(stage_truth_root(paths), "true_states.dat")))",
        "legacy_bottom_state_table $legacy_table_hash",
        "no_sif_triplet_set $(no_sif_triplet_set_hash(states, no_sif))",
        "stokes_coefficients $(file_sha256(paths.stokes_coefficients))",
        "snr_coefficients $(file_sha256(paths.snr_coefficients))",
        "migration_script $(file_sha256(@__FILE__))",
        "instrument_processor $(file_sha256(INSTRUMENT_PROCESSOR))",
        "noise_processor $(file_sha256(NOISE_PROCESSOR))",
    ]
    for (key, value) in sort!(collect(no_sif); by=first)
        push!(lines, "no_sif $(first(key)) $(last(key)) $value")
    end
    for (index, value) in sort!(collect(source_info.source_hashes); by=first)
        push!(lines, "full_source $index $value")
    end
    for (index, value) in sort!(collect(source_info.wavelength_hashes); by=first)
        push!(lines, "wavelength_grid $index $value")
    end
    return bytes2hex(sha256(join(lines, '\n') * "\n"))
end

function write_receipt(path, paths, states, source_info, no_sif, triplet_hashes,
                       legacy_table_hash; complete=false)
    input_hash = input_set_hash(
        paths, states, source_info, no_sif, legacy_table_hash)
    table_hash = file_sha256(joinpath(stage_truth_root(paths), "true_states.dat"))
    write_atomic(path) do io
        println(io, complete ?
            "# corrected bottom-layer SIF-v2 publication complete" :
            "# corrected bottom-layer SIF-v2 staging validation")
        println(io, "# release_schema $RELEASE_SCHEMA")
        println(io, "# sif_definition_version ",
                RRSXCO2Common.SIF_DEFINITION_VERSION)
        println(io, "# full_column_release_receipt_sha256 ",
                source_info.receipt_hash)
        println(io, "# full_column_state_table_sha256 ",
                source_info.full_table_hash)
        println(io, "# bottom_state_table_sha256 ", table_hash)
        println(io, "# legacy_bottom_state_table_sha256 ",
                legacy_table_hash)
        println(io, "# legacy_bottom_state_table_archive_relative ",
                joinpath("truth", "true_states.dat"))
        println(io, "# legacy_archive_manifest_relative ",
                basename(archive_manifest(paths)))
        println(io, "# legacy_archive_manifest_sha256 ",
                complete ? file_sha256(archive_manifest(paths)) : "pending")
        println(io, "# no_sif_byte_preservation_policy ",
                "preserve_truth_measurement_noise_bytes_and_legacy_state_table_hash_attributes")
        println(io, "# no_sif_triplet_set_sha256 ",
                no_sif_triplet_set_hash(states, no_sif))
        println(io, "# input_set_sha256 ", input_hash)
        println(io, "# completed_utc ",
                Dates.format(now(UTC), dateformat"yyyy-mm-ddTHH:MM:SS"), "Z")
        println(io, "# state truth_sha256 measurement_sha256 noise_sha256")
        for state in filter(state -> state.sif_index == 2, states)
            values = triplet_hashes[state.index]
            @printf(io, "%03d %s %s %s\n", state.index,
                    values.truth, values.measurement, values.noise)
        end
    end
    return path
end

function validate_release(paths::MigrationPaths=MigrationPaths();
                          write_staging_receipt=true)
    ispath(publication_marker(paths)) && error(
        "bottom SIF publication marker exists; restore explicitly before validation")
    states = compare_state_tables(paths)
    source_info = validate_full_column_release(paths,
        filter(state -> state.sif_index == 2, states))
    no_sif = no_sif_hashes(paths, states)
    legacy_table_hash = file_sha256(
        joinpath(truth_root(paths), "true_states.dat"))

    # Staged no-SIF files are dependencies/manifest inputs only. They must be
    # exact copies of the active accepted products.
    for aerosol_index in 1:2
        staged_wavelength = joinpath(
            scene_root(stage_truth_root(paths), aerosol_index),
            "sim_wavelength.nc")
        isfile(staged_wavelength) || error(
            "missing staged wavelength grid: $staged_wavelength")
        file_sha256(staged_wavelength) ==
            source_info.wavelength_hashes[aerosol_index] || error(
            "staged wavelength grid differs for aerosol index $aerosol_index")
    end
    for state in filter(state -> state.sif_index == 1, states)
        for (canonical, staged) in (
                (truth_scene(truth_root(paths), state),
                 truth_scene(stage_truth_root(paths), state)),
                (oco_scene(oco_root(paths), state.index),
                 oco_scene(stage_oco_root(paths), state.index)),
                (noise_scene(noise_root(paths), state.index),
                 noise_scene(stage_noise_root(paths), state.index)))
            isfile(staged) || error("missing staged no-SIF dependency: $staged")
            file_sha256(staged) == file_sha256(canonical) || error(
                "staged no-SIF dependency differs from canonical: $staged")
        end
    end

    triplets = Dict{Int,NamedTuple}()
    for state in filter(state -> state.sif_index == 2, states)
        truth = validate_truth_triplet(paths, state, source_info)
        measurement = validate_measurement(paths, state)
        noise = validate_noise(paths, state)
        triplets[state.index] = (
            truth=file_sha256(truth),
            measurement=file_sha256(measurement),
            noise=file_sha256(noise),
        )
    end
    length(triplets) == 40 || error("expected 40 staged SIF-v2 triplets")
    validate_no_sif_unchanged(paths, states, no_sif)
    receipt = write_staging_receipt ? write_receipt(
        stage_receipt(paths), paths, states, source_info, no_sif, triplets,
        legacy_table_hash) :
        nothing
    return (; paths, states, source_info, no_sif, triplets,
            legacy_table_hash, receipt)
end

function stage_release(paths::MigrationPaths=MigrationPaths())
    ispath(publication_marker(paths)) && error(
        "cannot stage while a bottom SIF publication marker exists")
    write_stage_tables!(paths)
    states = compare_state_tables(paths)
    validate_full_column_release(paths,
        filter(state -> state.sif_index == 2, states))
    no_sif_before = no_sif_hashes(paths, states)
    copy_no_sif_dependencies!(paths, states)
    table_hash = file_sha256(joinpath(stage_truth_root(paths), "true_states.dat"))
    for state in filter(state -> state.sif_index == 2, states)
        # Assembly is cheap relative to RT and must never reuse a composite
        # made by an older migration script or from earlier source bytes.
        assemble_truth_scene!(paths, state, table_hash; force=true)
    end
    build_instrument_and_noise!(paths, states)
    validate_no_sif_unchanged(paths, states, no_sif_before)
    return validate_release(paths)
end

function active_retrieval_products(paths, sif_states)
    result = String[]
    state_set = Set(sif_states)
    root = retrieval_root(paths)
    for class in ("corrected", "uncorrected")
        directory = joinpath(root, class)
        isdir(directory) || continue
        for name in readdir(directory)
            match_result = match(r"^retrieval_state(\d{3})_", name)
            isnothing(match_result) && continue
            parse(Int, only(match_result.captures)) in state_set || continue
            path = joinpath(directory, name)
            isfile(path) && push!(result, path)
        end
    end
    # Manifests may combine SIF and no-SIF rows, so retire them as a whole and
    # let the next retrieval launcher regenerate an authoritative manifest.
    isdir(root) && for name in readdir(root)
        startswith(name, "retrieval_manifest") && endswith(name, ".dat") || continue
        path = joinpath(root, name)
        isfile(path) && push!(result, path)
    end
    return sort!(unique(result))
end

function replacement_pairs(paths, states)
    pairs = Pair{String,String}[]
    for name in STATIC_TABLES
        push!(pairs, joinpath(stage_truth_root(paths), name) =>
                     joinpath(truth_root(paths), name))
    end
    for state in filter(state -> state.sif_index == 2, states)
        push!(pairs,
            truth_scene(stage_truth_root(paths), state) =>
                truth_scene(truth_root(paths), state),
            oco_scene(stage_oco_root(paths), state.index) =>
                oco_scene(oco_root(paths), state.index),
            noise_scene(stage_noise_root(paths), state.index) =>
                noise_scene(noise_root(paths), state.index),
        )
    end
    stage_manifest = joinpath(stage_noise_root(paths),
                              "noise_covariance_manifest.dat")
    isfile(stage_manifest) || error("missing complete staged noise manifest")
    push!(pairs, stage_manifest => joinpath(
        noise_root(paths), "noise_covariance_manifest.dat"))
    return pairs
end

function archive_relative(paths, canonical)
    relative = relpath(canonical, paths.campaign_root)
    startswith(relative, "..") && error(
        "refusing to archive a path outside the bottom campaign: $canonical")
    return joinpath(paths.archive_root, relative)
end

function archive_current!(paths, states, pairs, retrievals)
    archived = Dict{String,String}()
    for canonical in vcat(last.(pairs), retrievals)
        isfile(canonical) || error("missing canonical legacy product: $canonical")
        digest = file_sha256(canonical)
        destination = archive_relative(paths, canonical)
        copy_verified(canonical, destination; expected_hash=digest)
        archived[canonical] = digest
    end
    lines = sort(["$(relpath(path, paths.campaign_root)) $digest"
                  for (path, digest) in archived])
    manifest = archive_manifest(paths)
    temporary = manifest * ".tmp.$(getpid()).$(uuid4())"
    mkpath(dirname(manifest))
    open(temporary, "w") do io
        println(io, "# canonical_relative_path legacy_sha256")
        for line in lines
            println(io, line)
        end
    end
    if isfile(manifest)
        read(manifest) == read(temporary) || error(
            "existing bottom legacy archive manifest differs: $manifest")
        rm(temporary)
    else
        Base.Filesystem.rename(temporary, manifest)
    end
    for (canonical, digest) in archived
        file_sha256(archive_relative(paths, canonical)) == digest || error(
            "archive changed after verification: $canonical")
    end
    return archived
end

function prepare_sidecars(pairs, staged_hashes, release_id)
    sidecars = Pair{String,String}[]
    for (source, destination) in pairs
        expected = staged_hashes[source]
        file_sha256(source) == expected || error(
            "staged source changed after validation: $source")
        temporary = destination * ".sif-v2-pending.$release_id"
        ispath(temporary) && error("publication sidecar already exists: $temporary")
        copy_verified(source, temporary; expected_hash=expected)
        push!(sidecars, temporary => destination)
    end
    return sidecars
end

function restore_archived_pairs!(paths, archived)
    # Restore the legacy table first. This prevents a corrected table from
    # remaining visible while any product is being rolled back.
    paths_to_restore = sort!(collect(keys(archived)); by=path ->
        basename(path) == "true_states.dat" ? 0 : 1)
    for destination in paths_to_restore
        source = archive_relative(paths, destination)
        expected = archived[destination]
        temporary = destination * ".restore.$(uuid4())"
        copy_verified(source, temporary; expected_hash=expected)
        Base.Filesystem.rename(temporary, destination)
    end
end

function validate_published(paths, validation, pairs, staged_hashes)
    for (source, destination) in pairs
        file_sha256(destination) == staged_hashes[source] || error(
            "published destination differs from validated stage: $destination")
    end
    for state in filter(state -> state.sif_index == 2, validation.states)
        validate_truth_triplet(paths, state, validation.source_info;
                               root=truth_root(paths))
        validate_measurement(paths, state; root=oco_root(paths))
        validate_noise(paths, state; root=noise_root(paths),
                       measurement_root=oco_root(paths))
    end
    validate_no_sif_unchanged(
        paths, validation.states, validation.no_sif)
end

function publish_release(paths::MigrationPaths=MigrationPaths();
                         _test_fail_after_promotions::Union{Nothing,Int}=nothing)
    get(ENV, "CONFIRM_BOTTOM_SIF_V2_PUBLICATION", "") ==
        CONFIRM_PUBLICATION || error(
        "publication requires CONFIRM_BOTTOM_SIF_V2_PUBLICATION=$CONFIRM_PUBLICATION")
    marker = publication_marker(paths)
    lock = publication_lock(paths)
    ispath(marker) && error(
        "interrupted bottom SIF publication marker exists; restore explicitly")
    if ispath(lock)
        get(ENV, "BREAK_STALE_BOTTOM_SIF_V2_LOCK", "0") == "1" || error(
            "bottom SIF publication lock exists: $lock; verify no publisher " *
            "is active, then set BREAK_STALE_BOTTOM_SIF_V2_LOCK=1")
        rm(lock; recursive=true)
    end
    claims = joinpath(retrieval_root(paths), ".partition_claims")
    if isdir(claims) && any(name -> endswith(name, ".claim"), readdir(claims))
        error("active/stale retrieval claims exist; audit and retire them before publication: $claims")
    end
    ownership = joinpath(paths.campaign_root, EXTERNAL_OWNERSHIP_RELATIVE)
    ispath(ownership) || error(
        "external SIF ownership marker is required before bottom publication: $ownership")

    mkdir(lock)
    validation = nothing
    pairs = Pair{String,String}[]
    sidecars = Pair{String,String}[]
    parked = Pair{String,String}[]
    archived = Dict{String,String}()
    promotions_started = false
    committed = false
    try
        validation = validate_release(paths)
        pairs = replacement_pairs(paths, validation.states)
        staged_hashes = Dict(source => file_sha256(source)
                             for (source, _) in pairs)
        retrievals = active_retrieval_products(
            paths, getfield.(filter(state -> state.sif_index == 2,
                                    validation.states), :index))
        archived = archive_current!(paths, validation.states, pairs, retrievals)
        release_id = file_sha256(validation.receipt)[1:16]
        write_atomic(marker) do io
            println(io, "release_id $release_id")
            println(io, "staging_receipt_sha256 ", file_sha256(validation.receipt))
            println(io, "started_utc ",
                    Dates.format(now(UTC), dateformat"yyyy-mm-ddTHH:MM:SS"), "Z")
        end
        sidecars = prepare_sidecars(pairs, staged_hashes, release_id)

        # Park legacy retrievals in their current directories so a caught
        # failure can restore them by one rename. Verified archive copies are
        # retained permanently after a successful publication.
        for path in retrievals
            temporary = path * ".sif-v1-pending.$release_id"
            ispath(temporary) && error("legacy retrieval sidecar exists: $temporary")
            Base.Filesystem.rename(path, temporary)
            push!(parked, temporary => path)
        end

        # State table is the final commit point.
        sort!(sidecars; by=pair -> basename(last(pair)) == "true_states.dat" ? 1 : 0)
        promotion_count = 0
        for (temporary, destination) in sidecars
            promotions_started = true
            Base.Filesystem.rename(temporary, destination)
            promotion_count += 1
            promotion_count == _test_fail_after_promotions && error(
                "injected bottom SIF publication failure after $promotion_count promotions")
        end
        validate_published(paths, validation, pairs, staged_hashes)
        published_triplets = Dict(state.index => (
                truth=file_sha256(truth_scene(truth_root(paths), state)),
                measurement=file_sha256(oco_scene(oco_root(paths), state.index)),
                noise=file_sha256(noise_scene(noise_root(paths), state.index)),
            ) for state in filter(state -> state.sif_index == 2,
                                  validation.states))
        write_receipt(canonical_receipt(paths), paths, validation.states,
                      validation.source_info, validation.no_sif,
                      published_triplets, validation.legacy_table_hash;
                      complete=true)
        for (temporary, _) in parked
            isfile(temporary) && rm(temporary)
        end
        rm(marker)
        committed = true
        return validation
    catch
        if promotions_started && !isempty(archived)
            try
                restore_archived_pairs!(paths, archived)
                for (temporary, destination) in parked
                    if isfile(temporary)
                        file_sha256(temporary) == archived[destination] || error(
                            "parked legacy retrieval changed during rollback: $temporary")
                        file_sha256(destination) == archived[destination] || error(
                            "restored legacy retrieval checksum mismatch: $destination")
                        rm(temporary)
                    end
                end
            catch rollback_error
                @error("automatic bottom SIF rollback failed; immutable archive and publication marker remain",
                       exception=(rollback_error, catch_backtrace()))
            end
        elseif !isempty(parked)
            for (temporary, destination) in parked
                isfile(temporary) && Base.Filesystem.rename(temporary, destination)
            end
        end
        rethrow()
    finally
        for (temporary, _) in sidecars
            isfile(temporary) && rm(temporary)
        end
        isdir(lock) && rm(lock; recursive=true)
        committed && isfile(marker) && rm(marker)
    end
end

function parse_archive_manifest(paths)
    manifest = archive_manifest(paths)
    isfile(manifest) || error("missing bottom legacy archive manifest: $manifest")
    result = Dict{String,String}()
    for line in eachline(manifest)
        fields = split(strip(line))
        (isempty(fields) || first(fields) == "#") && continue
        length(fields) == 2 || error("malformed archive manifest line: $line")
        destination = joinpath(paths.campaign_root, fields[1])
        result[destination] = fields[2]
        file_sha256(archive_relative(paths, destination)) == fields[2] || error(
            "archived checksum mismatch: $(fields[1])")
    end
    isempty(result) && error("bottom legacy archive manifest is empty")
    return result
end

function restore_legacy(paths::MigrationPaths=MigrationPaths())
    get(ENV, "CONFIRM_BOTTOM_SIF_V2_RESTORE", "") == CONFIRM_RESTORE || error(
        "restore requires CONFIRM_BOTTOM_SIF_V2_RESTORE=$CONFIRM_RESTORE")
    marker = publication_marker(paths)
    ispath(marker) || error("restore requires an interrupted publication marker")
    isfile(marker) || error(
        "publication marker exists but is not a regular file: $marker")
    marker_values, _ = read_keyed_receipt(marker)
    release_id = get(marker_values, "release_id", "")
    occursin(r"^[0-9a-f]{16}$", release_id) || error(
        "publication marker has no valid release_id; refusing ambiguous cleanup")
    lock = publication_lock(paths)
    ispath(lock) && get(ENV, "BREAK_STALE_BOTTOM_SIF_V2_LOCK", "0") != "1" &&
        error("publication lock exists; verify no publisher is active, then set BREAK_STALE_BOTTOM_SIF_V2_LOCK=1")
    ispath(lock) && rm(lock; recursive=true)
    mkdir(lock)
    try
        archived = parse_archive_manifest(paths)
        restore_archived_pairs!(paths, archived)
        # Remove only sidecars carrying this marker's release ID, and only
        # after every canonical legacy file has been restored and verified.
        # A hard kill can leave either corrected-product sidecars prepared
        # before promotion or parked legacy-retrieval sidecars.
        for (destination, digest) in archived
            file_sha256(destination) == digest || error(
                "legacy restore validation failed: $destination")
            for suffix in (".sif-v1-pending.$release_id",
                           ".sif-v2-pending.$release_id")
                sidecar = destination * suffix
                isfile(sidecar) && rm(sidecar)
            end
        end
        isfile(canonical_receipt(paths)) && rm(canonical_receipt(paths))
        rm(marker)
    finally
        isdir(lock) && rm(lock; recursive=true)
    end
end

function main()
    action = lowercase(get(ENV, "BOTTOM_SIF_V2_ACTION", "validate"))
    paths = MigrationPaths()
    if action == "stage"
        result = stage_release(paths)
        @info("staged and validated corrected bottom-layer SIF-v2 chain without canonical mutation",
              scenes=40, receipt=result.receipt)
    elseif action == "validate"
        result = validate_release(paths)
        @info("bottom-layer SIF-v2 staging gate passed without canonical mutation",
              scenes=40, receipt=result.receipt)
    elseif action == "publish"
        publish_release(paths)
        @info("published corrected bottom-layer SIF-v2 truth/instrument/noise chain",
              scenes=40, receipt=canonical_receipt(paths))
    elseif action == "restore"
        restore_legacy(paths)
        @info("restored archived legacy bottom-layer SIF products",
              archive=paths.archive_root)
    else
        error("BOTTOM_SIF_V2_ACTION must be stage, validate, publish, or restore")
    end
end

end # module BottomLayerSIFV2Migration

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    BottomLayerSIFV2Migration.main()
end
