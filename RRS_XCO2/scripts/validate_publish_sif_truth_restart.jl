#!/usr/bin/env julia

"""
Fail-closed release gate for the directly recomputed full-column SIF-on truth.

The default action is validation only.  It checks all 32 corrected SIF-on
staging scenes before any canonical file is touched:

* all three O2 A-band Stokes arrays are finite, complete direct-RT products;
* every SIF attribute implements definition version 2,
  `2pi*L_lambda(760 nm)=0.5`;
* each four-XCO2 group has bit-identical A-band arrays;
* both CO2-band arrays are bit-identical to the immutable legacy archive and
  to the corresponding canonical no-SIF scene;
* the corrected state table differs from the archived table only in the three
  SIF fields and the SIF-on label; and
* the archive manifest and every pre-publication canonical checksum agree.

Validation writes a checksum receipt only inside the isolated restart root.
Publication is deliberately difficult to invoke: set

    SIF_RELEASE_ACTION=publish
    CONFIRM_SIF_V2_PUBLICATION=publish_32_direct_sif_scenes

The publisher prepares and verifies same-filesystem sidecars for every file
before promotion, retains the immutable archive, publishes the corrected state
table last, validates the canonical result, and automatically rolls back on a
caught error.  An interruption marker makes an unclean process death visible;
resume only after confirming that no publisher is active, with
`SIF_RELEASE_RESUME=1 SIF_RELEASE_BREAK_STALE_LOCK=1` and the same confirmation
token.  A resumed publication accepts only canonical files whose checksums are
either the archived or the validated staged versions.

Environment paths:

* `FULL_COLUMN_TRUTH_ROOT` (default `RRS_XCO2/truth_map`)
* `SIF_RESTART_ROOT` (default `<truth root>/.sif_v2_restart`)
* `SIF_OBSOLETE_ROOT` (default
  `RRS_XCO2/obsolete/sif_wavelength_integral_0p5_20260904`)
* `SIF_EXTERNAL_INPUT_MANIFEST` must name the transferred
  `external_inputs.sha256` used by the eight Gattaca producer tasks (the
  default is `<restart root>/external_inputs.sha256`).
"""

module SIFTruthReleaseGate

using Dates
using DelimitedFiles
using JLD2
using NCDatasets
using Printf
using SHA
using UUIDs
using vSmartMOM: sif_reference_state

include(joinpath(@__DIR__, "common.jl"))
using .RRSXCO2Common

export ReleasePaths, validate_release, publish_release, main

const O2_VARIABLES = (
    "radiance_rayleigh_o2a",
    "radiance_cabannes_o2a",
    "radiance_rrs_o2a",
)
const CO2_VARIABLES = (
    "radiance_rayleigh_weak_co2",
    "radiance_rayleigh_strong_co2",
)
const LEGACY_CASE = "total_0p5"
const CONFIRMATION = "publish_32_direct_sif_scenes"
const RELEASE_SCHEMA = 1
const DEFAULT_APPROVED_PRODUCER_SHA =
    "a2968e26e472d202983b8e53b1c68853bbf885ae"
const O2_REGEN_VERSION = 3
const PRODUCER_PARTITIONS = (
    (task=0, class="clear", filter="none", first=1, last=16,
     chunks=1, surface=1),
    (task=1, class="clear", filter="none", first=17, last=32,
     chunks=1, surface=2),
    (task=2, class="clear", filter="none", first=33, last=48,
     chunks=1, surface=3),
    (task=3, class="clear", filter="none", first=49, last=64,
     chunks=1, surface=4),
    (task=4, class="aerosol", filter="aerosol", first=1, last=16,
     chunks=43, surface=1),
    (task=5, class="aerosol", filter="aerosol", first=17, last=32,
     chunks=43, surface=2),
    (task=6, class="aerosol", filter="aerosol", first=33, last=48,
     chunks=43, surface=3),
    (task=7, class="aerosol", filter="aerosol", first=49, last=64,
     chunks=43, surface=4),
)
const PHYSICAL_SCALE = let corrected = RRSXCO2Common.campaign_sif_state()
    legacy = sif_reference_state(
        total_sif=0.5, reference_wavelength_nm=760)
    corrected.SIF760 / legacy.SIF760
end

struct ReleasePaths
    truth_root::String
    restart_root::String
    obsolete_root::String
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
           !(relative == ".." || startswith(relative, "..$(Base.Filesystem.path_separator)"))
end

function ReleasePaths(; truth_root=nothing, restart_root=nothing,
                        obsolete_root=nothing)
    rrs_root = normpath(joinpath(@__DIR__, ".."))
    truth = resolved_target(something(truth_root, get(
        ENV, "FULL_COLUMN_TRUTH_ROOT", joinpath(rrs_root, "truth_map"))))
    restart = resolved_target(something(restart_root, get(
        ENV, "SIF_RESTART_ROOT", joinpath(truth, ".sif_v2_restart"))))
    obsolete = resolved_target(something(obsolete_root, get(
        ENV,
        "SIF_OBSOLETE_ROOT",
        joinpath(rrs_root, "obsolete", "sif_wavelength_integral_0p5_20260904"),
    )))
    contains_path(restart, truth) && error(
        "restart root must not equal or contain canonical truth root")
    (contains_path(obsolete, truth) || contains_path(truth, obsolete)) && error(
        "obsolete archive and canonical truth root must not overlap")
    (contains_path(obsolete, restart) || contains_path(restart, obsolete)) && error(
        "obsolete archive and restart root must not overlap")
    return ReleasePaths(truth, restart, obsolete)
end

file_sha256(path::AbstractString) = open(path, "r") do io
    bytes2hex(sha256(io))
end

struct TableState
    index::Int
    surface_index::Int
    surface::String
    aerosol_index::Int
    aerosol_case::String
    sif_index::Int
    sif_case::String
    xco2_index::Int
    xco2_ppm::Int
    sif_control::Float64
    sif760::Float64
    msif::Float64
end

function read_table(path::AbstractString)
    isfile(path) || error("missing state table: $path")
    raw = readdlm(path; comments=true)
    size(raw) == (64, 27) || error("state table is not 64x27: $path")
    states = TableState[]
    for row in eachrow(raw)
        push!(states, TableState(
            Int(row[1]), Int(row[2]), String(row[3]), Int(row[4]),
            String(row[5]), Int(row[6]), String(row[7]), Int(row[8]),
            Int(row[9]), Float64(row[16]), Float64(row[17]), Float64(row[18])))
    end
    getfield.(states, :index) == collect(1:64) || error(
        "state table indices are not exactly 1:64: $path")
    return states, raw
end

aerosol_on(state::TableState) = state.aerosol_index == 2
sif_on(state::TableState) = state.sif_index == 2

function canonical_scene(paths::ReleasePaths, state::TableState)
    directory = aerosol_on(state) ?
        joinpath(paths.truth_root, "aerosol_chunked") : paths.truth_root
    return joinpath(directory, @sprintf("hiressim_%03d.nc", state.index))
end

function staged_scene(paths::ReleasePaths, state::TableState)
    directory = joinpath(paths.restart_root,
                         aerosol_on(state) ? "aerosol" : "clear")
    return joinpath(directory, @sprintf("hiressim_%03d.nc", state.index))
end

function archived_scene(paths::ReleasePaths, state::TableState)
    directory = aerosol_on(state) ?
        joinpath(paths.obsolete_root, "truth_map", "aerosol_chunked") :
        joinpath(paths.obsolete_root, "truth_map")
    return joinpath(directory, @sprintf("hiressim_%03d.nc", state.index))
end

corrected_table(paths::ReleasePaths) = joinpath(paths.restart_root, "true_states.dat")
canonical_table(paths::ReleasePaths) = joinpath(paths.truth_root, "true_states.dat")
archived_table(paths::ReleasePaths) =
    joinpath(paths.obsolete_root, "truth_map", "true_states.dat")
archive_manifest(paths::ReleasePaths) =
    joinpath(paths.obsolete_root, "full_column_input_manifest.dat")
release_receipt(paths::ReleasePaths) =
    joinpath(paths.restart_root, "sif_v2_release_validation.dat")
canonical_receipt(paths::ReleasePaths) =
    joinpath(paths.truth_root, "sif_v2_release_complete.dat")
publication_marker(paths::ReleasePaths) =
    joinpath(paths.truth_root, ".sif_v2_publication_in_progress")
publication_lock(paths::ReleasePaths) =
    joinpath(paths.restart_root, ".sif_v2_publication_lock")

producer_checkpoint(paths::ReleasePaths, partition) = joinpath(
    paths.restart_root, partition.class,
    "o2_exact_grid_$(partition.filter)_$(partition.first)_$(partition.last)_checkpoint.jld2")

function parse_producer_log(path)
    values = Dict{String,String}()
    hashes = Dict{String,String}()
    for line in eachline(path)
        stripped = strip(line)
        isempty(stripped) && continue
        hash_match = match(r"^([0-9a-f]{64})\s+(.+)$", stripped)
        if hash_match !== nothing
            digest, source = hash_match.captures
            hashes[basename(source)] = digest
            continue
        end
        fields = split(stripped, '='; limit=2)
        length(fields) == 2 && (values[fields[1]] = fields[2])
    end
    return (; values, hashes)
end

"Verify every transferred external input relative to its manifest directory."
function validate_external_input_manifest(path)
    isfile(path) || error("missing external-input manifest: $path")
    root = realpath(dirname(path))
    entries = Dict{String,String}()
    for line in eachline(path)
        isempty(strip(line)) && continue
        matched = match(r"^([0-9a-f]{64}) [ *](.+)$", line)
        matched === nothing && error(
            "malformed external-input manifest line in $path: $line")
        digest, relative = matched.captures
        isabspath(relative) && error(
            "external-input manifest path must be relative: $relative")
        normalized = normpath(relative)
        (normalized == ".." || startswith(
            normalized, "..$(Base.Filesystem.path_separator)")) && error(
            "external-input manifest path escapes its directory: $relative")
        haskey(entries, normalized) && error(
            "duplicate external-input manifest path: $relative")
        target = joinpath(root, normalized)
        isfile(target) || error(
            "missing external input listed by $path: $target")
        resolved = realpath(target)
        contains_path(root, resolved) || error(
            "external-input manifest target escapes through a symlink: $relative")
        file_sha256(resolved) == digest || error(
            "external input checksum mismatch: $target")
        entries[normalized] = digest
    end
    isempty(entries) && error("external-input manifest is empty: $path")
    return entries
end

function require_checkpoint_value(checkpoint, key, expected, path)
    actual = get(checkpoint, key, nothing)
    actual == expected || error(
        "$path has checkpoint $key=$(repr(actual)); expected $(repr(expected))")
end

"""Bind staged products to all eight completed Gattaca producer partitions."""
function validate_producer_provenance(paths::ReleasePaths;
                                      approved_sha=get(
                                          ENV,
                                          "APPROVED_SIF_TRUTH_CHECKPOINT_SHA",
                                          DEFAULT_APPROVED_PRODUCER_SHA))
    occursin(r"^[0-9a-f]{40}$", approved_sha) || error(
        "APPROVED_SIF_TRUTH_CHECKPOINT_SHA must be a full lowercase Git SHA")
    logs_root = joinpath(paths.restart_root, "logs")
    isdir(logs_root) || error("missing Gattaca producer log directory: $logs_root")
    log_paths = sort(filter(
        path -> endswith(path, ".provenance.txt"),
        joinpath.(logs_root, readdir(logs_root))))
    isempty(log_paths) && error("no Gattaca producer provenance logs were found")
    parsed_logs = Dict(path => parse_producer_log(path) for path in log_paths)
    selected_logs = String[]
    checkpoint_hashes = Dict{String,String}()
    common_hashes = Dict{String,Set{String}}(
        name => Set{String}() for name in
        ("external_inputs.sha256", "Manifest.toml",
         "lambertian_legendre_inputs.dat", "true_states.dat",
         "sif-spectra.csv"))

    current_sif = RRSXCO2Common.campaign_sif_state()
    for partition in PRODUCER_PARTITIONS
        candidates = filter(log_paths) do path
            record = parsed_logs[path]
            values = record.values
            get(values, "completed_utc", "") != "" &&
            get(values, "source_sha", "") == approved_sha &&
            get(values, "slurm_array_task_id", "") == string(partition.task) &&
            get(values, "aerosol_filter", "") == partition.filter &&
            get(values, "state_range", "") ==
                "$(partition.first):$(partition.last)" &&
            get(values, "force", "") == "0"
        end
        length(candidates) == 1 || error(
            "expected exactly one completed producer provenance log for task " *
            "$(partition.task), found $(length(candidates))")
        log = only(candidates)
        push!(selected_logs, log)
        for name in keys(common_hashes)
            digest = get(parsed_logs[log].hashes, name, "")
            length(digest) == 64 || error("$log lacks a SHA-256 record for $name")
            push!(common_hashes[name], digest)
        end

        checkpoint_path = producer_checkpoint(paths, partition)
        isfile(checkpoint_path) || error(
            "missing version-3 producer checkpoint: $checkpoint_path")
        checkpoint = JLD2.load(checkpoint_path)
        required = (
            ("version", O2_REGEN_VERSION),
            ("aerosol_case_filter", partition.filter),
            ("first_state", partition.first),
            ("last_state", partition.last),
            ("o2_chunk_points", partition.filter == "aerosol" ? 64 : 2735),
            ("float_type", "Float32"),
            ("nstreams", partition.filter == "aerosol" ? 9 : 5),
            ("psurf_hpa", 1000.0),
            ("nlayers", 16),
            ("sza_deg", Float32(30)),
            ("vza_deg", Float32(0)),
            ("relative_azimuth_deg", Float32(0)),
            ("sif_definition_version", RRSXCO2Common.SIF_DEFINITION_VERSION),
            ("sif_case_on", RRSXCO2Common.SIF_CASE_ON),
            ("sif_angular_integral_760", RRSXCO2Common.SIF_ANGULAR_INTEGRAL_760),
            ("sif_radiance_760", RRSXCO2Common.SIF_RADIANCE_760),
            ("sif760_native", current_sif.SIF760),
            ("msif_native", current_sif.mSIF),
            ("sif_template_wavelength_integral",
             current_sif.wavelength_integral),
            ("raman_shoulder_cm", 234.0),
            ("o2_core_grid_version", 2),
            ("surface_coordinate_version", 1),
            ("absco_version", RRSXCO2Common.ABSCO_VERSION),
        )
        for (key, expected) in required
            require_checkpoint_value(checkpoint, key, expected, checkpoint_path)
        end
        expected_tags = Set(
            "o2_surface$(partition.surface)_sif2_chunk$chunk"
            for chunk in 1:partition.chunks)
        Set(String.(get(checkpoint, "completed", String[]))) == expected_tags ||
            error("$checkpoint_path does not contain its exact completed v3 tag set")
        checkpoint_hashes[checkpoint_path] = file_sha256(checkpoint_path)
    end

    for (name, digests) in common_hashes
        length(digests) == 1 || error(
            "completed producer tasks disagree on the $name checksum")
    end
    local_inputs = Dict(
        "Manifest.toml" => normpath(joinpath(@__DIR__, "..", "..", "Manifest.toml")),
        "lambertian_legendre_inputs.dat" => normpath(joinpath(
            @__DIR__, "..", "surface_albedos", "lambertian_legendre_inputs.dat")),
        "true_states.dat" => corrected_table(paths),
        "sif-spectra.csv" => normpath(joinpath(
            @__DIR__, "..", "..", "src", "SIF_emission", "sif-spectra.csv")),
    )
    for (name, local_path) in local_inputs
        isfile(local_path) || error("missing local producer input: $local_path")
        only(common_hashes[name]) == file_sha256(local_path) || error(
            "producer provenance hash for $name differs from the local input")
    end
    external = get(ENV, "SIF_EXTERNAL_INPUT_MANIFEST",
                   joinpath(paths.restart_root, "external_inputs.sha256"))
    isfile(external) || error(
        "missing transferred Gattaca external-input manifest: $external; " *
        "set SIF_EXTERNAL_INPUT_MANIFEST to its verified transfer path")
    only(common_hashes["external_inputs.sha256"]) == file_sha256(external) ||
        error("external-input manifest differs from Gattaca provenance")
    external_entries = validate_external_input_manifest(external)

    provenance_lines = String["approved_sha $approved_sha"]
    append!(provenance_lines,
        "log $(basename(path)) $(file_sha256(path))" for path in selected_logs)
    append!(provenance_lines,
        "checkpoint $(relpath(path, paths.restart_root)) $digest"
        for (path, digest) in sort!(collect(checkpoint_hashes); by=first))
    for name in sort!(collect(keys(common_hashes)))
        push!(provenance_lines, "input $name $(only(common_hashes[name]))")
    end
    for (relative, digest) in sort!(collect(external_entries); by=first)
        push!(provenance_lines, "external $relative $digest")
    end
    provenance_digest = bytes2hex(sha256(
        join(sort!(provenance_lines), '\n') * "\n"))
    return (; approved_sha, selected_logs, checkpoint_hashes,
            input_hashes=Dict(name => only(values)
                              for (name, values) in common_hashes),
            external_entries,
            provenance_digest)
end

function require_attr(attributes, key, path)
    haskey(attributes, key) || error("$path lacks required attribute $key")
    return attributes[key]
end

function require_numeric_attr(attributes, key, expected, path;
                              atol=5e-15)
    value = Float64(require_attr(attributes, key, path))
    isapprox(value, Float64(expected); atol, rtol=0) || error(
        "$path has $key=$value; expected $expected")
end

function read_finite_variable(dataset, name, path; expected_points=nothing)
    haskey(dataset, name) || error("$path lacks $name")
    raw = dataset[name][:, :]
    any(ismissing, raw) && error("$path contains missing values in $name")
    values = Array(raw)
    size(values, 1) == 3 || error("$path:$name does not have three Stokes rows")
    expected_points === nothing || size(values, 2) == expected_points || error(
        "$path:$name has $(size(values, 2)) samples; expected $expected_points")
    all(isfinite, values) || error("$path contains non-finite values in $name")
    maximum(abs, values) < 1e30 || error(
        "$path:$name contains unwritten fill values")
    return values
end

function same_bits(left::AbstractArray, right::AbstractArray)
    size(left) == size(right) || return false
    eltype(left) == eltype(right) || return false
    isbitstype(eltype(left)) || return left == right
    return reinterpret(UInt8, vec(left)) == reinterpret(UInt8, vec(right))
end

function validate_corrected_table(paths::ReleasePaths)
    corrected, new = read_table(corrected_table(paths))
    legacy, old = read_table(archived_table(paths))
    current = RRSXCO2Common.campaign_sif_state()
    fixed_columns = vcat(collect(1:6), collect(8:15), collect(19:27))
    for index in 1:64
        all(old[index, column] == new[index, column] for column in fixed_columns) ||
            error("non-SIF state-table fields changed in state $index")
        state = corrected[index]
        if sif_on(state)
            legacy[index].sif_case == LEGACY_CASE || error(
                "archived state $index lacks the legacy SIF-on label")
            state.sif_case == RRSXCO2Common.SIF_CASE_ON || error(
                "corrected state $index has SIF label $(state.sif_case)")
            isapprox(state.sif_control,
                     RRSXCO2Common.SIF_ANGULAR_INTEGRAL_760;
                     atol=1e-14, rtol=0) || error(
                "corrected state $index has wrong 760-nm angular integral")
            isapprox(state.sif760, current.SIF760; atol=5e-15, rtol=0) ||
                error("corrected state $index has wrong SIF760")
            isapprox(state.msif, current.mSIF; atol=5e-16, rtol=0) ||
                error("corrected state $index has wrong mSIF")
        else
            legacy[index].sif_case == "off" && state.sif_case == "off" ||
                error("state $index has inconsistent no-SIF labels")
            all(iszero, (state.sif_control, state.sif760, state.msif)) ||
                error("corrected no-SIF state $index has nonzero SIF fields")
            all(old[index, column] == new[index, column] for column in 16:18) ||
                error("no-SIF fields changed in state $index")
        end
    end
    selected = filter(sif_on, corrected)
    length(selected) == 32 || error(
        "corrected table contains $(length(selected)) SIF-on states, expected 32")
    count(state -> !aerosol_on(state), selected) == 16 || error(
        "expected 16 clear SIF-on states")
    count(aerosol_on, selected) == 16 || error("expected 16 aerosol SIF-on states")
    return corrected, legacy, selected
end

function parse_archive_manifest(paths::ReleasePaths, selected)
    path = archive_manifest(paths)
    isfile(path) || error("missing immutable archive manifest: $path")
    scene_hashes = Dict{Int,String}()
    scalar_hashes = Dict{String,String}()
    for line in eachline(path)
        fields = split(strip(line))
        isempty(fields) && continue
        if first(fields) == "#"
            length(fields) >= 3 && endswith(fields[2], "sha256") &&
                (scalar_hashes[fields[2]] = fields[3])
            continue
        end
        length(fields) == 3 || error("malformed archive manifest line: $line")
        on = parse(Int, fields[2])
        haskey(scene_hashes, on) && error(
            "duplicate archived state $on in $path")
        scene_hashes[on] = fields[3]
    end
    expected = Set(state.index for state in selected)
    Set(keys(scene_hashes)) == expected || error(
        "archive manifest does not enumerate exactly the 32 SIF-on states")
    for state in selected
        archive = archived_scene(paths, state)
        isfile(archive) || error("missing archived scene: $archive")
        file_sha256(archive) == scene_hashes[state.index] || error(
            "immutable archive checksum mismatch for state $(state.index)")
    end
    get(scalar_hashes, "legacy_true_states_sha256", "") ==
        file_sha256(archived_table(paths)) || error(
            "archived state-table checksum differs from its manifest")
    get(scalar_hashes, "corrected_staging_true_states_sha256", "") ==
        file_sha256(corrected_table(paths)) || error(
            "corrected state-table checksum differs from the preparation manifest")
    return scene_hashes
end

function validate_sif_provenance(dataset, state::TableState, path)
    attributes = dataset.attrib
    Int(require_attr(attributes, "sif_definition_version", path)) ==
        RRSXCO2Common.SIF_DEFINITION_VERSION || error(
            "$path is not SIF definition version 2")
    String(require_attr(attributes, "sif_case", path)) ==
        RRSXCO2Common.SIF_CASE_ON || error("$path has stale SIF case metadata")
    String(require_attr(attributes, "sif_case_on_label", path)) ==
        RRSXCO2Common.SIF_CASE_ON || error("$path has stale SIF-on label")
    String(require_attr(attributes, "sif_definition", path)) ==
        "isotropic BOA radiance normalized by 2pi*L_lambda(760 nm)=0.5" ||
        error("$path has the wrong SIF definition")
    haskey(attributes, "sif_total_mW_m-2_sr-1") && error(
        "$path retains the superseded wavelength-integrated SIF attribute")
    current = RRSXCO2Common.campaign_sif_state()
    expected = (
        ("sif_reference_wavelength_nm", 760.0),
        ("sif_upwelling_solid_angle_sr", 2pi),
        ("sif_angular_integral_760_mW_m-2_nm-1", 0.5),
        ("sif_radiance_760_mW_m-2_sr-1_nm-1", current.radiance_760),
        ("sif_cosine_weighted_irradiance_760_mW_m-2_nm-1",
         pi * current.radiance_760),
        ("sif_SIF760_mW_m-2_sr-1_per_cm-1", current.SIF760),
        ("sif_mSIF_mW_m-2_sr-1_per_cm-2", current.mSIF),
        ("sif_template_wavelength_integral_mW_m-2_sr-1",
         current.wavelength_integral),
    )
    for (name, value) in expected
        require_numeric_attr(attributes, name, value, path)
    end
    Int(require_attr(attributes, "state_index", path)) == state.index || error(
        "$path has the wrong state index")
    String(require_attr(attributes, "surface", path)) == state.surface || error(
        "$path has the wrong surface label")
    String(require_attr(attributes, "aerosol_case", path)) ==
        state.aerosol_case || error("$path has the wrong aerosol label")
    Float64(require_attr(attributes, "xco2_ppm", path)) == state.xco2_ppm ||
        error("$path has the wrong XCO2 value")
    return nothing
end

# The producer's absolute source-table path may differ after a cross-host data
# transfer.  Its basename and content checksum, rather than its host path, are
# therefore the portable provenance contract.
function validate_direct_scene(paths, state; expected_o2_points=2735,
                               path=staged_scene(paths, state))
    isfile(path) || error("missing corrected direct-RT scene: $path")
    arrays = NCDataset(path, "r") do dataset
        validate_sif_provenance(dataset, state, path)
        source = String(require_attr(dataset.attrib, "source_state_table", path))
        basename(source) == "true_states.dat" || error(
            "$path has an unexpected source state-table name: $source")
        String(require_attr(dataset.attrib,
                            "source_state_table_sha256", path)) ==
            file_sha256(corrected_table(paths)) || error(
                "$path was not generated from the validated corrected state table")
        Int(require_attr(dataset.attrib, "simulation_complete", path)) == 1 ||
            error("$path is not complete")
        Int(require_attr(dataset.attrib,
                         "o2_absco_regeneration_complete", path)) == 1 ||
            error("$path has incomplete O2 regeneration")
        Int(require_attr(dataset.attrib, "o2_truth_regenerated", path)) == 1 ||
            error("$path is not marked as directly regenerated")
        Int(require_attr(dataset.attrib, "o2_truth_reused", path)) == 0 ||
            error("$path is still marked as reused O2 truth")
        String(require_attr(dataset.attrib, "o2_truth_reuse_source", path)) ==
            "none; regenerated from current model" || error(
                "$path does not carry direct-RT source provenance")
        Int(require_attr(dataset.attrib, "o2_core_grid_version", path)) == 2 ||
            error("$path does not use the exact retained O2 core grid")
        Int(require_attr(dataset.attrib, "o2_nstreams", path)) ==
            (aerosol_on(state) ? 9 : 5) || error("$path has wrong stream count")
        Int(require_attr(dataset.attrib, "o2_chunk_points", path)) ==
            (aerosol_on(state) ? 64 : expected_o2_points) || error(
                "$path has wrong retained chunk size")
        require_numeric_attr(dataset.attrib, "o2_raman_shoulder_cm-1",
                             234.0, path; atol=1e-12)
        String(require_attr(dataset.attrib, "spectroscopy_database", path)) ==
            "ABSCO" || error("$path does not use ABSCO")
        String(require_attr(dataset.attrib, "spectroscopy_version", path)) ==
            "5.2" || error("$path does not use ABSCO 5.2")
        Float64(require_attr(dataset.attrib, "o2_vmr", path)) == 0.21 ||
            error("$path has wrong O2 VMR")
        Float64(require_attr(dataset.attrib, "psurf_hpa", path)) == 1000.0 ||
            error("$path has wrong surface pressure")
        Int(require_attr(dataset.attrib, "atmospheric_layers", path)) == 16 ||
            error("$path has wrong layer count")
        aerosol_on(state) &&
            Int(require_attr(dataset.attrib,
                             "chunked_simulation_complete", path)) != 1 &&
            error("$path has incomplete aerosol chunking")
        Tuple(read_finite_variable(dataset, name, path;
                                   expected_points=expected_o2_points)
              for name in O2_VARIABLES)
    end
    return arrays
end

function co2_arrays(path)
    isfile(path) || error("missing CO2 comparison scene: $path")
    return NCDataset(path, "r") do dataset
        Tuple(read_finite_variable(dataset, name, path) for name in CO2_VARIABLES)
    end
end

function o2_arrays(path; expected_o2_points=2735)
    return NCDataset(path, "r") do dataset
        Tuple(read_finite_variable(dataset, name, path;
                                   expected_points=expected_o2_points)
              for name in O2_VARIABLES)
    end
end

function validate_preserved_co2(paths, state, off_state)
    stage = co2_arrays(staged_scene(paths, state))
    archive = co2_arrays(archived_scene(paths, state))
    off = co2_arrays(canonical_scene(paths, off_state))
    for (name, staged, archived, no_sif) in
            zip(CO2_VARIABLES, stage, archive, off)
        same_bits(staged, archived) || error(
            "state $(state.index) changed archived $name during O2 regeneration")
        same_bits(staged, no_sif) || error(
            "$name differs between SIF-on state $(state.index) and matching " *
            "no-SIF state $(off_state.index)")
    end
end

function validate_physical_rescaling(paths, state, off_state;
                                     expected_o2_points=2735,
                                     relative_tolerance=0.05)
    corrected = o2_arrays(staged_scene(paths, state);
                          expected_o2_points)
    legacy = o2_arrays(archived_scene(paths, state);
                       expected_o2_points)
    off = o2_arrays(canonical_scene(paths, off_state);
                    expected_o2_points)
    for (name, new, old, zero_sif) in zip(O2_VARIABLES, corrected, legacy, off)
        old_delta = Float64.(old) .- Float64.(zero_sif)
        new_delta = Float64.(new) .- Float64.(zero_sif)
        denominator = sum(abs2, old_delta)
        denominator > 0 || error(
            "legacy $name has no measurable SIF response in state $(state.index)")
        ratio = sum(new_delta .* old_delta) / denominator
        isapprox(ratio, PHYSICAL_SCALE; rtol=relative_tolerance, atol=0) ||
            error("state $(state.index) $name SIF-response scale is $ratio; " *
                  "expected approximately $PHYSICAL_SCALE")
        residual = new_delta .- PHYSICAL_SCALE .* old_delta
        relative_l2 = sqrt(sum(abs2, residual) / sum(abs2, new_delta))
        relative_l2 <= relative_tolerance || error(
            "state $(state.index) $name fails SIF linearity cross-check: " *
            "relative L2=$relative_l2")
    end
end

function validate_canonical_archive_preflight(paths, selected, scene_hashes;
                                              allow_mixed=false)
    old_table_hash = file_sha256(archived_table(paths))
    new_table_hash = file_sha256(corrected_table(paths))
    canonical_hash = file_sha256(canonical_table(paths))
    valid_table_hashes = allow_mixed ? (old_table_hash, new_table_hash) :
                                      (old_table_hash,)
    canonical_hash in valid_table_hashes || error(
        "canonical state table is neither the archived nor validated corrected table")
    for state in selected
        destination = canonical_scene(paths, state)
        isfile(destination) || error("missing canonical scene: $destination")
        canonical = file_sha256(destination)
        if allow_mixed
            canonical in (scene_hashes[state.index],
                          file_sha256(staged_scene(paths, state))) || error(
                "canonical state $(state.index) has an unknown checksum")
        else
            canonical == scene_hashes[state.index] || error(
                "canonical state $(state.index) differs from its immutable archive; " *
                "refusing publication")
        end
    end
end

function write_atomic(writer::Function, path)
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

function write_validation_receipt(paths, selected, scene_hashes, staged_hashes,
                                  corrected_table_hash, legacy_table_hash,
                                  producer)
    path = release_receipt(paths)
    write_atomic(path) do io
        println(io, "# corrected SIF truth release validation")
        println(io, "# release_schema $RELEASE_SCHEMA")
        println(io, "# validated_utc ",
                Dates.format(now(UTC), dateformat"yyyy-mm-ddTHH:MM:SS"), "Z")
        println(io, "# sif_definition_version ", RRSXCO2Common.SIF_DEFINITION_VERSION)
        println(io, "# sif_case ", RRSXCO2Common.SIF_CASE_ON)
        println(io, "# corrected_state_table_sha256 ", corrected_table_hash)
        println(io, "# archived_state_table_sha256 ", legacy_table_hash)
        println(io, "# approved_producer_git_sha ", producer.approved_sha)
        println(io, "# producer_provenance_sha256 ",
                producer.provenance_digest)
        println(io, "# index class staged_sha256 archived_sha256 canonical_destination")
        for state in selected
            @printf(io, "%03d %s %s %s %s\n",
                    state.index, aerosol_on(state) ? "aerosol" : "clear",
                    staged_hashes[state.index],
                    scene_hashes[state.index], canonical_scene(paths, state))
        end
    end
    return path
end

function validate_release(paths::ReleasePaths=ReleasePaths();
                          expected_o2_points=2735,
                          require_legacy_canonical=true,
                          write_receipt=true,
                          require_producer_provenance=true)
    corrected, legacy, selected = validate_corrected_table(paths)
    scene_hashes = parse_archive_manifest(paths, selected)
    corrected_table_hash = file_sha256(corrected_table(paths))
    legacy_table_hash = file_sha256(archived_table(paths))
    producer = require_producer_provenance ?
        validate_producer_provenance(paths) :
        (approved_sha="synthetic-test-bypass",
         provenance_digest="synthetic-test-bypass")
    staged_hashes = Dict(
        state.index => file_sha256(staged_scene(paths, state))
        for state in selected)
    arrays = Dict{Int,Any}()
    for state in selected
        arrays[state.index] = validate_direct_scene(
            paths, state; expected_o2_points)
        off_state = legacy[state.index - 4]
        off_state.sif_index == 1 || error(
            "state $(state.index) lacks its expected no-SIF partner")
        validate_preserved_co2(paths, state, off_state)
        validate_physical_rescaling(
            paths, state, off_state; expected_o2_points)
    end
    groups = Dict{Tuple{Int,Int},Vector{TableState}}()
    for state in selected
        push!(get!(groups, (state.surface_index, state.aerosol_index),
                   TableState[]), state)
    end
    length(groups) == 8 || error("expected eight surface/aerosol SIF-on groups")
    for (key, members) in groups
        sort!(members; by=state -> state.xco2_index)
        getfield.(members, :xco2_index) == collect(1:4) || error(
            "group $key does not contain all four XCO2 cases")
        reference = arrays[first(members).index]
        for member in members[2:end], (name, left, right) in
                zip(O2_VARIABLES, reference, arrays[member.index])
            same_bits(left, right) || error(
                "$name is not bit-identical across XCO2-only group $key")
        end
    end
    require_legacy_canonical && validate_canonical_archive_preflight(
        paths, selected, scene_hashes)
    file_sha256(corrected_table(paths)) == corrected_table_hash || error(
        "corrected state table changed while the release gate was running")
    for state in selected
        file_sha256(staged_scene(paths, state)) == staged_hashes[state.index] ||
            error("staged state $(state.index) changed while the release gate was running")
    end
    receipt = write_receipt ?
        write_validation_receipt(paths, selected, scene_hashes, staged_hashes,
                                 corrected_table_hash, legacy_table_hash,
                                 producer) : nothing
    return (; corrected, legacy, selected, scene_hashes, staged_hashes,
            corrected_table_hash, legacy_table_hash, producer, receipt)
end

function verified_sidecar(source, destination, release_id;
                          expected_hash=file_sha256(source))
    file_sha256(source) == expected_hash || error(
        "source changed after validation: $source")
    temporary = destination * ".sif-v2-pending.$release_id"
    ispath(temporary) && error("publication sidecar already exists: $temporary")
    cp(source, temporary)
    file_sha256(temporary) == expected_hash || error(
        "sidecar checksum mismatch: $source -> $temporary")
    return temporary
end

function atomic_replace(source, destination)
    dirname(source) == dirname(destination) || error(
        "atomic replacement requires a same-directory sidecar")
    Base.Filesystem.rename(source, destination)
    return destination
end

function restore_archive!(paths, selected, scene_hashes, legacy_table_hash)
    release_id = "rollback-$(getpid())-$(uuid4())"
    pairs = Tuple{String,String}[]
    try
        # Rollback reverses the publication commit order: restore the legacy
        # table first, so a corrected table is never visible while scenes are
        # being returned to their archived versions.
        push!(pairs, (verified_sidecar(
            archived_table(paths), canonical_table(paths), release_id;
            expected_hash=legacy_table_hash),
            canonical_table(paths)))
        for state in selected
            push!(pairs, (verified_sidecar(
                archived_scene(paths, state), canonical_scene(paths, state),
                release_id; expected_hash=scene_hashes[state.index]),
                canonical_scene(paths, state)))
        end
        for (temporary, destination) in pairs
            atomic_replace(temporary, destination)
        end
    finally
        for (temporary, _) in pairs
            isfile(temporary) && rm(temporary)
        end
    end
end

function validate_published(paths, validation; expected_o2_points=2735)
    file_sha256(canonical_table(paths)) == validation.corrected_table_hash ||
        error("published canonical state table differs from validated table")
    for state in validation.selected
        destination = canonical_scene(paths, state)
        file_sha256(destination) == validation.staged_hashes[state.index] ||
            error("published state $(state.index) differs from staged source")
        validate_direct_scene(
            paths, state; expected_o2_points, path=destination)
    end
end

function publish_release(paths::ReleasePaths=ReleasePaths();
                         expected_o2_points=2735,
                         _test_fail_after_scenes::Union{Nothing,Int}=nothing)
    get(ENV, "CONFIRM_SIF_V2_PUBLICATION", "") == CONFIRMATION || error(
        "publication requires CONFIRM_SIF_V2_PUBLICATION=$CONFIRMATION")
    marker = publication_marker(paths)
    lock = publication_lock(paths)
    resume = get(ENV, "SIF_RELEASE_RESUME", "0") == "1"
    resume && !ispath(marker) && error(
        "SIF_RELEASE_RESUME=1 requires an existing interruption marker")
    if isdir(lock)
        resume && get(ENV, "SIF_RELEASE_BREAK_STALE_LOCK", "0") == "1" ||
            error("publication lock exists: $lock")
        rm(lock; recursive=true)
    end
    mkdir(lock)
    validation = nothing
    sidecars = Tuple{String,String}[]
    committed = false
    promotions_started = false
    try
        if ispath(marker)
            isfile(marker) || error(
                "publication marker exists but is not a regular file: $marker")
            resume || error(
                "an interrupted publication marker exists: $marker; audit and resume explicitly")
            # No unknown canonical content is accepted during recovery.
            corrected, _, selected = validate_corrected_table(paths)
            scene_hashes = parse_archive_manifest(paths, selected)
            validate_canonical_archive_preflight(
                paths, selected, scene_hashes; allow_mixed=true)
        end
        # Staging is revalidated from the data, never merely trusted because a
        # previous receipt exists. During resume canonical files may be mixed.
        validation = validate_release(
            paths; expected_o2_points,
            require_legacy_canonical=!resume, write_receipt=true)
        release_id = file_sha256(validation.receipt)[1:16]
        write_atomic(marker) do io
            println(io, "release_id $release_id")
            println(io, "validation_receipt_sha256 ", file_sha256(validation.receipt))
            println(io, "started_utc ",
                    Dates.format(now(UTC), dateformat"yyyy-mm-ddTHH:MM:SS"), "Z")
        end
        for state in validation.selected
            destination = canonical_scene(paths, state)
            push!(sidecars, (verified_sidecar(
                staged_scene(paths, state), destination, release_id;
                expected_hash=validation.staged_hashes[state.index]), destination))
        end
        push!(sidecars, (verified_sidecar(
            corrected_table(paths), canonical_table(paths), release_id;
            expected_hash=validation.corrected_table_hash),
            canonical_table(paths)))

        # Promote scenes first and the authoritative table last. The marker
        # remains present throughout; downstream work must not begin while it
        # exists. Every individual replacement is a same-filesystem rename.
        validate_canonical_archive_preflight(
            paths, validation.selected, validation.scene_hashes;
            allow_mixed=resume)
        scene_promotions = 0
        for (temporary, destination) in sidecars
            promotions_started = true
            atomic_replace(temporary, destination)
            if destination != canonical_table(paths)
                scene_promotions += 1
                scene_promotions == _test_fail_after_scenes && error(
                    "injected release failure after $scene_promotions scenes")
            end
        end
        validate_published(paths, validation; expected_o2_points)
        write_atomic(canonical_receipt(paths)) do io
            println(io, "# corrected SIF truth publication complete")
            println(io, "release_schema $RELEASE_SCHEMA")
            println(io, "sif_definition_version ",
                    RRSXCO2Common.SIF_DEFINITION_VERSION)
            println(io, "release_id $release_id")
            println(io, "completed_utc ",
                    Dates.format(now(UTC), dateformat"yyyy-mm-ddTHH:MM:SS"), "Z")
            println(io, "validation_receipt_sha256 ", file_sha256(validation.receipt))
            println(io, "corrected_state_table_sha256 ",
                    validation.corrected_table_hash)
            println(io, "approved_producer_git_sha ",
                    validation.producer.approved_sha)
            println(io, "producer_provenance_sha256 ",
                    validation.producer.provenance_digest)
        end
        rm(marker)
        committed = true
        return validation
    catch exception
        if validation !== nothing && promotions_started
            try
                restore_archive!(paths, validation.selected,
                                 validation.scene_hashes,
                                 validation.legacy_table_hash)
            catch rollback_exception
                @error("automatic rollback failed; immutable archive remains intact",
                       rollback_exception)
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

function main()
    action = lowercase(get(ENV, "SIF_RELEASE_ACTION", "validate"))
    paths = ReleasePaths()
    if action == "validate"
        result = validate_release(paths)
        @info("corrected SIF truth release gate passed without canonical mutation",
              scenes=length(result.selected), receipt=result.receipt)
    elseif action == "publish"
        result = publish_release(paths)
        @info("published corrected SIF truth after complete release gate",
              scenes=length(result.selected), receipt=canonical_receipt(paths))
    else
        error("SIF_RELEASE_ACTION must be validate or publish")
    end
end

end # module SIFTruthReleaseGate

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    SIFTruthReleaseGate.main()
end
