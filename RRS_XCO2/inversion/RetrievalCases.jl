module RetrievalCases

using NCDatasets
using Printf
using Random
using SHA

export TruthCase,
       RetrievalExperiment,
       MeasurementRealization,
       default_truth_table,
       default_measurement_directory,
       default_noise_directory,
       read_truth_cases,
       read_no_sif_truth_cases,
       select_sif_truth_cases,
       build_experiments,
       external_sif_ownership_marker,
       enforce_sif_ownership,
       require_sif_release_barrier,
       paired_uniform_draw,
       validated_sif_provenance,
       matching_sif_provenance,
       SIF_PROVENANCE_ATTRIBUTES,
       load_measurement_realization,
       write_experiment_manifest,
       PERTURBED_REALIZATIONS,
       UNPERTURBED_INDEX

const INVERSION_ROOT = @__DIR__
const RRS_ROOT = normpath(joinpath(INVERSION_ROOT, ".."))
const DEFAULT_MASTER_SEED = UInt64(0x4f434f5f525253)
const MEASUREMENT_CLASSES = (:corrected, :uncorrected)
const PERTURBED_REALIZATIONS = 10
const UNPERTURBED_INDEX = 11
const EXTERNAL_SIF_OWNERSHIP_MARKER = joinpath(
    ".control", "sif_owned_externally")
const REQUIRED_SIF_CASE_ON = :angular_integral760_0p5
const REQUIRED_SIF_DEFINITION_VERSION = 2
const REQUIRED_SIF_REFERENCE_WAVELENGTH_NM = 760.0
const REQUIRED_SIF_UPWELLING_SOLID_ANGLE_SR = 2π
const REQUIRED_SIF_ANGULAR_INTEGRAL_760 = 0.5
const REQUIRED_SIF_DEFINITION =
    "isotropic BOA radiance normalized by 2pi*L_lambda(760 nm)=0.5"
const REQUIRED_SIF_RADIANCE_760 =
    REQUIRED_SIF_ANGULAR_INTEGRAL_760 /
    REQUIRED_SIF_UPWELLING_SOLID_ANGLE_SR
const FULL_SIF_RELEASE_MARKER = ".sif_v2_publication_in_progress"
const FULL_SIF_RELEASE_RECEIPT = "sif_v2_release_complete.dat"
const BOTTOM_SIF_RELEASE_MARKER =
    ".bottom_layer_sif_v2_publication_in_progress"
const BOTTOM_SIF_RELEASE_RECEIPT =
    "bottom_layer_sif_v2_release_complete.dat"
const SIF_PROVENANCE_ATTRIBUTES = (
    "sif_definition_version",
    "sif_definition",
    "sif_case_on_label",
    "sif_reference_wavelength_nm",
    "sif_upwelling_solid_angle_sr",
    "sif_angular_integral_760_mW_m-2_nm-1",
    "sif_radiance_760_mW_m-2_sr-1_nm-1",
    "sif_cosine_weighted_irradiance_760_mW_m-2_nm-1",
    "sif_SIF760_mW_m-2_sr-1_per_cm-1",
    "sif_mSIF_mW_m-2_sr-1_per_cm-2",
    "sif_template_wavelength_integral_mW_m-2_sr-1",
)

default_truth_table() = joinpath(RRS_ROOT, "truth_map", "true_states.dat")
default_measurement_directory() = joinpath(
    RRS_ROOT, "truth_map", "OCO_radiances")
default_noise_directory() = joinpath(
    default_measurement_directory(), "noise_covariances")

"""Return the marker that delegates SIF-on work away from an output tree."""
function external_sif_ownership_marker(output_root::AbstractString)
    configured = get(
        ENV, "RETRIEVAL_EXTERNAL_SIF_OWNERSHIP_MARKER",
        joinpath(output_root, EXTERNAL_SIF_OWNERSHIP_MARKER))
    return abspath(configured)
end

"""
    enforce_sif_ownership(output_root, experiments)

Refuse a selected SIF-on retrieval when the output tree contains the durable
external-ownership marker.  Host launchers start a fresh Julia runner for each
block, so this guard safely disarms their later SIF phase without interrupting
an already-running no-SIF block.  An external worker remains unaffected when
it uses its own isolated output root.
"""
function enforce_sif_ownership(output_root::AbstractString, experiments)
    any(experiment -> experiment.truth.sif_case != :off, experiments) ||
        return nothing
    marker = external_sif_ownership_marker(output_root)
    ispath(marker) || return nothing
    isfile(marker) || error(
        "SIF ownership marker exists but is not a regular file: $marker")
    owner = strip(read(marker, String))
    isempty(owner) && (owner = "external owner not recorded")
    error(
        "SIF-on retrievals are delegated outside output root " *
        "$(abspath(output_root)); ownership marker=$marker; $owner")
end

_sha256_file(path::AbstractString) = open(path, "r") do io
    bytes2hex(sha256(io))
end

function _required_receipt_hash(values, key, receipt)
    value = get(values, key, "")
    occursin(r"^[0-9a-f]{64}$", value) || error(
        "$receipt lacks a valid $key")
    return value
end

function _release_receipt(path)
    isfile(path) || error("missing required corrected-SIF release receipt: $path")
    values = Dict{String,String}()
    rows = Vector{Vector{String}}()
    for line in eachline(path)
        fields = split(strip(line))
        isempty(fields) && continue
        if first(fields) == "#"
            length(fields) == 3 && (values[fields[2]] = fields[3])
        elseif length(fields) == 2
            values[fields[1]] = fields[2]
        else
            push!(rows, fields)
        end
    end
    return values, rows
end

"""
    require_sif_release_barrier(truth_table, cases,
                                measurement_directory, noise_directory)

Refuse every SIF-on retrieval unless the corresponding all-scene publication
transaction is complete and its receipt still matches the authoritative state
table and current truth/measurement/noise files.  This is independent of the
external-ownership guard: ownership prevents two writers, while this barrier
prevents any writer from consuming a partial or superseded SIF campaign.
"""
function require_sif_release_barrier(
        truth_table::AbstractString, cases,
        measurement_directory::AbstractString,
        noise_directory::AbstractString)
    selected = filter(case -> case.sif_case != :off, cases)
    isempty(selected) && return nothing
    modes = unique(case.campaign for case in selected)
    length(modes) == 1 || error(
        "one SIF retrieval selection cannot mix truth campaigns")
    truth_root = dirname(abspath(truth_table))

    if only(modes) == :bottom_layer_XCO2
        marker = joinpath(truth_root, BOTTOM_SIF_RELEASE_MARKER)
        ispath(marker) && error(
            "bottom-layer corrected-SIF publication is in progress: $marker")
        receipt = joinpath(truth_root, BOTTOM_SIF_RELEASE_RECEIPT)
        values, rows = _release_receipt(receipt)
        get(values, "release_schema", "") == "1" || error(
            "unsupported bottom-layer SIF release schema in $receipt")
        get(values, "sif_definition_version", "") ==
            string(REQUIRED_SIF_DEFINITION_VERSION) || error(
            "bottom-layer receipt is not corrected SIF definition version 2")
        get(values, "bottom_state_table_sha256", "") ==
            _sha256_file(truth_table) || error(
            "bottom-layer truth table differs from its SIF release receipt")
        _required_receipt_hash(
            values, "legacy_bottom_state_table_sha256", receipt)
        _required_receipt_hash(
            values, "legacy_archive_manifest_sha256", receipt)
        _required_receipt_hash(
            values, "no_sif_triplet_set_sha256", receipt)
        get(values, "legacy_bottom_state_table_archive_relative", "") ==
            joinpath("truth", "true_states.dat") || error(
            "bottom-layer receipt has an unexpected legacy-table archive identity")
        get(values, "no_sif_byte_preservation_policy", "") ==
            "preserve_truth_measurement_noise_bytes_and_legacy_state_table_hash_attributes" ||
            error("bottom-layer receipt lacks the no-SIF byte-preservation policy")
        _required_receipt_hash(values, "input_set_sha256", receipt)
        full_truth_root = abspath(get(
            ENV, "FULL_COLUMN_TRUTH_ROOT",
            joinpath(@__DIR__, "..", "truth_map")))
        full_marker = joinpath(full_truth_root, FULL_SIF_RELEASE_MARKER)
        ispath(full_marker) && error(
            "full-column corrected-SIF publication is in progress: $full_marker")
        full_receipt = joinpath(full_truth_root, FULL_SIF_RELEASE_RECEIPT)
        isfile(full_receipt) || error(
            "missing published full-column SIF source receipt: $full_receipt")
        _sha256_file(full_receipt) == _required_receipt_hash(
            values, "full_column_release_receipt_sha256", receipt) || error(
            "bottom-layer receipt no longer matches its full-column source release")
        hashes = Dict{Int,NTuple{3,String}}()
        for fields in rows
            length(fields) == 4 || error(
                "malformed bottom-layer SIF release row in $receipt")
            index = parse(Int, fields[1])
            haskey(hashes, index) && error(
                "duplicate state $index in bottom-layer SIF release receipt")
            hashes[index] = (fields[2], fields[3], fields[4])
        end
        expected = sort([case.state_index for case in read_truth_cases(truth_table)
                         if case.sif_case != :off])
        sort!(collect(keys(hashes))) == expected || error(
            "bottom-layer release receipt does not enumerate every SIF-on state")
        all_cases = Dict(case.state_index => case
                         for case in read_truth_cases(truth_table))
        for index in expected
            truth = all_cases[index]
            truth_directory = truth.aerosol_case == :none ? truth_root :
                joinpath(truth_root, "aerosol_chunked")
            truth_path = joinpath(
                truth_directory, @sprintf("hiressim_%03d.nc", index))
            measurement_path = joinpath(
                measurement_directory, @sprintf("OCO2sims_%03d.nc", index))
            noise_path = joinpath(
                noise_directory, @sprintf("OCO2noise_%03d.nc", index))
            expected_hashes = hashes[index]
            for (kind, path, expected_hash) in (
                    ("truth", truth_path, expected_hashes[1]),
                    ("measurement", measurement_path, expected_hashes[2]),
                    ("noise", noise_path, expected_hashes[3]))
                isfile(path) || error("missing released SIF $kind product: $path")
                _sha256_file(path) == expected_hash || error(
                    "released SIF $kind state $index differs from its receipt")
            end
        end
        return receipt
    end

    only(modes) == :full_column_XCO2 || error(
        "unsupported SIF retrieval campaign $(only(modes))")
    marker = joinpath(truth_root, FULL_SIF_RELEASE_MARKER)
    ispath(marker) && error(
        "full-column corrected-SIF publication is in progress: $marker")
    receipt = joinpath(truth_root, FULL_SIF_RELEASE_RECEIPT)
    values, _ = _release_receipt(receipt)
    get(values, "release_schema", "") == "1" || error(
        "unsupported full-column SIF release schema in $receipt")
    get(values, "sif_definition_version", "") ==
        string(REQUIRED_SIF_DEFINITION_VERSION) || error(
        "full-column receipt is not corrected SIF definition version 2")
    get(values, "corrected_state_table_sha256", "") ==
        _sha256_file(truth_table) || error(
        "full-column truth table differs from its SIF release receipt")
    validation_receipt = abspath(get(
        ENV, "FULL_COLUMN_SIF_VALIDATION_RECEIPT",
        joinpath(truth_root, ".sif_v2_restart",
                 "sif_v2_release_validation.dat")))
    isfile(validation_receipt) || error(
        "missing full-column SIF validation receipt: $validation_receipt")
    _sha256_file(validation_receipt) == _required_receipt_hash(
        values, "validation_receipt_sha256", receipt) || error(
        "full-column validation receipt differs from publication receipt")
    _, validation_rows = _release_receipt(validation_receipt)
    hashes = Dict{Int,String}()
    for fields in validation_rows
        length(fields) == 5 || error(
            "malformed full-column SIF validation row in $validation_receipt")
        index = parse(Int, fields[1])
        haskey(hashes, index) && error(
            "duplicate state $index in full-column SIF validation receipt")
        hashes[index] = fields[3]
    end
    all_cases = read_truth_cases(truth_table)
    expected = sort([case.state_index for case in all_cases
                     if case.sif_case != :off])
    sort!(collect(keys(hashes))) == expected || error(
        "full-column validation receipt does not enumerate every SIF-on state")
    cases_by_index = Dict(case.state_index => case for case in all_cases)
    for index in expected
        truth = cases_by_index[index]
        directory = truth.aerosol_case == :none ? truth_root :
            joinpath(truth_root, "aerosol_chunked")
        path = joinpath(directory, @sprintf("hiressim_%03d.nc", index))
        isfile(path) || error("missing released full-column SIF truth: $path")
        _sha256_file(path) == hashes[index] || error(
            "released full-column SIF truth state $index differs from its receipt")
    end
    return receipt
end

"""The identifying metadata for one selected truth scene."""
struct TruthCase
    state_index::Int
    surface_index::Int
    surface::Symbol
    aerosol_index::Int
    aerosol_case::Symbol
    sif_case::Symbol
    xco2_index::Int
    xco2_ppm::Float64
    campaign::Symbol
    co2_profile_mode::Symbol
    fixed_upper_co2_ppm::Float64
    background_co2_ppm::Float64
    bottom_layer_index::Int
    bottom_co2_ppm::Float64
end

# Preserve the original public constructor for callers describing a uniform
# full-column truth scene. Bottom-layer scenes are constructed by the table
# parser below because they require additional profile metadata.
TruthCase(state_index::Int, surface_index::Int, surface::Symbol,
          aerosol_index::Int, aerosol_case::Symbol, sif_case::Symbol,
          xco2_index::Int, xco2_ppm::Real) = TruthCase(
    state_index, surface_index, surface, aerosol_index, aerosol_case,
    sif_case, xco2_index, Float64(xco2_ppm), :full_column_XCO2,
    :uniform_column, Float64(xco2_ppm), Float64(xco2_ppm), 0,
    Float64(xco2_ppm))

"""One retrieval solve from a paired corrected/uncorrected experiment."""
struct RetrievalExperiment
    retrieval_index::Int
    pair_index::Int
    truth::TruthCase
    noise_index::Int
    measurement_class::Symbol
    random_seed::UInt64
    measurement_path::String
    noise_path::String
end

"""Frozen measurement, covariance, and paired perturbation for one experiment."""
struct MeasurementRealization
    noiseless::Vector{Float64}
    perturbed::Vector{Float64}
    noise_std::Vector{Float64}
    variance::Vector{Float64}
    normalized_draw::Vector{Float64}
    wavelength_nm::Vector{Float64}
    band_ranges::Vector{UnitRange{Int}}
    provenance::Dict{String,Any}
end

# Preserve the original constructor for no-SIF callers and lightweight tests.
# A SIF-on retrieval is rejected later unless the versioned provenance is
# supplied by `load_measurement_realization`.
MeasurementRealization(noiseless, perturbed, noise_std, variance,
                       normalized_draw, wavelength_nm, band_ranges) =
    MeasurementRealization(
        noiseless, perturbed, noise_std, variance, normalized_draw,
        wavelength_nm, band_ranges, Dict{String,Any}())

function _required_attribute(attributes, key, source)
    haskey(attributes, key) || error(
        "$source is missing required corrected-SIF provenance attribute $key")
    return attributes[key]
end

function _require_close(value, expected, key, source;
                        atol=1e-14, rtol=1e-12)
    numeric = Float64(value)
    isfinite(numeric) && isapprox(numeric, expected; atol, rtol) || error(
        "$source has $key=$numeric; expected $expected for corrected SIF")
    return numeric
end

"""
    validated_sif_provenance(attributes, sif_case; source="input")

Return a portable copy of the version-2 SIF normalization record.  Every
SIF-on retrieval input must state the campaign convention
`2π Lλ(760 nm) = 0.5 mW m⁻² nm⁻¹` and must be internally
consistent.  SIF-off products intentionally return an empty dictionary so
older no-SIF products that predate this metadata remain valid.
"""
function validated_sif_provenance(attributes, sif_case;
                                  source::AbstractString="input")
    case = Symbol(sif_case)
    case == :off && return Dict{String,Any}()
    case == REQUIRED_SIF_CASE_ON || error(
        "$source uses stale or unsupported SIF case '$case'; expected " *
        "'$REQUIRED_SIF_CASE_ON'")

    provenance = Dict{String,Any}(
        key => _required_attribute(attributes, key, source)
        for key in SIF_PROVENANCE_ATTRIBUTES)
    Int(provenance["sif_definition_version"]) ==
        REQUIRED_SIF_DEFINITION_VERSION || error(
        "$source does not use sif_definition_version=" *
        "$REQUIRED_SIF_DEFINITION_VERSION")
    String(provenance["sif_definition"]) == REQUIRED_SIF_DEFINITION || error(
        "$source has an unexpected sif_definition for version " *
        "$REQUIRED_SIF_DEFINITION_VERSION")
    String(provenance["sif_case_on_label"]) ==
        String(REQUIRED_SIF_CASE_ON) || error(
        "$source has a SIF case-label/provenance mismatch")

    reference = _require_close(
        provenance["sif_reference_wavelength_nm"],
        REQUIRED_SIF_REFERENCE_WAVELENGTH_NM,
        "sif_reference_wavelength_nm", source)
    solid_angle = _require_close(
        provenance["sif_upwelling_solid_angle_sr"],
        REQUIRED_SIF_UPWELLING_SOLID_ANGLE_SR,
        "sif_upwelling_solid_angle_sr", source)
    angular_integral = _require_close(
        provenance["sif_angular_integral_760_mW_m-2_nm-1"],
        REQUIRED_SIF_ANGULAR_INTEGRAL_760,
        "sif_angular_integral_760_mW_m-2_nm-1", source)
    radiance = _require_close(
        provenance["sif_radiance_760_mW_m-2_sr-1_nm-1"],
        REQUIRED_SIF_RADIANCE_760,
        "sif_radiance_760_mW_m-2_sr-1_nm-1", source)
    cosine_irradiance = Float64(
        provenance["sif_cosine_weighted_irradiance_760_mW_m-2_nm-1"])
    isapprox(cosine_irradiance, π * radiance; atol=1e-14, rtol=1e-12) ||
        error("$source has inconsistent cosine-weighted SIF irradiance")
    isapprox(angular_integral, solid_angle * radiance;
             atol=1e-14, rtol=1e-12) ||
        error("$source has inconsistent SIF angular integral and radiance")

    native_sif760 = Float64(
        provenance["sif_SIF760_mW_m-2_sr-1_per_cm-1"])
    isapprox(native_sif760, radiance * reference^2 / 1e7;
             atol=1e-14, rtol=1e-12) ||
        error("$source has inconsistent wavelength/wavenumber SIF760 units")
    for key in ("sif_mSIF_mW_m-2_sr-1_per_cm-2",
                "sif_template_wavelength_integral_mW_m-2_sr-1")
        isfinite(Float64(provenance[key])) || error(
            "$source has a non-finite corrected-SIF provenance value $key")
    end
    return provenance
end

"""Require two versioned SIF records to describe the identical truth source."""
function matching_sif_provenance(first, second;
                                 first_source::AbstractString="first input",
                                 second_source::AbstractString="second input")
    Set(keys(first)) == Set(keys(second)) || error(
        "$first_source and $second_source carry different SIF provenance fields")
    for key in keys(first)
        left, right = first[key], second[key]
        matches = left isa Real && right isa Real ?
            isapprox(Float64(left), Float64(right); atol=1e-14, rtol=1e-12) :
            string(left) == string(right)
        matches || error(
            "SIF provenance attribute $key differs between " *
            "$first_source and $second_source")
    end
    return Dict{String,Any}(first)
end

function _truth_rows(path)
    isfile(path) || throw(ArgumentError("missing truth-state table: $path"))
    names = String[]
    rows = Vector{Dict{String,String}}()
    for line in eachline(path)
        stripped = strip(line)
        isempty(stripped) && continue
        if startswith(stripped, "# index ")
            names = split(stripped[3:end])
        elseif !startswith(stripped, '#')
            isempty(names) && error("truth-state column header was not found in $path")
            values = split(stripped)
            length(values) == length(names) || error(
                "truth-state row has $(length(values)) fields; expected $(length(names))")
            push!(rows, Dict(zip(names, values)))
        end
    end
    isempty(rows) && error("truth-state table contains no data rows: $path")
    return rows
end

function _required_float(row, name)
    haskey(row, name) || error("truth-state table is missing required column $name")
    value = parse(Float64, row[name])
    isfinite(value) || error("truth-state column $name contains a non-finite value")
    return value
end

function _required_int(row, name)
    haskey(row, name) || error("truth-state table is missing required column $name")
    return parse(Int, row[name])
end

"""Infer the uniform or bottom-layer CO2 profile metadata for one table row."""
function _co2_metadata(row)
    is_bottom_layer = haskey(row, "bottom_co2_ppm") ||
        haskey(row, "bottom_co2_index") || haskey(row, "background_co2_ppm")
    if is_bottom_layer
        bottom_index = _required_int(row, "bottom_co2_index")
        background = _required_float(row, "background_co2_ppm")
        bottom_layer = _required_int(row, "bottom_layer_index")
        bottom = _required_float(row, "bottom_co2_ppm")
        campaign = Symbol(get(row, "campaign", "bottom_layer_XCO2"))
        return (; xco2_index=bottom_index, campaign,
                co2_profile_mode=:bottom_layer,
                fixed_upper_co2_ppm=background,
                background_co2_ppm=background,
                bottom_layer_index=bottom_layer,
                bottom_co2_ppm=bottom)
    end

    xco2_index = _required_int(row, "xco2_index")
    xco2 = _required_float(row, "xco2_ppm")
    campaign = Symbol(get(row, "campaign", "full_column_XCO2"))
    return (; xco2_index, campaign, co2_profile_mode=:uniform_column,
            fixed_upper_co2_ppm=xco2, background_co2_ppm=xco2,
            bottom_layer_index=0, bottom_co2_ppm=xco2)
end

function _validate_truth_cases(cases, path)
    length(unique(case.state_index for case in cases)) == length(cases) || error(
        "truth-state indices are not unique in $path")
    sort([case.state_index for case in cases]) == collect(1:length(cases)) || error(
        "truth-state indices in $path must be contiguous from 1")
    length(unique(case.surface for case in cases)) == 4 || error(
        "truth table does not contain four surfaces")
    length(unique(case.aerosol_index for case in cases)) == 2 || error(
        "truth table does not contain aerosol/no-aerosol cases")
    length(unique(case.sif_case for case in cases)) == 2 || error(
        "truth table does not contain two SIF cases")

    modes = unique(case.co2_profile_mode for case in cases)
    length(modes) == 1 || error(
        "truth table mixes incompatible CO2 profile modes: $(join(modes, ", "))")
    mode = only(modes)
    expected_co2_cases = mode == :uniform_column ? 4 :
        mode == :bottom_layer ? 5 : error("unsupported CO2 profile mode $mode")
    length(unique(case.xco2_index for case in cases)) == expected_co2_cases || error(
        "$mode truth table does not contain $expected_co2_cases CO2 cases")
    expected_count = 4 * 2 * 2 * expected_co2_cases
    length(cases) == expected_count || error(
        "expected $expected_count $mode truth cases; found $(length(cases))")
    combinations = Set((case.surface_index, case.aerosol_index,
                        case.sif_case, case.xco2_index) for case in cases)
    length(combinations) == expected_count || error(
        "truth table does not contain a complete surface/aerosol/SIF/CO2 grid")

    if mode == :bottom_layer
        length(unique(case.background_co2_ppm for case in cases)) == 1 || error(
            "bottom-layer truth table has inconsistent background CO2")
        length(unique(case.bottom_layer_index for case in cases)) == 1 || error(
            "bottom-layer truth table has inconsistent perturbed-layer indices")
        all(case -> case.fixed_upper_co2_ppm == case.background_co2_ppm,
            cases) || error("bottom-layer fixed-upper and background CO2 differ")
        length(unique(case.bottom_co2_ppm for case in cases)) == 5 || error(
            "bottom-layer truth table does not contain five bottom-layer VMRs")
    else
        all(case -> case.fixed_upper_co2_ppm == case.xco2_ppm &&
                    case.background_co2_ppm == case.xco2_ppm &&
                    case.bottom_co2_ppm == case.xco2_ppm, cases) || error(
            "uniform truth metadata is internally inconsistent")
    end
    return cases
end

"""
Read a complete full-column (64-state) or bottom-layer (80-state) truth table.

The original defaults and `xco2_index` field remain unchanged for the
full-column campaign. A bottom-layer table supplies `bottom_co2_index`
instead; it is exposed through the same generic case-index field and carries
the 400 ppm background/fixed-upper value separately from column XCO2.
"""
function read_truth_cases(path::AbstractString=default_truth_table())
    cases = TruthCase[]
    for row in _truth_rows(path)
        co2 = _co2_metadata(row)
        push!(cases, TruthCase(
            parse(Int, row["index"]),
            parse(Int, row["surface_index"]),
            Symbol(row["surface"]),
            parse(Int, row["aerosol_index"]),
            Symbol(row["aerosol_case"]),
            Symbol(row["sif_case"]),
            co2.xco2_index,
            _required_float(row, "xco2_ppm"),
            co2.campaign,
            co2.co2_profile_mode,
            co2.fixed_upper_co2_ppm,
            co2.background_co2_ppm,
            co2.bottom_layer_index,
            co2.bottom_co2_ppm,
        ))
    end
    sort!(cases; by=case -> case.state_index)
    return _validate_truth_cases(cases, path)
end

"""
Select truth scenes by SIF state.

`sif_filter` accepts `off`, `on`, or `all` as either a string or symbol. `on`
means any named SIF case other than `off`, so the selector remains usable if
additional nonzero-SIF truth amplitudes are introduced later.
"""
function select_sif_truth_cases(cases, sif_filter)
    normalized = Symbol(lowercase(String(sif_filter)))
    normalized in (:off, :on, :all) || throw(ArgumentError(
        "SIF case filter must be off, on, or all"))
    normalized == :all && return collect(cases)
    want_sif = normalized == :on
    return filter(case -> (case.sif_case != :off) == want_sif, cases)
end

"""Read and validate the no-SIF half of either supported truth campaign."""
function read_no_sif_truth_cases(path::AbstractString=default_truth_table())
    cases = read_truth_cases(path)
    selected = select_sif_truth_cases(cases, :off)
    length(selected) == length(cases) ÷ 2 || error(
        "expected $(length(cases) ÷ 2) no-SIF truth cases; " *
        "found $(length(selected))")
    length(unique(case.surface for case in selected)) == 4 || error(
        "no-SIF truth subset does not contain four surfaces")
    length(unique(case.aerosol_index for case in selected)) == 2 || error(
        "no-SIF truth subset does not contain aerosol/no-aerosol cases")
    expected_co2_cases = only(unique(case.co2_profile_mode for case in selected)) ==
        :bottom_layer ? 5 : 4
    length(unique(case.xco2_index for case in selected)) == expected_co2_cases ||
        error("no-SIF truth subset does not contain $expected_co2_cases CO2 cases")
    return selected
end

@inline function _pair_seed(master_seed::UInt64, state_index::Int, noise_index::Int)
    # Fixed integer arithmetic avoids Julia's session-randomized hash and gives
    # the corrected/uncorrected pair exactly the same standardized noise draw.
    return master_seed + UInt64(10_000state_index + noise_index)
end

"""
Build the canonical experiment sequence with adjacent class pairs.

Indices 1:10 are independent unit-variance uniform noise draws. Index 11 is
the unperturbed experiment and has an exact zero injected-noise vector. When
it is included, index 11 is scheduled first within each truth state, followed
by indices 1:10; stored indices and random seeds are unchanged.
"""
function build_experiments(cases=read_no_sif_truth_cases();
                           noise_realizations::Int=PERTURBED_REALIZATIONS,
                           include_unperturbed::Bool=true,
                           master_seed::UInt64=DEFAULT_MASTER_SEED,
                           measurement_directory::AbstractString=
                               default_measurement_directory(),
                           noise_directory::AbstractString=
                               default_noise_directory(),
                           validate_inputs::Bool=true)
    noise_realizations == PERTURBED_REALIZATIONS || throw(ArgumentError(
        "production suite requires exactly $PERTURBED_REALIZATIONS perturbed realizations"))
    perturbation_indices = collect(1:noise_realizations)
    include_unperturbed && pushfirst!(perturbation_indices, UNPERTURBED_INDEX)
    experiments = RetrievalExperiment[]
    pairs_per_case = length(perturbation_indices)
    classes_per_pair = length(MEASUREMENT_CLASSES)
    for (case_position, truth) in enumerate(cases),
            noise_index in perturbation_indices
        # IDs retain their canonical numeric-perturbation mapping even though
        # execution begins with perturbation 11. This keeps existing products
        # and manifests comparable across the scheduling change.
        pair_index = (case_position - 1) * pairs_per_case + noise_index
        # The seed is deliberately zero for the unperturbed member so its
        # manifest cannot be misread as identifying a random realization.
        seed = noise_index == UNPERTURBED_INDEX ? UInt64(0) :
            _pair_seed(master_seed, truth.state_index, noise_index)
        measurement_path = joinpath(
            measurement_directory, @sprintf("OCO2sims_%03d.nc", truth.state_index))
        noise_path = joinpath(
            noise_directory, @sprintf("OCO2noise_%03d.nc", truth.state_index))
        if validate_inputs
            isfile(measurement_path) || error("missing measurement file: $measurement_path")
            isfile(noise_path) || error("missing noise file: $noise_path")
        end
        for (class_position, measurement_class) in enumerate(MEASUREMENT_CLASSES)
            retrieval_index = (pair_index - 1) * classes_per_pair + class_position
            push!(experiments, RetrievalExperiment(
                retrieval_index, pair_index, truth, noise_index,
                measurement_class, seed, measurement_path, noise_path))
        end
    end
    expected = length(cases) * length(perturbation_indices) *
        classes_per_pair
    length(experiments) == expected || error("experiment-count construction failed")
    return experiments
end

"""Generate a reproducible unit-variance uniform draw in [-sqrt(3),sqrt(3)]."""
function paired_uniform_draw(seed::UInt64, count::Integer)
    count > 0 || throw(ArgumentError("draw length must be positive"))
    rng = Xoshiro(seed)
    return (2sqrt(3.0)) .* rand(rng, count) .- sqrt(3.0)
end

function _finite_vector(dataset, name)
    haskey(dataset, name) || error("NetCDF input is missing $name")
    values = Float64.(nomissing(dataset[name][:], NaN))
    all(isfinite, values) || error("$name contains non-finite or missing values")
    return values
end

"""Load one frozen `S_epsilon` product and apply its deterministic pair draw."""
function load_measurement_realization(experiment::RetrievalExperiment)
    label = String(experiment.measurement_class)
    truth = experiment.truth
    measurement_provenance = NCDataset(experiment.measurement_path) do dataset
        get(dataset.attrib, "instrument_processing_complete", 0) == 1 || error(
            "measurement file is not marked complete: " *
            experiment.measurement_path)
        Int(dataset.attrib["state_index"]) == truth.state_index || error(
            "measurement-file state metadata disagrees with experiment")
        measurement_sif = Symbol(get(dataset.attrib, "sif_case", "off"))
        measurement_sif == truth.sif_case || error(
            "measurement-file SIF case disagrees with experiment")
        validated_sif_provenance(
            dataset.attrib, measurement_sif; source=experiment.measurement_path)
    end
    return NCDataset(experiment.noise_path) do dataset
        get(dataset.attrib, "noise_covariance_complete", 0) == 1 || error(
            "noise file is not marked complete: $(experiment.noise_path)")
        Int(dataset.attrib["state_index"]) == experiment.truth.state_index || error(
            "noise-file state metadata disagrees with experiment")
        noise_sif = Symbol(get(dataset.attrib, "sif_case", "off"))
        noise_sif == truth.sif_case || error(
            "noise-file SIF case disagrees with experiment")
        noise_provenance = validated_sif_provenance(
            dataset.attrib, noise_sif; source=experiment.noise_path)
        provenance = matching_sif_provenance(
            measurement_provenance, noise_provenance;
            first_source=experiment.measurement_path,
            second_source=experiment.noise_path)
        if truth.co2_profile_mode == :bottom_layer
            get(dataset.attrib, "campaign", "") == String(truth.campaign) || error(
                "noise-file campaign metadata disagrees with bottom-layer experiment")
            Float64(get(dataset.attrib, "background_co2_ppm", NaN)) ==
                truth.background_co2_ppm || error(
                "noise-file background CO2 disagrees with bottom-layer experiment")
            Int(get(dataset.attrib, "bottom_co2_layer_index", -1)) ==
                truth.bottom_layer_index || error(
                "noise-file bottom-layer index disagrees with experiment")
            Float64(get(dataset.attrib, "bottom_co2_ppm", NaN)) ==
                truth.bottom_co2_ppm || error(
                "noise-file bottom-layer VMR disagrees with experiment")
            isapprox(Float64(get(dataset.attrib, "xco2_ppm", NaN)),
                     truth.xco2_ppm; atol=1e-12, rtol=0) || error(
                "noise-file column XCO2 disagrees with bottom-layer experiment")
        end
        noiseless = _finite_vector(dataset, "measurement_$label")
        noise_std = _finite_vector(dataset, "noise_std_$label")
        variance = _finite_vector(dataset, "Se_diagonal_$label")
        variance ≈ noise_std .^ 2 || error(
            "stored variance is inconsistent with noise standard deviation")
        all(>(0), variance) || error("measurement variances must be positive")
        wavelength = _finite_vector(dataset, "wavelength")
        starts = Int.(dataset["band_start_index"][:])
        stops = Int.(dataset["band_end_index"][:])
        ranges = UnitRange{Int}[first:last for (first, last) in zip(starts, stops)]
        reduce(vcat, collect.(ranges)) == collect(eachindex(noiseless)) || error(
            "stored band ranges do not partition the measurement vector")
        draw = experiment.noise_index == UNPERTURBED_INDEX ?
            zeros(Float64, length(noiseless)) :
            paired_uniform_draw(experiment.random_seed, length(noiseless))
        perturbed = noiseless + noise_std .* draw
        return MeasurementRealization(
            noiseless, perturbed, noise_std, variance, draw, wavelength, ranges,
            provenance)
    end
end

"""Write a human-readable, whitespace-separated catalog of all solves."""
function write_experiment_manifest(
        experiments;
        output_path::Union{Nothing,AbstractString}=nothing,
        inversion_root::AbstractString=INVERSION_ROOT)
    output_path = isnothing(output_path) ?
        joinpath(inversion_root, "retrieval_manifest.dat") : String(output_path)
    mkpath(dirname(output_path))
    open(output_path, "w") do io
        println(io, "# RRS-XCO2 retrieval suite: paired corrected/uncorrected solves")
        println(io, "# Perturbations 01:10: noise_i=sigma_i*u_i, " *
                    "u_i uniform on [-sqrt(3),sqrt(3)]; the same u is used within each pair.")
        println(io, "# Perturbation 11: unperturbed measurement; u_i=noise_i=0 exactly.")
        println(io, "# Execution order within each state: 11, then 01:10; canonical IDs are unchanged.")
        println(io, "# retrieval pair state surface aerosol sif campaign co2_profile " *
                    "co2_case xco2_ppm fixed_upper_co2_ppm background_co2_ppm " *
                    "bottom_layer bottom_co2_ppm perturbation noise_injected inputs_ready class " *
                    "seed output_file measurement_file noise_file")
        for experiment in experiments
            truth = experiment.truth
            output_file = joinpath(
                inversion_root, String(experiment.measurement_class),
                @sprintf("retrieval_state%03d_perturbation%02d.nc",
                         truth.state_index, experiment.noise_index))
            noise_injected = experiment.noise_index != UNPERTURBED_INDEX
            inputs_ready = isfile(experiment.measurement_path) &&
                isfile(experiment.noise_path)
            @printf(io, "%04d %03d %03d %-7s %-12s %-10s %-18s %-12s %02d %.12f %.12f %.12f %02d %.12f %02d %d %d %-11s %020u %s %s %s\n",
                    experiment.retrieval_index,
                    experiment.pair_index,
                    truth.state_index,
                    String(truth.surface),
                    String(truth.aerosol_case),
                    String(truth.sif_case),
                    String(truth.campaign),
                    String(truth.co2_profile_mode),
                    truth.xco2_index,
                    truth.xco2_ppm,
                    truth.fixed_upper_co2_ppm,
                    truth.background_co2_ppm,
                    truth.bottom_layer_index,
                    truth.bottom_co2_ppm,
                    experiment.noise_index,
                    noise_injected,
                    inputs_ready,
                    String(experiment.measurement_class),
                    experiment.random_seed,
                    output_file,
                    experiment.measurement_path,
                    experiment.noise_path)
        end
    end
    return output_path
end

end # module RetrievalCases
