#!/usr/bin/env julia

"""
Build one disjoint partition of the 80-scene bottom-layer CO2 truth campaign.

This is deliberately a CO2-only producer.  The accepted uniform-400-ppm truth
map already contains the exact O2 A-band result needed for every bottom-layer
CO2 case.  The `nosif` stage copies that O2 result and recomputes only weak and
strong CO2 for 360, 380, 420, and 440 ppm in layer 16.  The 400-ppm control is
an exact component copy.  The later `sif` stage combines the accepted SIF-on
O2 result with the already completed no-SIF CO2 bands and performs no RT.

Required environment controls:

- `AEROSOL_CASE_FILTER=none|aerosol` selects the Curry or Wurst partition.
- `BOTTOM_TRUTH_STAGE=nosif|sif` enforces the agreed no-SIF-first sequence.

Optional controls:

- `BOTTOM_XCO2_CAMPAIGN_ROOT` (default
  `RRS_XCO2/bottom_layer_XCO2_retrievals`)
- `FULL_COLUMN_TRUTH_ROOT` (default `RRS_XCO2/truth_map`)
- `CO2_CHUNK_POINTS` (default 2048; both current bands fit in one chunk)
- `CUDA_DEVICE` (used only by the `nosif` stage). The partition launcher
  selects a physical GPU with `BOTTOM_TRUTH_PHYSICAL_GPU`, masks all other
  devices, and sets this logical index to 0.
- `BOTTOM_TRUTH_PLAN_ONLY=1` validates and prints the plan without writing.

Final files are never opened in place.  Every scene is protected by an atomic
claim directory, assembled in a unique staging directory on the same file
system, validated, and then published by rename.  A stale claim is never
deleted automatically; its owner file identifies the interrupted producer.
"""

include(joinpath(@__DIR__, "bottom_layer_truth_common.jl"))
import .BottomLayerTruthCommon

const BL_FILTER = lowercase(get(ENV, "AEROSOL_CASE_FILTER", ""))
BL_FILTER in ("none", "aerosol") || error(
    "AEROSOL_CASE_FILTER must be explicitly set to none or aerosol")

const BL_STAGE = lowercase(get(ENV, "BOTTOM_TRUTH_STAGE", ""))
BL_STAGE in ("nosif", "sif") || error(
    "BOTTOM_TRUTH_STAGE must be explicitly set to nosif or sif")

const BL_CAMPAIGN_ROOT = BottomLayerTruthCommon.CAMPAIGN_ROOT
const BL_TRUTH_ROOT = BottomLayerTruthCommon.TRUTH_ROOT
const BL_OUT = BL_FILTER == "aerosol" ?
    joinpath(BL_TRUTH_ROOT, "aerosol_chunked") : BL_TRUTH_ROOT
ENV["TRUTH_OUT"] = BL_OUT

# Reuse the exact production surface, aerosol, solar, ABSCO, and CO2-band
# forward helpers.  `prepare_truth_co2_profile!` is specialized below for the
# new state type; the existing 64-scene producer retains its uniform method.
include(joinpath(@__DIR__, "generate_truth_map.jl"))

using Dates
using NCDatasets
using Printf
using SHA
using Sockets
using UUIDs

const BLState = BottomLayerTruthCommon.BottomTruthState
const FULL_COLUMN_TRUTH_ROOT = BottomLayerTruthCommon.FULL_COLUMN_TRUTH_ROOT
const FULL_COLUMN_STATE_TABLE = BottomLayerTruthCommon.FULL_COLUMN_STATE_FILE
const CO2_CHUNK_POINTS_BOTTOM = parse(Int, get(ENV, "CO2_CHUNK_POINTS", "2048"))
const PLAN_ONLY = lowercase(get(ENV, "BOTTOM_TRUTH_PLAN_ONLY", "0")) in
                  ("1", "true", "yes", "on")
const CLAIMS_ROOT = joinpath(BL_TRUTH_ROOT, ".claims")
const STAGING_ROOT = joinpath(BL_OUT, ".staging")
const STATE_TABLE = BottomLayerTruthCommon.STATE_FILE
const REQUIRED_VARIABLES = (
    "radiance_rayleigh_o2a",
    "radiance_cabannes_o2a",
    "radiance_rrs_o2a",
    "radiance_rayleigh_weak_co2",
    "radiance_rayleigh_strong_co2",
)
const O2_VARIABLES = REQUIRED_VARIABLES[1:3]
const CO2_VARIABLES = REQUIRED_VARIABLES[4:5]
const VALIDATED_SOURCE_GRIDS = Set{String}()
const FULL_COLUMN_STATE_ROWS = Ref{Union{Nothing,Dict{Int,Vector{Any}}}}(nothing)

CO2_CHUNK_POINTS_BOTTOM > 0 || error("CO2_CHUNK_POINTS must be positive")
NLAYERS == 16 || error("bottom-layer truth requires NLAYERS=16")
FT === Float32 || error("bottom-layer production is fixed to Float32")
SZA_DEG == FT(30) && VZA_DEG == zero(FT) &&
    RELATIVE_AZIMUTH_DEG == zero(FT) || error(
        "bottom-layer truth requires SZA=30, VZA=0, relative azimuth=0")
AEROSOL_NSTREAMS == 9 || error(
    "bottom-layer aerosol truth requires AEROSOL_NSTREAMS=9")

"Materialize the common grid first, then change only TOA-to-BOA layer 16."
function prepare_truth_co2_profile!(params, state::BLState)
    state.bottom_layer_index == NLAYERS || error(
        "bottom-layer index $(state.bottom_layer_index) does not match $NLAYERS layers")
    params.absorption_params.vmr["CO2"] =
        FT(state.background_co2_ppm * 1e-6)
    RRSXCO2Common.prepare_shared_profile!(
        params; psurf=state.psurf_hpa, nlayers=NLAYERS)
    co2 = fill(FT(state.background_co2_ppm * 1e-6), NLAYERS)
    co2[state.bottom_layer_index] = FT(state.bottom_co2_ppm * 1e-6)
    params.absorption_params.vmr["CO2"] = co2
    return params
end

scene_path(state::BLState) =
    joinpath(BL_OUT, @sprintf("hiressim_%03d.nc", state.index))

function source_root(state::BLState)
    return state.aerosol_index == 2 ?
        joinpath(FULL_COLUMN_TRUTH_ROOT, "aerosol_chunked") :
        FULL_COLUMN_TRUTH_ROOT
end

source_scene_path(state::BLState) = joinpath(
    source_root(state),
    @sprintf("hiressim_%03d.nc",
             BottomLayerTruthCommon.old_control_index(state)))

function no_sif_scene_path(state::BLState)
    index = BottomLayerTruthCommon.paired_no_sif_index(state)
    return joinpath(BL_OUT, @sprintf("hiressim_%03d.nc", index))
end

function campaign_scene_path(state::BLState)
    root = state.aerosol_index == 2 ?
        joinpath(BL_TRUTH_ROOT, "aerosol_chunked") : BL_TRUTH_ROOT
    return joinpath(root, @sprintf("hiressim_%03d.nc", state.index))
end

function state_at_index(index::Integer)
    states = BottomLayerTruthCommon.read_bottom_states(STATE_TABLE)
    1 <= index <= length(states) || throw(BoundsError(states, index))
    states[index].index == index || error("state table is not index ordered")
    return states[index]
end

_attribute(ds, name, default=nothing) =
    haskey(ds.attrib, name) ? ds.attrib[name] : default

file_sha256(path::AbstractString) = bytes2hex(sha256(read(path)))

function _relative_to_experiment(path::AbstractString)
    return relpath(abspath(path), abspath(BottomLayerTruthCommon.ROOT))
end

"""Validate corrected SIF provenance while accepting untouched legacy off files."""
function validate_scene_sif_provenance(dataset, state::BLState, path)
    if state.sif_index == 1
        state.sif_case == "off" || error("invalid SIF-off state metadata: $path")
        # No-SIF radiances are numerically independent of the normalization
        # convention. Existing accepted files may predate versioned SIF attrs;
        # if those attrs are present, however, they must describe zero SIF.
        if haskey(dataset.attrib, "sif_definition_version")
            Int(dataset.attrib["sif_definition_version"]) ==
                RRSXCO2Common.SIF_DEFINITION_VERSION || error(
                    "unsupported SIF definition version on no-SIF scene: $path")
            iszero(Float64(_attribute(dataset,
                "sif_angular_integral_760_mW_m-2_nm-1", NaN))) || error(
                    "no-SIF scene carries a nonzero SIF angular integral: $path")
        end
        return nothing
    end

    Int(_attribute(dataset, "sif_definition_version", 0)) ==
        RRSXCO2Common.SIF_DEFINITION_VERSION || error(
        "SIF-on scene predates the corrected 760-nm SIF definition: $path")
    isapprox(Float64(_attribute(
                 dataset, "sif_angular_integral_760_mW_m-2_nm-1", NaN)),
             state.sif_angular_integral760; atol=1e-12, rtol=0) || error(
        "SIF angular integral differs: $path")
    return nothing
end

function expected_sizes(grids)
    return Dict(
        "radiance_rayleigh_o2a" => (3, length(grids[1])),
        "radiance_cabannes_o2a" => (3, length(grids[1])),
        "radiance_rrs_o2a" => (3, length(grids[1])),
        "radiance_rayleigh_weak_co2" => (3, length(grids[2])),
        "radiance_rayleigh_strong_co2" => (3, length(grids[3])),
    )
end

function full_column_state_rows()
    cached = FULL_COLUMN_STATE_ROWS[]
    cached === nothing || return cached
    isfile(FULL_COLUMN_STATE_TABLE) || error(
        "missing accepted full-column state table: $FULL_COLUMN_STATE_TABLE")
    raw = readdlm(FULL_COLUMN_STATE_TABLE; comments=true)
    size(raw) == (64, 27) || error(
        "expected a 64x27 full-column state table, found $(size(raw))")
    rows = Dict{Int,Vector{Any}}()
    for row in eachrow(raw)
        values = Any[item for item in row]
        index = Int(values[1])
        haskey(rows, index) && error(
            "duplicate state $index in $FULL_COLUMN_STATE_TABLE")
        rows[index] = values
    end
    length(rows) == 64 || error("full-column state table is not complete")
    FULL_COLUMN_STATE_ROWS[] = rows
    return rows
end

function validate_full_column_state_record(state::BLState)
    old_index = BottomLayerTruthCommon.old_control_index(state)
    row = get(full_column_state_rows(), old_index, nothing)
    row === nothing && error("full-column state table lacks state $old_index")
    Int(row[2]) == state.surface_index && String(row[3]) == state.surface ||
        error("full-column state-table surface differs for state $old_index")
    Int(row[4]) == state.aerosol_index &&
        String(row[5]) == state.aerosol_case || error(
            "full-column state-table aerosol differs for state $old_index")
    Int(row[6]) == state.sif_index && String(row[7]) == state.sif_case ||
        error("full-column state-table SIF case differs for state $old_index")
    Int(row[9]) == 400 || error(
        "full-column reuse source $old_index is not the 400-ppm control")
    Float64(row[10]) == state.psurf_hpa &&
        Float64(row[11]) == state.sza_deg &&
        Float64(row[12]) == state.vza_deg || error(
            "full-column state-table pressure/geometry differs for state $old_index")
    old_aod = FT.(Float64.(row[13:15]))
    old_aod == FT.(collect(state.aod550)) || error(
        "full-column state-table AOD differs for state $old_index")
    FT(Float64(row[16])) == FT(state.sif_angular_integral760) || error(
        "full-column state-table SIF angular integral differs for state $old_index")
    FT(Float64(row[17])) == FT(state.sif760) &&
        FT(Float64(row[18])) == FT(state.msif) || error(
            "full-column state-table SIF shape differs for state $old_index")
    old_coefficients = FT.(Float64.(row[19:27]))
    new_coefficients = FT.(collect(Iterators.flatten(state.coeff)))
    old_coefficients == new_coefficients || error(
        "full-column state-table surface coefficients differ for state $old_index")
    return row
end

function validate_source_grid(root::AbstractString, grids)
    path = joinpath(root, "sim_wavelength.nc")
    path in VALIDATED_SOURCE_GRIDS && return path
    isfile(path) || error("missing accepted wavelength file: $path")
    NCDataset(path, "r") do ds
        for (name, grid) in zip(
                ("o2a", "weak_co2", "strong_co2"), grids)
            wn_name = "$(name)_wavenumber"
            wl_name = "$(name)_wavelength"
            haskey(ds, wn_name) && haskey(ds, wl_name) || error(
                "accepted wavelength file lacks $name coordinates: $path")
            stored_wn = Float64.(ds[wn_name][:])
            stored_wl = Float64.(ds[wl_name][:])
            stored_wn == Float64.(grid) || error(
                "$name wavenumbers differ from the canonical production grid: $path")
            stored_wl == Float64.(1e7 ./ grid) || error(
                "$name wavelengths differ from the canonical production grid: $path")
        end
    end
    push!(VALIDATED_SOURCE_GRIDS, path)
    return path
end

function validate_source_scene(state::BLState, path::AbstractString, grids)
    isfile(path) || error("missing accepted full-column source: $path")
    validate_full_column_state_record(state)
    validate_source_grid(dirname(path), grids)
    NCDataset(path, "r") do ds
        Int(_attribute(ds, "simulation_complete", 0)) == 1 || error(
            "source is not marked complete: $path")
        Int(_attribute(ds, "state_index", -1)) ==
            BottomLayerTruthCommon.old_control_index(state) || error(
                "source state-index mismatch: $path")
        String(_attribute(ds, "surface", "")) == state.surface || error(
            "source surface mismatch: $path")
        String(_attribute(ds, "aerosol_case", "")) == state.aerosol_case ||
            error("source aerosol mismatch: $path")
        String(_attribute(ds, "sif_case", "")) == state.sif_case || error(
            "source SIF mismatch: $path")
        Float64(_attribute(ds, "xco2_ppm", NaN)) == 400.0 || error(
            "source is not the uniform-400 control: $path")
        Int(_attribute(ds, "atmospheric_layers", -1)) == 16 || error(
            "source does not use 16 layers: $path")
        Float64(_attribute(ds, "psurf_hpa", NaN)) == state.psurf_hpa ||
            error("source surface pressure differs: $path")
        Float64(_attribute(ds, "sza_deg", NaN)) == state.sza_deg || error(
            "source solar zenith angle differs: $path")
        Float64(_attribute(ds, "vza_deg", NaN)) == state.vza_deg || error(
            "source viewing zenith angle differs: $path")
        Float64(_attribute(ds, "relative_azimuth_deg", 0.0)) ==
            state.relative_azimuth_deg || error(
                "source relative azimuth differs: $path")
        source_aod = Float64[
            _attribute(ds, "aod550_sulfate", NaN),
            _attribute(ds, "aod550_organic_carbon", NaN),
            _attribute(ds, "aod550_stratospheric", NaN),
        ]
        all(isapprox.(source_aod, collect(state.aod550);
                      atol=5e-7, rtol=0)) || error(
            "source aerosol optical depths differ: $path")
        validate_scene_sif_provenance(ds, state, path)
        isapprox(Float64(_attribute(ds, "o2_vmr", NaN)), 0.21;
                 atol=1e-12, rtol=0) || error("source O2 VMR differs: $path")
        String(_attribute(ds, "spectroscopy_database", "")) == "ABSCO" ||
            error("source does not record ABSCO spectroscopy: $path")
        String(_attribute(ds, "spectroscopy_version", "")) ==
            RRSXCO2Common.ABSCO_VERSION || error(
                "source ABSCO version differs: $path")
        sizes = expected_sizes(grids)
        for variable in REQUIRED_VARIABLES
            haskey(ds, variable) || error("source lacks $variable: $path")
            size(ds[variable]) == sizes[variable] || error(
                "$variable has size $(size(ds[variable])) in $path; expected $(sizes[variable])")
            all(isfinite, ds[variable][:, :]) || error(
                "$variable contains non-finite values in $path")
        end
    end
    return path
end

function selected_states()
    states = BottomLayerTruthCommon.read_bottom_states(STATE_TABLE)
    aerosol_index = BL_FILTER == "aerosol" ? 2 : 1
    sif_index = BL_STAGE == "sif" ? 2 : 1
    selected = filter(
        state -> state.aerosol_index == aerosol_index &&
                 state.sif_index == sif_index,
        states)
    length(selected) == 20 || error(
        "expected 20 states for $BL_FILTER/$BL_STAGE, found $(length(selected))")
    return selected
end

function ensure_wavelength_file!()
    source = joinpath(BL_FILTER == "aerosol" ?
        joinpath(FULL_COLUMN_TRUTH_ROOT, "aerosol_chunked") :
        FULL_COLUMN_TRUTH_ROOT, "sim_wavelength.nc")
    destination = joinpath(BL_OUT, "sim_wavelength.nc")
    isfile(source) || error("missing accepted wavelength file: $source")
    if isfile(destination)
        file_sha256(destination) == file_sha256(source) || error(
            "existing wavelength file differs from accepted source: $destination")
        return destination
    end
    mkpath(BL_OUT)
    temporary = destination * ".tmp.$(getpid()).$(uuid4())"
    cp(source, temporary)
    mv(temporary, destination)
    return destination
end

function _profile_ppm(state::BLState)
    profile = fill(state.background_co2_ppm, NLAYERS)
    profile[state.bottom_layer_index] = state.bottom_co2_ppm
    return profile
end

function mark_staging!(path::AbstractString, state::BLState,
                       source::AbstractString, mode::AbstractString;
                       co2_source::Union{Nothing,AbstractString}=nothing)
    NCDataset(path, "a") do ds
        ds.attrib["campaign"] = "bottom_layer_XCO2"
        ds.attrib["state_index"] = Int32(state.index)
        ds.attrib["surface"] = state.surface
        ds.attrib["aerosol_case"] = state.aerosol_case
        ds.attrib["sif_case"] = state.sif_case
        RRSXCO2Common.write_sif_provenance!(
            ds.attrib, state.sif_angular_integral760 > 0)
        ds.attrib["xco2_ppm"] = state.xco2_ppm
        ds.attrib["column_xco2_ppm"] = state.xco2_ppm
        ds.attrib["background_co2_ppm"] = state.background_co2_ppm
        ds.attrib["bottom_co2_ppm"] = state.bottom_co2_ppm
        ds.attrib["bottom_co2_index"] = Int32(state.bottom_co2_index)
        ds.attrib["bottom_co2_layer_index"] = Int32(state.bottom_layer_index)
        ds.attrib["co2_profile_order"] = "TOA-to-BOA"
        ds.attrib["co2_profile_definition"] =
            "layers 1:15 fixed at background_co2_ppm; layer 16 set to bottom_co2_ppm"
        ds.attrib["co2_profile_ppm"] = _profile_ppm(state)
        ds.attrib["bottom_layer_dry_column_fraction"] =
            BottomLayerTruthCommon.BOTTOM_DRY_COLUMN_FRACTION
        ds.attrib["source_state_table"] =
            BL_FILTER == "aerosol" ? "../true_states.dat" : "true_states.dat"
        ds.attrib["state_table_sha256"] =
            BottomLayerTruthCommon.state_table_sha256(STATE_TABLE)
        ds.attrib["full_column_source_state"] = Int32(
            BottomLayerTruthCommon.old_control_index(state))
        ds.attrib["full_column_source_scene"] = _relative_to_experiment(source)
        ds.attrib["full_column_source_sha256"] = file_sha256(source)
        ds.attrib["full_column_state_table_sha256"] =
            file_sha256(FULL_COLUMN_STATE_TABLE)
        ds.attrib["o2_truth_reused"] = Int32(1)
        ds.attrib["o2_truth_reuse_source"] = _relative_to_experiment(source)
        ds.attrib["bottom_co2_truth_mode"] = mode
        if co2_source !== nothing
            ds.attrib["co2_truth_reuse_source"] =
                _relative_to_experiment(co2_source)
            ds.attrib["co2_truth_reuse_sha256"] = file_sha256(co2_source)
        else
            ds.attrib["co2_truth_reuse_source"] = "none"
        end
        ds.attrib["producer_script"] = _relative_to_experiment(@__FILE__)
        ds.attrib["producer_script_sha256"] = file_sha256(@__FILE__)
        ds.attrib["strong_co2_shoulder_source"] =
            "canonical 995-point strong-band grid carried in this scene"
        ds.attrib["co2_absco_regeneration_complete"] = Int32(0)
        ds.attrib["co2_absco_completed"] = "pending"
        ds.attrib["simulation_complete"] = Int32(0)
        ds.attrib["bottom_layer_truth_complete"] = Int32(0)
        if haskey(ds.attrib, "chunked_simulation_complete")
            ds.attrib["chunked_simulation_complete"] = Int32(0)
        end
        ds.attrib["created"] = string(now(UTC))
        ds.attrib["completed"] = "pending"
    end
    return path
end

function mark_complete!(path::AbstractString)
    NCDataset(path, "a") do ds
        ds.attrib["simulation_complete"] = Int32(1)
        ds.attrib["bottom_layer_truth_complete"] = Int32(1)
        if haskey(ds.attrib, "chunked_simulation_complete")
            ds.attrib["chunked_simulation_complete"] = Int32(1)
        end
        completed = string(now(UTC))
        ds.attrib["co2_absco_regeneration_complete"] = Int32(1)
        ds.attrib["co2_absco_completed"] = completed
        ds.attrib["strong_co2_shoulder_merged_utc"] = completed
        ds.attrib["completed"] = completed
    end
    return path
end

function _arrays_equal(path_a, path_b, variables)
    NCDataset(path_a, "r") do a
        NCDataset(path_b, "r") do b
            return all(variable -> a[variable][:, :] == b[variable][:, :],
                       variables)
        end
    end
end

function validate_bottom_scene(state::BLState, path::AbstractString, grids;
                               require_complete::Bool=true,
                               co2_source::Union{Nothing,AbstractString}=nothing)
    isfile(path) || error("missing bottom-layer truth scene: $path")
    NCDataset(path, "r") do ds
        require_complete &&
            Int(_attribute(ds, "bottom_layer_truth_complete", 0)) != 1 &&
            error("bottom-layer scene is not marked complete: $path")
        String(_attribute(ds, "campaign", "")) == "bottom_layer_XCO2" ||
            error("campaign metadata mismatch: $path")
        Int(_attribute(ds, "state_index", -1)) == state.index || error(
            "state-index metadata mismatch: $path")
        String(_attribute(ds, "surface", "")) == state.surface || error(
            "surface metadata mismatch: $path")
        String(_attribute(ds, "aerosol_case", "")) == state.aerosol_case ||
            error("aerosol metadata mismatch: $path")
        String(_attribute(ds, "sif_case", "")) == state.sif_case || error(
            "SIF metadata mismatch: $path")
        validate_scene_sif_provenance(ds, state, path)
        isapprox(Float64(_attribute(ds, "xco2_ppm", NaN)), state.xco2_ppm;
                 atol=5e-10, rtol=0) || error("XCO2 metadata mismatch: $path")
        Float64(_attribute(ds, "bottom_co2_ppm", NaN)) ==
            state.bottom_co2_ppm || error("bottom CO2 metadata mismatch: $path")
        profile = Float64.(_attribute(ds, "co2_profile_ppm", Float64[]))
        profile == _profile_ppm(state) || error("CO2 profile metadata mismatch: $path")
        sizes = expected_sizes(grids)
        for variable in REQUIRED_VARIABLES
            haskey(ds, variable) || error("scene lacks $variable: $path")
            size(ds[variable]) == sizes[variable] || error(
                "$variable has incorrect size in $path")
            all(isfinite, ds[variable][:, :]) || error(
                "$variable contains non-finite values in $path")
        end
    end

    source = source_scene_path(state)
    _arrays_equal(path, source, O2_VARIABLES) || error(
        "O2 arrays are not bit-identical to their accepted source: $path")
    if co2_source !== nothing
        _arrays_equal(path, co2_source, CO2_VARIABLES) || error(
            "SIF-on CO2 arrays differ from the paired no-SIF scene: $path")
    elseif state.bottom_co2_ppm == state.background_co2_ppm
        _arrays_equal(path, source, CO2_VARIABLES) || error(
            "400-ppm control CO2 arrays differ from their accepted source: $path")
    else
        for variable in CO2_VARIABLES
            _arrays_equal(path, source, (variable,)) && error(
                "$variable is unchanged from the 400-ppm source: $path")
        end
    end
    return path
end

function completed_scene_is_valid(state::BLState, grids)
    path = scene_path(state)
    isfile(path) || return false
    co2_source = state.sif_index == 2 ? no_sif_scene_path(state) : nothing
    validate_bottom_scene(state, path, grids; co2_source)
    return true
end

function acquire_claim(state::BLState)
    mkpath(CLAIMS_ROOT)
    claim = joinpath(CLAIMS_ROOT, @sprintf("state_%03d.claim", state.index))
    try
        mkdir(claim)
    catch error_value
        isdir(claim) && error(
            "state $(state.index) is already claimed at $claim; inspect owner.txt and never delete a live claim")
        rethrow(error_value)
    end
    token = string(uuid4())
    stage_dir = joinpath(STAGING_ROOT,
        @sprintf("state_%03d_%s", state.index, token))
    try
        mkpath(stage_dir)
        open(joinpath(claim, "owner.txt"), "w") do io
            println(io, "host=$(gethostname())")
            println(io, "pid=$(getpid())")
            println(io, "started_utc=$(now(UTC))")
            println(io, "state=$(state.index)")
            println(io, "partition=$BL_FILTER")
            println(io, "stage=$BL_STAGE")
            println(io, "token=$token")
            println(io, "staging=$stage_dir")
        end
    catch
        isdir(stage_dir) && rm(stage_dir; recursive=true)
        isdir(claim) && rm(claim; recursive=true)
        rethrow()
    end
    return (; claim, stage_dir,
            stage_file=joinpath(stage_dir,
                @sprintf("hiressim_%03d.nc", state.index)))
end

function record_failure!(claim, error_value)
    open(joinpath(claim.claim, "FAILED.txt"), "w") do io
        println(io, "failed_utc=$(now(UTC))")
        println(io, "error=$(sprint(showerror, error_value))")
    end
end

function release_claim!(claim; attempts::Int=3)
    attempts > 0 || throw(ArgumentError("claim-release attempts must be positive"))
    for path in (claim.stage_dir, claim.claim)
        for attempt in 1:attempts
            !ispath(path) && break
            try
                # Both paths are private to this producer.  `recursive=true`
                # also tolerates delayed directory entries left by a shared
                # filesystem while `force=true` makes an already-removed path
                # an idempotent success.
                rm(path; recursive=true, force=true)
                break
            catch
                attempt == attempts && rethrow()
                yield()
                sleep(0.05 * attempt)
            end
        end
    end
    return nothing
end

"""Release bookkeeping after a validated scene has already been published.

Claim cleanup is deliberately best-effort at this point.  Once the staging
file has passed validation and `mv` has committed it at the final path, a
transient filesystem cleanup error must not relabel that valid scene as a
failed simulation or abort the remaining queue.  The warning retains both
paths so a stale claim can be audited before it is retired manually.
"""
function release_published_claim!(state::BLState, claim;
                                  cleanup::Function=release_claim!)
    try
        cleanup(claim)
        return true
    catch error_value
        @warn("published bottom-layer truth scene, but claim cleanup failed; final output remains valid and the stale claim requires an owner audit",
              state=state.index,
              output=scene_path(state),
              claim=claim.claim,
              staging=claim.stage_dir,
              exception=(error_value, catch_backtrace()))
        return false
    end
end

function chunk_ranges_bottom(n::Int)
    return [first:min(first + CO2_CHUNK_POINTS_BOTTOM - 1, n)
            for first in 1:CO2_CHUNK_POINTS_BOTTOM:n]
end

function recompute_co2!(state::BLState, path::AbstractString, grids, solar_T)
    for (iband, band, variable) in (
            (2, "weak_co2", CO2_VARIABLES[1]),
            (3, "strong_co2", CO2_VARIABLES[2]))
        band_grid = grids[iband]
        ranges = chunk_ranges_bottom(length(band_grid))
        for (ichunk, indices) in enumerate(ranges)
            solve_grid = band_grid[indices]
            @info("bottom-layer CO2 truth",
                  state=state.index,
                  partition=BL_FILTER,
                  band=band,
                  bottom_co2_ppm=state.bottom_co2_ppm,
                  ichunk=ichunk,
                  nchunks=length(ranges),
                  nsolve=length(solve_grid))
            values = simulate_co2(
                state, iband, band_grid, solve_grid, solar_T)
            NCDataset(path, "a") do ds
                ds[variable][:, indices] = values
            end
            CUDA.synchronize()
            GC.gc()
            CUDA.reclaim()
        end
    end
    return path
end

function copy_co2_from_no_sif!(state::BLState, path::AbstractString, grids)
    source = no_sif_scene_path(state)
    paired_state = state_at_index(
        BottomLayerTruthCommon.paired_no_sif_index(state))
    validate_bottom_scene(paired_state, source, grids)
    NCDataset(source, "r") do src
        NCDataset(path, "a") do dst
            for variable in CO2_VARIABLES
                dst[variable][:, :] = src[variable][:, :]
            end
        end
    end
    return path
end

function publish_state!(state::BLState, grids, solar_T)
    if isfile(scene_path(state))
        completed_scene_is_valid(state, grids)
        @info("skip validated bottom-layer truth scene",
              state=state.index,
              path=scene_path(state))
        return scene_path(state)
    end

    source = source_scene_path(state)
    validate_source_scene(state, source, grids)
    claim = acquire_claim(state)
    try
        # Recheck after the atomic claim in case another producer published
        # between the first existence test and claim acquisition.
        isfile(scene_path(state)) && error(
            "final scene appeared after claim acquisition: $(scene_path(state))")
        cp(source, claim.stage_file)
        if BL_STAGE == "nosif"
            is_control = state.bottom_co2_ppm == state.background_co2_ppm
            mode = is_control ?
                "full accepted uniform-400 control reused" :
                "accepted O2 reused; weak/strong CO2 recomputed for layer-16 perturbation"
            co2_source = is_control ? source : nothing
            mark_staging!(claim.stage_file, state, source, mode;
                          co2_source=co2_source)
            is_control ||
                recompute_co2!(state, claim.stage_file, grids, solar_T)
        else
            paired = no_sif_scene_path(state)
            mark_staging!(claim.stage_file, state, source,
                "accepted SIF-on O2 reused; weak/strong CO2 reused from paired no-SIF bottom-layer scene";
                co2_source=paired)
            copy_co2_from_no_sif!(state, claim.stage_file, grids)
            co2_source = paired
        end
        mark_complete!(claim.stage_file)
        validate_bottom_scene(
            state, claim.stage_file, grids; co2_source=co2_source)
        mv(claim.stage_file, scene_path(state))
        release_published_claim!(state, claim)
        @info("published bottom-layer truth scene",
              state=state.index,
              path=scene_path(state))
        return scene_path(state)
    catch error_value
        record_failure!(claim, error_value)
        @error("bottom-layer truth state failed; claim and staging retained",
               state=state.index,
               claim=claim.claim,
               staging=claim.stage_dir,
               exception=(error_value, catch_backtrace()))
        rethrow(error_value)
    end
end

function print_plan(states, grids)
    compute_states = count(state -> BL_STAGE == "nosif" &&
        state.bottom_co2_ppm != state.background_co2_ppm, states)
    control_copies = count(state -> state.bottom_co2_ppm ==
        state.background_co2_ppm, states)
    no_sif_dependencies = filter(
        state -> state.sif_index == 1,
        BottomLayerTruthCommon.read_bottom_states(STATE_TABLE))
    dependencies_present = BL_STAGE == "sif" ? count(
        state -> isfile(campaign_scene_path(state)), no_sif_dependencies) : 0
    dependencies_required = BL_STAGE == "sif" ? 40 : 0
    @info("bottom-layer truth plan (source inputs validated below)",
          partition=BL_FILTER,
          stage=BL_STAGE,
          outputs=length(states),
          control_copies=control_copies,
          compute_states=compute_states,
          new_band_solves=2 * compute_states,
          global_no_sif_dependencies_present=dependencies_present,
          global_no_sif_dependencies_required=dependencies_required,
          execution_ready=(dependencies_present == dependencies_required),
          output=BL_OUT)
    for state in states
        source = source_scene_path(state)
        validate_source_scene(state, source, grids)
        action = BL_STAGE == "sif" ? "assemble from SIF-on O2 + no-SIF CO2" :
            state.bottom_co2_ppm == state.background_co2_ppm ?
                "copy uniform-400 control" : "recompute weak + strong CO2"
        @info("planned state",
              state=state.index,
              surface=state.surface,
              bottom_co2_ppm=state.bottom_co2_ppm,
              action=action,
              source_state=BottomLayerTruthCommon.old_control_index(state))
    end
end

function main_bottom_truth()
    states = selected_states()
    grids = output_grids()
    if PLAN_ONLY
        print_plan(states, grids)
        return
    end

    mkpath(BL_OUT)
    ensure_wavelength_file!()
    solar_T = nothing
    if BL_STAGE == "nosif" &&
            any(state -> !isfile(scene_path(state)) &&
                         state.bottom_co2_ppm != state.background_co2_ppm,
                states)
        CUDA.functional() || error("CUDA is not functional")
        CUDA.device!(DEVICE)
        solar_T = solar_interpolator()
    end
    if BL_STAGE == "sif"
        # This is a global campaign barrier, not merely a per-host check: all
        # 40 no-SIF clear+aerosol scenes must validate before either host may
        # begin assembling SIF-on outputs.
        no_sif_states = filter(
            state -> state.sif_index == 1,
            BottomLayerTruthCommon.read_bottom_states(STATE_TABLE))
        length(no_sif_states) == 40 || error(
            "expected 40 no-SIF states at the campaign barrier")
        for state in no_sif_states
            validate_bottom_scene(state, campaign_scene_path(state), grids)
        end
    end
    for state in states
        publish_state!(state, grids, solar_T)
    end
    all(state -> completed_scene_is_valid(state, grids), states) || error(
        "partition validation failed after publication")
    @info("bottom-layer truth partition complete",
          partition=BL_FILTER,
          stage=BL_STAGE,
          scenes=length(states),
          output=BL_OUT)
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main_bottom_truth()
