#!/usr/bin/env julia

"""
Regenerate only the O₂ A-band truth spectra while preserving the completed
weak and strong CO₂ arrays in every active scene file.

The retained A-band core is always sliced from the single canonical full-band
grid.  Each spectral chunk is then surrounded by independently constructed
Raman shoulders, so `solve_grid[keep] == canonical_core` exactly in both
Float32 and Float64.  O₂ is independent of XCO₂ in this experiment: one solve
for `(surface, aerosol case, SIF case)` is written to all four matching XCO₂
states.

Environment controls:

- `AEROSOL_CASE_FILTER`: `none` or `aerosol` (required)
- `TRUTH_OUT`: active destination directory
- `FIRST_STATE`, `LAST_STATE`: inclusive state interval (default 1:64)
- `O2_CHUNK_POINTS`: retained core points per Raman solve (default 64 for
  aerosol scenes and the full 2735-point band for clear scenes)
- `CUDA_DEVICE`: CUDA device index (default 1)
- `FORCE=1`: clear this producer's checkpoint and recompute selected chunks;
  it never recreates a scene file or writes either CO₂ variable
"""

include(joinpath(@__DIR__, "generate_truth_map.jl"))
using JLD2
using SHA

# Version 2 invalidated checkpoints made before the mandatory m>=3 Fourier
# convergence guard. Version 3 additionally invalidates every checkpoint made
# with the superseded wavelength-integrated-SIF normalization. Without this
# bump, a corrected SIF-on run could accept all old chunk tags and publish the
# old radiances under new state-table metadata without executing any RT.
const O2_REGEN_VERSION = 3
const O2_CORE_GRID_VERSION = 2
const SURFACE_COORDINATE_VERSION = 1
const CAMPAIGN_SIF_STATE = RRSXCO2Common.campaign_sif_state()
const DEFAULT_O2_CHUNK_POINTS = AEROSOL_CASE_FILTER == "aerosol" ? 64 : 2735
const O2_CHUNK_POINTS = parse(Int, get(
    ENV, "O2_CHUNK_POINTS", string(DEFAULT_O2_CHUNK_POINTS)))
const O2_CHECKPOINT = joinpath(
    OUT, "o2_exact_grid_$(AEROSOL_CASE_FILTER)_$(FIRST_STATE)_$(LAST_STATE)_checkpoint.jld2")
const O2_VARIABLES = (
    "radiance_rayleigh_o2a",
    "radiance_cabannes_o2a",
    "radiance_rrs_o2a",
)
const CO2_VARIABLES = (
    "radiance_rayleigh_weak_co2",
    "radiance_rayleigh_strong_co2",
)

AEROSOL_CASE_FILTER in ("none", "aerosol") || error(
    "AEROSOL_CASE_FILTER must be either none or aerosol")
O2_CHUNK_POINTS > 0 || error("O2_CHUNK_POINTS must be positive")

scene_path(state) = joinpath(OUT, @sprintf("hiressim_%03d.nc", state.index))
chunk_ranges(n::Int) = [i:min(i + O2_CHUNK_POINTS - 1, n)
                        for i in 1:O2_CHUNK_POINTS:n]

function o2_groups(states)
    groups = Dict{Tuple{Int,Int},Vector{TruthState}}()
    for state in states
        # There is no CO₂ opacity in the A-band configuration, so the four
        # XCO₂ states share exactly one O₂ calculation.
        push!(get!(groups, (state.surface_index, state.sif_index), TruthState[]),
              state)
    end
    return groups
end

function save_checkpoint!(completed::Set{String})
    tmp = O2_CHECKPOINT * ".tmp"
    jldsave(tmp;
        completed=sort!(collect(completed)),
        version=O2_REGEN_VERSION,
        aerosol_case_filter=AEROSOL_CASE_FILTER,
        first_state=FIRST_STATE,
        last_state=LAST_STATE,
        o2_chunk_points=O2_CHUNK_POINTS,
        float_type=string(FT),
        nstreams=AEROSOL_CASE_FILTER == "aerosol" ? AEROSOL_NSTREAMS : 5,
        psurf_hpa=1000.0,
        nlayers=NLAYERS,
        sza_deg=SZA_DEG,
        vza_deg=VZA_DEG,
        relative_azimuth_deg=RELATIVE_AZIMUTH_DEG,
        sif_definition_version=RRSXCO2Common.SIF_DEFINITION_VERSION,
        sif_case_on=RRSXCO2Common.SIF_CASE_ON,
        sif_angular_integral_760=
            RRSXCO2Common.SIF_ANGULAR_INTEGRAL_760,
        sif_radiance_760=RRSXCO2Common.SIF_RADIANCE_760,
        sif760_native=CAMPAIGN_SIF_STATE.SIF760,
        msif_native=CAMPAIGN_SIF_STATE.mSIF,
        sif_template_wavelength_integral=
            CAMPAIGN_SIF_STATE.wavelength_integral,
        raman_shoulder_cm=SHOULDER_CM,
        o2_core_grid_version=O2_CORE_GRID_VERSION,
        surface_coordinate_version=SURFACE_COORDINATE_VERSION,
        absco_version=RRSXCO2Common.ABSCO_VERSION)
    mv(tmp, O2_CHECKPOINT; force=true)
end

function load_checkpoint()
    isfile(O2_CHECKPOINT) || return Set{String}()
    saved = load(O2_CHECKPOINT)
    checks = (
        ("version", O2_REGEN_VERSION),
        ("aerosol_case_filter", AEROSOL_CASE_FILTER),
        ("first_state", FIRST_STATE),
        ("last_state", LAST_STATE),
        ("o2_chunk_points", O2_CHUNK_POINTS),
        ("float_type", string(FT)),
        ("psurf_hpa", 1000.0),
        ("nlayers", NLAYERS),
        ("sza_deg", SZA_DEG),
        ("vza_deg", VZA_DEG),
        ("relative_azimuth_deg", RELATIVE_AZIMUTH_DEG),
        ("sif_definition_version", RRSXCO2Common.SIF_DEFINITION_VERSION),
        ("sif_case_on", RRSXCO2Common.SIF_CASE_ON),
        ("sif_angular_integral_760",
         RRSXCO2Common.SIF_ANGULAR_INTEGRAL_760),
        ("sif_radiance_760", RRSXCO2Common.SIF_RADIANCE_760),
        ("sif760_native", CAMPAIGN_SIF_STATE.SIF760),
        ("msif_native", CAMPAIGN_SIF_STATE.mSIF),
        ("sif_template_wavelength_integral",
         CAMPAIGN_SIF_STATE.wavelength_integral),
        ("raman_shoulder_cm", SHOULDER_CM),
        ("o2_core_grid_version", O2_CORE_GRID_VERSION),
        ("surface_coordinate_version", SURFACE_COORDINATE_VERSION),
        ("absco_version", RRSXCO2Common.ABSCO_VERSION),
    )
    for (name, expected) in checks
        get(saved, name, nothing) == expected || error(
            "O2 checkpoint field $name differs from this run; use FORCE=1")
    end
    return Set{String}(saved["completed"])
end

function finite_complete_variable(dataset, name)
    haskey(dataset, name) || error("active truth scene is missing $name")
    values = Float64.(nomissing(dataset[name][:, :], NaN))
    all(isfinite, values) || error("$name contains non-finite values")
    maximum(abs, values) < 1e30 || error("$name contains unwritten fill values")
    return nothing
end

"Validate that a state row carries exactly the current campaign SIF definition."
function validate_state_sif(state)
    if state.sif_index == 1
        state.sif_case == "off" || error(
            "state $(state.index) has SIF index 1 but case $(state.sif_case)")
        all(iszero, (state.sif_angular_integral760,
                     state.sif760, state.msif)) || error(
            "state $(state.index) has nonzero SIF coordinates in the off case")
        return false
    elseif state.sif_index == 2
        state.sif_case == RRSXCO2Common.SIF_CASE_ON || error(
            "state $(state.index) uses stale/unknown SIF case $(state.sif_case); " *
            "expected $(RRSXCO2Common.SIF_CASE_ON)")
        state.sif_angular_integral760 ==
            RRSXCO2Common.SIF_ANGULAR_INTEGRAL_760 || error(
            "state $(state.index) has the wrong 760-nm SIF angular integral")
        isapprox(state.sif760, CAMPAIGN_SIF_STATE.SIF760;
                 atol=5e-15, rtol=0) || error(
            "state $(state.index) has stale SIF760=$(state.sif760)")
        isapprox(state.msif, CAMPAIGN_SIF_STATE.mSIF;
                 atol=5e-16, rtol=0) || error(
            "state $(state.index) has stale mSIF=$(state.msif)")
        return true
    end
    error("state $(state.index) has unsupported SIF index $(state.sif_index)")
end

"Rewrite, then verify, every version-2 campaign SIF attribute on a scene."
function write_and_validate_sif_provenance!(attributes, state, path)
    enabled = validate_state_sif(state)
    attributes["sif_case"] = state.sif_case
    RRSXCO2Common.write_sif_provenance!(attributes, enabled)

    # This legacy attribute encoded the superseded wavelength-integrated 0.5
    # control variable. Leaving it in an otherwise corrected file would make
    # the metadata internally contradictory.
    haskey(attributes, "sif_total_mW_m-2_sr-1") &&
        delete!(attributes, "sif_total_mW_m-2_sr-1")

    String(attributes["sif_case"]) == state.sif_case || error(
        "failed to rewrite SIF case metadata in $path")
    Int(attributes["sif_definition_version"]) ==
        RRSXCO2Common.SIF_DEFINITION_VERSION || error(
        "failed to rewrite SIF definition version in $path")
    String(attributes["sif_case_on_label"]) ==
        RRSXCO2Common.SIF_CASE_ON || error(
        "failed to rewrite the SIF-on label in $path")
    String(attributes["sif_definition"]) ==
        "isotropic BOA radiance normalized by 2pi*L_lambda(760 nm)=0.5" ||
        error("failed to rewrite the SIF definition in $path")

    active = enabled ? 1.0 : 0.0
    expected = (
        ("sif_reference_wavelength_nm",
         RRSXCO2Common.SIF_REFERENCE_WAVELENGTH_NM),
        ("sif_upwelling_solid_angle_sr",
         RRSXCO2Common.SIF_UPWELLING_SOLID_ANGLE_SR),
        ("sif_angular_integral_760_mW_m-2_nm-1",
         active * RRSXCO2Common.SIF_ANGULAR_INTEGRAL_760),
        ("sif_radiance_760_mW_m-2_sr-1_nm-1",
         active * CAMPAIGN_SIF_STATE.radiance_760),
        ("sif_cosine_weighted_irradiance_760_mW_m-2_nm-1",
         active * pi * CAMPAIGN_SIF_STATE.radiance_760),
        ("sif_SIF760_mW_m-2_sr-1_per_cm-1",
         active * CAMPAIGN_SIF_STATE.SIF760),
        ("sif_mSIF_mW_m-2_sr-1_per_cm-2",
         active * CAMPAIGN_SIF_STATE.mSIF),
        ("sif_template_wavelength_integral_mW_m-2_sr-1",
         active * CAMPAIGN_SIF_STATE.wavelength_integral),
    )
    for (name, value) in expected
        isapprox(Float64(attributes[name]), Float64(value);
                 atol=5e-15, rtol=0) || error(
            "failed to rewrite SIF provenance attribute $name in $path")
    end
    return enabled
end

function prepare_scene!(state)
    path = scene_path(state)
    isfile(path) || error(
        "missing active scene $path; this producer will not recreate files " *
        "because doing so could destroy the completed CO2 bands")
    validate_state_sif(state)
    NCDataset(path, "a") do ds
        Int(ds.attrib["state_index"]) == state.index || error(
            "state-index mismatch in $path")
        get(ds.attrib, "co2_absco_regeneration_complete", 0) == 1 || error(
            "completed CO2 regeneration is required before O2 replacement: $path")
        for variable in CO2_VARIABLES
            finite_complete_variable(ds, variable)
        end
        for variable in O2_VARIABLES
            haskey(ds, variable) || error("active truth scene is missing $variable")
        end

        # This is set before any replacement write. Downstream consumers can
        # never mistake an interrupted partially regenerated file for truth.
        ds.attrib["o2_absco_regeneration_complete"] = Int32(0)
        ds.attrib["o2_truth_regenerated"] = Int32(0)
        ds.attrib["simulation_complete"] = Int32(0)
        if haskey(ds.attrib, "chunked_simulation_complete")
            ds.attrib["chunked_simulation_complete"] = Int32(0)
        end
        write_and_validate_sif_provenance!(ds.attrib, state, path)
    end
    return path
end

function write_o2_chunk!(states, indices, result)
    data = (result.rayleigh, result.cabannes, result.rrs)
    for state in states
        NCDataset(scene_path(state), "a") do ds
            for (variable, values) in zip(O2_VARIABLES, data)
                ds[variable][:, indices] = values
            end
        end
    end
end

function finalize_scene!(state)
    NCDataset(scene_path(state), "a") do ds
        for variable in (O2_VARIABLES..., CO2_VARIABLES...)
            finite_complete_variable(ds, variable)
        end
        RRSXCO2Common.write_absco_provenance!(ds.attrib)
        RRSXCO2Common.write_fourier_convergence_provenance!(ds.attrib)
        write_and_validate_sif_provenance!(ds.attrib, state, scene_path(state))
        # Seed files were copied from the superseded campaign and therefore
        # still point at its state table.  A corrected scene must identify the
        # exact version-2 table that drove this direct solve; the release gate
        # rejects stale source-table provenance.  The checksum is portable
        # across hosts even when the absolute producer path is not.
        ds.attrib["source_state_table"] = abspath(STATE_FILE)
        ds.attrib["source_state_table_sha256"] = open(STATE_FILE, "r") do io
            bytes2hex(sha256(io))
        end
        ds.attrib["o2_truth_reused"] = Int32(0)
        ds.attrib["o2_truth_reuse_source"] = "none; regenerated from current model"
        ds.attrib["o2_truth_regenerated"] = Int32(1)
        ds.attrib["o2_absco_regeneration_complete"] = Int32(1)
        ds.attrib["o2_truth_spectroscopy"] =
            "ABSCO v5.2 O2 + rebuilt HITRAN H2O; O2 VMR=0.21"
        ds.attrib["o2_core_grid_version"] = Int32(O2_CORE_GRID_VERSION)
        ds.attrib["o2_core_grid_construction"] =
            "canonical retained output nodes inserted verbatim between independently constructed Raman shoulders"
        ds.attrib["o2_truth_grid_note"] =
            "retained nodes are exactly identical to the canonical 0.1 cm-1 output grid"
        ds.attrib["o2_chunk_points"] = Int32(O2_CHUNK_POINTS)
        ds.attrib["o2_raman_shoulder_cm-1"] = Float64(SHOULDER_CM)
        ds.attrib["o2_nstreams"] = Int32(
            AEROSOL_CASE_FILTER == "aerosol" ? AEROSOL_NSTREAMS : 5)
        ds.attrib["surface_coordinate"] =
            "Legendre x is defined once over each complete output band before chunking"
        ds.attrib["o2_absco_completed"] = string(now())
        ds.attrib["simulation_complete"] = Int32(1)
        if haskey(ds.attrib, "chunked_simulation_complete")
            ds.attrib["chunked_simulation_complete"] = Int32(1)
        end
        ds.attrib["completed"] = string(now())
    end
end

function run_o2_regeneration()
    mkpath(OUT)
    CUDA.functional() || error("CUDA is not functional")
    CUDA.device!(DEVICE)
    states = selected_states()
    expected_aerosol = AEROSOL_CASE_FILTER == "aerosol"
    all(state -> any(>(0), state.aod550) == expected_aerosol, states) || error(
        "selected states do not match AEROSOL_CASE_FILTER")

    if FORCE && isfile(O2_CHECKPOINT)
        mv(O2_CHECKPOINT, O2_CHECKPOINT * ".previous"; force=true)
    end
    completed = FORCE ? Set{String}() : load_checkpoint()
    for state in states
        prepare_scene!(state)
    end

    output_grid = output_grids()[1]
    ranges = chunk_ranges(length(output_grid))
    solar_T = solar_interpolator()
    for (key, members) in sort!(collect(o2_groups(states)); by=first)
        representative = first(members)
        for (ichunk, indices) in enumerate(ranges)
            tag = "o2_surface$(key[1])_sif$(key[2])_chunk$(ichunk)"
            tag in completed && continue
            core_grid = output_grid[indices]
            solve_grid, keep = o2_solve_grid(core_grid)
            solve_grid[keep] == core_grid || error(
                "Raman solve grid changed canonical retained nodes for $tag")
            @info "exact-grid O2 truth chunk" AEROSOL_CASE_FILTER key ichunk nchunks=length(ranges) ncore=length(core_grid) nsolve=length(solve_grid)
            # The complete canonical band remains the surface/aerosol spectral
            # anchor even though this solve contains only one retained chunk.
            result = simulate_o2(
                representative, (output_grid, solve_grid, keep), solar_T)
            write_o2_chunk!(members, indices, result)
            push!(completed, tag)
            save_checkpoint!(completed)
            CUDA.synchronize()
            GC.gc()
            CUDA.reclaim()
        end
    end

    for state in states
        finalize_scene!(state)
    end
    @info "completed exact-grid O2 truth regeneration" AEROSOL_CASE_FILTER scenes=length(states) groups=length(o2_groups(states)) chunks_per_group=length(ranges) output=OUT checkpoint=O2_CHECKPOINT
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && run_o2_regeneration()
