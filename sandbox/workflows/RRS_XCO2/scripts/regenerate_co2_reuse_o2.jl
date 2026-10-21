#!/usr/bin/env julia

"""
Regenerate only the weak and strong CO₂ bands while retaining the archived
O₂ A-band truth spectra.

The archived A-band calculations already use ABSCO v5.2 O₂, O₂ VMR 0.21,
the common p/T/q profile, and the required Raman solve. This driver copies a
scene into the active truth directory, leaves its three A-band variables
unchanged, and replaces only the two CO₂-band variables using the current
band-specific ABSCO v5.2 tables. Files remain marked incomplete until both
CO₂ bands have been written.

Environment controls:

- `AEROSOL_CASE_FILTER`: `none` or `aerosol` (required)
- `TRUTH_OUT`: destination directory
- `REUSE_TRUTH_ROOT`: archived truth root
- `FIRST_STATE`, `LAST_STATE`: inclusive state range
- `CO2_CHUNK_POINTS`: retained points per solve (default 2048, normally one
  full-band solve)
- `WRITE_WAVELENGTHS=0` skips the shared wavelength file when several
  disjoint state ranges run concurrently (default 1)
- `CUDA_DEVICE`: CUDA device index (default 1)
- `FORCE=1`: restore selected scenes from the archive and clear this driver's
  checkpoint before recomputing
"""

include(joinpath(@__DIR__, "generate_truth_map.jl"))
using JLD2

const REUSE_TRUTH_ROOT = get(ENV, "REUSE_TRUTH_ROOT", joinpath(
    RRSXCO2Common.ROOT, "obsolete", "pre_absco_closure_20260829", "truth_map"))
const CO2_CHUNK_POINTS = parse(Int, get(ENV, "CO2_CHUNK_POINTS", "2048"))
const WRITE_WAVELENGTHS = lowercase(get(ENV, "WRITE_WAVELENGTHS", "1")) in
                          ("1", "true", "yes", "on")
const CO2_REUSE_VERSION = 1

AEROSOL_CASE_FILTER in ("none", "aerosol") || error(
    "AEROSOL_CASE_FILTER must be either none or aerosol")
CO2_CHUNK_POINTS > 0 || error("CO2_CHUNK_POINTS must be positive")

const REUSE_SOURCE = AEROSOL_CASE_FILTER == "aerosol" ?
    joinpath(REUSE_TRUTH_ROOT, "aerosol_chunked") : REUSE_TRUTH_ROOT
const CHECKPOINT = joinpath(OUT,
    "co2_absco_$(AEROSOL_CASE_FILTER)_$(FIRST_STATE)_$(LAST_STATE)_checkpoint.jld2")

chunk_ranges(n::Int) = [i:min(i + CO2_CHUNK_POINTS - 1, n)
                        for i in 1:CO2_CHUNK_POINTS:n]
scene_path(state) = joinpath(OUT, @sprintf("hiressim_%03d.nc", state.index))
source_scene_path(state) = joinpath(
    REUSE_SOURCE, @sprintf("hiressim_%03d.nc", state.index))

function co2_groups(states)
    groups = Dict{Tuple{Int,Int},Vector{TruthState}}()
    for state in states
        # SIF is absent from both CO2 bands, so its two truth states share one
        # calculation for a given surface and XCO2 value.
        push!(get!(groups, (state.surface_index, state.xco2_index), TruthState[]),
              state)
    end
    return groups
end

function save_checkpoint!(completed::Set{String})
    tmp = CHECKPOINT * ".tmp"
    jldsave(tmp;
        completed=sort!(collect(completed)),
        version=CO2_REUSE_VERSION,
        aerosol_case_filter=AEROSOL_CASE_FILTER,
        first_state=FIRST_STATE,
        last_state=LAST_STATE,
        chunk_points=CO2_CHUNK_POINTS,
        absco_version=RRSXCO2Common.ABSCO_VERSION)
    mv(tmp, CHECKPOINT; force=true)
end

function load_checkpoint()
    isfile(CHECKPOINT) || return Set{String}()
    saved = load(CHECKPOINT)
    saved["version"] == CO2_REUSE_VERSION || error("stale CO2-only checkpoint")
    saved["aerosol_case_filter"] == AEROSOL_CASE_FILTER || error(
        "checkpoint aerosol selection differs")
    saved["first_state"] == FIRST_STATE || error("checkpoint first state differs")
    saved["last_state"] == LAST_STATE || error("checkpoint last state differs")
    saved["chunk_points"] == CO2_CHUNK_POINTS || error("checkpoint chunk size differs")
    saved["absco_version"] == RRSXCO2Common.ABSCO_VERSION || error(
        "checkpoint ABSCO version differs")
    return Set{String}(saved["completed"])
end

function prepare_scene!(state)
    source = source_scene_path(state)
    destination = scene_path(state)
    isfile(source) || error("missing archived truth scene: $source")
    if FORCE || !isfile(destination)
        cp(source, destination; force=true)
    end
    NCDataset(destination, "a") do ds
        # Preserve the archived A-band values while recording their exact
        # spectroscopy separately from the regenerated CO2 bands.
        RRSXCO2Common.write_absco_provenance!(ds.attrib)
        ds.attrib["o2_truth_reused"] = Int32(1)
        ds.attrib["o2_truth_reuse_source"] = abspath(source)
        ds.attrib["o2_truth_spectroscopy"] =
            "ABSCO v5.2 O2 + rebuilt HITRAN H2O; O2 VMR=0.21"
        ds.attrib["o2_truth_grid_note"] = AEROSOL_CASE_FILTER == "aerosol" ?
            "archived 256-point Raman chunks; retained Float32 nodes may differ from the nominal grid by <=9.765625e-4 cm-1" :
            "archived full-band Raman solve; retained nodes equal the nominal grid"
        ds.attrib["co2_absco_regeneration_complete"] = Int32(0)
        ds.attrib["simulation_complete"] = Int32(0)
        if haskey(ds.attrib, "chunked_simulation_complete")
            ds.attrib["chunked_simulation_complete"] = Int32(0)
        end
    end
    return destination
end

function write_chunk!(states, band::String, indices, values)
    variable = "radiance_rayleigh_$band"
    for state in states
        NCDataset(scene_path(state), "a") do ds
            ds[variable][:, indices] = values
        end
    end
end

function run_band!(states, iband::Int, band::String, band_grid, solar_T,
                   completed::Set{String})
    for (key, members) in sort!(collect(co2_groups(states)); by=first)
        representative = first(members)
        for (ichunk, indices) in enumerate(chunk_ranges(length(band_grid)))
            tag = "$(band)_surface$(key[1])_xco2$(key[2])_chunk$(ichunk)"
            tag in completed && continue
            solve_grid = band_grid[indices]
            @info "CO2-only ABSCO truth" AEROSOL_CASE_FILTER band key ichunk nchunks=length(chunk_ranges(length(band_grid))) nsolve=length(solve_grid)
            values = simulate_co2(
                representative, iband, band_grid, solve_grid, solar_T)
            write_chunk!(members, band, indices, values)
            push!(completed, tag)
            save_checkpoint!(completed)
            CUDA.synchronize()
            GC.gc()
            CUDA.reclaim()
        end
    end
end

function main_co2_reuse()
    mkpath(OUT)
    CUDA.functional() || error("CUDA is not functional")
    CUDA.device!(DEVICE)
    states = selected_states()
    expected_aerosol = AEROSOL_CASE_FILTER == "aerosol"
    all(state -> any(>(0), state.aod550) == expected_aerosol, states) || error(
        "selected state set does not match AEROSOL_CASE_FILTER")

    if FORCE && isfile(CHECKPOINT)
        mv(CHECKPOINT, CHECKPOINT * ".previous"; force=true)
    end
    completed = FORCE ? Set{String}() : load_checkpoint()
    for state in states
        prepare_scene!(state)
    end
    grids = output_grids()
    WRITE_WAVELENGTHS && write_wavelengths(grids)
    solar_T = solar_interpolator()
    run_band!(states, 2, "weak_co2", grids[2], solar_T, completed)
    run_band!(states, 3, "strong_co2", grids[3], solar_T, completed)

    for state in states
        NCDataset(scene_path(state), "a") do ds
            ds.attrib["co2_absco_regeneration_complete"] = Int32(1)
            ds.attrib["simulation_complete"] = Int32(1)
            if haskey(ds.attrib, "chunked_simulation_complete")
                ds.attrib["chunked_simulation_complete"] = Int32(1)
            end
            ds.attrib["co2_absco_completed"] = string(now())
        end
    end
    @info "completed CO2-only ABSCO truth regeneration" AEROSOL_CASE_FILTER scenes=length(states) output=OUT
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main_co2_reuse()
