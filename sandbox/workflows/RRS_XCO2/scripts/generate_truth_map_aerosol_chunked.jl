#!/usr/bin/env julia

"""
Generate the aerosol-on half of the RRS–XCO2 truth map in spectral chunks.

The O2 A-band output grid is divided into core chunks. Every core is solved on
an expanded grid with `SHOULDER_CM` on both sides so Raman source wavelengths
outside the retained interval are present; only the core indices are written.
The weak and strong CO2 bands contain no RRS and are chunked without shoulders.

Results are written directly into full-size `hiressim_NNN.nc` files. A JLD2
checkpoint records completed `(band, physical-state key, chunk)` units, making
the calculation resumable without retaining spectra in memory.

Environment controls:

- `TRUTH_OUT` (default `RRS_XCO2/truth_map/aerosol_chunked`)
- `O2_CHUNK_POINTS` (default 256 output points)
- `CO2_CHUNK_POINTS` (default 512 output points)
- `RAMAN_SHOULDER_CM` (default 250 cm-1 on each side)
- `CUDA_DEVICE` (default 1)
- `FIRST_STATE`, `LAST_STATE` (default 1, 64; aerosol-off states are skipped)
- `FORCE=1` recreates output files and clears the checkpoint
"""

const AEROSOL_SCRIPT_ROOT = normpath(joinpath(@__DIR__, ".."))
ENV["TRUTH_OUT"] = get(ENV, "TRUTH_OUT",
    joinpath(AEROSOL_SCRIPT_ROOT, "truth_map", "aerosol_chunked"))
include(joinpath(@__DIR__, "generate_truth_map.jl"))
using JLD2

const AEROSOL_OUT = OUT
const AEROSOL_CHECKPOINT = joinpath(AEROSOL_OUT, "aerosol_chunks_checkpoint.jld2")
const O2_CHUNK_POINTS = parse(Int, get(ENV, "O2_CHUNK_POINTS", "256"))
const CO2_CHUNK_POINTS = parse(Int, get(ENV, "CO2_CHUNK_POINTS", "512"))

O2_CHUNK_POINTS > 0 || error("O2_CHUNK_POINTS must be positive")
CO2_CHUNK_POINTS > 0 || error("CO2_CHUNK_POINTS must be positive")

chunk_ranges(n::Int, width::Int) =
    [i:min(i + width - 1, n) for i in 1:width:n]

function checkpoint!(completed::Set{String})
    tmp = AEROSOL_CHECKPOINT * ".tmp"
    jldsave(tmp; completed=sort!(collect(completed)),
            o2_chunk_points=O2_CHUNK_POINTS,
            co2_chunk_points=CO2_CHUNK_POINTS,
            psurf_hpa=1000.0, nlayers=NLAYERS,
            float_type=string(FT))
    mv(tmp, AEROSOL_CHECKPOINT; force=true)
end

function load_checkpoint()
    isfile(AEROSOL_CHECKPOINT) || return Set{String}()
    saved = load(AEROSOL_CHECKPOINT)
    saved["o2_chunk_points"] == O2_CHUNK_POINTS ||
        error("checkpoint O2 chunk size differs from requested size")
    saved["co2_chunk_points"] == CO2_CHUNK_POINTS ||
        error("checkpoint CO2 chunk size differs from requested size")
    saved["psurf_hpa"] == 1000.0 || error("checkpoint surface pressure is stale")
    saved["nlayers"] == NLAYERS || error("checkpoint layer count differs")
    saved["float_type"] == string(FT) || error("checkpoint float type differs")
    return Set{String}(saved["completed"])
end

scene_path(state) = joinpath(AEROSOL_OUT, @sprintf("hiressim_%03d.nc", state.index))

function initialize_scene!(state, grids)
    path = scene_path(state)
    if isfile(path) && !FORCE
        return path
    end
    NCDataset(path, "c") do ds
        defDim(ds, "stokes", 3)
        for (name, ν) in zip(("o2a", "weak_co2", "strong_co2"), grids)
            defDim(ds, name, length(ν))
        end
        for name in ("radiance_rayleigh_o2a", "radiance_cabannes_o2a",
                     "radiance_rrs_o2a")
            v = defVar(ds, name, Float32, ("stokes", "o2a"))
            v.attrib["units"] = "mW m-2 sr-1 (cm-1)-1"
        end
        for (name, band) in (("radiance_rayleigh_weak_co2", "weak_co2"),
                             ("radiance_rayleigh_strong_co2", "strong_co2"))
            v = defVar(ds, name, Float32, ("stokes", band))
            v.attrib["units"] = "mW m-2 sr-1 (cm-1)-1"
        end
        ds.attrib["state_index"] = Int32(state.index)
        ds.attrib["surface"] = state.surface
        ds.attrib["aerosol_case"] = state.aerosol_case
        ds.attrib["sif_case"] = state.sif_case
        ds.attrib["xco2_ppm"] = state.xco2_ppm
        ds.attrib["psurf_hpa"] = 1000.0
        ds.attrib["sza_deg"] = 30.0
        ds.attrib["vza_deg"] = 0.0
        ds.attrib["atmospheric_layers"] = Int32(NLAYERS)
        ds.attrib["aod550_sulfate"] = state.aod550[1]
        ds.attrib["aod550_organic_carbon"] = state.aod550[2]
        ds.attrib["aod550_stratospheric"] = state.aod550[3]
        ds.attrib["sif_total_mW_m-2_sr-1"] = state.sif_total
        ds.attrib["source_state_table"] = "true_states.dat"
        ds.attrib["spectral_chunking"] =
            "O2 cores carry ±$(SHOULDER_CM) cm-1 Raman shoulders; CO2 has no shoulders"
        ds.attrib["created"] = string(now())
    end
    return path
end

function write_o2_chunk!(states, irange, result)
    for state in states
        NCDataset(scene_path(state), "a") do ds
            ds["radiance_rayleigh_o2a"][:, irange] = result.rayleigh
            ds["radiance_cabannes_o2a"][:, irange] = result.cabannes
            ds["radiance_rrs_o2a"][:, irange] = result.rrs
        end
    end
end

function write_co2_chunk!(states, band::String, irange, result)
    varname = "radiance_rayleigh_$(band)"
    for state in states
        NCDataset(scene_path(state), "a") do ds
            ds[varname][:, irange] = result
        end
    end
end

function o2_groups(states)
    groups = Dict{Tuple{Int,Int},Vector{TruthState}}()
    for state in states
        # O2 contains no CO2 absorption, so XCO2 is not part of this key.
        push!(get!(groups, (state.surface_index, state.sif_index), TruthState[]), state)
    end
    groups
end

function co2_groups(states)
    groups = Dict{Tuple{Int,Int},Vector{TruthState}}()
    for state in states
        # SIF is absent from both CO2 bands.
        push!(get!(groups, (state.surface_index, state.xco2_index), TruthState[]), state)
    end
    groups
end

function run_o2_chunks!(states, output_ν, solar_T, completed)
    ranges = chunk_ranges(length(output_ν), O2_CHUNK_POINTS)
    for (key, members) in sort!(collect(o2_groups(states)); by=first)
        representative = first(members)
        for (ichunk, irange) in enumerate(ranges)
            tag = "o2_s$(key[1])_sif$(key[2])_chunk$(ichunk)"
            tag in completed && continue
            coreν = output_ν[irange]
            solveν, keep = o2_solve_grid(coreν)
            @info "aerosol O2 Raman chunk" key ichunk nchunks=length(ranges) ncore=length(coreν) nsolve=length(solveν)
            result = simulate_o2(representative, (coreν, solveν, keep), solar_T)
            write_o2_chunk!(members, irange, result)
            push!(completed, tag)
            checkpoint!(completed)
            CUDA.synchronize(); GC.gc(); CUDA.reclaim()
        end
    end
end

function run_co2_chunks!(states, iband, band, output_ν, solar_T, completed)
    ranges = chunk_ranges(length(output_ν), CO2_CHUNK_POINTS)
    for (key, members) in sort!(collect(co2_groups(states)); by=first)
        representative = first(members)
        for (ichunk, irange) in enumerate(ranges)
            tag = "$(band)_s$(key[1])_xco2$(key[2])_chunk$(ichunk)"
            tag in completed && continue
            ν = output_ν[irange]
            @info "aerosol CO2 chunk" band key ichunk nchunks=length(ranges) nsolve=length(ν)
            result = simulate_co2(representative, iband, ν, solar_T)
            write_co2_chunk!(members, band, irange, result)
            push!(completed, tag)
            checkpoint!(completed)
            CUDA.synchronize(); GC.gc(); CUDA.reclaim()
        end
    end
end

function main_aerosol_chunked()
    mkpath(AEROSOL_OUT)
    CUDA.functional() || error("CUDA is not functional")
    CUDA.device!(DEVICE)
    all_states = read_states()
    selected = all_states[FIRST_STATE:LAST_STATE]
    states = filter(s -> any(>(0), s.aod550), selected)
    isempty(states) && error("selected state range contains no aerosol-on scenes")

    grids = output_grids()
    if FORCE
        isfile(AEROSOL_CHECKPOINT) && mv(AEROSOL_CHECKPOINT,
            AEROSOL_CHECKPOINT * ".previous"; force=true)
    end
    completed = FORCE ? Set{String}() : load_checkpoint()
    if !isempty(completed)
        missing = filter(s -> !isfile(scene_path(s)), states)
        isempty(missing) || error(
            "checkpoint exists but scene files are missing: " *
            join((string(s.index) for s in missing), ", "))
    end
    for state in states
        initialize_scene!(state, grids)
    end
    # Keep a local wavelength file beside this independently runnable subset.
    write_wavelengths(grids)

    solar_T = solar_interpolator()
    run_o2_chunks!(states, grids[1], solar_T, completed)
    run_co2_chunks!(states, 2, "weak_co2", grids[2], solar_T, completed)
    run_co2_chunks!(states, 3, "strong_co2", grids[3], solar_T, completed)
    for state in states
        NCDataset(scene_path(state), "a") do ds
            ds.attrib["chunked_simulation_complete"] = Int32(1)
            ds.attrib["completed"] = string(now())
        end
    end
    @info "completed chunked aerosol truth map" scenes=length(states) output=AEROSOL_OUT
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main_aerosol_chunked()
