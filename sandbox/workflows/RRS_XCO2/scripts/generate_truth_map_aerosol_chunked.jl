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
- `RAMAN_SHOULDER_CM` (default 234 cm-1 on each side)
- `CUDA_DEVICE` (default 1)
- `FIRST_STATE`, `LAST_STATE` (default 1, 64; aerosol-off states are skipped)
- `SIF_CASE_FILTER` (`all`, `off`, or `on`; default `all`)
- `TRUTH_SZA_DEG`, `TRUTH_VZA_DEG`, `TRUTH_RELATIVE_AZIMUTH_DEG`
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
const SURFACE_COORDINATE_VERSION = 1
# Version 2 means that the retained O2 core is inserted verbatim between
# independently constructed Raman shoulders.  In particular, a Float32
# range that starts at the left shoulder is no longer allowed to reconstruct
# (and slightly shift) the output-band nodes.
const O2_CORE_GRID_VERSION = 2
const SIF_DEFINITION_VERSION = RRSXCO2Common.SIF_DEFINITION_VERSION

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
            float_type=string(FT),
            sza_deg=SZA_DEG, vza_deg=VZA_DEG,
            relative_azimuth_deg=RELATIVE_AZIMUTH_DEG,
            sif_case_filter=SIF_CASE_FILTER,
            sif_definition_version=SIF_DEFINITION_VERSION,
            sif_case_on=RRSXCO2Common.SIF_CASE_ON,
            sif_angular_integral_760=
                RRSXCO2Common.SIF_ANGULAR_INTEGRAL_760,
            sif_radiance_760=RRSXCO2Common.SIF_RADIANCE_760,
            surface_coordinate_version=SURFACE_COORDINATE_VERSION,
            o2_core_grid_version=O2_CORE_GRID_VERSION,
            raman_shoulder_cm=SHOULDER_CM,
            absco_version=RRSXCO2Common.ABSCO_VERSION)
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
    get(saved, "sza_deg", FT(30)) == SZA_DEG || error("checkpoint SZA differs")
    get(saved, "vza_deg", zero(FT)) == VZA_DEG || error("checkpoint VZA differs")
    get(saved, "relative_azimuth_deg", zero(FT)) == RELATIVE_AZIMUTH_DEG ||
        error("checkpoint relative azimuth differs")
    get(saved, "sif_case_filter", "all") == SIF_CASE_FILTER ||
        error("checkpoint SIF selection differs")
    get(saved, "sif_definition_version", 0) == SIF_DEFINITION_VERSION ||
        error("checkpoint predates the corrected 760-nm SIF normalization; restart with FORCE=1")
    get(saved, "sif_case_on", "") == RRSXCO2Common.SIF_CASE_ON ||
        error("checkpoint SIF case label differs")
    get(saved, "sif_angular_integral_760", NaN) ==
        RRSXCO2Common.SIF_ANGULAR_INTEGRAL_760 ||
        error("checkpoint SIF angular integral differs")
    get(saved, "sif_radiance_760", NaN) ==
        RRSXCO2Common.SIF_RADIANCE_760 ||
        error("checkpoint SIF stream radiance differs")
    get(saved, "surface_coordinate_version", 0) == SURFACE_COORDINATE_VERSION ||
        error("checkpoint predates full-band surface normalization; restart with FORCE=1")
    get(saved, "o2_core_grid_version", 0) == O2_CORE_GRID_VERSION ||
        error("checkpoint predates exact retained O2 core nodes; restart with FORCE=1")
    get(saved, "raman_shoulder_cm", nothing) == SHOULDER_CM ||
        error("checkpoint Raman shoulder differs from requested width")
    get(saved, "absco_version", "") == RRSXCO2Common.ABSCO_VERSION ||
        error("checkpoint spectroscopy differs from ABSCO $(RRSXCO2Common.ABSCO_VERSION)")
    return Set{String}(saved["completed"])
end

scene_path(state) = joinpath(AEROSOL_OUT, @sprintf("hiressim_%03d.nc", state.index))

function initialize_scene!(state, grids)
    path = scene_path(state)
    if isfile(path) && !FORCE
        return path
    end
    NCDataset(path, "c") do ds
        RRSXCO2Common.write_absco_provenance!(ds.attrib)
        RRSXCO2Common.write_fourier_convergence_provenance!(ds.attrib)
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
        ds.attrib["sza_deg"] = SZA_DEG
        ds.attrib["vza_deg"] = VZA_DEG
        ds.attrib["relative_azimuth_deg"] = RELATIVE_AZIMUTH_DEG
        ds.attrib["strong_co2_short_shoulder_points"] = Int32(8)
        ds.attrib["strong_co2_convolution_support_sigma"] = 6.0
        ds.attrib["atmospheric_layers"] = Int32(NLAYERS)
        ds.attrib["aod550_sulfate"] = state.aod550[1]
        ds.attrib["aod550_organic_carbon"] = state.aod550[2]
        ds.attrib["aod550_stratospheric"] = state.aod550[3]
        RRSXCO2Common.write_sif_provenance!(
            ds.attrib, state.sif_angular_integral760 > 0)
        ds.attrib["source_state_table"] = "true_states.dat"
        ds.attrib["spectral_chunking"] =
            "O2 cores carry ±$(SHOULDER_CM) cm-1 Raman shoulders; CO2 has no shoulders"
        ds.attrib["o2_core_grid_version"] = Int32(O2_CORE_GRID_VERSION)
        ds.attrib["o2_core_grid_construction"] =
            "retained output nodes inserted verbatim between independently built shoulders"
        ds.attrib["surface_coordinate"] =
            "Legendre x is defined once over each complete output band before chunking"
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
            # Surface normalization is defined by the complete output band;
            # only its already-defined spectrum is partitioned into chunks.
            result = simulate_o2(representative, (output_ν, solveν, keep), solar_T)
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
            result = simulate_co2(representative, iband, output_ν, ν, solar_T)
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
    selected = selected_states()
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
            # Also stamp files created before Fourier convergence was enabled
            # when a checkpointed calculation is resumed.
            RRSXCO2Common.write_fourier_convergence_provenance!(ds.attrib)
            ds.attrib["simulation_complete"] = Int32(1)
            ds.attrib["chunked_simulation_complete"] = Int32(1)
            ds.attrib["completed"] = string(now())
        end
    end
    @info "completed chunked aerosol truth map" scenes=length(states) output=AEROSOL_OUT
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main_aerosol_chunked()
