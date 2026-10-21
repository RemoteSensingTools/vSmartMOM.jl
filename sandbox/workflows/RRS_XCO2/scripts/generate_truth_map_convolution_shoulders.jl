#!/usr/bin/env julia

"""
Generate only the high-resolution edge wavelengths missing from the completed
truth map but required by the synthetic Gaussian OCO-2 convolution.

This script is intentionally separate from the main truth files: completed
spectra are never resized or rewritten. The retrieval measurement operator
combines each scene's base spectrum with its supplemental file in wavelength
space before convolution.

The required interval is derived from the Table-1 synthetic OCO grid and a
default +/-6 sigma Gaussian support. With the current 0.1 cm^-1 truth grids,
only eight extra samples on the short-wavelength side of the strong CO2 band
are required. The calculation is nevertheless derived rather than
hard-coded, so a later FWHM/grid change fails transparently or produces the
new required nodes.

Run only after the corrected aerosol truth production has completed:

    CUDA_DEVICE=1 julia --project=. \
        RRS_XCO2/scripts/generate_truth_map_convolution_shoulders.jl

Environment variables:

- `CONVOLUTION_SHOULDER_OUT` (default
  `RRS_XCO2/truth_map/convolution_shoulders`)
- `CURRENT_AEROSOL_TRUTH` (default
  `RRS_XCO2/truth_map/aerosol_chunked`)
- `GAUSSIAN_SUPPORT_SIGMA` (default 6)
- `ALLOW_BEFORE_TRUTH_COMPLETE=1` bypasses the completion gate for testing
- the standard truth generator controls (`CUDA_DEVICE`, `RRS_XCO2_FLOAT_TYPE`,
  `AEROSOL_NSTREAMS`, `FIRST_STATE`, `LAST_STATE`, and `FORCE`)
"""

const SHOULDER_SCRIPT_ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(@__DIR__, "generate_truth_map.jl"))
include(joinpath(SHOULDER_SCRIPT_ROOT, "inversion", "instrument", "SyntheticOCO2.jl"))
using .SyntheticOCO2

const CONVOLUTION_SHOULDER_OUT = get(
    ENV, "CONVOLUTION_SHOULDER_OUT",
    joinpath(SHOULDER_SCRIPT_ROOT, "truth_map", "convolution_shoulders"))
const CURRENT_AEROSOL_TRUTH = get(
    ENV, "CURRENT_AEROSOL_TRUTH",
    joinpath(SHOULDER_SCRIPT_ROOT, "truth_map", "aerosol_chunked"))
const GAUSSIAN_SUPPORT_SIGMA = parse(
    Float64, get(ENV, "GAUSSIAN_SUPPORT_SIGMA", "6"))
const ALLOW_BEFORE_TRUTH_COMPLETE =
    lowercase(get(ENV, "ALLOW_BEFORE_TRUTH_COMPLETE", "0")) in
    ("1", "true", "yes", "on")

GAUSSIAN_SUPPORT_SIGMA > 0 || error("GAUSSIAN_SUPPORT_SIGMA must be positive")

function current_truth_complete()
    paths = sort(filter(name -> occursin(r"^hiressim_\d+\.nc$", name),
                        isdir(CURRENT_AEROSOL_TRUTH) ? readdir(CURRENT_AEROSOL_TRUTH) : String[]))
    length(paths) == 32 || return false
    for name in paths
        complete = NCDataset(joinpath(CURRENT_AEROSOL_TRUTH, name)) do dataset
            haskey(dataset.attrib, "chunked_simulation_complete") &&
                Int(dataset.attrib["chunked_simulation_complete"]) == 1
        end
        complete || return false
    end
    return true
end

function extension_nodes(base_wavenumber, spec)
    extension = required_wavenumber_extensions(
        base_wavenumber, spec;
        step_cm=Float64(STEP),
        support_sigma=GAUSSIAN_SUPPORT_SIGMA)
    return (short=FT.(extension.short),
            long=FT.(extension.long),
            required_short_nm=extension.required_short_nm,
            required_long_nm=extension.required_long_nm)
end

function required_segments(grids)
    segments = NamedTuple[]
    for (band_index, spec) in pairs(BAND_SPECS)
        extension = extension_nodes(grids[band_index], spec)
        for side in (:short, :long)
            nodes = getfield(extension, side)
            isempty(nodes) && continue
            band_index == 1 && error(
                "O2 A-band convolution support is missing on the $side side; " *
                "this elastic-only supplemental script cannot synthesize RRS shoulders")
            push!(segments, (
                band_index=band_index,
                band=spec.name,
                side=side,
                wavenumber=nodes,
                required_short_nm=extension.required_short_nm,
                required_long_nm=extension.required_long_nm,
            ))
        end
    end
    return segments
end

function physical_groups(states)
    groups = Dict{Tuple{Int,Int,Int},Vector{TruthState}}()
    for state in states
        key = (state.surface_index, state.aerosol_index, state.xco2_index)
        push!(get!(groups, key, TruthState[]), state)
    end
    return groups
end

supplement_path(state) = joinpath(
    CONVOLUTION_SHOULDER_OUT, @sprintf("hires_shoulders_%03d.nc", state.index))

function write_segment!(state, segment, radiance)
    path = supplement_path(state)
    mode = isfile(path) ? "a" : "c"
    dimension = "$(segment.band)_$(segment.side)_shoulder"
    wavelength_name = "$(dimension)_wavelength"
    wavenumber_name = "$(dimension)_wavenumber"
    radiance_name = "radiance_rayleigh_$(dimension)"

    NCDataset(path, mode) do dataset
        if !haskey(dataset.dim, dimension)
            defDim(dataset, dimension, length(segment.wavenumber))
            haskey(dataset.dim, "stokes") || defDim(dataset, "stokes", 3)
            wavelength = defVar(dataset, wavelength_name, Float64, (dimension,))
            wavelength.attrib["units"] = "nm"
            wavelength[:] = 1e7 ./ Float64.(segment.wavenumber)
            wavenumber = defVar(dataset, wavenumber_name, Float64, (dimension,))
            wavenumber.attrib["units"] = "cm-1"
            wavenumber[:] = Float64.(segment.wavenumber)
            variable = defVar(dataset, radiance_name, Float32, ("stokes", dimension))
            variable.attrib["units"] = "mW m-2 sr-1 (cm-1)-1"
            variable[:, :] = Float32.(radiance)
        else
            dataset[radiance_name][:, :] = Float32.(radiance)
        end

        dataset.attrib["state_index"] = Int32(state.index)
        dataset.attrib["surface"] = state.surface
        dataset.attrib["aerosol_case"] = state.aerosol_case
        dataset.attrib["sif_case"] = state.sif_case
        dataset.attrib["xco2_ppm"] = state.xco2_ppm
        dataset.attrib["psurf_hpa"] = 1000.0
        dataset.attrib["atmospheric_layers"] = Int32(NLAYERS)
        dataset.attrib["gaussian_support_sigma"] = GAUSSIAN_SUPPORT_SIGMA
        dataset.attrib["required_wavelength_interval_nm"] =
            "$(segment.required_short_nm) $(segment.required_long_nm)"
        dataset.attrib["surface_coordinate"] =
            "same complete base-band Legendre coordinate as the main truth run"
        dataset.attrib["purpose"] =
            "supplemental edge basis for synthetic Gaussian OCO convolution"
        dataset.attrib["updated_utc"] = string(now(UTC))
    end
end

function mark_complete!(states, segments)
    description = join(("$(s.band):$(s.side):$(length(s.wavenumber))" for s in segments), ";")
    for state in states
        NCDataset(supplement_path(state), "a") do dataset
            dataset.attrib["supplemental_convolution_shoulders_complete"] = Int32(1)
            dataset.attrib["segments"] = description
            dataset.attrib["completed_utc"] = string(now(UTC))
        end
    end
end

function main_shoulders()
    ALLOW_BEFORE_TRUTH_COMPLETE || current_truth_complete() || error(
        "corrected aerosol truth map is not complete; wait for the production run " *
        "or set ALLOW_BEFORE_TRUTH_COMPLETE=1 only for a development test")
    CUDA.functional() || error("CUDA is not functional")
    CUDA.device!(DEVICE)
    mkpath(CONVOLUTION_SHOULDER_OUT)

    states = read_states()[FIRST_STATE:LAST_STATE]
    grids = output_grids()
    segments = required_segments(grids)
    isempty(segments) && (@info "truth grids already cover every Gaussian support interval"; return)

    if FORCE
        for state in states
            path = supplement_path(state)
            isfile(path) && rm(path)
        end
    end

    solar_T = solar_interpolator()
    groups = physical_groups(states)
    for segment in segments
        base_wavenumber = grids[segment.band_index]
        @info "supplemental convolution segment" band=segment.band side=segment.side npoints=length(segment.wavenumber) wavelength_nm=extrema(1e7 ./ Float64.(segment.wavenumber))
        for (key, members) in sort!(collect(groups); by=first)
            if all(state -> begin
                    path = supplement_path(state)
                    isfile(path) && NCDataset(path) do dataset
                        haskey(dataset, "radiance_rayleigh_$(segment.band)_$(segment.side)_shoulder")
                    end
                end, members) && !FORCE
                @info "skip existing supplemental group" segment.band segment.side key
                continue
            end
            representative = first(members)
            @info "simulate supplemental convolution shoulder" segment.band segment.side key states=getfield.(members, :index)
            radiance = simulate_co2(
                representative,
                segment.band_index,
                base_wavenumber,
                segment.wavenumber,
                solar_T,
            )
            for state in members
                write_segment!(state, segment, radiance)
            end
            CUDA.synchronize(); GC.gc(); CUDA.reclaim()
        end
    end
    mark_complete!(states, segments)
    @info "completed supplemental convolution shoulders" output=CONVOLUTION_SHOULDER_OUT scenes=length(states)
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main_shoulders()
