#!/usr/bin/env julia

"""
Merge the separately simulated strong-CO2 convolution shoulder into the
canonical high-resolution truth files.

Each truth scene's `strong_co2` dimension is expanded from 987 to 995 samples
by appending the eight short-wavelength (high-wavenumber) samples.  The local
`sim_wavelength.nc` files beside the no-aerosol and aerosol scene collections
are rewritten consistently.  Files are written to a temporary path, validated,
and atomically renamed over their target only after the new file is complete.

The source shoulder directory is deliberately retained.  Delete it only after
running an independent full-dataset validation.

Environment variables:

- `TRUTH_ROOT` (default `RRS_XCO2/truth_map`)
- `AEROSOL_TRUTH` (default `\$TRUTH_ROOT/aerosol_chunked`)
- `CONVOLUTION_SHOULDER_DIR` (default `\$TRUTH_ROOT/convolution_shoulders`)
- `MERGE_FIRST_STATE`, `MERGE_LAST_STATE` (defaults 1 and 64)
- `MERGE_WAVELENGTH_FILES=0` skips the two shared wavelength files (testing)
"""

using Dates
using NCDatasets
using Printf

const SCRIPT_ROOT = normpath(joinpath(@__DIR__, ".."))
const TRUTH_ROOT = get(ENV, "TRUTH_ROOT", joinpath(SCRIPT_ROOT, "truth_map"))
const AEROSOL_TRUTH = get(
    ENV, "AEROSOL_TRUTH", joinpath(TRUTH_ROOT, "aerosol_chunked"))
const SHOULDER_DIR = get(
    ENV, "CONVOLUTION_SHOULDER_DIR",
    joinpath(TRUTH_ROOT, "convolution_shoulders"))
const FIRST_STATE = parse(Int, get(ENV, "MERGE_FIRST_STATE", "1"))
const LAST_STATE = parse(Int, get(ENV, "MERGE_LAST_STATE", "64"))
const MERGE_WAVELENGTH_FILES =
    lowercase(get(ENV, "MERGE_WAVELENGTH_FILES", "1")) in
    ("1", "true", "yes", "on")

const SHOULDER_STEM = "strong_co2_short_shoulder"
const SHOULDER_RADIANCE = "radiance_rayleigh_$(SHOULDER_STEM)"
const STRONG_RADIANCE = "radiance_rayleigh_strong_co2"
const SCENE_VARIABLES = (
    ("radiance_rayleigh_o2a", ("stokes", "o2a")),
    ("radiance_cabannes_o2a", ("stokes", "o2a")),
    ("radiance_rrs_o2a", ("stokes", "o2a")),
    ("radiance_rayleigh_weak_co2", ("stokes", "weak_co2")),
    (STRONG_RADIANCE, ("stokes", "strong_co2")),
)

1 <= FIRST_STATE <= LAST_STATE <= 64 || error(
    "invalid merge state interval $FIRST_STATE:$LAST_STATE")

attributes(proxy) = Dict(String(key) => proxy[key] for key in keys(proxy))

function write_attributes!(proxy, values)
    for (name, value) in values
        proxy[name] = value
    end
end

function finite_array(variable, label)
    indices = ntuple(_ -> Colon(), ndims(variable))
    raw = variable[indices...]
    any(ismissing, raw) && error("$label contains missing/fill values")
    values = collect(skipmissing(raw))
    all(isfinite, values) || error("$label contains non-finite values")
    maximum(abs, values) < 1e30 || error("$label contains NetCDF fill values")
    return reshape(values, size(raw))
end

shoulder_path(index) = joinpath(
    SHOULDER_DIR, @sprintf("hires_shoulders_%03d.nc", index))

function scene_path(index)
    name = @sprintf("hiressim_%03d.nc", index)
    candidates = (joinpath(TRUTH_ROOT, name), joinpath(AEROSOL_TRUTH, name))
    existing = filter(isfile, collect(candidates))
    length(existing) == 1 || error(
        "state $index must have exactly one active truth file; found " *
        string(existing))
    return only(existing)
end

function read_shoulder(index)
    path = shoulder_path(index)
    isfile(path) || error("missing convolution shoulder: $path")
    return NCDataset(path, "r") do dataset
        get(dataset.attrib, "state_index", -1) == index || error(
            "$path has the wrong state_index")
        get(dataset.attrib, "supplemental_convolution_shoulders_complete", 0) == 1 ||
            error("$path is not marked complete")
        get(dataset.attrib, "segments", "") == "strong_co2:short:8" ||
            error("$path has unexpected segment metadata")
        wavenumber = Float64.(finite_array(
            dataset["$(SHOULDER_STEM)_wavenumber"], "$path wavenumber"))
        wavelength = Float64.(finite_array(
            dataset["$(SHOULDER_STEM)_wavelength"], "$path wavelength"))
        radiance = Float32.(finite_array(
            dataset[SHOULDER_RADIANCE], "$path radiance"))
        length(wavenumber) == 8 || error("$path does not contain eight points")
        size(radiance) == (3, 8) || error(
            "$path radiance has shape $(size(radiance)); expected (3, 8)")
        maximum(abs.(wavelength .- 1e7 ./ wavenumber)) <= 1e-10 || error(
            "$path wavelength and wavenumber coordinates disagree")
        issorted(wavenumber) || error("$path wavenumbers are not increasing")
        return (; wavenumber, wavelength, radiance)
    end
end

function validate_join(base_wavenumber, shoulder_wavenumber, label)
    length(base_wavenumber) == 987 || error(
        "$label has $(length(base_wavenumber)) base strong-band points; expected 987")
    base_wavenumber[end] < shoulder_wavenumber[1] || error(
        "$label shoulder does not follow the base grid")
    all_step = diff(vcat(base_wavenumber, shoulder_wavenumber))
    maximum(abs.(all_step .- 0.1)) <= 2e-3 || error(
        "$label is not contiguous at 0.1 cm^-1")
end

function validate_merged_scene(path, base_strong, shoulder_radiance)
    NCDataset(path, "r") do dataset
        size(dataset[STRONG_RADIANCE]) == (3, 995) || error(
            "$path merged strong-band shape is $(size(dataset[STRONG_RADIANCE]))")
        merged = Float32.(finite_array(dataset[STRONG_RADIANCE], path))
        merged[:, 1:987] == base_strong || error(
            "$path changed the original strong-band radiance")
        merged[:, 988:995] == shoulder_radiance || error(
            "$path did not append the shoulder exactly")
        get(dataset.attrib, "strong_co2_short_shoulder_merged", 0) == 1 ||
            error("$path lacks the merge completion attribute")
    end
end

function merge_scene!(index, canonical_wavenumber)
    path = scene_path(index)
    shoulder = read_shoulder(index)
    shoulder.wavenumber == canonical_wavenumber || error(
        "state $index uses a different shoulder grid")

    scene = NCDataset(path, "r") do dataset
        strong_size = size(dataset[STRONG_RADIANCE], 2)
        if strong_size == 995
            get(dataset.attrib, "strong_co2_short_shoulder_merged", 0) == 1 ||
                error("$path has 995 points but is not marked as a completed merge")
            return nothing
        end
        strong_size == 987 || error(
            "$path has an unexpected strong-band length $strong_size")
        variable_data = Dict{String,Array}()
        variable_attributes = Dict{String,Dict{String,Any}}()
        for (name, _) in SCENE_VARIABLES
            variable_data[name] = finite_array(dataset[name], "$path:$name")
            variable_attributes[name] = attributes(dataset[name].attrib)
        end
        return (
            global_attributes=attributes(dataset.attrib),
            variable_data,
            variable_attributes,
        )
    end
    isnothing(scene) && (@info "skip already merged scene" index path; return)

    base_strong = Float32.(scene.variable_data[STRONG_RADIANCE])
    merged_strong = hcat(base_strong, shoulder.radiance)
    temporary = path * ".merge_tmp"
    isfile(temporary) && rm(temporary)
    try
        NCDataset(temporary, "c") do output
            defDim(output, "stokes", 3)
            defDim(output, "o2a", 2735)
            defDim(output, "weak_co2", 1281)
            defDim(output, "strong_co2", 995)
            for (name, dimensions) in SCENE_VARIABLES
                variable = defVar(output, name, Float32, dimensions)
                write_attributes!(variable.attrib, scene.variable_attributes[name])
                data = name == STRONG_RADIANCE ? merged_strong :
                    Float32.(scene.variable_data[name])
                variable[:, :] = data
            end
            write_attributes!(output.attrib, scene.global_attributes)
            output.attrib["simulation_complete"] = Int32(1)
            output.attrib["strong_co2_short_shoulder_merged"] = Int32(1)
            output.attrib["strong_co2_short_shoulder_points"] = Int32(8)
            output.attrib["strong_co2_convolution_support_sigma"] = 6.0
            output.attrib["strong_co2_shoulder_source"] = basename(shoulder_path(index))
            output.attrib["strong_co2_shoulder_merged_utc"] = string(now(UTC))
        end
        validate_merged_scene(temporary, base_strong, shoulder.radiance)
        mv(temporary, path; force=true)
    finally
        isfile(temporary) && rm(temporary)
    end
    @info "merged strong-CO2 shoulder into scene" index path
end

function merge_wavelength_file!(path, shoulder)
    isfile(path) || error("missing wavelength file: $path")
    contents = NCDataset(path, "r") do dataset
        strong_wavenumber = Float64.(finite_array(
            dataset["strong_co2_wavenumber"], "$path strong wavenumber"))
        if length(strong_wavenumber) == 995
            get(dataset.attrib, "strong_co2_short_shoulder_merged", 0) == 1 ||
                error("$path has 995 points but lacks merge metadata")
            return nothing
        end
        validate_join(strong_wavenumber, shoulder.wavenumber, path)
        names = (
            "o2a_wavelength", "o2a_wavenumber",
            "weak_co2_wavelength", "weak_co2_wavenumber",
            "strong_co2_wavelength", "strong_co2_wavenumber",
        )
        data = Dict(name => Float64.(finite_array(dataset[name], "$path:$name"))
                    for name in names)
        variable_attributes = Dict(
            name => attributes(dataset[name].attrib) for name in names)
        return (; data, variable_attributes,
                global_attributes=attributes(dataset.attrib))
    end
    isnothing(contents) && (@info "skip already merged wavelength file" path; return)

    contents.data["strong_co2_wavelength"] = vcat(
        contents.data["strong_co2_wavelength"], shoulder.wavelength)
    contents.data["strong_co2_wavenumber"] = vcat(
        contents.data["strong_co2_wavenumber"], shoulder.wavenumber)
    temporary = path * ".merge_tmp"
    isfile(temporary) && rm(temporary)
    try
        NCDataset(temporary, "c") do output
            defDim(output, "o2a", 2735)
            defDim(output, "weak_co2", 1281)
            defDim(output, "strong_co2", 995)
            for band in ("o2a", "weak_co2", "strong_co2")
                for coordinate in ("wavelength", "wavenumber")
                    name = "$(band)_$(coordinate)"
                    variable = defVar(output, name, Float64, (band,))
                    write_attributes!(variable.attrib,
                                      contents.variable_attributes[name])
                    variable[:] = contents.data[name]
                end
            end
            write_attributes!(output.attrib, contents.global_attributes)
            output.attrib["strong_co2_short_shoulder_merged"] = Int32(1)
            output.attrib["strong_co2_short_shoulder_points"] = Int32(8)
            output.attrib["strong_co2_shoulder_merged_utc"] = string(now(UTC))
        end
        NCDataset(temporary, "r") do dataset
            merged_wavenumber = Float64.(finite_array(
                dataset["strong_co2_wavenumber"], temporary))
            merged_wavelength = Float64.(finite_array(
                dataset["strong_co2_wavelength"], temporary))
            length(merged_wavenumber) == 995 || error(
                "$temporary does not have 995 strong-band points")
            merged_wavenumber[988:995] == shoulder.wavenumber || error(
                "$temporary wavenumber shoulder differs from its source")
            merged_wavelength[988:995] == shoulder.wavelength || error(
                "$temporary wavelength shoulder differs from its source")
        end
        mv(temporary, path; force=true)
    finally
        isfile(temporary) && rm(temporary)
    end
    @info "merged strong-CO2 shoulder into wavelength file" path
end

function main()
    canonical = read_shoulder(FIRST_STATE)
    for index in FIRST_STATE:LAST_STATE
        merge_scene!(index, canonical.wavenumber)
    end
    if MERGE_WAVELENGTH_FILES
        FIRST_STATE == 1 && LAST_STATE == 64 || error(
            "shared wavelength files may be merged only for the complete 1:64 run")
        merge_wavelength_file!(joinpath(TRUTH_ROOT, "sim_wavelength.nc"), canonical)
        merge_wavelength_file!(joinpath(AEROSOL_TRUTH, "sim_wavelength.nc"), canonical)
    end
    @info "completed strong-CO2 shoulder merge" states=FIRST_STATE:LAST_STATE
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
