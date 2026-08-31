#!/usr/bin/env julia

"""Validate every production synthetic OCO radiance file."""

using NCDatasets
using Printf

include(joinpath(@__DIR__, "SyntheticOCO2.jl"))
using .SyntheticOCO2

const WORKFLOW_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const RRS_ROOT = normpath(get(ENV, "RRS_XCO2_DATA_ROOT", WORKFLOW_ROOT))
const DEFAULT_OUTPUT = joinpath(RRS_ROOT, "truth_map", "OCO_radiances")

function output_files(directory)
    isdir(directory) || error("missing OCO radiance directory: $directory")
    files = filter(name -> occursin(r"^OCO2sims_\d{3}\.nc$", name),
                   readdir(directory))
    return sort(joinpath.(directory, files))
end

function state_from_name(path)
    match_result = match(r"OCO2sims_(\d{3})\.nc$", basename(path))
    isnothing(match_result) && error("cannot parse state index from $path")
    return parse(Int, only(match_result.captures))
end

function finite_vector(dataset, name, expected_length)
    haskey(dataset, name) || error("NetCDF output is missing $name")
    values = Float64.(nomissing(dataset[name][:], NaN))
    length(values) == expected_length || error(
        "$name has $(length(values)) samples; expected $expected_length")
    all(isfinite, values) || error("$name contains missing or non-finite values")
    maximum(abs, values) < 1e30 || error("$name contains NetCDF fill values")
    return values
end

function validate_file(path, expected_state)
    extrema_by_band = Dict{Symbol,Tuple{Float64,Float64}}()
    NCDataset(path) do dataset
        get(dataset.attrib, "instrument_processing_complete", 0) == 1 ||
            error("$path is not marked instrument_processing_complete=1")
        Int(dataset.attrib["state_index"]) == expected_state ||
            error("state metadata does not match filename for $path")
        source_path = String(dataset.attrib["source_truth_scene"])
        isfile(source_path) || error("source truth scene is missing for $path")

        for spec in BAND_SPECS
            grid = synthetic_grid(spec)
            coordinate_name = "$(spec.name)_wavelength"
            coordinate = finite_vector(dataset, coordinate_name, length(grid))
            coordinate == grid || error("$path has the wrong $(spec.name) grid")
            coordinate_variable = dataset[coordinate_name]
            coordinate_variable.attrib["sampling_interval_nm"] ==
                spec.sampling_interval_nm || error("wrong sampling interval in $path")
            coordinate_variable.attrib["gaussian_fwhm_nm"] == spec.fwhm_nm ||
                error("wrong Gaussian FWHM in $path")
            coordinate_variable.attrib["gaussian_sigma_nm"] ==
                fwhm_to_sigma(spec.fwhm_nm) || error("wrong Gaussian sigma in $path")

            components = spec.name == :o2a ?
                (:rayleigh, :cabannes, :rrs, :corrected, :uncorrected) :
                (:rayleigh, :corrected, :uncorrected)
            values = Dict(component => finite_vector(
                dataset, "I_OCO_$(component)_$(spec.name)", length(grid))
                for component in components)

            values[:corrected] == values[:rayleigh] || error(
                "corrected $(spec.name) is not identical to Rayleigh in $path")
            if spec.name == :o2a
                values[:uncorrected] == values[:cabannes] .+ values[:rrs] ||
                    error("uncorrected O2 A-band is not Cabannes + RRS in $path")
            else
                values[:uncorrected] == values[:rayleigh] || error(
                    "uncorrected $(spec.name) is not identical to Rayleigh in $path")
            end
            extrema_by_band[spec.name] = extrema(values[:uncorrected])
        end
    end
    return extrema_by_band
end

function main()
    directory = get(ENV, "SYNTHETIC_OCO_OUT", DEFAULT_OUTPUT)
    files = output_files(directory)
    states = state_from_name.(files)
    states == collect(1:64) || error(
        "expected states 001:064, found $(join(states, ','))")

    global_extrema = Dict(spec.name => (Inf, -Inf) for spec in BAND_SPECS)
    for (path, state) in zip(files, states)
        extrema_by_band = validate_file(path, state)
        for (band, (minimum_value, maximum_value)) in extrema_by_band
            current_minimum, current_maximum = global_extrema[band]
            global_extrema[band] = (
                min(current_minimum, minimum_value),
                max(current_maximum, maximum_value),
            )
        end
    end

    @info "validated synthetic OCO radiances" directory files=length(files)
    for spec in BAND_SPECS
        minimum_value, maximum_value = global_extrema[spec.name]
        @printf("%-10s uncorrected I_OCO range: %.9g to %.9g mW m-2 sr-1 nm-1\n",
                spec.name, minimum_value, maximum_value)
    end
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
