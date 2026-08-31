#!/usr/bin/env julia

"""
Generate state-dependent diagonal OCO-2 noise covariance products.

For each `OCO2sims_NNN.nc`, the corrected and uncorrected measurement vectors
are converted from mW m^-2 sr^-1 nm^-1 to photon radiance, evaluated with
OCO-2 L1B ATBD Eq. (3-8), and converted back to energy-radiance noise. The
stored `Se_diagonal_*` vectors are NEN^2 and exactly define the requested
diagonal covariance matrices without allocating dense 2742 x 2742 arrays.

Environment variables:

- `SYNTHETIC_OCO_DIR`: input directory containing `OCO2sims_NNN.nc`
- `OCO_NOISE_OUT`: output directory
- `OCO_SNR_COEFFICIENTS`: representative SNR coefficient NetCDF
- `BANDS`: comma-separated subset of `o2a,weak_co2,strong_co2` (default all)
- `FIRST_STATE`, `LAST_STATE`: inclusive state limits (default 1 and 64)
- `FORCE=1`: replace existing noise products
"""

using Dates
using NCDatasets
using Printf
using Statistics

include(joinpath(@__DIR__, "SyntheticOCO2.jl"))
include(joinpath(@__DIR__, "OCO2Noise.jl"))
using .SyntheticOCO2
using .OCO2Noise

const WORKFLOW_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const RRS_ROOT = normpath(get(ENV, "RRS_XCO2_DATA_ROOT", WORKFLOW_ROOT))
const DEFAULT_SYNTHETIC_DIR = joinpath(RRS_ROOT, "truth_map", "OCO_radiances")
const DEFAULT_OUTPUT_DIR = joinpath(DEFAULT_SYNTHETIC_DIR, "noise_covariances")
const DEFAULT_COEFFICIENTS = joinpath(
    RRS_ROOT, "inversion", "instrument", "representative_snr_coefficients.nc")
const CLASSES = (:corrected, :uncorrected)
const RADIANCE_UNITS = "mW m-2 sr-1 nm-1"
const VARIANCE_UNITS = "(mW m-2 sr-1 nm-1)^2"

env_flag(name) = lowercase(get(ENV, name, "0")) in ("1", "true", "yes", "on")

function selected_band_specs()
    requested = Symbol.(strip.(split(
        get(ENV, "BANDS", "o2a,weak_co2,strong_co2"), ',')))
    isempty(requested) && error("BANDS selected no bands")
    length(unique(requested)) == length(requested) ||
        error("BANDS contains a duplicate band")
    return [band_spec(name) for name in requested]
end

function state_index(path)
    match_result = match(r"OCO2sims_(\d{3})\.nc$", basename(path))
    isnothing(match_result) && error("cannot parse state index from $path")
    return parse(Int, only(match_result.captures))
end

function discover_inputs(directory, first_state, last_state)
    isdir(directory) || error("missing synthetic OCO directory: $directory")
    files = sort(filter(path -> begin
        occursin(r"^OCO2sims_\d{3}\.nc$", basename(path)) || return false
        first_state <= state_index(path) <= last_state
    end, readdir(directory; join=true)); by=state_index)
    isempty(files) && error("no selected synthetic OCO files in $directory")
    return files
end

function finite_vector(dataset, name)
    haskey(dataset, name) || error("input is missing $name")
    values = Float64.(nomissing(dataset[name][:], NaN))
    all(isfinite, values) || error("$name contains missing or non-finite values")
    return values
end

function copy_scene_attributes(dataset)
    attributes = Dict{String,Any}()
    for key in ("state_index", "surface", "aerosol_case", "sif_case",
                "xco2_ppm", "psurf_hpa", "sza_deg", "vza_deg",
                "atmospheric_layers")
        haskey(dataset.attrib, key) && (attributes[key] = dataset.attrib[key])
    end
    return attributes
end

function define_vector(output, name, values, units, long_name)
    variable = defVar(output, name, Float64, ("measurement",))
    variable.attrib["units"] = units
    variable.attrib["long_name"] = long_name
    variable[:] = values
    return variable
end

function process_noise_file(input_path, output_path, coefficient_path, coefficients,
                            specs=BAND_SPECS)
    index = state_index(input_path)
    band_wavelength = Dict{Symbol,Vector{Float64}}()
    band_coefficients = Dict{Symbol,Any}()
    attributes = Dict{String,Any}()
    class_data = Dict(class => Dict{Symbol,Vector{Float64}}(
        :measurement => Float64[],
        :signal_photon => Float64[],
        :noise_std => Float64[],
        :variance => Float64[],
        :snr => Float64[],
        :below_min => Float64[],
        :above_max => Float64[],
    ) for class in CLASSES)
    band_counts = Dict(class => Dict(
        :below_min => Int[], :above_max => Int[]) for class in CLASSES)

    NCDataset(input_path) do input
        get(input.attrib, "instrument_processing_complete", 0) == 1 ||
            error("input is not marked instrument_processing_complete=1: $input_path")
        Int(input.attrib["state_index"]) == index ||
            error("state metadata disagrees with filename: $input_path")
        attributes = copy_scene_attributes(input)

        for spec in specs
            name = spec.name
            name_string = String(name)
            wavelength = finite_vector(input, "$(name_string)_wavelength")
            coefficient = coefficients[name]
            wavelength == coefficient.wavelength || error(
                "SNR coefficient grid does not match $name grid in $input_path")
            band_wavelength[name] = wavelength
            band_coefficients[name] = coefficient
            noise_spec = noise_band_spec(name)

            for class in CLASSES
                radiance = finite_vector(
                    input, "I_OCO_$(class)_$(name_string)")
                length(radiance) == length(wavelength) || error(
                    "$class $name radiance length does not match wavelength")
                statistics = noise_statistics(
                    wavelength, radiance,
                    coefficient.c_photon, coefficient.c_background,
                    noise_spec)
                append!(class_data[class][:measurement], radiance)
                append!(class_data[class][:signal_photon], statistics.signal_photon)
                append!(class_data[class][:noise_std], statistics.nen_energy)
                append!(class_data[class][:variance], statistics.variance_energy)
                append!(class_data[class][:snr], statistics.snr)
                append!(class_data[class][:below_min], Float64.(statistics.below_min_ms))
                append!(class_data[class][:above_max], Float64.(statistics.above_max_ms))
                push!(band_counts[class][:below_min], count(statistics.below_min_ms))
                push!(band_counts[class][:above_max], count(statistics.above_max_ms))
            end
        end
    end

    wavelength = vcat((band_wavelength[spec.name] for spec in specs)...)
    function canonical_band_index(spec)
        index = findfirst(candidate -> candidate.name == spec.name, BAND_SPECS)
        isnothing(index) && error("$(spec.name) is not a canonical OCO band")
        return index
    end
    band_index = vcat((fill(Int8(canonical_band_index(spec)),
                            length(band_wavelength[spec.name]))
                       for spec in specs)...)
    c_photon = vcat((band_coefficients[spec.name].c_photon
                     for spec in specs)...)
    c_background = vcat((band_coefficients[spec.name].c_background
                         for spec in specs)...)
    extrapolated_count = vcat((
        band_coefficients[spec.name].extrapolated_source_count
        for spec in specs)...)
    band_lengths = [length(band_wavelength[spec.name]) for spec in specs]
    band_end = cumsum(band_lengths)
    band_start = vcat(1, band_end[1:end-1] .+ 1)

    mkpath(dirname(output_path))
    isfile(output_path) && rm(output_path)
    NCDataset(output_path, "c") do output
        defDim(output, "measurement", length(wavelength))
        defDim(output, "band", length(specs))

        coordinate = defVar(output, "wavelength", Float64, ("measurement",))
        coordinate.attrib["units"] = "nm"
        coordinate.attrib["long_name"] =
            "concatenated synthetic OCO sample-center wavelength"
        coordinate[:] = wavelength
        band_variable = defVar(output, "band_index", Int8, ("measurement",))
        band_variable.attrib["flag_values"] = Int8[1, 2, 3]
        band_variable.attrib["flag_meanings"] = "o2a weak_co2 strong_co2"
        band_variable[:] = band_index
        start_variable = defVar(output, "band_start_index", Int32, ("band",))
        start_variable.attrib["index_convention"] = "one-based inclusive"
        start_variable[:] = Int32.(band_start)
        end_variable = defVar(output, "band_end_index", Int32, ("band",))
        end_variable.attrib["index_convention"] = "one-based inclusive"
        end_variable[:] = Int32.(band_end)

        max_ms = defVar(output, "max_ms", Float64, ("band",))
        max_ms.attrib["units"] = "photons s-1 m-2 sr-1 um-1"
        max_ms.attrib["source"] = "OCO L1B ATBD Table 3-5"
        max_ms[:] = [noise_band_spec(spec.name).max_ms for spec in specs]
        min_ms = defVar(output, "min_ms", Float64, ("band",))
        min_ms.attrib["units"] = "photons s-1 m-2 sr-1 um-1"
        min_ms.attrib["source"] = "OCO L1B ATBD Table 3-6, OCO-2 column"
        min_ms.attrib["usage"] = "dynamic-range validation only; not an Eq. 3-8 floor"
        min_ms[:] = [noise_band_spec(spec.name).min_ms for spec in specs]

        define_vector(output, "c_photon", c_photon, "1",
                      "representative L1B photon-noise coefficient")
        define_vector(output, "c_background", c_background, "1",
                      "representative L1B background-noise coefficient")
        extrapolated = defVar(
            output, "snr_coefficient_extrapolated_source_count", Int16,
            ("measurement",))
        extrapolated.attrib["long_name"] =
            "number of 32 L1B coefficient spectra using nearest-edge extrapolation"
        extrapolated[:] = extrapolated_count

        for class in CLASSES
            label = String(class)
            data = class_data[class]
            define_vector(output, "measurement_$(label)", data[:measurement],
                          RADIANCE_UNITS,
                          "$label synthetic OCO measurement vector")
            define_vector(output, "signal_photon_$(label)", data[:signal_photon],
                          "photons s-1 m-2 sr-1 um-1",
                          "$label measurement converted to photon radiance")
            define_vector(output, "noise_std_$(label)", data[:noise_std],
                          RADIANCE_UNITS,
                          "$label noise-equivalent energy radiance from Eq. 3-8")
            variance = define_vector(
                output, "Se_diagonal_$(label)", data[:variance],
                VARIANCE_UNITS,
                "diagonal of the $label measurement-noise covariance matrix")
            variance.attrib["matrix_reconstruction"] =
                "LinearAlgebra.Diagonal(Se_diagonal_$(label))"
            define_vector(output, "snr_$(label)", data[:snr], "1",
                          "$label signal-to-noise ratio")

            below = defVar(output, "below_min_ms_$(label)", Int8,
                           ("measurement",))
            below.attrib["flag_values"] = Int8[0, 1]
            below.attrib["flag_meanings"] = "within_or_above_MinMS below_MinMS"
            below[:] = Int8.(data[:below_min])
            above = defVar(output, "above_max_ms_$(label)", Int8,
                           ("measurement",))
            above.attrib["flag_values"] = Int8[0, 1]
            above.attrib["flag_meanings"] = "within_or_below_MaxMS above_MaxMS"
            above[:] = Int8.(data[:above_max])

            below_count = defVar(output, "below_min_ms_count_$(label)", Int32,
                                 ("band",))
            below_count[:] = Int32.(band_counts[class][:below_min])
            above_count = defVar(output, "above_max_ms_count_$(label)", Int32,
                                 ("band",))
            above_count[:] = Int32.(band_counts[class][:above_max])
        end

        for (key, value) in attributes
            output.attrib[key] = value
        end
        output.attrib["source_synthetic_measurement"] = abspath(input_path)
        output.attrib["representative_snr_coefficients"] = abspath(coefficient_path)
        output.attrib["band_order"] = join(string.(getfield.(specs, :name)), " ")
        output.attrib["measurement_vector_order"] =
            "o2a samples, weak_co2 samples, strong_co2 samples"
        output.attrib["measurement_classes"] = "corrected uncorrected"
        output.attrib["noise_equation"] =
            "NEN=(MaxMS/100)*sqrt(abs(100*N/MaxMS)*Cphoton^2+Cbackground^2)"
        output.attrib["noise_equation_source"] =
            "OCO-2 Level 1B ATBD v3.0 rev0 Eq. 3-8"
        output.attrib["photon_conversion"] =
            "mW/nm is numerically W/um; N_photon=L_energy/(h*c/lambda)"
        output.attrib["covariance_storage"] =
            "diagonal vectors; off-diagonal elements are exactly zero"
        output.attrib["correlation_assumption"] =
            "independent spectral-sample noise as requested; Eq. 3-8 supplies per-sample NEN"
        output.attrib["min_ms_usage"] =
            "Table 3-6 validation flag only; no clipping or noise-floor substitution"
        output.attrib["created_utc"] = string(now(UTC))
        output.attrib["noise_covariance_complete"] = 1
    end

    return [(state=index,
             class=class,
             minimum_snr=minimum(class_data[class][:snr]),
             median_snr=median(class_data[class][:snr]),
             maximum_snr=maximum(class_data[class][:snr]),
             below_min=sum(band_counts[class][:below_min]),
             above_max=sum(band_counts[class][:above_max]),
             attributes=attributes) for class in CLASSES]
end

function summarize_existing(path)
    NCDataset(path) do dataset
        get(dataset.attrib, "noise_covariance_complete", 0) == 1 ||
            error("existing output is incomplete: $path")
        index = Int(dataset.attrib["state_index"])
        attributes = copy_scene_attributes(dataset)
        return [(state=index,
                 class=class,
                 minimum_snr=minimum(finite_vector(dataset, "snr_$(class)")),
                 median_snr=median(finite_vector(dataset, "snr_$(class)")),
                 maximum_snr=maximum(finite_vector(dataset, "snr_$(class)")),
                 below_min=sum(Int.(dataset["below_min_ms_count_$(class)"][:])),
                 above_max=sum(Int.(dataset["above_max_ms_count_$(class)"][:])),
                 attributes=attributes) for class in CLASSES]
    end
end

function write_manifest(path, summaries)
    open(path, "w") do io
        println(io, "# state surface aerosol sif xco2_ppm class min_snr median_snr max_snr below_MinMS above_MaxMS")
        for summary in sort(summaries; by=value -> (value.state, value.class))
            attributes = summary.attributes
            @printf(io, "%03d %s %s %s %d %s %.10g %.10g %.10g %d %d\n",
                    summary.state,
                    get(attributes, "surface", "unknown"),
                    get(attributes, "aerosol_case", "unknown"),
                    get(attributes, "sif_case", "unknown"),
                    Int(get(attributes, "xco2_ppm", -1)),
                    summary.class,
                    summary.minimum_snr, summary.median_snr,
                    summary.maximum_snr, summary.below_min,
                    summary.above_max)
        end
    end
end

function main()
    input_directory = get(ENV, "SYNTHETIC_OCO_DIR", DEFAULT_SYNTHETIC_DIR)
    output_directory = get(ENV, "OCO_NOISE_OUT", DEFAULT_OUTPUT_DIR)
    coefficient_path = get(
        ENV, "OCO_SNR_COEFFICIENTS", DEFAULT_COEFFICIENTS)
    first_state = parse(Int, get(ENV, "FIRST_STATE", "1"))
    last_state = parse(Int, get(ENV, "LAST_STATE", "64"))
    force = env_flag("FORCE")
    specs = selected_band_specs()
    coefficients = read_representative_snr_coefficients(coefficient_path)
    inputs = discover_inputs(input_directory, first_state, last_state)
    mkpath(output_directory)

    summaries = NamedTuple[]
    for input_path in inputs
        index = state_index(input_path)
        output_path = joinpath(
            output_directory, @sprintf("OCO2noise_%03d.nc", index))
        if isfile(output_path) && !force
            @info "skip existing OCO noise covariance" index output_path
            append!(summaries, summarize_existing(output_path))
        else
            @info "generate OCO noise covariance" index input_path output_path
            append!(summaries, process_noise_file(
                input_path, output_path, coefficient_path, coefficients, specs))
        end
    end
    manifest = joinpath(output_directory, "noise_covariance_manifest.dat")
    write_manifest(manifest, summaries)
    @info "wrote OCO noise covariance products" output_directory files=length(inputs) manifest
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
