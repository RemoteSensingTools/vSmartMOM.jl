#!/usr/bin/env julia

"""Validate selected diagonal OCO-2 noise covariance products.

Set `EXPECTED_STATES` to a comma-separated list of indices and/or inclusive
ranges (for example `1-5,21-25`) when validating a campaign subset. The
historical default remains `1-64`.
"""

using NCDatasets
using Printf

include(joinpath(@__DIR__, "SyntheticOCO2.jl"))
include(joinpath(@__DIR__, "OCO2Noise.jl"))
include(joinpath(@__DIR__, "..", "RetrievalCases.jl"))
using .SyntheticOCO2
using .OCO2Noise
using .RetrievalCases: validated_sif_provenance,
                       matching_sif_provenance

const RRS_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const DEFAULT_SYNTHETIC_DIR = joinpath(RRS_ROOT, "truth_map", "OCO_radiances")
const DEFAULT_NOISE_DIR = joinpath(DEFAULT_SYNTHETIC_DIR, "noise_covariances")
const DEFAULT_COEFFICIENTS = joinpath(@__DIR__, "representative_snr_coefficients.nc")
const CLASSES = (:corrected, :uncorrected)
const BOTTOM_LAYER_PROVENANCE_ATTRIBUTES = (
    "campaign", "source_state_table", "column_xco2_ppm", "background_co2_ppm",
    "bottom_co2_ppm", "bottom_co2_index", "bottom_co2_layer_index",
    "co2_profile_order", "co2_profile_definition", "co2_profile_ppm",
    "bottom_layer_dry_column_fraction", "state_table_sha256",
    "full_column_source_state", "full_column_source_scene",
    "full_column_source_sha256", "full_column_state_table_sha256",
    "bottom_co2_truth_mode", "co2_truth_reuse_source", "producer_script",
    "producer_script_sha256", "bottom_layer_truth_complete",
)

function expected_states(value=get(ENV, "EXPECTED_STATES", "1-64"))
    result = Int[]
    for token in split(value, ',')
        fields = split(strip(token), '-'; limit=2)
        if length(fields) == 1
            push!(result, parse(Int, only(fields)))
        else
            append!(result, parse(Int, fields[1]):parse(Int, fields[2]))
        end
    end
    isempty(result) && error("EXPECTED_STATES selected no states")
    length(unique(result)) == length(result) || error(
        "EXPECTED_STATES contains duplicate indices")
    return sort!(result)
end

function finite_vector(dataset, name)
    haskey(dataset, name) || error("missing $name")
    values = Float64.(nomissing(dataset[name][:], NaN))
    all(isfinite, values) || error("$name contains non-finite values")
    return values
end

function validate_file(path, synthetic_path, coefficients, expected_state)
    band_lengths = collect(length.(synthetic_grid.(BAND_SPECS)))
    band_end = cumsum(band_lengths)
    band_start = vcat(1, band_end[1:end-1] .+ 1)
    total_below = Dict(class => 0 for class in CLASSES)
    total_above = Dict(class => 0 for class in CLASSES)
    snr_extrema = Dict(class => (Inf, -Inf) for class in CLASSES)

    NCDataset(synthetic_path) do synthetic
        NCDataset(path) do noise
            get(noise.attrib, "noise_covariance_complete", 0) == 1 ||
                error("$path is not marked noise_covariance_complete=1")
            Int(noise.attrib["state_index"]) == expected_state ||
                error("state metadata mismatch in $path")
            synthetic_sif = Symbol(get(synthetic.attrib, "sif_case", "off"))
            noise_sif = Symbol(get(noise.attrib, "sif_case", "off"))
            synthetic_provenance = validated_sif_provenance(
                synthetic.attrib, synthetic_sif; source=synthetic_path)
            noise_provenance = validated_sif_provenance(
                noise.attrib, noise_sif; source=path)
            matching_sif_provenance(
                synthetic_provenance, noise_provenance;
                first_source=synthetic_path, second_source=path)
            if get(synthetic.attrib, "campaign", "") == "bottom_layer_XCO2"
                for key in BOTTOM_LAYER_PROVENANCE_ATTRIBUTES
                    haskey(noise.attrib, key) || error(
                        "$path is missing bottom-layer provenance attribute $key")
                    haskey(synthetic.attrib, key) || error(
                        "$synthetic_path is missing bottom-layer provenance attribute $key")
                    noise.attrib[key] == synthetic.attrib[key] || error(
                        "$key differs between $path and its synthetic-radiance source")
                end
            end
            Int.(noise["band_start_index"][:]) == band_start ||
                error("wrong band starts in $path")
            Int.(noise["band_end_index"][:]) == band_end ||
                error("wrong band ends in $path")

            for class in CLASSES
                reconstructed_measurement = Float64[]
                reconstructed_std = Float64[]
                reconstructed_variance = Float64[]
                reconstructed_snr = Float64[]
                reconstructed_below = Bool[]
                reconstructed_above = Bool[]
                for spec in BAND_SPECS
                    name = String(spec.name)
                    wavelength = finite_vector(synthetic, "$(name)_wavelength")
                    radiance = finite_vector(
                        synthetic, "I_OCO_$(class)_$(name)")
                    coefficient = coefficients[spec.name]
                    statistics = noise_statistics(
                        wavelength, radiance,
                        coefficient.c_photon, coefficient.c_background,
                        noise_band_spec(spec.name))
                    append!(reconstructed_measurement, radiance)
                    append!(reconstructed_std, statistics.nen_energy)
                    append!(reconstructed_variance, statistics.variance_energy)
                    append!(reconstructed_snr, statistics.snr)
                    append!(reconstructed_below, statistics.below_min_ms)
                    append!(reconstructed_above, statistics.above_max_ms)
                end

                measurement = finite_vector(noise, "measurement_$(class)")
                std = finite_vector(noise, "noise_std_$(class)")
                variance = finite_vector(noise, "Se_diagonal_$(class)")
                snr = finite_vector(noise, "snr_$(class)")
                measurement == reconstructed_measurement ||
                    error("stored $class measurement mismatch in $path")
                std == reconstructed_std ||
                    error("stored $class NEN mismatch in $path")
                variance == reconstructed_variance ||
                    error("stored $class covariance mismatch in $path")
                snr == reconstructed_snr ||
                    error("stored $class SNR mismatch in $path")
                variance == std .^ 2 ||
                    error("$class covariance is not NEN^2 in $path")
                all(variance .> 0) || error("non-positive covariance in $path")

                below = Bool.(noise["below_min_ms_$(class)"][:])
                above = Bool.(noise["above_max_ms_$(class)"][:])
                below == reconstructed_below ||
                    error("stored $class MinMS flags mismatch in $path")
                above == reconstructed_above ||
                    error("stored $class MaxMS flags mismatch in $path")
                total_below[class] += count(below)
                total_above[class] += count(above)
                current_min, current_max = snr_extrema[class]
                snr_extrema[class] = (min(current_min, minimum(snr)),
                                      max(current_max, maximum(snr)))
            end

            weak_strong = (band_start[2]:band_end[3])
            finite_vector(noise, "Se_diagonal_corrected")[weak_strong] ==
                finite_vector(noise, "Se_diagonal_uncorrected")[weak_strong] ||
                error("corrected/uncorrected CO2 covariance differs in $path")
        end
    end
    return total_below, total_above, snr_extrema
end

function main()
    synthetic_directory = get(
        ENV, "SYNTHETIC_OCO_DIR", DEFAULT_SYNTHETIC_DIR)
    noise_directory = get(ENV, "OCO_NOISE_OUT", DEFAULT_NOISE_DIR)
    coefficient_path = get(
        ENV, "OCO_SNR_COEFFICIENTS", DEFAULT_COEFFICIENTS)
    coefficients = read_representative_snr_coefficients(coefficient_path)
    global_below = Dict(class => 0 for class in CLASSES)
    global_above = Dict(class => 0 for class in CLASSES)
    global_snr = Dict(class => (Inf, -Inf) for class in CLASSES)

    states = expected_states()
    for state in states
        noise_path = joinpath(noise_directory, @sprintf("OCO2noise_%03d.nc", state))
        synthetic_path = joinpath(
            synthetic_directory, @sprintf("OCO2sims_%03d.nc", state))
        isfile(noise_path) || error("missing $noise_path")
        isfile(synthetic_path) || error("missing $synthetic_path")
        below, above, extrema_by_class = validate_file(
            noise_path, synthetic_path, coefficients, state)
        for class in CLASSES
            global_below[class] += below[class]
            global_above[class] += above[class]
            current_min, current_max = global_snr[class]
            file_min, file_max = extrema_by_class[class]
            global_snr[class] = (min(current_min, file_min),
                                 max(current_max, file_max))
        end
    end

    @info "validated diagonal OCO-2 noise covariances" directory=noise_directory files=length(states) states
    for class in CLASSES
        @printf("%-11s SNR range %.8g to %.8g; below MinMS=%d; above MaxMS=%d\n",
                class, global_snr[class]..., global_below[class],
                global_above[class])
    end
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
