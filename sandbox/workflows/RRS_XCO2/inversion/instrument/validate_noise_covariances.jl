#!/usr/bin/env julia

"""Validate all 64 diagonal OCO-2 noise covariance products."""

using NCDatasets
using Printf

include(joinpath(@__DIR__, "SyntheticOCO2.jl"))
include(joinpath(@__DIR__, "OCO2Noise.jl"))
using .SyntheticOCO2
using .OCO2Noise

const WORKFLOW_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const RRS_ROOT = normpath(get(ENV, "RRS_XCO2_DATA_ROOT", WORKFLOW_ROOT))
const DEFAULT_SYNTHETIC_DIR = joinpath(RRS_ROOT, "truth_map", "OCO_radiances")
const DEFAULT_NOISE_DIR = joinpath(DEFAULT_SYNTHETIC_DIR, "noise_covariances")
const DEFAULT_COEFFICIENTS = joinpath(
    RRS_ROOT, "inversion", "instrument", "representative_snr_coefficients.nc")
const CLASSES = (:corrected, :uncorrected)

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

    for state in 1:64
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

    @info "validated diagonal OCO-2 noise covariances" directory=noise_directory files=64
    for class in CLASSES
        @printf("%-11s SNR range %.8g to %.8g; below MinMS=%d; above MaxMS=%d\n",
                class, global_snr[class]..., global_below[class],
                global_above[class])
    end
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
