#!/usr/bin/env julia

"""Compare reference and optimized Round-5 forward/Jacobian benchmark products."""

using JLD2
using LinearAlgebra
using NCDatasets
using Printf
using Statistics

required_file(name) = begin
    path = get(ENV, name, "")
    !isempty(path) && isfile(path) || error("$name must identify a file")
    realpath(path)
end

function load_noise(path::AbstractString)
    return NCDataset(path, "r") do dataset
        Float64.(dataset["noise_std_corrected"][:])
    end
end

function relative_l2(delta, reference)
    denominator = norm(reference)
    return iszero(denominator) ? norm(delta) : norm(delta) / denominator
end

function report_line(io, text="")
    println(io, text)
    println(text)
end

function main()
    legacy_path = required_file("ROUND5_BENCH_LEGACY")
    optimized_path = required_file("ROUND5_BENCH_OPTIMIZED")
    noise_path = required_file("ROUND5_BENCH_NOISE")
    report_path = abspath(get(
        ENV, "ROUND5_BENCH_REPORT", "round5_migration_benchmark.md"))

    legacy = load(legacy_path)
    optimized = load(optimized_path)
    legacy["implementation"] in ("legacy", "candidate_reference") || error(
        "reference product mislabeled")
    optimized["implementation"] == "optimized" || error(
        "optimized product mislabeled")

    identity_keys = (
        "seed", "state_index", "surface", "aerosol_case", "sif_case",
        "bottom_co2_ppm", "column_xco2_ppm", "fixed_upper_co2_ppm",
        "truth_table_sha256", "prior_sha256", "coefficient_sha256",
        "component_sha256", "config_sha256", "active_parameter_names",
        "active_to_full", "band_starts", "band_stops", "active_state")
    for key in identity_keys
        legacy[key] == optimized[key] || error(
            "benchmark identity mismatch for $key")
    end

    y_legacy = Float64.(legacy["measurement"])
    y_optimized = Float64.(optimized["measurement"])
    K_legacy = Float64.(legacy["jacobian"])
    K_optimized = Float64.(optimized["jacobian"])
    size(y_legacy) == size(y_optimized) || error("measurement shapes differ")
    size(K_legacy) == size(K_optimized) || error("Jacobian shapes differ")
    size(K_legacy, 1) == length(y_legacy) || error("invalid Jacobian row count")
    noise = load_noise(noise_path)
    length(noise) == length(y_legacy) || error("noise vector length differs")
    all(>(0), noise) || error("noise standard deviation must be positive")

    delta_y = y_optimized .- y_legacy
    delta_K = K_optimized .- K_legacy
    forward_max_abs = maximum(abs, delta_y)
    forward_relative_l2 = relative_l2(delta_y, y_legacy)
    forward_max_noise = maximum(abs.(delta_y) ./ noise)

    names = String.(legacy["active_parameter_names"])
    jacobian_relative_l2 = [
        relative_l2(@view(delta_K[:, column]), @view(K_legacy[:, column]))
        for column in axes(K_legacy, 2)]
    jacobian_noise_weighted_relative_l2 = [
        relative_l2(@view(delta_K[:, column]) ./ noise,
                    @view(K_legacy[:, column]) ./ noise)
        for column in axes(K_legacy, 2)]
    worst_column = argmax(jacobian_relative_l2)
    worst_noise_weighted_column = argmax(jacobian_noise_weighted_relative_l2)

    forward_gate_noise = 0.05
    jacobian_gate_relative_l2 = 0.005
    forward_pass = forward_max_noise <= forward_gate_noise
    jacobian_pass = jacobian_relative_l2[worst_column] <=
        jacobian_gate_relative_l2
    parity_pass = forward_pass && jacobian_pass

    legacy_time = median(Float64.(legacy["wall_seconds"]))
    optimized_time = median(Float64.(optimized["wall_seconds"]))
    speedup = legacy_time / optimized_time
    starts = Int.(legacy["band_starts"])
    stops = Int.(legacy["band_stops"])
    band_names = ("O2 A", "weak CO2", "strong CO2")
    seed = legacy["seed"]
    state_index = legacy["state_index"]
    surface = legacy["surface"]
    aerosol_case = legacy["aerosol_case"]
    sif_case = legacy["sif_case"]
    bottom_co2_ppm = legacy["bottom_co2_ppm"]
    legacy_checkpoint = legacy["checkpoint"]
    optimized_checkpoint = optimized["checkpoint"]
    reference_implementation = legacy["implementation"]
    acceptance = parity_pass ? "PASS" : "FAIL"

    mkpath(dirname(report_path))
    open(report_path, "w") do io
        report_line(io, "# Round-5 reference versus optimized closure and timing")
        report_line(io)
        report_line(io, "- Random seed: `$seed`")
        report_line(io, "- Selected state: `$(lpad(state_index, 3, '0'))`")
        report_line(io, "- Scene: `$surface`, `$aerosol_case`, `$sif_case`, " *
                        "bottom-layer CO2 `$bottom_co2_ppm ppm`")
        report_line(io, "- Reference implementation: `$reference_implementation`")
        report_line(io, "- Reference commit: `$legacy_checkpoint`")
        report_line(io, "- Optimized base commit: `$optimized_checkpoint`")
        report_line(io, "- Compared output: complete Mueller-processed, convolved, " *
                        "resampled three-band forward vector and all 28 analytical " *
                        "Round-5 Jacobian columns.")
        report_line(io)
        report_line(io, "## Numerical closure")
        report_line(io)
        report_line(io, @sprintf(
            "- Forward max absolute difference: `%.9e`", forward_max_abs))
        report_line(io, @sprintf(
            "- Forward relative L2 difference: `%.9e`", forward_relative_l2))
        report_line(io, @sprintf(
            "- Forward max difference/noise sigma: `%.9e` (gate `%.3g`)",
            forward_max_noise, forward_gate_noise))
        report_line(io, @sprintf(
            "- Worst Jacobian-column relative L2: `%.9e` for `%s` (gate `%.3g`)",
            jacobian_relative_l2[worst_column], names[worst_column],
            jacobian_gate_relative_l2))
        report_line(io, @sprintf(
            "- Worst noise-weighted Jacobian relative L2: `%.9e` for `%s`",
            jacobian_noise_weighted_relative_l2[worst_noise_weighted_column],
            names[worst_noise_weighted_column]))
        report_line(io, "- Acceptance: **$acceptance**")
        report_line(io)
        report_line(io, "| Band | max abs forward difference | max difference/noise | relative L2 |")
        report_line(io, "|---|---:|---:|---:|")
        for (band, first_index, last_index) in zip(band_names, starts, stops)
            indices = first_index:last_index
            report_line(io, @sprintf(
                "| %s | %.9e | %.9e | %.9e |", band,
                maximum(abs, @view(delta_y[indices])),
                maximum(abs.(@view(delta_y[indices])) ./ @view(noise[indices])),
                relative_l2(@view(delta_y[indices]), @view(y_legacy[indices]))))
        end
        report_line(io)
        report_line(io, "## Warmed timing")
        report_line(io)
        report_line(io, @sprintf("- Reference median: `%.6f s`", legacy_time))
        report_line(io, @sprintf("- Optimized median: `%.6f s`", optimized_time))
        report_line(io, @sprintf("- Speedup: **`%.3fx`**", speedup))
        report_line(io)
        report_line(io, "Compilation and the first warm-up evaluation are excluded. " *
                        "Each timed evaluation includes parameter copying, model " *
                        "construction, all three linearized RT bands, Mueller " *
                        "processing, convolution, and resampling.")
    end
    parity_pass || error("Round-5 migration parity gate failed; see $report_path")
    println("Report saved to $report_path")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
