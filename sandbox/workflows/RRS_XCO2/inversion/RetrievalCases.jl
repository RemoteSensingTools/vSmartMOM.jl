module RetrievalCases

using NCDatasets
using Printf
using Random

export TruthCase,
       RetrievalExperiment,
       MeasurementRealization,
       default_truth_table,
       default_measurement_directory,
       default_noise_directory,
       read_no_sif_truth_cases,
       build_experiments,
       paired_uniform_draw,
       load_measurement_realization,
       write_experiment_manifest

const WORKFLOW_ROOT = normpath(joinpath(@__DIR__, ".."))
const RRS_ROOT = normpath(get(ENV, "RRS_XCO2_DATA_ROOT", WORKFLOW_ROOT))
const INVERSION_ROOT = joinpath(RRS_ROOT, "inversion")
const DEFAULT_MASTER_SEED = UInt64(0x4f434f5f525253)
const MEASUREMENT_CLASSES = (:corrected, :uncorrected)

default_truth_table() = joinpath(RRS_ROOT, "truth_map", "true_states.dat")
default_measurement_directory() = joinpath(
    RRS_ROOT, "truth_map", "OCO_radiances")
default_noise_directory() = joinpath(
    default_measurement_directory(), "noise_covariances")

"""The identifying metadata for one selected no-SIF truth scene."""
struct TruthCase
    state_index::Int
    surface_index::Int
    surface::Symbol
    aerosol_index::Int
    aerosol_case::Symbol
    xco2_index::Int
    xco2_ppm::Float64
end

"""One of 640 retrieval solves (320 noise draws times two measurement classes)."""
struct RetrievalExperiment
    retrieval_index::Int
    pair_index::Int
    truth::TruthCase
    noise_index::Int
    measurement_class::Symbol
    random_seed::UInt64
    measurement_path::String
    noise_path::String
end

"""Frozen measurement, covariance, and paired perturbation for one experiment."""
struct MeasurementRealization
    noiseless::Vector{Float64}
    perturbed::Vector{Float64}
    noise_std::Vector{Float64}
    variance::Vector{Float64}
    normalized_draw::Vector{Float64}
    wavelength_nm::Vector{Float64}
    band_ranges::Vector{UnitRange{Int}}
end

function _truth_rows(path)
    isfile(path) || throw(ArgumentError("missing truth-state table: $path"))
    names = String[]
    rows = Vector{Dict{String,String}}()
    for line in eachline(path)
        stripped = strip(line)
        isempty(stripped) && continue
        if startswith(stripped, "# index ")
            names = split(stripped[3:end])
        elseif !startswith(stripped, '#')
            isempty(names) && error("truth-state column header was not found in $path")
            values = split(stripped)
            length(values) == length(names) || error(
                "truth-state row has $(length(values)) fields; expected $(length(names))")
            push!(rows, Dict(zip(names, values)))
        end
    end
    isempty(rows) && error("truth-state table contains no data rows: $path")
    return rows
end

"""Read and validate the 4 surface x 2 aerosol x 4 CO2 no-SIF truth subset."""
function read_no_sif_truth_cases(path::AbstractString=default_truth_table())
    selected = TruthCase[]
    for row in _truth_rows(path)
        row["sif_case"] == "off" || continue
        push!(selected, TruthCase(
            parse(Int, row["index"]),
            parse(Int, row["surface_index"]),
            Symbol(row["surface"]),
            parse(Int, row["aerosol_index"]),
            Symbol(row["aerosol_case"]),
            parse(Int, row["xco2_index"]),
            parse(Float64, row["xco2_ppm"]),
        ))
    end
    sort!(selected; by=case -> (
        case.surface_index, case.aerosol_index, case.xco2_index))
    length(selected) == 32 || error(
        "expected 32 no-SIF truth cases; found $(length(selected))")
    length(unique(case.surface for case in selected)) == 4 || error(
        "no-SIF truth subset does not contain four surfaces")
    length(unique(case.aerosol_index for case in selected)) == 2 || error(
        "no-SIF truth subset does not contain aerosol/no-aerosol cases")
    length(unique(case.xco2_ppm for case in selected)) == 4 || error(
        "no-SIF truth subset does not contain four CO2 values")
    return selected
end

@inline function _pair_seed(master_seed::UInt64, state_index::Int, noise_index::Int)
    # Fixed integer arithmetic avoids Julia's session-randomized hash and gives
    # the corrected/uncorrected pair exactly the same standardized noise draw.
    return master_seed + UInt64(10_000state_index + noise_index)
end

"""Build the canonical 640-solve experiment sequence with adjacent pairs."""
function build_experiments(cases=read_no_sif_truth_cases();
                           noise_realizations::Int=10,
                           master_seed::UInt64=DEFAULT_MASTER_SEED,
                           measurement_directory::AbstractString=
                               default_measurement_directory(),
                           noise_directory::AbstractString=
                               default_noise_directory())
    noise_realizations > 0 || throw(ArgumentError(
        "noise_realizations must be positive"))
    experiments = RetrievalExperiment[]
    pair_index = 0
    retrieval_index = 0
    for truth in cases, noise_index in 1:noise_realizations
        pair_index += 1
        seed = _pair_seed(master_seed, truth.state_index, noise_index)
        measurement_path = joinpath(
            measurement_directory, @sprintf("OCO2sims_%03d.nc", truth.state_index))
        noise_path = joinpath(
            noise_directory, @sprintf("OCO2noise_%03d.nc", truth.state_index))
        isfile(measurement_path) || error("missing measurement file: $measurement_path")
        isfile(noise_path) || error("missing noise file: $noise_path")
        for measurement_class in MEASUREMENT_CLASSES
            retrieval_index += 1
            push!(experiments, RetrievalExperiment(
                retrieval_index, pair_index, truth, noise_index,
                measurement_class, seed, measurement_path, noise_path))
        end
    end
    expected = length(cases) * noise_realizations * length(MEASUREMENT_CLASSES)
    length(experiments) == expected || error("experiment-count construction failed")
    return experiments
end

"""Generate a reproducible unit-variance uniform draw in [-sqrt(3),sqrt(3)]."""
function paired_uniform_draw(seed::UInt64, count::Integer)
    count > 0 || throw(ArgumentError("draw length must be positive"))
    rng = Xoshiro(seed)
    return (2sqrt(3.0)) .* rand(rng, count) .- sqrt(3.0)
end

function _finite_vector(dataset, name)
    haskey(dataset, name) || error("NetCDF input is missing $name")
    values = Float64.(nomissing(dataset[name][:], NaN))
    all(isfinite, values) || error("$name contains non-finite or missing values")
    return values
end

"""Load one frozen `S_epsilon` product and apply its deterministic pair draw."""
function load_measurement_realization(experiment::RetrievalExperiment)
    label = String(experiment.measurement_class)
    return NCDataset(experiment.noise_path) do dataset
        get(dataset.attrib, "noise_covariance_complete", 0) == 1 || error(
            "noise file is not marked complete: $(experiment.noise_path)")
        Int(dataset.attrib["state_index"]) == experiment.truth.state_index || error(
            "noise-file state metadata disagrees with experiment")
        noiseless = _finite_vector(dataset, "measurement_$label")
        noise_std = _finite_vector(dataset, "noise_std_$label")
        variance = _finite_vector(dataset, "Se_diagonal_$label")
        variance ≈ noise_std .^ 2 || error(
            "stored variance is inconsistent with noise standard deviation")
        all(>(0), variance) || error("measurement variances must be positive")
        wavelength = _finite_vector(dataset, "wavelength")
        starts = Int.(dataset["band_start_index"][:])
        stops = Int.(dataset["band_end_index"][:])
        ranges = UnitRange{Int}[first:last for (first, last) in zip(starts, stops)]
        reduce(vcat, collect.(ranges)) == collect(eachindex(noiseless)) || error(
            "stored band ranges do not partition the measurement vector")
        draw = paired_uniform_draw(experiment.random_seed, length(noiseless))
        perturbed = noiseless + noise_std .* draw
        return MeasurementRealization(
            noiseless, perturbed, noise_std, variance, draw, wavelength, ranges)
    end
end

"""Write a human-readable, whitespace-separated catalog of all solves."""
function write_experiment_manifest(experiments;
                                   output_path=joinpath(
                                       INVERSION_ROOT, "retrieval_manifest.dat"))
    mkpath(dirname(output_path))
    open(output_path, "w") do io
        println(io, "# RRS-XCO2 first retrieval suite: paired corrected/uncorrected solves")
        println(io, "# Noise: sigma_i*u_i, u_i uniform on [-sqrt(3),sqrt(3)]; " *
                    "the same u is used within each pair.")
        println(io, "# retrieval pair state surface aerosol xco2_ppm perturbation class " *
                    "seed output_file measurement_file noise_file")
        for experiment in experiments
            truth = experiment.truth
            output_file = joinpath(
                INVERSION_ROOT, String(experiment.measurement_class),
                @sprintf("retrieval_state%03d_perturbation%02d.nc",
                         truth.state_index, experiment.noise_index))
            @printf(io, "%04d %03d %03d %-7s %-12s %6.1f %02d %-11s %020u %s %s %s\n",
                    experiment.retrieval_index,
                    experiment.pair_index,
                    truth.state_index,
                    String(truth.surface),
                    String(truth.aerosol_case),
                    truth.xco2_ppm,
                    experiment.noise_index,
                    String(experiment.measurement_class),
                    experiment.random_seed,
                    output_file,
                    experiment.measurement_path,
                    experiment.noise_path)
        end
    end
    return output_path
end

end # module RetrievalCases
