#!/usr/bin/env julia

"""
Add the standard XCO2 diagnostics to completed retrieval products written by
an older retrieval process. No radiative-transfer calculation is repeated.

With no arguments, all completed `retrieval_state*_perturbation*.nc` products
under `corrected/` and `uncorrected/` are considered. Explicit file paths may
instead be supplied on the command line. Products that already contain
`XCO2` are left unchanged.
"""

using NCDatasets

include(joinpath(@__DIR__, "OptimalEstimation.jl"))
include(joinpath(@__DIR__, "RetrievalCases.jl"))
include(joinpath(@__DIR__, "VSmartMOMForward.jl"))
include(joinpath(@__DIR__, "RetrievalOutput.jl"))
using .VSmartMOMForward
using .RetrievalOutput

function default_products()
    paths = String[]
    for class in ("corrected", "uncorrected")
        directory = joinpath(@__DIR__, class)
        isdir(directory) || continue
        append!(paths, filter(
            path -> occursin(r"retrieval_state\d+_perturbation\d+\.nc$", path),
            readdir(directory; join=true)))
    end
    return sort!(paths)
end

function backfill!(evaluator::OCOForwardEvaluator, path::AbstractString)
    isfile(path) || throw(ArgumentError("retrieval product does not exist: $path"))
    return NCDataset(path, "a") do dataset
        get(dataset.attrib, "retrieval_complete", 0) == 1 || begin
            @info "skip incomplete retrieval" path
            return false
        end
        haskey(dataset, "XCO2") && begin
            @info "skip existing XCO2" path
            return false
        end

        fixed_upper_ppm = Float64(dataset.attrib["fixed_upper_co2_ppm"])
        set_fixed_upper_co2_ppm!(evaluator, fixed_upper_ppm)
        prior = Float64.(dataset["a_priori_state"][:])
        final = Float64.(dataset["final_state"][:])
        states = Float64.(dataset["state_at_trial"][:, :])
        ntrial = size(states, 2)
        diagnostics = (
            a_priori_ppm=column_averaged_co2_ppm(evaluator, prior),
            trial_ppm=[column_averaged_co2_ppm(evaluator, @view(states[:, trial]))
                       for trial in 1:ntrial],
            final_ppm=column_averaged_co2_ppm(evaluator, final),
        )
        write_xco2_diagnostics!(dataset, diagnostics, ntrial)
        dataset.attrib["xco2_diagnostic"] =
            "dry-air-column-weighted VMR including scene-fixed upper layers"
        @info "added XCO2" path XCO2_ppm=diagnostics.final_ppm
        return true
    end
end

function main(args=ARGS)
    paths = isempty(args) ? default_products() : abspath.(args)
    isempty(paths) && error("no retrieval products were found")
    evaluator = OCOForwardEvaluator(;
        architecture=:CPU, float_type=Float32, nstreams=9)
    updated = count(path -> backfill!(evaluator, path), paths)
    println("XCO2 backfill complete: updated=$updated considered=$(length(paths))")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
