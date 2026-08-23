#!/usr/bin/env julia

"""Combine external-solar CPU shards into a plotting-friendly ASCII table."""

using JLD2
using Printf

const ROOT = normpath(joinpath(@__DIR__, ".."))
const INPUTS = [
    joinpath(ROOT, "truth_map_aerosols", "untruncated_convergence_cpu.jld2"),
    joinpath(ROOT, "truth_map_aerosols", "untruncated_convergence_cpu_curry.jld2"),
]
const OUTPUT = joinpath(ROOT, "truth_map_aerosols", "truncation_iqU.dat")
const CAPS = (8, 12, 16, 24, 32)
const VZAS = 0:5:70

function combined_results()
    out = Dict{String,Any}()
    for path in INPUTS
        isfile(path) || error("Missing benchmark shard: $path")
        merge!(out, load(path, "results"))
    end
    return out
end

function main()
    results = combined_results()
    open(OUTPUT, "w") do io
        println(io, "# vza_deg vaz_deg case I Q U")
        for vza in VZAS
            tag = lpad(string(vza), 2, '0')
            for case in (string.(CAPS)..., "untruncated")
                key = "band1_vza$(tag)_" * (case == "untruncated" ? case : "cap$(case)")
                haskey(results, key) || error("Missing result $key")
                r = results[key]
                size(r.radiance) == (2, 3, 1) || error(
                    "Unexpected radiance shape for $key: $(size(r.radiance))")
                for iview in 1:2
                    @printf(io, "%5.1f %5.1f %12s %.9e %.9e %.9e\n",
                            r.vza[iview], r.vaz[iview], case,
                            r.radiance[iview, 1, 1], r.radiance[iview, 2, 1],
                            r.radiance[iview, 3, 1])
                end
            end
        end
    end
    println(OUTPUT)
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
