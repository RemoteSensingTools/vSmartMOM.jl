#!/usr/bin/env julia

"""
Run native-angle land-glint truth-map phases sequentially on one visible GPU.

Environment controls:

- `GLINT_AEROSOL_STREAM_INDICES`: comma-separated indices, default `6,7,8,9`.
- `GLINT_NO_AEROSOL_STREAM_INDICES`: comma-separated indices, default empty.
- `WAIT_FOR_PID`: optional process ID that must exit before the sequence starts.

The stream indices address the 5th through 9th nodes of the nine-stream
half-space Gauss-Legendre rule: 5 is 60 degrees and 9 is 10.237292 degrees.
Each child inherits the CUDA visibility and truth-map precision settings of
this process. The land-glint driver itself fixes aerosol calculations to nine
streams and 256-point O2 cores (11 chunks).
"""

using Dates

const WORKFLOW_ROOT = normpath(joinpath(@__DIR__, ".."))
const REPO = normpath(joinpath(WORKFLOW_ROOT, "..", "..", ".."))
const DRIVER = joinpath(@__DIR__, "generate_truth_map_land_glint.jl")
const DATA_ROOT = normpath(get(ENV, "RRS_XCO2_DATA_ROOT", WORKFLOW_ROOT))
const OUT_ROOT = joinpath(DATA_ROOT, "truth_map", "land_glint")
const STREAM_TAG = Dict(
    5 => "60",
    6 => "48p537727",
    7 => "36p226627",
    8 => "23p362328",
    9 => "10p237292",
)

function parse_indices(name, default)
    value = strip(get(ENV, name, default))
    isempty(value) && return Int[]
    indices = parse.(Int, strip.(split(value, ',')))
    all(haskey(STREAM_TAG, index) for index in indices) ||
        error("$name may contain only indices 5, 6, 7, 8, and 9")
    length(unique(indices)) == length(indices) ||
        error("$name contains a duplicate stream index")
    return indices
end

function process_exists(pid::Int)
    command = pipeline(ignorestatus(`kill -0 $pid`), stdout=devnull, stderr=devnull)
    return run(command).exitcode == 0
end

function wait_for_predecessor()
    haskey(ENV, "WAIT_FOR_PID") || return
    pid = parse(Int, ENV["WAIT_FOR_PID"])
    println("[$(now())] waiting for PID $pid")
    flush(stdout)
    while process_exists(pid)
        sleep(30)
    end
    println("[$(now())] predecessor PID $pid exited")
    flush(stdout)
end

function run_phase(index::Int, phase::String)
    tag = STREAM_TAG[index]
    geometry = "sza$(tag)_vza$(tag)_relaz00"
    phase_dir = phase == "aerosol" ? "aerosol_chunked" : "no_aerosol"
    output_dir = joinpath(OUT_ROOT, geometry, phase_dir)
    mkpath(output_dir)
    log_path = joinpath(output_dir, "run.log")

    child_env = copy(ENV)
    delete!(child_env, "WAIT_FOR_PID")
    child_env["GLINT_STREAM_INDEX"] = string(index)
    child_env["GLINT_PHASE"] = phase
    child_env["RRS_XCO2_FLOAT_TYPE"] = "Float32"
    child_env["NLAYERS"] = "16"
    child_env["RAMAN_SHOULDER_CM"] = "234"

    command = Cmd(`$(Base.julia_cmd()) --project=$REPO $DRIVER`;
                  dir=REPO, env=child_env)
    println("[$(now())] starting stream $index ($tag deg), phase $phase")
    flush(stdout)
    open(log_path, "w") do log
        run(pipeline(command, stdout=log, stderr=log))
    end
    println("[$(now())] completed stream $index ($tag deg), phase $phase")
    flush(stdout)
end

function main()
    aerosol = parse_indices("GLINT_AEROSOL_STREAM_INDICES", "6,7,8,9")
    no_aerosol = parse_indices("GLINT_NO_AEROSOL_STREAM_INDICES", "")
    wait_for_predecessor()
    for index in aerosol
        run_phase(index, "aerosol")
    end
    for index in no_aerosol
        run_phase(index, "no_aerosol")
    end
    println("[$(now())] native-angle sequence complete")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
