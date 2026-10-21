#!/usr/bin/env julia

"""
Recompute one aerosol O₂ truth-map chunk for a controlled Raman-shoulder test.

Environment variables:
- `CHUNK_INDEX` (default `43`) selects the 64-point production chunk.
- `RAMAN_SHOULDER_CM` is consumed by `generate_truth_map.jl`.
- `DIAGNOSTIC_OUT` selects the JLD2 result path.
- `STATE_INDEX` (default `57`) selects a representative aerosol scene.

This diagnostic deliberately uses the same production construction so that
two runs differing only in `RAMAN_SHOULDER_CM` isolate shoulder sensitivity.
It does not modify truth-map NetCDF files or checkpoints.
"""

ENV["RRS_XCO2_FLOAT_TYPE"] = "Float32"
ENV["AEROSOL_NSTREAMS"] = get(ENV, "AEROSOL_NSTREAMS", "9")

include(joinpath(@__DIR__, "generate_truth_map.jl"))
using JLD2

chunk_ranges(n::Int, width::Int) =
    [i:min(i + width - 1, n) for i in 1:width:n]

function main()
    CUDA.functional() || error("CUDA is not functional")
    CUDA.device!(DEVICE)

    width = 64
    ichunk = parse(Int, get(ENV, "CHUNK_INDEX", "43"))
    state_index = parse(Int, get(ENV, "STATE_INDEX", "57"))
    output = get(ENV, "DIAGNOSTIC_OUT",
        joinpath(TRUTH_ROOT, "aerosol_chunked",
                 "shoulder_diagnostic_chunk$(ichunk).jld2"))

    full_ν = output_grids()[1]
    ranges = chunk_ranges(length(full_ν), width)
    1 <= ichunk <= length(ranges) || error("invalid CHUNK_INDEX=$ichunk")
    irange = ranges[ichunk]
    core_ν = full_ν[irange]
    solve_ν, keep = o2_solve_grid(core_ν)
    state = only(filter(s -> s.index == state_index, read_states()))
    solar_T = solar_interpolator()

    t0 = time_ns()
    result = simulate_o2(state, (full_ν, solve_ν, keep), solar_T)
    CUDA.synchronize()
    elapsed_s = (time_ns() - t0) / 1e9

    jldsave(output;
        state_index, ichunk, irange=collect(irange),
        shoulder_cm=SHOULDER_CM, core_ν, solve_ν,
        rayleigh=Array(result.rayleigh),
        cabannes=Array(result.cabannes), rrs=Array(result.rrs), elapsed_s)
    println("SHOULDER_DIAGNOSTIC output=$output shoulder_cm=$SHOULDER_CM " *
            "chunk=$ichunk ncore=$(length(core_ν)) nsolve=$(length(solve_ν)) " *
            "elapsed_s=$elapsed_s")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
