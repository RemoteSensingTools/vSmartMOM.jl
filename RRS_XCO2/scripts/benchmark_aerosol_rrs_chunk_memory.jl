#!/usr/bin/env julia

"""
Benchmark one Float32 aerosol-on O2 A-band Raman chunk.

The requested `CORE_POINTS` are retained output wavelengths. The actual solve
also contains the standard +/-250 cm^-1 Raman shoulders. Set
`AEROSOL_NSTREAMS` (default 16 here) independently of the production default.
Run each core size in a fresh Julia process so CUDA-pool high-water marks do
not leak between cases.
"""

ENV["RRS_XCO2_FLOAT_TYPE"] = "Float32"
ENV["AEROSOL_NSTREAMS"] = get(ENV, "AEROSOL_NSTREAMS", "16")

include(joinpath(@__DIR__, "generate_truth_map.jl"))

function main()
    CUDA.functional() || error("CUDA is not functional")
    CUDA.device!(DEVICE)
    grids = output_grids()
    ncore = min(parse(Int, get(ENV, "CORE_POINTS", "256")), length(grids[1]))
    first_index = max(1, (length(grids[1]) - ncore) ÷ 2 + 1)
    coreν = grids[1][first_index:first_index+ncore-1]
    solveν, keep = o2_solve_grid(coreν)
    state = first(filter(s -> any(>(0), s.aod550), read_states()))
    solar_T = solar_interpolator()

    CUDA.reclaim()
    free_before = CUDA.available_memory()
    t0 = time_ns()
    result = simulate_o2(state, (coreν, solveν, keep), solar_T)
    CUDA.synchronize()
    elapsed = (time_ns() - t0) / 1e9
    free_after = CUDA.available_memory()

    println("BENCHMARK_RESULT",
            " core_points=", ncore,
            " solve_points=", length(solveν),
            " nstreams=", AEROSOL_NSTREAMS,
            " float_type=", FT,
            " elapsed_s=", elapsed,
            " gpu_used_delta_mib=", (free_before-free_after)/2.0^20,
            " checksum=", sum(abs, result.rayleigh) +
                          sum(abs, result.cabannes) + sum(abs, result.rrs))
    CUDA.pool_status()
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
