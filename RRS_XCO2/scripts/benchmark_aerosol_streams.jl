#!/usr/bin/env julia

"""Test δBGE angular caps for a representative 16-layer aerosol-on scene."""

include(joinpath(@__DIR__, "generate_truth_map.jl"))
using JLD2

const STREAM_OUT = joinpath(ROOT, "truth_map_aerosols", "stream_convergence.jld2")
const CAPS = (8, 12, 16, 24, 32, 40, 48, 64)
const NSPEC_TEST = parse(Int, get(ENV, "NSPEC_TEST", "33"))

function representative_grid(iband)
    full = output_grids()[iband]
    idx = unique(round.(Int, range(1, length(full); length=NSPEC_TEST)))
    return full[idx]
end

function build_test_model(state, iband, ν, cap)
    p = load_parameters()
    set_common!(p, state)
    p.nstreams = cld(cap + 1, 2)
    p.l_trunc = cap
    p.max_m = cap + 1
    p.truncation = vSmartMOM.Scattering.δBGE{FT}(cap, zero(FT))
    surface = CoreRT.LambertianSurfaceLegendre(FT.(collect(state.coeff[iband])))
    select_band!(p, iband, ν, surface)
    return model_from_parameters(p; external_solar=false)
end

function run_case(state, iband, cap, solar_T)
    ν = representative_grid(iband)
    model = build_test_model(state, iband, ν, cap)
    F0, sources = source_set(ν, false, solar_T)
    elapsed = @elapsed result = CoreRT.rt_run_test(make_nors(F0), model, 1; sources)
    radiance = toa3(result[1])
    return (; ν=Array(ν), radiance, elapsed,
            nstreams=model.quad_points.Nstreams,
            nquad=model.quad_points.Nquad,
            m_max=model.solver.m_max_bands[1],
            l_max=model.solver.l_max[1])
end

function main_streams()
    CUDA.functional() || error("CUDA is not functional")
    CUDA.device!(DEVICE)
    state = read_states()[9] # urban, aerosol on, no SIF, 380 ppm CO2
    solar_T = solar_interpolator()
    results = Dict{String,Any}()
    for iband in 1:3, cap in CAPS
        @info "stream convergence" iband cap
        results["band$(iband)_cap$(cap)"] = run_case(state, iband, cap, solar_T)
        CUDA.synchronize(); GC.gc(); CUDA.reclaim()
    end
    jldsave(STREAM_OUT; results, caps=collect(CAPS), nspec_test=NSPEC_TEST,
            state_index=state.index)
    println(STREAM_OUT)
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main_streams()
