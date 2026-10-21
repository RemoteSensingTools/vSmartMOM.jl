#!/usr/bin/env julia

"""Diagnose full-Mie VZA smoothness with paired views in one operator."""

include(joinpath(@__DIR__, "benchmark_aerosol_untruncated.jl"))
using JLD2

const DIAG_VZAS = FT.(parse.(Float64, split(get(ENV, "VZAS_DIAG", "20,25,30"), ',')))
const DIAG_FULL = get(ENV, "FULL_DIAG", "0") == "1"
const DIAG_SS = get(ENV, "SS_DIAG", "0") == "1"
const DIAG_OUT = get(ENV, "DIAG_OUT", joinpath(
    ROOT, "truth_map_aerosols", "untruncated_vza_diagnostic.jld2"))

function diagnostic_model(state, ν)
    p = load_parameters()
    set_common!(p, state)
    p.architecture = CPU()
    p.vza = repeat(DIAG_VZAS, inner=2)
    p.vaz = repeat(TEST_VAZ, outer=length(DIAG_VZAS))
    req = full_mie_requirements(ν)
    p.nstreams = req.nstreams
    p.l_trunc = req.ncoeff
    p.max_m = req.max_order + 1
    p.truncation = vSmartMOM.Scattering.NoTruncation()
    surface = CoreRT.LambertianSurfaceScalar(FT(sum(state.coeff[1])))
    select_band!(p, 1, ν, surface)
    return model_from_parameters(p; external_solar=false), req
end

function main_diag()
    state = read_states()[9]
    ν = sparse_grid(1)
    model, requirements = diagnostic_model(state, ν)
    F0, sources = source_set(ν, false, solar_interpolator())

    ss_elapsed = NaN
    ss_radiance = nothing
    if DIAG_SS
        ss_elapsed = @elapsed ss = CoreRT.rt_run_ss(model; i_band=1, sources)
        ss_radiance = Array(ss[1])[:, 1:3, :]
    end

    full_elapsed = NaN
    full_radiance = nothing
    if DIAG_FULL
        full_elapsed = @elapsed full = CoreRT.rt_run(model; i_band=1, sources)
        full_radiance = Array(full.toa)[:, 1:3, :]
    end
    jldsave(DIAG_OUT; vza=Array(model.obs_geom.vza), vaz=Array(model.obs_geom.vaz),
            ν=Array(ν), requirements, nquad=model.quad_points.Nquad,
            nstreams=model.quad_points.Nstreams, ss_elapsed, ss_radiance,
            full_elapsed, full_radiance)
    println(DIAG_OUT)
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main_diag()
