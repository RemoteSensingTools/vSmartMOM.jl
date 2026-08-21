#!/usr/bin/env julia

"""Testing-only extraction of m=0 for the existing delta-BGE cap cases."""

include(joinpath(@__DIR__, "benchmark_aerosol_untruncated.jl"))
using JLD2
using DelimitedFiles

const M0_OUT = joinpath(ROOT, "truth_map_aerosols", "truncated_m0_stokes.jld2")

function truncated_m0_model(state, ν, cap)
    p = load_parameters()
    set_common!(p, state)
    p.architecture = CPU()
    p.vza = repeat(TEST_VZA, inner=2)
    p.vaz = repeat(TEST_VAZ, outer=length(TEST_VZA))
    p.nstreams = cld(cap, 2)
    p.l_trunc = cap
    p.max_m = 1                 # testing only: solve exactly m=0
    p.truncation = vSmartMOM.Scattering.δBGE{FT}(cap, zero(FT))
    surface = CoreRT.LambertianSurfaceScalar(FT(sum(state.coeff[1])))
    select_band!(p, 1, ν, surface)
    return model_from_parameters(p; external_solar=false)
end

function main_truncated_m0()
    state = read_states()[9]
    ν = sparse_grid(1)
    _, sources = source_set(ν, false, solar_interpolator())
    m0_values = Dict{Int,Array{FT,3}}()
    elapsed = Dict{Int,Float64}()

    for cap in TEST_CAPS
        @info "truncated m=0 diagnostic" cap
        model = truncated_m0_model(state, ν, cap)
        out = zeros(FT, length(TEST_VZA), 2, 3)
        callback = function (step)
            J = Array(step.composite_layer.J₀⁻)
            for slot in values(step.composite_layer.J₀_by_src)
                J .+= Array(slot.J₀⁻)
            end
            for iv in eachindex(TEST_VZA)
                iμ = CoreRT.nearest_point(step.qp_μ, cosd(TEST_VZA[iv]))
                _, first_stokes, last_stokes = CoreRT.get_indices(iμ, step.pol_type)
                raw = @view J[first_stokes:last_stokes, 1, 1]
                for iaz in eachindex(TEST_VAZ)
                    out[iv, iaz, 1] = step.weight * raw[1]
                    out[iv, iaz, 2] = step.weight * raw[2]
                    out[iv, iaz, 3] = zero(FT)
                end
            end
            nothing
        end
        FTtype = CoreRT.float_type(model)
        elapsed[cap] = @elapsed CoreRT._rt_run_column(
            vSmartMOM.InelasticScattering.noRS{FTtype}(), model, [1];
            sources, streams_callback=callback)
        m0_values[cap] = out
        jldsave(M0_OUT; caps=collect(TEST_CAPS), vza=TEST_VZA, vaz=TEST_VAZ,
                m0_values, elapsed, wavelength_nm=1e7 / ν[1])
    end
    println(M0_OUT)
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main_truncated_m0()
