#!/usr/bin/env julia

"""Compare δBGE caps with a full-Mie, untruncated angular reference."""

include(joinpath(@__DIR__, "generate_truth_map.jl"))
using JLD2

const UNTRUNC_OUT = get(ENV, "UNTRUNC_OUT",
    joinpath(ROOT, "truth_map_aerosols", "untruncated_convergence_cpu.jld2"))
const TEST_CAPS = (8, 12, 16, 24, 32)
const NSPEC_UNTRUNC = parse(Int, get(ENV, "NSPEC_UNTRUNC", "1"))
const BANDS_UNTRUNC = Tuple(parse.(Int, split(get(ENV, "BANDS_UNTRUNC", "1"), ',')))
const TEST_VZA = FT.(parse.(Float64,
    split(get(ENV, "VZAS_UNTRUNC", join(0:5:70, ',')), ',')))
const TEST_VAZ = FT[0, 180]

function sparse_grid(iband)
    full = output_grids()[iband]
    # The CPU backend supports a genuine singleton spectral dimension.  Keep
    # only the shortest wavelength (largest wavenumber) for this experiment.
    NSPEC_UNTRUNC == 1 && return FT[maximum(full)]
    idx = unique(round.(Int, range(1, length(full); length=NSPEC_UNTRUNC)))
    return full[idx]
end

function full_mie_requirements(ν)
    λmin = 1e4 / maximum(ν) # micron
    x_max = 2π * 10.0 / λmin
    mie_nmax = vSmartMOM.Scattering.get_n_max(x_max)
    ncoeff = 2 * mie_nmax - 1
    max_order = ncoeff - 1
    return (; λmin, x_max, mie_nmax, ncoeff, max_order,
            nstreams=cld(max_order + 1, 2))
end

function model_for_case(state, iband, ν, cap, vza; external_solar=true)
    p = load_parameters()
    set_common!(p, state)
    p.architecture = CPU()
    # Two paired look geometries share one distinct camera zenith node.
    p.vza = FT[vza, vza]
    p.vaz = copy(TEST_VAZ)
    if cap === nothing
        req = full_mie_requirements(ν)
        p.nstreams = req.nstreams
        p.l_trunc = req.ncoeff
        p.max_m = req.max_order + 1
        p.truncation = vSmartMOM.Scattering.NoTruncation()
    else
        p.nstreams = cld(cap, 2)
        p.l_trunc = cap
        p.max_m = cap
        p.truncation = vSmartMOM.Scattering.δBGE{FT}(cap, zero(FT))
    end
    # The singleton is the maximum-wavenumber (shortest-wavelength) endpoint,
    # hence x=+1 in the band's Legendre coordinate.  P₀(1)=P₁(1)=P₂(1)=1.
    surface = CoreRT.LambertianSurfaceScalar(FT(sum(state.coeff[iband])))
    select_band!(p, iband, ν, surface)
    return model_from_parameters(p; external_solar)
end

function execute(state, iband, ν, cap, solar_T, vza)
    model = model_for_case(state, iband, ν, cap, vza)
    F0, sources = source_set(ν, false, solar_T)
    elapsed = @elapsed result = CoreRT.rt_run_toa(model; i_band=1, sources)
    radiance = FT.(Array(result)[:, 1:3, :])
    return (; radiance, ν=Array(ν), vza=FT[vza, vza], vaz=copy(TEST_VAZ), elapsed,
            nstreams=model.quad_points.Nstreams, nquad=model.quad_points.Nquad,
            m_max=model.solver.m_max_bands[1], l_max=model.solver.l_max[1])
end

function main_untruncated()
    state = read_states()[9]
    solar_T = solar_interpolator()
    if isfile(UNTRUNC_OUT)
        saved = load(UNTRUNC_OUT)
        results = Dict{String,Any}(saved["results"])
        requirements = Dict{String,Any}(saved["requirements"])
        @info "resuming checkpoint" UNTRUNC_OUT completed=length(results)
    else
        results = Dict{String,Any}()
        requirements = Dict{String,Any}()
    end
    for iband in BANDS_UNTRUNC
        ν = sparse_grid(iband)
        requirements["band$(iband)"] = full_mie_requirements(ν)
        for vza in TEST_VZA
            vtag = lpad(string(Int(vza)), 2, '0')
            refkey = "band$(iband)_vza$(vtag)_untruncated"
            if haskey(results, refkey)
                @info "skipping completed VZA" iband vza
                continue
            end
            for cap in TEST_CAPS
                @info "δBGE reference comparison" iband vza cap
                results["band$(iband)_vza$(vtag)_cap$(cap)"] =
                    execute(state, iband, ν, cap, solar_T, vza)
                GC.gc()
            end
            @info "full untruncated reference" iband vza requirements["band$(iband)"]
            results[refkey] =
                execute(state, iband, ν, nothing, solar_T, vza)
            GC.gc()
            # Each VZA includes both requested relative azimuths.  Checkpoint
            # here because a full 111-stream reference is intentionally slow.
            jldsave(UNTRUNC_OUT; results, requirements,
                    caps=collect(TEST_CAPS), nspec=NSPEC_UNTRUNC,
                    state_index=state.index, vza=TEST_VZA, vaz=TEST_VAZ,
                    architecture="CPU")
        end
    end
    println(UNTRUNC_OUT)
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main_untruncated()
