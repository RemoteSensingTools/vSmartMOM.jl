#!/usr/bin/env julia
# Round 7_altgrd validation driver (see ../ROUND7_ALTGRD.md sections 4-5).
# Run with --project=<worktree> and the Round 7 environment exported
# (configure_round7_altgrd, or the same variables by hand). Physical GPU 1 only:
#   CUDA_VISIBLE_DEVICES=<GPU1 UUID> CUDA_DEVICE=0
# Central finite difference of the full instrument-space measurement in psurf
# against the analytic psurf Jacobian column, scenes 011 and 036 truth states.
#   OUT=<file> GRID=pressure|altitude [FTYPE=Float32|Float64] [STEPS=2.0,0.5]
using Serialization, LinearAlgebra
const W = dirname(Base.active_project())
const INV = joinpath(W, "sandbox", "workflows", "RRS_XCO2", "inversion")
ENV["PRIOR_OFF"] = get(ENV, "PRIOR_OFF", replace(ENV["RETRIEVAL_PRIOR_PATH"], "_on_" => "_off_"))
for name in ("OptimalEstimation.jl", "RetrievalCases.jl", "RetrievalState.jl",
             "RetrievalOutput.jl", "Round4SIFTruthConvention.jl",
             "Round5FixedSIF.jl", "Round5RetrievalCampaign.jl")
    include(joinpath(INV, name))
end
using .OptimalEstimation, .RetrievalCases, .RetrievalState, .Round5FixedSIF
import .Round5RetrievalCampaign as R5
include(joinpath(INV, "VSmartMOMForward.jl"))
using .VSmartMOMForward
using vSmartMOM, vSmartMOM.CoreRT

arch = Symbol(get(ENV, "ARCH", "GPU"))
FT = get(ENV, "FTYPE", "Float32") == "Float64" ? Float64 : Float32
grid_kw = haskey(ENV, "GRID") ?
    (; vertical_grid = ENV["GRID"] == "altitude" ? FixedAltitudeGrid() : FixedPressureGrid()) : (;)
core = OCOForwardEvaluator(; architecture=arch, float_type=FT, nstreams=9,
    coefficient_path=ENV["RETRIEVAL_STOKES_COEFFICIENT_PATH"],
    component_path=ENV["RETRIEVAL_SCENE_COMPONENTS_PATH"], grid_kw...)
truths = read_truth_cases(ENV["RETRIEVAL_TRUTH_TABLE"])
bytruth(i) = only(filter(t -> t.state_index == i, truths))

const ROWS = Dict(parse(Int, r["index"]) => r for r in RetrievalCases._truth_rows(ENV["RETRIEVAL_TRUTH_TABLE"]))
function truth_core_state(truth, params, map)
    row = ROWS[truth.state_index]
    f(name) = parse(Float64, row[name])
    x = zeros(Float64, 30)
    x[1] = 1000.0
    x[VSmartMOMForward.CO2_RANGE] .= 400e-6
    x[last(VSmartMOMForward.CO2_RANGE)] = truth.bottom_co2_ppm * 1e-6
    scale = VSmartMOMForward.tau_ref_per_aod760(ENV["RETRIEVAL_SCENE_COMPONENTS_PATH"])
    x[VSmartMOMForward.LOG_AOD_RANGE] .= log.([f("sulfate_aod550"), f("organic_aod550"), f("stratospheric_aod550")] ./ scale)
    x[VSmartMOMForward.LOG_HEIGHT_RANGE] .= [aer.profile.μ for aer in params.scattering_params.rt_aerosols]
    for iband in 1:3
        i0 = first(VSmartMOMForward.SURFACE_RANGE) + 3 * (iband - 1)
        prefix = ("o2a", "weak", "strong")[iband]
        x[i0:i0+2] .= [f("$(prefix)_P$k") for k in 0:2]
    end
    x[29:30] .= expand_round5_state(map, zeros(28))[29:30]
    return x
end


using NCDatasets
const OBS = ENV["ROUND7_OBSERVATION_DIR"]
results = Dict{String,Any}()
for index in (11, 36)
    truth = bytruth(index)
    map = R5.resolve_round5_sif_map([truth])
    x = truth_core_state(truth, core.base_parameters, map)
    noise = NCDataset(joinpath(OBS, "OCO2round7_$(lpad(index, 3, '0')).nc")) do ds
        Float64.(ds["noise_standard_deviation"][:])
    end
    set_fixed_upper_co2_ppm!(core, 400.0)
    e0 = core(x)
    K1 = e0.jacobian[:, 1]
    for h in parse.(Float64, split(get(ENV, "STEPS", "2.0,0.5"), ","))
        xp = copy(x); xp[1] += h
        xm = copy(x); xm[1] -= h
        fp = core(xp).measurement; fm = core(xm).measurement
        fd = (fp .- fm) ./ (2h)
        rel = norm(K1 .- fd) / norm(fd)
        relw = norm((K1 .- fd) ./ noise) / norm(fd ./ noise)
        bands = [norm(K1[r] .- fd[r]) / norm(fd[r]) for r in e0.band_ranges]
        println("state=$index h=$h rel_l2=$rel noise_weighted=$relw bands=$bands")
        results["s$(index)_h$(h)"] = (; K1, fd, rel, relw, bands, noise)
    end
end
serialize(ENV["OUT"], results)
