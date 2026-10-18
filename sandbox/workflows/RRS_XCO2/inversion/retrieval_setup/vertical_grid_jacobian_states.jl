#!/usr/bin/env julia
# Round 7_altgrd validation driver (see ../ROUND7_ALTGRD.md sections 4-5).
# Run with --project=<worktree> and the Round 7 environment exported
# (configure_round7_altgrd, or the same variables by hand). Physical GPU 1 only:
#   CUDA_VISIBLE_DEVICES=<GPU1 UUID> CUDA_DEVICE=0
# Evaluate the core OCO forward model + Jacobian at fixed test states and
# serialize them, for bit-identity comparison between checkouts/conventions.
#   OUT=<file> [GRID=pressure|altitude] [ARCH=GPU|CPU] [FTYPE=Float32|Float64]
# Without GRID the evaluator uses its default (also works on 7acab570).
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

cases = Pair{String,Vector{Float64}}[]
for (label, index, dp) in (("s011_truth", 11, 0.0), ("s011_truth_p+3.7", 11, 3.7),
                           ("s036_truth", 36, 0.0), ("s036_truth_p-6.2", 36, -6.2))
    truth = bytruth(index)
    map = R5.resolve_round5_sif_map([truth])
    x = truth_core_state(truth, core.base_parameters, map)
    x[1] += dp
    push!(cases, label => x)
end
let truth = bytruth(1), map = R5.resolve_round5_sif_map([truth])
    prior = load_retrieval_prior(truth.surface; path=ENV["PRIOR_OFF"])
    push!(cases, "s001_prior" => expand_round5_state(map, prior.xa))
end
results = Dict{String,Any}()
for (label, x) in cases
    set_fixed_upper_co2_ppm!(core, 400.0)
    e = core(x)
    results[label] = (; state=x, measurement=e.measurement, jacobian=e.jacobian,
        xco2=column_averaged_co2_ppm(core, x))
    println(label, " xco2=", results[label].xco2, " |K1|=", norm(e.jacobian[:, 1]))
end
# Repeat one state to measure run-to-run determinism of the backend.
e = core(cases[1].second)
results["repeat_s011_truth"] = (; measurement=e.measurement, jacobian=e.jacobian)
println("repeat identical: ", e.measurement == results["s011_truth"].measurement,
        " ", e.jacobian == results["s011_truth"].jacobian)
serialize(ENV["OUT"], results)
