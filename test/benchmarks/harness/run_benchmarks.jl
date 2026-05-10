# Phase 2 benchmark harness — scenario runner
# =========================================================================
#
# Entry point: loads scenarios.toml, builds each model, invokes the
# designated driver (rt_run / rt_run_ss), collects metrics via metrics.jl,
# and writes per-scenario JSON + output-array JLD2 files under
# <save_dir>/<scenario_name>_metrics.json + <scenario_name>_output.jld2.
#
# Usage:
#   julia --project=test test/benchmarks/harness/run_benchmarks.jl \
#         --save-dir test/benchmarks/baseline_output/sanghavi-unified_<sha>/
#
# Or programmatically:
#   include("test/benchmarks/harness/run_benchmarks.jl")
#   run_benchmarks("test/benchmarks/harness/scenarios.toml",
#                  "test/benchmarks/baseline_output/sanghavi-unified_<sha>/")

using vSmartMOM
using vSmartMOM.CoreRT
using vSmartMOM.InelasticScattering
using TOML
using JSON
using JLD2
using Statistics

include(joinpath(@__DIR__, "metrics.jl"))

"""
    _build_rs_type(rs_type_str, model, iBand)

Build an `AbstractRamanType` for the scenario. For `"noRS"` returns the
float-parametric `noRS{FT}()`; for `"RRS"` constructs a minimal placeholder
RRS from the model's spectral band (matching the pattern in
test_forward_raman_phase1b.jl).
"""
function _build_rs_type(rs_type_str::String, model, iBand::Integer)
    FT = CoreRT.float_type(model)
    if rs_type_str == "noRS"
        return noRS{FT}()
    elseif rs_type_str == "RRS"
        nPol = CoreRT.polarization_type(model).n
        ν    = model.atmosphere.spec_bands[iBand]
        nSpec = length(ν)
        ν̃    = (ν[1] + ν[end]) / 2
        effT = (model.profile.vcd_dry' * model.profile.T) / sum(model.profile.vcd_dry)
        n2, o2 = InelasticScattering.getRamanAtmoConstants(ν̃, effT)

        F₀ = zeros(FT, nPol, nSpec); F₀[1, :] .= 1
        SIF₀ = zeros(FT, nPol, nSpec)

        rs = InelasticScattering.RRS(
            n2 = n2, o2 = o2,
            greek_raman = InelasticScattering.GreekCoefs(
                [FT(1)], [FT(1)], [FT(1)], [FT(1)], [FT(1)], [FT(1)]),
            fscattRayl  = [FT(1)], ϖ_Cabannes = [FT(1)],
            ϖ_λ₁λ₀ = zeros(FT, 1), i_λ₁λ₀ = zeros(Int, 1),
            Z⁻⁺_λ₁λ₀ = zeros(FT, 1, 1), Z⁺⁺_λ₁λ₀ = zeros(FT, 1, 1),
            i_ref = argmin(abs.(ν .- ν̃)),
            n_Raman = 0,
            F₀ = F₀, SIF₀ = SIF₀)
        CoreRT.getRamanSSProp!(rs, 1e7/ν̃, ν)
        return rs
    else
        error("Unknown rs_type: $rs_type_str (expected \"noRS\" or \"RRS\")")
    end
end

function _unpack_driver_output(result, driver::String)
    if driver == "rt_run"
        # rt_run SFI=true return is (R_SFI, T_SFI, ieR_SFI, ieT_SFI, hdr, bhr_uw, bhr_dw)
        # or (R_SFI, T_SFI, ieR_SFI, ieT_SFI) for the 4-tuple path.
        if length(result) >= 4
            return (R_SFI=result[1], T_SFI=result[2],
                    ieR_SFI=result[3], ieT_SFI=result[4])
        end
    elseif driver == "rt_run_ss"
        # rt_run_ss SFI=true return is 6-tuple with hem_R, hem_T tail.
        if length(result) == 6
            return (R_SFI=result[1], T_SFI=result[2],
                    ieR_SFI=result[3], ieT_SFI=result[4],
                    hem_R=result[5], hem_T=result[6])
        elseif length(result) == 4
            return (R_SFI=result[1], T_SFI=result[2],
                    ieR_SFI=result[3], ieT_SFI=result[4])
        end
    end
    error("Unrecognized driver result shape (driver=$driver, len=$(length(result)))")
end

"""
    run_scenario(scen::Dict, save_dir::AbstractString)

Run one scenario and return the metrics NamedTuple.
"""
function run_scenario(scen::Dict, save_dir::AbstractString)
    name    = scen["name"]
    yaml    = scen["yaml"]
    arch    = get(scen, "architecture", "cpu")
    rs_type = scen["rs_type"]
    driver  = scen["driver"]
    n_samp  = get(scen, "n_samples", 5)

    @info "Running scenario" name=name yaml=yaml arch=arch rs_type=rs_type driver=driver

    yaml_full = isabspath(yaml) ? yaml : joinpath(@__DIR__, "..", "..", yaml)
    params = parameters_from_yaml(yaml_full)
    params.architecture = arch == "gpu" ? vSmartMOM.Architectures.GPU() :
                                          vSmartMOM.Architectures.CPU()
    model = model_from_parameters(params)

    iBand = 1

    # Closure that returns the driver output as a NamedTuple. RS_type must be
    # rebuilt per-call because some RS_type fields (ϖ_λ₁λ₀, Z_λ₁λ₀) get
    # mutated inside rt_run.
    run_fn = () -> begin
        rs = _build_rs_type(rs_type, model, iBand)
        raw = driver == "rt_run_ss" ?
                CoreRT.rt_run_test_ss(rs, model, iBand) :
                CoreRT.rt_run_test(rs, model, iBand)
        _unpack_driver_output(raw, driver)
    end

    return capture_run(name, run_fn, save_dir; n=n_samp)
end

function _to_json_serializable(v)
    v === nothing && return nothing
    v isa AbstractString && return v
    v isa Symbol         && return string(v)
    v isa AbstractVector && return [_to_json_serializable(x) for x in v]
    v isa Tuple          && return [_to_json_serializable(x) for x in v]
    v isa NamedTuple     && return Dict(string(k) => _to_json_serializable(getfield(v,k))
                                         for k in keys(v))
    v isa AbstractDict   && return Dict(string(k) => _to_json_serializable(val)
                                         for (k,val) in v)
    return v
end

"""
    run_benchmarks(toml_path, save_dir)

Load all scenarios from `toml_path`, run each, and write both:
  - `<save_dir>/<name>_metrics.json` — per-scenario metrics + hw info
  - `<save_dir>/<name>_output.jld2`  — output arrays from one run
  - `<save_dir>/_summary.json`       — index of all scenarios + hw + git sha
"""
function run_benchmarks(toml_path::AbstractString, save_dir::AbstractString)
    cfg = TOML.parsefile(toml_path)
    scenarios = cfg["scenario"]
    mkpath(save_dir)

    # Record git sha of the current branch for traceability.
    git_sha = try
        chomp(read(`git -C $(dirname(dirname(dirname(@__DIR__)))) rev-parse HEAD`, String))
    catch
        "unknown"
    end

    summary = Dict{String,Any}(
        "git_sha"  => git_sha,
        "toml"     => toml_path,
        "hw"       => _to_json_serializable(hw_identifier()),
        "scenarios"=> String[],
    )

    for scen in scenarios
        metrics = run_scenario(scen, save_dir)
        metrics_path = joinpath(save_dir, scen["name"] * "_metrics.json")
        open(metrics_path, "w") do io
            JSON.print(io, _to_json_serializable(metrics), 2)
        end
        push!(summary["scenarios"], scen["name"])
        @info "Captured" scenario=scen["name"] metrics_path=metrics_path
    end

    summary_path = joinpath(save_dir, "_summary.json")
    open(summary_path, "w") do io
        JSON.print(io, summary, 2)
    end
    @info "Harness complete" summary=summary_path
    return summary
end

# When run as a script via `julia --project=test run_benchmarks.jl`,
# parse CLI args and dispatch.
if abspath(PROGRAM_FILE) == @__FILE__
    save_dir = "test/benchmarks/baseline_output/current"
    toml_path = joinpath(@__DIR__, "scenarios.toml")
    for (i, arg) in enumerate(ARGS)
        if arg == "--save-dir" && i < length(ARGS)
            save_dir = ARGS[i+1]
        elseif arg == "--scenarios" && i < length(ARGS)
            toml_path = ARGS[i+1]
        end
    end
    run_benchmarks(toml_path, save_dir)
end
