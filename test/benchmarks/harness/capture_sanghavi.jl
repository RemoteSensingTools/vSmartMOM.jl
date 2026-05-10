# Phase 2b — run the harness against the sanghavi worktree.
# =========================================================================
#
# Must be launched from the sanghavi worktree so that `using vSmartMOM`
# resolves to the sanghavi source tree:
#
#   cd /home/sanghavi/code/github/vSmartMOM.jl
#   julia --project=. /home/sanghavi/code/github/uni_vSmartMOM/test/benchmarks/harness/capture_sanghavi.jl
#
# Writes metrics + output JLD2 into uni_vSmartMOM/test/benchmarks/baseline_output/
# sanghavi_9ee9a75/, overwriting the Phase 2a output-only version.
#
# Sanghavi-side patches (no edits to sanghavi source, per read-only policy):
# 1. Monkey-patch the 21-arg `postprocessing_vza!(::RRS, ...)` call that
#    sanghavi rt_run invokes but whose 18-arg RRS implementation hasn't
#    been updated for the deferred hem_R/hem_T API. Forwards 21 → 18 args.
# 2. The scenario runner here uses `model.params.polarization_type` and
#    `model.params.spec_bands[iBand]` instead of unified's `RTModel`
#    sub-struct accessors.

const UNI_ROOT = "/home/sanghavi/code/github/uni_vSmartMOM"
const SAVE_DIR = joinpath(UNI_ROOT, "test/benchmarks/baseline_output/sanghavi_9ee9a75")

using vSmartMOM
using vSmartMOM.CoreRT
using vSmartMOM.InelasticScattering
using TOML
using JSON
using JLD2
using Statistics

include(joinpath(UNI_ROOT, "test/benchmarks/harness/metrics.jl"))

# --- Patch 1: 21-arg postprocessing_vza! forwarder for sanghavi RRS path ---
@eval vSmartMOM.CoreRT begin
    function postprocessing_vza!(
        RS_type::Union{InelasticScattering.RRS, InelasticScattering.VS_0to1_plus, InelasticScattering.VS_1to0_plus},
        iμ₀, pol_type, composite_layer, vza, qp_μ, m, vaz, μ₀, weight,
        nSpec, SFI, R, R_SFI, T, T_SFI, ieR_SFI, ieT_SFI,
        hem_R, hem_T, wt_μ)
        # Forward to the 18-arg RRS implementation; hem_R/hem_T/wt_μ dropped
        # per amendments §2.3 deferral on sanghavi.
        postprocessing_vza!(RS_type, iμ₀, pol_type, composite_layer,
            vza, qp_μ, m, vaz, μ₀, weight,
            nSpec, SFI, R, R_SFI, T, T_SFI, ieR_SFI, ieT_SFI)
    end
end

# --- RS_type builder, sanghavi flavor ---
function _build_rs_type_sanghavi(rs_type_str::String, model, iBand::Integer)
    FT = model.params.float_type
    pol  = model.params.polarization_type
    nPol = pol.n
    ν    = model.params.spec_bands[iBand]
    nSpec = length(ν)
    F₀_init = zeros(FT, nPol, nSpec); F₀_init[1, :] .= 1
    SIF₀_init = zeros(FT, nPol, nSpec)
    if rs_type_str == "noRS"
        return noRS{FT}(F₀ = F₀_init, SIF₀ = SIF₀_init)
    elseif rs_type_str == "RRS"
        ν̃    = (ν[1] + ν[end]) / 2
        effT = (model.profile.vcd_dry' * model.profile.T) / sum(model.profile.vcd_dry)
        n2, o2 = InelasticScattering.getRamanAtmoConstants(ν̃, effT)

        rs = InelasticScattering.RRS(
            n2 = n2, o2 = o2,
            greek_raman = InelasticScattering.GreekCoefs(
                [FT(1)], [FT(1)], [FT(1)], [FT(1)], [FT(1)], [FT(1)]),
            fscattRayl  = [FT(1)], ϖ_Cabannes = [FT(1)],
            ϖ_λ₁λ₀ = zeros(FT, 1), i_λ₁λ₀ = zeros(Int, 1),
            Z⁻⁺_λ₁λ₀ = zeros(FT, 1, 1), Z⁺⁺_λ₁λ₀ = zeros(FT, 1, 1),
            i_ref = argmin(abs.(ν .- ν̃)),
            n_Raman = 0,
            F₀ = F₀_init, SIF₀ = SIF₀_init)
        CoreRT.getRamanSSProp!(rs, 1e7/ν̃, ν)
        return rs
    else
        error("Unknown rs_type: $rs_type_str")
    end
end

# --- Unpack driver output into a NamedTuple ---
function _unpack(result, driver::String)
    if driver == "rt_run"
        # sanghavi rt_run with SFI returns 6-tuple (R_SFI, T_SFI, ieR_SFI, ieT_SFI, hem_R, hem_T)
        # for RRS, or (R_SFI, T_SFI, ieR_SFI, ieT_SFI) for noRS.
        n = length(result)
        if n == 6
            return (R_SFI=result[1], T_SFI=result[2],
                    ieR_SFI=result[3], ieT_SFI=result[4],
                    hem_R=result[5], hem_T=result[6])
        elseif n == 4
            return (R_SFI=result[1], T_SFI=result[2],
                    ieR_SFI=result[3], ieT_SFI=result[4])
        end
    elseif driver == "rt_run_ss"
        n = length(result)
        if n == 6
            return (R_SFI=result[1], T_SFI=result[2],
                    ieR_SFI=result[3], ieT_SFI=result[4],
                    hem_R=result[5], hem_T=result[6])
        elseif n == 4
            return (R_SFI=result[1], T_SFI=result[2],
                    ieR_SFI=result[3], ieT_SFI=result[4])
        end
    end
    error("Unrecognized driver result shape (driver=$driver, len=$(length(result)))")
end

function run_scenario_sanghavi(scen::Dict, save_dir::AbstractString)
    name    = scen["name"]
    yaml    = scen["yaml"]
    arch    = get(scen, "architecture", "cpu")
    rs_type = scen["rs_type"]
    driver  = scen["driver"]
    n_samp  = get(scen, "n_samples", 5)

    @info "Running (sanghavi)" name=name yaml=yaml rs_type=rs_type driver=driver

    yaml_full = isabspath(yaml) ? yaml : joinpath(UNI_ROOT, "test", yaml)
    params = parameters_from_yaml(yaml_full)
    params.architecture = arch == "gpu" ? vSmartMOM.Architectures.GPU() :
                                          vSmartMOM.Architectures.CPU()
    model = model_from_parameters(params)

    iBand = 1

    run_fn = () -> begin
        rs = _build_rs_type_sanghavi(rs_type, model, iBand)
        raw = driver == "rt_run_ss" ?
                CoreRT.rt_run_ss(rs, model, iBand) :
                CoreRT.rt_run(rs, model, iBand)
        _unpack(raw, driver)
    end

    return capture_run(name, run_fn, save_dir; n=n_samp)
end

function _to_json(v)
    v === nothing && return nothing
    v isa AbstractString && return v
    v isa Symbol         && return string(v)
    v isa AbstractVector && return [_to_json(x) for x in v]
    v isa Tuple          && return [_to_json(x) for x in v]
    v isa NamedTuple     && return Dict(string(k) => _to_json(getfield(v,k))
                                         for k in keys(v))
    v isa AbstractDict   && return Dict(string(k) => _to_json(val)
                                         for (k,val) in v)
    return v
end

function main()
    toml_path = joinpath(UNI_ROOT, "test/benchmarks/harness/scenarios.toml")
    cfg = TOML.parsefile(toml_path)
    scenarios = cfg["scenario"]
    mkpath(SAVE_DIR)

    git_sha = try
        chomp(read(`git -C /home/sanghavi/code/github/vSmartMOM.jl rev-parse HEAD`, String))
    catch
        "9ee9a75"
    end

    summary = Dict{String,Any}(
        "git_sha"  => git_sha * " (sanghavi tip at Phase 2b freeze)",
        "toml"     => toml_path,
        "hw"       => _to_json(hw_identifier()),
        "scenarios"=> String[],
        "note"     => "Phase 2b sanghavi-side baseline captured via capture_sanghavi.jl with 21-arg postprocessing_vza! monkey-patch."
    )

    for scen in scenarios
        m = run_scenario_sanghavi(scen, SAVE_DIR)
        path = joinpath(SAVE_DIR, scen["name"] * "_metrics.json")
        open(path, "w") do io
            JSON.print(io, _to_json(m), 2)
        end
        push!(summary["scenarios"], scen["name"])
        @info "Captured" scenario=scen["name"] metrics=path
    end

    open(joinpath(SAVE_DIR, "_summary.json"), "w") do io
        JSON.print(io, summary, 2)
    end
    @info "Sanghavi baseline capture complete" save_dir=SAVE_DIR
end

main()
