# Phase 2 benchmark harness — side-by-side comparison + tolerance gate
# =========================================================================
#
# Loads two baseline_output/<branch>_<sha>/ directories (e.g., sanghavi_9ee9a75/
# and sanghavi-unified_b4a5ab1/) and produces:
#   * Per-scenario per-quantity relative-difference summary (max |u-s|, max
#     relative error, median relative error, ratio-stats).
#   * Pass/fail gate against a tolerance table (per-quantity atol/rtol).
#   * Wall-clock and GPU-memory comparison (primary gate per v2 plan Phase 2).
#
# Tolerance table is user-specified before Phase 2b baseline freeze. Defaults
# are permissive so the first diagnostic run (Phase 2a) surfaces deltas
# without blocking.
#
# Usage (programmatic):
#   include("test/benchmarks/harness/report.jl")
#   r = compare_baselines(
#         "test/benchmarks/baseline_output/sanghavi_9ee9a75/",
#         "test/benchmarks/baseline_output/sanghavi-unified_b4a5ab1/")

using JSON
using JLD2
using Statistics
using Printf

# Default tolerance table (LOOSE — used only when `compare_baselines` is
# called without explicit tolerances). Replace with user-provided numbers
# before Phase 2b freeze.
const DEFAULT_TOLERANCES = Dict{String, NamedTuple}(
    # quantity => (atol, rtol)
    "R_SFI"   => (atol=1e-6, rtol=0.05),
    "T_SFI"   => (atol=1e-6, rtol=0.05),
    "ieR_SFI" => (atol=1e-6, rtol=0.05),
    "ieT_SFI" => (atol=1e-6, rtol=0.05),
    "hem_R"   => (atol=1e-6, rtol=0.05),
    "hem_T"   => (atol=1e-6, rtol=0.05),
)

"""
    load_baseline(dir)

Scan a baseline directory for `_summary.json`, then load each
`<name>_metrics.json` + `<name>_output.jld2` referenced inside.
Returns a Dict `scenario_name => (metrics, outputs)` plus the summary.
"""
function load_baseline(dir::AbstractString)
    summary_path = joinpath(dir, "_summary.json")
    isfile(summary_path) || error("No _summary.json under $dir")
    summary = JSON.parsefile(summary_path)

    scenarios = Dict{String,NamedTuple}()
    for name in summary["scenarios"]
        metrics_path = joinpath(dir, name * "_metrics.json")
        output_path  = joinpath(dir, name * "_output.jld2")
        metrics = isfile(metrics_path) ? JSON.parsefile(metrics_path) : nothing
        outputs = isfile(output_path)  ? load(output_path)             : Dict{String,Any}()
        scenarios[name] = (metrics=metrics, outputs=outputs)
    end
    return (summary=summary, scenarios=scenarios)
end

"""
    array_diff_summary(a, b)

Elementwise delta summary between two arrays. Returns
NamedTuple of (max_abs_diff, max_rel_diff, median_rel_diff, ratio_min,
ratio_max, ratio_median, shape_a, shape_b).
"""
function array_diff_summary(a::AbstractArray, b::AbstractArray)
    if size(a) != size(b)
        return (max_abs_diff=NaN, max_rel_diff=NaN, median_rel_diff=NaN,
                ratio_min=NaN, ratio_max=NaN, ratio_median=NaN,
                shape_a=size(a), shape_b=size(b),
                note="shape mismatch")
    end
    diff = a .- b
    ref  = max.(abs.(b), eps(eltype(b)))
    rel  = abs.(diff) ./ ref
    # Ratio stats — uninformative when b has zeros, so mask.
    nonzero = b .!= 0
    ratios = any(nonzero) ? (a[nonzero] ./ b[nonzero]) : [NaN]
    return (max_abs_diff=maximum(abs.(diff)),
            max_rel_diff=maximum(rel),
            median_rel_diff=median(rel),
            ratio_min=minimum(ratios), ratio_max=maximum(ratios),
            ratio_median=median(ratios),
            shape_a=size(a), shape_b=size(b))
end

"""
    within_tol(diff_summary, atol, rtol)

True iff max_abs_diff ≤ atol or max_rel_diff ≤ rtol.
"""
function within_tol(d::NamedTuple, atol::Real, rtol::Real)
    isnan(d.max_abs_diff) && return false
    return d.max_abs_diff <= atol || d.max_rel_diff <= rtol
end

"""
    compare_baselines(a_dir, b_dir; tolerances=DEFAULT_TOLERANCES, verbose=true)

Load two baselines, iterate shared scenarios + shared output keys, print a
formatted report to stdout, and return a Dict of (scenario, quantity) =>
per-quantity diff summary. `a_dir` is the "target" (e.g., sanghavi —
physics truth), `b_dir` is the "candidate" (e.g., sanghavi-unified).
"""
function compare_baselines(a_dir::AbstractString, b_dir::AbstractString;
                           tolerances=DEFAULT_TOLERANCES, verbose=true)
    A = load_baseline(a_dir)
    B = load_baseline(b_dir)

    shared = intersect(keys(A.scenarios), keys(B.scenarios))
    verbose && println("="^75)
    verbose && println("Baseline comparison")
    verbose && println("  A (target):    $a_dir  [git=$(get(A.summary, "git_sha", "?"))]")
    verbose && println("  B (candidate): $b_dir  [git=$(get(B.summary, "git_sha", "?"))]")
    verbose && println("  Shared scenarios: $(length(shared))")
    verbose && println("="^75)

    results = Dict{Tuple{String,String},Any}()
    all_pass = true

    for scen in sort(collect(shared))
        ma = A.scenarios[scen].metrics
        mb = B.scenarios[scen].metrics
        oa = A.scenarios[scen].outputs
        ob = B.scenarios[scen].outputs

        verbose && println()
        verbose && println("--- $scen ---")
        if ma !== nothing && mb !== nothing
            wa = get(ma, "wall_median_s", nothing); wb = get(mb, "wall_median_s", nothing)
            aa = get(ma, "cpu_alloc_bytes", nothing); ab = get(mb, "cpu_alloc_bytes", nothing)
            ga = get(ma, "gpu_used_delta", nothing); gb = get(mb, "gpu_used_delta", nothing)
            if !isnothing(wa) && !isnothing(wb)
                verbose && @printf "  wall (s)  : %8.3f -> %8.3f   (ratio %.2fx)\n" wa wb wb/wa
            elseif !isnothing(wb)
                verbose && @printf "  wall (s)  :   (target unmeasured) -> %8.3f\n" wb
            end
            if !isnothing(aa) && !isnothing(ab)
                verbose && @printf "  cpu allocs: %10d -> %10d   (ratio %.2fx)\n" aa ab ab/max(aa,1)
            elseif !isnothing(ab)
                verbose && @printf "  cpu allocs:   (target unmeasured) -> %10d\n" ab
            end
            if !isnothing(ga) && !isnothing(gb)
                verbose && @printf "  gpu Δmem  : %10d -> %10d   bytes\n" ga gb
            end
        end

        shared_keys = intersect(keys(oa), keys(ob))
        for k in sort(collect(shared_keys))
            a_arr = oa[k]; b_arr = ob[k]
            if a_arr isa AbstractArray && b_arr isa AbstractArray
                d = array_diff_summary(b_arr, a_arr)  # note: b is candidate, a is target
                results[(scen, k)] = d
                tol = get(tolerances, k, (atol=Inf, rtol=Inf))
                pass = within_tol(d, tol.atol, tol.rtol)
                all_pass &= pass
                verbose && @printf("  %-8s: max|Δ|=%-10.3g max_rel=%-10.3g median_rel=%-10.3g  ratio[%.4f..%.4f]  %s\n",
                                   k, d.max_abs_diff, d.max_rel_diff, d.median_rel_diff,
                                   d.ratio_min, d.ratio_max,
                                   pass ? "PASS" : "FAIL")
            end
        end
    end

    verbose && println("="^75)
    verbose && println("Overall: $(all_pass ? "PASS" : "FAIL")")
    verbose && println("="^75)

    return (results=results, all_pass=all_pass, shared=shared, A=A, B=B)
end
