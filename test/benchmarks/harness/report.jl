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

# Per-Stokes-component tolerances, Phase 2b (user-specified 2026-04-21).
# I component: relative-only gate (rtol = 1e-4).
# Q, U, V components: absolute-only gate (atol = 1e-8).
# A tolerance of (atol=Inf, rtol=r) means "rtol gate only"; (atol=a, rtol=0)
# means "atol gate only". `within_tol` takes atol-OR-rtol per v2 plan.
const STOKES_TOL_I   = (atol=0.0,  rtol=1e-4)
const STOKES_TOL_QUV = (atol=1e-8, rtol=0.0)

# Mapping from (quantity, stokes_index) => tolerance. Stokes indexing is
# 1=I, 2=Q, 3=U, 4=V as stored along axis 2 of the (Nvza, Nstokes, Nspec)
# SFI arrays.
const DEFAULT_TOLERANCES = Dict{Tuple{String,Int}, NamedTuple}(
    ("R_SFI",   1) => STOKES_TOL_I,
    ("R_SFI",   2) => STOKES_TOL_QUV,
    ("R_SFI",   3) => STOKES_TOL_QUV,
    ("R_SFI",   4) => STOKES_TOL_QUV,
    ("T_SFI",   1) => STOKES_TOL_I,
    ("T_SFI",   2) => STOKES_TOL_QUV,
    ("T_SFI",   3) => STOKES_TOL_QUV,
    ("T_SFI",   4) => STOKES_TOL_QUV,
    ("ieR_SFI", 1) => STOKES_TOL_I,
    ("ieR_SFI", 2) => STOKES_TOL_QUV,
    ("ieR_SFI", 3) => STOKES_TOL_QUV,
    ("ieR_SFI", 4) => STOKES_TOL_QUV,
    ("ieT_SFI", 1) => STOKES_TOL_I,
    ("ieT_SFI", 2) => STOKES_TOL_QUV,
    ("ieT_SFI", 3) => STOKES_TOL_QUV,
    ("ieT_SFI", 4) => STOKES_TOL_QUV,
    # Hemispheric — scalar-like outputs; use I tolerance.
    ("hem_R",   1) => STOKES_TOL_I,
    ("hem_T",   1) => STOKES_TOL_I,
)

const STOKES_LABELS = ("I", "Q", "U", "V")

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

        # Per-Stokes-component comparison. Arrays are (Nvza, Nstokes, Nspec).
        shared_keys = intersect(keys(oa), keys(ob))
        for k in sort(collect(shared_keys))
            a_arr = oa[k]; b_arr = ob[k]
            if !(a_arr isa AbstractArray && b_arr isa AbstractArray)
                continue
            end
            if size(a_arr) != size(b_arr)
                verbose && @printf("  %-8s: SHAPE MISMATCH %s vs %s\n", k,
                                   size(a_arr), size(b_arr))
                all_pass = false
                continue
            end
            nstokes = ndims(a_arr) >= 2 ? size(a_arr, 2) : 1
            for s in 1:nstokes
                slice_a = ndims(a_arr) >= 3 ? selectdim(a_arr, 2, s) :
                          (nstokes == 1 ? a_arr : selectdim(a_arr, 2, s))
                slice_b = ndims(b_arr) >= 3 ? selectdim(b_arr, 2, s) :
                          (nstokes == 1 ? b_arr : selectdim(b_arr, 2, s))
                d = array_diff_summary(slice_b, slice_a)
                results[(scen, k * "_" * STOKES_LABELS[s])] = d
                tol = get(tolerances, (k, s), (atol=Inf, rtol=Inf))
                pass = within_tol(d, tol.atol, tol.rtol)
                all_pass &= pass
                gate = tol.atol == 0    ? @sprintf("rtol≤%.1e",  tol.rtol) :
                       tol.rtol == 0    ? @sprintf("atol≤%.1e",  tol.atol) :
                       tol.atol == Inf  ? "(no gate)" :
                                          @sprintf("atol≤%.1e|rtol≤%.1e", tol.atol, tol.rtol)
                verbose && @printf("  %-10s: max|Δ|=%-10.3g max_rel=%-10.3g  %s  %s\n",
                                   k * "[" * STOKES_LABELS[s] * "]",
                                   d.max_abs_diff, d.max_rel_diff,
                                   gate, pass ? "PASS" : "FAIL")
            end
        end

        # Perf-ceiling flags. Each dimension is "≤ sanghavi" per user policy
        # 2026-04-21 — if the candidate exceeds the target, we FLAG but don't
        # necessarily gate the whole comparison.
        if ma !== nothing && mb !== nothing
            flag_perf!(verbose, ma, mb)
        end
    end

    verbose && println("="^75)
    verbose && println("Overall (physics): $(all_pass ? "PASS" : "FAIL")")
    verbose && println("="^75)

    return (results=results, all_pass=all_pass, shared=shared, A=A, B=B)
end

"""
    flag_perf!(verbose, ma, mb)

Print perf-ceiling flags per user policy 2026-04-21: sanghavi (a) is the
ceiling; any candidate (b) metric that exceeds it is flagged with
`FLAG` (not `FAIL` — perf regressions are tracked, not auto-blocking).
"""
function flag_perf!(verbose::Bool, ma::Dict, mb::Dict)
    function _cmp(label::String, a, b, unit::String="")
        (a === nothing || b === nothing) && return
        # For wall-clock, we use (b > a) as the flag condition.
        # Same for cpu_alloc_bytes and gpu_used_delta.
        a, b = float(a), float(b)
        ratio = b / max(abs(a), eps())
        verdict = b > a ? "FLAG (> sanghavi)" : "ok (≤ sanghavi)"
        verbose && @printf("  perf %-18s: sanghavi=%-12.4g candidate=%-12.4g  ratio=%5.2fx  %s\n",
                           label * unit, a, b, ratio, verdict)
    end
    _cmp("wall_median_s",    get(ma, "wall_median_s",   nothing), get(mb, "wall_median_s",   nothing), " s")
    _cmp("cpu_alloc_bytes",  get(ma, "cpu_alloc_bytes", nothing), get(mb, "cpu_alloc_bytes", nothing), " B")
    _cmp("gpu_used_delta",   get(ma, "gpu_used_delta",  nothing), get(mb, "gpu_used_delta",  nothing), " B")
end
