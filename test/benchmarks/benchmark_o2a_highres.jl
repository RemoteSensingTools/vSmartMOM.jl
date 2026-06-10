# =============================================================================
# benchmark_o2a_highres.jl
#
# Speed / memory benchmark for vSmartMOM.jl on a high-resolution O2 A-band
# forward run *with an aerosol layer*. Designed to be run unchanged on two
# revisions and have the outputs diffed:
#
#   * origin/main            — legacy in-tree HITRAN LBL engine, Z matrices
#                              repeat-expanded to (NquadN, NquadN, nSpec).
#   * externalAbsorption     — AtmosphericAbsorption.jl LBL engine, GPU Mie,
#                              flat (singleton) Z fast path.
#
# The aerosol exercises the Mie setup, δBGE truncation, aerosol layers in the
# RT, and the GPU Mie path on the branch. The case is built so that the *only*
# thing that varies between revisions is the source code under test, never the
# config: we start from test/test_parameters/O2Parameters.yaml (byte-identical
# on both revisions) and override only the spectral grid (N_SPEC points) in
# memory. WITH_AEROSOL can strip the aerosol for an A/B comparison; the default
# is true (the user's hard requirement).
#
# Measures, per (revision x CPU/GPU):
#   * model_from_parameters wall time (LBL absorption + Mie live here)
#   * rt_run wall time
#   * Julia allocation totals + GC time fraction (via @timed)
#   * peak process RSS (Sys.maxrss())
#   * GPU: CUDA pool / free-memory before & after, pool peak
#   * max |R| (radiance) for a downstream branch-vs-main radiance diff
#
# USAGE
#   cd test
#   # GPU:
#   JULIA_NUM_THREADS=16 julia --project=. benchmarks/benchmark_o2a_highres.jl
#   # override via env vars (all optional):
#   BENCH_N_SPEC=40000 BENCH_ARCH=GPU BENCH_N_REPS=3 BENCH_WITH_AEROSOL=1 \
#       BENCH_BLAS_THREADS=4 BENCH_TAG=branch \
#       julia --project=. benchmarks/benchmark_o2a_highres.jl
#
# The script writes a one-line-per-field summary to stdout and, if BENCH_OUT is
# set, appends a machine-readable key=value block to that file so a wrapper can
# collect both revisions' numbers and diff the saved radiance vector.
# =============================================================================

# ----------------------------- parameters (top) ------------------------------
# Edit here, or override with the BENCH_* environment variables shown above.
const N_SPEC        = parse(Int,  get(ENV, "BENCH_N_SPEC",     "40000"))  # spectral points
const ARCH          =            get(ENV, "BENCH_ARCH",        "GPU")     # "GPU" or "CPU"
const N_REPS        = parse(Int,  get(ENV, "BENCH_N_REPS",     ARCH == "GPU" ? "3" : "1"))
const WITH_AEROSOL  =            get(ENV, "BENCH_WITH_AEROSOL", "1") in ("1", "true", "yes")
const BLAS_THREADS  = parse(Int,  get(ENV, "BENCH_BLAS_THREADS", "4"))
const TAG           =            get(ENV, "BENCH_TAG",         ARCH)      # label for output
const OUT_FILE      =            get(ENV, "BENCH_OUT",         "")        # optional kv dump
const YAML_PATH     =            get(ENV, "BENCH_YAML",
    joinpath(@__DIR__, "..", "test_parameters", "O2Parameters.yaml"))
# Save the band-1 radiance vector here for a branch-vs-main radiance delta.
const RAD_FILE      =            get(ENV, "BENCH_RAD_OUT",     "")        # optional .txt
# -----------------------------------------------------------------------------

using vSmartMOM
using vSmartMOM.CoreRT
using LinearAlgebra
using Printf
using DelimitedFiles
using TimerOutputs   # both revisions' model_from_parameters use @timeit (default timer)

# CUDA is a weak dependency; load it only if we need the GPU and it is present.
const _USE_GPU = (ARCH == "GPU")
const _HAS_CUDA = if _USE_GPU
    try
        @eval using CUDA
        CUDA.functional()
    catch err
        @warn "CUDA requested but not loadable" err
        false
    end
else
    false
end

if _USE_GPU && !_HAS_CUDA
    error("BENCH_ARCH=GPU but CUDA is not functional in this environment.")
end

# --------------------------- helpers -----------------------------------------

gib(bytes) = bytes / 2^30

"Build the params object: load YAML, override the spectral grid to N_SPEC pts, " *
"force the architecture, and optionally drop the aerosol. Identical inputs on " *
"both revisions because the YAML is byte-identical and the overrides are explicit."
function build_params()
    p = parameters_from_yaml(YAML_PATH)

    # --- spectral grid override: keep [ν_start, ν_end] from the YAML, resample
    #     to exactly N_SPEC points. This is the only knob that changes the
    #     workload size; everything else (geometry, profile, surface, aerosol,
    #     quadrature, polarization) comes straight from the shared YAML.
    FT = p.float_type
    νref = p.spec_bands[1]
    ν0, ν1 = FT(first(νref)), FT(last(νref))
    p.spec_bands[1] = collect(range(ν0, ν1, length = N_SPEC))

    # --- architecture override (so one YAML drives both CPU and GPU runs).
    #     `GPU`/`CPU` are exported by vSmartMOM (from the Architectures submodule).
    p.architecture = (ARCH == "GPU") ? GPU() : CPU()

    # --- BLAS threads kept identical on both sides (recorded below)
    if p.numerics !== nothing
        try
            p.numerics.blas_threads = BLAS_THREADS
        catch
            # older main may not have numerics.blas_threads; we also set BLAS
            # globally below, so this is best-effort.
        end
    end

    # --- aerosol toggle
    if !WITH_AEROSOL
        empty!(p.scattering_params.rt_aerosols)
    end

    return p
end

"@timed wrapper that also captures peak RSS delta and (optionally) GPU pool stats."
function timed_stage(f; gpu::Bool)
    GC.gc(true)
    rss0 = Sys.maxrss()
    gpu_free0 = gpu_used0 = -1
    if gpu
        CUDA.reclaim()
        gpu_free0, gpu_tot = CUDA.MemoryInfo().free_bytes, CUDA.MemoryInfo().total_bytes
        gpu_used0 = gpu_tot - gpu_free0
    end
    # Synchronize INSIDE the timed block: CUDA kernels are launched
    # asynchronously, so timing must include queued device work, not just
    # the host-side launch latency (codex review finding).
    stat = if gpu
        @timed begin
            r = f()
            CUDA.synchronize()
            r
        end
    else
        @timed f()
    end
    rss1 = Sys.maxrss()
    gpu_free1 = gpu_used1 = gpu_pool = -1
    if gpu
        gpu_free1 = CUDA.MemoryInfo().free_bytes
        gpu_used1 = CUDA.MemoryInfo().total_bytes - gpu_free1
        gpu_pool = try CUDA.pool_status_bytes() catch; -1 end
    end
    return (value = stat.value, time = stat.time, bytes = stat.bytes,
            gctime = stat.gctime, rss0 = rss0, rss1 = rss1,
            gpu_free0 = gpu_free0, gpu_used0 = gpu_used0,
            gpu_free1 = gpu_free1, gpu_used1 = gpu_used1, gpu_pool = gpu_pool)
end

# CUDA.pool_status_bytes may not exist across CUDA.jl versions — provide a probe.
if _HAS_CUDA
    if !isdefined(CUDA, :pool_status_bytes)
        @eval CUDA pool_status_bytes() = -1
    end
end

# --------------------------- run ---------------------------------------------

# Keep BLAS thread count identical on both revisions.
BLAS.set_num_threads(BLAS_THREADS)

println("="^78)
println("benchmark_o2a_highres.jl  —  tag=$(TAG)")
println("  N_SPEC=$(N_SPEC)  ARCH=$(ARCH)  N_REPS=$(N_REPS)  WITH_AEROSOL=$(WITH_AEROSOL)")
println("  JULIA_NUM_THREADS=$(Threads.nthreads())  BLAS_THREADS=$(BLAS.get_num_threads())")
println("  YAML=$(YAML_PATH)")
if _HAS_CUDA
    dev = CUDA.device()
    mi = CUDA.MemoryInfo()
    println("  GPU=$(CUDA.name(dev))  total=$(round(gib(mi.total_bytes); digits=2)) GiB  " *
            "free@start=$(round(gib(mi.free_bytes); digits=2)) GiB")
end
println("="^78)

p = build_params()
println("  configured aerosols = $(length(p.scattering_params.rt_aerosols))")
println("  spectral points     = $(length(p.spec_bands[1]))  " *
        "[$(round(first(p.spec_bands[1]); digits=2)) … $(round(last(p.spec_bands[1]); digits=2)) cm⁻¹]")

# ---- model_from_parameters (LBL absorption + Mie live here) ----
# 1 warmup (compile + first-touch) then 1 timed build. model build is cheap
# relative to rt_run and not the headline, so a single timed build suffices.
print("\n[model_from_parameters] warmup … "); flush(stdout)
# Cold (first-call) time = Julia/JIT compile + first-touch; reported separately
# below as cold_model_build_s — never folded into the warm medians.
cold_model_build_s = let s = timed_stage(() -> model_from_parameters(p); gpu = _HAS_CUDA)
    @printf("done (%.2f s, warmup)\n", s.time)
    s.time
end

print("[model_from_parameters] timed  … "); flush(stdout)
reset_timer!()   # isolate the timed build in the TimerOutputs default timer
model_stat = timed_stage(() -> model_from_parameters(p); gpu = _HAS_CUDA)
model = model_stat.value
@printf("%.3f s | alloc %.3f GiB | gc %.1f%% | ΔRSS %.3f GiB\n",
        model_stat.time, gib(model_stat.bytes),
        100 * model_stat.gctime / max(model_stat.time, eps()),
        gib(model_stat.rss1 - model_stat.rss0))

# LBL-vs-Mie split: both revisions wrap "Read HITRAN"/"Absorption Coeff"/"Mie calc"
# in @timeit against the TimerOutputs default timer.
println("  [model_from_parameters breakdown — TimerOutputs default timer]")
try
    show(stdout, TimerOutputs.get_defaulttimer(); allocations = true, compact = true)
    println()
catch err
    @warn "could not print TimerOutputs breakdown" err
end

ns = CoreRT.polarization_type(model).n
NquadN = model.quad_points.Nquad * ns
@printf("  NquadN = %d (Nquad=%d × n_stokes=%d) ; n_aerosols(model)=%d\n",
        NquadN, model.quad_points.Nquad, ns, CoreRT.n_aerosols(model))
# main repeat-expands Z to (NquadN, NquadN, nSpec); estimate one expanded Z array.
zbytes = NquadN^2 * N_SPEC * sizeof(p.float_type)
@printf("  one repeat-expanded Z array @ nSpec=%d ≈ %.3f GiB (main); flat on branch\n",
        N_SPEC, gib(zbytes))

# ---- rt_run (the RT solver) ----
print("\n[rt_run] warmup … "); flush(stdout)
# Cold rt_run is dominated by kernel compilation (KA/CUBLAS on GPU) — at 40k
# points typically ~10x the warm time. Reported as cold_rt_run_s.
cold_rt_run_s = let
    s = timed_stage(() -> rt_run(model), gpu = _HAS_CUDA)
    @printf("done (%.3f s, warmup)\n", s.time)
    s.time
end

rt_times  = Float64[]
rt_allocs = Float64[]
rt_gc     = Float64[]
R = nothing
last_rt_stat = nothing
for rep in 1:N_REPS
    global R, last_rt_stat
    print("[rt_run] rep $rep/$N_REPS … "); flush(stdout)
    s = timed_stage(() -> rt_run(model), gpu = _HAS_CUDA)
    R, _ = s.value
    last_rt_stat = s
    push!(rt_times,  s.time)
    push!(rt_allocs, gib(s.bytes))
    push!(rt_gc,     100 * s.gctime / max(s.time, eps()))
    @printf("%.3f s | alloc %.3f GiB | gc %.1f%%", s.time, gib(s.bytes),
            100 * s.gctime / max(s.time, eps()))
    if _HAS_CUDA
        @printf(" | GPU used Δ %.3f GiB (now %.2f GiB) | pool %.2f GiB",
                gib(s.gpu_used1 - s.gpu_used0), gib(s.gpu_used1),
                s.gpu_pool >= 0 ? gib(s.gpu_pool) : NaN)
    end
    println()
end

# pull host copy of R for the radiance-diff sanity check
Rh = Array(R)
maxabsR = maximum(abs, Rh)
peak_rss = gib(Sys.maxrss())

best_rt = minimum(rt_times)
mean_rt = sum(rt_times) / length(rt_times)

println("\n" * "-"^78)
@printf("SUMMARY  tag=%s\n", TAG)
@printf("  cold_model_build_s = %.3f   (first call: compile + first-touch)\n", cold_model_build_s)
@printf("  cold_rt_run_s      = %.3f   (first call: kernel compilation)\n", cold_rt_run_s)
@printf("  model_build_s      = %.3f\n", model_stat.time)
@printf("  model_alloc_GiB    = %.3f\n", gib(model_stat.bytes))
@printf("  rt_run_best_s      = %.3f   (reps: %s)\n", best_rt,
        join((@sprintf("%.3f", t) for t in rt_times), ", "))
@printf("  rt_run_mean_s      = %.3f\n", mean_rt)
@printf("  rt_run_alloc_GiB   = %.3f\n", minimum(rt_allocs))
@printf("  rt_run_gc_pct      = %.2f\n", minimum(rt_gc))
@printf("  peak_RSS_GiB       = %.3f\n", peak_rss)
@printf("  NquadN             = %d\n", NquadN)
@printf("  expandedZ_est_GiB  = %.3f\n", gib(zbytes))
@printf("  max_abs_R          = %.10e\n", maxabsR)
if _HAS_CUDA
    mi = CUDA.MemoryInfo()
    @printf("  gpu_used_now_GiB   = %.3f\n", gib(mi.total_bytes - mi.free_bytes))
    @printf("  gpu_pool_GiB       = %.3f\n",
            last_rt_stat.gpu_pool >= 0 ? gib(last_rt_stat.gpu_pool) : NaN)
end
println("-"^78)

# Optional machine-readable dump (one block per invocation; safe to append).
if !isempty(OUT_FILE)
    reps_csv = join(rt_times, ",")
    open(OUT_FILE, "a") do io
        println(io, "### tag=$(TAG) arch=$(ARCH) n_spec=$(N_SPEC) with_aerosol=$(WITH_AEROSOL)")
        println(io, "cold_model_build_s=$(cold_model_build_s)")
        println(io, "cold_rt_run_s=$(cold_rt_run_s)")
        println(io, "model_build_s=$(model_stat.time)")
        println(io, "model_alloc_bytes=$(model_stat.bytes)")
        println(io, "model_gc_s=$(model_stat.gctime)")
        println(io, "rt_run_best_s=$(best_rt)")
        println(io, "rt_run_mean_s=$(mean_rt)")
        println(io, "rt_run_reps_s=$(reps_csv)")
        println(io, "rt_run_alloc_bytes=$(last_rt_stat.bytes)")
        println(io, "rt_run_gc_pct=$(minimum(rt_gc))")
        println(io, "peak_rss_bytes=$(Sys.maxrss())")
        println(io, "nquadn=$(NquadN)")
        println(io, "expanded_z_bytes=$(zbytes)")
        println(io, "max_abs_R=$(maxabsR)")
        if _HAS_CUDA
            mi = CUDA.MemoryInfo()
            println(io, "gpu_used_now_bytes=$(mi.total_bytes - mi.free_bytes)")
            println(io, "gpu_pool_bytes=$(last_rt_stat.gpu_pool)")
        end
        println(io, "")
    end
    println("appended key=value summary to $(OUT_FILE)")
end

# Optional: dump the radiance vector for a branch-vs-main delta.
if !isempty(RAD_FILE)
    writedlm(RAD_FILE, vec(Rh))
    println("wrote radiance vector ($(length(vec(Rh))) values) to $(RAD_FILE)")
end

println("DONE tag=$(TAG)")
