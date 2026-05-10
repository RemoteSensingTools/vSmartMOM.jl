# Phase 2 benchmark harness — metrics collection helpers
# =========================================================================
#
# Wraps a single forward-RT run (or any zero-arg callable) in the
# instrumentation required by the v2 plan's Phase 2 regression harness:
#
#  * Wall-clock distribution (N runs, median + MAD, outlier trim).
#  * CPU allocation count (per-run total allocs + bytes).
#  * GPU used-memory delta, when CUDA is functional.
#  * GPU peak during the run, via CUDA's memory statistics.
#  * Output arrays (I/Q/U and ieR/ieT) saved to disk for
#    physics-parity reporting in report.jl.
#
# All timings are done AFTER a warm-up call so first-run compilation
# doesn't dominate the median. `capture_run` returns a NamedTuple
# that `run_benchmarks.jl` serializes to JSON.

using Statistics
using StatsBase: mad
using JLD2

"""
    warmup_then_time(f; n=5)

Discard the first call result (compile / GPU-init amortization), then run
`f` `n` times and return a vector of per-run wall-clock in seconds.
"""
function warmup_then_time(f; n::Integer=5)
    f()                              # warm-up; discard
    times = Float64[]
    for _ in 1:n
        t0 = time()
        f()
        push!(times, time() - t0)
    end
    return times
end

"""
    trim_and_summarize(times)

Return a NamedTuple with median, MAD, min, max after dropping the single
longest timing (simple outlier trim for MCMC-style GC-spike resilience).
"""
function trim_and_summarize(times::AbstractVector{<:Real})
    length(times) < 3 && return (median=median(times), mad=mad(times, normalize=false),
                                 min=minimum(times), max=maximum(times), n=length(times))
    trimmed = sort(times)[1:end-1]
    return (median=median(trimmed), mad=mad(trimmed, normalize=false),
            min=minimum(trimmed), max=maximum(trimmed), n=length(trimmed))
end

"""
    capture_alloc_count(f)

Single call to `f`, capturing (allocs, bytes_allocated, gctime_s) via @timed.
Returns a NamedTuple with the three fields.
"""
function capture_alloc_count(f)
    stats = @timed f()
    return (value=stats.value, allocs=stats.bytes, gctime=stats.gctime)
end

"""
    capture_gpu_memory_delta(f)

Run `f` once and measure the GPU used-memory delta (after - before) via
`CUDA.Mem.info()`. When CUDA is unavailable or non-functional, returns
`(used_before=nothing, used_after=nothing, delta=nothing, gpu="none")`.
"""
function capture_gpu_memory_delta(f)
    gpu_info = try
        @eval using CUDA
        if Main.CUDA.functional()
            Main.CUDA.synchronize()
            free_before, total = Main.CUDA.Mem.info()
            used_before = total - free_before
            f()
            Main.CUDA.synchronize()
            free_after, _ = Main.CUDA.Mem.info()
            used_after = total - free_after
            return (used_before=used_before, used_after=used_after,
                    delta=used_after - used_before,
                    gpu=string(Main.CUDA.name(Main.CUDA.device())))
        end
    catch
        nothing
    end
    if gpu_info === nothing
        f()  # still run once so semantics match
        return (used_before=nothing, used_after=nothing, delta=nothing, gpu="none")
    end
    return gpu_info
end

"""
    hw_identifier()

Capture a terse hardware string (CPU model fragment + GPU name if any)
for the output JSON. Re-capture flagged when hardware changes (v2 plan Phase 2).
"""
function hw_identifier()
    cpu = try
        Sys.cpu_info()[1].model
    catch
        "unknown-cpu"
    end
    gpu = try
        @eval using CUDA
        Main.CUDA.functional() ? string(Main.CUDA.name(Main.CUDA.device())) : "no-cuda"
    catch
        "no-cuda"
    end
    return (cpu=cpu, gpu=gpu, julia_version=string(VERSION))
end

"""
    capture_run(name::AbstractString, run_fn::Function, save_dir::AbstractString;
                n::Int=5, output_keys=(:R_SFI, :T_SFI, :ieR_SFI, :ieT_SFI))

Top-level per-scenario capture. Calls `run_fn` once to obtain output arrays
(saved to `<save_dir>/<name>_output.jld2`), then times it with a warm-up +
`n` timed runs, captures allocation counts and GPU memory delta.

`run_fn` must be a no-argument function that returns a NamedTuple or Dict
mapping symbol names (e.g. `:R_SFI`) to arrays. Keys in `output_keys` are
the subset extracted into the saved JLD2. Missing keys are silently skipped.

Returns a NamedTuple suitable for JSON serialization (all numeric fields).
"""
function capture_run(name::AbstractString, run_fn::Function, save_dir::AbstractString;
                     n::Integer=5, output_keys=(:R_SFI, :T_SFI, :ieR_SFI, :ieT_SFI,
                                                :hem_R, :hem_T))
    mkpath(save_dir)

    # Combined capture run: GPU memory delta + outputs + allocation stats
    # from a single invocation. Avoids redundant ~3-min runs on RRS configs.
    local gpu_info
    local alloc_stats
    begin
        before_free = nothing
        total = nothing
        has_cuda = try
            @eval using CUDA
            Main.CUDA.functional()
        catch
            false
        end
        if has_cuda
            Main.CUDA.synchronize()
            before_free, total = Main.CUDA.Mem.info()
        end
        alloc_stats = capture_alloc_count(run_fn)
        if has_cuda
            Main.CUDA.synchronize()
            after_free, _ = Main.CUDA.Mem.info()
            used_before = total - before_free
            used_after  = total - after_free
            gpu_info = (used_before=used_before, used_after=used_after,
                        delta=used_after - used_before,
                        gpu=string(Main.CUDA.name(Main.CUDA.device())))
        else
            gpu_info = (used_before=nothing, used_after=nothing,
                        delta=nothing, gpu="none")
        end
    end
    outputs = alloc_stats.value

    # Persist output arrays for downstream report.jl comparison.
    out_path = joinpath(save_dir, name * "_output.jld2")
    outdict = Dict{String,Any}()
    for k in output_keys
        outputs isa NamedTuple && haskey(outputs, k) &&
            (outdict[string(k)] = outputs[k])
        outputs isa AbstractDict && haskey(outputs, k) &&
            (outdict[string(k)] = outputs[k])
    end
    isempty(outdict) || save(out_path, outdict)

    # Wall-clock distribution (warmed — first call above already compiled).
    times = warmup_then_time(run_fn; n=n)
    wc = trim_and_summarize(times)

    return (
        name = name,
        wall_median_s = wc.median,
        wall_mad_s    = wc.mad,
        wall_min_s    = wc.min,
        wall_max_s    = wc.max,
        wall_n        = wc.n,
        cpu_alloc_bytes = alloc_stats.allocs,
        cpu_gctime_s    = alloc_stats.gctime,
        gpu_used_before = gpu_info.used_before,
        gpu_used_after  = gpu_info.used_after,
        gpu_used_delta  = gpu_info.delta,
        gpu_name        = gpu_info.gpu,
        output_path     = isempty(outdict) ? nothing : out_path,
        output_keys     = collect(keys(outdict)),
        hw              = hw_identifier(),
    )
end
