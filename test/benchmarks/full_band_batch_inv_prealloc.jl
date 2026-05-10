#!/usr/bin/env julia
#
# Full-band O2-A CUDA benchmark for regular vs preallocated batched inverses.
#
# Run from repo/test:
#   julia --project=. benchmarks/full_band_batch_inv_prealloc.jl
#
# Useful overrides:
#   VSMARTMOM_BENCH_NSPEC=full|2048
#   VSMARTMOM_BENCH_NRUNS=2
#   VSMARTMOM_BENCH_VERBOSE_RT=1

using Logging
using Printf
using Statistics

using vSmartMOM
using vSmartMOM.CoreRT

mutable struct ModeSamples
    times::Vector{Float64}
    bytes::Vector{Int}
    result::Any
end

const CUDA_OK = try
    @eval using CUDA
    CUDA.functional()
catch
    false
end

const YAML_PATH = get(ENV, "VSMARTMOM_BENCH_YAML",
    normpath(joinpath(@__DIR__, "..", "test_parameters", "O2Parameters.yaml")))

@inline function env_enabled(name::AbstractString)
    return lowercase(get(ENV, name, "0")) in ("1", "true", "yes", "on")
end

function format_bytes(bytes)
    bytes < 1024 && return string(bytes, " B")
    units = ("KiB", "MiB", "GiB", "TiB")
    value = Float64(bytes)
    unit = "B"
    for u in units
        value /= 1024
        unit = u
        value < 1024 && break
    end
    return @sprintf("%.2f %s", value, unit)
end

function requested_arch()
    arch = uppercase(get(ENV, "VSMARTMOM_BENCH_ARCH", "GPU"))
    arch == "GPU" || error("This benchmark is GPU-only; CPU preallocation is intentionally not benchmarked.")
    CUDA_OK || error("CUDA is not functional, so the GPU batch_inv! benchmark cannot run")
    return vSmartMOM.Architectures.GPU(), "GPU"
end

function trim_to_first_band!(params)
    length(params.spec_bands) == 1 && return params

    params.spec_bands = [params.spec_bands[1]]
    params.brdf = [params.brdf[1]]
    if params.absorption_params !== nothing
        ap = params.absorption_params
        ap.fixed_molecules = ap.fixed_molecules[1:1]
        ap.variable_molecules = ap.variable_molecules[1:1]
        ap.luts = ap.luts[1:1]
        ap.h2o_lut = ap.h2o_lut[1:1]
    end
    return params
end

function subset_first_band!(params)
    full_band = collect(params.spec_bands[1])
    raw = lowercase(get(ENV, "VSMARTMOM_BENCH_NSPEC", "full"))
    nspec = raw in ("full", "all", "0") ? length(full_band) : parse(Int, raw)
    nspec = min(nspec, length(full_band))
    nspec <= 0 && error("VSMARTMOM_BENCH_NSPEC must be positive, full, all, or 0")

    if nspec < length(full_band)
        start = max(1, cld(length(full_band) - nspec, 2))
        stop = start + nspec - 1
        params.spec_bands = [full_band[start:stop]]
    else
        params.spec_bands = [full_band]
    end
    return params
end

function benchmark_params(arch)
    params = parameters_from_yaml(YAML_PATH)
    trim_to_first_band!(params)
    subset_first_band!(params)
    params.architecture = arch
    return params
end

function set_prealloc_mode!(arch_label::AbstractString, enabled::Bool)
    ENV["VSMARTMOM_CUDA_BATCH_INV_PREALLOC"] = arch_label == "GPU" && enabled ? "1" : "0"
end

function maybe_quiet(f)
    env_enabled("VSMARTMOM_BENCH_VERBOSE_RT") && return f()
    devnull_path = Sys.iswindows() ? "NUL" : "/dev/null"
    open(devnull_path, "w") do io
        redirect_stdout(io) do
            redirect_stderr(io) do
                Logging.with_logger(NullLogger()) do
                    return f()
                end
            end
        end
    end
end

function sync_if_needed(arch_label)
    if arch_label == "GPU"
        CUDA.synchronize()
    end
    return nothing
end

function timed_rt_run(model, arch_label)
    GC.gc()
    sync_if_needed(arch_label)
    stats = @timed maybe_quiet(() -> rt_run(model))
    sync_if_needed(arch_label)
    return stats
end

function timed_mode!(store, mode::Symbol, model, arch_label::String)
    set_prealloc_mode!(arch_label, mode == :prealloc)
    stats = timed_rt_run(model, arch_label)
    samples = store[mode]
    push!(samples.times, stats.time)
    push!(samples.bytes, stats.bytes)
    samples.result = stats.value
    return nothing
end

function run_comparison(arch, arch_label::String; nruns::Int)
    set_prealloc_mode!(arch_label, false)

    params = benchmark_params(arch)
    model_stats = @timed maybe_quiet(() -> model_from_parameters(params))
    model = model_stats.value
    nspec = length(model.atmosphere.spec_bands[1])

    # Warm both modes on the same model. The prealloc warm-up creates the
    # dimension-specific scratch closure before timed samples are collected.
    set_prealloc_mode!(arch_label, false)
    regular_warmup = timed_rt_run(model, arch_label)
    set_prealloc_mode!(arch_label, true)
    prealloc_warmup = timed_rt_run(model, arch_label)
    regular_ref = regular_warmup.value
    prealloc_ref = prealloc_warmup.value
    regular_warmup = nothing
    prealloc_warmup = nothing

    store = Dict{Symbol, ModeSamples}(
        :regular => ModeSamples(Float64[], Int[], regular_ref),
        :prealloc => ModeSamples(Float64[], Int[], prealloc_ref),
    )

    for i in 1:nruns
        order = isodd(i) ? (:regular, :prealloc) : (:prealloc, :regular)
        timed_mode!(store, order[1], model, arch_label)
        timed_mode!(store, order[2], model, arch_label)
    end

    function summarize(mode)
        times = store[mode].times
        bytes = store[mode].bytes
        result = store[mode].result
        return (
            mode = mode,
            nspec = nspec,
            median_rt = median(times),
            best_rt = minimum(times),
            median_bytes = round(Int, median(bytes)),
            best_bytes = minimum(bytes),
            R = Array(result[1]),
            T = Array(result[2]),
        )
    end

    return summarize(:regular), summarize(:prealloc), (
        nspec = nspec,
        model_time = model_stats.time,
        model_bytes = model_stats.bytes,
    )
end

function max_drift(test, ref)
    absdiff = abs.(test .- ref)
    scale = max.(abs.(ref), eps(eltype(ref)))
    return maximum(absdiff), maximum(absdiff ./ scale)
end

function print_result(row)
    @printf("%-10s  rt median %8.3fs  rt best %8.3fs  julia alloc median %10s  julia alloc best %10s\n",
        String(row.mode),
        row.median_rt,
        row.best_rt,
        format_bytes(row.median_bytes),
        format_bytes(row.best_bytes))
end

function main()
    arch, arch_label = requested_arch()
    nruns = parse(Int, get(ENV, "VSMARTMOM_BENCH_NRUNS",
        get(ENV, "VSMARTMOM_BENCH_NRUNC", "2")))

    println("="^86)
    println("vSmartMOM full-band batch_inv! preallocation benchmark")
    println("Julia:  $(VERSION)")
    println("YAML:   $(YAML_PATH)")
    println("Arch:   $(arch_label)")
    println("CUDA:   $(CUDA_OK)")
    println("Runs:   $(nruns) timed runs after one warmup per mode")
    println("NSpec:  $(get(ENV, "VSMARTMOM_BENCH_NSPEC", "full"))")
    println("="^86)

    regular, prealloc, build = run_comparison(arch, arch_label; nruns)

    @printf("\nModel build: %.3fs | Julia alloc: %s | nSpec: %d\n",
        build.model_time, format_bytes(build.model_bytes), regular.nspec)
    println("\nResults")
    print_result(regular)
    print_result(prealloc)

    r_abs, r_rel = max_drift(prealloc.R, regular.R)
    t_abs, t_rel = max_drift(prealloc.T, regular.T)
    speedup = regular.median_rt / prealloc.median_rt
    alloc_ratio = regular.median_bytes / max(prealloc.median_bytes, 1)

    println("\nComparison: prealloc vs regular")
    @printf("  nSpec:             %d\n", regular.nspec)
    @printf("  median speedup:    %.3fx\n", speedup)
    @printf("  Julia alloc ratio: %.3fx\n", alloc_ratio)
    @printf("  R max abs/rel:     %.6e / %.6e\n", r_abs, r_rel)
    @printf("  T max abs/rel:     %.6e / %.6e\n", t_abs, t_rel)
    println("  Note: @timed allocation counts are Julia host allocations, not CUDA device allocation totals.")
end

main()
