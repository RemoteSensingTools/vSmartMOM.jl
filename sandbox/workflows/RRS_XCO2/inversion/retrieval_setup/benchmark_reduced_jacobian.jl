#!/usr/bin/env julia

"""
Benchmark the retrieval-selected `OCO_RRS_synth` Jacobian against the
historical full physical Jacobian on the same three-band, 16-layer scene.

The reported model-build time includes absorption, aerosol optics, and their
upstream tangents. The per-band RT time includes compact/full tangent
propagation through elemental, doubling, interaction, and surface/source
coupling. Julia compilation is excluded by one warm-up build and one warm-up
RT call of each path.

Environment controls:

- `JACOBIAN_BENCH_ARCH=GPU|CPU` (default `GPU`)
- `CUDA_DEVICE` (default `1`)
- `JACOBIAN_BENCH_FLOAT_TYPE=Float32|Float64` (default `Float32`)
- `JACOBIAN_BENCH_NSPEC` (default `256`; `0` keeps every 0.1 cm^-1 point)
- `JACOBIAN_BENCH_MAX_M` (optional diagnostic Fourier-order cap)
- `JACOBIAN_BENCH_BUILD_REPS` (default `2`)
- `JACOBIAN_BENCH_RT_REPS` (default `3`)
- `JACOBIAN_BENCH_OUT` (default `jacobian_timing_latest.dat` beside this file)
"""

include(joinpath(@__DIR__, "..", "..", "scripts", "common.jl"))
using .RRSXCO2Common
using Dates
using Printf
using Sockets: gethostname
using Statistics
using vSmartMOM
using vSmartMOM.CoreRT
using YAML

const USE_GPU = uppercase(get(ENV, "JACOBIAN_BENCH_ARCH", "GPU")) == "GPU"
const DEVICE = parse(Int, get(ENV, "CUDA_DEVICE", "1"))
const FLOAT_NAME = get(ENV, "JACOBIAN_BENCH_FLOAT_TYPE", "Float32")
const NSPEC = parse(Int, get(ENV, "JACOBIAN_BENCH_NSPEC", "256"))
const MAX_M = haskey(ENV, "JACOBIAN_BENCH_MAX_M") ?
    parse(Int, ENV["JACOBIAN_BENCH_MAX_M"]) : nothing
const BUILD_REPS = parse(Int, get(ENV, "JACOBIAN_BENCH_BUILD_REPS", "2"))
const RT_REPS = parse(Int, get(ENV, "JACOBIAN_BENCH_RT_REPS", "3"))
const OUT = get(ENV, "JACOBIAN_BENCH_OUT",
    joinpath(@__DIR__, "jacobian_timing_latest.dat"))

if USE_GPU
    @eval using CUDA
    CUDA.functional() || error("CUDA is not functional on this host")
    CUDA.device!(DEVICE)
end

synchronize_backend() = USE_GPU ? CUDA.synchronize() : nothing
function release_backend!()
    GC.gc(true)
    USE_GPU && CUDA.reclaim()
    return nothing
end

function timed(f)
    synchronize_backend()
    value = nothing
    elapsed = @elapsed begin
        value = f()
        synchronize_backend()
    end
    return value, elapsed
end

function truncate_profile!(params, psurf=1000.0)
    p0, T0, q0 = Float64.(params.p), Float64.(params.T), Float64.(params.q)
    k = searchsortedlast(p0, psurf - eps(psurf))
    pcenter = (p0[1:end-1] .+ p0[2:end]) ./ 2
    function at_pressure(values, pressure)
        hi = searchsortedfirst(pcenter, pressure)
        hi <= 1 && return values[1]
        hi > length(values) && return values[end]
        lo = hi - 1
        weight = (log(pressure) - log(pcenter[lo])) /
                 (log(pcenter[hi]) - log(pcenter[lo]))
        return values[lo] + weight * (values[hi] - values[lo])
    end
    FT = params.float_type
    params.p = FT.(vcat(p0[1:k], psurf))
    params.T = FT.(vcat(T0[1:k-1], at_pressure(T0, psurf)))
    params.q = FT.(vcat(q0[1:k-1], at_pressure(q0, psurf)))
    return params
end

function center_subset(values, n)
    all_values = collect(values)
    (n <= 0 || n >= length(all_values)) && return all_values
    first_index = (length(all_values) - n) ÷ 2 + 1
    return all_values[first_index:(first_index + n - 1)]
end

function benchmark_parameters()
    # Select Float32 before parsing: vSmartMOM_Parameters is parameterized by
    # the parsed precision, so assigning `params.float_type` afterwards would
    # leave nested aerosol/surface objects in the old precision.
    configure_luts!()
    input = YAML.load_file(CONFIG)
    rt = input["radiative_transfer"]
    FLOAT_NAME in ("Float32", "Float64") || error(
        "JACOBIAN_BENCH_FLOAT_TYPE must be Float32 or Float64")
    rt["float_type"] = FLOAT_NAME
    rt["architecture"] = USE_GPU ? "GPU()" : "CPU()"
    params = parameters_from_dict(input)
    FT = params.float_type
    params.sza = FT(30)
    params.vza = FT[0]
    params.vaz = FT[0]
    params.profile_reduction_n = 16
    truncate_profile!(params)
    params.spec_bands = [FT.(center_subset(band, NSPEC))
                         for band in params.spec_bands]
    MAX_M === nothing || (params.m_max_override = MAX_M)
    return params
end

full_build(params) = model_from_parameters(
    LinMode(), deepcopy(params); external_solar=true)
reduced_build(params) = model_from_parameters(
    OCO_RRS_synth(), deepcopy(params); external_solar=true)

function full_rt(model, lin_model, params, i_band)
    sources = sources_for_band(params, i_band)
    return rt_run_lin(model, lin_model,
        length(params.scattering_params.rt_aerosols),
        size(lin_model.τ̇_abs[i_band], 1), 3;
        i_band, sources)
end

function reduced_rt(model, lin_model, params, i_band)
    return rt_run_lin(model, lin_model; i_band,
        sources=sources_for_band(params, i_band))
end

function selected_full_columns(full_layout, active_layout)
    columns = copy(native_layer_columns(active_layout))
    append!(columns, surface_range(full_layout))
    if !isempty(sif_range(active_layout))
        append!(columns, sif_range(full_layout))
    end
    return columns
end

function comparison_metrics(full, reduced)
    columns = selected_full_columns(full.layout, reduced.layout)
    full_forward = Array(full.toa)
    reduced_forward = Array(reduced.toa)
    full_selected = Array(full.toa_jacobian[:, :, :, columns])
    reduced_jacobian = Array(reduced.toa_jacobian)
    forward_finite = isfinite.(full_forward) .& isfinite.(reduced_forward)
    jacobian_finite = isfinite.(full_selected) .& isfinite.(reduced_jacobian)
    forward_abs = all(forward_finite) ?
        maximum(abs.(full_forward .- reduced_forward)) : Inf
    jacobian_abs = all(jacobian_finite) ?
        maximum(abs.(full_selected .- reduced_jacobian)) : Inf
    scale = all(jacobian_finite) ? maximum(abs, full_selected) : Inf
    jacobian_rel = iszero(scale) ? jacobian_abs : jacobian_abs / scale
    full_nonfinite_by_column = [count(!isfinite, full_selected[:, :, :, j])
                                for j in axes(full_selected, 4)]
    reduced_nonfinite_by_column = [count(!isfinite, reduced_jacobian[:, :, :, j])
                                   for j in axes(reduced_jacobian, 4)]
    return (; columns, forward_abs, jacobian_abs, jacobian_rel,
            full_nonfinite_by_column, reduced_nonfinite_by_column)
end

function main()
    NSPEC >= 0 || error("JACOBIAN_BENCH_NSPEC must be nonnegative")
    BUILD_REPS > 0 || error("JACOBIAN_BENCH_BUILD_REPS must be positive")
    RT_REPS > 0 || error("JACOBIAN_BENCH_RT_REPS must be positive")

    params = benchmark_parameters()
    architecture = USE_GPU ? "GPU" : "CPU"
    hardware = USE_GPU ? string(CUDA.name(CUDA.device())) : Sys.CPU_NAME
    println("="^78)
    println("Full vs OCO_RRS_synth Jacobian benchmark")
    println("host=$(gethostname()) architecture=$architecture hardware=$hardware")
    println("layers=16 spectra=$(length.(params.spec_bands)) float_type=$(params.float_type)")
    println("build_repetitions=$BUILD_REPS rt_repetitions=$RT_REPS")
    println("="^78)

    # Compile and populate file/driver caches before collecting build times.
    println("Warming model construction ...")
    warm_full, _ = timed(() -> full_build(params))
    warm_reduced, _ = timed(() -> reduced_build(params))
    warm_full = warm_reduced = nothing
    release_backend!()

    full_build_times = Float64[]
    reduced_build_times = Float64[]
    full_pair = reduced_pair = nothing
    for repetition in 1:BUILD_REPS
        full_pair, full_time = timed(() -> full_build(params))
        push!(full_build_times, full_time)
        full_pair = nothing
        release_backend!()
        reduced_pair, reduced_time = timed(() -> reduced_build(params))
        push!(reduced_build_times, reduced_time)
        reduced_pair = nothing
        release_backend!()
        @printf("build repetition %d: full=%.3f s reduced=%.3f s speedup=%.2fx\n",
                repetition, full_time, reduced_time, full_time / reduced_time)
    end

    # Retain one pair of models for all per-band timings.
    full_pair, _ = timed(() -> full_build(params))
    reduced_pair, _ = timed(() -> reduced_build(params))
    full_model, full_lin = full_pair
    reduced_model, reduced_lin = reduced_pair
    n_global(reduced_lin.plan) == 30 || error(
        "Expected the 16-layer OCO plan to contain 30 global parameters")

    rows = NamedTuple[]
    for (i_band, band_name) in enumerate(BAND_NAMES)
        println("Warming $band_name RT paths ...")
        full_result, _ = timed(() -> full_rt(
            full_model, full_lin, params, i_band))
        reduced_result, _ = timed(() -> reduced_rt(
            reduced_model, reduced_lin, params, i_band))

        full_times = Float64[]
        reduced_times = Float64[]
        for repetition in 1:RT_REPS
            full_result, full_time = timed(() -> full_rt(
                full_model, full_lin, params, i_band))
            reduced_result, reduced_time = timed(() -> reduced_rt(
                reduced_model, reduced_lin, params, i_band))
            push!(full_times, full_time)
            push!(reduced_times, reduced_time)
            @printf("RT %-10s repetition %d: full=%.3f s reduced=%.3f s speedup=%.2fx\n",
                    band_name, repetition, full_time, reduced_time,
                    full_time / reduced_time)
        end

        metrics = comparison_metrics(full_result, reduced_result)
        if any(>(0), metrics.full_nonfinite_by_column) ||
           any(>(0), metrics.reduced_nonfinite_by_column)
            println(stderr, "$band_name non-finite selected Jacobian entries")
            for j in eachindex(metrics.columns)
                nf = metrics.full_nonfinite_by_column[j]
                nr = metrics.reduced_nonfinite_by_column[j]
                (nf > 0 || nr > 0) && println(stderr,
                    "  local=$j full_column=$(metrics.columns[j]) " *
                    "full_nonfinite=$nf reduced_nonfinite=$nr")
            end
            error("$band_name contains non-finite selected Jacobian entries")
        end
        metrics.forward_abs <= 2f-6 || error(
            "$band_name forward mismatch: $(metrics.forward_abs)")
        metrics.jacobian_abs <= 2f-5 || error(
            "$band_name selected Jacobian mismatch: $(metrics.jacobian_abs)")
        full_count = n_total(full_result.layout)
        reduced_count = n_total(reduced_result.layout)
        full_median = median(full_times)
        reduced_median = median(reduced_times)
        push!(rows, (band=band_name,
            nspec=length(params.spec_bands[i_band]),
            full_parameters=full_count,
            reduced_parameters=reduced_count,
            full_seconds=full_median,
            reduced_seconds=reduced_median,
            speedup=full_median / reduced_median,
            forward_max_abs=metrics.forward_abs,
            jacobian_max_abs=metrics.jacobian_abs,
            jacobian_scaled_error=metrics.jacobian_rel))
    end

    full_build_median = median(full_build_times)
    reduced_build_median = median(reduced_build_times)
    mkpath(dirname(OUT))
    open(OUT, "w") do io
        println(io, "# Full physical Jacobian versus OCO_RRS_synth compact Jacobian")
        println(io, "# created $(now())")
        println(io, "# host $(gethostname())")
        println(io, "# architecture $architecture")
        println(io, "# hardware $hardware")
        println(io, "# float_type $(params.float_type)")
        println(io, "# layers 16")
        println(io, "# nstreams $(params.nstreams)")
        println(io, "# truncation $(params.truncation)")
        println(io, "# spectral_points $(join(length.(params.spec_bands), ','))")
        @printf(io, "# full_build_median_s %.9f\n", full_build_median)
        @printf(io, "# reduced_build_median_s %.9f\n", reduced_build_median)
        @printf(io, "# build_speedup %.9f\n",
                full_build_median / reduced_build_median)
        println(io, "band nspec full_parameters reduced_parameters " *
                    "full_rt_s reduced_rt_s speedup forward_max_abs " *
                    "jacobian_max_abs jacobian_scaled_error")
        for row in rows
            @printf(io, "%s %d %d %d %.9f %.9f %.9f %.9e %.9e %.9e\n",
                row.band, row.nspec, row.full_parameters,
                row.reduced_parameters, row.full_seconds,
                row.reduced_seconds, row.speedup, row.forward_max_abs,
                row.jacobian_max_abs, row.jacobian_scaled_error)
        end
    end

    println("-"^78)
    @printf("model build median: full=%.3f s reduced=%.3f s speedup=%.2fx\n",
            full_build_median, reduced_build_median,
            full_build_median / reduced_build_median)
    for row in rows
        @printf("RT %-10s: %2d -> %2d params, full=%.3f s reduced=%.3f s speedup=%.2fx\n",
                row.band, row.full_parameters, row.reduced_parameters,
                row.full_seconds, row.reduced_seconds, row.speedup)
    end
    println("results: $OUT")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
