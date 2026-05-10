#!/usr/bin/env julia
# ==========================================================================
# Mie + CoreRT Performance Profiling Script
# ==========================================================================
# Compares forward, ForwardDiff AD, and analytic linearized Mie code, then
# benchmarks forward vs linearized RT on CPU (and GPU if available).
#
# Usage:
#   julia --project=. sandbox/emit/benchmark_mie_and_rt.jl
# ==========================================================================

using Revise
using vSmartMOM, vSmartMOM.CoreRT, vSmartMOM.Scattering
using vSmartMOM.InelasticScattering
using Distributions
using LinearAlgebra
using Printf
using TimerOutputs
using ForwardDiff, DiffResults

CUDA_OK = false
try
    using CUDA
    global CUDA_OK = CUDA.functional()
catch
end

# ─────────────────────────────────────────────────────────────────
# Section 1: Mie Benchmarks
# ─────────────────────────────────────────────────────────────────

const YAML_PATH = joinpath(@__DIR__, "ParamsEMIT.yaml")

function build_mie_model(; λ=0.55, nquad_radius=2500, r_max=50.0,
                          rₘ=0.1, σ_g=1.75, nᵣ=1.5, nᵢ=0.001)
    aerosol = Aerosol(LogNormal(log(rₘ), log(σ_g)), nᵣ, nᵢ)
    return make_mie_model(NAI2(), aerosol, λ, Stokes_IQU(), δBGE(5, 2.0),
                          r_max, nquad_radius)
end

function bench_mie_forward(model; nruns=3)
    compute_aerosol_optical_properties(model)  # warmup
    times = zeros(nruns)
    allocs = zeros(Int, nruns)
    for i in 1:nruns
        stats = @timed compute_aerosol_optical_properties(model)
        times[i] = stats.time
        allocs[i] = stats.bytes
    end
    return (time=minimum(times), alloc=minimum(allocs))
end

function bench_mie_ad(model; nruns=3)
    compute_aerosol_optical_properties(model; autodiff=true)  # warmup
    times = zeros(nruns)
    allocs = zeros(Int, nruns)
    for i in 1:nruns
        stats = @timed compute_aerosol_optical_properties(model; autodiff=true)
        times[i] = stats.time
        allocs[i] = stats.bytes
    end
    return (time=minimum(times), alloc=minimum(allocs))
end

function bench_mie_linearized(model; nruns=3)
    compute_aerosol_optical_properties(LinMode(), model)  # warmup
    times = zeros(nruns)
    allocs = zeros(Int, nruns)
    for i in 1:nruns
        stats = @timed compute_aerosol_optical_properties(LinMode(), model)
        times[i] = stats.time
        allocs[i] = stats.bytes
    end
    return (time=minimum(times), alloc=minimum(allocs))
end

function bench_mie_ref_extinction_lin(model; nruns=3)
    compute_ref_aerosol_extinction(LinMode(), model)  # warmup
    times = zeros(nruns)
    allocs = zeros(Int, nruns)
    for i in 1:nruns
        stats = @timed compute_ref_aerosol_extinction(LinMode(), model)
        times[i] = stats.time
        allocs[i] = stats.bytes
    end
    return (time=minimum(times), alloc=minimum(allocs))
end

function format_bytes(b)
    if b >= 1e9
        @sprintf("%.2f GiB", b / 1024^3)
    elseif b >= 1e6
        @sprintf("%.1f MiB", b / 1024^2)
    elseif b >= 1e3
        @sprintf("%.1f KiB", b / 1024)
    else
        @sprintf("%d B", b)
    end
end

# ─────────────────────────────────────────────────────────────────
# Section 2: RT Benchmarks
# ─────────────────────────────────────────────────────────────────

function bench_rt_forward(arch; yaml=YAML_PATH, nruns=2)
    params = parameters_from_yaml(yaml)
    params.architecture = arch
    params.spec_bands = [params.spec_bands[2]]
    params.brdf = [params.brdf[2]]
    params.absorption_params.molecules = [params.absorption_params.molecules[2]]
    model = model_from_parameters(params)
    reset_timer!()
    rt_run(model)  # warmup

    times = zeros(nruns)
    for i in 1:nruns
        reset_timer!()
        GC.gc()
        times[i] = @elapsed rt_run(model)
    end
    return (rt=minimum(times),)
end

function bench_rt_linearized(arch; yaml=YAML_PATH, nruns=2)
    params = parameters_from_yaml(yaml)
    params.architecture = arch
    params.spec_bands = [params.spec_bands[2]]
    params.brdf = [params.brdf[2]]
    params.absorption_params.molecules = [params.absorption_params.molecules[2]]
    model, lin_model = model_from_parameters(LinMode(), params)
    NAer  = length(params.scattering_params.rt_aerosols)
    NGas  = size(lin_model.τ̇_abs[1], 1)
    NSurf = 1
    Nparams = 7*NAer + NGas + NSurf
    reset_timer!()
    rt_run(model, lin_model, NAer, NGas, NSurf)  # warmup

    times = zeros(nruns)
    for i in 1:nruns
        reset_timer!()
        GC.gc()
        times[i] = @elapsed rt_run(model, lin_model, NAer, NGas, NSurf)
    end
    print_timer()
    return (rt=minimum(times), Nparams=Nparams, NAer=NAer, NGas=NGas, NSurf=NSurf)
end

# ─────────────────────────────────────────────────────────────────
# Section 3: Run everything
# ─────────────────────────────────────────────────────────────────

println("="^72)
println("  vSmartMOM Mie + RT Performance Profiling")
println("  Julia $(VERSION)")
println("  CUDA available: $CUDA_OK")
println("="^72)

# --- Mie Benchmarks ---
println("\n" * "─"^72)
println("  PART 1: Mie Scattering Performance (λ=0.55µm)")
println("─"^72)

mie_model = build_mie_model(nquad_radius=2500)

println("\nBenchmarking forward-only Mie...")
r_fwd = bench_mie_forward(mie_model)
@printf("  Forward:         %8.3fs | %s\n", r_fwd.time, format_bytes(r_fwd.alloc))

println("Benchmarking ForwardDiff AD Mie...")
r_ad = bench_mie_ad(mie_model)
@printf("  ForwardDiff AD:  %8.3fs | %s\n", r_ad.time, format_bytes(r_ad.alloc))

println("Benchmarking analytic linearized Mie...")
r_lin = bench_mie_linearized(mie_model)
@printf("  Analytic Lin:    %8.3fs | %s\n", r_lin.time, format_bytes(r_lin.alloc))

println("Benchmarking linearized ref extinction...")
r_ext = bench_mie_ref_extinction_lin(mie_model)
@printf("  Lin RefExtinct:  %8.3fs | %s\n", r_ext.time, format_bytes(r_ext.alloc))

println("\n  Overhead ratios (relative to forward):")
@printf("    ForwardDiff AD / Forward:  %6.1fx time, %6.1fx alloc\n",
        r_ad.time / r_fwd.time, r_ad.alloc / r_fwd.alloc)
@printf("    Analytic Lin   / Forward:  %6.1fx time, %6.1fx alloc\n",
        r_lin.time / r_fwd.time, r_lin.alloc / r_fwd.alloc)
@printf("    Analytic Lin   / AD:       %6.1fx time, %6.1fx alloc\n",
        r_lin.time / r_ad.time, r_lin.alloc / r_ad.alloc)

# Also test a smaller model to see scaling
println("\n  Scaling test (nquad_radius=500):")
mie_small = build_mie_model(nquad_radius=500)
r_fwd_s = bench_mie_forward(mie_small)
r_ad_s  = bench_mie_ad(mie_small)
r_lin_s = bench_mie_linearized(mie_small)
@printf("    Forward:       %8.3fs | %s\n", r_fwd_s.time, format_bytes(r_fwd_s.alloc))
@printf("    ForwardDiff:   %8.3fs | %s\n", r_ad_s.time, format_bytes(r_ad_s.alloc))
@printf("    Analytic Lin:  %8.3fs | %s\n", r_lin_s.time, format_bytes(r_lin_s.alloc))
@printf("    Lin/Fwd ratio: %6.1fx time, %6.1fx alloc\n",
        r_lin_s.time / r_fwd_s.time, r_lin_s.alloc / r_fwd_s.alloc)

# --- RT Benchmarks ---
println("\n" * "─"^72)
println("  PART 2: RT Kernel Performance (single band, O2-B)")
println("─"^72)

println("\nBenchmarking forward RT (CPU)...")
t_fwd_cpu = bench_rt_forward(vSmartMOM.Architectures.CPU())
@printf("  Forward  (CPU): %8.3fs\n", t_fwd_cpu.rt)

println("\nBenchmarking linearized RT (CPU)...")
t_lin_cpu = bench_rt_linearized(vSmartMOM.Architectures.CPU())
@printf("  Linearized (CPU): %8.3fs  (Nparams=%d)\n", t_lin_cpu.rt, t_lin_cpu.Nparams)
@printf("  Lin/Fwd overhead: %6.1fx\n", t_lin_cpu.rt / t_fwd_cpu.rt)

if CUDA_OK
    println("\nBenchmarking forward RT (GPU)...")
    t_fwd_gpu = bench_rt_forward(vSmartMOM.Architectures.GPU())
    @printf("  Forward  (GPU): %8.3fs\n", t_fwd_gpu.rt)

    println("\nBenchmarking linearized RT (GPU)...")
    t_lin_gpu = bench_rt_linearized(vSmartMOM.Architectures.GPU())
    @printf("  Linearized (GPU): %8.3fs  (Nparams=%d)\n", t_lin_gpu.rt, t_lin_gpu.Nparams)
    @printf("  Lin/Fwd overhead: %6.1fx\n", t_lin_gpu.rt / t_fwd_gpu.rt)

    println("\n  GPU speedup:")
    @printf("    Forward:    %5.1fx\n", t_fwd_cpu.rt / t_fwd_gpu.rt)
    @printf("    Linearized: %5.1fx\n", t_lin_cpu.rt / t_lin_gpu.rt)
end

# --- Summary ---
println("\n" * "="^72)
println("  SUMMARY TABLE")
println("="^72)
@printf("\n  %-25s %10s %15s %10s\n", "Component", "Time (s)", "Allocations", "Overhead")
println("  " * "─"^62)
@printf("  %-25s %10.3f %15s %10s\n", "Mie Forward", r_fwd.time, format_bytes(r_fwd.alloc), "1.0x")
@printf("  %-25s %10.3f %15s %9.1fx\n", "Mie ForwardDiff AD", r_ad.time, format_bytes(r_ad.alloc), r_ad.time/r_fwd.time)
@printf("  %-25s %10.3f %15s %9.1fx\n", "Mie Analytic Lin", r_lin.time, format_bytes(r_lin.alloc), r_lin.time/r_fwd.time)
println("  " * "─"^62)
@printf("  %-25s %10.3f %15s %10s\n", "RT Forward (CPU)", t_fwd_cpu.rt, "—", "1.0x")
@printf("  %-25s %10.3f %15s %9.1fx\n", "RT Linearized (CPU)", t_lin_cpu.rt, "—", t_lin_cpu.rt/t_fwd_cpu.rt)
if CUDA_OK
    @printf("  %-25s %10.3f %15s %10s\n", "RT Forward (GPU)", t_fwd_gpu.rt, "—", "1.0x")
    @printf("  %-25s %10.3f %15s %9.1fx\n", "RT Linearized (GPU)", t_lin_gpu.rt, "—", t_lin_gpu.rt/t_fwd_gpu.rt)
end
println("  " * "─"^62)

println("\n  Diagnosis:")
if r_lin.time / r_fwd.time > 10
    @printf("  [FIXABLE] Mie Lin/Fwd = %.0fx — dominated by temporary allocations\n",
            r_lin.time / r_fwd.time)
    @printf("            Mie Lin allocates %s vs Forward %s (%.0fx more)\n",
            format_bytes(r_lin.alloc), format_bytes(r_fwd.alloc),
            r_lin.alloc / r_fwd.alloc)
end
if t_lin_cpu.rt / t_fwd_cpu.rt > 3
    overhead = t_lin_cpu.rt / t_fwd_cpu.rt
    theoretical = 1 + t_lin_cpu.Nparams * 1.8
    if overhead > theoretical * 1.5
        @printf("  [FIXABLE] RT Lin/Fwd = %.1fx — exceeds theoretical %.1fx for Nparams=%d\n",
                overhead, theoretical, t_lin_cpu.Nparams)
    else
        @printf("  [EXPECTED] RT Lin/Fwd = %.1fx — consistent with theoretical %.1fx for Nparams=%d\n",
                overhead, theoretical, t_lin_cpu.Nparams)
    end
end

println("\n" * "="^72)
println("  Profiling complete.")
println("="^72)
