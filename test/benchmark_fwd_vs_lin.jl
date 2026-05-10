#!/usr/bin/env julia
# ==========================================================================
# Benchmark: Forward-Only vs Linearized RT — Full Pipeline
# ==========================================================================
# Times the complete end-to-end pipeline:
#   model_from_parameters() + rt_run()
# for both forward and linearized modes, on CPU and GPU.
#
# O2 A-band config: ~6800 spectral points, IQU polarization, 1 aerosol.
#
# Usage:  julia --project=. test/benchmark_fwd_vs_lin.jl
# ==========================================================================

using vSmartMOM, vSmartMOM.CoreRT
using Statistics, Printf, TimerOutputs

CUDA_OK = false
try
    using CUDA
    global CUDA_OK = CUDA.functional()
catch
end

const YAML = joinpath(@__DIR__, "test_parameters", "O2Parameters.yaml")
const NRUNS = 5

println("="^72)
println("  Forward vs Linearized RT — Full Pipeline Benchmark")
println("  Julia ", VERSION, " | ", Sys.CPU_NAME)
println("  CUDA available: ", CUDA_OK)
println("  Config: O2Parameters.yaml | ", NRUNS, " timed runs (after warmup)")
println("="^72)

# ─────────────────────────────────────────────────────────────────────────
function run_benchmark(arch; yaml=YAML, nruns=NRUNS)
    arch_name = arch isa vSmartMOM.Architectures.GPU ? "GPU" : "CPU"

    println("\n", "="^72)
    println("  Architecture: ", arch_name)
    println("="^72)

    # ── Determine problem size ────────────────────────────────────────
    params0 = parameters_from_yaml(yaml)
    params0.architecture = arch
    m0 = model_from_parameters(params0)
    nSpec   = length(m0.params.spec_bands[1])
    nLayers = size(m0.τ_aer[1], 3)
    nStokes = length(params0.polarization_type.I₀)
    nQuad   = nStokes * (params0.l_trunc + 1)

    m0_lin, lm0 = model_from_parameters(LinMode(), params0)
    NAer  = length(params0.scattering_params.rt_aerosols)
    NGas  = size(lm0.τ̇_abs[1], 1)
    NSurf = 1
    Nparams = NAer * 7 + NGas + NSurf

    @printf("  nSpec=%d  nLayers=%d  nQuad=%d (nStokes=%d)  Nparams=%d\n",
            nSpec, nLayers, nQuad, nStokes, Nparams)
    @printf("  NAer=%d  NGas=%d  NSurf=%d  Matrix: %d x %d per lambda\n",
            NAer, NGas, NSurf, nQuad, nQuad)

    # ── Warmup (full pipeline, 2 runs each) ───────────────────────────
    println("\n  Warming up...")
    for _ in 1:2
        p = parameters_from_yaml(yaml); p.architecture = arch
        m = model_from_parameters(p); reset_timer!(); rt_run(m)
    end
    for _ in 1:2
        p = parameters_from_yaml(yaml); p.architecture = arch
        m, lm = model_from_parameters(LinMode(), p)
        na = length(p.scattering_params.rt_aerosols)
        ng = size(lm.τ̇_abs[1], 1)
        reset_timer!(); rt_run(m, lm, na, ng, 1)
    end

    # ── Forward-only: full pipeline ───────────────────────────────────
    println("\n", "-"^72)
    println("  Forward-only (", arch_name, ") -- full pipeline")
    println("-"^72)

    t_fwd_total = zeros(nruns)
    t_fwd_model = zeros(nruns)
    t_fwd_rt    = zeros(nruns)
    for i in 1:nruns
        p = parameters_from_yaml(yaml); p.architecture = arch
        reset_timer!()
        t_fwd_model[i] = @elapsed begin m = model_from_parameters(p) end
        t_fwd_rt[i]    = @elapsed begin rt_run(m) end
        t_fwd_total[i] = t_fwd_model[i] + t_fwd_rt[i]
    end

    # Show TimerOutputs for last forward run
    p = parameters_from_yaml(yaml); p.architecture = arch
    m = model_from_parameters(p)
    reset_timer!()
    rt_run(m)
    println("\n  Forward TimerOutputs (RT kernel only, last run):")
    print_timer()

    # ── Linearized: full pipeline ─────────────────────────────────────
    println("\n", "-"^72)
    println("  Linearized (", arch_name, ", Nparams=", Nparams, ") -- full pipeline")
    println("-"^72)

    t_lin_total = zeros(nruns)
    t_lin_model = zeros(nruns)
    t_lin_rt    = zeros(nruns)
    for i in 1:nruns
        p = parameters_from_yaml(yaml); p.architecture = arch
        reset_timer!()
        t_lin_model[i] = @elapsed begin m, lm = model_from_parameters(LinMode(), p) end
        na = length(p.scattering_params.rt_aerosols)
        ng = size(lm.τ̇_abs[1], 1)
        t_lin_rt[i]    = @elapsed begin rt_run(m, lm, na, ng, 1) end
        t_lin_total[i] = t_lin_model[i] + t_lin_rt[i]
    end

    # Show TimerOutputs for last linearized run
    p = parameters_from_yaml(yaml); p.architecture = arch
    m, lm = model_from_parameters(LinMode(), p)
    na = length(p.scattering_params.rt_aerosols)
    ng = size(lm.τ̇_abs[1], 1)
    reset_timer!()
    rt_run(m, lm, na, ng, 1)
    println("\n  Linearized TimerOutputs (RT kernel only, last run):")
    print_timer()

    # ── Per-architecture results ──────────────────────────────────────
    fwd_med   = median(t_fwd_total)
    lin_med   = median(t_lin_total)
    fwd_rt    = median(t_fwd_rt)
    lin_rt    = median(t_lin_rt)
    fwd_model = median(t_fwd_model)
    lin_model = median(t_lin_model)

    println("\n", "-"^72)
    println("  ", arch_name, " Results (median of ", nruns, " runs)")
    println("-"^72)
    @printf("                   model_build     rt_run       total\n")
    @printf("  Forward:        %8.3fs    %8.3fs    %8.3fs\n",
            fwd_model, fwd_rt, fwd_med)
    @printf("  Linearized:     %8.3fs    %8.3fs    %8.3fs\n",
            lin_model, lin_rt, lin_med)
    @printf("  Lin/Fwd ratio:  %8.2fx    %8.2fx    %8.2fx\n",
            lin_model/fwd_model, lin_rt/fwd_rt, lin_med/fwd_med)
    println()
    println("  Per-run totals (ms):")
    @printf("    Forward:    %s\n",
            join([@sprintf("%.1f", t*1000) for t in t_fwd_total], ", "))
    @printf("    Linearized: %s\n",
            join([@sprintf("%.1f", t*1000) for t in t_lin_total], ", "))

    return (arch=arch_name,
            fwd_total=fwd_med, lin_total=lin_med,
            fwd_rt=fwd_rt, lin_rt=lin_rt,
            fwd_model=fwd_model, lin_model=lin_model,
            overhead_total=lin_med/fwd_med, overhead_rt=lin_rt/fwd_rt,
            nSpec=nSpec, Nparams=Nparams)
end

# ─────────────────────────────────────────────────────────────────────────
# Run benchmarks
# ─────────────────────────────────────────────────────────────────────────
results = []
push!(results, run_benchmark(vSmartMOM.Architectures.CPU()))
if CUDA_OK
    push!(results, run_benchmark(vSmartMOM.Architectures.GPU()))
end

# ─────────────────────────────────────────────────────────────────────────
# Final comparison
# ─────────────────────────────────────────────────────────────────────────
println("\n", "="^72)
println("  FINAL COMPARISON  (O2 A-band, nSpec=", results[1].nSpec,
        ", Nparams=", results[1].Nparams, ")")
println("="^72)

@printf("\n  %-6s  %12s  %12s  %12s  %12s\n",
        "", "Fwd total", "Lin total", "Overhead", "RT overhead")
@printf("  %-6s  %12s  %12s  %12s  %12s\n",
        "------", "------------", "------------", "------------", "------------")
for r in results
    @printf("  %-6s  %9.1f ms  %9.1f ms  %9.2fx      %9.2fx\n",
            r.arch, r.fwd_total*1000, r.lin_total*1000,
            r.overhead_total, r.overhead_rt)
end

if length(results) == 2
    cpu, gpu = results[1], results[2]
    @printf("\n  GPU speedup -- total:  %.1fx (fwd), %.1fx (lin)\n",
            cpu.fwd_total/gpu.fwd_total, cpu.lin_total/gpu.lin_total)
    @printf("  GPU speedup -- RT:     %.1fx (fwd), %.1fx (lin)\n",
            cpu.fwd_rt/gpu.fwd_rt, cpu.lin_rt/gpu.lin_rt)
end
println("="^72)
