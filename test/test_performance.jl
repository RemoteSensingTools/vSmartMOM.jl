#!/usr/bin/env julia
# ==========================================================================
# vSmartMOM Performance Report
# ==========================================================================
# Compares forward vs linearized RT on CPU and GPU.
# Reports wall-clock timings, linearization overhead, GPU speedup,
# and per-component breakdown from TimerOutputs.
#
# Usage:
#   julia --project=. test/test_performance.jl
#
# This is a diagnostic script, not a @testset — not included in runtests.jl.
# ==========================================================================

using vSmartMOM, vSmartMOM.CoreRT
using Statistics
using Printf
using TimerOutputs

CUDA_OK = false
try
    using CUDA
    global CUDA_OK = CUDA.functional()
catch
end

const YAML_FAST = "test_parameters/JacobianTestFast.yaml"
const YAML_RAYLEIGH = "test_parameters/PureRayleighParameters.yaml"

# ---------------------------------------------------------------
# Benchmark helpers
# ---------------------------------------------------------------
function bench_forward(arch; yaml=YAML_RAYLEIGH, nruns=3)
    params = parameters_from_yaml(yaml)
    params.architecture = arch
    model = model_from_parameters(params)
    reset_timer!()
    rt_run(model)  # warmup

    t_build = zeros(nruns)
    t_rt    = zeros(nruns)
    for i in 1:nruns
        reset_timer!()
        t_build[i] = @elapsed begin
            model = model_from_parameters(params)
        end
        t_rt[i] = @elapsed rt_run(model)
    end
    return (model=median(t_build), rt=median(t_rt), total=median(t_build .+ t_rt))
end

function bench_linearized(arch; yaml=YAML_FAST, nruns=3)
    params = parameters_from_yaml(yaml)
    params.architecture = arch
    model, lin_model = model_from_parameters(LinMode(), params)
    NAer  = length(params.scattering_params.rt_aerosols)
    NGas  = size(lin_model.τ̇_abs[1], 1)
    NSurf = 1
    reset_timer!()
    rt_run(model, lin_model, NAer, NGas, NSurf)  # warmup

    t_build = zeros(nruns)
    t_rt    = zeros(nruns)
    for i in 1:nruns
        reset_timer!()
        t_build[i] = @elapsed begin
            model, lin_model = model_from_parameters(LinMode(), params)
        end
        t_rt[i] = @elapsed rt_run(model, lin_model, NAer, NGas, NSurf)
    end
    # Capture the last run's timer breakdown
    return (model=median(t_build), rt=median(t_rt), total=median(t_build .+ t_rt),
            NAer=NAer, NGas=NGas, NSurf=NSurf)
end

function bench_linearized_breakdown(arch; yaml=YAML_FAST)
    params = parameters_from_yaml(yaml)
    params.architecture = arch
    model, lin_model = model_from_parameters(LinMode(), params)
    NAer  = length(params.scattering_params.rt_aerosols)
    NGas  = size(lin_model.τ̇_abs[1], 1)
    NSurf = 1

    reset_timer!()
    rt_run(model, lin_model, NAer, NGas, NSurf)
    println("\n--- TimerOutputs Breakdown (Linearized, $arch) ---")
    print_timer()
    reset_timer!()
end

# ---------------------------------------------------------------
# Run benchmarks
# ---------------------------------------------------------------
nruns = 2

println("="^70)
println("  vSmartMOM Performance Report")
println("  Julia $(VERSION) | $(Sys.CPU_NAME)")
println("  CUDA available: $CUDA_OK")
println("  Each timing: median of $nruns runs (after 1 warmup)")
println("="^70)

# --- CPU Forward ---
print("\nForward  (CPU)... ")
t_fwd_cpu = bench_forward(vSmartMOM.Architectures.CPU(); nruns=nruns)
@printf("model: %.3fs | rt: %.3fs | total: %.3fs\n",
        t_fwd_cpu.model, t_fwd_cpu.rt, t_fwd_cpu.total)

# --- CPU Linearized ---
print("Linearized (CPU)... ")
t_lin_cpu = bench_linearized(vSmartMOM.Architectures.CPU(); nruns=nruns)
@printf("model: %.3fs | rt: %.3fs | total: %.3fs  (NAer=%d NGas=%d NSurf=%d)\n",
        t_lin_cpu.model, t_lin_cpu.rt, t_lin_cpu.total,
        t_lin_cpu.NAer, t_lin_cpu.NGas, t_lin_cpu.NSurf)

# --- CPU Breakdown ---
bench_linearized_breakdown(vSmartMOM.Architectures.CPU())

# --- GPU ---
if CUDA_OK
    print("\nForward  (GPU)... ")
    t_fwd_gpu = bench_forward(vSmartMOM.Architectures.GPU(); nruns=nruns)
    @printf("model: %.3fs | rt: %.3fs | total: %.3fs\n",
            t_fwd_gpu.model, t_fwd_gpu.rt, t_fwd_gpu.total)

    print("Linearized (GPU)... ")
    t_lin_gpu = bench_linearized(vSmartMOM.Architectures.GPU(); nruns=nruns)
    @printf("model: %.3fs | rt: %.3fs | total: %.3fs\n",
            t_lin_gpu.model, t_lin_gpu.rt, t_lin_gpu.total)

    bench_linearized_breakdown(vSmartMOM.Architectures.GPU())

    # --- Summary ---
    println("\n" * "="^70)
    println("  Summary")
    println("="^70)
    @printf("                      CPU          GPU        GPU Speedup\n")
    @printf("Forward  (rt only)  %7.3fs     %7.3fs       %5.1fx\n",
            t_fwd_cpu.rt, t_fwd_gpu.rt, t_fwd_cpu.rt / t_fwd_gpu.rt)
    @printf("Linearized (rt)     %7.3fs     %7.3fs       %5.1fx\n",
            t_lin_cpu.rt, t_lin_gpu.rt, t_lin_cpu.rt / t_lin_gpu.rt)
    @printf("\nLin/Fwd overhead    %7.1fx     %7.1fx\n",
            t_lin_cpu.rt / t_fwd_cpu.rt, t_lin_gpu.rt / t_fwd_gpu.rt)

    println("\nBottleneck analysis:")
    if t_lin_cpu.model > t_lin_cpu.rt
        println("  CPU: Model build (Mie) dominates — $(round(t_lin_cpu.model/t_lin_cpu.total*100, digits=1))% of total")
    else
        println("  CPU: RT kernel dominates — $(round(t_lin_cpu.rt/t_lin_cpu.total*100, digits=1))% of total")
    end
    if t_lin_gpu.model > t_lin_gpu.rt
        println("  GPU: Model build (Mie) dominates — $(round(t_lin_gpu.model/t_lin_gpu.total*100, digits=1))% of total")
    else
        println("  GPU: RT kernel dominates — $(round(t_lin_gpu.rt/t_lin_gpu.total*100, digits=1))% of total")
    end
    println("  Note: In inversions, Mie is cached — RT kernel time matters most.")
else
    println("\n" * "="^70)
    @printf("\nLin/Fwd overhead (CPU): %.1fx\n", t_lin_cpu.rt / t_fwd_cpu.rt)
    println("  No GPU available for comparison.")
end

println("\n" * "="^70)
println("  Performance report complete.")
println("="^70)
