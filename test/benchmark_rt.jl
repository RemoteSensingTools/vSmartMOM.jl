#!/usr/bin/env julia
# =================================================================
# vSmartMOM RT Benchmark Script
#
# Measures wall-clock time for the main RT code paths:
#   1. Forward RT (noRS)        — Pure Rayleigh, CPU
#   2. Forward RT (Raman/RRS)   — O2-A band, CPU
#   3. Linearized RT            — Forward + Jacobians, CPU
#   4. GPU variants             — conditional on CUDA
#
# Usage:
#   julia --project=. test/benchmark_rt.jl
# =================================================================

using vSmartMOM, vSmartMOM.CoreRT
using vSmartMOM.InelasticScattering
using Statistics
using Printf

const CAN_USE_GPU = try
    using CUDA
    CUDA.functional()
catch
    false
end

function benchmark_forward_noRS(arch; yaml="test/test_parameters/PureRayleighParameters.yaml", nruns=3)
    params = parameters_from_yaml(yaml)
    params.architecture = arch

    # Warmup
    model = model_from_parameters(params)
    rt_run(model)

    # Timed runs
    t_model = zeros(nruns)
    t_rt    = zeros(nruns)
    for i in 1:nruns
        t_model[i] = @elapsed model = model_from_parameters(params)
        t_rt[i]    = @elapsed rt_run(model)
    end
    return (model_build=median(t_model), rt_run=median(t_rt))
end

function benchmark_raman(arch; yaml="test/test_parameters/O2Parameters_GPU.yaml", nruns=3)
    params = parameters_from_yaml(yaml)
    params.architecture = arch
    FT = Float64

    model = model_from_parameters(params)
    iBand = 1
    ν = model.params.spec_bands[iBand]
    nSpec = length(ν)
    ν̃ = mean(ν)

    effT = (model.profile.vcd_dry' * model.profile.T) / sum(model.profile.vcd_dry)
    n2, o2 = InelasticScattering.getRamanAtmoConstants(ν̃, effT)

    nPol = model.params.polarization_type.n
    F₀ = zeros(FT, nPol, nSpec); F₀[1, :] .= 1.0
    SIF₀ = zeros(FT, nPol, nSpec)

    RS_type = InelasticScattering.RRS(
        n2=n2, o2=o2,
        greek_raman=InelasticScattering.GreekCoefs(
            [FT(1)], [FT(1)], [FT(1)], [FT(1)], [FT(1)], [FT(1)]),
        fscattRayl=[FT(1)], ϖ_Cabannes=[FT(1)],
        ϖ_λ₁λ₀=zeros(FT,1), i_λ₁λ₀=zeros(Int,1),
        Z⁻⁺_λ₁λ₀=zeros(FT,1,1), Z⁺⁺_λ₁λ₀=zeros(FT,1,1),
        i_ref=argmin(abs.(ν .- ν̃)), n_Raman=0,
        F₀=F₀, SIF₀=SIF₀)
    CoreRT.getRamanSSProp!(RS_type, 1e7/ν̃, ν)

    # Warmup
    CoreRT.rt_run_test(RS_type, model, iBand)

    t_rt = zeros(nruns)
    for i in 1:nruns
        t_rt[i] = @elapsed CoreRT.rt_run_test(RS_type, model, iBand)
    end
    return (rt_run=median(t_rt), nSpec=nSpec)
end

function benchmark_linearized(arch; yaml="test/test_parameters/ParamsEMIT_fast.yaml", nruns=3)
    params = parameters_from_yaml(yaml)
    params.architecture = arch

    # Warmup
    model, lin_model = model_from_parameters(LinMode(), params)
    NAer = length(params.scattering_params.rt_aerosols)
    NGas = size(lin_model.τ̇_abs[1], 1)
    NSurf = 1
    rt_run(model, lin_model, NAer, NGas, NSurf)

    t_model = zeros(nruns)
    t_rt    = zeros(nruns)
    for i in 1:nruns
        t_model[i] = @elapsed begin
            model, lin_model = model_from_parameters(LinMode(), params)
        end
        t_rt[i] = @elapsed rt_run(model, lin_model, NAer, NGas, NSurf)
    end
    return (model_build=median(t_model), rt_run=median(t_rt),
            NAer=NAer, NGas=NGas, NSurf=NSurf)
end

# =================================================================
# Run benchmarks
# =================================================================
println("="^70)
println("  vSmartMOM RT Benchmarks")
println("  Julia $(VERSION) | $(Sys.CPU_NAME)")
println("  CUDA available: $CAN_USE_GPU")
println("="^70)

nruns = 3
println("\nEach timing is the median of $nruns runs (after 1 warmup).\n")

# --- Forward noRS (CPU) ---
print("Forward noRS (CPU)... ")
t = benchmark_forward_noRS(vSmartMOM.Architectures.CPU(); nruns=nruns)
@printf("model: %.3fs  |  rt_run: %.3fs  |  total: %.3fs\n",
        t.model_build, t.rt_run, t.model_build + t.rt_run)

# --- Raman RRS (CPU) ---
print("Raman RRS    (CPU)... ")
t_r = benchmark_raman(vSmartMOM.Architectures.CPU(); nruns=nruns)
@printf("rt_run: %.3fs  (nSpec=%d)\n", t_r.rt_run, t_r.nSpec)

# --- Linearized RT (CPU) ---
print("Linearized   (CPU)... ")
t_l = benchmark_linearized(vSmartMOM.Architectures.CPU(); nruns=nruns)
@printf("model: %.3fs  |  rt_run: %.3fs  |  total: %.3fs  (NAer=%d, NGas=%d, NSurf=%d)\n",
        t_l.model_build, t_l.rt_run, t_l.model_build + t_l.rt_run,
        t_l.NAer, t_l.NGas, t_l.NSurf)

# --- GPU variants ---
if CAN_USE_GPU
    println()
    print("Forward noRS (GPU)... ")
    t_g = benchmark_forward_noRS(vSmartMOM.Architectures.GPU(); nruns=nruns)
    @printf("model: %.3fs  |  rt_run: %.3fs  |  total: %.3fs\n",
            t_g.model_build, t_g.rt_run, t_g.model_build + t_g.rt_run)

    print("Raman RRS    (GPU)... ")
    t_rg = benchmark_raman(vSmartMOM.Architectures.GPU(); nruns=nruns)
    @printf("rt_run: %.3fs  (nSpec=%d)\n", t_rg.rt_run, t_rg.nSpec)

    print("Linearized   (GPU)... ")
    t_lg = benchmark_linearized(vSmartMOM.Architectures.GPU(); nruns=nruns)
    @printf("model: %.3fs  |  rt_run: %.3fs  |  total: %.3fs\n",
            t_lg.model_build, t_lg.rt_run, t_lg.model_build + t_lg.rt_run)

    # Speedup summary
    println("\n--- GPU Speedup ---")
    @printf("Forward noRS: %.1fx\n", t.rt_run / t_g.rt_run)
    @printf("Raman RRS:    %.1fx\n", t_r.rt_run / t_rg.rt_run)
    @printf("Linearized:   %.1fx\n", t_l.rt_run / t_lg.rt_run)
end

println("\n" * "="^70)
println("  Benchmark complete.")
println("="^70)
