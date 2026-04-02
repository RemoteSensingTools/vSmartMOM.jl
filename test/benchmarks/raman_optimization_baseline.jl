#=
Raman Optimization Baseline Benchmark
======================================
Runs a single-band Raman RT computation and saves:
  1. Reference output arrays (R, T, ieR, ieT) to JLD2 for bit-exact comparison
  2. Timing and allocation statistics for performance comparison

Output directory: test/benchmarks/raman_opttest_output/
This script NEVER writes to /home/sanghavi/data/ to avoid overwriting production results.

Usage:
  julia --project=/home/sanghavi/code/github/vSmartMOM.jl \
        /home/sanghavi/code/github/vSmartMOM.jl/test/benchmarks/raman_optimization_baseline.jl
=#

using CUDA
device!(1)
using vSmartMOM, vSmartMOM.CoreRT, vSmartMOM.SolarModel
using vSmartMOM.InelasticScattering
using Statistics
using Interpolations
using DelimitedFiles
using JLD2

const PROJECT_ROOT = "/home/sanghavi/code/github/vSmartMOM.jl"
const OUTPUT_DIR   = joinpath(PROJECT_ROOT, "test", "benchmarks", "raman_opttest_output")
mkpath(OUTPUT_DIR)

println("="^60)
println("Raman Optimization Baseline Benchmark")
println("="^60)
println("Output directory: $OUTPUT_DIR")
println("GPU: ", CUDA.device())
println("Free GPU memory: ", round(CUDA.available_memory() / 1e9, digits=2), " GB")
println()

# ── Load parameters (1-band config) ──
yaml_path = joinpath(PROJECT_ROOT, "test", "test_parameters", "O2_parameters2_1band_opttest.yaml")
println("Loading parameters from: $yaml_path")
parameters = parameters_from_yaml(yaml_path)
FT = parameters.float_type

# ── Create model ──
println("Creating model...")
t_model = @elapsed model = model_from_parameters(parameters)
println("  Model creation: $(round(t_model, digits=2)) s")

# ── Setup geometry (fixed: sza=19, vza=0, albedo=0) ──
model.obs_geom = CoreRT.ObsGeometry(parameters.sza, parameters.vza, parameters.vaz, parameters.obs_alt)
model.quad_points = CoreRT.rt_set_streams(parameters.quadrature_type,
                                          parameters.l_trunc,
                                          model.obs_geom,
                                          parameters.polarization_type,
                                          array_type(parameters.architecture))

# ── Setup Raman types ──
ν = model.params.spec_bands[1]
ν̃ = mean(ν)
effT = (model.profile.vcd_dry' * model.profile.T) / sum(model.profile.vcd_dry)
n2, o2 = InelasticScattering.getRamanAtmoConstants(ν̃, effT)

RS_type0 = InelasticScattering.noRS(
    fscattRayl  = [FT(1)],
    ϖ_Cabannes  = [FT(1)],
    bandSpecLim = [],
    iBand       = [1],
    F₀          = zeros(FT, 1, 1),
    SIF₀        = zeros(FT, 1, 1))

RS_type1 = InelasticScattering.RRS(
    n2=n2, o2=o2,
    greek_raman = InelasticScattering.GreekCoefs([FT(1)], [FT(1)], [FT(1)], [FT(1)], [FT(1)], [FT(1)]),
    fscattRayl  = [FT(1)],
    ϖ_Cabannes  = [FT(1)],
    ϖ_λ₁λ₀      = zeros(FT, 1),
    i_λ₁λ₀      = zeros(Int, 1),
    Z⁻⁺_λ₁λ₀    = zeros(FT, 1, 1),
    Z⁺⁺_λ₁λ₀    = zeros(FT, 1, 1),
    i_ref       = 0,
    n_Raman     = 0,
    F₀          = zeros(FT, 1, 1),
    SIF₀        = zeros(FT, 1, 1))

# ── Compute Raman properties ──
println("Computing Raman SS properties...")
CoreRT.getRamanSSProp!(RS_type1, 1e7 / mean(ν), ν)
println("  n_Raman = $(RS_type1.n_Raman)")
println("  nSpec   = $(length(ν))")

# ── Setup solar spectrum ──
Tsolar = solar_transmission_from_file(joinpath(PROJECT_ROOT, "src", "SolarModel", "solar.out"))
Tsolar_interp = LinearInterpolation(Tsolar[4:end, 1], Tsolar[4:end, 2])
T_sun = 5777.0
P = planck_spectrum_wn(T_sun, ν) * 2.1629e-05 * π
F₀ = Tsolar_interp.(ν) .* P

RS_type0.F₀ = zeros(model.params.polarization_type.n, length(P))
RS_type1.F₀ = zeros(model.params.polarization_type.n, length(P))
RS_type0.SIF₀ = zeros(model.params.polarization_type.n, length(ν))
RS_type1.SIF₀ = zeros(model.params.polarization_type.n, length(ν))
RS_type0.F₀[1, :] = F₀
RS_type1.F₀[1, :] = F₀

# ── Warmup run (JIT compilation) ──
println("\nWarmup run (JIT compilation)...")
t_warmup = @elapsed begin
    R1_w, T1_w, ieR1_w, ieT1_w = CoreRT.rt_run_test(RS_type1, model, 1)
end
println("  Warmup: $(round(t_warmup, digits=2)) s")
CUDA.reclaim()

# ── Timed benchmark run (Raman) ──
println("\nBenchmark run (RRS, 1 band)...")
CUDA.reclaim()
gpu_mem_before = CUDA.available_memory()

t_rrs = @elapsed begin
    R1, T1, ieR1, ieT1 = CoreRT.rt_run_test(RS_type1, model, 1)
end
gpu_mem_after = CUDA.available_memory()

println("  RRS time:      $(round(t_rrs, digits=3)) s")
println("  GPU mem used:  ~$(round((gpu_mem_before - gpu_mem_after) / 1e9, digits=2)) GB")

# ── Timed benchmark run (noRS) ──
println("\nBenchmark run (noRS, 1 band)...")
t_nors = @elapsed begin
    R0, T0, _, _ = CoreRT.rt_run_test(RS_type0, model, 1)
end
println("  noRS time:     $(round(t_nors, digits=3)) s")
println("  Speedup ratio: $(round(t_rrs / t_nors, digits=1))x slower with Raman")

# ── Save reference results ──
ref_file = joinpath(OUTPUT_DIR, "raman_reference.jld2")
println("\nSaving reference results to: $ref_file")
@save ref_file R1 T1 ieR1 ieT1 R0 T0 ν F₀

# ── Save timing info ──
timing_file = joinpath(OUTPUT_DIR, "timing_baseline.txt")
open(timing_file, "w") do f
    println(f, "Raman Optimization Baseline Timing")
    println(f, "Date: ", Dates.now())
    println(f, "Git commit: ", strip(read(`git -C $PROJECT_ROOT rev-parse --short HEAD`, String)))
    println(f, "")
    println(f, "Config: O2_parameters2_1band_opttest.yaml")
    println(f, "nSpec = $(length(ν))")
    println(f, "n_Raman = $(RS_type1.n_Raman)")
    println(f, "NquadN = $(size(R1, 1) > 0 ? "see output" : "N/A")")
    println(f, "")
    println(f, "RRS time (s):  $(round(t_rrs, digits=3))")
    println(f, "noRS time (s): $(round(t_nors, digits=3))")
    println(f, "Ratio:         $(round(t_rrs / t_nors, digits=1))x")
    println(f, "")
    println(f, "R1 size:  $(size(R1))")
    println(f, "ieR1 size: $(size(ieR1))")
end

println("\nBaseline complete. Files saved:")
println("  Reference: $ref_file")
println("  Timing:    $timing_file")
println("="^60)
