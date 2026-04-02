#=
Raman Optimization Regression Test
====================================
Runs the same 1-band computation as the baseline and compares results
against the saved reference for bit-exact match.

Reports: timing improvement, allocation reduction, and correctness.

Usage:
  julia --project=/home/sanghavi/code/github/vSmartMOM.jl \
        /home/sanghavi/code/github/vSmartMOM.jl/test/benchmarks/raman_optimization_compare.jl
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

# ── Load reference ──
ref_file = joinpath(OUTPUT_DIR, "raman_reference.jld2")
if !isfile(ref_file)
    error("Reference file not found: $ref_file\nRun raman_optimization_baseline.jl first!")
end
println("Loading reference from: $ref_file")
@load ref_file R1 T1 ieR1 ieT1 R0 T0 ν F₀
R1_ref, T1_ref, ieR1_ref, ieT1_ref = R1, T1, ieR1, ieT1
R0_ref, T0_ref = R0, T0

println("="^60)
println("Raman Optimization Regression Test")
println("="^60)
println("GPU: ", CUDA.device())
println("Free GPU memory: ", round(CUDA.available_memory() / 1e9, digits=2), " GB")
println()

# ── Load parameters (1-band config) ──
yaml_path = joinpath(PROJECT_ROOT, "test", "test_parameters", "O2_parameters2_1band_opttest.yaml")
parameters = parameters_from_yaml(yaml_path)
FT = parameters.float_type
model = model_from_parameters(parameters)

model.obs_geom = CoreRT.ObsGeometry(parameters.sza, parameters.vza, parameters.vaz, parameters.obs_alt)
model.quad_points = CoreRT.rt_set_streams(parameters.quadrature_type,
                                          parameters.l_trunc,
                                          model.obs_geom,
                                          parameters.polarization_type,
                                          array_type(parameters.architecture))

# ── Setup Raman types (identical to baseline) ──
ν_model = model.params.spec_bands[1]
ν̃ = mean(ν_model)
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

CoreRT.getRamanSSProp!(RS_type1, 1e7 / mean(ν_model), ν_model)

Tsolar = solar_transmission_from_file(joinpath(PROJECT_ROOT, "src", "SolarModel", "solar.out"))
Tsolar_interp = LinearInterpolation(Tsolar[4:end, 1], Tsolar[4:end, 2])
T_sun = 5777.0
P = planck_spectrum_wn(T_sun, ν_model) * 2.1629e-05 * π
F₀_new = Tsolar_interp.(ν_model) .* P

RS_type0.F₀ = zeros(model.params.polarization_type.n, length(P))
RS_type1.F₀ = zeros(model.params.polarization_type.n, length(P))
RS_type0.SIF₀ = zeros(model.params.polarization_type.n, length(ν_model))
RS_type1.SIF₀ = zeros(model.params.polarization_type.n, length(ν_model))
RS_type0.F₀[1, :] = F₀_new
RS_type1.F₀[1, :] = F₀_new

# ── Warmup ──
println("Warmup run...")
CoreRT.rt_run_test(RS_type1, model, 1)
CUDA.reclaim()

# ── Timed run ──
println("Benchmark run (RRS)...")
CUDA.reclaim()
t_rrs = @elapsed begin
    R1_new, T1_new, ieR1_new, ieT1_new = CoreRT.rt_run_test(RS_type1, model, 1)
end

println("Benchmark run (noRS)...")
t_nors = @elapsed begin
    R0_new, T0_new, _, _ = CoreRT.rt_run_test(RS_type0, model, 1)
end

# ── Correctness checks ──
println("\n" * "="^60)
println("CORRECTNESS CHECKS")
println("="^60)

pass = true

function check_match(name, ref, new)
    if size(ref) != size(new)
        println("  FAIL $name: size mismatch $(size(ref)) vs $(size(new))")
        return false
    end
    exact = all(ref .== new)
    if exact
        println("  PASS $name: bit-exact match")
        return true
    else
        maxdiff = maximum(abs.(ref .- new))
        reldiff = maximum(abs.(ref .- new) ./ max.(abs.(ref), 1e-30))
        nmismatch = sum(ref .!= new)
        println("  FAIL $name: $nmismatch / $(length(ref)) elements differ")
        println("       max abs diff: $maxdiff")
        println("       max rel diff: $reldiff")
        return false
    end
end

pass &= check_match("R1 (Raman reflectance)", R1_ref, R1_new)
pass &= check_match("T1 (Raman transmittance)", T1_ref, T1_new)
pass &= check_match("ieR1 (inelastic reflectance)", ieR1_ref, ieR1_new)
pass &= check_match("ieT1 (inelastic transmittance)", ieT1_ref, ieT1_new)
pass &= check_match("R0 (noRS reflectance)", R0_ref, R0_new)
pass &= check_match("T0 (noRS transmittance)", T0_ref, T0_new)

# ── Timing comparison ──
println("\n" * "="^60)
println("TIMING")
println("="^60)

# Read baseline timing
timing_file = joinpath(OUTPUT_DIR, "timing_baseline.txt")
t_rrs_baseline = NaN
if isfile(timing_file)
    for line in readlines(timing_file)
        m = match(r"RRS time \(s\):\s+([\d.]+)", line)
        if m !== nothing
            t_rrs_baseline = parse(Float64, m.captures[1])
        end
    end
end

println("  RRS time (now):      $(round(t_rrs, digits=3)) s")
println("  noRS time (now):     $(round(t_nors, digits=3)) s")
if !isnan(t_rrs_baseline)
    speedup = t_rrs_baseline / t_rrs
    println("  RRS baseline:        $(round(t_rrs_baseline, digits=3)) s")
    println("  Speedup:             $(round(speedup, digits=2))x")
end

# ── Final verdict ──
println("\n" * "="^60)
if pass
    println("ALL CHECKS PASSED")
else
    println("SOME CHECKS FAILED - results differ from reference!")
end
println("="^60)
