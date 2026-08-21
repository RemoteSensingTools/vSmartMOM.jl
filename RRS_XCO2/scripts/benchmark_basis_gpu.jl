#!/usr/bin/env julia

include(joinpath(@__DIR__, "common.jl"))
using .RRSXCO2Common
using vSmartMOM
using vSmartMOM.CoreRT
using CUDA
using JLD2
using Printf
using Statistics
using Dates
using DelimitedFiles

const GPU_DEVICE = parse(Int, get(ENV, "CUDA_DEVICE", "1"))
const NREPS = parse(Int, get(ENV, "NREPS", "2"))
const RUN_LINEAR = get(ENV, "RUN_LINEAR", "1") == "1"
const BAND_INDICES = parse.(Int, split(get(ENV, "BAND_INDICES", "1,2,3"), ','))
const OUTDIR = get(ENV, "RRS_XCO2_OUTPUT_DIR", joinpath(ROOT, "output"))

CUDA.functional() || error("CUDA is not functional on this host")
CUDA.device!(GPU_DEVICE)
mkpath(OUTDIR)

function timed_cuda(f)
    CUDA.synchronize()
    t = @elapsed value = f()
    CUDA.synchronize()
    return value, t
end

function labels(layout, band)
    result = ["surface_pressure"]
    aerosol_names = ("sulfate", "organic_carbon", "stratospheric_sulfate")
    aerosol_fields = ("tau_ref", "n_real", "n_imag", "median_radius",
                      "geometric_width", "profile_location", "profile_width")
    for aerosol in aerosol_names, field in aerosol_fields
        push!(result, "$(aerosol)_$(field)")
    end
    n_layers = layout.n_gases ÷ 2
    for species in ("h2o", "co2"), iz in 1:n_layers
        push!(result, "$(species)_layer_$(lpad(iz, 2, '0'))")
    end
    append!(result, ("surface_P0", "surface_P1", "surface_P2"))
    layout.n_sif == 2 && append!(result, ("SIF760", "mSIF"))
    length(result) == CoreRT.n_total(layout) ||
        error("Label count $(length(result)) does not match layout $(CoreRT.n_total(layout))")
    return result
end

println("="^78)
println("OCO three-band basis benchmark")
println("timestamp: ", now())
println("GPU device index: ", GPU_DEVICE, "  name: ", CUDA.name(CUDA.device()))
println("repetitions: ", NREPS)
println("output: ", OUTDIR)
println("="^78)

params = load_parameters()
@assert all(all(isapprox.(diff(band), 0.1; atol=1e-10, rtol=0))
            for band in params.spec_bands)
println("spectral points: ", length.(params.spec_bands))

model, forward_build_s = timed_cuda(() ->
    model_from_parameters(params; external_solar=true))
n_layers = size(model.τ_rayl[1], 2)
n_layers == 12 || error("Expected 12 layers, got $n_layers")

forward_rt_s = Dict{String,Vector{Float64}}()
for iband in BAND_INDICES
    band = BAND_NAMES[iband]
    # Compile/first-touch separately from reported warm timings.
    band_sources = sources_for_band(params, iband)
    warm, cold_s = timed_cuda(() -> rt_run_toa(model; i_band=iband, sources=band_sources))
    times = Float64[]
    result = warm
    for _ in 1:NREPS
        result, t = timed_cuda(() -> rt_run_toa(model; i_band=iband, sources=band_sources))
        push!(times, t)
    end
    forward_rt_s[band] = times
    ν = collect(params.spec_bands[iband])
    jldsave(joinpath(OUTDIR, "forward_$(band).jld2");
            band, wavenumber_cm1=ν, wavelength_nm=wavelengths_nm(ν),
            toa=Array(result), n_layers, cold_rt_s=cold_s, warm_rt_s=times,
            operator_size=length(model.quad_points.qp_μN))
    writedlm(joinpath(OUTDIR, "forward_$(band).csv"),
             hcat(wavelengths_nm(ν), vec(Array(result[1, 1, :]))), ',')
    @printf("forward %-10s cold=%8.3f s warm median=%8.3f s\n",
            band, cold_s, median(times))
end

# Measure linearized feasibility independently of retained forward workspaces.
model = nothing
GC.gc(true)
CUDA.reclaim()

linear_status = "not requested"
linear_build_s = NaN
linear_rt_s = Dict{String,Vector{Float64}}()
if RUN_LINEAR
    try
        lin_model_pair, build_time = timed_cuda(() ->
            model_from_parameters(LinMode(), params; external_solar=true))
        global linear_build_s = build_time
        lin_forward_model, lin_model = lin_model_pair
        NAer = length(params.scattering_params.rt_aerosols)
        for iband in BAND_INDICES
            band = BAND_NAMES[iband]
            NGas = size(lin_model.τ̇_abs[iband], 1)
            NSurf = 3
            band_sources = sources_for_band(params, iband)
            warm, cold_s = timed_cuda(() ->
                rt_run(lin_forward_model, lin_model, NAer, NGas, NSurf;
                       i_band=iband, sources=band_sources))
            times = Float64[]
            result = warm
            for _ in 1:NREPS
                result, t = timed_cuda(() ->
                    rt_run(lin_forward_model, lin_model, NAer, NGas, NSurf;
                           i_band=iband, sources=band_sources))
                push!(times, t)
            end
            linear_rt_s[band] = times
            ν = collect(params.spec_bands[iband])
            parameter_labels = labels(result.layout, band)
            jldsave(joinpath(OUTDIR, "linearized_$(band).jld2");
                    band, wavenumber_cm1=ν, wavelength_nm=wavelengths_nm(ν),
                    toa=Array(result.toa), jacobian=Array(result.toa_jacobian),
                    parameter_labels, n_layers, cold_rt_s=cold_s, warm_rt_s=times)
            jac_I = Array(result.toa_jacobian[1, 1, :, :])
            csv_path = joinpath(OUTDIR, "linearized_$(band).csv")
            open(csv_path, "w") do io
                println(io, join(vcat("wavelength_nm", parameter_labels), ','))
                writedlm(io, hcat(wavelengths_nm(ν), jac_I), ',')
            end
            @printf("linear  %-10s cold=%8.3f s warm median=%8.3f s\n",
                    band, cold_s, median(times))
        end
        global linear_status = "completed"
    catch err
        global linear_status = sprint(showerror, err, catch_backtrace())
        println(stderr, "Linearized benchmark failed:\n", linear_status)
    end
end

summary_path = joinpath(OUTDIR, "basis_benchmark_summary.jld2")
jldsave(summary_path;
        timestamp=string(now()), gpu_device=GPU_DEVICE,
        gpu_name=string(CUDA.name(CUDA.device())), n_layers,
        spectral_points=length.(params.spec_bands),
        forward_build_s, forward_rt_s,
        linear_build_s, linear_rt_s, linear_status)
println("summary: ", summary_path)
