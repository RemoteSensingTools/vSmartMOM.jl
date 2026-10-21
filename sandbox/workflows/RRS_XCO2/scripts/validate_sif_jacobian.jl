#!/usr/bin/env julia

include(joinpath(@__DIR__, "common.jl"))
using .RRSXCO2Common
using vSmartMOM
using vSmartMOM.CoreRT
using Test

params = load_parameters()
params.architecture = CPU()
# Preserve the band endpoints and reference neighbourhood while keeping this
# finite-difference regression inexpensive.
params.spec_bands = [band[unique(round.(Int, range(1, length(band); length=9)))]
                     for band in params.spec_bands]
model, lin_model = model_from_parameters(LinMode(), params; external_solar=true)
NAer = length(params.scattering_params.rt_aerosols)
NGas = size(lin_model.τ̇_abs[1], 1)
NSurf = length(params.brdf[1].legendre_coeff)

state = RRSXCO2Common.campaign_sif_state()
base_ref, base_slope = state.SIF760, state.mSIF
sources = sources_for_band(params, 1; SIF760=base_ref, mSIF=base_slope)
analytic = rt_run(model, lin_model, NAer, NGas, NSurf;
                  i_band=1, sources)
sif_cols = CoreRT.sif_range(analytic.layout)

function forward(sif760, msif)
    src = sources_for_band(params, 1; SIF760=sif760, mSIF=msif)
    return Array(rt_run_toa(model; i_band=1, sources=src))
end

@testset "SIF760/mSIF analytic Jacobian" begin
    for (k, h) in ((1, 1e-7), (2, 1e-9))
        plus = k == 1 ? forward(base_ref + h, base_slope) :
                        forward(base_ref, base_slope + h)
        minus = k == 1 ? forward(base_ref - h, base_slope) :
                         forward(base_ref, base_slope - h)
        fd = (plus .- minus) ./ (2h)
        ad = Array(@view analytic.toa_jacobian[:, :, :, sif_cols[k]])
        @test isapprox(ad, fd; rtol=2e-5, atol=2e-8)
    end
    @test CoreRT.surface_sif_parameter_count(sources_for_band(params, 2)) == 0
    @test CoreRT.surface_sif_parameter_count(sources_for_band(params, 3)) == 0
end

@testset "SIF amplitude linear superposition" begin
    legacy = sif_reference_state(total_sif=0.5, reference_wavelength_nm=760)
    scale = base_ref / legacy.SIF760
    no_sif = forward(0.0, 0.0)
    legacy_sif = forward(legacy.SIF760, legacy.mSIF)
    corrected_direct = forward(base_ref, base_slope)
    corrected_reconstructed =
        no_sif .+ scale .* (legacy_sif .- no_sif)
    @test isapprox(corrected_reconstructed, corrected_direct;
                   rtol=2e-13, atol=2e-14)
end
