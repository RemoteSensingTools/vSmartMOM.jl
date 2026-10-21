#!/usr/bin/env julia

using DataInterpolations
using LinearAlgebra
using Test
using vSmartMOM.SolarModel: planck_spectrum_wn, solar_transmission_from_file

include(joinpath(@__DIR__, "OptimalEstimation.jl"))
include(joinpath(@__DIR__, "RetrievalState.jl"))
include(joinpath(@__DIR__, "VSmartMOMForward.jl"))
using .RetrievalState
using .VSmartMOMForward

@testset "Truth/retrieval O2 retained grid identity" begin
    for FT in (Float32, Float64)
        core = VSmartMOMForward.RRSXCO2Common.basis_grids(FT)[1]
        solve, keep =
            VSmartMOMForward.RRSXCO2Common.raman_solve_grid(core)
        nshoulder = round(Int, FT(234) / FT(0.1))
        @test solve[keep] == core
        @test length(solve) == length(core) + 2nshoulder
        @test first(keep) == nshoulder + 1
    end
end

@testset "scene-specific fixed upper CO2" begin
    evaluator = OCOForwardEvaluator(;
        architecture=:CPU, float_type=Float32, nstreams=9)
    prior = load_retrieval_prior(:urban)

    for truth_ppm in (380.0, 400.0, 420.0, 440.0)
        set_fixed_upper_co2_ppm!(evaluator, truth_ppm)
        params = deepcopy(evaluator.base_parameters)
        apply_retrieval_state!(
            params, prior.xa, evaluator.tau_ref_scale;
            fixed_upper_co2_vmr=evaluator.fixed_upper_co2_vmr)
        co2 = params.absorption_params.vmr["CO2"]

        @test length(co2) == 16
        @test co2[1:4] == fill(Float32(truth_ppm * 1e-6), 4)
        @test co2[5:16] == Float32.(prior.xa[2:13])

        uniform_state = copy(prior.xa)
        uniform_state[2:13] .= truth_ppm * 1e-6
        @test column_averaged_co2_ppm(evaluator, uniform_state) ≈ truth_ppm rtol=2e-7
    end


    # Independent hydrostatic weighting check for a vertically varying state.
    set_fixed_upper_co2_ppm!(evaluator, 380.0)
    varying_state = copy(prior.xa)
    varying_state[1] = 975.0
    varying_state[2:13] .= range(385e-6, 440e-6, length=12)
    p_half = Float64.(evaluator.base_parameters.p)
    p_half[end] = varying_state[1]
    q = Float64.(evaluator.base_parameters.q)
    dry_mass = 28.9644e-3
    wet_mass = 18.01534e-3
    h2o_dry_ratio = q ./ (1 .- q) .* (dry_mass / wet_mass)
    x_dry = 1.0 ./ (1.0 .+ h2o_dry_ratio)
    mixture_mass = x_dry .* dry_mass .+ (1.0 .- x_dry) .* wet_mass
    dry_weights = diff(p_half) .* x_dry ./ mixture_mass
    co2 = vcat(fill(380e-6, 4), varying_state[2:13])
    expected_xco2 = 1e6 * dot(co2, dry_weights) / sum(dry_weights)
    @test column_averaged_co2_ppm(evaluator, varying_state) ≈
          expected_xco2 rtol=3e-7

    @test_throws ArgumentError set_fixed_upper_co2_ppm!(evaluator, 0.0)
    @test_throws ArgumentError set_fixed_upper_co2_ppm!(evaluator, Inf)
end

@testset "truth-consistent solar source" begin
    evaluator = OCOForwardEvaluator(;
        architecture=:CPU, float_type=Float32, nstreams=9)
    params = evaluator.base_parameters
    raw = solar_transmission_from_file(VSmartMOMForward.RRSXCO2Common.SOLAR_OUT)
    truth_solar = LinearInterpolation(
        Float32.(raw[4:end, 2]), Float32.(raw[4:end, 1]);
        extrapolation=ExtrapolationType.Linear)

    for iband in 1:3
        sources = VSmartMOMForward.RRSXCO2Common.sources_for_band(
            params, iband; SIF760=0, mSIF=0,
            solar_T=evaluator.solar_transmission)
        beam = iband == 1 ? sources[1] : sources
        ν = Float32.(params.spec_bands[iband])
        planck = Float32.(planck_spectrum_wn(Float32(5777), collect(ν)) .*
                          Float32(2.1629e-5 * π))
        expected = Float32.(truth_solar.(ν)) .* planck
        @test beam.F₀[1, :] == expected
        @test all(iszero, beam.F₀[2:3, :])
    end
end

@testset "truth-coordinate surface coefficients" begin
    evaluator = OCOForwardEvaluator(;
        architecture=:CPU, float_type=Float32, nstreams=9)
    params = evaluator.base_parameters
    coefficients = Float32[0.24, -0.035, 0.004]

    for iband in 1:3
        transform = surface_coefficient_transform(params, iband)
        native = transform * coefficients
        solve = params.spec_bands[iband]
        x_native = collect(range(Float32(-1), Float32(1), length=length(solve)))
        actual = native[1] .+ native[2] .* x_native .+
                 native[3] .* (3 .* x_native.^2 .- 1) ./ 2

        base = VSmartMOMForward.RRSXCO2Common.surface_basis_grids(Float32)[iband]
        x_truth = 2 .* (solve .- first(base)) ./ (last(base) - first(base)) .- 1
        expected = coefficients[1] .+ coefficients[2] .* x_truth .+
                   coefficients[3] .* (3 .* x_truth.^2 .- 1) ./ 2
        @test actual ≈ expected rtol=3e-6 atol=3e-7
    end

    @test surface_coefficient_transform(params, 1) ≈ Matrix{Float32}(I, 3, 3)
    @test surface_coefficient_transform(params, 2) ≈ Matrix{Float32}(I, 3, 3)
    @test surface_coefficient_transform(params, 3) != Matrix{Float32}(I, 3, 3)
end
