#!/usr/bin/env julia

using LinearAlgebra
using Random
using Statistics
using Test
using NCDatasets

include(joinpath(@__DIR__, "OptimalEstimation.jl"))
include(joinpath(@__DIR__, "RetrievalCases.jl"))
include(joinpath(@__DIR__, "RetrievalOutput.jl"))
using .OptimalEstimation
using .RetrievalCases
using .RetrievalOutput

@testset "Connor-style optimal estimation" begin
    K = [1.0 0.2;
         0.7 1.1;
         0.1 1.5;
         1.3 -0.2]
    offset = [0.2, -0.1, 0.4, 0.7]
    truth = [1.2, -0.4]
    measurement = offset + K * truth
    variance = fill(0.05^2, length(measurement))
    xa = [0.0, 0.0]
    Sa = Matrix(Diagonal([1.5^2, 1.0^2]))
    ranges = [1:2, 3:4]

    evaluate(x) = ForwardEvaluation(offset + K * x, K, ranges)
    result = solve_optimal_estimation(
        evaluate, measurement, variance, xa, Sa;
        settings=OESettings(convergence_threshold=1e-8,
                            maximum_iterations=7))

    W = Diagonal(1.0 ./ variance)
    expected = (K' * W * K + inv(Sa)) \
               (K' * W * (measurement - offset) + inv(Sa) * xa)
    @test result.converged
    @test result.outcome == 1
    @test result.final_state ≈ expected rtol=2e-8 atol=2e-10
    @test result.final_jacobian == K
    @test result.posterior_covariance ≈ inv(K' * W * K + inv(Sa)) rtol=1e-12
    @test result.gain_matrix ≈ result.posterior_covariance * K' * W rtol=1e-12
    @test result.averaging_kernel ≈ result.gain_matrix * K rtol=1e-12
    @test length(result.final_band_chi_squared) == 2
    @test all(record -> length(record.state) == 2, result.records)
end

@testset "class-specific NetCDF output schema" begin
    K = [1.0 0.0; 0.5 1.0; 0.2 0.7; 1.1 -0.1]
    xa = zeros(2)
    Sa = Matrix{Float64}(I, 2, 2)
    variance = fill(0.1^2, 4)
    truth_state = [0.3, -0.2]
    measurement = K * truth_state
    ranges = [1:2, 3:3, 4:4]
    evaluate(x) = ForwardEvaluation(K * x, K, ranges)
    result = solve_optimal_estimation(evaluate, measurement, variance, xa, Sa)
    # Output-schema validation must not depend on the production campaign
    # having completed all 32 no-SIF measurements. Construct one synthetic
    # experiment record directly; build_experiments deliberately enforces
    # completeness and is tested separately against production inputs.
    truth = TruthCase(1, 1, :urban, 1, :none, 1, 380.0)
    experiment = RetrievalExperiment(
        1, 1, truth, 1, :corrected, UInt64(1),
        "synthetic_measurement.nc", "synthetic_noise.nc")
    realization = MeasurementRealization(
        measurement, measurement, sqrt.(variance), variance, zeros(4),
        [760.0, 760.1, 1600.0, 2050.0], ranges)
    xco2_diagnostics = (
        a_priori_ppm=400.0,
        trial_ppm=collect(range(400.0, 405.0, length=length(result.records))),
        final_ppm=404.75,
    )

    mktempdir() do directory
        path = joinpath(directory, "corrected",
                        "retrieval_state001_perturbation01.nc")
        write_retrieval_result(
            experiment, realization, result, xa, Sa, ["x1", "x2"];
            output_path=path, xco2_diagnostics,
            provenance=Dict("spectroscopy_database" => "ABSCO"))
        @test isfile(path)
        NCDataset(path) do dataset
            @test dataset.attrib["retrieval_complete"] == 1
            @test dataset.attrib["truth_state_index"] == 1
            @test dataset.attrib["perturbation_index"] == 1
            @test dataset.attrib["spectroscopy_database"] == "ABSCO"
            @test size(dataset["final_jacobian"]) == (4, 2)
            @test size(dataset["gain_matrix"]) == (2, 4)
            @test length(dataset["trial_index"]) == length(result.records)
            @test dataset["a_priori_XCO2"][] == 400.0
            @test dataset["XCO2_at_trial"][:] == xco2_diagnostics.trial_ppm
            @test dataset["XCO2"][] == 404.75
            @test dataset["XCO2"].attrib["units"] == "ppm"
        end
    end
end

@testset "uniform perturbation has covariance scale one" begin
    rng = MersenneTwister(0x4f434f525253)
    draws = (2sqrt(3)) .* rand(rng, 1_000_000) .- sqrt(3)
    @test abs(mean(draws)) < 4e-3
    @test abs(var(draws; corrected=true) - 1) < 4e-3
    @test extrema(draws)[1] >= -sqrt(3)
    @test extrema(draws)[2] <= sqrt(3)
end
