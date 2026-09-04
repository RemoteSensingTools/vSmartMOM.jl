#!/usr/bin/env julia

using LinearAlgebra
using Random
using Statistics
using Test
using NCDatasets

include(joinpath(@__DIR__, "OptimalEstimation.jl"))
include(joinpath(@__DIR__, "RetrievalCases.jl"))
include(joinpath(@__DIR__, "RetrievalOutput.jl"))
include(joinpath(@__DIR__, "RetrievalState.jl"))
using .OptimalEstimation
using .RetrievalCases
using .RetrievalOutput
using .RetrievalState

@testset "truth-case SIF selection" begin
    all_cases = read_truth_cases()
    sif_off = select_sif_truth_cases(all_cases, :off)
    sif_on = select_sif_truth_cases(all_cases, "on")

    @test length(all_cases) == 64
    @test length(sif_off) == 32
    @test length(sif_on) == 32
    @test select_sif_truth_cases(all_cases, :all) == all_cases
    @test all(case -> case.sif_case == :off, sif_off)
    @test all(case -> case.sif_case != :off, sif_on)
    @test isempty(intersect(
        Set(case.state_index for case in sif_off),
        Set(case.state_index for case in sif_on)))
    @test union(
        Set(case.state_index for case in sif_off),
        Set(case.state_index for case in sif_on)) == Set(1:64)
    @test [case.state_index for case in read_no_sif_truth_cases()] ==
          [case.state_index for case in sif_off]
    @test_throws ArgumentError select_sif_truth_cases(all_cases, :enabled)

    off_experiments = build_experiments(
        [first(sif_off)]; validate_inputs=false)
    on_experiments = build_experiments(
        [first(sif_on)]; validate_inputs=false)
    off_perturbation_one = only(filter(experiment ->
        experiment.noise_index == 1 &&
        experiment.measurement_class == :corrected, off_experiments))
    on_perturbation_one = only(filter(experiment ->
        experiment.noise_index == 1 &&
        experiment.measurement_class == :corrected, on_experiments))
    @test off_perturbation_one.random_seed ==
          RetrievalCases.DEFAULT_MASTER_SEED + UInt64(10_001)
    @test on_perturbation_one.random_seed ==
          RetrievalCases.DEFAULT_MASTER_SEED + UInt64(50_001)
    @test basename(retrieval_output_path(off_perturbation_one)) ==
          "retrieval_state001_perturbation01.nc"
    @test basename(retrieval_output_path(on_perturbation_one)) ==
          "retrieval_state005_perturbation01.nc"
end

@testset "mapped ACOS CO2 prior covariance" begin
    prior = load_retrieval_prior(:urban)
    @test length(prior.xa) == 30
    @test isposdef(Symmetric(prior.Sa))
    # Active order is psurf, then CO2 layers 5:16.
    co2 = prior.Sa[2:13, 2:13]
    @test sqrt(co2[1, 1]) * 1e6 ≈ 6.839421605168 rtol=1e-11
    @test sqrt(co2[end, end]) * 1e6 ≈ 43.52698551448 rtol=1e-11
    @test maximum(abs.(co2 - Diagonal(diag(co2)))) * 1e12 > 1400
end

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
    experiment = only(filter(
        experiment -> experiment.noise_index == 1 &&
                      experiment.measurement_class == :corrected,
        build_experiments([first(read_no_sif_truth_cases())])))
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
            @test dataset.attrib["sif_case"] == "off"
            @test dataset.attrib["spectroscopy_database"] == "ABSCO"
            @test size(dataset["final_jacobian"]) == (4, 2)
            @test size(dataset["gain_matrix"]) == (2, 4)
            @test length(dataset["trial_index"]) == length(result.records)
            @test dataset["a_priori_XCO2"][] == 400.0
            @test dataset["XCO2_at_trial"][:] == xco2_diagnostics.trial_ppm
            @test dataset["XCO2"][] == 404.75
            @test dataset["XCO2"].attrib["units"] == "ppm"
            @test dataset["injected_measurement_noise"][:] == zeros(4)
            @test dataset.attrib["noise_injected"] == 1
        end
    end
end

@testset "perturbation 11 is unperturbed and scheduled first" begin
    truth = first(read_no_sif_truth_cases())
    experiments = build_experiments([truth])
    @test length(experiments) == 22
    @test [experiment.noise_index for experiment in experiments] ==
          repeat(vcat(11, collect(1:10)), inner=2)
    @test [experiment.measurement_class for experiment in experiments[1:4]] ==
          [:corrected, :uncorrected, :corrected, :uncorrected]
    @test [experiment.pair_index for experiment in experiments[1:4]] ==
          [11, 11, 1, 1]
    @test [experiment.retrieval_index for experiment in experiments[1:4]] ==
          [21, 22, 1, 2]
    perturbation_one = filter(experiment -> experiment.noise_index == 1,
                              experiments)
    @test length(perturbation_one) == 2
    @test all(experiment -> experiment.random_seed ==
              RetrievalCases.DEFAULT_MASTER_SEED + UInt64(10_001),
              perturbation_one)
    unperturbed = filter(
        experiment -> experiment.noise_index == UNPERTURBED_INDEX,
        experiments)
    @test length(unperturbed) == 2
    @test all(experiment -> experiment.random_seed == 0, unperturbed)
    @test Set(experiment.measurement_class for experiment in unperturbed) ==
          Set((:corrected, :uncorrected))
    for experiment in unperturbed
        realization = load_measurement_realization(experiment)
        @test all(iszero, realization.normalized_draw)
        @test realization.perturbed == realization.noiseless
    end
end

@testset "noise draws are paired only within one truth state" begin
    cases = read_no_sif_truth_cases()
    state001 = only(filter(case -> case.state_index == 1, cases))
    state009 = only(filter(case -> case.state_index == 9, cases))
    experiments = build_experiments(
        [state001, state009]; include_unperturbed=false,
        validate_inputs=false)

    function perturbation_pair(state_index, noise_index)
        selected = filter(experiment ->
            experiment.truth.state_index == state_index &&
            experiment.noise_index == noise_index, experiments)
        @test length(selected) == 2
        return selected
    end

    pair001 = perturbation_pair(1, 1)
    pair009 = perturbation_pair(9, 1)
    @test pair001[1].random_seed == pair001[2].random_seed
    @test pair009[1].random_seed == pair009[2].random_seed
    @test pair001[1].random_seed != pair009[1].random_seed

    draw001 = paired_uniform_draw(pair001[1].random_seed, 64)
    draw009 = paired_uniform_draw(pair009[1].random_seed, 64)
    @test draw001 == paired_uniform_draw(pair001[2].random_seed, 64)
    @test draw009 == paired_uniform_draw(pair009[2].random_seed, 64)
    @test draw001 != draw009
end

@testset "uniform perturbation has covariance scale one" begin
    rng = MersenneTwister(0x4f434f525253)
    draws = (2sqrt(3)) .* rand(rng, 1_000_000) .- sqrt(3)
    @test abs(mean(draws)) < 4e-3
    @test abs(var(draws; corrected=true) - 1) < 4e-3
    @test extrema(draws)[1] >= -sqrt(3)
    @test extrema(draws)[2] <= sqrt(3)
end
