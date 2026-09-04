#!/usr/bin/env julia

using Test
using vSmartMOM
using vSmartMOM.CoreRT

include(joinpath(@__DIR__, "validate_bottom_layer_retrievals.jl"))
include(joinpath(@__DIR__, "..", "scripts", "common.jl"))
using .RRSXCO2Common

@testset "760-nm angular-integral SIF normalization" begin
    state = RRSXCO2Common.campaign_sif_state()
    L760 = state.SIF760 * 1.0e7 /
        RRSXCO2Common.SIF_REFERENCE_WAVELENGTH_NM^2
    @test RRSXCO2Common.SIF_CASE_ON == "angular_integral760_0p5"
    @test isapprox(L760, 0.5 / (2π); atol=2e-16, rtol=0)
    @test isapprox(2π * L760, 0.5; atol=2e-15, rtol=0)
    @test isapprox(state.SIF760, 0.004596394756493938;
                   atol=2e-18, rtol=0)
    @test isapprox(state.mSIF, 1.2291230681458325e-5;
                   atol=2e-19, rtol=0)

    source = SurfaceSIF(
        SIF760=state.SIF760,
        mSIF=state.mSIF,
        wavenumber_cm1=[state.ν_ref],
    )
    prepared = CoreRT.prepare_source(source, Float64, 3, 1, Array)
    @test prepared.SIF₀[1, 1] == π * state.SIF760
    @test (0.5 / π) * (2 * prepared.SIF₀[1, 1]) == state.SIF760
    @test all(iszero, prepared.SIF₀[2:3, :])
end

function outcome_attributes(outcome, converged, fit_quality_ok)
    return Dict{String,Any}(
        "retrieval_complete" => 1,
        "outcome" => outcome,
        "converged" => converged,
        "fit_quality_ok" => fit_quality_ok,
    )
end

function test_prior(; names=["parameter_$index" for index in 1:30],
                    state=collect(Float64, 1:30),
                    covariance=Matrix{Float64}(I, 30, 30))
    names = copy(names)
    names[ACTIVE_LOG_AOD_INDICES] .= ACTIVE_LOG_AOD_NAMES
    covariance = copy(covariance)
    for index in ACTIVE_LOG_AOD_INDICES
        covariance[index, index] = REQUIRED_LOG_AOD_VARIANCE
    end
    return RetrievalPrior(:urban, copy(state), covariance,
                          collect(1:30), names)
end

@testset "exact embedded campaign prior validation" begin
    expected = test_prior()
    @test isnothing(validate_embedded_prior(
        copy(expected.parameter_names), copy(expected.xa), copy(expected.Sa),
        expected, "retrieval.nc"))

    changed_names = copy(expected.parameter_names)
    changed_names[14] = "wrong_parameter"
    @test_throws ErrorException validate_embedded_prior(
        changed_names, copy(expected.xa), copy(expected.Sa), expected,
        "retrieval.nc")

    mislabeled_prior = test_prior()
    mislabeled_prior.parameter_names[14] = "not_ln_aod"
    @test_throws ErrorException validate_embedded_prior(
        copy(mislabeled_prior.parameter_names), copy(mislabeled_prior.xa),
        copy(mislabeled_prior.Sa), mislabeled_prior, "retrieval.nc")

    changed_state = copy(expected.xa)
    changed_state[1] = nextfloat(changed_state[1])
    @test_throws ErrorException validate_embedded_prior(
        copy(expected.parameter_names), changed_state, copy(expected.Sa),
        expected, "retrieval.nc")

    changed_covariance = copy(expected.Sa)
    changed_covariance[1, 2] = nextfloat(changed_covariance[1, 2])
    @test_throws ErrorException validate_embedded_prior(
        copy(expected.parameter_names), copy(expected.xa), changed_covariance,
        expected, "retrieval.nc")

    stale_covariance = copy(expected.Sa)
    for index in ACTIVE_LOG_AOD_INDICES
        stale_covariance[index, index] = 2.0^2
    end
    stale_prior = RetrievalPrior(
        :urban, copy(expected.xa), stale_covariance, collect(1:30),
        copy(expected.parameter_names))
    @test_throws ErrorException validate_embedded_prior(
        copy(stale_prior.parameter_names), copy(stale_prior.xa),
        copy(stale_prior.Sa), stale_prior, "retrieval.nc")
end

@testset "aerosol SIF-stage smoke selection" begin
    aerosol_off = (
        surface=:urban,
        aerosol_case=:aod760_0p28,
        sif_case=:off,
        bottom_layer_index=16,
        bottom_co2_ppm=400.0,
        fixed_upper_co2_ppm=400.0,
    )
    aerosol_on = merge(
        aerosol_off, (sif_case=:angular_integral760_0p5,))
    truth = Dict(13 => aerosol_off, 18 => aerosol_on)
    @test isnothing(validate_aerosol_controls_selection([13], [11], truth))
    @test isnothing(validate_aerosol_controls_selection([18], [11], truth))
    @test isnothing(validate_aerosol_controls_selection([13, 18], [11], truth))
    @test_throws ErrorException validate_aerosol_controls_selection(
        [13, 18], [1], truth)
    @test_throws ErrorException validate_aerosol_controls_selection(
        [43], [11], truth)
    bad_truth = Dict(13 => merge(aerosol_off, (aerosol_case=:none,)),
                     18 => aerosol_on)
    @test_throws ErrorException validate_aerosol_controls_selection(
        [13], [11], bad_truth)
end

@testset "class-aware bottom-layer outcome validation" begin
    closure = outcome_attributes(1, 1, 1)
    @test validate_outcome_attributes(
        closure, "corrected.nc";
        require_scientific_closure=true).outcome == 1

    converged_bad_fit = outcome_attributes(2, 1, 0)
    @test validate_outcome_attributes(
        converged_bad_fit, "uncorrected.nc").outcome == 2
    @test_throws ErrorException validate_outcome_attributes(
        converged_bad_fit, "corrected.nc";
        require_scientific_closure=true)

    maximum_iterations = outcome_attributes(3, 0, 1)
    @test validate_outcome_attributes(
        maximum_iterations, "uncorrected.nc").outcome == 3
    @test_throws ErrorException validate_outcome_attributes(
        maximum_iterations, "corrected.nc";
        require_scientific_closure=true)

    @test_throws ErrorException validate_outcome_attributes(
        outcome_attributes(5, 0, 0), "illegal.nc")
    @test_throws ErrorException validate_outcome_attributes(
        outcome_attributes(1, 0, 1), "inconsistent.nc")
    @test_throws ErrorException validate_outcome_attributes(
        outcome_attributes(3, 1, 1), "inconsistent.nc")
end

@testset "desert control smoke selection" begin
    desert_control = (
        surface=:desert,
        aerosol_case=:none,
        sif_case=:off,
        bottom_layer_index=16,
        bottom_co2_ppm=400.0,
        fixed_upper_co2_ppm=400.0,
    )
    truth = Dict(43 => desert_control)
    @test isnothing(validate_desert_control_selection([43], [11], truth))
    @test_throws ErrorException validate_desert_control_selection(
        [3], [11], truth)
    @test_throws ErrorException validate_desert_control_selection(
        [43], [1], truth)
    bad_truth = Dict(43 => merge(desert_control, (surface=:urban,)))
    @test_throws ErrorException validate_desert_control_selection(
        [43], [11], bad_truth)
end

@testset "cross-host structural barrier scene selection" begin
    clear_nosif = (aerosol_case=:none, sif_case=:off)
    clear_sif = (
        aerosol_case=:none, sif_case=:angular_integral760_0p5)
    aerosol_nosif = (aerosol_case=:aod760_0p28, sif_case=:off)
    aerosol_sif = (
        aerosol_case=:aod760_0p28,
        sif_case=:angular_integral760_0p5,
    )

    for (state, truth) in enumerate(
            (clear_nosif, clear_sif, aerosol_nosif, aerosol_sif))
        @test isnothing(validate_scene_membership(truth, state, "all"))
    end
    @test isnothing(validate_scene_membership(clear_nosif, 1, "clear_nosif"))
    @test_throws ErrorException validate_scene_membership(
        clear_sif, 6, "clear_nosif")
    @test isnothing(validate_scene_membership(
        aerosol_sif, 16, "aerosol_all_sif"))
    @test_throws ErrorException validate_scene_membership(
        clear_sif, 6, "aerosol_all_sif")
end
