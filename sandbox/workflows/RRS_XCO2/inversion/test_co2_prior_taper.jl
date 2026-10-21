#!/usr/bin/env julia

using LinearAlgebra
using NCDatasets
using Test

include(joinpath(@__DIR__, "retrieval_setup", "build_apriori.jl"))

const EXPECTED_TAPERED_ADJACENT_CORRELATION = [
    0.9599841058533383,
    0.9199699332712215,
    0.8724550511563143,
    0.8248973950782327,
    0.7881857693079238,
    0.7454300983633975,
    0.7008696412340173,
    0.6547195361487076,
    0.6024021854175606,
    0.5617510691168982,
    0.5013049539297996,
]

@testset "opt-in tapered mapped-ACOS CO2 covariance" begin
    default_prior = build_prior(:urban, 0.1)
    explicit_default = build_prior(
        :urban, 0.1;
        co2_covariance_model=ACOS_MAPPED_CO2_COVARIANCE_MODEL)
    tapered_prior = build_prior(
        :urban, 0.1;
        co2_covariance_model=TAPERED_CO2_COVARIANCE_MODEL)

    # The active campaign remains bit-for-bit on the original default unless
    # the new model is selected explicitly.
    @test DEFAULT_CO2_COVARIANCE_MODEL == "acos_mapped"
    @test CO2_COVARIANCE_MODEL == DEFAULT_CO2_COVARIANCE_MODEL
    @test default_prior.co2_covariance_model == "acos_mapped"
    @test default_prior.xa == explicit_default.xa
    @test default_prior.Sa == explicit_default.Sa

    @test tapered_prior.co2_covariance_model ==
        "acos_mapped_tapered_vertical_correlation"
    @test tapered_prior.co2_covariance_base_model == "acos_mapped"
    @test tapered_prior.co2_correlation_taper_type ==
        "nonstationary_ar1_schur"
    @test tapered_prior.xa == default_prior.xa
    @test tapered_prior.active == default_prior.active

    expected_retention = collect(range(0.98, 0.55; length=11))
    @test tapered_prior.co2_correlation_taper_adjacent_retention ≈
        expected_retention atol=2e-16 rtol=0
    @test tapered_prior.co2_selected_adjacent_correlation ≈
        EXPECTED_TAPERED_ADJACENT_CORRELATION atol=2e-14 rtol=0
    approved_rounded_correlations =
        [0.960, 0.920, 0.872, 0.825, 0.788, 0.745,
         0.701, 0.655, 0.602, 0.562, 0.501]
    @test maximum(abs.(
        tapered_prior.co2_selected_adjacent_correlation .-
        approved_rounded_correlations)) < 5e-4
    @test tapered_prior.co2_selected_adjacent_correlation ≈
        tapered_prior.co2_base_adjacent_correlation .* expected_retention
    @test all(diff(tapered_prior.co2_selected_adjacent_correlation) .< 0)
    @test all(tapered_prior.co2_selected_adjacent_correlation .> 0)

    base_co2 = default_prior.Sa[2:17, 2:17]
    tapered_co2 = tapered_prior.Sa[2:17, 2:17]
    @test all(iszero, tapered_co2[1:4, :])
    @test all(iszero, tapered_co2[:, 1:4])
    @test diag(tapered_co2) == diag(base_co2)
    @test isposdef(Symmetric(tapered_co2[5:16, 5:16]))
    @test isposdef(Symmetric(tapered_prior.Sa[tapered_prior.active,
                                               tapered_prior.active]))

    # Only the CO2 covariance block changes; every state mean and every other
    # covariance block remains identical to the current prior.
    covariance_difference = tapered_prior.Sa - default_prior.Sa
    covariance_difference[2:17, 2:17] .= 0
    @test all(iszero, covariance_difference)
    @test isapprox(co2_xco2_sigma_ppm(base_co2), 13.715633108;
                   atol=1e-9, rtol=0)
    @test isapprox(co2_xco2_sigma_ppm(tapered_co2), 9.309754891;
                   atol=1e-9, rtol=0)

    @test_throws ErrorException build_prior(
        :urban, 0.1; co2_covariance_model="misspelled_model")
    @test_throws ArgumentError nonstationary_ar1_correlation([0.5, 1.0])
    @test_throws ArgumentError nonstationary_ar1_correlation([0.5, NaN])

    mktempdir() do directory
        netcdf_path = joinpath(directory, "tapered_apriori.nc")
        summary_path = joinpath(directory, "tapered_apriori.dat")
        priors = Dict(
            surface => build_prior(
                surface, 0.1;
                co2_covariance_model=TAPERED_CO2_COVARIANCE_MODEL)
            for surface in SURFACES)
        write_netcdf(priors; output_path=netcdf_path)
        write_summary(priors; output_path=summary_path)

        NCDataset(netcdf_path) do dataset
            @test dataset.attrib["co2_covariance_model"] ==
                TAPERED_CO2_COVARIANCE_MODEL
            @test dataset.attrib["co2_covariance_base_model"] ==
                ACOS_MAPPED_CO2_COVARIANCE_MODEL
            @test dataset.attrib["co2_correlation_taper_type"] ==
                "nonstationary_ar1_schur"
            @test occursin("R_selected=R.*T",
                           dataset.attrib["co2_covariance_construction"])
            @test Float64.(dataset[
                "co2_correlation_taper_adjacent_retention"][:]) ≈
                expected_retention atol=2e-16 rtol=0
            @test Float64.(dataset[
                "co2_selected_adjacent_correlation"][:]) ≈
                EXPECTED_TAPERED_ADJACENT_CORRELATION atol=2e-14 rtol=0
            @test isapprox(
                Float64(dataset.attrib[
                    "co2_prior_xco2_sigma_ppm_at_1000hpa"]),
                9.309754891; atol=1e-9, rtol=0)
            @test isapprox(
                Float64(dataset.attrib[
                    "co2_covariance_base_xco2_sigma_ppm_at_1000hpa"]),
                13.715633108; atol=1e-9, rtol=0)
            @test isapprox(
                Float64(dataset.attrib[
                    "co2_prior_xco2_variance_ppm2_at_1000hpa"]),
                9.309754891308744^2; atol=1e-8, rtol=0)
        end

        summary = read(summary_path, String)
        @test occursin(
            "CO2 covariance model: acos_mapped_tapered_vertical_correlation",
            summary)
        @test occursin("CO2-only sigma(XCO2) at 1000 hPa = 9.309754891 ppm",
                       summary)
        @test occursin(
            "CO2-only sigma(XCO2) consequence: base=13.715633108 selected=9.309754891 ppm",
            summary)
        @test occursin("0.980000000000 0.937000000000", summary)
        @test occursin("0.959984105853 0.919969933271", summary)
    end
end
