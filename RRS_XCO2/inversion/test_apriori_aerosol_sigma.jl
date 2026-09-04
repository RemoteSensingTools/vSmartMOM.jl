#!/usr/bin/env julia

using LinearAlgebra
using NCDatasets
using Test

include(joinpath(@__DIR__, "retrieval_setup", "build_apriori.jl"))
include(joinpath(@__DIR__, "RetrievalState.jl"))
include(joinpath(@__DIR__, "preflight_bottom_layer_aerosol_retrievals.jl"))
using .RetrievalState

@testset "isolated aerosol ln(AOD) prior width" begin
    legacy = build_prior(:urban, 0.1; aerosol_ln_aod_sigma=2.0)
    current = build_prior(:urban, 0.1)

    @test DEFAULT_AEROSOL_LN_AOD_SIGMA == 0.75
    @test current.xa == legacy.xa
    @test sqrt.(diag(current.Sa))[18:20] == fill(0.75, 3)

    expected_difference = zeros(34, 34)
    for index in 18:20
        expected_difference[index, index] = 0.75^2 - 2.0^2
    end
    @test current.Sa - legacy.Sa == expected_difference
    @test_throws ArgumentError build_prior(
        :urban, 0.1; aerosol_ln_aod_sigma=0.0)
    @test_throws ArgumentError build_prior(
        :urban, 0.1; aerosol_ln_aod_sigma=Inf)

    mktempdir() do directory
        netcdf_path = joinpath(directory, "apriori_states.nc")
        summary_path = joinpath(directory, "apriori_states.dat")
        priors = Dict(
            surface => build_prior(
                surface, 0.1;
                sif_wavelength_slope_mean=0.0,
                sif_wavelength_slope_sigma=0.002625,
                surface_p1_sigmas=ntuple(_ -> 0.002, 3),
                surface_p2_sigmas=ntuple(_ -> 0.002, 3))
            for surface in SURFACES)
        write_netcdf(priors; output_path=netcdf_path)
        write_summary(priors; output_path=summary_path)

        prior = load_retrieval_prior(:urban; path=netcdf_path)
        @test prior.active_to_full[14:16] == collect(18:20)
        @test prior.parameter_names[14:16] == [
            "ln_sulfate_aod760",
            "ln_organic_carbon_aod760",
            "ln_utls_sulfate_aod760",
        ]
        @test diag(prior.Sa)[14:16] == fill(0.75^2, 3)

        NCDataset(netcdf_path) do dataset
            @test Float64(dataset.attrib["aerosol_ln_aod_sigma"]) == 0.75
            @test Float64.(dataset["prior_sigma"][18:20, :]) ==
                fill(0.75, 3, 4)
            @test Float64.(dataset["Sa_active"][14:16, 14:16, 1]) ==
                Diagonal(fill(0.75^2, 3))
        end
        @test isnothing(validate_prior(netcdf_path))

        NCDataset(netcdf_path, "a") do dataset
            dataset.attrib["aerosol_ln_aod_sigma"] = 2.0
        end
        @test_throws ErrorException validate_prior(netcdf_path)
        NCDataset(netcdf_path, "a") do dataset
            dataset.attrib["aerosol_ln_aod_sigma"] = 0.75
            dataset["Sa_active"][14, 14, 4] = 4.0
        end
        @test_throws ErrorException validate_prior(netcdf_path)

        @test occursin(
            "Aerosol ln(AOD760) sigma: 7.500000000000e-01",
            read(summary_path, String))
    end
end
