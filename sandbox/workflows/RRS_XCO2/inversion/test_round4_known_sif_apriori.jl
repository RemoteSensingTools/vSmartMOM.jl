#!/usr/bin/env julia

using LinearAlgebra
using NCDatasets
using Test

module Round3PriorFixture
include(joinpath(@__DIR__, "retrieval_setup", "build_apriori.jl"))
end

include(joinpath(
    @__DIR__, "retrieval_setup", "build_round4_known_sif_apriori.jl"))
using .Round4KnownSIFPriorBuilder

include(joinpath(@__DIR__, "RetrievalState.jl"))
using .RetrievalState: load_retrieval_prior

function write_round3_fixture(path)
    priors = Dict(
        surface => Round3PriorFixture.build_prior(
            surface, 0.1;
            co2_covariance_model=
                Round3PriorFixture.TAPERED_CO2_COVARIANCE_MODEL,
            sif_wavelength_slope_mean=0.0,
            sif_wavelength_slope_sigma=0.002625,
            aerosol_ln_aod_sigma=0.75,
            surface_p1_sigmas=(0.002, 0.002, 0.002),
            surface_p2_sigmas=(0.002, 0.002, 0.002))
        for surface in Round3PriorFixture.SURFACES)
    Round3PriorFixture.write_netcdf(priors; output_path=path)
    return path
end

@testset "round-4 known-SIF prior builder" begin
    mktempdir() do temporary
        source_path = write_round3_fixture(joinpath(temporary, "round3.nc"))
        source_sha = file_sha256(source_path)
        NCDataset(source_path, "r") do source
            source_xa = Float64.(source["xa"][:, :])
            source_Sa = Float64.(source["Sa"][:, :, :])

            for mode in (:off, :on)
                prior = build_round4_prior(source_path, mode)
                expected_full = mode == :on ?
                    ROUND4_ON_ACTIVE_TO_FULL : ROUND4_OFF_ACTIVE_TO_FULL
                expected_core = mode == :on ?
                    vcat(collect(1:28), 30) : collect(1:28)

                @test prior.active_to_full == expected_full
                @test prior.active_to_core == expected_core
                @test length(prior.active_to_full) == (mode == :on ? 29 : 28)
                @test prior.xa[1:32, :] == source_xa[1:32, :]
                @test prior.Sa[1:32, 1:32, :] == source_Sa[1:32, 1:32, :]
                @test all(surface -> isposdef(Symmetric(
                    prior.Sa_active[:, :, surface])), axes(prior.Sa_active, 3))

                output_path = joinpath(
                    temporary, output_filename(mode; extension="nc"))
                summary_path = joinpath(
                    temporary, output_filename(mode; extension="dat"))
                write_round4_prior(prior; output_path)
                write_round4_summary(prior; output_path=summary_path)
                @test_throws ArgumentError write_round4_prior(
                    prior; output_path)
                @test occursin("Source prior SHA-256: $source_sha",
                               read(summary_path, String))

                NCDataset(output_path, "r") do output
                    @test output.attrib["retrieval_state_model"] ==
                        "round4_known_sif759"
                    @test output.attrib["round4_prior_model"] == ROUND4_MODEL
                    @test output.attrib["round4_sif_case"] == String(mode)
                    @test output.attrib["source_prior_sha256"] == source_sha
                    @test occursin("do not Gaussian-condition",
                        String(output.attrib[
                            "round4_production_prior_policy"]))
                    @test output.attrib["active_state_count"] ==
                        length(expected_full)
                    @test Int.(output["active_parameter_index"][:]) ==
                        expected_full
                    @test Int.(output["active_core_parameter_index"][:]) ==
                        expected_core
                    @test Float64.(output["xa"][1:32, :]) ==
                        source_xa[1:32, :]
                    @test Float64.(output["Sa"][1:32, 1:32, :]) ==
                        source_Sa[1:32, 1:32, :]
                    @test size(output["Sa_active"]) ==
                        (length(expected_full), length(expected_full), 4)
                    @test output.attrib[
                        "sif_wavelength_slope_sigma_mw_m2_sr_nm2"] ==
                        (mode == :on ? 0.002625 : 0.0)
                end

                for surface in (:urban, :rural, :desert, :forest)
                    loaded = load_retrieval_prior(surface; path=output_path)
                    @test length(loaded.xa) == length(expected_full)
                    @test loaded.active_to_full == expected_full
                    @test isposdef(Symmetric(loaded.Sa))
                end
            end
        end
    end
end

@testset "round-4 SIF slope semantics" begin
    convention =
        Round4KnownSIFPriorBuilder.Round4SIFTruthConvention.
            validate_round4_sif_truth_convention()
    Lnu759 = convention.diagnostic.Lnu
    slope = native_slope_prior(Lnu759, 0.0, 0.002625)

    @test isapprox(slope.mean, -7.34275712494799e-7;
                   atol=2e-21, rtol=2e-14)
    @test isapprox(slope.sigma, 8.780708772523119e-6;
                   atol=2e-20, rtol=2e-14)
    @test isapprox(slope.SIF760_mean, 0.004830761266413151;
                   atol=3e-18, rtol=2e-14)

    # Transform the native one-dimensional prior back to the physical
    # wavelength-space slope at 760 nm. It must reproduce the requested
    # current constraint exactly, rather than the old uncertain-amplitude
    # marginal or a conditionally shifted slope.
    recovered_mean = slope.wavelength_slope_alpha * slope.mean +
                     slope.wavelength_slope_beta
    recovered_sigma = abs(slope.wavelength_slope_alpha) * slope.sigma
    @test isapprox(recovered_mean, 0.0; atol=2e-19, rtol=0)
    @test isapprox(recovered_sigma, 0.002625;
                   atol=2e-18, rtol=2e-14)

    mktempdir() do temporary
        source_path = write_round3_fixture(joinpath(temporary, "round3.nc"))
        on = build_round4_prior(source_path, :on)
        delta_nu = on.native_slope.delta_nu_760_minus_759
        @test isapprox(on.xa[33, 1] - delta_nu * on.xa[34, 1],
                       on.known_Lnu759; atol=3e-18, rtol=0)
        @test on.Sa[33:34, 33:34, 1] ==
            on.native_slope.variance .* [delta_nu^2 delta_nu; delta_nu 1.0]

        # Conditioning the old uncertain-amplitude two-coordinate prior is a
        # different scientific assumption. Confirm it does not silently
        # reproduce the dedicated production construction.
        NCDataset(source_path, "r") do source
            old_mean = Float64.(source["xa"][33:34, 1])
            old_covariance = Float64.(source["Sa"][33:34, 33:34, 1])
            nu_759 = 1.0e7 / 759.0
            nu_760 = 1.0e7 / 760.0
            H = [1.0, nu_759 - nu_760]
            cross = old_covariance * H
            conditioned_mean = old_mean + cross *
                ((on.known_Lnu759 - dot(H, old_mean)) /
                 dot(H, cross))
            conditioned_variance = old_covariance -
                (cross * transpose(cross)) / dot(H, cross)
            @test !isapprox(conditioned_mean[2], on.native_slope.mean;
                            atol=0, rtol=1e-6)
            @test !isapprox(sqrt(conditioned_variance[2, 2]),
                            on.native_slope.sigma; atol=0, rtol=1e-6)
        end

        off = build_round4_prior(source_path, :off)
        @test all(iszero, off.xa[33:34, :])
        @test all(iszero, off.Sa[33:34, :, :])
        @test all(iszero, off.Sa[:, 33:34, :])
        @test off.wavelength_slope_mean == 0.0
        @test off.wavelength_slope_sigma == 0.0

        coupled_path = write_round3_fixture(
            joinpath(temporary, "round3_with_unexpected_sif_coupling.nc"))
        NCDataset(coupled_path, "a") do coupled
            coupled["Sa"][33, 1, 1] = 1e-12
            coupled["Sa"][1, 33, 1] = 1e-12
        end
        @test_throws ErrorException build_round4_prior(coupled_path, :on)
    end
end
