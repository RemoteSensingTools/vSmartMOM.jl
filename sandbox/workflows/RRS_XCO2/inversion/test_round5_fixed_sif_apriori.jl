#!/usr/bin/env julia

using LinearAlgebra
using NCDatasets
using Test

include(joinpath(
    @__DIR__, "retrieval_setup", "build_round5_fixed_sif_apriori.jl"))
using .Round5FixedSIFPriorBuilder

include(joinpath(@__DIR__, "RetrievalState.jl"))
using .RetrievalState: load_retrieval_prior

const SOURCE_PRIOR = joinpath(
    @__DIR__, "..", "bottom_layer_XCO2_retrievals", "retrieval_setup",
    "apriori_states_acos_mapped_tapered_vertical_correlation.nc")

@testset "round-5 fixed-SIF and tight-UTLS prior" begin
    isfile(SOURCE_PRIOR) || error("missing test source prior: $SOURCE_PRIOR")
    NCDataset(SOURCE_PRIOR, "r") do source_dataset
        source_xa = Float64.(source_dataset["xa"][:, :])
        source_Sa = Float64.(source_dataset["Sa"][:, :, :])

        for mode in (:off, :on)
            prior = build_round5_prior(SOURCE_PRIOR, mode)
            @test prior.active_to_full == ROUND5_ACTIVE_TO_FULL
            @test prior.active_to_core == ROUND5_ACTIVE_TO_CORE
            @test size(prior.Sa_active) == (28, 28, 4)
            @test all(surface -> isposdef(Symmetric(
                prior.Sa_active[:, :, surface])), 1:4)

            # Apples-to-apples round-4 comparison: the complete CO2 mean and
            # covariance—including fixed upper layers—must be bit-identical.
            @test prior.xa[2:17, :] == source_xa[2:17, :]
            @test prior.Sa[2:17, 2:17, :] == source_Sa[2:17, 2:17, :]

            # Only the third aerosol species is tightened. Sigma / 10 means
            # variance / 100; both tropospheric species remain bit-identical.
            @test prior.Sa[18:19, 18:19, :] ==
                source_Sa[18:19, 18:19, :]
            @test prior.Sa[21:22, 21:22, :] ==
                source_Sa[21:22, 21:22, :]
            @test prior.aerosol_aod_sigma ≈ [0.75, 0.75, 0.075]
            @test prior.aerosol_z0_sigma ≈ [0.10, 0.10, 0.01]
            @test prior.Sa[20, 20, :] ≈
                source_Sa[20, 20, :] .* 0.01
            @test prior.Sa[23, 23, :] ≈
                source_Sa[23, 23, :] .* 0.01

            @test all(iszero, prior.Sa[33:34, :, :])
            @test all(iszero, prior.Sa[:, 33:34, :])
            if mode == :off
                @test all(iszero, prior.xa[33:34, :])
            else
                @test all(prior.xa[33, :] .== prior.fixed_sif.SIF760)
                @test all(prior.xa[34, :] .== prior.fixed_sif.mSIF)
                @test prior.fixed_sif.mSIF == 1.2291230681458325e-5
            end

            mktempdir() do temporary
                nc = joinpath(temporary, output_filename(mode))
                dat = joinpath(temporary,
                    output_filename(mode; extension="dat"))
                write_round5_prior(prior; output_path=nc)
                write_round5_summary(prior; output_path=dat)
                @test_throws ArgumentError write_round5_prior(
                    prior; output_path=nc)
                @test occursin(
                    "CO2 prior: mean and full 16x16 covariance copied exactly",
                    read(dat, String))

                NCDataset(nc, "r") do output
                    @test output.attrib["retrieval_state_model"] ==
                        "round5_fixed_sif"
                    @test output.attrib["round5_prior_model"] == ROUND5_MODEL
                    @test output.attrib["round5_sif_case"] == String(mode)
                    @test output.attrib["active_state_count"] == 28
                    @test Int.(output["active_parameter_index"][:]) ==
                        ROUND5_ACTIVE_TO_FULL
                    @test Int.(output["active_core_parameter_index"][:]) ==
                        ROUND5_ACTIVE_TO_CORE
                    @test Float64.(output["xa"][2:17, :]) ==
                        source_xa[2:17, :]
                    @test Float64.(output["Sa"][2:17, 2:17, :]) ==
                        source_Sa[2:17, 2:17, :]
                end

                for surface in (:urban, :rural, :desert, :forest)
                    loaded = load_retrieval_prior(surface; path=nc)
                    @test length(loaded.xa) == 28
                    @test loaded.active_to_full == ROUND5_ACTIVE_TO_FULL
                    @test isposdef(Symmetric(loaded.Sa))
                    @test all(name -> name ∉ ("SIF760", "mSIF"),
                              loaded.parameter_names)
                end
            end
        end
    end
end
