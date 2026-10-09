#!/usr/bin/env julia

using LinearAlgebra
using NCDatasets
using Test

include(joinpath(@__DIR__, "OptimalEstimation.jl"))
include(joinpath(@__DIR__, "Round4SIFTruthConvention.jl"))
include(joinpath(@__DIR__, "Round5FixedSIF.jl"))
include(joinpath(@__DIR__, "Round5RetrievalCampaign.jl"))
include(joinpath(
    @__DIR__, "retrieval_setup", "build_round5_fixed_sif_apriori.jl"))
include(joinpath(@__DIR__, "RetrievalState.jl"))

using .Round4SIFTruthConvention
using .Round5FixedSIF
using .Round5RetrievalCampaign
using .Round5FixedSIFPriorBuilder
using .RetrievalState: load_retrieval_prior

const SOURCE_PRIOR = joinpath(
    @__DIR__, "..", "bottom_layer_XCO2_retrievals", "retrieval_setup",
    "apriori_states_acos_mapped_tapered_vertical_correlation.nc")

@testset "round-5 campaign mode and prior contracts" begin
    off_scene = (sif_case=:off,)
    on_scene = (sif_case=:angular_integral760_0p5,)
    off = resolve_round5_sif_map([off_scene])
    on = resolve_round5_sif_map([on_scene]; requested="on")
    @test !off.sif_on
    @test on.sif_on
    @test active_state_count(off) == active_state_count(on) == 28
    @test active_core_indices(off) == active_core_indices(on) == collect(1:28)
    @test expected_active_to_full(off) == vcat(1, collect(6:32))
    @test expected_active_to_full(on) == expected_active_to_full(off)
    @test on.Lnu759 == truth_sif_at_wavelength_nm(759.0).Lnu
    @test on.mSIF == validate_round4_sif_truth_convention().mSIF
    @test_throws ArgumentError resolve_round5_sif_map([off_scene, on_scene])
    @test_throws ArgumentError resolve_round5_sif_map(
        [off_scene]; requested="on")
    @test_throws ArgumentError resolve_round5_sif_map(
        [on_scene]; known_wavelength_nm=760)

    mktempdir() do directory
        for (mode, map) in ((:off, off), (:on, on))
            prior_record = build_round5_prior(SOURCE_PRIOR, mode)
            path = joinpath(directory, output_filename(mode))
            write_round5_prior(prior_record; output_path=path)
            for surface in (:urban, :rural, :desert, :forest)
                prior = load_retrieval_prior(surface; path)
                @test validate_round5_prior(prior, map; path) === prior
            end

            identity = Dict(
                "round5_campaign_identity_status" => "complete",
                "round5_code_checkpoint_sha" => repeat("a", 40),
                "round5_codeset_sha256" => repeat("b", 64),
                "round5_input_set_sha256" => repeat("c", 64),
                "round5_campaign_identity_sha256" => repeat("d", 64),
            )
            provenance = round5_output_provenance(
                map; prior_path=path, campaign_identity=identity)
            provenance["jacobian_flavor"] = ROUND5_JACOBIAN_FLAVOR
            provenance["state_dimension"] = 28
            @test validate_round5_output_provenance(
                provenance, map; prior_path=path,
                campaign_identity=identity)
        end
    end
end

@testset "round-5 identity and realization contracts" begin
    environment = Dict(
        "ROUND5_CODE_CHECKPOINT_SHA" => repeat("a", 40),
        "ROUND5_CODESET_SHA256" => repeat("b", 64),
        "ROUND5_INPUT_SET_SHA256" => repeat("c", 64),
        "ROUND5_CAMPAIGN_IDENTITY_SHA256" => repeat("d", 64),
    )
    identity = round5_campaign_identity(environment)
    @test identity["round5_campaign_identity_status"] == "complete"
    @test_throws ErrorException round5_campaign_identity(
        Dict("ROUND5_CODE_CHECKPOINT_SHA" => repeat("a", 40)))
    @test round5_campaign_identity(Dict(); required=false)[
        "round5_campaign_identity_status"] == "not_provided_nonproduction"

    off = Round5SIFMap(false)
    @test validate_round5_realization(
        off, (sif_case=:off,), (provenance=Dict{String,Any}(),)) !== nothing
    convention = validate_round4_sif_truth_convention()
    on = Round5SIFMap(true;
        Lnu759=convention.diagnostic.Lnu, mSIF=convention.mSIF)
    realization = (provenance=
        expected_round4_sif_provenance(; enabled=true),)
    @test validate_round5_realization(
        on, (sif_case=:angular_integral760_0p5,), realization) === realization
    @test_throws ErrorException validate_round5_realization(
        on, (sif_case=:off,), realization)

    @test endswith(validate_round5_output_root(
        "/tmp/round5_fixed_sif/results"), "/tmp/round5_fixed_sif/results")
    @test_throws ArgumentError validate_round5_output_root(
        "/tmp/round4_known_sif759/results")
end
