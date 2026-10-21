#!/usr/bin/env julia

using LinearAlgebra
using NCDatasets
using Test

include(joinpath(@__DIR__, "OptimalEstimation.jl"))
include(joinpath(@__DIR__, "Round4KnownSIF.jl"))
include(joinpath(@__DIR__, "Round4SIFTruthConvention.jl"))
include(joinpath(@__DIR__, "Round4RetrievalCampaign.jl"))

using .Round4KnownSIF
using .Round4SIFTruthConvention
using .Round4RetrievalCampaign

function write_prior_identity(path, map)
    NCDataset(path, "c") do dataset
        defDim(dataset, "active", active_state_count(map))
        core_index = defVar(dataset, "active_core_parameter_index", Int16,
                            ("active",))
        core_index[:] = Int16.(active_core_indices(map))
        dataset.attrib["retrieval_state_model"] = ROUND4_STATE_MODEL
        dataset.attrib["round4_prior_model"] =
            "round4_known_sif759_acos_mapped_tapered_vertical_correlation"
        dataset.attrib["round4_prior_definition_version"] = Int32(1)
        dataset.attrib["round4_sif_case"] = map.sif_on ? "on" : "off"
        dataset.attrib["known_sif_wavelength_nm"] = 759.0
        dataset.attrib["known_sif_Lnu_mW_m-2_sr-1_per_cm-1"] = map.Lν759
        dataset.attrib["known_sif_Llambda_mW_m-2_sr-1_nm-1"] =
            map.sif_on ? lnu_to_llambda(map.Lν759, 759.0) : 0.0
        dataset.attrib["active_state_count"] = active_state_count(map)
        dataset.attrib["co2_covariance_model"] =
            "acos_mapped_tapered_vertical_correlation"
        dataset.attrib["round4_active_to_full"] =
            join(expected_active_to_full(map), " ")
        dataset.attrib["round4_active_to_core"] =
            join(active_core_indices(map), " ")
        dataset.attrib["source_prior_sha256"] = repeat("a", 64)
        dataset.attrib["source_prior_path"] = "/source/round3_prior.nc"
    end
end

function fake_prior(map)
    n = active_state_count(map)
    names = ["parameter_$index" for index in 1:n]
    map.sif_on && (names[end] = "mSIF")
    return (;
        xa=zeros(n),
        Sa=Matrix{Float64}(I, n, n),
        active_to_full=expected_active_to_full(map),
        parameter_names=names,
    )
end

@testset "round-4 retrieval mode selection" begin
    off_scene = (sif_case=:off,)
    on_scene = (sif_case=:angular_integral760_0p5,)

    off = resolve_round4_sif_map([off_scene])
    on = resolve_round4_sif_map([on_scene]; requested="on")
    @test !off.sif_on
    @test iszero(off.Lν759)
    @test on.sif_on
    @test on.Lν759 == truth_sif_at_wavelength_nm(759.0).Lnu
    @test active_state_count(off) == 28
    @test active_state_count(on) == 29
    @test expected_active_to_full(off) == vcat(1, collect(6:32))
    @test expected_active_to_full(on) == vcat(1, collect(6:32), 34)

    @test_throws ArgumentError resolve_round4_sif_map(
        [off_scene, on_scene])
    @test_throws ArgumentError resolve_round4_sif_map(
        [off_scene]; requested="on")
    @test_throws ArgumentError resolve_round4_sif_map(
        [on_scene]; known_wavelength_nm=760)
end

@testset "round-4 prior and output identity" begin
    identity_environment = Dict(
        "ROUND4_CODE_CHECKPOINT_SHA" => repeat("a", 40),
        "ROUND4_CODESET_SHA256" => repeat("b", 64),
        "ROUND4_INPUT_SET_SHA256" => repeat("c", 64),
        "ROUND4_CAMPAIGN_IDENTITY_SHA256" => repeat("d", 64),
    )
    identity = round4_campaign_identity(identity_environment)
    @test identity["round4_campaign_identity_status"] == "complete"
    @test identity["round4_code_checkpoint_sha"] == repeat("a", 40)
    @test_throws ErrorException round4_campaign_identity(
        Dict("ROUND4_CODE_CHECKPOINT_SHA" => repeat("a", 40)))
    @test_throws ErrorException round4_campaign_identity(merge(
        identity_environment,
        Dict("ROUND4_CODE_CHECKPOINT_SHA" => repeat("A", 40))))
    @test round4_campaign_identity(Dict(); required=false)[
        "round4_campaign_identity_status"] == "not_provided_nonproduction"

    mktempdir() do directory
        for map in (Round4SIFMap(false),
                    Round4SIFMap(true;
                        Lν759=truth_sif_at_wavelength_nm(759.0).Lnu))
            path = joinpath(directory, map.sif_on ? "on.nc" : "off.nc")
            write_prior_identity(path, map)
            prior = fake_prior(map)
            @test validate_round4_prior(prior, map; path) === prior

            provenance = round4_output_provenance(
                map; prior_path=path, campaign_identity=identity)
            if !map.sif_on
                independent = round4_output_provenance(
                    map; prior_path=path,
                    campaign_identity=identity,
                    convention_loader=() -> error(
                        "SIF truth resource must not be read"))
                @test independent == provenance
                @test provenance[
                    "round4_known_sif_Lnu_mW_m-2_sr-1_per_cm-1"] == 0.0
                @test provenance["round4_canonical_mSIF_fixed_mean"] == 0.0
                @test provenance[
                    "round4_canonical_mSIF_fixed_variance"] == 0.0
                @test provenance[
                    "round4_canonical_mSIF_excluded_from_active_solve"] == 1
            end
            provenance["jacobian_flavor"] = ROUND4_JACOBIAN_FLAVOR
            provenance["state_dimension"] = active_state_count(map)
            @test validate_round4_output_provenance(
                provenance, map; prior_path=path,
                campaign_identity=identity)

            wrong_prior = merge(prior, (
                active_to_full=collect(1:active_state_count(map)),))
            @test_throws ErrorException validate_round4_prior(
                wrong_prior, map; path)
            wrong_output = copy(provenance)
            wrong_output["round4_known_sif_wavelength_nm"] = 760.0
            @test_throws ErrorException validate_round4_output_provenance(
                wrong_output, map; prior_path=path,
                campaign_identity=identity)
        end
    end

    @test endswith(
        validate_round4_output_root("/tmp/round4_known_sif759/results"),
        "/tmp/round4_known_sif759/results")
    @test_throws ArgumentError validate_round4_output_root(
        "/tmp/round3/results")
end

@testset "round-4 truth-realization consistency" begin
    off = Round4SIFMap(false)
    off_truth = (sif_case=:off,)
    @test validate_round4_realization(
        off, off_truth, (provenance=Dict{String,Any}(),)) !== nothing

    convention = validate_round4_sif_truth_convention()
    on = Round4SIFMap(true; Lν759=convention.diagnostic.Lnu)
    on_truth = (sif_case=:angular_integral760_0p5,)
    on_realization = (
        provenance=expected_round4_sif_provenance(; enabled=true),)
    @test validate_round4_realization(
        on, on_truth, on_realization) === on_realization
    @test_throws ErrorException validate_round4_realization(
        on, off_truth, on_realization)
end
