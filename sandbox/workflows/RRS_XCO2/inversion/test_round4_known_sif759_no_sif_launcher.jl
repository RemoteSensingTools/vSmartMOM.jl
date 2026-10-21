#!/usr/bin/env julia

using LinearAlgebra
using NCDatasets
using Test

const ROUND4_LAUNCHER = joinpath(
    @__DIR__, "run_round4_known_sif759_no_sif_partition.sh")
const ROUND4_PREFLIGHT = joinpath(
    @__DIR__, "preflight_round4_known_sif759_no_sif_retrievals.jl")

include(ROUND4_PREFLIGHT)

module Round3PriorFixture
include(joinpath(@__DIR__, "retrieval_setup", "build_apriori.jl"))
end

include(joinpath(
    @__DIR__, "retrieval_setup", "build_round4_known_sif_apriori.jl"))
using .Round4KnownSIFPriorBuilder

function launcher_plan(worker; extra=Dict{String,String}())
    environment = merge(copy(ENV), Dict(
        "ROUND4_NOSIF_PLAN_ONLY" => "1",
        "ROUND4_NOSIF_EXECUTE" => "0",
    ), extra)
    return read(setenv(`bash $ROUND4_LAUNCHER $worker`, environment), String)
end

function failed_launcher(arguments...)
    environment = merge(copy(ENV), Dict(
        "ROUND4_NOSIF_PLAN_ONLY" => "1",
        "ROUND4_NOSIF_EXECUTE" => "0",
    ))
    command = setenv(`bash $ROUND4_LAUNCHER $arguments`, environment)
    return run(pipeline(ignorestatus(command), stdout=devnull,
                        stderr=devnull)).exitcode
end

function plan_value(plan, name)
    prefix = name * "="
    line = only(filter(line -> startswith(line, prefix), split(plan, '\n')))
    return line[length(prefix) + 1:end]
end

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

function write_round4_off_fixture(source_path, path)
    prior = build_round4_prior(source_path, :off)
    write_round4_prior(prior; output_path=path)
    return path
end

@testset "round-4 local schedule is exact and disjoint" begin
    source = read(ROUND4_LAUNCHER, String)
    preflight_source = read(ROUND4_PREFLIGHT, String)

    @test occursin("ROUND4_NOSIF_EXECUTE:-0", source)
    @test occursin("ROUND4_SIF_MODE=off", source)
    @test occursin("ROUND4_KNOWN_WAVELENGTH_NM=759", source)
    @test occursin("SIF_CASE_FILTER=off", source)
    @test occursin("RETRIEVAL_CLASS=paired", source)
    @test occursin("FIRST_PERTURBATION=1", source)
    @test occursin("LAST_PERTURBATION=11", source)
    @test occursin("FORCE=0", source)
    @test occursin("RETRIEVAL_WRITE_MANIFEST=0", source)
    @test occursin(".state_claims", source)
    @test occursin("owner.dat", source)
    @test occursin("require_assigned_gpu_idle", source)
    @test occursin("for block in \"\${STATE_BLOCKS[@]}\"", source)
    @test occursin("ROUND4_CODESET_SHA256", source)
    @test occursin("ROUND4_INPUT_SET_SHA256", source)
    @test occursin("campaign_identity.dat", preflight_source)
    @test occursin("round4_input_set_sha256", preflight_source)

    plans = Dict(worker => launcher_plan(worker) for worker in
                 ("curry1", "wurst0", "wurst1"))
    @test plan_value(plans["curry1"], "physical_gpu") == "1"
    @test plan_value(plans["curry1"], "scene_class") == "none"
    @test plan_value(plans["curry1"], "blocks") ==
        "1:5,21:25,41:45,61:65"
    @test plan_value(plans["wurst0"], "physical_gpu") == "0"
    @test plan_value(plans["wurst0"], "scene_class") == "aerosol"
    @test plan_value(plans["wurst0"], "blocks") == "11:15,51:55"
    @test plan_value(plans["wurst1"], "physical_gpu") == "1"
    @test plan_value(plans["wurst1"], "scene_class") == "aerosol"
    @test plan_value(plans["wurst1"], "blocks") == "31:35,71:75"

    assignments = Dict(worker => parse.(Int,
        split(plan_value(plan, "states"), ',')) for (worker, plan) in plans)
    @test all(length(intersect(assignments[left], assignments[right])) == 0
              for (left, right) in (("curry1", "wurst0"),
                                     ("curry1", "wurst1"),
                                     ("wurst0", "wurst1")))
    @test sort(vcat(values(assignments)...)) ==
        parse_state_spec(ROUND4_GLOBAL_STATE_SPEC)
    @test all(length(values) == length(unique(values))
              for values in values(assignments))

    @test failed_launcher("curry0") != 0
    @test failed_launcher("curry1", "1-5") != 0
    @test failed_launcher("wurst0", "31-35") != 0
    @test failed_launcher("unknown") != 0
end

@testset "round-4 namespaces and campaign-aware truth paths" begin
    campaign = "/campaign"
    round4 = joinpath(campaign, "round4_known_sif759")
    prior = joinpath(round4, "retrieval_setup", "prior.nc")
    source = joinpath(campaign, "retrieval_setup", "round3.nc")
    output = joinpath(round4, "retrievals_nosif")
    @test isnothing(validate_path_isolation(
        campaign, prior, source, output))
    @test_throws ErrorException validate_path_isolation(
        campaign, source, source, output)
    @test_throws ErrorException validate_path_isolation(
        campaign, prior, source, joinpath(campaign, "retrievals"))
    @test_throws ErrorException validate_path_isolation(
        campaign, prior, source,
        joinpath(campaign,
            "retrievals_acos_mapped_tapered_vertical_correlation_nosif"))
    @test_throws ErrorException validate_path_isolation(
        campaign, prior, source, joinpath(round4, "retrievals"))

    clear = (state_index=1, aerosol_case=:none)
    aerosol = (state_index=11, aerosol_case=:aod760_0p28)
    @test truth_scene_path(clear, "/campaign/truth") ==
        "/campaign/truth/hiressim_001.nc"
    @test truth_scene_path(aerosol, "/campaign/truth") ==
        "/campaign/truth/aerosol_chunked/hiressim_011.nc"
end

@testset "round-4 no-SIF prior is pinned and invariant" begin
    mktempdir() do temporary
        source_path = write_round3_fixture(joinpath(temporary, "round3.nc"))
        round4_path = write_round4_off_fixture(
            source_path, joinpath(temporary, "round4_off.nc"))
        source_sha = file_sha256(source_path)
        round4_sha = file_sha256(round4_path)

        identity = validate_round4_nosif_prior(
            round4_path; expected_sha256=round4_sha,
            source_prior_path=source_path,
            source_prior_sha256=source_sha)
        @test identity.prior_sha256 == round4_sha
        @test identity.source_prior_sha256 == source_sha
        @test_throws ErrorException validate_round4_nosif_prior(
            round4_path; expected_sha256=repeat("0", 64),
            source_prior_path=source_path,
            source_prior_sha256=source_sha)
        @test_throws ErrorException validate_round4_nosif_prior(
            round4_path; expected_sha256=round4_sha,
            source_prior_path=source_path,
            source_prior_sha256=repeat("0", 64))

        # Non-adjacent CO2/non-SIF corruption is rejected even if Sa_active is
        # changed in lockstep to preserve internal indexing consistency.
        NCDataset(round4_path, "a") do dataset
            value = Float64(dataset["Sa"][6, 8, 1]) + 1e-14
            dataset["Sa"][6, 8, 1] = value
            dataset["Sa"][8, 6, 1] = value
            dataset["Sa_active"][2, 4, 1] = value
            dataset["Sa_active"][4, 2, 1] = value
        end
        @test_throws ErrorException validate_round4_nosif_prior(
            round4_path; expected_sha256=file_sha256(round4_path),
            source_prior_path=source_path,
            source_prior_sha256=source_sha)

        round4_path = joinpath(temporary, "round4_off_clean.nc")
        write_round4_off_fixture(source_path, round4_path)
        NCDataset(round4_path, "a") do dataset
            dataset["Sa"][33, 33, 1] = 1e-10
        end
        @test_throws ErrorException validate_round4_nosif_prior(
            round4_path; expected_sha256=file_sha256(round4_path),
            source_prior_path=source_path,
            source_prior_sha256=source_sha)

        output_root = joinpath(temporary, "round4", "retrievals_nosif")
        expected = "identity_schema=1\ncampaign_id=round4-a\n"
        path = initialize_campaign_identity!(output_root, expected)
        @test read(path, String) == expected
        @test initialize_campaign_identity!(output_root, expected) == path
        @test_throws ErrorException initialize_campaign_identity!(
            output_root, "identity_schema=1\ncampaign_id=round4-b\n")

        good_output_identity = Dict(
            "round4_campaign_identity_status" => "complete",
            "round4_code_checkpoint_sha" => repeat("a", 40),
            "round4_codeset_sha256" => repeat("b", 64),
            "round4_input_set_sha256" => repeat("c", 64),
            "round4_campaign_identity_sha256" => repeat("d", 64),
        )
        @test validate_output_campaign_identity(good_output_identity;
            code_checkpoint=repeat("a", 40),
            codeset_sha256=repeat("b", 64),
            input_set_sha256=repeat("c", 64),
            campaign_identity_sha256=repeat("d", 64))
        mixed = copy(good_output_identity)
        mixed["round4_input_set_sha256"] = repeat("e", 64)
        @test_throws ErrorException validate_output_campaign_identity(mixed;
            code_checkpoint=repeat("a", 40),
            codeset_sha256=repeat("b", 64),
            input_set_sha256=repeat("c", 64),
            campaign_identity_sha256=repeat("d", 64))
    end
end
