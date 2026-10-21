#!/usr/bin/env julia

using Test
using NCDatasets

const LAUNCHER = joinpath(
    @__DIR__, "run_tapered_no_sif_retrieval_partition.sh")
const PREFLIGHT = joinpath(
    @__DIR__, "preflight_tapered_no_sif_retrievals.jl")

include(PREFLIGHT)

module PriorBuilder
include(joinpath(@__DIR__, "retrieval_setup", "build_apriori.jl"))
end

function plan(worker, states; extra=Dict{String,String}())
    environment = copy(ENV)
    merge!(environment, Dict(
        "TAPERED_NOSIF_PLAN_ONLY" => "1",
        "TAPERED_NOSIF_EXECUTE" => "0",
    ))
    merge!(environment, extra)
    command = setenv(`bash $LAUNCHER $worker $states`, environment)
    return read(command, String)
end

function failed_plan(worker, states)
    command = setenv(
        `bash $LAUNCHER $worker $states`,
        merge(copy(ENV), Dict("TAPERED_NOSIF_PLAN_ONLY" => "1")))
    process = run(pipeline(ignorestatus(command), stdout=devnull,
                           stderr=devnull))
    return process.exitcode
end

@testset "isolated tapered no-SIF launcher contract" begin
    source = read(LAUNCHER, String)
    preflight_source = read(PREFLIGHT, String)

    @test occursin("TAPERED_NOSIF_EXECUTE:-0", source)
    @test occursin("SIF_CASE_FILTER=off", source)
    @test !occursin("SIF_CASE_FILTER=\"\${", source)
    @test occursin("FORCE=0", source)
    @test occursin("RETRIEVAL_WRITE_MANIFEST=0", source)
    @test occursin(".state_claims", source)
    @test occursin("TAPERED_PRIOR_SHA256", source)
    @test occursin("CURRENT_PRIOR_SHA256", source)
    @test occursin("campaign_identity.dat", preflight_source)
    @test occursin("apriori_states_\${MODEL_NAME}.nc", source)
    @test occursin("retrievals_\${MODEL_NAME}_nosif", source)
    @test occursin("validate_existing_outputs", preflight_source)
    @test occursin("acos_mapped_tapered_vertical_correlation",
                   preflight_source)

    curry = plan("curry0", "1-5,21-25")
    @test occursin("expected_host=curry", curry)
    @test occursin("physical_gpu=0", curry)
    @test occursin("states=1,2,3,4,5,21,22,23,24,25", curry)
    @test occursin("sif_case_filter=off", curry)
    @test occursin("execute=0", curry)

    wurst = plan("wurst1", "11-15,31-35,41")
    @test occursin("expected_host=wurst", wurst)
    @test occursin("physical_gpu=1", wurst)
    @test occursin("states=11,12,13,14,15,31,32,33,34,35,41", wurst)

    # Every sixth-through-tenth state in a decade is SIF-on. Duplicate and
    # malformed assignments must also fail before any input or GPU is touched.
    @test failed_plan("curry0", "6") != 0
    @test failed_plan("wurst0", "1-5,5") != 0
    @test failed_plan("curry1", "5-1") != 0
    @test failed_plan("wurst1", "not-a-state") != 0

    @test parse_state_spec("1-5,21,31-32") ==
        [1, 2, 3, 4, 5, 21, 31, 32]
    @test_throws ErrorException parse_state_spec("1-5,5")
    @test_throws ErrorException parse_state_spec("5-1")
    @test isnothing(validate_path_isolation(
        "/campaign", "/campaign/tapered.nc", "/campaign/tapered_outputs"))
    @test_throws ErrorException validate_path_isolation(
        "/campaign", "/campaign/retrieval_setup/apriori_states.nc",
        "/campaign/tapered_outputs")
    @test_throws ErrorException validate_path_isolation(
        "/campaign", "/campaign/tapered.nc", "/campaign/retrievals/round3")

    clear = (state_index=1, aerosol_case=:none)
    aerosol = (state_index=11, aerosol_case=:aod760_0p28)
    @test truth_scene_path(clear, "/campaign/truth") ==
        "/campaign/truth/hiressim_001.nc"
    @test truth_scene_path(aerosol, "/campaign/truth") ==
        "/campaign/truth/aerosol_chunked/hiressim_011.nc"
end

@testset "tapered no-SIF prior retains every non-CO2 setting" begin
    mktempdir() do directory
        reference_path = joinpath(directory, "apriori_states.nc")
        tapered_path = joinpath(
            directory,
            "apriori_states_acos_mapped_tapered_vertical_correlation.nc")
        common = (
            sif_wavelength_slope_mean=0.0,
            sif_wavelength_slope_sigma=0.002625,
            aerosol_ln_aod_sigma=0.75,
            surface_p1_sigmas=ntuple(_ -> 0.002, 3),
            surface_p2_sigmas=ntuple(_ -> 0.002, 3),
        )
        reference = Dict(
            surface => PriorBuilder.build_prior(
                surface, 0.1;
                co2_covariance_model=
                    PriorBuilder.ACOS_MAPPED_CO2_COVARIANCE_MODEL,
                common...)
            for surface in PriorBuilder.SURFACES)
        tapered = Dict(
            surface => PriorBuilder.build_prior(
                surface, 0.1;
                co2_covariance_model=PriorBuilder.TAPERED_CO2_COVARIANCE_MODEL,
                common...)
            for surface in PriorBuilder.SURFACES)
        PriorBuilder.write_netcdf(reference; output_path=reference_path)
        PriorBuilder.write_netcdf(tapered; output_path=tapered_path)

        reference_sha256 = file_sha256(reference_path)
        tapered_sha256 = file_sha256(tapered_path)

        identity = validate_tapered_prior(
            tapered_path; expected_sha256=tapered_sha256,
            reference_prior_path=reference_path,
            reference_sha256)
        @test identity.prior_sha256 == tapered_sha256
        @test identity.reference_prior_sha256 == reference_sha256
        @test_throws ErrorException validate_tapered_prior(
            tapered_path; expected_sha256=repeat("0", 64),
            reference_prior_path=reference_path,
            reference_sha256)
        @test_throws ErrorException validate_tapered_prior(
            tapered_path; expected_sha256=tapered_sha256,
            reference_prior_path=reference_path,
            reference_sha256=repeat("0", 64))

        # Corrupt a non-adjacent active CO2 covariance without touching the
        # stored adjacent-correlation diagnostic. Keep Sa and Sa_active
        # internally consistent so only the complete taper check can catch it.
        NCDataset(tapered_path, "a") do dataset
            value = Float64(dataset["Sa"][6, 8, 1]) + 1e-14
            dataset["Sa"][6, 8, 1] = value
            dataset["Sa"][8, 6, 1] = value
            dataset["Sa_active"][2, 4, 1] = value
            dataset["Sa_active"][4, 2, 1] = value
        end
        @test_throws ErrorException validate_tapered_prior(
            tapered_path; expected_sha256=file_sha256(tapered_path),
            reference_prior_path=reference_path,
            reference_sha256)

        # Restore a pristine target for the remaining adversarial checks.
        PriorBuilder.write_netcdf(tapered; output_path=tapered_path)

        NCDataset(tapered_path, "a") do dataset
            dataset.attrib["surface_p2_sigma_o2a"] = 0.001
        end
        @test_throws ErrorException validate_tapered_prior(
            tapered_path; expected_sha256=file_sha256(tapered_path),
            reference_prior_path=reference_path,
            reference_sha256)

        output_root = joinpath(directory, "isolated_output")
        expected_identity = "identity_schema=1\ncampaign_id=test-a\n"
        identity_path = initialize_campaign_identity!(
            output_root, expected_identity)
        @test read(identity_path, String) == expected_identity
        @test initialize_campaign_identity!(
            output_root, expected_identity) == identity_path
        @test_throws ErrorException initialize_campaign_identity!(
            output_root, "identity_schema=1\ncampaign_id=test-b\n")
    end
end
