#!/usr/bin/env julia

using LinearAlgebra
using NCDatasets
using Test

module Round3PriorFixture
include(joinpath(@__DIR__, "retrieval_setup", "build_apriori.jl"))
end

include(joinpath(
    @__DIR__, "retrieval_setup", "build_round4_known_sif_apriori.jl"))
include(joinpath(@__DIR__, "GattacaRound4SIFReadiness.jl"))

using .Round4KnownSIFPriorBuilder
using .GattacaRound4SIFReadiness

const LAUNCHER = joinpath(
    @__DIR__, "gattaca_round4_known_sif759_retrievals.sbatch")
const SUBMITTER = joinpath(
    @__DIR__, "submit_gattaca_round4_known_sif759_retrievals.sh")
const RUNNER = joinpath(
    @__DIR__, "run_round4_known_sif_retrievals.jl")

function write_source_prior(path)
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

function fixture_paths(root)
    repo = joinpath(root, "round4_repo")
    legacy_repo = joinpath(root, "round3_data_repo")
    private = joinpath(root, "private")
    full = joinpath(legacy_repo, "RRS_XCO2", "truth_map")
    restart = joinpath(private, "results", ".sif_v2_restart")
    bottom = joinpath(
        legacy_repo, "RRS_XCO2", "bottom_layer_XCO2_retrievals")
    stokes = joinpath(legacy_repo, "RRS_XCO2", "inversion", "instrument",
                       "representative_stokes_coefficients.nc")
    components = joinpath(bottom, "truth", "scene_components.dat")
    sif_template = joinpath(
        legacy_repo, "src", "SIF_emission", "sif-spectra.csv")
    for path in (repo, legacy_repo, private, full, restart, bottom)
        mkpath(path)
    end
    mkpath.(dirname.((stokes, components, sif_template)))
    write(stokes, "representative-stokes-fixture\n")
    write(components,
          "aod760_0p28 0.224 0.0504 0.0056 0.28 0.224 0.0504 0.0056 0.28\n")
    cp(joinpath(@__DIR__, "..", "..", "src", "SIF_emission",
                "sif-spectra.csv"), sif_template)
    run(`git -C $repo init -q`)
    run(`git -C $legacy_repo init -q`)
    return Round4ReadinessPaths(
        repo_root=repo,
        private_root=private,
        full_truth_root=full,
        restart_root=restart,
        bottom_campaign_root=bottom,
        stokes_coefficient_path=stokes,
        scene_components_path=components,
        sif_template_path=sif_template)
end

function install_prior_fixture(paths)
    mkpath(dirname(paths.source_prior_path))
    write_source_prior(paths.source_prior_path)
    prior = build_round4_prior(paths.source_prior_path, :on)
    write_round4_prior(prior; output_path=paths.prior_path)
    write_round4_summary(prior; output_path=paths.summary_path)
    return (
        prior_sha256=GattacaRound4SIFReadiness.file_sha256(paths.prior_path),
        summary_sha256=GattacaRound4SIFReadiness.file_sha256(paths.summary_path),
        source_prior_sha256=GattacaRound4SIFReadiness.file_sha256(
            paths.source_prior_path),
    )
end

function fake_barrier_keywords(; states=EXPECTED_SIF_STATES)
    full_validator = _ -> (
        validation_receipt_sha256=repeat("1", 64),)
    bottom_validator = (_, _) -> (
        receipt_sha256=repeat("2", 64),
        input_set_sha256=repeat("3", 64),)
    input_validator = (_, _, _) -> (
        states=copy(states),
        input_set_sha256=repeat("4", 64),
        table_sha256=repeat("5", 64),)
    return (; full_validator, bottom_validator, input_validator)
end

function write_output_fixture(path, paths, readiness, measurement_class)
    prior_xa, prior_Sa = NCDataset(paths.prior_path, "r") do prior
        (Float64.(Array(prior["xa"])[vcat(1, collect(6:32), 34), 1]),
         Float64.(Array(prior["Sa_active"])[:, :, 1]))
    end
    mkpath(dirname(path))
    NCDataset(path, "c") do output
        defDim(output, "state", 29)
        defDim(output, "state_2", 29)
        defDim(output, "measurement", 3)
        active = defVar(output, "active_core_parameter_index", Int16,
                        ("state",))
        active[:] = Int16.(vcat(collect(1:28), 30))
        xa = defVar(output, "a_priori_state", Float64, ("state",))
        xa[:] = prior_xa
        Sa = defVar(output, "a_priori_covariance", Float64,
                    ("state", "state_2"))
        Sa[:, :] = prior_Sa
        noise = defVar(output, "injected_measurement_noise", Float64,
                       ("measurement",))
        noise[:] = zeros(3)
        output.attrib["retrieval_complete"] = 1
        output.attrib["truth_state_index"] = 18
        output.attrib["perturbation_index"] = 11
        output.attrib["measurement_class"] = measurement_class
        output.attrib["surface"] = "urban"
        output.attrib["retrieval_state_model"] = "round4_known_sif759"
        output.attrib["round4_sif_case"] = "on"
        output.attrib["round4_known_sif_wavelength_nm"] = 759.0
        output.attrib["round4_known_sif_Lnu_mW_m-2_sr-1_per_cm-1"] =
            readiness.prior.known_Lnu
        output.attrib["state_dimension"] = 29
        output.attrib["jacobian_flavor"] =
            "OCO_RRS_synth_round4_known_sif759_boundary_chain"
        output.attrib["round4_prior_sha256"] = readiness.prior.sha256
        output.attrib["round4_code_checkpoint_sha"] = repeat("a", 40)
        output.attrib["round4_codeset_sha256"] = readiness.codeset
        output.attrib["round4_input_set_sha256"] = readiness.input_set
        output.attrib["round4_campaign_identity_sha256"] =
            readiness.identity_sha256
    end
    return path
end

@testset "Gattaca round-4 array and submission contract" begin
    @test success(`bash -n $LAUNCHER`)
    @test success(`bash -n $SUBMITTER`)
    states = Int[]
    for task in 0:39
        command = setenv(`bash $LAUNCHER`,
            "SLURM_ARRAY_TASK_ID" => string(task),
            "GATTACA_ROUND4_PLAN_ONLY" => "1")
        output = read(command, String)
        matched = match(r"state=(\d{3})", output)
        @test matched !== nothing
        push!(states, parse(Int, only(matched.captures)))
        @test occursin("sif_filter=on", output)
        @test occursin("active=29", output)
        @test occursin(CAMPAIGN_ID, output)
    end
    @test states == EXPECTED_SIF_STATES
    @test length(unique(states)) == 40
    @test states[8] == 18
    bad = setenv(`bash $LAUNCHER`,
        "SLURM_ARRAY_TASK_ID" => "40",
        "GATTACA_ROUND4_PLAN_ONLY" => "1")
    @test !success(pipeline(bad; stdout=devnull, stderr=devnull))

    launcher = read(LAUNCHER, String)
    submitter = read(SUBMITTER, String)
    runner = read(RUNNER, String)
    @test occursin("#SBATCH --array=0-39%2", launcher)
    @test occursin("run_round4_known_sif_retrievals.jl", launcher)
    @test occursin("ROUND4_SIF_MODE=on", launcher)
    @test occursin("FIRST_PERTURBATION", launcher)
    @test occursin("source \"\${identity_env}\"", launcher)
    @test occursin("set -a\nsource \"\${identity_env}\"\nset +a", launcher)
    @test occursin("symbolic-ref -q HEAD", launcher)
    @test occursin("CUDA_VISIBLE_DEVICES", launcher)
    @test occursin("state018_perturbation11_smoke.complete", launcher)
    @test occursin(
        ": \"\${RRS_REPO:?export the separate round-4 source checkout}\"",
        launcher)
    @test occursin("bottom_campaign=\"\${BOTTOM_XCO2_CAMPAIGN_ROOT}\"",
                   launcher)
    @test occursin("full_truth_root=\"\${FULL_COLUMN_TRUTH_ROOT}\"",
                   launcher)
    @test occursin("RETRIEVAL_STOKES_COEFFICIENT_PATH", launcher)
    @test occursin("RETRIEVAL_SCENE_COMPONENTS_PATH", launcher)
    @test occursin("RRS_XCO2_SIF_TEMPLATE_PATH", launcher)
    @test occursin(
        "source_apriori_states_acos_mapped_tapered_vertical_correlation.nc",
        launcher)
    @test !occursin(
        "results/bottom_layer_sif_acos_mapped_tapered_vertical_correlation_v1/retrieval_setup/apriori_states",
        launcher)
    @test !occursin(
        "bottom_campaign=\"\${repo_root}/RRS_XCO2/bottom_layer_XCO2_retrievals\"",
        launcher)
    @test occursin(
        "path_is_within \"\${repo_root}\" \"\${legacy_git_root}\"",
        launcher)
    @test occursin("--dependency=\"afterok:\${smoke_job}\"", submitter)
    @test occursin("--array=0-39%2", submitter)
    @test occursin("RRS_REPO=\${repo_root}", submitter)
    @test occursin("BOTTOM_XCO2_CAMPAIGN_ROOT=\${bottom_campaign}",
                   submitter)
    @test occursin("FULL_COLUMN_TRUTH_ROOT=\${full_truth_root}", submitter)
    @test occursin("LEGACY_INPUT_REPO=%q", submitter)
    @test occursin("[[ -f \"\${repo_root}/Manifest.toml\" ]]", submitter)
    @test occursin("component_path=scene_components_path", runner)
    @test occursin("coefficient_path=stokes_coefficient_path", runner)
    @test !occursin("run_retrievals.jl\"", launcher)
    for forbidden in ("#SBATCH --account", "#SBATCH --qos",
                      "#SBATCH --output", "#SBATCH --error")
        @test !occursin(forbidden, launcher)
    end
    @test isnothing(match(r"/home/[A-Za-z0-9_.-]+", launcher))
end

@testset "Gattaca round-4 requires two canonical checkouts" begin
    mktempdir() do root
        paths = fixture_paths(root)
        separation = validate_checkout_separation(paths)
        @test separation.code_repo == paths.repo_root
        @test separation.legacy_repo != separation.code_repo
        @test separation.expected_bottom == paths.bottom_campaign_root
        @test separation.expected_full == paths.full_truth_root

        same_checkout = Round4ReadinessPaths(
            repo_root=separation.legacy_repo,
            private_root=paths.private_root,
            full_truth_root=paths.full_truth_root,
            restart_root=paths.restart_root,
            bottom_campaign_root=paths.bottom_campaign_root,
            stokes_coefficient_path=paths.stokes_coefficient_path,
            scene_components_path=paths.scene_components_path,
            sif_template_path=paths.sif_template_path)
        @test_throws ErrorException validate_checkout_separation(same_checkout)

        wrong_full = joinpath(separation.legacy_repo, "wrong_truth")
        mkpath(wrong_full)
        mismatched = Round4ReadinessPaths(
            repo_root=paths.repo_root,
            private_root=paths.private_root,
            full_truth_root=wrong_full,
            restart_root=paths.restart_root,
            bottom_campaign_root=paths.bottom_campaign_root,
            stokes_coefficient_path=paths.stokes_coefficient_path,
            scene_components_path=paths.scene_components_path,
            sif_template_path=paths.sif_template_path)
        @test_throws ErrorException validate_checkout_separation(mismatched)
    end

    withenv("RRS_REPO" => nothing,
            "BOTTOM_XCO2_CAMPAIGN_ROOT" => nothing,
            "FULL_COLUMN_TRUTH_ROOT" => nothing,
            "RETRIEVAL_STOKES_COEFFICIENT_PATH" => nothing,
            "RETRIEVAL_SCENE_COMPONENTS_PATH" => nothing,
            "RRS_XCO2_SIF_TEMPLATE_PATH" => nothing) do
        @test_throws ErrorException Round4ReadinessPaths()
    end
end

@testset "Gattaca round-4 29-D prior identity" begin
    mktempdir() do root
        paths = fixture_paths(root)
        @test dirname(paths.source_prior_path) == dirname(paths.prior_path)
        @test basename(paths.source_prior_path) ==
            "source_apriori_states_acos_mapped_tapered_vertical_correlation.nc"
        @test !GattacaRound4SIFReadiness.contains_path(
            legacy_input_repo_root(paths), paths.source_prior_path)
        @test paths.source_prior_path !=
            GattacaRound4SIFReadiness.GattacaTaperedSIFReadiness.
                required_prior_path(paths.private_root)
        hashes = install_prior_fixture(paths)
        identity = validate_prior_identity(paths; hashes...)
        @test identity.active == vcat(1, collect(6:32), 34)
        @test identity.active_core == vcat(collect(1:28), 30)
        @test size(identity.Sa_active) == (29, 29, 4)
        @test identity.known_Lnu ==
            Round4KnownSIFPriorBuilder.Round4SIFTruthConvention.
                truth_sif_at_wavelength_nm(759.0).Lnu

        @test_throws ErrorException validate_prior_identity(paths;
            prior_sha256=repeat("0", 64),
            summary_sha256=hashes.summary_sha256,
            source_prior_sha256=hashes.source_prior_sha256)

        NCDataset(paths.prior_path, "a") do prior
            prior.attrib["known_sif_wavelength_nm"] = 760.0
        end
        @test_throws ErrorException validate_prior_identity(paths;
            prior_sha256=GattacaRound4SIFReadiness.file_sha256(
                paths.prior_path),
            summary_sha256=hashes.summary_sha256,
            source_prior_sha256=hashes.source_prior_sha256)
    end
end

@testset "Gattaca round-4 publication barrier is mandatory" begin
    mktempdir() do root
        paths = fixture_paths(root)
        calls = Symbol[]
        expected_legacy = legacy_input_repo_root(paths)
        full = legacy -> begin
            @test legacy.repo_root == expected_legacy
            @test legacy.full_truth_root == paths.full_truth_root
            push!(calls, :full)
            (validation_receipt_sha256=repeat("1", 64),)
        end
        bottom = (legacy, _) -> begin
            @test legacy.repo_root == expected_legacy
            @test legacy.bottom_campaign_root == paths.bottom_campaign_root
            push!(calls, :bottom)
            (receipt_sha256=repeat("2", 64),
             input_set_sha256=repeat("3", 64))
        end
        inputs = (legacy, _, _) -> begin
            @test legacy.repo_root == expected_legacy
            push!(calls, :inputs)
            (states=copy(EXPECTED_SIF_STATES),
             input_set_sha256=repeat("4", 64),
             table_sha256=repeat("5", 64))
        end
        release = validate_release_barrier(paths;
            full_validator=full, bottom_validator=bottom,
            input_validator=inputs)
        @test calls == [:full, :bottom, :inputs]
        @test release.inputs.states == EXPECTED_SIF_STATES

        broken = fake_barrier_keywords(states=EXPECTED_SIF_STATES[1:39])
        @test_throws ErrorException validate_release_barrier(paths; broken...)
        @test_throws ErrorException validate_release_barrier(paths;
            full_validator=_ -> error("publication incomplete"),
            bottom_validator=bottom, input_validator=inputs)
    end
end

@testset "Gattaca round-4 identity and output isolation" begin
    mktempdir() do root
        # Use the real checkout only for the read-only code-set digest.
        paths0 = fixture_paths(root)
        paths = Round4ReadinessPaths(
            repo_root=normpath(joinpath(@__DIR__, "..", "..")),
            private_root=paths0.private_root,
            full_truth_root=paths0.full_truth_root,
            restart_root=paths0.restart_root,
            bottom_campaign_root=paths0.bottom_campaign_root,
            stokes_coefficient_path=paths0.stokes_coefficient_path,
            scene_components_path=paths0.scene_components_path,
            sif_template_path=paths0.sif_template_path)
        hashes = install_prior_fixture(paths)
        barrier = fake_barrier_keywords()
        readiness = prepare_readiness!(paths;
            hashes...,
            checkpoint=repeat("a", 40), barrier...)
        @test isfile(readiness.identity)
        @test isfile(readiness.environment)
        identity_text = read(readiness.identity, String)
        @test occursin(
            "representative_stokes_coefficients_sha256=" *
            GattacaRound4SIFReadiness.file_sha256(
                paths.stokes_coefficient_path), identity_text)
        @test occursin(
            "bottom_scene_components_sha256=" *
            GattacaRound4SIFReadiness.file_sha256(
                paths.scene_components_path), identity_text)
        @test occursin(
            "sif_spectral_template_sha256=" *
            GattacaRound4SIFReadiness.file_sha256(
                paths.sif_template_path), identity_text)
        @test readiness == prepare_readiness!(paths;
            hashes..., checkpoint=repeat("a", 40), barrier...)
        original_stokes = read(paths.stokes_coefficient_path)
        write(paths.stokes_coefficient_path,
              vcat(original_stokes, codeunits("changed\n")))
        @test_throws ErrorException prepare_readiness!(paths;
            hashes..., checkpoint=repeat("a", 40), barrier...)
        write(paths.stokes_coefficient_path, original_stokes)
        @test validate_checkout_separation(paths).legacy_repo ==
            GattacaRound4SIFReadiness.legacy_input_repo_root(paths)
        @test validate_output_isolation(paths) == paths.output_root

        for measurement_class in ("corrected", "uncorrected")
            output = joinpath(paths.output_root, measurement_class,
                "retrieval_state018_perturbation11.nc")
            write_output_fixture(
                output, paths, readiness, measurement_class)
        end
        @test validate_state_outputs(paths, 18; perturbations=[11]) == 2
        corrupted = joinpath(paths.output_root, "corrected",
            "retrieval_state018_perturbation11.nc")
        NCDataset(corrupted, "a") do output
            output["a_priori_state"][1] += 1
        end
        @test_throws ErrorException validate_state_outputs(
            paths, 18; perturbations=[11])

        wrong = Round4ReadinessPaths(
            paths.repo_root, paths.private_root, paths.full_truth_root,
            paths.restart_root, paths.bottom_campaign_root,
            paths.bottom_truth_table, paths.measurement_directory,
            paths.noise_directory, paths.source_prior_path,
            paths.prior_path, paths.summary_path,
            GattacaRound4SIFReadiness.GattacaTaperedSIFReadiness.
                required_output_root(paths.private_root),
            paths.stokes_coefficient_path, paths.scene_components_path,
            paths.sif_template_path)
        @test_throws ErrorException validate_output_isolation(wrong)
    end
end
