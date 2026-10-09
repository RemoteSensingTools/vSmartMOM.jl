#!/usr/bin/env julia

# Run from test/: julia --project=/path/to/uni_vSmartMOM_round5_jacobian \
#   ../RRS_XCO2/inversion/test_round7_retrieval_campaign.jl
# No GPU or truth-observation generation is performed. Synthetic observations
# are confined to mktempdir; the released v2 table and Round-5 priors are read-only.
using Test
using LinearAlgebra

const TEST_INPUTS = get(ENV, "ROUND7_TEST_INPUT_DIR", normpath(joinpath(@__DIR__,
    "..", "bottom_layer_XCO2_retrievals", "round6_fixed_sif", "deployment",
    "round6_fixed_sif_standard_utls_v1", "inputs")))
const TEST_PRIORS = get(ENV, "ROUND7_TEST_PRIOR_DIRECTORY", normpath(joinpath(@__DIR__,
    "..", "bottom_layer_XCO2_retrievals", "round5_fixed_sif", "retrieval_setup")))
ENV["ROUND7_PREFLIGHT_ONLY"] = "1"
ENV["RRS_XCO2_SIF_TEMPLATE_PATH"] = joinpath(TEST_INPUTS, "sif-spectra.csv")
include(joinpath(@__DIR__, "run_round7_imperfect_correction_retrievals.jl"))
using .Round4SIFTruthConvention

const TEST_TRUTH = joinpath(TEST_INPUTS, "true_states_corrected_sif_v2.dat")
const TEST_SOURCE_PRIOR = joinpath(TEST_INPUTS, "source_prior.nc")
const TEST_TRUTHS = read_truth_cases(TEST_TRUTH)
const TEST_TRUTH_HASH = file_sha256(TEST_TRUTH)
const TEST_OFF = R5.resolve_round5_sif_map([TEST_TRUTHS[1]])
const TEST_ON = R5.resolve_round5_sif_map([TEST_TRUTHS[6]])
test_prior_path(mode) = joinpath(TEST_PRIORS,
    "apriori_states_round5_fixed_sif_$(mode)_tight_utls_acos_mapped_tapered_vertical_correlation.nc")

function fixture_observation(path, truth; reversed_dimensions=false)
    n = 2742
    wavelength = vcat(collect(758.0:0.015:772.0), collect(1594.0:0.031:1619.0),
                      collect(2042.0:0.04:2082.0))
    @assert length(wavelength) == n
    uncorrected = 0.2 .+ collect(1:n) .* 1e-6
    correction = vcat(fill(0.0001, 934), zeros(n - 934))
    imperfect = uncorrected - correction
    noise_std = fill(0.002, n)
    draw = [sin(i + j) for i in 1:n, j in 1:11]
    draw[:, 11] .= 0
    injected = noise_std .* draw
    NCDataset(path, "c") do ds
        for (name, length) in (("measurement", n), ("perturbation", 11),
                               ("band", 3), ("surface_coefficient", 3))
            defDim(ds, name, length)
        end
        for (key, value) in Round7RetrievalCampaign.DEFINITION
            ds.attrib[key] = value
        end
        ds.attrib["round7_observation_complete"] = 1
        for key in (:state_index, :surface, :aerosol_case, :sif_case)
            value = getproperty(truth, key)
            ds.attrib[string(key)] = value isa Symbol ? string(value) : value
        end
        ds.attrib["source_truth_table_sha256"] = TEST_TRUTH_HASH
        ds.attrib["round5_prior_sha256"] = file_sha256(test_prior_path(truth.sif_case == :off ? "off" : "on"))
        ds.attrib["source_analyzer_sha256"] = file_sha256(joinpath(TEST_INPUTS, "representative_stokes_coefficients.nc"))
        ds.attrib["source_generator_sha256"] = repeat("a", 64)
        ds.attrib["source_grid_sha256"] = repeat("b", 64)
        for (key, value) in ("psurf_hpa" => 1000.0, "sza_deg" => 30.0,
                             "vza_deg" => 0.0, "relative_azimuth_deg" => 0.0)
            ds.attrib[key] = value
        end
        if truth.sif_case != :off
            for (key, value) in expected_round4_sif_provenance(; enabled=true)
                ds.attrib[key] = value
            end
        end
        for (name, values) in (
                "wavelength" => wavelength, "measurement_uncorrected" => uncorrected,
                "measurement_ideal_corrected" => uncorrected .- 2correction,
                "measurement_imperfectly_corrected" => imperfect,
                "correction_observation" => correction,
                "noise_standard_deviation" => noise_std,
                "Se_diagonal" => noise_std .^ 2)
            defVar(ds, name, Float64, ("measurement",); deflatelevel=1)[:] = values
        end
        for (name, values, dimension) in (
                ("band_start_index", [1, 935, 1742], "band"),
                ("band_end_index", [934, 1741, 2742], "band"),
                ("perturbation_index", collect(1:11), "perturbation"))
            defVar(ds, name, Int32, (dimension,))[:] = values
        end
        defVar(ds, "random_seed_uint64", UInt64, ("perturbation",))[:] =
            vcat(UInt64(2)^60 .+ UInt64.(1:10), UInt64(0))
        defVar(ds, "mean_o2a_surface_coefficients", Float64,
               ("surface_coefficient",))[:] = [0.2 + truth.state_index * 1e-4, 0.001, -0.002]
        dims = reversed_dimensions ? ("perturbation", "measurement") : ("measurement", "perturbation")
        for (name, values) in (
                "uncorrected_perturbed" => uncorrected .+ injected,
                "imperfectly_corrected_perturbed" => imperfect .+ injected,
                "normalized_noise_draw" => draw, "injected_measurement_noise" => injected)
            defVar(ds, name, Float64, dims; deflatelevel=1)[:, :] =
                reversed_dimensions ? permutedims(values) : values
        end
    end
    return path
end

function fixture_manifest(directory; count=80)
    path = joinpath(directory, "SHA256SUMS")
    isfile(joinpath(directory, "provenance.json")) || write(joinpath(directory, "provenance.json"), "{}\n")
    open(path, "w") do io
        for index in 1:count
            name = @sprintf("OCO2round7_%03d.nc", index)
            println(io, file_sha256(joinpath(directory, name)), "  ", name)
        end
        println(io, file_sha256(joinpath(directory, "provenance.json")), "  provenance.json")
    end
    return path
end

function fixture_load(path, truth=TEST_TRUTHS[1], map=TEST_OFF)
    # A singleton handle is used only in isolated schema tests. Production uses
    # validate_round7_manifest, which always checks exactly all 80 files.
    manifest = (; directory=dirname(path), entries=Dict(basename(path) => file_sha256(path)))
    load_round7_observation(manifest, truth, map; truth_table_sha256=TEST_TRUTH_HASH,
        prior_sha256=file_sha256(test_prior_path(map.sif_on ? "on" : "off")),
        analyzer_sha256=file_sha256(joinpath(TEST_INPUTS, "representative_stokes_coefficients.nc")))
end

@testset "pinned Round-5 map/prior and CPU-only imports" begin
    @test !isdefined(Main, :VSmartMOMForward)
    @test !isdefined(Main, :CUDA)
    runtime = round7_runtime_provenance()
    @test runtime["round7_vsmartmom_source_path"] == realpath(joinpath(BASE_PROJECT_DIRECTORY, "src", "vSmartMOM.jl"))
    @test runtime["round7_project_sha256"] == file_sha256(Base.active_project())
    @test ROUND7_BASE_CHECKPOINT == strip(read(`git -C $BASE_PROJECT_DIRECTORY rev-parse HEAD`, String))
    @test active_core_indices(TEST_OFF) == active_core_indices(TEST_ON) == collect(1:28)
    @test R5.expected_active_to_full(TEST_OFF) == vcat(1, collect(6:32))
    @test all(iszero, expand_round5_state(TEST_OFF, zeros(28))[29:30])
    for (mode, map) in (("off", TEST_OFF), ("on", TEST_ON)), surface in (:urban, :rural, :desert, :forest)
        path = test_prior_path(mode)
        prior = load_retrieval_prior(surface; path)
        @test validate_round7_prior(prior, map; path,
            source_prior_path=TEST_SOURCE_PRIOR) === prior
    end
    mktempdir() do directory
        path = joinpath(directory, "prior.nc")
        cp(test_prior_path("off"), path)
        NCDataset(path, "a") do ds
            ds["xa"][33, 1] = 0.1
        end
        prior = load_retrieval_prior(:urban; path)
        @test_throws ErrorException validate_round7_prior(prior, TEST_OFF; path,
            source_prior_path=TEST_SOURCE_PRIOR)
    end
end

@testset "checksummed 80-scene manifest and portable identity" begin
    mktempdir() do directory
        for index in 1:80
            write(joinpath(directory, @sprintf("OCO2round7_%03d.nc", index)), string(index))
        end
        path = fixture_manifest(directory)
        digest = file_sha256(path)
        manifest = validate_round7_manifest(directory; expected_sha256=digest)
        @test length(manifest.entries) == 81
        @test_throws ErrorException validate_round7_manifest(directory; expected_sha256=repeat("0", 64))
        identities = round7_identity_values(; checkpoint=ROUND7_BASE_CHECKPOINT,
            base_tree_sha256=repeat("a", 64), code_paths=[path], input_paths=[path], observation_manifest=path)
        @test round7_campaign_identity(identities; expected=identities)["round7_campaign_identity_status"] == "complete"
        @test_throws ErrorException round7_campaign_identity(Dict())
        for key in keys(identities)
            changed = copy(identities)
            changed[key] = repeat("0", length(changed[key]))
            @test_throws ErrorException round7_campaign_identity(changed; expected=identities)
        end
        alias = copy(identities)
        alias["OBSERVATION_MANIFEST_SHA256"] = pop!(alias, "ROUND7_OBSERVATION_MANIFEST_SHA256")
        @test round7_campaign_identity(alias; expected=identities)["round7_observation_manifest_sha256"] == digest
        # Tampering even with an unselected scene invalidates the whole set.
        write(joinpath(directory, "OCO2round7_080.nc"), "tampered")
        @test_throws ErrorException validate_round7_manifest(directory; expected_sha256=digest)
        fixture_manifest(directory; count=79)
        @test_throws ErrorException validate_round7_manifest(directory; expected_sha256=file_sha256(path))
        fixture_manifest(directory)
        open(path, "a") do io
            println(io, repeat("0", 64), "  ../outside.nc")
        end
        @test_throws ErrorException validate_round7_manifest(directory; expected_sha256=file_sha256(path))
        fixture_manifest(directory)
        open(path, "a") do io
            println(io, file_sha256(joinpath(directory, "OCO2round7_001.nc")), "  OCO2round7_001.nc")
        end
        @test_throws ErrorException validate_round7_manifest(directory; expected_sha256=file_sha256(path))
    end
end

@testset "observation schema, copied SIF, and exact stored noise" begin
    mktempdir() do directory
        path = joinpath(directory, "OCO2round7_001.nc")
        fixture_observation(path, TEST_TRUTHS[1])
        observation = fixture_load(path)
        @test size(observation.perturbed) == (2742, 11)
        @test observation.seeds[1] == UInt64(2)^60 + 1
        @test observation.seeds[11] == 0
        for index in 1:11
            realization = round7_realization(observation, index)
            @test realization.perturbed == realization.noiseless + realization.noise_std .* realization.normalized_draw
            @test realization.perturbed == observation.perturbed[:, index]
            @test realization.variance == observation.variance
        end
        @test round7_realization(observation, 11).perturbed == observation.imperfect
        @test_throws ArgumentError round7_realization(observation, 0)
        @test_throws ArgumentError round7_realization(observation, 12)

        mutations = [
            ds -> (ds.attrib["round7_observation_complete"] = 0),
            ds -> (ds.attrib["round7_definition_version"] = 2),
            ds -> (ds.attrib["round7_correction_operation"] = "add"),
            ds -> (ds.attrib["round7_pressure_geometry_source"] = "retrieved"),
            ds -> (ds.attrib["round7_noise_policy"] = "redraw"),
            ds -> (ds.attrib["source_truth_table_sha256"] = repeat("0", 64)),
            ds -> (ds.attrib["round5_prior_sha256"] = repeat("0", 64)),
            ds -> (ds.attrib["source_analyzer_sha256"] = repeat("0", 64)),
            ds -> (ds.attrib["psurf_hpa"] = 999.0),
            ds -> (ds.attrib["sza_deg"] = 31.0),
            ds -> (ds.attrib["vza_deg"] = 1.0),
            ds -> (ds.attrib["relative_azimuth_deg"] = 1.0),
            ds -> (ds.attrib["state_index"] = 2),
            ds -> (ds.attrib["surface"] = "forest"),
            ds -> (ds.attrib["aerosol_case"] = "unexpected"),
            ds -> (ds["band_start_index"][2] = 936),
            ds -> (ds["perturbation_index"][1] = 0),
            ds -> (ds["measurement_ideal_corrected"][1] = NaN),
            ds -> (ds["noise_standard_deviation"][1] = -0.002),
            ds -> (ds["Se_diagonal"][1] = 0),
            ds -> (ds["correction_observation"][1] = 0.1),
            ds -> (ds["correction_observation"][935] = eps(Float64)),
            ds -> (ds["random_seed_uint64"][11] = UInt64(1)),
            ds -> (ds["injected_measurement_noise"][1, 1] = 1),
            ds -> (ds["uncorrected_perturbed"][1, 1] = 1),
            ds -> (ds["imperfectly_corrected_perturbed"][1, 1] = 1),
            ds -> (ds["normalized_noise_draw"][1, 11] = 0.01),
            # Within arithmetic tolerance, but not exact for the base writer.
            ds -> (ds["imperfectly_corrected_perturbed"][1, 1] = nextfloat(ds["imperfectly_corrected_perturbed"][1, 1]))]
        for mutate in mutations
            fixture_observation(path, TEST_TRUTHS[1])
            NCDataset(mutate, path, "a")
            @test_throws ErrorException fixture_load(path)
        end
        fixture_observation(path, TEST_TRUTHS[1]; reversed_dimensions=true)
        @test_throws ErrorException fixture_load(path)
        on_path = joinpath(directory, "OCO2round7_006.nc")
        fixture_observation(on_path, TEST_TRUTHS[6])
        @test fixture_load(on_path, TEST_TRUTHS[6], TEST_ON).provenance["sif_definition_version"] == 2
        NCDataset(on_path, "a") do ds
            ds.attrib["sif_definition_version"] = 1
        end
        @test_throws ErrorException fixture_load(on_path, TEST_TRUTHS[6], TEST_ON)
    end
end

@testset "unchanged writer, strict smoke and immutable resume" begin
    mktempdir() do directory
        path = fixture_observation(joinpath(directory, "OCO2round7_001.nc"), TEST_TRUTHS[1])
        observation = fixture_load(path)
        experiments = round7_experiments(TEST_TRUTHS, directory)
        @test length(experiments) == 880
        @test all(e -> e.measurement_class == :corrected, experiments)
        @test first(experiments).noise_index == 11
        experiment = first(experiments)
        realization = round7_realization(observation, 11)
        prior_path = test_prior_path("off")
        prior = load_retrieval_prior(:urban; path=prior_path)
        provenance = round7_output_provenance(TEST_OFF; prior_path, observation,
            campaign_identity=Dict("round7_campaign_identity_sha256" => repeat("b", 64)),
            truth_table_sha256=TEST_TRUTH_HASH)
        @test provenance["retrieval_state_model"] == ROUND7_STATE_MODEL
        @test provenance["round7_underlying_state_model"] == "round5_fixed_sif"
        @test provenance["round5_stratospheric_aerosol_sigma_scale"] == 0.1
        @test provenance["round7_observation_kind"] == "imperfectly_corrected"
        jacobian = zeros(2742, 28)
        for i in 1:28
            jacobian[i, i] = 1e-5
        end
        evaluate(x) = ForwardEvaluation(realization.noiseless + jacobian * (x - prior.xa),
                                        jacobian, realization.band_ranges)
        result = solve_optimal_estimation(evaluate, realization.perturbed,
            realization.variance, prior.xa, prior.Sa)
        @test result.converged && result.fit_quality_ok
        output = joinpath(directory, "retrieval.nc")
        @test !validate_round7_existing_output(output, experiment, realization, prior, TEST_OFF; provenance)
        write_retrieval_result(experiment, realization, result, prior.xa, prior.Sa,
            prior.parameter_names; output_path=output, provenance,
            jacobian_flavor=R5.ROUND5_JACOBIAN_FLAVOR, active_to_core=collect(1:28))
        @test validate_round7_existing_output(output, experiment, realization, prior, TEST_OFF; provenance, smoke=true)
        original_hash = file_sha256(output)
        link = joinpath(directory, "published.nc")
        hardlink(output, link)
        @test file_sha256(link) == original_hash
        @test_throws Base.IOError hardlink(output, link)
        bad = copy(provenance)
        bad["round7_source_observation_sha256"] = repeat("0", 64)
        @test_throws ErrorException validate_round7_existing_output(output, experiment,
            realization, prior, TEST_OFF; provenance=bad)
        @test file_sha256(output) == original_hash
        @test_throws ArgumentError write_retrieval_result(experiment, realization, result,
            prior.xa, prior.Sa, prior.parameter_names; output_path=output)
        NCDataset(output, "a") do ds
            ds.attrib["fit_quality_ok"] = 0
        end
        @test validate_round7_existing_output(output, experiment, realization, prior, TEST_OFF; provenance)
        @test_throws ErrorException validate_round7_existing_output(output, experiment,
            realization, prior, TEST_OFF; provenance, smoke=true)
        NCDataset(output, "a") do ds
            ds.attrib["retrieval_complete"] = 0
        end
        @test_throws ErrorException validate_round7_existing_output(output, experiment,
            realization, prior, TEST_OFF; provenance)
        @test_throws ErrorException round7_smoke_gate(true, false)
        @test_throws ErrorException round7_smoke_gate(false, true)
        @test round7_smoke_gate(true, true)
        @test endswith(validate_round7_output_root(joinpath(directory, "round7", "retrievals_nosif");
            sif_on=false), "retrievals_nosif")
        @test_throws ArgumentError validate_round7_output_root(joinpath(directory, "round5", "retrievals_nosif"); sif_on=false)
        @test_throws ArgumentError validate_round7_output_root(joinpath(directory, "round7", "retrievals_sif"); sif_on=false)
    end
end

@testset "full no-GPU preflight across both SIF modes and filters" begin
    mktempdir() do directory
        observations = joinpath(directory, "observations")
        mkdir(observations)
        for truth in TEST_TRUTHS
            fixture_observation(joinpath(observations, @sprintf("OCO2round7_%03d.nc", truth.state_index)), truth)
        end
        fixture_manifest(observations)
        base_environment = Dict(
            "ROUND7_OBSERVATION_DIR" => observations,
            "ROUND7_SOURCE_PRIOR" => TEST_SOURCE_PRIOR,
            "RETRIEVAL_TRUTH_TABLE" => TEST_TRUTH,
            "RETRIEVAL_STOKES_COEFFICIENT_PATH" => joinpath(TEST_INPUTS, "representative_stokes_coefficients.nc"),
            "RETRIEVAL_SCENE_COMPONENTS_PATH" => joinpath(TEST_INPUTS, "scene_components.dat"),
            "RRS_XCO2_SIF_TEMPLATE_PATH" => joinpath(TEST_INPUTS, "sif-spectra.csv"),
            # Preflight hashes these fixtures; it never executes RT or loads solar data.
            "SOLAR_OUT" => joinpath(TEST_INPUTS, "sif-spectra.csv"),
            "RRS_XCO2_CONFIG" => joinpath(BASE_PROJECT_DIRECTORY, "sandbox", "workflows", "RRS_XCO2", "config", "oco_grass_3aerosol.yaml"),
            "ROUND7_PREFLIGHT_ONLY" => "1", "ROUND7_PRINT_IDENTITY" => "0",
            "ROUND7_SMOKE" => "0", "RETRIEVAL_CLASS" => "corrected", "FORCE" => "0",
            "RETRIEVAL_ARCH" => "GPU", "RETRIEVAL_FLOAT_TYPE" => "Float32",
            "RETRIEVAL_NSTREAMS" => "9", "FIRST_STATE" => "1", "LAST_STATE" => "80",
            "FIRST_PERTURBATION" => "1", "LAST_PERTURBATION" => "11")
        withenv(base_environment...) do
            for mode in ("off", "on"), aerosol in ("none", "aerosol")
                output = joinpath(directory, "round7", mode == "on" ? "retrievals_sif" : "retrievals_nosif")
                withenv("SIF_CASE_FILTER" => mode, "AEROSOL_CASE_FILTER" => aerosol,
                        "RETRIEVAL_OUTPUT_ROOT" => output, "RETRIEVAL_PRIOR_PATH" => test_prior_path(mode)) do
                    identity = checked_round7_identity(round7_inputs())
                    withenv(identity...) do
                        campaign = main()
                        @test length(campaign.experiments) == 220
                        @test length(campaign.observations) == 20
                        @test !ispath(output)
                        @test !isdefined(Main, :VSmartMOMForward)
                        @test campaign.experiments[2].random_seed == UInt64(2)^60 + 1
                        withenv("FORCE" => "1") do
                            @test_throws ErrorException prepare_round7()
                        end
                        withenv("ROUND7_SMOKE" => "1") do
                            @test_throws ErrorException prepare_round7()
                        end
                        withenv("RETRIEVAL_CLASS" => "uncorrected") do
                            @test_throws ErrorException prepare_round7()
                        end
                    end
                end
            end
        end
    end
end
