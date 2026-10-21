using Test
using vSmartMOM: sif_data_path

const WORKFLOW_TEST_ROOT = normpath(joinpath(@__DIR__, "..", "..", "..",
    "sandbox", "workflows", "RRS_XCO2", "inversion"))

module WorkflowFixtures
    using NCDatasets, Printf
    include(joinpath(@__DIR__, "..", "..", "..", "sandbox", "workflows",
                     "RRS_XCO2", "inversion", "RetrievalCases.jl"))
    using .RetrievalCases
    include(joinpath(@__DIR__, "..", "..", "..", "sandbox", "workflows",
                     "RRS_XCO2", "inversion", "test_support", "synthetic_campaign.jl"))
end

mktempdir() do fixture_root
    table = WorkflowFixtures.write_test_truth_table(
        joinpath(fixture_root, "truth_map", "true_states.dat"))
    for truth in WorkflowFixtures.read_truth_cases(table)
        number = lpad(truth.state_index, 3, '0')
        measurement_dir = joinpath(fixture_root, "truth_map", "OCO_radiances")
        provenance = truth.sif_case == :off ? Dict{String,Any}() :
                     WorkflowFixtures.corrected_sif_provenance()
        WorkflowFixtures.write_test_measurement(
            joinpath(measurement_dir, "OCO2sims_$number.nc"), truth; provenance)
        WorkflowFixtures.write_test_noise(
            joinpath(measurement_dir, "noise_covariances", "OCO2noise_$number.nc"),
            truth; provenance)
    end
    albedo = joinpath(fixture_root, "lambertian_legendre_inputs.dat")
    write(albedo, "# Synthetic surface coefficients for checksum contract tests\n")
    manifest = joinpath(fixture_root, "Manifest.toml")
    write(manifest, "# Synthetic producer manifest for provenance contract tests\n")
    withenv("RRS_XCO2_DATA_ROOT" => fixture_root,
            "SIF_PRODUCER_MANIFEST" => manifest, "SIF_PRODUCER_ALBEDO" => albedo) do
        covariance = get(ENV, "CO2_COVARIANCE_FILE", "")
        has_covariance = isfile(covariance)
        if has_covariance
            prior_module = Module(gensym(:WorkflowPrior))
            Base.include(prior_module, joinpath(WORKFLOW_TEST_ROOT,
                "retrieval_setup", "build_apriori.jl"))
            Core.eval(prior_module, quote
                mkpath(OUTPUT_ROOT)
                priors = Dict(surface => build_prior(surface, 0.1) for surface in SURFACES)
                write_netcdf(priors)
            end)
        end
        @testset "Portable retrieval workflow contracts" begin
            for path in sort(readdir(WORKFLOW_TEST_ROOT; join=true))
                name = basename(path)
                startswith(name, "test_") && endswith(name, ".jl") || continue
                # This separate integration test constructs the real OCO
                # evaluator and requires external spectroscopy/solar/prior data.
                name == "test_forward_state_mapping.jl" && continue
                if !isfile(sif_data_path("sif-spectra.csv")) && name ∉ (
                    "test_optimal_estimation.jl", "test_retrieval_cases_campaigns.jl")
                    @test_skip "External SIF template required: $name"
                    continue
                end
                if !has_covariance && name in (
                    "test_apriori_aerosol_sigma.jl", "test_co2_prior_taper.jl",
                    "test_round4_known_sif_apriori.jl", "test_tapered_no_sif_launcher.jl")
                    @test_skip "External CO2 covariance required: $name"
                    continue
                end
                println("WORKFLOW_TEST ", name); flush(stdout)
                @testset "$name" begin
                    mod = Module(gensym(:WorkflowTests))
                    Core.eval(mod, :(include(path) = Base.include($mod, path)))
                    Base.include(mod, path)
                end
            end
        end
    end
end
