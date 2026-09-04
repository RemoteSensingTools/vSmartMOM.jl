#!/usr/bin/env julia

using Test
using NCDatasets
using Printf
using SHA

include(joinpath(@__DIR__, "RetrievalCases.jl"))
using .RetrievalCases

const TEST_SURFACES = (:urban, :rural, :desert, :forest)

function corrected_sif_provenance()
    radiance = 0.5 / (2π)
    return Dict{String,Any}(
        "sif_definition_version" => Int32(2),
        "sif_definition" =>
            "isotropic BOA radiance normalized by 2pi*L_lambda(760 nm)=0.5",
        "sif_case_on_label" => "angular_integral760_0p5",
        "sif_reference_wavelength_nm" => 760.0,
        "sif_upwelling_solid_angle_sr" => 2π,
        "sif_angular_integral_760_mW_m-2_nm-1" => 0.5,
        "sif_radiance_760_mW_m-2_sr-1_nm-1" => radiance,
        "sif_cosine_weighted_irradiance_760_mW_m-2_nm-1" => π * radiance,
        "sif_SIF760_mW_m-2_sr-1_per_cm-1" =>
            radiance * 760.0^2 / 1e7,
        "sif_mSIF_mW_m-2_sr-1_per_cm-2" => 1.2291230681458325e-5,
        "sif_template_wavelength_integral_mW_m-2_sr-1" =>
            15.368806005166872,
    )
end

function write_test_measurement(path, truth; provenance=Dict{String,Any}())
    mkpath(dirname(path))
    NCDataset(path, "c") do dataset
        dataset.attrib["instrument_processing_complete"] = 1
        dataset.attrib["state_index"] = truth.state_index
        dataset.attrib["sif_case"] = String(truth.sif_case)
        for (key, value) in provenance
            dataset.attrib[key] = value
        end
    end
    return path
end

function write_test_noise(path, truth; provenance=Dict{String,Any}())
    mkpath(dirname(path))
    NCDataset(path, "c") do dataset
        defDim(dataset, "measurement", 3)
        defDim(dataset, "band", 3)
        for class in ("corrected", "uncorrected")
            defVar(dataset, "measurement_$class", Float64,
                   ("measurement",))[:] = [1.0, 2.0, 3.0]
            defVar(dataset, "noise_std_$class", Float64,
                   ("measurement",))[:] = fill(0.1, 3)
            defVar(dataset, "Se_diagonal_$class", Float64,
                   ("measurement",))[:] = fill(0.01, 3)
        end
        defVar(dataset, "wavelength", Float64,
               ("measurement",))[:] = [760.0, 1600.0, 2050.0]
        defVar(dataset, "band_start_index", Int32,
               ("band",))[:] = Int32[1, 2, 3]
        defVar(dataset, "band_end_index", Int32,
               ("band",))[:] = Int32[1, 2, 3]
        dataset.attrib["noise_covariance_complete"] = 1
        dataset.attrib["state_index"] = truth.state_index
        dataset.attrib["sif_case"] = String(truth.sif_case)
        dataset.attrib["campaign"] = String(truth.campaign)
        dataset.attrib["background_co2_ppm"] = truth.background_co2_ppm
        dataset.attrib["bottom_co2_layer_index"] = truth.bottom_layer_index
        dataset.attrib["bottom_co2_ppm"] = truth.bottom_co2_ppm
        dataset.attrib["xco2_ppm"] = truth.xco2_ppm
        for (key, value) in provenance
            dataset.attrib[key] = value
        end
    end
    return path
end

function test_experiment(truth, measurement_path, noise_path)
    return RetrievalExperiment(
        1, 1, truth, UNPERTURBED_INDEX, :corrected, UInt64(0),
        measurement_path, noise_path)
end

test_sha256(path) = open(path, "r") do io
    bytes2hex(sha256(io))
end

function write_bottom_release_fixture(directory, table, cases)
    measurement_directory = joinpath(directory, "OCO_radiances")
    noise_directory = joinpath(measurement_directory, "noise_covariances")
    rows = String[]
    for truth in filter(case -> case.sif_case != :off, cases)
        truth_directory = truth.aerosol_case == :none ? directory :
            joinpath(directory, "aerosol_chunked")
        truth_path = joinpath(
            truth_directory, @sprintf("hiressim_%03d.nc", truth.state_index))
        measurement_path = joinpath(
            measurement_directory, @sprintf("OCO2sims_%03d.nc", truth.state_index))
        noise_path = joinpath(
            noise_directory, @sprintf("OCO2noise_%03d.nc", truth.state_index))
        for (path, label) in ((truth_path, "truth"),
                              (measurement_path, "measurement"),
                              (noise_path, "noise"))
            mkpath(dirname(path))
            write(path, "$label $(truth.state_index)\n")
        end
        push!(rows, @sprintf("%03d %s %s %s", truth.state_index,
                            test_sha256(truth_path),
                            test_sha256(measurement_path),
                            test_sha256(noise_path)))
    end
    full_truth_root = joinpath(directory, "full_truth")
    mkpath(full_truth_root)
    full_receipt = joinpath(full_truth_root, "sif_v2_release_complete.dat")
    write(full_receipt, "synthetic full-column release receipt\n")
    receipt = joinpath(directory, "bottom_layer_sif_v2_release_complete.dat")
    open(receipt, "w") do io
        println(io, "# release_schema 1")
        println(io, "# sif_definition_version 2")
        println(io, "# full_column_release_receipt_sha256 ",
                test_sha256(full_receipt))
        println(io, "# bottom_state_table_sha256 ", test_sha256(table))
        println(io, "# legacy_bottom_state_table_sha256 ", repeat("b", 64))
        println(io, "# legacy_bottom_state_table_archive_relative truth/true_states.dat")
        println(io, "# legacy_archive_manifest_sha256 ", repeat("c", 64))
        println(io, "# no_sif_byte_preservation_policy preserve_truth_measurement_noise_bytes_and_legacy_state_table_hash_attributes")
        println(io, "# no_sif_triplet_set_sha256 ", repeat("d", 64))
        println(io, "# input_set_sha256 ", repeat("a", 64))
        println(io, "# state truth_sha256 measurement_sha256 noise_sha256")
        foreach(line -> println(io, line), rows)
    end
    return (; receipt, measurement_directory, noise_directory,
            full_truth_root)
end

function write_test_truth_table(path; bottom_layer=false)
    mkpath(dirname(path))
    open(path, "w") do io
        if bottom_layer
            println(io, "# index surface_index surface aerosol_index aerosol_case " *
                        "sif_index sif_case bottom_co2_index background_co2_ppm " *
                        "bottom_layer_index bottom_co2_ppm xco2_ppm")
        else
            println(io, "# index surface_index surface aerosol_index aerosol_case " *
                        "sif_index sif_case xco2_index xco2_ppm")
        end
        index = 0
        co2_values = bottom_layer ? (360, 380, 400, 420, 440) :
                                    (380, 400, 420, 440)
        for (surface_index, surface) in enumerate(TEST_SURFACES),
                (aerosol_index, aerosol) in enumerate((:none, :aod760_0p28)),
                (sif_index, sif) in enumerate((
                    :off, :angular_integral760_0p5)),
                (co2_index, co2) in enumerate(co2_values)
            index += 1
            if bottom_layer
                xco2 = 400 + 0.06229780735737046 * (co2 - 400)
                println(io, "$index $surface_index $surface $aerosol_index " *
                            "$aerosol $sif_index $sif $co2_index 400 16 " *
                            "$co2 $xco2")
            else
                println(io, "$index $surface_index $surface $aerosol_index " *
                            "$aerosol $sif_index $sif $co2_index $co2")
            end
        end
    end
    return path
end

@testset "versioned SIF provenance gates retrieval inputs" begin
    mktempdir() do directory
        table = write_test_truth_table(
            joinpath(directory, "bottom.dat"); bottom_layer=true)
        truths = read_truth_cases(table)
        sif_on = truths[6]
        no_sif = truths[1]

        measurement = joinpath(directory, "OCO2sims_006.nc")
        noise = joinpath(directory, "OCO2noise_006.nc")
        provenance = corrected_sif_provenance()
        write_test_measurement(measurement, sif_on; provenance)
        write_test_noise(noise, sif_on; provenance)
        realization = load_measurement_realization(
            test_experiment(sif_on, measurement, noise))
        @test realization.provenance["sif_definition_version"] == 2
        @test realization.provenance[
            "sif_angular_integral_760_mW_m-2_nm-1"] == 0.5

        stale = copy(provenance)
        stale["sif_definition_version"] = Int32(1)
        write_test_measurement(measurement, sif_on; provenance=stale)
        @test_throws ErrorException load_measurement_realization(
            test_experiment(sif_on, measurement, noise))

        mismatched = copy(provenance)
        mismatched["sif_mSIF_mW_m-2_sr-1_per_cm-2"] *= 2
        write_test_measurement(measurement, sif_on; provenance)
        write_test_noise(noise, sif_on; provenance=mismatched)
        @test_throws ErrorException load_measurement_realization(
            test_experiment(sif_on, measurement, noise))

        # No-SIF products do not need the metadata introduced by definition 2.
        clear_measurement = joinpath(directory, "OCO2sims_001.nc")
        clear_noise = joinpath(directory, "OCO2noise_001.nc")
        write_test_measurement(clear_measurement, no_sif)
        write_test_noise(clear_noise, no_sif)
        clear = load_measurement_realization(
            test_experiment(no_sif, clear_measurement, clear_noise))
        @test isempty(clear.provenance)
    end
end

@testset "full-column and bottom-layer truth case schemas" begin
    mktempdir() do directory
        full_path = write_test_truth_table(joinpath(directory, "full.dat"))
        bottom_path = write_test_truth_table(
            joinpath(directory, "bottom.dat"); bottom_layer=true)

        full = read_truth_cases(full_path)
        @test length(full) == 64
        @test length(read_no_sif_truth_cases(full_path)) == 32
        @test all(case -> case.campaign == :full_column_XCO2, full)
        @test all(case -> case.co2_profile_mode == :uniform_column, full)
        @test full[1].fixed_upper_co2_ppm == 380
        @test full[1].background_co2_ppm == 380
        @test full[1].bottom_layer_index == 0
        @test full[1].bottom_co2_ppm == 380

        bottom = read_truth_cases(bottom_path)
        @test length(bottom) == 80
        @test length(read_no_sif_truth_cases(bottom_path)) == 40
        @test all(case -> case.campaign == :bottom_layer_XCO2, bottom)
        @test all(case -> case.co2_profile_mode == :bottom_layer, bottom)
        @test all(case -> case.fixed_upper_co2_ppm == 400, bottom)
        @test all(case -> case.background_co2_ppm == 400, bottom)
        @test all(case -> case.bottom_layer_index == 16, bottom)
        @test bottom[1].xco2_index == 1
        @test bottom[1].bottom_co2_ppm == 360
        @test bottom[1].xco2_ppm ≈ 397.508087705705

        # The legacy positional constructor remains a uniform-column shortcut.
        legacy = TruthCase(1, 1, :urban, 1, :none, :off, 1, 380)
        @test legacy.campaign == :full_column_XCO2
        @test legacy.co2_profile_mode == :uniform_column
        @test legacy.fixed_upper_co2_ppm == 380
        @test legacy.bottom_layer_index == 0
    end
end

@testset "durable external ownership of SIF-on retrievals" begin
    mktempdir() do directory
        table = write_test_truth_table(
            joinpath(directory, "bottom.dat"); bottom_layer=true)
        truth = read_truth_cases(table)
        no_sif = build_experiments(
            [truth[1]]; measurement_directory=directory,
            noise_directory=directory, validate_inputs=false)
        with_sif = build_experiments(
            [truth[6]]; measurement_directory=directory,
            noise_directory=directory, validate_inputs=false)
        output_root = joinpath(directory, "retrievals")
        marker = external_sif_ownership_marker(output_root)

        @test endswith(marker, joinpath(".control", "sif_owned_externally"))
        @test isnothing(enforce_sif_ownership(output_root, no_sif))
        @test isnothing(enforce_sif_ownership(output_root, with_sif))

        mkpath(dirname(marker))
        write(marker, "owner=external-worker\n")
        @test isnothing(enforce_sif_ownership(output_root, no_sif))
        error = try
            enforce_sif_ownership(output_root, with_sif)
            nothing
        catch exception
            exception
        end
        @test error isa ErrorException
        @test occursin("owner=external-worker", sprint(showerror, error))

        alternate = joinpath(directory, "alternate.marker")
        withenv("RETRIEVAL_EXTERNAL_SIF_OWNERSHIP_MARKER" => alternate) do
            @test external_sif_ownership_marker(output_root) == alternate
            @test isnothing(enforce_sif_ownership(output_root, with_sif))
        end
    end
end

@testset "SIF retrievals require an intact all-scene release receipt" begin
    mktempdir() do directory
        table = write_test_truth_table(
            joinpath(directory, "true_states.dat"); bottom_layer=true)
        cases = read_truth_cases(table)
        no_sif = select_sif_truth_cases(cases, :off)
        with_sif = select_sif_truth_cases(cases, :on)
        measurement_directory = joinpath(directory, "OCO_radiances")
        noise_directory = joinpath(measurement_directory, "noise_covariances")

        # No-SIF work remains independent of the corrected-SIF publication.
        @test isnothing(require_sif_release_barrier(
            table, no_sif, measurement_directory, noise_directory))
        @test_throws ErrorException require_sif_release_barrier(
            table, with_sif, measurement_directory, noise_directory)

        fixture = write_bottom_release_fixture(directory, table, cases)
        withenv("FULL_COLUMN_TRUTH_ROOT" => fixture.full_truth_root) do
            @test require_sif_release_barrier(
                table, with_sif, fixture.measurement_directory,
                fixture.noise_directory) == fixture.receipt
        end

        marker = joinpath(
            directory, ".bottom_layer_sif_v2_publication_in_progress")
        write(marker, "test interrupted publication\n")
        withenv("FULL_COLUMN_TRUTH_ROOT" => fixture.full_truth_root) do
            @test_throws ErrorException require_sif_release_barrier(
                table, with_sif, fixture.measurement_directory,
                fixture.noise_directory)
        end
        rm(marker)

        tampered = joinpath(
            fixture.measurement_directory, "OCO2sims_006.nc")
        write(tampered, "tampered\n")
        withenv("FULL_COLUMN_TRUTH_ROOT" => fixture.full_truth_root) do
            @test_throws ErrorException require_sif_release_barrier(
                table, with_sif, fixture.measurement_directory,
                fixture.noise_directory)
        end
    end
end

@testset "campaign-local experiment inputs and manifest outputs" begin
    mktempdir() do directory
        table = write_test_truth_table(
            joinpath(directory, "bottom.dat"); bottom_layer=true)
        truth = first(read_no_sif_truth_cases(table))
        measurement_directory = joinpath(directory, "measurements")
        noise_directory = joinpath(directory, "noise")
        experiments = build_experiments(
            [truth]; measurement_directory, noise_directory,
            validate_inputs=false)
        @test length(experiments) == 22
        @test all(experiment -> startswith(
            experiment.measurement_path, measurement_directory), experiments)
        @test all(experiment -> startswith(
            experiment.noise_path, noise_directory), experiments)

        output_root = joinpath(directory, "retrievals")
        manifest = write_experiment_manifest(
            experiments; inversion_root=output_root)
        @test manifest == joinpath(output_root, "retrieval_manifest.dat")
        text = read(manifest, String)
        @test occursin("bottom_layer_XCO2", text)
        @test occursin("bottom_layer", text)
        @test occursin(joinpath(output_root, "corrected",
                                "retrieval_state001_perturbation11.nc"), text)
        @test occursin(measurement_directory, text)
        @test occursin(noise_directory, text)
    end
end
