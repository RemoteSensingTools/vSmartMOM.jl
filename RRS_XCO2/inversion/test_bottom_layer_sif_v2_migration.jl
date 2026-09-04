#!/usr/bin/env julia

using NCDatasets
using Printf
using SHA
using Test

include(joinpath(@__DIR__, "..", "scripts",
                 "migrate_bottom_layer_sif_v2.jl"))
using .BottomLayerSIFV2Migration
const Migration = BottomLayerSIFV2Migration

file_hash(path) = open(path, "r") do io
    bytes2hex(sha256(io))
end

function bottom_fields(index)
    surface_index = fld(index - 1, 20) + 1
    local_index = mod(index - 1, 20) + 1
    aerosol_index = fld(local_index - 1, 10) + 1
    within_aerosol = mod(local_index - 1, 10) + 1
    sif_index = fld(within_aerosol - 1, 5) + 1
    bottom_co2_index = mod(within_aerosol - 1, 5) + 1
    surfaces = ("urban", "rural", "desert", "forest")
    aerosols = ("none", "aod760_0p28")
    bottom_values = (360.0, 380.0, 400.0, 420.0, 440.0)
    bottom = bottom_values[bottom_co2_index]
    return (; index, surface_index, surface=surfaces[surface_index],
            aerosol_index, aerosol=aerosols[aerosol_index], sif_index,
            bottom_co2_index, bottom,
            xco2=Migration.BottomLayerTruthCommon.column_xco2_ppm(bottom))
end

function write_bottom_table(path; corrected)
    mkpath(dirname(path))
    sif = Migration.RRSXCO2Common.campaign_sif_state()
    open(path, "w") do io
        println(io, "# synthetic 80-state bottom-layer migration fixture")
        println(io, "# index surface_index surface aerosol_index aerosol_case " *
                    "sif_index sif_case bottom_co2_index background_co2_ppm " *
                    "bottom_layer_index bottom_co2_ppm xco2_ppm psurf_hpa " *
                    "sza_deg vza_deg relative_azimuth_deg sulfate_aod550 " *
                    "organic_aod550 stratospheric_aod550 " *
                    "sif_angular_integral760 SIF760 mSIF o2a_P0 o2a_P1 " *
                    "o2a_P2 weak_P0 weak_P1 weak_P2 strong_P0 strong_P1 " *
                    "strong_P2")
        for index in 1:80
            state = bottom_fields(index)
            on = state.sif_index == 2
            label = on ? (corrected ?
                Migration.RRSXCO2Common.SIF_CASE_ON : "total_0p5") : "off"
            angular = on ? 0.5 : 0.0
            sif760 = on ? (corrected ? sif.SIF760 : 0.5) : 0.0
            msif = on ? (corrected ? sif.mSIF : 0.01) : 0.0
            aod = state.aerosol_index == 2 ? (0.20, 0.07, 0.01) :
                                                   (0.0, 0.0, 0.0)
            values = Any[
                index, state.surface_index, state.surface,
                state.aerosol_index, state.aerosol, state.sif_index, label,
                state.bottom_co2_index, 400.0, 16, state.bottom, state.xco2,
                1000.0, 30.0, 0.0, 0.0, aod..., angular, sif760, msif,
                0.2, 0.001, 0.0001,
                0.3, 0.002, 0.0002,
                0.4, 0.003, 0.0003,
            ]
            println(io, join(values, ' '))
        end
    end
    return path
end

function truth_arrays(seed)
    o2 = [Float32(seed + 0.01r + 0.001j) for r in 1:3, j in 1:4]
    weak = [Float32(seed + 0.02r + 0.002j) for r in 1:3, j in 1:3]
    strong = [Float32(seed + 0.03r + 0.003j) for r in 1:3, j in 1:2]
    return o2, o2 .+ 0.1f0, o2 .+ 0.2f0, weak, strong
end

function write_truth(path, seed; corrected_sif=false, state_index=0)
    mkpath(dirname(path))
    values = truth_arrays(seed)
    NCDataset(path, "c") do dataset
        defDim(dataset, "stokes", 3)
        defDim(dataset, "o2a", 4)
        defDim(dataset, "weak", 3)
        defDim(dataset, "strong", 2)
        dimensions = ("o2a", "o2a", "o2a", "weak", "strong")
        for (name, value, dimension) in zip(
                Migration.ALL_TRUTH_VARIABLES, values, dimensions)
            defVar(dataset, name, Float32,
                   ("stokes", dimension))[:, :] = value
        end
        dataset.attrib["simulation_complete"] = Int32(1)
        dataset.attrib["state_index"] = Int32(state_index)
        dataset.attrib["co2_absco_completed"] = "source-$seed"
        dataset.attrib["co2_absco_regeneration_complete"] = Int32(1)
        dataset.attrib["strong_co2_short_shoulder_merged"] = Int32(1)
        dataset.attrib["strong_co2_short_shoulder_points"] = Int32(seed)
        dataset.attrib["strong_co2_shoulder_source"] = "source-$seed.nc"
        if corrected_sif
            dataset.attrib["o2_truth_regenerated"] = Int32(1)
            dataset.attrib["o2_truth_reused"] = Int32(0)
            dataset.attrib["sif_case"] =
                Migration.RRSXCO2Common.SIF_CASE_ON
            Migration.RRSXCO2Common.write_sif_provenance!(
                dataset.attrib, true)
        end
    end
    return path
end

function write_full_release!(paths, states)
    full_table = joinpath(paths.full_truth_root, "true_states.dat")
    mkpath(dirname(full_table))
    write(full_table, "synthetic corrected full-column state table\n")
    for aerosol_index in 1:2
        directory = Migration.scene_root(
            paths.full_truth_root, aerosol_index)
        mkpath(directory)
        write(joinpath(directory, "sim_wavelength.nc"),
              "shared synthetic wavelength grid\n")
    end

    source_hashes = Dict{Int,String}()
    for state in states
        old = Migration.BottomLayerTruthCommon.old_control_index(state)
        haskey(source_hashes, old) && continue
        source = Migration.full_source_scene(paths, state)
        write_truth(source, old; corrected_sif=true, state_index=old)
        source_hashes[old] = file_hash(source)
    end

    validation = paths.full_validation_receipt
    mkpath(dirname(validation))
    open(validation, "w") do io
        println(io, "# release_schema 1")
        println(io, "# sif_definition_version 2")
        println(io, "# corrected_state_table_sha256 ", file_hash(full_table))
        println(io, "# index class staged_sha256 archived_sha256 canonical_destination")
        for index in Migration.EXPECTED_FULL_COLUMN_SIF_STATES
            digest = get(source_hashes, index, repeat("0", 64))
            @printf(io, "%03d clear %s %s synthetic\n",
                    index, digest, repeat("1", 64))
        end
    end
    receipt = joinpath(paths.full_truth_root,
                       Migration.FULL_RELEASE_RECEIPT_NAME)
    open(receipt, "w") do io
        println(io, "release_schema 1")
        println(io, "sif_definition_version 2")
        println(io, "validation_receipt_sha256 ", file_hash(validation))
        println(io, "corrected_state_table_sha256 ", file_hash(full_table))
    end
    return source_hashes
end

function write_measurement(path, paths, state)
    mkpath(dirname(path))
    NCDataset(path, "c") do dataset
        for (band, n) in (("o2a", 4), ("weak_co2", 3),
                          ("strong_co2", 2))
            defDim(dataset, band, n)
            rayleigh = collect(range(1.0, 1.0 + 0.01(n - 1); length=n))
            corrected = copy(rayleigh)
            if band == "o2a"
                cabannes = rayleigh .+ 0.2
                rrs = fill(0.03, n)
                uncorrected = cabannes .+ rrs
                defVar(dataset, "I_OCO_cabannes_o2a", Float64,
                       (band,))[:] = cabannes
                defVar(dataset, "I_OCO_rrs_o2a", Float64,
                       (band,))[:] = rrs
            else
                uncorrected = copy(rayleigh)
            end
            defVar(dataset, "I_OCO_rayleigh_$band", Float64,
                   (band,))[:] = rayleigh
            defVar(dataset, "I_OCO_corrected_$band", Float64,
                   (band,))[:] = corrected
            defVar(dataset, "I_OCO_uncorrected_$band", Float64,
                   (band,))[:] = uncorrected
        end
        dataset.attrib["instrument_processing_complete"] = Int32(1)
        dataset.attrib["state_index"] = Int32(state.index)
        dataset.attrib["source_truth_scene"] = abspath(
            Migration.truth_scene(Migration.truth_root(paths), state))
        dataset.attrib["source_truth_sha256"] = file_hash(
            Migration.truth_scene(Migration.stage_truth_root(paths), state))
        dataset.attrib["instrument_processor_script"] =
            Migration.INSTRUMENT_PROCESSOR
        dataset.attrib["instrument_processor_script_sha256"] =
            file_hash(Migration.INSTRUMENT_PROCESSOR)
        dataset.attrib["representative_stokes_coefficients"] =
            paths.stokes_coefficients
        dataset.attrib["sif_case"] = state.sif_case
        Migration.RRSXCO2Common.write_sif_provenance!(dataset.attrib, true)
    end
end

function measurement_vector(path, class)
    NCDataset(path) do dataset
        return vcat((Array(dataset["I_OCO_$(class)_$band"][:])
                     for band in ("o2a", "weak_co2", "strong_co2"))...)
    end
end

function write_noise(path, paths, state, measurement_path)
    mkpath(dirname(path))
    corrected = measurement_vector(measurement_path, "corrected")
    uncorrected = measurement_vector(measurement_path, "uncorrected")
    NCDataset(path, "c") do dataset
        defDim(dataset, "measurement", length(corrected))
        for (class, measurement) in (("corrected", corrected),
                                     ("uncorrected", uncorrected))
            sigma = fill(0.1, length(measurement))
            variance = sigma .^ 2
            defVar(dataset, "measurement_$class", Float64,
                   ("measurement",))[:] = measurement
            defVar(dataset, "noise_std_$class", Float64,
                   ("measurement",))[:] = sigma
            defVar(dataset, "Se_diagonal_$class", Float64,
                   ("measurement",))[:] = variance
        end
        dataset.attrib["noise_covariance_complete"] = Int32(1)
        dataset.attrib["state_index"] = Int32(state.index)
        dataset.attrib["source_synthetic_measurement"] = abspath(
            Migration.oco_scene(Migration.oco_root(paths), state.index))
        dataset.attrib["source_synthetic_measurement_sha256"] =
            file_hash(measurement_path)
        dataset.attrib["noise_processor_script"] = Migration.NOISE_PROCESSOR
        dataset.attrib["noise_processor_script_sha256"] =
            file_hash(Migration.NOISE_PROCESSOR)
        dataset.attrib["representative_snr_coefficients"] =
            paths.snr_coefficients
        dataset.attrib["sif_case"] = state.sif_case
        Migration.RRSXCO2Common.write_sif_provenance!(dataset.attrib, true)
    end
end

function fixture(root)
    campaign = joinpath(root, "bottom")
    full = joinpath(root, "full")
    restart = joinpath(root, "stage")
    archive = joinpath(root, "archive")
    stokes = joinpath(root, "stokes.nc")
    snr = joinpath(root, "snr.nc")
    write(stokes, "synthetic stokes coefficients\n")
    write(snr, "synthetic snr coefficients\n")
    paths = MigrationPaths(;
        campaign_root=campaign, full_truth_root=full,
        restart_root=restart, archive_root=archive,
        full_validation_receipt=joinpath(
            full, ".sif_v2_restart", "sif_v2_release_validation.dat"),
        stokes_coefficients=stokes, snr_coefficients=snr)

    canonical_truth = Migration.truth_root(paths)
    staged_truth = Migration.stage_truth_root(paths)
    write_bottom_table(joinpath(canonical_truth, "true_states.dat");
                       corrected=false)
    write_bottom_table(joinpath(staged_truth, "true_states.dat");
                       corrected=true)
    for name in ("control_reuse_map.dat", "scene_components.dat")
        mkpath(canonical_truth)
        write(joinpath(canonical_truth, name), "legacy $name\n")
        mkpath(staged_truth)
        write(joinpath(staged_truth, name), "corrected $name\n")
    end
    states = Migration.BottomLayerTruthCommon.read_bottom_states(
        joinpath(staged_truth, "true_states.dat"))
    write_full_release!(paths,
        filter(state -> state.sif_index == 2, states))

    for aerosol_index in 1:2
        directory = Migration.scene_root(canonical_truth, aerosol_index)
        mkpath(directory)
        write(joinpath(directory, "sim_wavelength.nc"),
              "shared synthetic wavelength grid\n")
    end

    for state in states
        canonical_scene = Migration.truth_scene(canonical_truth, state)
        canonical_measurement = Migration.oco_scene(
            Migration.oco_root(paths), state.index)
        canonical_noise = Migration.noise_scene(
            Migration.noise_root(paths), state.index)
        if state.sif_index == 1
            write_truth(canonical_scene, state.index;
                        state_index=state.index)
            mkpath(dirname(canonical_measurement))
            write(canonical_measurement,
                  "legacy no-SIF measurement $(state.index)\n")
            mkpath(dirname(canonical_noise))
            write(canonical_noise, "legacy no-SIF noise $(state.index)\n")
        else
            mkpath(dirname(canonical_scene))
            write(canonical_scene, "legacy SIF truth $(state.index)\n")
            mkpath(dirname(canonical_measurement))
            write(canonical_measurement,
                  "legacy SIF measurement $(state.index)\n")
            mkpath(dirname(canonical_noise))
            write(canonical_noise, "legacy SIF noise $(state.index)\n")
        end
    end
    write(joinpath(Migration.noise_root(paths),
                   "noise_covariance_manifest.dat"), "legacy manifest\n")

    Migration.copy_no_sif_dependencies!(paths, states)
    table_hash = file_hash(joinpath(staged_truth, "true_states.dat"))
    for state in filter(state -> state.sif_index == 2, states)
        Migration.assemble_truth_scene!(paths, state, table_hash)
        measurement = Migration.oco_scene(
            Migration.stage_oco_root(paths), state.index)
        write_measurement(measurement, paths, state)
        write_noise(Migration.noise_scene(
            Migration.stage_noise_root(paths), state.index),
            paths, state, measurement)
    end
    write(joinpath(Migration.stage_noise_root(paths),
                   "noise_covariance_manifest.dat"), "corrected manifest\n")

    control = joinpath(Migration.retrieval_root(paths), ".control")
    mkpath(control)
    write(joinpath(control, "sif_owned_externally"),
          "owner=synthetic-external-compute\n")
    sif_retrieval = joinpath(Migration.retrieval_root(paths), "corrected",
                             "retrieval_state006_perturbation11.nc")
    no_sif_retrieval = joinpath(Migration.retrieval_root(paths), "corrected",
                                "retrieval_state001_perturbation11.nc")
    mkpath(dirname(sif_retrieval))
    write(sif_retrieval, "legacy SIF retrieval\n")
    write(no_sif_retrieval, "accepted no-SIF retrieval\n")
    write(joinpath(Migration.retrieval_root(paths), "retrieval_manifest.dat"),
          "legacy mixed retrieval manifest\n")
    return (; paths, states, sif_retrieval, no_sif_retrieval)
end

@testset "bottom-layer SIF-v2 validation and publication transaction" begin
    mktempdir() do root
        built = fixture(root)
        paths = built.paths
        validation = validate_release(paths)
        @test length(validation.triplets) == 40
        @test isfile(validation.receipt)
        @test !isfile(Migration.publication_marker(paths))

        # A same-label stale convolution/noise product cannot be laundered by
        # updating only its source pathname.  Stage production always uses
        # FORCE=1; these source-byte bindings independently make stale bytes
        # fail validation before publication.
        state = first(filter(state -> state.sif_index == 2,
                             validation.states))
        measurement = Migration.oco_scene(
            Migration.stage_oco_root(paths), state.index)
        truth_digest = file_hash(Migration.truth_scene(
            Migration.stage_truth_root(paths), state))
        NCDataset(measurement, "a") do dataset
            dataset.attrib["source_truth_sha256"] = repeat("0", 64)
        end
        @test_throws ErrorException validate_release(
            paths; write_staging_receipt=false)
        NCDataset(measurement, "a") do dataset
            dataset.attrib["source_truth_sha256"] = truth_digest
        end
        noise = Migration.noise_scene(
            Migration.stage_noise_root(paths), state.index)
        measurement_digest = file_hash(measurement)
        NCDataset(noise, "a") do dataset
            dataset.attrib["source_synthetic_measurement_sha256"] =
                repeat("0", 64)
        end
        @test_throws ErrorException validate_release(
            paths; write_staging_receipt=false)
        NCDataset(noise, "a") do dataset
            dataset.attrib["source_synthetic_measurement_sha256"] =
                measurement_digest
        end
        validation = validate_release(paths)

        pairs = Migration.replacement_pairs(paths, validation.states)
        legacy_hashes = Dict(destination => file_hash(destination)
                             for (_, destination) in pairs)
        no_sif_hash = file_hash(built.no_sif_retrieval)
        sif_retrieval_hash = file_hash(built.sif_retrieval)

        withenv("CONFIRM_BOTTOM_SIF_V2_PUBLICATION" =>
                Migration.CONFIRM_PUBLICATION) do
            mkdir(Migration.publication_lock(paths))
            @test_throws ErrorException publish_release(paths)
            withenv("BREAK_STALE_BOTTOM_SIF_V2_LOCK" => "1") do
                @test_throws ErrorException publish_release(
                    paths; _test_fail_after_promotions=7)
            end
        end
        @test isfile(Migration.publication_marker(paths))
        @test all(file_hash(path) == digest
                  for (path, digest) in legacy_hashes)
        @test file_hash(built.sif_retrieval) == sif_retrieval_hash
        @test file_hash(built.no_sif_retrieval) == no_sif_hash

        marker_values, _ = Migration.read_keyed_receipt(
            Migration.publication_marker(paths))
        release_id = marker_values["release_id"]
        corrected_sidecar = last(first(pairs)) *
            ".sif-v2-pending.$release_id"
        retrieval_sidecar = built.sif_retrieval *
            ".sif-v1-pending.$release_id"
        write(corrected_sidecar, "hard-kill corrected sidecar\n")
        write(retrieval_sidecar, "hard-kill retrieval sidecar\n")
        withenv("CONFIRM_BOTTOM_SIF_V2_RESTORE" =>
                Migration.CONFIRM_RESTORE) do
            restore_legacy(paths)
        end
        @test !isfile(Migration.publication_marker(paths))
        @test !isfile(corrected_sidecar)
        @test !isfile(retrieval_sidecar)
        @test all(file_hash(path) == digest
                  for (path, digest) in legacy_hashes)

        withenv("CONFIRM_BOTTOM_SIF_V2_PUBLICATION" =>
                Migration.CONFIRM_PUBLICATION) do
            publish_release(paths)
        end
        @test !isfile(Migration.publication_marker(paths))
        @test isfile(Migration.canonical_receipt(paths))
        @test file_hash(joinpath(Migration.truth_root(paths),
                                 "true_states.dat")) ==
              file_hash(joinpath(Migration.stage_truth_root(paths),
                                 "true_states.dat"))
        @test !isfile(built.sif_retrieval)
        @test file_hash(built.no_sif_retrieval) == no_sif_hash
        @test isfile(Migration.archive_relative(paths, built.sif_retrieval))

        receipt_values, receipt_rows = Migration.read_keyed_receipt(
            Migration.canonical_receipt(paths))
        @test receipt_values["release_schema"] == "1"
        @test receipt_values["sif_definition_version"] == "2"
        @test receipt_values["legacy_bottom_state_table_sha256"] ==
              legacy_hashes[joinpath(Migration.truth_root(paths),
                                     "true_states.dat")]
        @test receipt_values["legacy_archive_manifest_sha256"] ==
              file_hash(Migration.archive_manifest(paths))
        @test receipt_values["no_sif_byte_preservation_policy"] ==
              "preserve_truth_measurement_noise_bytes_and_legacy_state_table_hash_attributes"
        @test length(receipt_rows) == 40

        cases = Migration.RetrievalCases.read_truth_cases(
            joinpath(Migration.truth_root(paths), "true_states.dat"))
        selected = Migration.RetrievalCases.select_sif_truth_cases(cases, :on)
        withenv("FULL_COLUMN_TRUTH_ROOT" => paths.full_truth_root) do
            @test Migration.RetrievalCases.require_sif_release_barrier(
                joinpath(Migration.truth_root(paths), "true_states.dat"),
                selected, Migration.oco_root(paths),
                Migration.noise_root(paths)) ==
                  Migration.canonical_receipt(paths)
        end
    end
end
