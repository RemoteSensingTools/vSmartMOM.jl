#!/usr/bin/env julia

using Test

include(joinpath(@__DIR__, "RetrievalCases.jl"))
using .RetrievalCases

const TEST_SURFACES = (:urban, :rural, :desert, :forest)

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
                (sif_index, sif) in enumerate((:off, :total_0p5)),
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
