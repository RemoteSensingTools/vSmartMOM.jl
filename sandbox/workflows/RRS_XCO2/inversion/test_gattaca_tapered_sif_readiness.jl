#!/usr/bin/env julia

using NCDatasets
using LinearAlgebra
using Test

include(joinpath(@__DIR__, "GattacaTaperedSIFReadiness.jl"))
using .GattacaTaperedSIFReadiness

const LAUNCHER = joinpath(
    @__DIR__, "gattaca_tapered_sif_retrievals.sbatch")
const STATIC_PACKAGER = normpath(joinpath(
    @__DIR__, "..", "scripts",
    "package_gattaca_sif_release_dependencies.sh"))
const SOURCE_SHA = repeat("a", 40)

function test_paths(root; output_root=nothing, prior_path=nothing)
    repo = joinpath(root, "repo")
    private = joinpath(root, "private")
    full = joinpath(repo, "RRS_XCO2", "truth_map")
    restart = joinpath(private, "results", ".sif_v2_restart")
    bottom = joinpath(repo, "RRS_XCO2", "bottom_layer_XCO2_retrievals")
    for path in (repo, private, full, restart, bottom)
        mkpath(path)
    end
    table = joinpath(bottom, "truth", "true_states.dat")
    measurement = joinpath(bottom, "truth", "OCO_radiances")
    noise = joinpath(measurement, "noise_covariances")
    prior = something(prior_path, required_prior_path(private))
    output = something(output_root, required_output_root(private))
    return ReadinessPaths(
        realpath(repo), realpath(private), realpath(full), realpath(restart),
        realpath(bottom), abspath(table), abspath(measurement), abspath(noise),
        abspath(prior), abspath(output))
end

function synthetic_co2_correlation(adjacent)
    result = Matrix{Float64}(I, 12, 12)
    for row in 1:11
        value = 1.0
        for column in row + 1:12
            value *= adjacent[column - 1]
            result[row, column] = result[column, row] = value
        end
    end
    return result
end

function write_prior(path; model=EXPECTED_PRIOR_MODEL,
                     base_model=EXPECTED_PRIOR_BASE_MODEL,
                     alter_nonco2=false, alter_xa=false,
                     alter_nonadjacent_co2=false,
                     alter_co2_diagonal=false,
                     alter_active_mapping=false,
                     alter_prior_sigma=false,
                     aerosol_ln_aod_sigma=0.75,
                     surface_p1_sigmas=(0.002, 0.002, 0.002),
                     surface_p2_sigmas=(0.002, 0.002, 0.002),
                     sif_slope_mean=0.0,
                     sif_slope_sigma=0.002625)
    mkpath(dirname(path))
    active = vcat(1, collect(6:34))
    alter_active_mapping && (active[2] = 5)
    xa_values = repeat(reshape(collect(1.0:34.0), 34, 1), 1, 4)
    alter_xa && (xa_values[18, 1] += 1)
    covariance = zeros(Float64, 34, 34, 4)
    base_adjacent = fill(0.8, 11)
    retention = model == EXPECTED_PRIOR_MODEL ?
        copy(GattacaTaperedSIFReadiness.EXPECTED_TAPER_RETENTION) : ones(11)
    selected_adjacent = base_adjacent .* retention
    base_correlation = synthetic_co2_correlation(base_adjacent)
    selected_correlation = synthetic_co2_correlation(selected_adjacent)
    co2_sigma = collect(range(10e-6, 21e-6; length=12))
    co2_covariance = Diagonal(co2_sigma) * selected_correlation *
        Diagonal(co2_sigma)
    for surface in 1:4
        covariance[1, 1, surface] = 2500.0
        covariance[6:17, 6:17, surface] .= co2_covariance
        for index in 18:34
            covariance[index, index, surface] = 0.01 * index
        end
        covariance[33, 34, surface] = covariance[34, 33, surface] = 0.001
    end
    alter_nonco2 && (covariance[18, 18, 1] += 1)
    alter_co2_diagonal && (covariance[6, 6, 1] += 1e-12)
    if alter_nonadjacent_co2
        covariance[6, 8, 1] += 1e-12
        covariance[8, 6, 1] = covariance[6, 8, 1]
    end
    NCDataset(path, "c") do dataset
        defDim(dataset, "parameter", 34)
        defDim(dataset, "parameter_2", 34)
        defDim(dataset, "active_parameter", 30)
        defDim(dataset, "active_parameter_2", 30)
        defDim(dataset, "surface", 4)
        defDim(dataset, "sif_wavelength_parameter", 2)
        defDim(dataset, "co2_adjacent_layer_pair", 11)
        active_variable = defVar(
            dataset, "active_parameter_index", Int32, ("active_parameter",))
        active_variable[:] = Int32.(active)
        xa = defVar(dataset, "xa", Float64, ("parameter", "surface"))
        xa[:, :] = xa_values
        Sa = defVar(dataset, "Sa", Float64,
                    ("parameter", "parameter_2", "surface"))
        Sa[:, :, :] = covariance
        Sa_active = defVar(dataset, "Sa_active", Float64,
            ("active_parameter", "active_parameter_2", "surface"))
        Sa_active[:, :, :] = covariance[active, active, :]
        sigma = defVar(
            dataset, "prior_sigma", Float64, ("parameter", "surface"))
        prior_sigma = sqrt.(max.(
            [covariance[index, index, surface]
             for index in 1:34, surface in 1:4], 0.0))
        alter_prior_sigma && (prior_sigma[18, 1] += 1e-6)
        sigma[:, :] = prior_sigma
        fixed = defVar(dataset, "active_mask", Int8, ("parameter",))
        fixed[:] = Int8.([index in active for index in 1:34])
        centers = defVar(
            dataset, "co2_layer_center_height", Float64, ("parameter",))
        centers[:] = vcat(NaN, collect(15.5:-1.0:0.5), fill(NaN, 17))
        pressures = defVar(
            dataset, "co2_layer_center_pressure", Float64, ("parameter",))
        pressures[:] = vcat(NaN, collect(100.0:50.0:850.0), fill(NaN, 17))
        sif_state = defVar(dataset, "sif_wavelength_state", Float64,
                           ("sif_wavelength_parameter", "surface"))
        sif_state[:, :] = repeat([0.1, sif_slope_mean], 1, 4)
        sif_sigma = defVar(dataset, "sif_wavelength_sigma", Float64,
                           ("sif_wavelength_parameter", "surface"))
        sif_sigma[:, :] = repeat([0.25, sif_slope_sigma], 1, 4)
        taper = defVar(dataset,
            "co2_correlation_taper_adjacent_retention", Float64,
            ("co2_adjacent_layer_pair",))
        taper[:] = retention
        base = defVar(dataset, "co2_base_adjacent_correlation", Float64,
                      ("co2_adjacent_layer_pair",))
        base[:] = base_adjacent
        selected = defVar(dataset, "co2_selected_adjacent_correlation", Float64,
                          ("co2_adjacent_layer_pair",))
        selected[:] = selected_adjacent
        dataset.attrib["apriori_complete"] = 1
        dataset.attrib["co2_covariance_model"] = model
        dataset.attrib["co2_covariance_base_model"] = base_model
        dataset.attrib["co2_correlation_taper_type"] =
            model == EXPECTED_PRIOR_MODEL ? "nonstationary_ar1_schur" : "none"
        dataset.attrib["aerosol_ln_aod_sigma"] = aerosol_ln_aod_sigma
        for (index, band) in enumerate(("o2a", "weak_co2", "strong_co2"))
            dataset.attrib["surface_p1_sigma_$band"] = surface_p1_sigmas[index]
            dataset.attrib["surface_p2_sigma_$band"] = surface_p2_sigmas[index]
        end
        dataset.attrib["sif_slope_reference_radiance_mw_m2_sr_nm"] = 0.1
        dataset.attrib["sif_wavelength_slope_prior_mw_m2_sr_nm2"] =
            sif_slope_mean
        dataset.attrib["sif_wavelength_slope_sigma_mw_m2_sr_nm2"] =
            sif_slope_sigma
        dataset.attrib["surface_order"] = "urban rural desert forest"
        dataset.attrib["parameter_names"] = join(
            ["parameter_$index" for index in 1:34], ' ')
        dataset.attrib["parameter_units"] = join(fill("1", 34), " | ")
        dataset.attrib["state_order"] = "synthetic test state"
        dataset.attrib["aerosol_coordinate_transform"] = "natural logarithm"
    end
    return path
end

function write_full_truth_table(path)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "# index surface_index surface aerosol_index aerosol_case " *
                    "sif_index sif_case xco2_index xco2_ppm")
        index = 0
        for (surface_index, surface) in enumerate(("urban", "rural", "desert", "forest")),
                (aerosol_index, aerosol) in enumerate(("none", "aod760_0p28")),
                (sif_index, sif) in enumerate(("off", "angular_integral760_0p5")),
                (xco2_index, xco2) in enumerate((380, 400, 420, 440))
            index += 1
            println(io, join((index, surface_index, surface, aerosol_index,
                              aerosol, sif_index, sif, xco2_index, xco2), ' '))
        end
    end
    return path
end

function synthetic_release!(paths; receipt_root=paths.full_truth_root)
    table = write_full_truth_table(joinpath(paths.full_truth_root, "true_states.dat"))
    validation = joinpath(paths.restart_root, "sif_v2_release_validation.dat")
    open(validation, "w") do io
        println(io, "# corrected SIF truth release validation")
        println(io, "# release_schema 1")
        println(io, "# sif_definition_version 2")
        println(io, "# sif_case angular_integral760_0p5")
        println(io, "# corrected_state_table_sha256 ", file_sha256(table))
        println(io, "# index class staged_sha256 archived_sha256 canonical_destination")
        for index in vcat(
                [collect(first_state:(first_state + 3))
                 for first_state in (5, 13, 21, 29, 37, 45, 53, 61)]...)
            block_position = mod(index - 1, 16) + 1
            class = block_position <= 8 ? "clear" : "aerosol"
            directory = class == "clear" ? paths.full_truth_root :
                joinpath(paths.full_truth_root, "aerosol_chunked")
            mkpath(directory)
            destination = joinpath(directory, "hiressim_$(lpad(index, 3, '0')).nc")
            write(destination, "published corrected SIF state $index\n")
            hash = file_sha256(destination)
            receipt_destination = joinpath(
                receipt_root, class == "clear" ? "" : "aerosol_chunked",
                basename(destination))
            println(io, "$index $class $hash $(repeat("b", 64)) " *
                        receipt_destination)
        end
    end
    validation_hash = file_sha256(validation)
    complete = joinpath(paths.full_truth_root, "sif_v2_release_complete.dat")
    open(complete, "w") do io
        println(io, "# corrected SIF truth publication complete")
        println(io, "release_id ", validation_hash[1:16])
        println(io, "validation_receipt_sha256 ", validation_hash)
    end
    return validation, complete
end

function write_bottom_truth_table(path)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "# index surface_index surface aerosol_index aerosol_case " *
                    "sif_index sif_case bottom_co2_index background_co2_ppm " *
                    "bottom_layer_index bottom_co2_ppm xco2_ppm campaign")
        index = 0
        for (surface_index, surface) in enumerate(("urban", "rural", "desert", "forest")),
                (aerosol_index, aerosol) in enumerate(("none", "aod760_0p28")),
                (sif_index, sif) in enumerate(("off", "angular_integral760_0p5")),
                (co2_index, bottom_co2) in enumerate((360, 380, 400, 420, 440))
            index += 1
            xco2 = 400 + 0.06229780735737046 * (bottom_co2 - 400)
            println(io, join((index, surface_index, surface, aerosol_index,
                              aerosol, sif_index, sif, co2_index, 400.0, 16,
                              bottom_co2, xco2, "bottom_layer_XCO2"), ' '))
        end
    end
    return path
end

function synthetic_bottom_release!(paths, full_release)
    table = write_bottom_truth_table(paths.bottom_truth_table)
    mkpath(paths.measurement_directory)
    mkpath(paths.noise_directory)
    receipt = joinpath(
        dirname(paths.bottom_truth_table),
        "bottom_layer_sif_v2_release_complete.dat")
    rows = Tuple[]
    for index in EXPECTED_SIF_STATES
        block_position = mod(index - 1, 20) + 1
        truth_directory = block_position <= 10 ?
            dirname(paths.bottom_truth_table) :
            joinpath(dirname(paths.bottom_truth_table), "aerosol_chunked")
        mkpath(truth_directory)
        truth = joinpath(truth_directory,
                         "hiressim_$(lpad(index, 3, '0')).nc")
        measurement = joinpath(paths.measurement_directory,
                               "OCO2sims_$(lpad(index, 3, '0')).nc")
        noise = joinpath(paths.noise_directory,
                         "OCO2noise_$(lpad(index, 3, '0')).nc")
        write(truth, "truth $index\n")
        write(measurement, "measurement $index\n")
        write(noise, "noise $index\n")
        push!(rows, (index, file_sha256(truth), file_sha256(measurement),
                     file_sha256(noise)))
    end
    bottom_cases = GattacaTaperedSIFReadiness.read_truth_cases(table)
    for truth_case in GattacaTaperedSIFReadiness.select_sif_truth_cases(
            bottom_cases, :off)
        index = truth_case.state_index
        truth = GattacaTaperedSIFReadiness._truth_scene(paths, truth_case)
        measurement = joinpath(paths.measurement_directory,
                               "OCO2sims_$(lpad(index, 3, '0')).nc")
        noise = joinpath(paths.noise_directory,
                         "OCO2noise_$(lpad(index, 3, '0')).nc")
        mkpath(dirname(truth))
        write(truth, "preserved no-SIF truth $index\n")
        write(measurement, "preserved no-SIF measurement $index\n")
        write(noise, "preserved no-SIF noise $index\n")
    end
    archive_root = joinpath(
        paths.private_root, "archive", "sif_wavelength_integral_0p5_20260904",
        "bottom_layer_XCO2_retrievals")
    legacy_table = joinpath(archive_root, "truth", "true_states.dat")
    archive_manifest = joinpath(
        archive_root, "bottom_layer_legacy_manifest.dat")
    mkpath(dirname(legacy_table))
    cp(table, legacy_table)
    write(archive_manifest, "# synthetic legacy manifest\n")
    no_sif_hash =
        GattacaTaperedSIFReadiness._no_sif_triplet_set_sha256(paths)
    open(receipt, "w") do io
        println(io, "# bottom-layer SIF-v2 release complete")
        println(io, "# release_schema 1")
        println(io, "# sif_definition_version 2")
        println(io, "# full_column_release_receipt_sha256 ",
                full_release.complete_receipt_sha256)
        println(io, "# full_column_state_table_sha256 ",
                full_release.table_sha256)
        println(io, "# bottom_state_table_sha256 ", file_sha256(table))
        println(io, "# legacy_bottom_state_table_sha256 ",
                file_sha256(legacy_table))
        println(io, "# legacy_bottom_state_table_archive_relative truth/true_states.dat")
        println(io, "# legacy_archive_manifest_relative bottom_layer_legacy_manifest.dat")
        println(io, "# legacy_archive_manifest_sha256 ",
                file_sha256(archive_manifest))
        println(io, "# no_sif_byte_preservation_policy " *
                    "preserve_truth_measurement_noise_bytes_and_legacy_state_table_hash_attributes")
        println(io, "# no_sif_triplet_set_sha256 ", no_sif_hash)
        println(io, "# input_set_sha256 ", repeat("c", 64))
        println(io, "# state truth_sha256 measurement_sha256 noise_sha256")
        for row in rows
            println(io, join(row, ' '))
        end
    end
    return receipt
end

@testset "Gattaca tapered-SIF array is SIF-only and complete" begin
    states = Int[]
    for task in 0:39
        command = setenv(`bash $LAUNCHER`,
            "SLURM_ARRAY_TASK_ID" => string(task),
            "GATTACA_TAPERED_PLAN_ONLY" => "1")
        output = read(command, String)
        matched = match(r"state=(\d{3})", output)
        @test matched !== nothing
        push!(states, parse(Int, only(matched.captures)))
        @test occursin("sif_filter=on", output)
        @test occursin("prior_model=$(EXPECTED_PRIOR_MODEL)", output)
        @test occursin("campaign=$(CAMPAIGN_ID)", output)
    end
    @test states == EXPECTED_SIF_STATES
    @test length(unique(states)) == 40

    command = setenv(`bash $LAUNCHER`,
        "SLURM_ARRAY_TASK_ID" => "40",
        "GATTACA_TAPERED_PLAN_ONLY" => "1")
    @test !success(pipeline(command; stdout=devnull, stderr=devnull))
end

@testset "tracked Gattaca helpers do not publish private scheduler identity" begin
    launcher_source = read(LAUNCHER, String)
    packager_source = read(STATIC_PACKAGER, String)
    @test success(`bash -n $LAUNCHER`)
    @test success(`bash -n $STATIC_PACKAGER`)
    @test isnothing(match(r"/home/[A-Za-z0-9_.-]+", launcher_source))
    for private_text in (
            "#SBATCH --account", "#SBATCH --qos", "#SBATCH --output",
            "#SBATCH --error")
        @test !occursin(private_text, launcher_source)
    end
    @test isnothing(match(r"/home/[A-Za-z0-9_.-]+", packager_source))
    @test occursin("EXPECTED_SIF_OWNERSHIP_SHA256", packager_source)
    @test occursin("require_packager_at_head", packager_source)
    @test occursin("rev-parse --is-inside-work-tree", packager_source)
    @test !occursin("[[ -d \"\$repo_root/.git\" ]]", packager_source)
    @test occursin("git -C \"\$repo_root\" hash-object", packager_source)
    @test occursin("verify_packaged_archives", packager_source)
    @test occursin("bundle_schema=2", packager_source)
    @test occursin("collect_archive_files", packager_source)
    @test occursin("outside its 35-member allowlist", packager_source)
    @test occursin("bundle output must be outside the Git checkout",
                   packager_source)
    for required_private_input in (
            "RRS_XCO2/inversion/retrieval_setup/co2_prior_covariances.dat",
            "RRS_XCO2/bottom_layer_XCO2_retrievals/retrieval_setup/apriori_states.nc",
            "RRS_XCO2/surface_albedos/lambertian_legendre_inputs.dat")
        @test occursin(required_private_input, packager_source)
    end

    package_start = first(findfirst("package_bundle()", packager_source))
    package_end = first(findnext(
        "safe_tar_members()", packager_source, package_start)) - 1
    package_source = packager_source[package_start:package_end]
    @test !occursin("tracked checkout has unstaged changes", package_source)
    @test !occursin("tracked checkout has staged changes", package_source)

    install_start = first(findfirst("install_bundle()", packager_source))
    install_source = packager_source[install_start:end]
    @test occursin("tracked checkout has unstaged changes", install_source)
    @test occursin("tracked checkout has staged changes", install_source)
end

@testset "cross-host provenance paths retain canonical semantic tails" begin
    @test GattacaTaperedSIFReadiness._has_path_tail(
        "/publisher/code/github/uni_vSmartMOM/RRS_XCO2/truth_map/" *
        "aerosol_chunked/hiressim_013.nc",
        "truth_map/aerosol_chunked/hiressim_013.nc")
    @test GattacaTaperedSIFReadiness._has_path_tail(
        "/publisher/code/uni_vSmartMOM/RRS_XCO2/" *
        "bottom_layer_XCO2_retrievals/truth/hiressim_006.nc",
        "bottom_layer_XCO2_retrievals/truth/hiressim_006.nc")
    @test !GattacaTaperedSIFReadiness._has_path_tail(
        "/publisher/other/hiressim_006.nc",
        "bottom_layer_XCO2_retrievals/truth/hiressim_006.nc")
end

@testset "tapered prior identity is path/model/hash pinned" begin
    mktempdir() do root
        paths = test_paths(root)
        reference_path = joinpath(
            paths.bottom_campaign_root, "retrieval_setup", "apriori_states.nc")
        write_prior(reference_path; model=EXPECTED_PRIOR_BASE_MODEL)
        reference_hash = file_sha256(reference_path)
        write_prior(paths.prior_path)
        hash = file_sha256(paths.prior_path)
        identity = validate_prior_identity(
            paths.prior_path, hash;
            reference_path, reference_sha256=reference_hash)
        @test identity.model == EXPECTED_PRIOR_MODEL
        @test identity.sha256 == hash
        @test_throws ErrorException validate_prior_identity(
            paths.prior_path, repeat("0", 64);
            reference_path, reference_sha256=reference_hash)
        @test_throws ErrorException validate_prior_identity(
            paths.prior_path, hash;
            reference_path, reference_sha256=repeat("0", 64))

        wrong_name = joinpath(root, "apriori_states_wrong_name.nc")
        cp(paths.prior_path, wrong_name)
        @test_throws ErrorException validate_prior_identity(
            wrong_name, file_sha256(wrong_name);
            reference_path, reference_sha256=reference_hash)

        rm(paths.prior_path)
        write_prior(paths.prior_path; model="acos_mapped")
        @test_throws ErrorException validate_prior_identity(
            paths.prior_path, file_sha256(paths.prior_path);
            reference_path, reference_sha256=reference_hash)

        rm(paths.prior_path)
        write_prior(paths.prior_path; alter_nonco2=true)
        @test_throws ErrorException validate_prior_identity(
            paths.prior_path, file_sha256(paths.prior_path);
            reference_path, reference_sha256=reference_hash)

        rm(paths.prior_path)
        write_prior(paths.prior_path; alter_xa=true)
        @test_throws ErrorException validate_prior_identity(
            paths.prior_path, file_sha256(paths.prior_path);
            reference_path, reference_sha256=reference_hash)

        rm(paths.prior_path)
        write_prior(paths.prior_path; alter_nonadjacent_co2=true)
        @test_throws ErrorException validate_prior_identity(
            paths.prior_path, file_sha256(paths.prior_path);
            reference_path, reference_sha256=reference_hash)

        for options in (
                (; alter_co2_diagonal=true),
                (; alter_active_mapping=true),
                (; alter_prior_sigma=true),
                (; aerosol_ln_aod_sigma=0.76),
                (; surface_p1_sigmas=(0.002, 0.0021, 0.002)),
                (; surface_p2_sigmas=(0.002, 0.002, 0.0021)),
                (; sif_slope_mean=1e-6),
                (; sif_slope_sigma=0.002626))
            rm(paths.prior_path)
            write_prior(paths.prior_path; options...)
            @test_throws ErrorException validate_prior_identity(
                paths.prior_path, file_sha256(paths.prior_path);
                reference_path, reference_sha256=reference_hash)
        end
    end
end

@testset "corrected full-column SIF release is cryptographically linked" begin
    mktempdir() do root
        paths = test_paths(root)
        synthetic_release!(
            paths;
            receipt_root="/publisher/checkout/RRS_XCO2/truth_map")
        release = validate_published_sif_release(paths)
        @test length(release.source_hashes) == 32
        @test length(release.release_id) == 16

        full_marker = joinpath(
            paths.full_truth_root, ".sif_v2_publication_in_progress")
        write(full_marker, "synthetic interrupted publication\n")
        @test_throws ErrorException validate_published_sif_release(paths)
        rm(full_marker)

        scene = joinpath(paths.full_truth_root, "hiressim_005.nc")
        write(scene, "changed after publication\n")
        @test_throws ErrorException validate_published_sif_release(paths)
    end
end

@testset "bottom-layer SIF truth, OCO, and noise have a release barrier" begin
    mktempdir() do root
        paths = test_paths(root)
        synthetic_release!(paths)
        full_release = validate_published_sif_release(paths)
        synthetic_bottom_release!(paths, full_release)
        bottom_release = validate_published_bottom_sif_release(
            paths, full_release)
        @test length(bottom_release.rows) == 40
        @test bottom_release.bottom_table_sha256 ==
              file_sha256(paths.bottom_truth_table)

        preserved = joinpath(
            paths.noise_directory, "OCO2noise_001.nc")
        preserved_bytes = read(preserved)
        write(preserved, "changed preserved no-SIF input\n")
        @test_throws ErrorException validate_published_bottom_sif_release(
            paths, full_release)
        write(preserved, preserved_bytes)

        bottom_marker = joinpath(
            dirname(paths.bottom_truth_table),
            ".bottom_layer_sif_v2_publication_in_progress")
        write(bottom_marker, "synthetic interrupted publication\n")
        @test_throws ErrorException validate_published_bottom_sif_release(
            paths, full_release)
        rm(bottom_marker)

        changed = joinpath(paths.noise_directory, "OCO2noise_006.nc")
        write(changed, "changed after bottom-layer release\n")
        @test_throws ErrorException validate_published_bottom_sif_release(
            paths, full_release)
    end
end

@testset "isolated output and immutable campaign identity" begin
    mktempdir() do root
        paths = test_paths(root)
        @test validate_output_isolation(paths) == paths.output_root
        @test validate_prior_isolation(paths) == paths.prior_path
        wrong = ReadinessPaths(
            paths.repo_root, paths.private_root, paths.full_truth_root,
            paths.restart_root, paths.bottom_campaign_root,
            paths.bottom_truth_table, paths.measurement_directory,
            paths.noise_directory, paths.prior_path,
            joinpath(paths.private_root, "results", "old", "retrievals"))
        @test_throws ErrorException validate_output_isolation(wrong)
        wrong_prior = ReadinessPaths(
            paths.repo_root, paths.private_root, paths.full_truth_root,
            paths.restart_root, paths.bottom_campaign_root,
            paths.bottom_truth_table, paths.measurement_directory,
            paths.noise_directory,
            joinpath(paths.bottom_campaign_root, "retrieval_setup",
                     basename(paths.prior_path)), paths.output_root)
        @test_throws ErrorException validate_prior_isolation(wrong_prior)

        readiness = (
            prior=(model=EXPECTED_PRIOR_MODEL, sha256=repeat("1", 64),
                   path=paths.prior_path,
                   reference_path=joinpath(root, "current.nc"),
                   reference_sha256=repeat("9", 64)),
            release=(release_id=repeat("2", 16),
                     validation_receipt_sha256=repeat("3", 64)),
            bottom_release=(receipt_sha256=repeat("7", 64),
                            input_set_sha256=repeat("8", 64),
                            no_sif_triplet_set_sha256=repeat("6", 64),
                            archive_manifest_sha256=repeat("a", 64)),
            inputs=(input_set_sha256=repeat("4", 64),
                    table_sha256=repeat("5", 64)),
        )
        identity = initialize_campaign_identity!(
            paths, readiness, SOURCE_SHA)
        @test isfile(identity)
        @test identity == initialize_campaign_identity!(
            paths, readiness, SOURCE_SHA)
        changed = merge(readiness, (inputs=(
            input_set_sha256=repeat("6", 64),
            table_sha256=repeat("5", 64)),))
        @test_throws ErrorException initialize_campaign_identity!(
            paths, changed, SOURCE_SHA)
    end

    mktempdir() do root
        paths = test_paths(root)
        directory = joinpath(paths.output_root, "corrected")
        mkpath(directory)
        write(joinpath(directory, "retrieval_state006_perturbation11.nc"),
              "unlabelled prior campaign")
        readiness = (
            prior=(model=EXPECTED_PRIOR_MODEL, sha256=repeat("1", 64),
                   path=paths.prior_path,
                   reference_path=joinpath(root, "current.nc"),
                   reference_sha256=repeat("9", 64)),
            release=(release_id=repeat("2", 16),
                     validation_receipt_sha256=repeat("3", 64)),
            bottom_release=(receipt_sha256=repeat("7", 64),
                            input_set_sha256=repeat("8", 64),
                            no_sif_triplet_set_sha256=repeat("6", 64),
                            archive_manifest_sha256=repeat("a", 64)),
            inputs=(input_set_sha256=repeat("4", 64),
                    table_sha256=repeat("5", 64)),
        )
        @test_throws ErrorException initialize_campaign_identity!(
            paths, readiness, SOURCE_SHA)
    end
end
