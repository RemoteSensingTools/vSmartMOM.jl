using Test
using NCDatasets
using Printf
using SHA
using vSmartMOM: sif_reference_state

include(joinpath(@__DIR__, "..", "scripts",
                 "validate_publish_sif_truth_restart.jl"))
using .SIFTruthReleaseGate

const Gate = SIFTruthReleaseGate
const N_O2_TEST = 11
const N_WEAK_TEST = 5
const N_STRONG_TEST = 6

function state_fields(index)
    surface_index = fld(index - 1, 16) + 1
    local_index = mod(index - 1, 16) + 1
    aerosol_index = local_index > 8 ? 2 : 1
    within_aerosol = mod(local_index - 1, 8) + 1
    sif_index = within_aerosol > 4 ? 2 : 1
    xco2_index = mod(index - 1, 4) + 1
    surfaces = ("urban", "rural", "desert", "forest")
    xco2 = (380, 400, 420, 440)
    return (; index, surface_index, surface=surfaces[surface_index],
            aerosol_index,
            aerosol_case=aerosol_index == 1 ? "none" : "aod760_0p28",
            sif_index, xco2_index, xco2_ppm=xco2[xco2_index])
end

function write_state_table(path; corrected)
    mkpath(dirname(path))
    current = Gate.RRSXCO2Common.campaign_sif_state()
    legacy = sif_reference_state(total_sif=0.5,
                                 reference_wavelength_nm=760)
    open(path, "w") do io
        println(io, "# synthetic 64-state release-gate fixture")
        for index in 1:64
            state = state_fields(index)
            on = state.sif_index == 2
            sif_case = on ? (corrected ?
                Gate.RRSXCO2Common.SIF_CASE_ON : "total_0p5") : "off"
            control = on ? 0.5 : 0.0
            sif760 = on ? (corrected ? current.SIF760 : legacy.SIF760) : 0.0
            msif = on ? (corrected ? current.mSIF : legacy.mSIF) : 0.0
            aods = state.aerosol_index == 2 ? (0.36, 0.08, 0.01) : (0.0, 0.0, 0.0)
            values = Any[
                index, state.surface_index, state.surface,
                state.aerosol_index, state.aerosol_case,
                state.sif_index, sif_case, state.xco2_index, state.xco2_ppm,
                1000.0, 30.0, 0.0, aods..., control, sif760, msif,
                0.2, 0.0, 0.0, 0.3, 0.0, 0.0, 0.4, 0.0, 0.0,
            ]
            println(io, join(values, ' '))
        end
    end
end

function o2_components(state; corrected=false)
    base = Float32(state.surface_index + 0.25state.aerosol_index)
    x = Float32.(1:N_O2_TEST)
    off = [base + Float32(0.01r) + Float32(0.001) * value
           for r in 1:3, value in x]
    delta0 = [Float32(2e-4r + 1e-5j)
              for r in 1:3, j in 1:N_O2_TEST]
    scale = state.sif_index == 1 ? 0.0f0 :
            (corrected ? Float32(Gate.PHYSICAL_SCALE) : 1.0f0)
    return Tuple(off .+ Float32(component) .* scale .* delta0
                 for component in (1, 2, 3))
end

function co2_components(state)
    weak = fill(Float32(state.xco2_ppm / 1000 + state.surface_index / 100),
                3, N_WEAK_TEST)
    strong = fill(Float32(state.xco2_ppm / 900 + state.aerosol_index / 100),
                  3, N_STRONG_TEST)
    return weak, strong
end

function create_scene(path, state; corrected=false, source_table=nothing)
    mkpath(dirname(path))
    o2 = o2_components(state; corrected)
    weak, strong = co2_components(state)
    NCDataset(path, "c") do dataset
        defDim(dataset, "stokes", 3)
        defDim(dataset, "o2a", N_O2_TEST)
        defDim(dataset, "weak", N_WEAK_TEST)
        defDim(dataset, "strong", N_STRONG_TEST)
        for (name, values) in zip(Gate.O2_VARIABLES, o2)
            defVar(dataset, name, Float32, ("stokes", "o2a"))[:, :] = values
        end
        defVar(dataset, Gate.CO2_VARIABLES[1], Float32,
               ("stokes", "weak"))[:, :] = weak
        defVar(dataset, Gate.CO2_VARIABLES[2], Float32,
               ("stokes", "strong"))[:, :] = strong
        attributes = dataset.attrib
        attributes["state_index"] = Int32(state.index)
        attributes["surface"] = state.surface
        attributes["aerosol_case"] = state.aerosol_case
        attributes["xco2_ppm"] = Float64(state.xco2_ppm)
        attributes["simulation_complete"] = Int32(1)
        attributes["co2_absco_regeneration_complete"] = Int32(1)
        attributes["o2_absco_regeneration_complete"] = Int32(1)
        attributes["o2_truth_regenerated"] = Int32(1)
        attributes["o2_truth_reused"] = Int32(0)
        attributes["o2_truth_reuse_source"] =
            "none; regenerated from current model"
        attributes["o2_core_grid_version"] = Int32(2)
        attributes["o2_nstreams"] = Int32(state.aerosol_index == 2 ? 9 : 5)
        attributes["o2_chunk_points"] =
            Int32(state.aerosol_index == 2 ? 64 : N_O2_TEST)
        attributes["o2_raman_shoulder_cm-1"] = 234.0
        attributes["spectroscopy_database"] = "ABSCO"
        attributes["spectroscopy_version"] = "5.2"
        attributes["o2_vmr"] = 0.21
        attributes["psurf_hpa"] = 1000.0
        attributes["atmospheric_layers"] = Int32(16)
        state.aerosol_index == 2 &&
            (attributes["chunked_simulation_complete"] = Int32(1))
        if corrected
            Gate.RRSXCO2Common.write_sif_provenance!(attributes, true)
            attributes["sif_case"] = Gate.RRSXCO2Common.SIF_CASE_ON
            attributes["source_state_table"] = abspath(source_table)
            attributes["source_state_table_sha256"] = Gate.file_sha256(source_table)
        else
            attributes["sif_case"] = state.sif_index == 2 ? "total_0p5" : "off"
            state.sif_index == 2 &&
                (attributes["sif_total_mW_m-2_sr-1"] = 0.5)
        end
    end
end

function fixture(root)
    truth = joinpath(root, "truth_map")
    restart = joinpath(truth, ".sif_v2_restart")
    obsolete = joinpath(root, "obsolete")
    paths = ReleasePaths(truth_root=truth, restart_root=restart,
                         obsolete_root=obsolete)
    write_state_table(Gate.canonical_table(paths); corrected=false)
    write_state_table(Gate.corrected_table(paths); corrected=true)
    mkpath(dirname(Gate.archived_table(paths)))
    cp(Gate.canonical_table(paths), Gate.archived_table(paths))

    selected = Int[]
    pairs = Tuple{Int,Int,String}[]
    for index in 1:64
        state = state_fields(index)
        canonical = Gate.canonical_scene(paths, Gate.TableState(
            state.index, state.surface_index, state.surface,
            state.aerosol_index, state.aerosol_case, state.sif_index,
            state.sif_index == 2 ? "total_0p5" : "off",
            state.xco2_index, state.xco2_ppm, 0.0, 0.0, 0.0))
        create_scene(canonical, state; corrected=false)
        if state.sif_index == 2
            push!(selected, index)
            table_state = first(filter(s -> s.index == index,
                first(Gate.read_table(Gate.corrected_table(paths)))))
            archive = Gate.archived_scene(paths, table_state)
            mkpath(dirname(archive))
            cp(canonical, archive)
            stage = Gate.staged_scene(paths, table_state)
            create_scene(stage, state; corrected=true,
                         source_table=Gate.corrected_table(paths))
            push!(pairs, (index - 4, index, Gate.file_sha256(archive)))
        end
    end
    mkpath(dirname(Gate.archive_manifest(paths)))
    open(Gate.archive_manifest(paths), "w") do io
        println(io, "# off_state legacy_sif_state legacy_scene_sha256")
        for (off, on, digest) in pairs
            @printf(io, "%03d %03d %s\n", off, on, digest)
        end
        println(io, "# legacy_true_states_sha256 ",
                Gate.file_sha256(Gate.archived_table(paths)))
        println(io, "# corrected_staging_true_states_sha256 ",
                Gate.file_sha256(Gate.corrected_table(paths)))
    end
    return paths, selected
end

@testset "corrected SIF truth release gate" begin
    mktempdir() do root
        truth_probe = joinpath(root, "path-guard", "truth")
        restart_probe = joinpath(truth_probe, "restart")
        @test_throws ErrorException ReleasePaths(
            truth_root=truth_probe, restart_root=restart_probe,
            obsolete_root=joinpath(root, "path-guard"))
        mkpath(truth_probe)
        truth_alias = joinpath(root, "truth-alias")
        symlink(truth_probe, truth_alias)
        @test_throws ErrorException ReleasePaths(
            truth_root=truth_probe, restart_root=restart_probe,
            obsolete_root=truth_alias)

        paths, selected = fixture(root)
        archive_hashes = Dict(index => Gate.file_sha256(
            Gate.archived_scene(paths, first(filter(
                state -> state.index == index,
                first(Gate.read_table(Gate.corrected_table(paths)))))))
            for index in selected)
        result = validate_release(paths; expected_o2_points=N_O2_TEST)
        @test length(result.selected) == 32
        @test isfile(result.receipt)

        # A changed CO2 element must fail before any canonical publication.
        state = first(result.selected)
        staged = Gate.staged_scene(paths, state)
        original = NCDataset(staged, "r") do dataset
            dataset[Gate.CO2_VARIABLES[1]][1, 1]
        end
        NCDataset(staged, "a") do dataset
            dataset[Gate.CO2_VARIABLES[1]][1, 1] = original + 1.0f0
        end
        canonical_before = Gate.file_sha256(Gate.canonical_scene(paths, state))
        @test_throws ErrorException validate_release(
            paths; expected_o2_points=N_O2_TEST, write_receipt=false)
        @test Gate.file_sha256(Gate.canonical_scene(paths, state)) ==
              canonical_before
        NCDataset(staged, "a") do dataset
            dataset[Gate.CO2_VARIABLES[1]][1, 1] = original
        end

        old_confirmation = get(ENV, "CONFIRM_SIF_V2_PUBLICATION", nothing)
        old_resume = get(ENV, "SIF_RELEASE_RESUME", nothing)
        ENV["CONFIRM_SIF_V2_PUBLICATION"] = Gate.CONFIRMATION
        try
            # A fault after several scene promotions rolls every canonical
            # scene back and leaves the legacy table plus a visible marker.
            @test_throws ErrorException publish_release(
                paths; expected_o2_points=N_O2_TEST,
                _test_fail_after_scenes=7)
            @test Gate.file_sha256(Gate.canonical_table(paths)) ==
                  Gate.file_sha256(Gate.archived_table(paths))
            @test isfile(Gate.publication_marker(paths))
            @test !isdir(Gate.publication_lock(paths))
            for state in result.selected
                @test Gate.file_sha256(Gate.canonical_scene(paths, state)) ==
                      archive_hashes[state.index]
            end

            # Explicit resume revalidates the mixed-state barrier, then
            # finishes the release and publishes the table last.
            ENV["SIF_RELEASE_RESUME"] = "1"
            publish_release(paths; expected_o2_points=N_O2_TEST)
        finally
            if old_confirmation === nothing
                delete!(ENV, "CONFIRM_SIF_V2_PUBLICATION")
            else
                ENV["CONFIRM_SIF_V2_PUBLICATION"] = old_confirmation
            end
            if old_resume === nothing
                delete!(ENV, "SIF_RELEASE_RESUME")
            else
                ENV["SIF_RELEASE_RESUME"] = old_resume
            end
        end
        @test Gate.file_sha256(Gate.canonical_table(paths)) ==
              Gate.file_sha256(Gate.corrected_table(paths))
        @test isfile(Gate.canonical_receipt(paths))
        @test !isfile(Gate.publication_marker(paths))
        @test !isdir(Gate.publication_lock(paths))
        for state in result.selected
            @test Gate.file_sha256(Gate.canonical_scene(paths, state)) ==
                  Gate.file_sha256(Gate.staged_scene(paths, state))
            @test Gate.file_sha256(Gate.archived_scene(paths, state)) ==
                  archive_hashes[state.index]
        end
    end
end
