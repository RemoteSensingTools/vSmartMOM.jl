#!/usr/bin/env julia

"""
Prepare an isolated, resumable direct-RT restart of the full-column SIF-on
O2 A-band truth calculations.

Run `generate_truth_states.jl` with `TRUTH_OUT` set to the same staging root
before this script. This producer then:

1. validates the legacy and corrected 64-state tables;
2. copies the legacy table and all 32 legacy SIF-on scenes into a clearly
   named obsolete archive, with SHA-256 checksums;
3. seeds disjoint clear/aerosol staging directories with copies of those
   scenes so `regenerate_o2_preserve_co2.jl` can replace only their O2 arrays
   while retaining the accepted CO2 arrays.

Canonical truth files are not modified here. Existing staging files are left
untouched so interrupted direct calculations can resume from their versioned
checkpoints.

Environment controls:

- `FULL_COLUMN_TRUTH_ROOT` (default `RRS_XCO2/truth_map`)
- `SIF_RESTART_ROOT` (default `<truth root>/.sif_v2_restart`)
- `SIF_OBSOLETE_ROOT` (default
  `RRS_XCO2/obsolete/sif_wavelength_integral_0p5_20260904`)
"""

using DelimitedFiles
using NCDatasets
using Printf
using SHA
using UUIDs
using vSmartMOM

include(joinpath(@__DIR__, "common.jl"))
using .RRSXCO2Common

const RRS_ROOT = normpath(joinpath(@__DIR__, ".."))
const TRUTH_ROOT = abspath(get(
    ENV, "FULL_COLUMN_TRUTH_ROOT", joinpath(RRS_ROOT, "truth_map")))
const RESTART_ROOT = abspath(get(
    ENV, "SIF_RESTART_ROOT", joinpath(TRUTH_ROOT, ".sif_v2_restart")))
const OBSOLETE_ROOT = abspath(get(
    ENV,
    "SIF_OBSOLETE_ROOT",
    joinpath(RRS_ROOT, "obsolete", "sif_wavelength_integral_0p5_20260904"),
))
const LEGACY_CASE = "total_0p5"
const LEGACY_TABLE = joinpath(TRUTH_ROOT, "true_states.dat")
const CORRECTED_TABLE = joinpath(RESTART_ROOT, "true_states.dat")

file_sha256(path::AbstractString) = open(path, "r") do io
    bytes2hex(sha256(io))
end

function sif_pairs()
    pairs = Tuple{Int,Int}[]
    for surface_offset in 0:16:48, aerosol_offset in (0, 8), xco2 in 1:4
        off = surface_offset + aerosol_offset + xco2
        push!(pairs, (off, off + 4))
    end
    return pairs
end

aerosol_on(index::Integer) = ((index - 1) % 16) >= 8

function canonical_scene(index::Integer)
    directory = aerosol_on(index) ? joinpath(TRUTH_ROOT, "aerosol_chunked") :
                                    TRUTH_ROOT
    return joinpath(directory, @sprintf("hiressim_%03d.nc", index))
end

function staged_scene(index::Integer)
    directory = joinpath(RESTART_ROOT, aerosol_on(index) ? "aerosol" : "clear")
    return joinpath(directory, @sprintf("hiressim_%03d.nc", index))
end

function archived_scene(index::Integer)
    directory = joinpath(OBSOLETE_ROOT, "truth_map",
                         aerosol_on(index) ? "aerosol_chunked" : "")
    return joinpath(directory, @sprintf("hiressim_%03d.nc", index))
end

function copy_verified(source, destination)
    isfile(source) || error("missing source file: $source")
    if isfile(destination)
        file_sha256(destination) == file_sha256(source) || error(
            "existing archive differs from source: $destination")
        return destination
    end
    mkpath(dirname(destination))
    temporary = destination * ".tmp.$(getpid()).$(uuid4())"
    cp(source, temporary)
    file_sha256(temporary) == file_sha256(source) || error(
        "copy verification failed: $source -> $temporary")
    mv(temporary, destination)
    return destination
end

function validate_tables()
    isfile(LEGACY_TABLE) || error("missing legacy state table: $LEGACY_TABLE")
    isfile(CORRECTED_TABLE) || error(
        "missing corrected staging table: $CORRECTED_TABLE; first run " *
        "TRUTH_OUT=$RESTART_ROOT generate_truth_states.jl")
    old = readdlm(LEGACY_TABLE; comments=true)
    new = readdlm(CORRECTED_TABLE; comments=true)
    size(old) == (64, 27) || error("legacy table is not 64x27")
    size(new) == (64, 27) || error("corrected table is not 64x27")
    corrected = RRSXCO2Common.campaign_sif_state()
    legacy = sif_reference_state(total_sif=0.5, reference_wavelength_nm=760)
    on_indices = Set(last(pair) for pair in sif_pairs())
    fixed_columns = vcat(collect(1:6), collect(8:15), collect(19:27))
    for (row_old, row_new) in zip(eachrow(old), eachrow(new))
        Int(row_old[1]) == Int(row_new[1]) || error("table indices differ")
        index = Int(row_new[1])
        is_on = index in on_indices
        all(row_old[column] == row_new[column] for column in fixed_columns) ||
            error("non-SIF truth fields changed in state $index")
        if is_on
            String(row_old[7]) == LEGACY_CASE || error(
                "state $index is not a legacy SIF-on row")
            String(row_new[7]) == RRSXCO2Common.SIF_CASE_ON || error(
                "state $index lacks the corrected SIF label")
            isapprox(Float64(row_old[17]), legacy.SIF760;
                     atol=5e-15, rtol=0) || error(
                "state $index has an unexpected legacy SIF760")
            Float64(row_new[16]) ==
                RRSXCO2Common.SIF_ANGULAR_INTEGRAL_760 || error(
                "state $index has the wrong corrected angular integral")
            isapprox(Float64(row_new[17]), corrected.SIF760;
                     atol=5e-15, rtol=0) || error(
                "state $index has the wrong corrected SIF760")
            isapprox(Float64(row_new[18]), corrected.mSIF;
                     atol=5e-16, rtol=0) || error(
                "state $index has the wrong corrected mSIF")
        else
            String(row_old[7]) == "off" && String(row_new[7]) == "off" ||
                error("state $index has inconsistent no-SIF labels")
            all(iszero, Float64.(row_new[16:18])) || error(
                "state $index has nonzero corrected-table SIF fields")
        end
    end
    return (; legacy, corrected)
end

function validate_legacy_scene(index, path)
    NCDataset(path, "r") do dataset
        Int(get(dataset.attrib, "simulation_complete", 0)) == 1 || error(
            "legacy truth scene is incomplete: $path")
        Int(get(dataset.attrib, "state_index", -1)) == index || error(
            "legacy truth state-index mismatch: $path")
        String(get(dataset.attrib, "sif_case", "")) == LEGACY_CASE || error(
            "legacy truth SIF label mismatch: $path")
        for variable in ("radiance_rayleigh_o2a", "radiance_cabannes_o2a",
                         "radiance_rrs_o2a", "radiance_rayleigh_weak_co2",
                         "radiance_rayleigh_strong_co2")
            haskey(dataset, variable) || error("$path lacks $variable")
            all(isfinite, dataset[variable][:, :]) || error(
                "$path contains non-finite $variable")
        end
    end
    return path
end

function write_records(pairs, states)
    mkpath(OBSOLETE_ROOT)
    readme = joinpath(OBSOLETE_ROOT, "README.md")
    if !isfile(readme)
        open(readme, "w") do io
            println(io, "# Superseded wavelength-integrated SIF campaign")
            println(io)
            println(io, "The archived SIF-on inputs used a wavelength-integrated ",
                        "template area of 0.5 mW m^-2 sr^-1. They are retained ",
                        "only for provenance and must not be used for retrievals.")
            println(io)
            println(io, "The replacement campaign directly recomputes O2 truth ",
                        "with 2pi*L_lambda(760 nm)=0.5 mW m^-2 nm^-1, so each ",
                        "isotropic upwelling BOA stream has L_lambda=0.5/(2pi) ",
                        "mW m^-2 sr^-1 nm^-1.")
            println(io)
            println(io, "No-SIF truth and retrieval products are not invalidated ",
                        "by this normalization correction.")
        end
    end
    manifest = joinpath(OBSOLETE_ROOT, "full_column_input_manifest.dat")
    open(manifest, "w") do io
        println(io, "# off_state legacy_sif_state legacy_scene_sha256")
        for (off, on) in pairs
            @printf(io, "%03d %03d %s\n", off, on,
                    file_sha256(archived_scene(on)))
        end
        println(io, "# legacy_true_states_sha256 ",
                file_sha256(joinpath(OBSOLETE_ROOT, "truth_map", "true_states.dat")))
        println(io, "# corrected_staging_true_states_sha256 ",
                file_sha256(CORRECTED_TABLE))
        println(io, "# legacy_SIF760 ", states.legacy.SIF760)
        println(io, "# corrected_SIF760 ", states.corrected.SIF760)
        println(io, "# corrected_mSIF ", states.corrected.mSIF)
    end
    return (; readme, manifest)
end

function main()
    states = validate_tables()
    pairs = sif_pairs()
    copy_verified(
        LEGACY_TABLE, joinpath(OBSOLETE_ROOT, "truth_map", "true_states.dat"))
    seeded = 0
    resumed = 0
    for (_, on) in pairs
        source = canonical_scene(on)
        validate_legacy_scene(on, source)
        copy_verified(source, archived_scene(on))
        stage = staged_scene(on)
        if isfile(stage)
            resumed += 1
        else
            mkpath(dirname(stage))
            temporary = stage * ".tmp.$(getpid()).$(uuid4())"
            cp(source, temporary)
            file_sha256(temporary) == file_sha256(source) || error(
                "staging copy verification failed for state $on")
            mv(temporary, stage)
            seeded += 1
        end
    end
    records = write_records(pairs, states)
    @info("prepared direct SIF truth restart",
          pairs=length(pairs), seeded, resumed,
          corrected_table=CORRECTED_TABLE,
          clear_stage=joinpath(RESTART_ROOT, "clear"),
          aerosol_stage=joinpath(RESTART_ROOT, "aerosol"),
          archive=OBSOLETE_ROOT,
          records)
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
