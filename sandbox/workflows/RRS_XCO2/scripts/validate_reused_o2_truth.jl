#!/usr/bin/env julia

"""Validate bit-exact A-band reuse and completed CO2 regeneration."""

using NCDatasets
using Printf

const ROOT = normpath(joinpath(@__DIR__, ".."))
const ACTIVE_ROOT = joinpath(ROOT, "truth_map")
const ARCHIVE_ROOT = joinpath(
    ROOT, "obsolete", "pre_absco_closure_20260829", "truth_map")
const O2_VARIABLES = (
    "radiance_rayleigh_o2a",
    "radiance_cabannes_o2a",
    "radiance_rrs_o2a",
)
const CO2_VARIABLES = (
    "radiance_rayleigh_weak_co2",
    "radiance_rayleigh_strong_co2",
)

function scene_files(root)
    isdir(root) || error("missing truth directory: $root")
    return sort!(filter(path -> occursin(r"hiressim_\d{3}\.nc$", basename(path)),
                        readdir(root; join=true)))
end

function validate_subset(active_dir, archive_dir, expected_count)
    active = scene_files(active_dir)
    length(active) == expected_count || error(
        "expected $expected_count scenes in $active_dir, found $(length(active))")
    worst = 0.0
    for destination in active
        source = joinpath(archive_dir, basename(destination))
        isfile(source) || error("missing A-band source: $source")
        NCDataset(destination) do current
            get(current.attrib, "simulation_complete", 0) == 1 || error(
                "scene is incomplete: $destination")
            get(current.attrib, "co2_absco_regeneration_complete", 0) == 1 ||
                error("CO2 regeneration is incomplete: $destination")
            get(current.attrib, "o2_truth_reused", 0) == 1 || error(
                "scene lacks A-band reuse provenance: $destination")
            abspath(String(current.attrib["o2_truth_reuse_source"])) == abspath(source) ||
                error("A-band reuse source mismatch: $destination")
            NCDataset(source) do archived
                for variable in O2_VARIABLES
                    old = archived[variable][:, :]
                    new = current[variable][:, :]
                    size(old) == size(new) || error(
                        "$variable shape changed in $destination")
                    old == new || error(
                        "$variable is not bit-identical to $source")
                    worst = max(worst, maximum(abs, Float64.(new .- old)))
                end
            end
            for variable in CO2_VARIABLES
                values = Float64.(nomissing(current[variable][:, :], NaN))
                all(isfinite, values) || error(
                    "$variable contains non-finite values in $destination")
                maximum(abs, values) < 1e30 || error(
                    "$variable contains unwritten fill values in $destination")
            end
        end
    end
    return worst
end

function main()
    clear_worst = validate_subset(ACTIVE_ROOT, ARCHIVE_ROOT, 32)
    aerosol_worst = validate_subset(
        joinpath(ACTIVE_ROOT, "aerosol_chunked"),
        joinpath(ARCHIVE_ROOT, "aerosol_chunked"), 32)
    @printf("clear scenes: 32; maximum A-band copy difference = %.9e\n", clear_worst)
    @printf("aerosol scenes: 32; maximum A-band copy difference = %.9e\n", aerosol_worst)
    println("CO2 regeneration and A-band reuse validation: PASS")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
