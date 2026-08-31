#!/usr/bin/env julia

"""Validate all 64 exact-grid O₂ / completed-ABSCO-CO₂ truth scenes."""

using NCDatasets
using Printf

const ROOT = normpath(joinpath(@__DIR__, ".."))
const TRUTH_ROOT = joinpath(ROOT, "truth_map")
const AEROSOL_ROOT = joinpath(TRUTH_ROOT, "aerosol_chunked")
const O2_VARIABLES = (
    "radiance_rayleigh_o2a",
    "radiance_cabannes_o2a",
    "radiance_rrs_o2a",
)
const CO2_VARIABLES = (
    "radiance_rayleigh_weak_co2",
    "radiance_rayleigh_strong_co2",
)
const NO_AEROSOL_STATES = vcat(collect(1:8), collect(17:24),
                               collect(33:40), collect(49:56))
const AEROSOL_STATES = vcat(collect(9:16), collect(25:32),
                            collect(41:48), collect(57:64))

scene_path(index::Integer, aerosol::Bool) = joinpath(
    aerosol ? AEROSOL_ROOT : TRUTH_ROOT, @sprintf("hiressim_%03d.nc", index))

function read_variable(dataset, name)
    haskey(dataset, name) || error("scene is missing $name")
    values = Float32.(nomissing(dataset[name][:, :], NaN32))
    size(values, 1) == 3 || error("$name does not contain I,Q,U rows")
    all(isfinite, values) || error("$name contains non-finite values")
    maximum(abs, values) < 1f30 || error("$name contains unwritten fill values")
    return values
end

function validate_scene(index; aerosol::Bool)
    path = scene_path(index, aerosol)
    isfile(path) || error("missing scene $path")
    return NCDataset(path) do ds
        Int(ds.attrib["state_index"]) == index || error("state mismatch in $path")
        get(ds.attrib, "simulation_complete", 0) == 1 || error(
            "scene is incomplete: $path")
        get(ds.attrib, "co2_absco_regeneration_complete", 0) == 1 || error(
            "CO2 regeneration is incomplete: $path")
        get(ds.attrib, "o2_absco_regeneration_complete", 0) == 1 || error(
            "O2 regeneration is incomplete: $path")
        get(ds.attrib, "o2_truth_regenerated", 0) == 1 || error(
            "scene does not carry regenerated-O2 provenance: $path")
        get(ds.attrib, "o2_truth_reused", 1) == 0 || error(
            "scene is still marked as reused O2 truth: $path")
        Int(ds.attrib["o2_core_grid_version"]) == 2 || error(
            "scene does not use exact retained O2 core nodes: $path")
        String(ds.attrib["spectroscopy_database"]) == "ABSCO" || error(
            "scene does not use ABSCO: $path")
        String(ds.attrib["spectroscopy_version"]) == "5.2" || error(
            "scene uses the wrong ABSCO version: $path")
        Float64(ds.attrib["o2_vmr"]) == 0.21 || error("wrong O2 VMR in $path")
        Float64(ds.attrib["psurf_hpa"]) == 1000.0 || error(
            "wrong surface pressure in $path")
        Int(ds.attrib["atmospheric_layers"]) == 16 || error(
            "wrong layer count in $path")
        if aerosol
            get(ds.attrib, "chunked_simulation_complete", 0) == 1 || error(
                "aerosol scene lacks chunked completion: $path")
            Int(ds.attrib["o2_nstreams"]) == 9 || error(
                "aerosol scene does not use nine weighted streams: $path")
        else
            Int(ds.attrib["o2_nstreams"]) == 5 || error(
                "clear scene does not use five weighted streams: $path")
        end
        o2 = Tuple(read_variable(ds, name) for name in O2_VARIABLES)
        for name in CO2_VARIABLES
            read_variable(ds, name)
        end
        return o2
    end
end

function validate_xco2_invariance(indices, aerosol)
    # State ordering is surface -> aerosol -> SIF -> XCO2. Consecutive blocks
    # of four therefore differ only in XCO2, which the A-band does not use.
    for first_index in indices[1:4:end]
        reference = validate_scene(first_index; aerosol)
        for index in first_index + 1:first_index + 3
            candidate = validate_scene(index; aerosol)
            for (name, reference_values, candidate_values) in
                    zip(O2_VARIABLES, reference, candidate)
                candidate_values == reference_values || error(
                    "$name is not bit-identical across XCO2-only states " *
                    "$first_index and $index")
            end
        end
    end
end

function validate_sif_independent_co2(indices, aerosol)
    # Within each surface/aerosol block, states 1:4 and 5:8 differ only in SIF.
    # Neither CO2 band carries SIF in this experiment.
    for block_start in indices[1:8:end]
        for offset in 0:3
            off_path = scene_path(block_start + offset, aerosol)
            on_path = scene_path(block_start + offset + 4, aerosol)
            NCDataset(off_path) do off
                NCDataset(on_path) do on
                    for name in CO2_VARIABLES
                        read_variable(off, name) == read_variable(on, name) || error(
                            "$name changed between SIF-only states " *
                            "$(block_start + offset) and $(block_start + offset + 4)")
                    end
                end
            end
        end
    end
end

function validate_wavelengths(directory)
    path = joinpath(directory, "sim_wavelength.nc")
    isfile(path) || error("missing $path")
    NCDataset(path) do ds
        nu = Float32.(ds["o2a_wavenumber"][:])
        length(nu) == 2735 || error("unexpected A-band grid length in $path")
        all(diff(nu) .> 0) || error("non-monotonic A-band grid in $path")
        step_error = maximum(abs.(diff(Float64.(nu)) .- 0.1))
        step_error <= 1e-3 || error("A-band grid is not the intended 0.1 cm-1 grid")
    end
end

function main()
    validate_wavelengths(TRUTH_ROOT)
    validate_wavelengths(AEROSOL_ROOT)
    validate_xco2_invariance(NO_AEROSOL_STATES, false)
    validate_xco2_invariance(AEROSOL_STATES, true)
    validate_sif_independent_co2(NO_AEROSOL_STATES, false)
    validate_sif_independent_co2(AEROSOL_STATES, true)
    println("clear scenes: 32; aerosol scenes: 32")
    println("exact-grid O2, completed CO2, XCO2 invariance, and CO2/SIF invariance: PASS")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
