#!/usr/bin/env julia

"""Export the old and corrected forest RRS spectra around O2 chunks 42/43."""

using DelimitedFiles
using JLD2
using NCDatasets

const ROOT = normpath(joinpath(@__DIR__, ".."))
const OBSOLETE_DIAGNOSTICS = joinpath(
    ROOT, "obsolete", "truth_map", "chunk_surface_error_diagnostics")
const OLD_SCENE = joinpath(
    ROOT, "obsolete", "truth_map", "aerosol_chunked", "old_w_error",
    "hiressim_057.nc")
const BEFORE_AFTER_DAT = joinpath(
    OBSOLETE_DIAGNOSTICS, "rrs_chunk42_43_before_after.dat")
const FIXED42 = get(ENV, "FIXED_CHUNK42", "/tmp/fixed_surface_n8_chunk42.jld2")
const FIXED43 = get(ENV, "FIXED_CHUNK43", "/tmp/fixed_surface_n8_chunk43.jld2")

function main()
    mkpath(OBSOLETE_DIAGNOSTICS)
    for path in (OLD_SCENE, FIXED42, FIXED43)
        isfile(path) || error("required validation input is missing: $path")
    end

    old_rrs = NCDataset(OLD_SCENE) do ds
        Array(ds["radiance_rrs_o2a"][:, 2625:2735])
    end
    d42, d43 = load(FIXED42), load(FIXED43)
    corrected_rrs = hcat(d42["rrs"], d43["rrs"])
    ν = vcat(d42["core_ν"], d43["core_ν"])
    size(old_rrs) == size(corrected_rrs) || error("before/after grids differ")

    open(BEFORE_AFTER_DAT, "w") do io
        println(io, "# Forest aerosol scene 57, O2 chunks 42 and 43")
        println(io, "# before: production Float32 nstreams=9; after: corrected Float32 nstreams=8")
        println(io, "# chunk boundary lies between rows 64 and 65")
        println(io, "# wavelength_nm wavenumber_cm-1 rrs_I_before rrs_I_after")
        writedlm(io, hcat(1e7 ./ ν, ν, old_rrs[1, :], corrected_rrs[1, :]))
    end
    println(BEFORE_AFTER_DAT)
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
