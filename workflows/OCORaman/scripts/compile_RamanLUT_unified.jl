#!/usr/bin/env julia

# Compile the psurf-level nadir Raman LUT files and the off-nadir VZA/RAz
# chunk files into one coordinate-complete NetCDF container.
#
# Output grid:
#   psurf x sif x albedo x sza x vza x raz x stokes x wn
#
# The vza=0 slice is taken from the trusted psurf-level nadir files and
# replicated over all RAz values.  Off-nadir slices are taken from the
# VZA-expanded chunks.  Any accidental vza=0 chunk is ignored.

using Dates
using NCDatasets
using Printf

const DATA_ROOT = get(ENV, "DATA_ROOT", "/home/sanghavi/data/RamanSIFgrid")
const OUT_NC = get(ENV, "OUT_NC", joinpath(DATA_ROOT, "o2a_raman_lut_unified.nc"))
const REPORT_TXT = get(ENV, "REPORT_TXT", replace(OUT_NC, r"\.nc$" => "_compile_report.txt"))
const DRY_RUN = lowercase(get(ENV, "DRY_RUN", "0")) in ("1", "true", "yes", "y", "on")
const OVERWRITE = lowercase(get(ENV, "OVERWRITE", "0")) in ("1", "true", "yes", "y", "on")
const FILL = Float32(9.96921e36)
const AXIS_ATOL = parse(Float64, get(ENV, "AXIS_ATOL", "5e-4"))
const WN_ATOL = parse(Float64, get(ENV, "WN_ATOL", "5e-3"))
const PSURFS = Float64[500, 750, 1000]
const STOKES_VARS = ("stokes_rayleigh", "stokes_cabannes", "stokes_rrs")
const OPTICAL_VARS = ("tau_rayl", "tau_abs", "tau_aer", "varpi_cabannes", "effective_temperature")

function psurf_tag(psurf::Real)
    @sprintf("%04.0f", Float64(psurf))
end

function psurf_plain(psurf::Real)
    @sprintf("%.0f", Float64(psurf))
end

function nadir_file(psurf::Real)
    plain = psurf_plain(psurf)
    tag = psurf_tag(psurf)
    candidates = [
        joinpath(DATA_ROOT, "o2a_raman_lut_psurf$(plain).nc"),
        joinpath(DATA_ROOT, "o2a_raman_lut_psurf$(tag).nc"),
    ]
    for path in candidates
        isfile(path) && return path
    end
    error("Missing nadir psurf-level LUT for psurf=$psurf; checked $(join(candidates, ", "))")
end

function vza_root(psurf::Real)
    tag = psurf_tag(psurf)
    candidates = [
        joinpath(DATA_ROOT, "o2a_raman_lut_vza_psurf$(tag)"),
        joinpath(DATA_ROOT, "o2a_raman_lut_vza_psurf$(psurf_plain(psurf))"),
    ]
    for path in candidates
        isdir(path) && return path
    end
    error("Missing VZA chunk root for psurf=$psurf; checked $(join(candidates, ", "))")
end

function nc_files_under(root::AbstractString)
    files = String[]
    for (dir, _, names) in walkdir(root)
        for name in names
            endswith(name, ".nc") || continue
            push!(files, joinpath(dir, name))
        end
    end
    return sort(files)
end

function read_axis(path, name)
    NCDataset(path, "r") do ds
        haskey(ds, name) || error("$path lacks axis $name")
        return Float64.(ds[name][:])
    end
end

function read_int_axis(ds, name)
    haskey(ds, name) || error("Dataset lacks integer axis $name")
    return Int.(ds[name][:])
end

function same_axis(a, b; atol = AXIS_ATOL)
    length(a) == length(b) && all(isapprox.(Float64.(a), Float64.(b); atol, rtol = 0))
end

function physical_array(vals)
    any(x -> x === missing, vals) && return false
    a = Float64.(vals)
    isempty(a) && return false
    all(isfinite, a) || return false
    maximum(abs.(a)) < Float64(FILL) / 10 || return false
    maximum(abs.(a)) > 0 || return false
    return true
end

function validate_static_axes(path, axes)
    NCDataset(path, "r") do ds
        same_axis(ds["wn"][:], axes.wn; atol = WN_ATOL) || error("wn axis mismatch in $path")
        same_axis(ds["albedo"][:], axes.albedo) || error("albedo axis mismatch in $path")
        haskey(ds, "sif_on") && Int.(ds["sif_on"][:]) == axes.sif_on || error("sif_on mismatch in $path")
        if haskey(ds, "F0")
            same_axis(ds["F0"][:], axes.F0; atol = AXIS_ATOL) || error("F0 mismatch in $path")
        end
    end
    return nothing
end

function collect_axes()
    base_nadir = nadir_file(PSURFS[1])
    base_vza_root = vza_root(PSURFS[1])
    vza_probe = first(nc_files_under(base_vza_root))
    NCDataset(base_nadir, "r") do nd
        NCDataset(vza_probe, "r") do vz
            return (;
                psurf = Float32.(PSURFS),
                sif_on = Int32.(nd["sif_on"][:]),
                albedo = Float32.(nd["albedo"][:]),
                sza = Float32.(nd["sza"][:]),
                mu0 = Float32.(cosd.(Float64.(nd["sza"][:]))),
                vza = Float32.(vz["full_vza"][:]),
                raz = Float32.(vz["raz"][:]),
                stokes_index = Int32.(nd["stokes_index"][:]),
                wn = Float32.(nd["wn"][:]),
                wavelength_nm = Float32.(nd["wavelength_nm"][:]),
                band_index = Int32.(nd["band_index"][:]),
                F0 = Float32.(nd["F0"][:]),
                SIF0 = Float32.(nd["SIF0"][:, :]),
                band_solve_start = Float32.(nd["band_solve_start"][:]),
                band_solve_stop = Float32.(nd["band_solve_stop"][:]),
                band_output_start = Float32.(nd["band_output_start"][:]),
                band_output_stop = Float32.(nd["band_output_stop"][:]),
            )
        end
    end
end

function collect_vza_sources(axes)
    sources = Dict{Tuple{Int, Int, Int}, String}()
    ignored_vza1 = String[]
    for psurf in PSURFS
        ip = findfirst(x -> isapprox(x, psurf; atol = AXIS_ATOL, rtol = 0), Float64.(axes.psurf))
        ip === nothing && error("Internal psurf axis missing $psurf")
        for path in nc_files_under(vza_root(psurf))
            NCDataset(path, "r") do ds
                validate_static_axes(path, axes)
                vinds = read_int_axis(ds, "vza_full_index")
                sinds = read_int_axis(ds, "sza_full_index")
                albinds = read_int_axis(ds, "albedo_full_index")
                albinds == collect(1:length(axes.albedo)) ||
                    error("Expected all albedos in VZA chunk $path, got $albinds")
                length(sinds) == 1 || error("Expected one SZA per VZA chunk $path, got $sinds")
                for vi in vinds
                    if vi == 1
                        push!(ignored_vza1, path)
                        continue
                    end
                    key = (ip, first(sinds), vi)
                    if haskey(sources, key)
                        error("Duplicate VZA chunk node $key:\n  $(sources[key])\n  $path")
                    end
                    sources[key] = path
                end
            end
        end
    end
    return sources, unique(ignored_vza1)
end

function validate_coverage(sources, axes)
    missing = Tuple{Int, Int, Int}[]
    for ip in eachindex(axes.psurf), isza in eachindex(axes.sza), ivza in 2:length(axes.vza)
        haskey(sources, (ip, isza, ivza)) || push!(missing, (ip, isza, ivza))
    end
    isempty(missing) || error("Missing $(length(missing)) off-nadir nodes, first few: $(missing[1:min(end, 10)])")
    return nothing
end

function define_unified_file(path, axes)
    isfile(path) && begin
        OVERWRITE || error("Output exists: $path. Set OVERWRITE=1 to replace.")
        rm(path)
    end
    mkpath(dirname(path))
    ds = NCDataset(path, "c")
    defDim(ds, "psurf", length(axes.psurf))
    defDim(ds, "sif", length(axes.sif_on))
    defDim(ds, "albedo", length(axes.albedo))
    defDim(ds, "sza", length(axes.sza))
    defDim(ds, "vza", length(axes.vza))
    defDim(ds, "raz", length(axes.raz))
    defDim(ds, "stokes", length(axes.stokes_index))
    defDim(ds, "wn", length(axes.wn))
    defDim(ds, "band", length(axes.band_index |> unique))
    defDim(ds, "profile_level", 50)
    defDim(ds, "profile_layer", 49)

    defVar(ds, "psurf", axes.psurf, ("psurf",); attrib = ["units" => "hPa"])
    defVar(ds, "sif_on", axes.sif_on, ("sif",); attrib = ["description" => "0 = SurfaceSIF off, 1 = SurfaceSIF on"])
    defVar(ds, "albedo", axes.albedo, ("albedo",))
    defVar(ds, "sza", axes.sza, ("sza",); attrib = ["units" => "degree"])
    defVar(ds, "mu0", axes.mu0, ("sza",); attrib = ["description" => "cos(sza), equally spaced by construction"])
    defVar(ds, "vza", axes.vza, ("vza",); attrib = ["units" => "degree"])
    defVar(ds, "raz", axes.raz, ("raz",); attrib = ["units" => "degree", "description" => "Unique relative azimuth grid; 360 deg is omitted"])
    defVar(ds, "stokes_index", axes.stokes_index, ("stokes",); attrib = ["description" => "1=I, 2=Q, 3=U"])
    defVar(ds, "wn", axes.wn, ("wn",); attrib = ["units" => "cm-1"])
    defVar(ds, "wavelength_nm", axes.wavelength_nm, ("wn",); attrib = ["units" => "nm"])
    defVar(ds, "band_index", axes.band_index, ("wn",); attrib = ["description" => "1-based solve band that produced this output wavenumber"])
    defVar(ds, "F0", axes.F0, ("wn",); attrib = ["description" => "Solar source flux used in SolarBeam, Stokes I component"])
    defVar(ds, "SIF0", axes.SIF0, ("sif", "wn"); attrib = ["description" => "SurfaceSIF source term, Stokes I component"])
    defVar(ds, "band_solve_start", axes.band_solve_start, ("band",); attrib = ["units" => "cm-1"])
    defVar(ds, "band_solve_stop", axes.band_solve_stop, ("band",); attrib = ["units" => "cm-1"])
    defVar(ds, "band_output_start", axes.band_output_start, ("band",); attrib = ["units" => "cm-1"])
    defVar(ds, "band_output_stop", axes.band_output_stop, ("band",); attrib = ["units" => "cm-1"])

    defVar(ds, "profile_p", Float32, ("psurf", "profile_level"); fillvalue = FILL, nofill = true, attrib = ["units" => "hPa"])
    defVar(ds, "profile_T", Float32, ("psurf", "profile_layer"); fillvalue = FILL, nofill = true, attrib = ["units" => "K"])
    defVar(ds, "profile_q", Float32, ("psurf", "profile_layer"); fillvalue = FILL, nofill = true, attrib = ["units" => "kg kg-1"])
    defVar(ds, "profile_nlevel", Int32, ("psurf",); nofill = true)
    defVar(ds, "profile_nlayer", Int32, ("psurf",); nofill = true)

    dims = ("psurf", "sif", "albedo", "sza", "vza", "raz", "stokes", "wn")
    chunk = (1, 1, length(axes.albedo), 1, 1, length(axes.raz), length(axes.stokes_index), 1024)
    for name in STOKES_VARS
        defVar(ds, name, Float32, dims; chunksizes = chunk, shuffle = false,
               deflatelevel = 0, fillvalue = FILL, nofill = true,
               attrib = ["description" => "$(name) unified over psurf/sza/vza/raz"])
    end

    defVar(ds, "tau_rayl", Float32, ("psurf", "wn"); fillvalue = FILL, nofill = true,
           attrib = ["description" => "Column Rayleigh optical depth, summed over layers"])
    defVar(ds, "tau_abs", Float32, ("psurf", "wn"); fillvalue = FILL, nofill = true,
           attrib = ["description" => "Column gas absorption optical depth, summed over layers"])
    defVar(ds, "tau_aer", Float32, ("psurf", "wn"); fillvalue = FILL, nofill = true,
           attrib = ["description" => "Column aerosol optical depth; zero in this no-aerosol LUT"])
    defVar(ds, "varpi_cabannes", Float32, ("psurf", "band"); fillvalue = FILL, nofill = true,
           attrib = ["description" => "Band-center Cabannes single-scattering albedo from model"])
    defVar(ds, "effective_temperature", Float32, ("psurf",); fillvalue = FILL, nofill = true,
           attrib = ["description" => "Dry-column weighted temperature after profile reduction"])
    defVar(ds, "scene_complete", Int8, ("psurf", "sza", "vza"); fillvalue = Int8(-1),
           attrib = ["description" => "1 = complete; vza=0 replicated from nadir file; off-nadir from VZA chunks"])
    ds["scene_complete"][:, :, :] = zeros(Int8, length(axes.psurf), length(axes.sza), length(axes.vza))

    ds.attrib["title"] = "Unified O2 A-band Rayleigh/Cabannes/RRS Float32 GPU LUT"
    ds.attrib["source"] = "vSmartMOM.jl"
    ds.attrib["created_utc"] = string(now(UTC))
    ds.attrib["script"] = abspath(@__FILE__)
    ds.attrib["nadir_vza0_note"] = "vza=0 slice is replicated over all relative azimuths from the psurf-level nadir LUT files"
    ds.attrib["compile_report"] = REPORT_TXT
    return ds
end

function copy_static_psurf!(dst, axes)
    for (ip, psurf) in enumerate(PSURFS)
        path = nadir_file(psurf)
        NCDataset(path, "r") do src
            validate_static_axes(path, axes)
            dst["profile_p"][ip:ip, :] = Float32.(src["profile_p"][:, :])
            dst["profile_T"][ip:ip, :] = Float32.(src["profile_T"][:, :])
            dst["profile_q"][ip:ip, :] = Float32.(src["profile_q"][:, :])
            dst["profile_nlevel"][ip] = Int32(src["profile_nlevel"][1])
            dst["profile_nlayer"][ip] = Int32(src["profile_nlayer"][1])
            for name in OPTICAL_VARS
                haskey(src, name) || continue
                if length(size(src[name])) == 2
                    dst[name][ip:ip, :] = Float32.(src[name][:, :])
                else
                    dst[name][ip] = Float32(src[name][1])
                end
            end
        end
    end
end

function copy_nadir!(dst, axes)
    for (ip, psurf) in enumerate(PSURFS)
        path = nadir_file(psurf)
        @info "Copying nadir psurf slice" psurf path
        NCDataset(path, "r") do src
            for isza in eachindex(axes.sza)
                for name in STOKES_VARS
                    vals = Float32.(src[name][1:1, 1:1, :, isza:isza, :, :])
                    physical_array(vals) || error("Incomplete nadir source $name psurf=$psurf sza_index=$isza")
                    out = Array{Float32}(undef, 1, 1, length(axes.albedo), 1, 1,
                                         length(axes.raz), length(axes.stokes_index), length(axes.wn))
                    for iraz in eachindex(axes.raz)
                        out[:, :, :, :, 1, iraz, :, :] = vals
                    end
                    dst[name][ip:ip, 1:1, :, isza:isza, 1:1, :, :, :] = out
                end
                dst["scene_complete"][ip, isza, 1] = Int8(1)
            end
        end
        NCDatasets.sync(dst)
    end
end

function copy_offnadir!(dst, axes, sources)
    copied = 0
    for (key, path) in sort(collect(sources); by = x -> x[1])
        ip, isza, ivza_first_key = key
        NCDataset(path, "r") do src
            vinds = read_int_axis(src, "vza_full_index")
            vinds = filter(!=(1), vinds)
            sinds = read_int_axis(src, "sza_full_index")
            length(sinds) == 1 || error("Expected single SZA in $path")
            isza_src = 1
            contiguous = vinds == collect(first(vinds):last(vinds))
            contiguous || error("Non-contiguous VZA indices $vinds in $path")
            vdst = first(vinds):last(vinds)
            sza_idx = first(sinds)
            # Copy each source file only when visiting its first key.
            (isza == sza_idx && ivza_first_key == first(vinds)) || return
            @info "Copying off-nadir chunk" psurf = axes.psurf[ip] sza_index = sza_idx vza_indices = vinds path
            for name in STOKES_VARS
                vals = Float32.(src[name][:, :, :, isza_src:isza_src, :, :, :, :])
                physical_array(vals) || error("Incomplete off-nadir source $name in $path")
                dst[name][ip:ip, 1:1, :, sza_idx:sza_idx, vdst, :, :, :] = vals
            end
            for ivza in vinds
                dst["scene_complete"][ip, sza_idx, ivza] = Int8(1)
            end
            copied += 1
        end
        copied % 10 == 0 && NCDatasets.sync(dst)
    end
    NCDatasets.sync(dst)
    return copied
end

function validate_output(path, axes)
    NCDataset(path, "r") do ds
        mask = Int.(ds["scene_complete"][:, :, :])
        expected = length(axes.psurf) * length(axes.sza) * length(axes.vza)
        complete = count(==(1), mask)
        complete == expected ||
            error("Unified LUT scene_complete has $complete/$expected complete scenes")
        for name in STOKES_VARS
            @info "Validating unified variable" name
            for ip in eachindex(axes.psurf), isza in eachindex(axes.sza), ivza in eachindex(axes.vza)
                vals = ds[name][ip:ip, 1:1, :, isza:isza, ivza:ivza, :, :, :]
                physical_array(vals) ||
                    error("Unified output incomplete: $name psurf_index=$ip sza_index=$isza vza_index=$ivza")
            end
        end
    end
end

function write_report(axes, sources, ignored_vza1, copied_chunks; status)
    mkpath(dirname(REPORT_TXT))
    open(REPORT_TXT, "w") do io
        println(io, "Unified Raman LUT compilation report")
        println(io, "created_utc: ", string(now(UTC)))
        println(io, "status: ", status)
        println(io, "output_nc: ", OUT_NC)
        println(io, "data_root: ", DATA_ROOT)
        println(io, "psurf: ", join(axes.psurf, ","))
        println(io, "albedo_count: ", length(axes.albedo))
        println(io, "sza_count: ", length(axes.sza))
        println(io, "vza_count: ", length(axes.vza))
        println(io, "raz_count: ", length(axes.raz))
        println(io, "wn_count: ", length(axes.wn))
        println(io, "nadir_files:")
        for psurf in PSURFS
            println(io, "  ", psurf, " hPa: ", nadir_file(psurf))
        end
        println(io, "offnadir_nodes: ", length(sources))
        println(io, "offnadir_chunks_copied: ", copied_chunks)
        println(io, "ignored_vza0_chunks: ", length(ignored_vza1))
        for path in ignored_vza1
            println(io, "  ignored: ", path)
        end
        println(io, "note: vza=0 was replicated across all RAz values from the psurf-level nadir files.")
    end
end

function main()
    println("Unified Raman LUT compilation")
    println("  DATA_ROOT: $DATA_ROOT")
    println("  OUT_NC: $OUT_NC")
    println("  DRY_RUN: $DRY_RUN")
    println("  OVERWRITE: $OVERWRITE")

    axes = collect_axes()
    for psurf in PSURFS
        validate_static_axes(nadir_file(psurf), axes)
    end
    sources, ignored_vza1 = collect_vza_sources(axes)
    validate_coverage(sources, axes)

    expected_scenes = length(axes.psurf) * length(axes.sza) * length(axes.vza)
    expected_offnadir = length(axes.psurf) * length(axes.sza) * (length(axes.vza) - 1)
    println("  expected scenes: $expected_scenes")
    println("  expected off-nadir nodes: $expected_offnadir")
    println("  available off-nadir nodes: $(length(sources))")
    println("  ignored vza=0 chunk files: $(length(ignored_vza1))")

    if DRY_RUN
        write_report(axes, sources, ignored_vza1, 0; status = "dry_run_ok")
        println("DRY_RUN complete; report: $REPORT_TXT")
        return
    end

    ds = define_unified_file(OUT_NC, axes)
    copied_chunks = 0
    try
        copy_static_psurf!(ds, axes)
        copy_nadir!(ds, axes)
        copied_chunks = copy_offnadir!(ds, axes, sources)
        ds.attrib["completed_utc"] = string(now(UTC))
        NCDatasets.sync(ds)
    finally
        close(ds)
    end

    validate_output(OUT_NC, axes)
    write_report(axes, sources, ignored_vza1, copied_chunks; status = "complete")
    println("Unified LUT complete: $OUT_NC")
    println("Report: $REPORT_TXT")
end

main()
