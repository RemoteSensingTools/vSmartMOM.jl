#!/usr/bin/env julia

# Compile Raman LUT components into psurf/SZA shard files instead of one huge
# unified container.
#
# Each output file contains one psurf and one SZA:
#   psurf x sif x albedo x sza x vza x raz x stokes x wn
#
# The vza=0 slice comes from the psurf-level nadir LUT and is replicated over
# all RAz values.  Off-nadir slices come from the VZA/RAz chunk files.

using Dates
using NCDatasets
using Printf

const DATA_ROOT = get(ENV, "DATA_ROOT", "/home/sanghavi/data/RamanSIFgrid")
const OUTDIR = get(ENV, "OUTDIR", joinpath(DATA_ROOT, "o2a_raman_lut_by_psurf_sza"))
const DRY_RUN = lowercase(get(ENV, "DRY_RUN", "0")) in ("1", "true", "yes", "y", "on")
const OVERWRITE = lowercase(get(ENV, "OVERWRITE", "0")) in ("1", "true", "yes", "y", "on")
const FILL = Float32(9.96921e36)
const AXIS_ATOL = parse(Float64, get(ENV, "AXIS_ATOL", "5e-4"))
const WN_ATOL = parse(Float64, get(ENV, "WN_ATOL", "5e-3"))
const STOKES_VARS = ("stokes_rayleigh", "stokes_cabannes", "stokes_rrs")
const OPTICAL_VARS = ("tau_rayl", "tau_abs", "tau_aer", "varpi_cabannes", "effective_temperature")

function parse_list(::Type{T}, text::AbstractString) where {T}
    vals = T[]
    for raw in split(strip(text), ',')
        part = strip(raw)
        isempty(part) && continue
        push!(vals, parse(T, part))
    end
    return vals
end

const PSURFS = parse_list(Float64, get(ENV, "PSURFS", "500,750,1000"))
const SZA_IDXS = parse_list(Int, get(ENV, "SZA_IDXS", join(1:14, ",")))

psurf_tag(psurf::Real) = @sprintf("%04.0f", Float64(psurf))
psurf_plain(psurf::Real) = @sprintf("%.0f", Float64(psurf))

function nadir_file(psurf::Real)
    plain = psurf_plain(psurf)
    tag = psurf_tag(psurf)
    for path in (joinpath(DATA_ROOT, "o2a_raman_lut_psurf$(plain).nc"),
                 joinpath(DATA_ROOT, "o2a_raman_lut_psurf$(tag).nc"))
        isfile(path) && return path
    end
    error("Missing nadir psurf-level LUT for psurf=$psurf")
end

function vza_root(psurf::Real)
    tag = psurf_tag(psurf)
    plain = psurf_plain(psurf)
    for path in (joinpath(DATA_ROOT, "o2a_raman_lut_vza_psurf$(tag)"),
                 joinpath(DATA_ROOT, "o2a_raman_lut_vza_psurf$(plain)"))
        isdir(path) && return path
    end
    error("Missing VZA chunk root for psurf=$psurf")
end

function nc_files_under(root::AbstractString)
    out = String[]
    for (dir, _, names) in walkdir(root)
        for name in names
            endswith(name, ".nc") && push!(out, joinpath(dir, name))
        end
    end
    return sort(out)
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

function collect_axes()
    base_nadir = nadir_file(first(PSURFS))
    probe = first(nc_files_under(vza_root(first(PSURFS))))
    NCDataset(base_nadir, "r") do nd
        NCDataset(probe, "r") do vz
            return (;
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

function validate_static_axes(path, axes)
    NCDataset(path, "r") do ds
        same_axis(ds["wn"][:], axes.wn; atol = WN_ATOL) || error("wn mismatch in $path")
        same_axis(ds["albedo"][:], axes.albedo) || error("albedo mismatch in $path")
        Int.(ds["sif_on"][:]) == Int.(axes.sif_on) || error("sif_on mismatch in $path")
        haskey(ds, "F0") && same_axis(ds["F0"][:], axes.F0; atol = AXIS_ATOL) || error("F0 mismatch in $path")
    end
end

function find_vza_chunks(psurf, sza_idx)
    files = String[]
    ignored = String[]
    for path in nc_files_under(vza_root(psurf))
        NCDataset(path, "r") do ds
            haskey(ds, "sza_full_index") || return
            sza_inds = Int.(ds["sza_full_index"][:])
            sza_idx in sza_inds || return
            vza_inds = Int.(ds["vza_full_index"][:])
            if all(==(1), vza_inds)
                push!(ignored, path)
            else
                push!(files, path)
            end
        end
    end
    # Sort by first VZA index so writes are deterministic.
    sort!(files; by = path -> NCDataset(path, "r") do ds
        minimum(filter(!=(1), Int.(ds["vza_full_index"][:])))
    end)
    return files, ignored
end

function coverage_ok(psurf, axes)
    missing = Tuple{Int, Int}[]
    ignored = String[]
    for sza_idx in SZA_IDXS
        files, ign = find_vza_chunks(psurf, sza_idx)
        append!(ignored, ign)
        got = Set{Int}()
        for path in files
            NCDataset(path, "r") do ds
                union!(got, filter(!=(1), Int.(ds["vza_full_index"][:])))
            end
        end
        for ivza in 2:length(axes.vza)
            ivza in got || push!(missing, (sza_idx, ivza))
        end
    end
    return missing, unique(ignored)
end

function shard_path(psurf, sza_idx)
    dir = joinpath(OUTDIR, "psurf$(psurf_plain(psurf))")
    return joinpath(dir, @sprintf("o2a_raman_lut_psurf%s_sza%03d.nc", psurf_plain(psurf), sza_idx))
end

function define_shard(path, axes, psurf, sza_idx)
    isfile(path) && begin
        OVERWRITE || error("Output exists: $path. Set OVERWRITE=1 to replace.")
        rm(path)
    end
    mkpath(dirname(path))
    ds = NCDataset(path, "c")
    defDim(ds, "psurf", 1)
    defDim(ds, "sif", length(axes.sif_on))
    defDim(ds, "albedo", length(axes.albedo))
    defDim(ds, "sza", 1)
    defDim(ds, "vza", length(axes.vza))
    defDim(ds, "raz", length(axes.raz))
    defDim(ds, "stokes", length(axes.stokes_index))
    defDim(ds, "wn", length(axes.wn))
    defDim(ds, "band", length(unique(axes.band_index)))
    defDim(ds, "full_sza", length(axes.sza))
    defDim(ds, "profile_level", 50)
    defDim(ds, "profile_layer", 49)

    defVar(ds, "psurf", Float32[psurf], ("psurf",); attrib = ["units" => "hPa"])
    defVar(ds, "sif_on", axes.sif_on, ("sif",); attrib = ["description" => "0 = SurfaceSIF off, 1 = SurfaceSIF on"])
    defVar(ds, "albedo", axes.albedo, ("albedo",))
    defVar(ds, "sza", Float32[axes.sza[sza_idx]], ("sza",); attrib = ["units" => "degree"])
    defVar(ds, "sza_full_index", Int32[sza_idx], ("sza",); attrib = ["description" => "1-based index in full SZA grid"])
    defVar(ds, "full_sza", axes.sza, ("full_sza",); attrib = ["units" => "degree", "description" => "Full intended SZA axis"])
    defVar(ds, "mu0", Float32[axes.mu0[sza_idx]], ("sza",); attrib = ["description" => "cos(sza)"])
    defVar(ds, "vza", axes.vza, ("vza",); attrib = ["units" => "degree"])
    defVar(ds, "raz", axes.raz, ("raz",); attrib = ["units" => "degree"])
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
        defVar(ds, name, Float32, dims; chunksizes = chunk, fillvalue = FILL, nofill = true,
               attrib = ["description" => "$(name), psurf/SZA shard"])
    end
    defVar(ds, "tau_rayl", Float32, ("psurf", "wn"); fillvalue = FILL, nofill = true)
    defVar(ds, "tau_abs", Float32, ("psurf", "wn"); fillvalue = FILL, nofill = true)
    defVar(ds, "tau_aer", Float32, ("psurf", "wn"); fillvalue = FILL, nofill = true)
    defVar(ds, "varpi_cabannes", Float32, ("psurf", "band"); fillvalue = FILL, nofill = true)
    defVar(ds, "effective_temperature", Float32, ("psurf",); fillvalue = FILL, nofill = true)
    defVar(ds, "scene_complete", Int8, ("vza",); fillvalue = Int8(-1),
           attrib = ["description" => "1 = complete for this psurf/SZA shard"])
    ds["scene_complete"][:] = zeros(Int8, length(axes.vza))

    ds.attrib["title"] = "O2 A-band Raman LUT psurf/SZA shard"
    ds.attrib["source"] = "vSmartMOM.jl"
    ds.attrib["created_utc"] = string(now(UTC))
    ds.attrib["script"] = abspath(@__FILE__)
    ds.attrib["nadir_vza0_note"] = "vza=0 slice is replicated over all RAz values from the psurf-level nadir LUT file"
    return ds
end

function copy_static!(dst, axes, psurf)
    NCDataset(nadir_file(psurf), "r") do src
        dst["profile_p"][1:1, :] = Float32.(src["profile_p"][:, :])
        dst["profile_T"][1:1, :] = Float32.(src["profile_T"][:, :])
        dst["profile_q"][1:1, :] = Float32.(src["profile_q"][:, :])
        dst["profile_nlevel"][1] = Int32(src["profile_nlevel"][1])
        dst["profile_nlayer"][1] = Int32(src["profile_nlayer"][1])
        for name in OPTICAL_VARS
            haskey(src, name) || continue
            if length(size(src[name])) == 2
                dst[name][1:1, :] = Float32.(src[name][:, :])
            else
                dst[name][1] = Float32(src[name][1])
            end
        end
    end
end

function copy_nadir!(dst, axes, psurf, sza_idx)
    NCDataset(nadir_file(psurf), "r") do src
        for name in STOKES_VARS
            vals = Float32.(src[name][1:1, 1:1, :, sza_idx:sza_idx, :, :])
            physical_array(vals) || error("Incomplete nadir source $name psurf=$psurf sza_idx=$sza_idx")
            out = Array{Float32}(undef, 1, 1, length(axes.albedo), 1, 1,
                                 length(axes.raz), length(axes.stokes_index), length(axes.wn))
            for iraz in eachindex(axes.raz)
                out[:, :, :, :, 1, iraz, :, :] = vals
            end
            dst[name][1:1, 1:1, :, 1:1, 1:1, :, :, :] = out
        end
    end
    dst["scene_complete"][1] = Int8(1)
end

function copy_offnadir!(dst, axes, psurf, sza_idx)
    files, ignored = find_vza_chunks(psurf, sza_idx)
    got = Set{Int}()
    for path in files
        NCDataset(path, "r") do src
            validate_static_axes(path, axes)
            vinds = filter(!=(1), Int.(src["vza_full_index"][:]))
            sza_inds = Int.(src["sza_full_index"][:])
            sza_idx in sza_inds || return
            src_sza_pos = findfirst(==(sza_idx), sza_inds)
            vinds == collect(first(vinds):last(vinds)) || error("Non-contiguous VZA indices $vinds in $path")
            for name in STOKES_VARS
                vals = Float32.(src[name][:, :, :, src_sza_pos:src_sza_pos, :, :, :, :])
                physical_array(vals) || error("Incomplete off-nadir source $name in $path")
                dst[name][1:1, 1:1, :, 1:1, first(vinds):last(vinds), :, :, :] = vals
            end
            for ivza in vinds
                push!(got, ivza)
                dst["scene_complete"][ivza] = Int8(1)
            end
        end
    end
    missing = setdiff(Set(2:length(axes.vza)), got)
    isempty(missing) || error("Missing VZA indices for psurf=$psurf sza_idx=$sza_idx: $(sort(collect(missing)))")
    return length(files), ignored
end

function validate_shard(path, axes)
    NCDataset(path, "r") do ds
        all(Int.(ds["scene_complete"][:]) .== 1) || error("Incomplete scene_complete in $path")
        for name in STOKES_VARS
            for ivza in eachindex(axes.vza)
                vals = ds[name][1:1, 1:1, :, 1:1, ivza:ivza, :, :, :]
                physical_array(vals) || error("Incomplete output $name vza_idx=$ivza in $path")
            end
        end
    end
end

function main()
    axes = collect_axes()
    println("Raman LUT psurf/SZA shard compilation")
    println("  DATA_ROOT: $DATA_ROOT")
    println("  OUTDIR: $OUTDIR")
    println("  DRY_RUN: $DRY_RUN")
    println("  OVERWRITE: $OVERWRITE")
    println("  PSURFS: $(join(PSURFS, ","))")
    println("  SZA_IDXS: $(join(SZA_IDXS, ","))")
    println("  shard count: $(length(PSURFS) * length(SZA_IDXS))")

    total_ignored = Set{String}()
    for psurf in PSURFS
        validate_static_axes(nadir_file(psurf), axes)
        missing, ignored = coverage_ok(psurf, axes)
        union!(total_ignored, ignored)
        isempty(missing) || error("Missing VZA coverage for psurf=$psurf: $(missing[1:min(end, 10)])")
    end
    println("  coverage: complete")
    println("  ignored vza=0 VZA chunk files: $(length(total_ignored))")

    if DRY_RUN
        mkpath(OUTDIR)
        report = joinpath(OUTDIR, "compile_psurf_sza_shards_dryrun.txt")
        open(report, "w") do io
            println(io, "dry_run_utc: ", string(now(UTC)))
            println(io, "psurfs: ", join(PSURFS, ","))
            println(io, "sza_idxs: ", join(SZA_IDXS, ","))
            println(io, "shard_count: ", length(PSURFS) * length(SZA_IDXS))
            println(io, "ignored_vza0_chunks: ", length(total_ignored))
            for path in sort(collect(total_ignored))
                println(io, "ignored: ", path)
            end
        end
        println("DRY_RUN complete: $report")
        return
    end

    for psurf in PSURFS, sza_idx in SZA_IDXS
        out = shard_path(psurf, sza_idx)
        @info "Writing shard" psurf sza_idx out
        ds = define_shard(out, axes, psurf, sza_idx)
        copied = 0
        try
            copy_static!(ds, axes, psurf)
            copy_nadir!(ds, axes, psurf, sza_idx)
            copied, _ = copy_offnadir!(ds, axes, psurf, sza_idx)
            ds.attrib["completed_utc"] = string(now(UTC))
            ds.attrib["offnadir_chunk_files_copied"] = copied
            NCDatasets.sync(ds)
        finally
            close(ds)
        end
        validate_shard(out, axes)
    end
    println("Shard compilation complete: $OUTDIR")
end

main()
