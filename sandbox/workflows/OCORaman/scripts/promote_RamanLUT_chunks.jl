#!/usr/bin/env julia

# Promote validated chunk/subset Raman LUT scenes into the psurf-level NetCDF
# containers.  This is intentionally conservative:
#
#   * unreadable/open/corrupt chunks are skipped;
#   * only scenes with complete finite Rayleigh, Cabannes, and RRS Stokes
#     arrays are copied;
#   * already complete target scenes are left untouched unless OVERWRITE=1.
#
# Chunk files are mapped into the full LUT by albedo_full_index and
# sza_full_index when present, otherwise by coordinate values.

using Dates
using NCDatasets
using Printf

function parse_bool(text::AbstractString)
    lowercase(strip(text)) in ("1", "true", "yes", "y", "on")
end

function parse_list(::Type{T}, text::AbstractString) where {T}
    vals = T[]
    for raw in split(strip(text), ',')
        part = strip(raw)
        isempty(part) && continue
        push!(vals, parse(T, part))
    end
    return vals
end

const DATA_ROOT = get(ENV, "DATA_ROOT", "/home/sanghavi/data/RamanSIFgrid/O2ABand")
const PSURFS = parse_list(Float64, get(ENV, "PSURFS", "1000,750"))
const DRY_RUN = parse_bool(get(ENV, "DRY_RUN", "0"))
const OVERWRITE = parse_bool(get(ENV, "OVERWRITE", "0"))
const MIN_AGE_SECONDS = parse(Float64, get(ENV, "MIN_AGE_SECONDS", "60"))
const FILL_LIMIT = Float64(parse(Float32, get(ENV, "FILL_LIMIT", "1e20")))
const WN_ATOL = parse(Float64, get(ENV, "WN_ATOL", "5e-3"))
const AXIS_ATOL = parse(Float64, get(ENV, "AXIS_ATOL", "5e-4"))
const COPY_OPTICAL = parse_bool(get(ENV, "COPY_OPTICAL", "1"))
const STOKES_VARS = ("stokes_rayleigh", "stokes_cabannes", "stokes_rrs")
const OPTICAL_VARS = ("tau_rayl", "tau_abs", "tau_aer", "varpi_cabannes", "effective_temperature")

function psurf_tag(psurf::Real)
    @sprintf("%04.0f", Float64(psurf))
end

function psurf_plain(psurf::Real)
    @sprintf("%.0f", Float64(psurf))
end

function target_file(psurf::Real)
    plain = psurf_plain(psurf)
    tag = psurf_tag(psurf)
    for name in ("o2a_raman_lut_psurf$(plain).nc", "o2a_raman_lut_psurf$(tag).nc")
        path = joinpath(DATA_ROOT, name)
        isfile(path) && return path
    end
    return joinpath(DATA_ROOT, "o2a_raman_lut_psurf$(plain).nc")
end

function chunk_dirs(psurf::Real)
    tag = psurf_tag(psurf)
    plain = psurf_plain(psurf)
    dirs = [
        joinpath(DATA_ROOT, "o2a_raman_lut_priority_psurf$(tag)_rgb"),
        joinpath(DATA_ROOT, "o2a_raman_lut_priority_psurf$(plain)_rgb"),
        joinpath(DATA_ROOT, "o2a_raman_lut_chunked_psurf$(tag)"),
        joinpath(DATA_ROOT, "o2a_raman_lut_chunked_psurf$(plain)"),
    ]
    extra = strip(get(ENV, "CHUNK_DIRS", ""))
    if !isempty(extra)
        append!(dirs, split(extra, ':'))
    end
    return unique(filter(isdir, dirs))
end

function chunk_files(psurf::Real)
    tag = psurf_tag(psurf)
    plain = psurf_plain(psurf)
    files = String[]
    for dir in chunk_dirs(psurf)
        for path in sort(readdir(dir; join = true))
            isfile(path) || continue
            endswith(path, ".nc") || continue
            base = basename(path)
            (occursin("psurf$(tag)", base) || occursin("psurf$(plain)", base)) || continue
            push!(files, path)
        end
    end
    return unique(files)
end

function file_age_seconds(path::AbstractString)
    return time() - stat(path).mtime
end

function hasvar(ds, name::AbstractString)
    haskey(ds, name)
end

function axis_values(ds, name::AbstractString)
    hasvar(ds, name) || return nothing
    return Float64.(Array(ds[name][:]))
end

function full_indices(ds, full_index_name::AbstractString, axis_name::AbstractString)
    if hasvar(ds, full_index_name)
        return Int.(Array(ds[full_index_name][:]))
    end
    hasvar(ds, axis_name) || error("Dataset lacks both '$full_index_name' and '$axis_name'")
    return collect(1:length(ds[axis_name]))
end

function find_axis_position(ds, axis_name::AbstractString, value::Real; atol = AXIS_ATOL)
    vals = axis_values(ds, axis_name)
    vals === nothing && return nothing
    hits = findall(x -> isapprox(x, Float64(value); atol = atol, rtol = 0), vals)
    return isempty(hits) ? nothing : first(hits)
end

function target_axis_position(ds, full_index_name::AbstractString, axis_name::AbstractString,
                              full_index::Integer, value::Real)
    if hasvar(ds, full_index_name)
        vals = Int.(Array(ds[full_index_name][:]))
        hits = findall(==(full_index), vals)
        return isempty(hits) ? nothing : first(hits)
    end
    pos = find_axis_position(ds, axis_name, value)
    pos !== nothing && return pos
    return 1 <= full_index <= length(ds[axis_name]) ? full_index : nothing
end

function physical_array(vals)
    any(x -> x === missing, vals) && return false
    a = Float64.(vals)
    isempty(a) && return false
    all(isfinite, a) || return false
    all(abs.(a) .< FILL_LIMIT) || return false
    maximum(abs.(a)) > 0 || return false
    return true
end

function scene_complete(ds, ip::Integer, isif::Integer, ialb::Integer, isza::Integer)
    for name in STOKES_VARS
        hasvar(ds, name) || return false
        vals = Array(ds[name][ip, isif, ialb, isza, :, :])
        physical_array(vals) || return false
    end
    return true
end

function compatible_static_axes(src, dst)
    for name in ("wn", "F0", "band_index")
        hasvar(src, name) && hasvar(dst, name) || continue
        a = Array(src[name][:])
        b = Array(dst[name][:])
        length(a) == length(b) || return false, "$name length mismatch"
        if name == "band_index"
            Int.(a) == Int.(b) || return false, "$name values mismatch"
        else
            all(isapprox.(Float64.(a), Float64.(b); atol = name == "wn" ? WN_ATOL : AXIS_ATOL, rtol = 0)) ||
                return false, "$name values mismatch"
        end
    end
    return true, ""
end

function copy_optical!(dst, src, tip::Integer, sip::Integer, source_path::AbstractString, stats)
    COPY_OPTICAL || return nothing
    for name in OPTICAL_VARS
        (hasvar(src, name) && hasvar(dst, name)) || continue
        vals = try
            length(size(src[name])) == 2 ? Array(src[name][sip, :]) : Array(src[name][sip])
        catch
            continue
        end
        physical_array(vals) || continue
        if DRY_RUN
            stats[:optical_would_write] += 1
        else
            if length(size(dst[name])) == 2
                dst[name][tip, :] = vals
            else
                dst[name][tip] = vals
            end
            stats[:optical_written] += 1
        end
    end
    return nothing
end

function promote_scene!(dst, src, tip, sip, tisif, sisif, tialb, sialb, tisza, sisza, source_path, stats)
    if scene_complete(dst, tip, tisif, tialb, tisza) && !OVERWRITE
        stats[:target_complete] += 1
        return :target_complete
    end

    scene_complete(src, sip, sisif, sialb, sisza) || begin
        stats[:source_incomplete] += 1
        return :source_incomplete
    end

    if DRY_RUN
        stats[:would_write] += 1
        return :would_write
    end

    for name in STOKES_VARS
        dst[name][tip, tisif, tialb, tisza, :, :] = Array(src[name][sip, sisif, sialb, sisza, :, :])
    end
    stats[:written] += 1
    return :written
end

function annotate_target!(ds)
    DRY_RUN && return nothing
    ds.attrib["last_chunk_promotion_utc"] = string(now(UTC))
    ds.attrib["chunk_promotion_script"] = abspath(@__FILE__)
    return nothing
end

function promote_file!(dst, source_path::AbstractString, stats)
    age = file_age_seconds(source_path)
    if age < MIN_AGE_SECONDS
        @info "Skipping very recent chunk" source_path age_seconds = round(age; digits = 1)
        stats[:too_young] += 1
        return nothing
    end

    src = try
        NCDataset(source_path, "r")
    catch err
        @warn "Skipping unreadable/open/corrupt chunk" source_path error = sprint(showerror, err)
        stats[:unreadable] += 1
        return nothing
    end

    try
        ok, why = compatible_static_axes(src, dst)
        if !ok
            @warn "Skipping chunk with incompatible axes" source_path reason = why
            stats[:incompatible] += 1
            return nothing
        end

        sp = axis_values(src, "psurf")
        sp === nothing && error("Chunk lacks psurf axis")
        salbedo = axis_values(src, "albedo")
        ssza = axis_values(src, "sza")
        ssif = hasvar(src, "sif_on") ? Int.(Array(src["sif_on"][:])) : collect(1:length(src["sif"]))
        alb_full = full_indices(src, "albedo_full_index", "albedo")
        sza_full = full_indices(src, "sza_full_index", "sza")

        for sip in eachindex(sp)
            tip = find_axis_position(dst, "psurf", sp[sip])
            tip === nothing && continue
            copy_optical!(dst, src, tip, sip, source_path, stats)

            for sisif in eachindex(ssif)
                tisif = find_axis_position(dst, "sif_on", ssif[sisif]; atol = 0)
                tisif === nothing && (tisif = ssif[sisif] + 1)
                1 <= tisif <= length(dst["sif_on"]) || continue

                for sialb in eachindex(salbedo)
                    tialb = target_axis_position(dst, "albedo_full_index", "albedo",
                                                 alb_full[sialb], salbedo[sialb])
                    tialb === nothing && continue
                    for sisza in eachindex(ssza)
                        tisza = target_axis_position(dst, "sza_full_index", "sza",
                                                     sza_full[sisza], ssza[sisza])
                        tisza === nothing && continue
                        status = promote_scene!(dst, src, tip, sip, tisif, sisif,
                                                tialb, sialb, tisza, sisza,
                                                source_path, stats)
                        status in (:written, :would_write) && @info "Promoted Raman LUT scene" status source_path psurf = sp[sip] sza = ssza[sisza] albedo = salbedo[sialb]
                    end
                end
            end
        end
    finally
        close(src)
    end
    return nothing
end

function new_stats()
    return Dict{Symbol, Int}(
        :files_seen => 0,
        :missing_target => 0,
        :too_young => 0,
        :unreadable => 0,
        :incompatible => 0,
        :source_incomplete => 0,
        :target_complete => 0,
        :would_write => 0,
        :written => 0,
        :optical_would_write => 0,
        :optical_written => 0,
    )
end

function main()
    println("Raman LUT chunk promotion")
    println("  time UTC: $(now(UTC))")
    println("  DATA_ROOT: $DATA_ROOT")
    println("  PSURFS: $(join(PSURFS, ", "))")
    println("  DRY_RUN: $DRY_RUN")
    println("  OVERWRITE: $OVERWRITE")
    println("  MIN_AGE_SECONDS: $MIN_AGE_SECONDS")

    stats = new_stats()
    for psurf in PSURFS
        target = target_file(psurf)
        if !isfile(target)
            @warn "Skipping psurf with missing target NetCDF" psurf target
            stats[:missing_target] += 1
            continue
        end

        files = chunk_files(psurf)
        stats[:files_seen] += length(files)
        println("\npsurf=$(psurf_plain(psurf)) target=$target")
        println("  candidate chunks: $(length(files))")
        isempty(files) && continue

        dst = DRY_RUN ? NCDataset(target, "r") : NCDataset(target, "a")
        try
            for source_path in files
                promote_file!(dst, source_path, stats)
            end
            annotate_target!(dst)
            !DRY_RUN && NCDatasets.sync(dst)
        finally
            close(dst)
        end
    end

    println("\nPromotion summary")
    for key in sort(collect(keys(stats)); by = string)
        println("  $(key): $(stats[key])")
    end
end

main()
