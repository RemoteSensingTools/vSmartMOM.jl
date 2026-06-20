#!/usr/bin/env julia

# Print the requested full-grid albedo indices that are still missing for one
# psurf/SZA pair.  Completion is checked across:
#   - the monolithic psurf LUT, when present
#   - dense chunked LUT files, when present
#   - priority chunk files, when present
#
# The shell priority launcher consumes stdout.  Diagnostic messages go to
# stderr so command substitution remains clean.

using NCDatasets
using Printf

const FILL_LIMIT = Float32(1e20)
const DATA_ROOT = get(ENV, "DATA_ROOT", "/home/sanghavi/data/RamanSIFgrid")
const OUTDIR = get(ENV, "OUTDIR", "")

function parse_selector(text::AbstractString)
    vals = Int[]
    for raw in split(strip(text), ',')
        part = strip(raw)
        isempty(part) && continue
        if occursin(':', part)
            bits = split(part, ':')
            length(bits) == 2 || error("Only a:b selectors are supported here; got $part")
            append!(vals, parse(Int, bits[1]):parse(Int, bits[2]))
        else
            push!(vals, parse(Int, part))
        end
    end
    return unique(vals)
end

function psurf_tag(psurf::Real)
    return @sprintf("%04.0f", Float64(psurf))
end

function psurf_plain(psurf::Real)
    return @sprintf("%.0f", Float64(psurf))
end

function candidate_files(psurf::Real)
    tag = psurf_tag(psurf)
    plain = psurf_plain(psurf)
    files = String[]

    for name in ("o2a_raman_lut_psurf$(plain).nc", "o2a_raman_lut_psurf$(tag).nc")
        path = joinpath(DATA_ROOT, name)
        isfile(path) && push!(files, path)
    end

    for d in (
        joinpath(DATA_ROOT, "o2a_raman_lut_chunked_psurf$(tag)"),
        joinpath(DATA_ROOT, "o2a_raman_lut_chunked_psurf$(plain)"),
        OUTDIR,
    )
        (isempty(d) || !isdir(d)) && continue
        append!(files, sort(String.(filter(isfile, readdir(d; join = true)))))
    end

    return unique(filter(p -> endswith(p, ".nc"), files))
end

function hasvar(ds, name::AbstractString)
    return haskey(ds, name)
end

function selected_position(ds, full_index_name::AbstractString, axis_name::AbstractString, target::Integer)
    if hasvar(ds, full_index_name)
        vals = Int.(Array(ds[full_index_name][:]))
        hits = findall(==(target), vals)
        return isempty(hits) ? nothing : first(hits)
    end

    if hasvar(ds, axis_name)
        n = length(ds[axis_name])
        return 1 <= target <= n ? target : nothing
    end

    return nothing
end

function scene_complete(ds, ialb::Integer, isza::Integer)
    ip = 1
    isif = 1
    for name in ("stokes_rayleigh", "stokes_cabannes", "stokes_rrs")
        hasvar(ds, name) || return false
        vals = Array(ds[name][ip, isif, ialb, isza, :, :])
        physical = isfinite.(vals) .& (abs.(vals) .< FILL_LIMIT)
        all(physical) || return false
        maximum(abs.(vals)) > 0f0 || return false
    end
    return true
end

function complete_in_file(path::AbstractString, alb_idx::Integer, sza_idx::Integer)
    ds = try
        NCDataset(path, "r")
    catch err
        @warn "Skipping unreadable candidate LUT" path exception = (err, catch_backtrace())
        return false
    end
    try
        ialb = selected_position(ds, "albedo_full_index", "albedo", alb_idx)
        isza = selected_position(ds, "sza_full_index", "sza", sza_idx)
        (ialb === nothing || isza === nothing) && return false
        return scene_complete(ds, ialb, isza)
    catch err
        @warn "Candidate LUT did not validate as complete" path alb_idx sza_idx exception = (err, catch_backtrace())
        return false
    finally
        close(ds)
    end
end

function complete_anywhere(files::Vector{String}, alb_idx::Integer, sza_idx::Integer)
    for path in files
        complete_in_file(path, alb_idx, sza_idx) && return true
    end
    return false
end

function main()
    psurf = parse(Float64, get(ENV, "PSURF", get(ENV, "PSURFS", "1000")))
    sza_idx = parse(Int, get(ENV, "SZA_IDX", get(ENV, "SZA_IDXS", "1")))
    requested = parse_selector(get(ENV, "ALBEDO_IDXS", "1,7,21"))
    files = candidate_files(psurf)

    missing = Int[]
    for alb_idx in requested
        complete_anywhere(files, alb_idx, sza_idx) || push!(missing, alb_idx)
    end

    completed = setdiff(requested, missing)
    println(join(missing, ","))
    @info "Priority preflight" psurf sza_idx requested completed missing candidates = length(files)
end

main()
