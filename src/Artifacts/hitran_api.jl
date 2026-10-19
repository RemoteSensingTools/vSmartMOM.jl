#=

HITRAN download client — Julia-native equivalent of HAPI's fetch()/fetch_by_ids().

Downloads line-by-line spectroscopic data directly from hitran.org and caches
it locally with provenance metadata for reproducibility.

A cached file is reused only when `_hitran_cache_status` verifies it: the
metadata must cover the requested band and isotopologues, and the recorded
SHA-256 must match the file. `fetch_hitran_by_ids` and `artifact` share this
rule. Downloads go to a temporary file, are validated, and only then replace
the cache.

The edition string (e.g. "HITRAN2024") is a cache/provenance label chosen by
the caller. hitran.org serves only its current database, so the label does not
prove which historical HITRAN edition a file came from.

=#

const HITRAN_API_URL = "https://hitran.org/lbl/api"

# Spectral range of a full molecule dataset; `artifact` requires this coverage.
const _HITRAN_FULL_NUMIN = 0
const _HITRAN_FULL_NUMAX = 150000

"""Return the persistent scratch directory for HITRAN data cache."""
_hitran_scratch_dir() = @get_scratch!("hitran_data")

"""Return the expected .par file path for a molecule and edition."""
function _hitran_par_path(molecule::AbstractString, edition::AbstractString)
    joinpath(_hitran_scratch_dir(), edition, "$(molecule).par")
end

"""Return the expected .meta.toml file path for a molecule and edition."""
function _hitran_meta_path(molecule::AbstractString, edition::AbstractString)
    joinpath(_hitran_scratch_dir(), edition, "$(molecule).meta.toml")
end

"""Return all HITRAN global isotopologue IDs of a molecule."""
function _molecule_iso_ids(molecule::AbstractString)
    mol_id = Absorption.mol_number(molecule)
    return Int[gid for gid in Absorption.global_ids[mol_id, :] if gid > 0]
end

"""
    $(FUNCTIONNAME)(molecule::AbstractString, edition=get_hitran_edition())

Return `true` when `artifact(molecule)` would use a downloaded file for
`edition` without downloading: the cached file must cover the full spectral
range and all isotopologues of the molecule, and its SHA-256 must match its
metadata. A narrower or unverified file returns `false`.
"""
function hitran_is_cached(molecule::AbstractString,
                          edition::AbstractString=get_hitran_edition())
    molecule in Absorption.mol_names || return false
    verified, _ = _hitran_cache_status(molecule, edition, _molecule_iso_ids(molecule),
                                       _HITRAN_FULL_NUMIN, _HITRAN_FULL_NUMAX)
    return verified
end

"""
    $(FUNCTIONNAME)(molecule::AbstractString; numin=0, numax=150000,
                    edition="HITRAN2024", force=false)

Download HITRAN line-by-line data for a molecule directly from hitran.org.

Analogous to HAPI's `fetch()`. Resolves the molecule name to all its isotopologue
global IDs and queries the HITRAN API. Data is cached locally; subsequent calls
return the cached path unless `force=true` or the cache cannot be verified for
this request (see `fetch_hitran_by_ids`).

# Arguments
- `molecule`: Molecule name (e.g., "CO2", "H2O", "O3")
- `numin`: Lower wavenumber bound in cm⁻¹ (default: 0)
- `numax`: Upper wavenumber bound in cm⁻¹ (default: 150000)
- `edition`: Cache/provenance label for this download (default: "HITRAN2024").
  hitran.org always serves its current database; the label names the local
  cache directory and does not select or prove a historical HITRAN edition.
- `force`: Re-download even if cached (default: false)

# Returns
- Path to the downloaded `.par` file

# Example
```julia
path = fetch_hitran("CO2")
path = fetch_hitran("O3"; numin=500, numax=3000, edition="HITRAN2024")
```
"""
function fetch_hitran(molecule::AbstractString;
                      numin::Real=_HITRAN_FULL_NUMIN, numax::Real=_HITRAN_FULL_NUMAX,
                      edition::AbstractString="HITRAN2024",
                      force::Bool=false)
    @assert molecule in Absorption.mol_names "Unknown molecule: $molecule. Use Absorption.show_molecules() to list available molecules."

    return fetch_hitran_by_ids(molecule, _molecule_iso_ids(molecule);
                               numin=numin, numax=numax,
                               edition=edition, force=force)
end

"""
    $(FUNCTIONNAME)(name::AbstractString, iso_ids::Vector{Int};
                    numin=0, numax=150000, edition="HITRAN2024", force=false)

Download HITRAN data using explicit global isotopologue IDs.

Lower-level variant of `fetch_hitran` for when you need specific isotopologues.
Global IDs can be looked up via `Absorption.mol_globalID(mol, iso)`.

A cached file under the same `name` and `edition` is reused only if its
metadata covers `[numin, numax]`, includes every requested isotopologue, and
records a SHA-256 that matches the file. A reused file may therefore contain a
wider band or additional isotopologues. Missing or inconsistent metadata
triggers a fresh download.

A download is written to a temporary file and checked (non-empty, well-formed
fixed-width HITRAN records, only requested isotopologues within the requested
band) before it replaces the cache. If the download or check fails, an
existing cache is left unchanged. After a replacement, the process-local
parsed spectroscopy cache (`clear_spectroscopy_cache!`) is cleared.

# Arguments
- `name`: Label for the file (typically molecule name)
- `iso_ids`: Vector of HITRAN global isotopologue IDs
- `numin`, `numax`, `edition`, `force`: Same as `fetch_hitran`

# Returns
- Path to the downloaded `.par` file
"""
function fetch_hitran_by_ids(name::AbstractString, iso_ids::Vector{Int};
                              numin::Real=_HITRAN_FULL_NUMIN, numax::Real=_HITRAN_FULL_NUMAX,
                              edition::AbstractString="HITRAN2024",
                              force::Bool=false)
    par_path  = _hitran_par_path(name, edition)
    meta_path = _hitran_meta_path(name, edition)

    # Reuse a cached file only if it verifiably covers this request. The cache
    # is keyed on (name, edition); without this check a narrower-band or
    # fewer-isotopologue fetch made earlier would be silently reused, producing
    # an LUT with lines missing (root cause of the 2026-04-23 H2O truncation).
    if !force
        verified, reason = _hitran_cache_status(name, edition, iso_ids, numin, numax)
        if verified
            @info "Using cached HITRAN data: $par_path"
            return par_path
        elseif isfile(par_path)
            @info "Cached HITRAN data for $name not reused ($reason); re-downloading" par_path
        end
    end

    # Construct API URL (same pattern as HAPI's queryHITRAN)
    ids_str = join(iso_ids, ",")
    url = "$(HITRAN_API_URL)?iso_ids_list=$(ids_str)&numin=$(numin)&numax=$(numax)"

    edition_dir = dirname(par_path)
    mkpath(edition_dir)
    _remove_stale_hitran_files(par_path, meta_path)

    @info "Downloading HITRAN data for $name from hitran.org..." url=url

    # Download and describe the new file in temporary files next to the cache
    # (same file system, so the final renames are atomic).
    tmp_par  = _hitran_sidecar_path(par_path, "tmp")
    tmp_meta = _hitran_sidecar_path(meta_path, "tmp")
    try
        _HITRAN_DOWNLOAD[](url, tmp_par)
        _validate_hitran_download(tmp_par, iso_ids, numin, numax)
        _write_meta_toml(tmp_meta;
            molecule=name,
            edition=edition,
            download_date=string(Dates.now()),
            numin=numin,
            numax=numax,
            iso_ids=ids_str,
            source_url=url,
            sha256=_sha256_hex(tmp_par),
            file_size_bytes=filesize(tmp_par)
        )
        # Parsed line databases are keyed by path, so drop them together with
        # the old file. Lock order matches `CoreRT._hitran_lines` (parsed
        # cache first, then file cache), which calls `artifact` under its lock.
        lock(CoreRT._HITRAN_LINES_LOCK) do
            lock(_HITRAN_CACHE_LOCK) do
                _replace_hitran_cache!(tmp_par, par_path, tmp_meta, meta_path)
            end
            CoreRT.clear_spectroscopy_cache!()
        end
    finally
        rm(tmp_par; force=true)
        rm(tmp_meta; force=true)
    end

    @info "Downloaded $(name).par ($(Base.format_bytes(filesize(par_path)))) to $edition_dir"
    return par_path
end

# --- Cache verification ---

"""
    _hitran_cache_status(name, edition, iso_ids, numin, numax) -> (verified, reason)

The single reuse rule shared by `fetch_hitran_by_ids` and `artifact`. The
cached `name.par` for `edition` is verified only if its metadata

- names the same molecule label and edition,
- covers `[numin, numax]`,
- lists every requested global isotopologue ID, and
- records a size and SHA-256 that match the file on disk.

Missing, unreadable, or inconsistent metadata means "not verified"; `reason`
says why. Errors while checking are reported the same way, never thrown. Old
caches written by earlier versions carry all of these fields and verify
unchanged. The check holds `_HITRAN_CACHE_LOCK`, so within one process it
never sees a half-published replacement.
"""
function _hitran_cache_status(name::AbstractString, edition::AbstractString,
                              iso_ids::AbstractVector{<:Integer},
                              numin::Real, numax::Real)
    return lock(_HITRAN_CACHE_LOCK) do
        try
            _check_hitran_cache(name, edition, iso_ids, numin, numax)
        catch err
            err isa InterruptException && rethrow()
            (false, "cache could not be checked: $(sprint(showerror, err))")
        end
    end
end

function _check_hitran_cache(name, edition, iso_ids, numin, numax)
    par_path  = _hitran_par_path(name, edition)
    meta_path = _hitran_meta_path(name, edition)
    isfile(par_path)  || return (false, "no cached file")
    isfile(meta_path) || return (false, "provenance metadata $meta_path is missing")
    meta = _read_meta_toml(meta_path)

    get(meta, "molecule", nothing) == name ||
        return (false, "metadata molecule $(repr(get(meta, "molecule", nothing))) is not $(repr(name))")
    get(meta, "edition", nothing) == edition ||
        return (false, "metadata edition $(repr(get(meta, "edition", nothing))) is not $(repr(edition))")

    cnumin, cnumax = get(meta, "numin", nothing), get(meta, "numax", nothing)
    (cnumin isa Real && cnumax isa Real) ||
        return (false, "metadata has no numeric numin/numax")
    cnumin <= numin && cnumax >= numax ||
        return (false, "cached band [$cnumin, $cnumax] cm⁻¹ does not cover [$numin, $numax] cm⁻¹")

    cached_ids = _parse_iso_ids(get(meta, "iso_ids", nothing))
    cached_ids === nothing && return (false, "metadata has no readable iso_ids list")
    missing_ids = setdiff(iso_ids, cached_ids)
    isempty(missing_ids) ||
        return (false, "cached file lacks requested isotopologue global IDs $(join(missing_ids, ","))")

    sha, nbytes = get(meta, "sha256", nothing), get(meta, "file_size_bytes", nothing)
    (sha isa AbstractString && nbytes isa Real) ||
        return (false, "metadata has no sha256/file_size_bytes")
    filesize(par_path) == nbytes ||
        return (false, "file size $(filesize(par_path)) differs from metadata ($nbytes)")
    _cached_sha256_hex(par_path) == sha ||
        return (false, "file SHA-256 differs from metadata")
    return (true, "verified")
end

"""Parse a metadata `iso_ids` entry ("7,8,9" or a single number); `nothing` if unreadable."""
function _parse_iso_ids(value)
    value isa Real && isinteger(value) && abs(value) < 2^31 && return [Int(value)]
    value isa AbstractString || return nothing
    ids = tryparse.(Int, strip.(split(value, ',')))
    return any(isnothing, ids) ? nothing : Int.(ids)
end

_sha256_hex(path::AbstractString) = open(io -> bytes2hex(sha256(io)), path)

# Serializes cache checks and cache replacement within this process.
const _HITRAN_CACHE_LOCK = ReentrantLock()

# `artifact` runs the cache check on every model construction, so file hashes
# are remembered per path (guarded by `_HITRAN_CACHE_LOCK`). Any rewrite or
# rename changes the file's inode, size, mtime or ctime, and the hash is then
# recomputed. File timestamps can be coarse, so a hash is remembered only if
# the file was already a few seconds old before it was read and is unchanged
# afterwards; a later rewrite then necessarily carries a different timestamp.
const _HITRAN_HASHES = Dict{String,Tuple{Tuple{UInt64,Int64,Float64,Float64},String}}()

_file_identity(st) = (UInt64(st.inode), Int64(st.size), st.mtime, st.ctime)

function _cached_sha256_hex(path::AbstractString)
    started = time()
    st = stat(path)
    key = _file_identity(st)
    cached = get(_HITRAN_HASHES, path, nothing)
    cached !== nothing && first(cached) == key && return last(cached)
    delete!(_HITRAN_HASHES, path)
    hash = _sha256_hex(path)
    if max(st.mtime, st.ctime) < started - 2 && _file_identity(stat(path)) == key
        _HITRAN_HASHES[path] = (key, hash)
    end
    return hash
end

# --- Download, validation and replacement ---

"""
    _download_hitran(url, path)

Download `url` into `path`, throwing unless the transfer completed with a 2xx
HTTP status. Transport errors are failures, because a body truncated at a
record boundary passes record validation. The one exception is hitran.org's
known length mismatch (see `_is_hitran_newline_shortfall`); the body is still
checked record by record by `_validate_hitran_download` before use.
"""
function _download_hitran(url::AbstractString, path::AbstractString)
    resp = open(path, "w") do output
        Downloads.request(url; output=output, throw=false)
    end
    _check_hitran_response(resp, path)
    return path
end

function _check_hitran_response(resp, path::AbstractString)
    status = resp isa Downloads.Response ? resp.status :
             resp isa Downloads.RequestError ? resp.response.status : 0
    status == 403 && error("HITRAN API rate limit exceeded. You have exceeded the daily limit of API queries. Try again tomorrow.")
    if resp isa Downloads.RequestError
        200 <= status < 300 && _is_hitran_newline_shortfall(resp, path) && return nothing
        error("HITRAN download failed: $(resp.message)" *
              (status > 0 ? " (after HTTP status $status)" : ""))
    end
    200 <= status < 300 || error("HITRAN API request failed with HTTP status $status")
    return nothing
end

"""
    _is_hitran_newline_shortfall(err::Downloads.RequestError, path)

hitran.org declares a Content-Length that counts two bytes (CRLF) per line but
sends LF-terminated lines, so curl ends a complete transfer with error 18
("N bytes missing"), N being the number of records. Accept exactly that case:
error code 18; every Content-Length header agreeing on a length equal to the
received size plus the number of newlines; LF-only line endings; and a body
ending in a newline. Each record then accounts for 162 declared and 161
received bytes, so the equality means every declared record arrived, provided
the server over-counts as observed. A response whose declared length were
correct for LF endings and that lost exactly one record in 162 at a record
boundary would also match; this cannot be told apart from the server defect.
Any other transport error remains a failure.
"""
function _is_hitran_newline_shortfall(err, path::AbstractString)
    err.code == 18 || return false
    lengths = [tryparse(Int, strip(v)) for (k, v) in err.response.headers
               if lowercase(k) == "content-length"]
    (isempty(lengths) || any(isnothing, lengths) || !allequal(lengths)) && return false
    declared = first(lengths)
    nbytes = filesize(path)
    nbytes == 0 && return false
    newlines, carriage_returns, last_byte = open(path) do io
        n, r = 0, 0
        while !eof(io)
            chunk = read(io, 1 << 20)
            n += count(==(UInt8('\n')), chunk)
            r += count(==(UInt8('\r')), chunk)
        end
        seekend(io); skip(io, -1)
        n, r, read(io, UInt8)
    end
    return carriage_returns == 0 && last_byte == UInt8('\n') && declared == nbytes + newlines
end

# Download callable used by `fetch_hitran_by_ids`; `(url, path) -> path`.
# Tests replace it with offline fixtures.
const _HITRAN_DOWNLOAD = Ref{Any}(_download_hitran)

"""
    _validate_hitran_download(path, iso_ids, numin, numax)

Check a downloaded file before it may replace the cache. It must be non-empty;
every line must be a fixed-width 160-column HITRAN record whose numeric fields
parse as in `Absorption.read_hitran`; and every line must belong to a
requested isotopologue and lie within `[numin, numax]`. Returns the number of
records.
"""
function _validate_hitran_download(path::AbstractString, iso_ids::AbstractVector{<:Integer},
                                   numin::Real, numax::Real)
    filesize(path) > 0 || error("HITRAN download produced an empty file. Check your query parameters (numin=$numin, numax=$numax).")
    A = Absorption
    requested = Set(iso_ids)
    # Read through an explicit handle: an `eachline(path)` left by an error
    # keeps the file open, and Windows cannot then delete the download.
    nrecords = open(path) do io
        nread = 0
        for (n, record) in enumerate(eachline(io))
            # `read_hitran` tolerates trailing bytes; a download must not, since
            # two records joined without a newline would otherwise pass as one.
            ncodeunits(record) == A._HITRAN_RECORD_LENGTH ||
                error("$path:$n: HITRAN record has $(ncodeunits(record)) bytes; " *
                      "expected exactly $(A._HITRAN_RECORD_LENGTH)")
            fields = A._hitran_fields(record, path, n)
            mol = A._parse_hitran_number(Int, fields[1], path, n, "molec_id")
            iso = A._parse_hitran_isotopologue(fields[2], path, n)
            ν = A._parse_hitran_number(Float64, fields[3], path, n, "nu")
            for (i, field) in zip(4:10, ("sw", "a", "gamma_air", "gamma_self",
                                         "elower", "n_air", "delta_air"))
                A._parse_hitran_number(Float64, fields[i], path, n, field)
            end
            A._parse_hitran_number(Float64, fields[18], path, n, "gp"; blank_is_zero=true)
            A._parse_hitran_number(Float64, fields[19], path, n, "gpp"; blank_is_zero=true)

            gid = checkbounds(Bool, A.global_ids, mol, iso) ? Int(A.global_ids[mol, iso]) : -1
            gid in requested || error("$path:$n: line of molecule $mol, isotopologue $iso " *
                                      "was not requested (global IDs $(join(iso_ids, ",")))")
            numin <= ν <= numax || error("$path:$n: line at $ν cm⁻¹ is outside the " *
                                         "requested range [$numin, $numax] cm⁻¹")
            nread += 1
        end
        nread
    end
    nrecords > 0 || error("HITRAN download contains no line records.")
    return nrecords
end

"""
    _replace_hitran_cache!(tmp_par, par_path, tmp_meta, meta_path; rename, remove)

Move a validated data file and its metadata into place; callers hold
`_HITRAN_CACHE_LOCK`. The previous data file is first hard-linked (or copied)
to a `.<name>.par.bak.*` backup, so `par_path` exists throughout. The data
file is then replaced in one rename, followed by the metadata. If either
rename fails, the backup is renamed back, so a failed replacement leaves the
previous pair intact. A crash in between can leave new data beside old
metadata; the cache check then reports a SHA-256 mismatch and the data are
downloaded again rather than reused. The same check catches a mismatched pair
left by two processes replacing the same file at once.
"""
function _replace_hitran_cache!(tmp_par, par_path, tmp_meta, meta_path;
                                rename=_atomic_rename, remove=rm)
    backup = _hitran_sidecar_path(par_path, "bak")
    had_previous = _backup_hitran_file(par_path, backup)
    try
        rename(tmp_par, par_path)
        rename(tmp_meta, meta_path)
    catch
        if !had_previous
            rm(par_path; force=true)
        elseif ispath(par_path) && samefile(backup, par_path)
            rm(backup; force=true)  # data file never replaced; rename(2) would be a no-op
        else
            _atomic_rename(backup, par_path)
        end
        rethrow()
    end
    # The new pair is in place; failing to delete the backup must not undo that.
    try
        had_previous && remove(backup; force=true)
    catch err
        @warn "Could not remove previous HITRAN file $backup" exception=err
    end
    return par_path
end

"""Link (or, where hard links are unsupported, copy) `path` to `backup`; `false` if `path` does not exist."""
function _backup_hitran_file(path, backup)
    try
        hardlink(path, backup)
    catch err
        err isa Base.IOError && err.code == Base.UV_ENOENT && return false
        isfile(path) || return false
        cp(path, backup)
    end
    return true
end

"""
Return a unique hidden temporary/backup path `.<file>.<kind>.<pid>-<random>`
next to `path`, recognisable by `_remove_stale_hitran_files`. The random part
uses `Libc.rand`, leaving the user's random-number streams untouched.
"""
function _hitran_sidecar_path(path, kind)
    chars = "0123456789abcdefghijklmnopqrstuvwxyz"
    suffix = String([chars[Libc.rand() % 36 + 1] for _ in 1:12])
    return joinpath(dirname(path), ".$(basename(path)).$(kind).$(getpid())-$(suffix)")
end

const _HITRAN_SIDECAR_SUFFIX = r"^\.(tmp|bak)\.[0-9]+-[0-9a-z]{12}$"

"""
Delete temporary and backup files of this cache entry left by an interrupted
download or replacement, once they are older than `max_age` seconds. Only
names produced by `_hitran_sidecar_path` for these two files are touched.
"""
function _remove_stale_hitran_files(par_path, meta_path; max_age=86400)
    dir = dirname(par_path)
    stems = (".$(basename(par_path))", ".$(basename(meta_path))")
    for file in readdir(dir)
        any(stem -> startswith(file, stem) &&
                    occursin(_HITRAN_SIDECAR_SUFFIX, file[ncodeunits(stem)+1:end]), stems) || continue
        path = joinpath(dir, file)
        try
            st = stat(path)
            max(st.mtime, st.ctime) < time() - max_age && rm(path; force=true)
        catch err
            @warn "Could not remove stale HITRAN file $path" exception=err
        end
    end
end

"""
Rename `src` to `dst`, replacing `dst` in one step or failing. Unlike
`Base.Filesystem.rename` on Julia 1.10/1.11, this never falls back to copying.
"""
function _atomic_rename(src::AbstractString, dst::AbstractString)
    err = ccall(:jl_fs_rename, Int32, (Cstring, Cstring), src, dst)
    err < 0 && Base.uv_error("rename($(repr(src)), $(repr(dst)))", err)
    return dst
end

# --- Metadata TOML helpers ---

function _write_meta_toml(path::AbstractString; kwargs...)
    open(path, "w") do io
        println(io, "# HITRAN download metadata")
        println(io, "# Auto-generated by vSmartMOM.jl — do not edit")
        println(io)
        for (k, v) in pairs(kwargs)
            if v isa AbstractString
                println(io, "$k = \"$v\"")
            else
                println(io, "$k = $v")
            end
        end
    end
end

function _read_meta_toml(path::AbstractString)
    meta = Dict{String,Any}()
    # Explicit handle, as in `_validate_hitran_download`: closed on errors too.
    open(path) do io
        for line in eachline(io)
            line = strip(line)
            (isempty(line) || startswith(line, '#')) && continue
            if occursin(" = ", line)
                key, val = split(line, " = "; limit=2)
                # Strip quotes from string values
                if length(val) >= 2 && startswith(val, '"') && endswith(val, '"')
                    meta[strip(key)] = chop(val; head=1, tail=1)
                else
                    # Try numeric parse
                    meta[strip(key)] = tryparse(Float64, val) !== nothing ? tryparse(Float64, val) : val
                end
            end
        end
    end
    return meta
end
