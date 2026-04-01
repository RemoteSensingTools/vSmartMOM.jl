#=

HITRAN download client — Julia-native equivalent of HAPI's fetch()/fetch_by_ids().

Downloads line-by-line spectroscopic data directly from hitran.org and caches
it locally with version metadata for reproducibility.

=#

const HITRAN_API_URL = "https://hitran.org/lbl/api"

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

"""Check whether HITRAN data for a molecule is already cached for the given edition."""
function hitran_is_cached(molecule::AbstractString,
                          edition::AbstractString=get_hitran_edition())
    isfile(_hitran_par_path(molecule, edition))
end

"""
    $(FUNCTIONNAME)(molecule::AbstractString; numin=0, numax=150000,
                    edition="HITRAN2024", force=false)

Download HITRAN line-by-line data for a molecule directly from hitran.org.

Analogous to HAPI's `fetch()`. Resolves the molecule name to all its isotopologue
global IDs and queries the HITRAN API. Data is cached locally; subsequent calls
return the cached path unless `force=true`.

# Arguments
- `molecule`: Molecule name (e.g., "CO2", "H2O", "O3")
- `numin`: Lower wavenumber bound in cm⁻¹ (default: 0)
- `numax`: Upper wavenumber bound in cm⁻¹ (default: 150000)
- `edition`: Label for this download (default: "HITRAN2024")
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
                      numin::Real=0, numax::Real=150000,
                      edition::AbstractString="HITRAN2024",
                      force::Bool=false)
    @assert molecule in Absorption.mol_names "Unknown molecule: $molecule. Use Absorption.show_molecules() to list available molecules."

    mol_id = Absorption.mol_number(molecule)

    # Gather all global isotopologue IDs for this molecule
    iso_ids = Int[]
    for iso in 1:size(Absorption.global_ids, 2)
        gid = Absorption.global_ids[mol_id, iso]
        if gid > 0
            push!(iso_ids, gid)
        end
    end

    return fetch_hitran_by_ids(molecule, iso_ids;
                               numin=numin, numax=numax,
                               edition=edition, force=force)
end

"""
    $(FUNCTIONNAME)(name::AbstractString, iso_ids::Vector{Int};
                    numin=0, numax=150000, edition="HITRAN2024", force=false)

Download HITRAN data using explicit global isotopologue IDs.

Lower-level variant of `fetch_hitran` for when you need specific isotopologues.
Global IDs can be looked up via `Absorption.mol_globalID(mol, iso)`.

# Arguments
- `name`: Label for the file (typically molecule name)
- `iso_ids`: Vector of HITRAN global isotopologue IDs
- `numin`, `numax`, `edition`, `force`: Same as `fetch_hitran`

# Returns
- Path to the downloaded `.par` file
"""
function fetch_hitran_by_ids(name::AbstractString, iso_ids::Vector{Int};
                              numin::Real=0, numax::Real=150000,
                              edition::AbstractString="HITRAN2024",
                              force::Bool=false)
    par_path = _hitran_par_path(name, edition)

    # Return cached file unless force re-download
    if isfile(par_path) && !force
        @info "Using cached HITRAN data: $par_path"
        return par_path
    end

    # Construct API URL (same pattern as HAPI's queryHITRAN)
    ids_str = join(iso_ids, ",")
    url = "$(HITRAN_API_URL)?iso_ids_list=$(ids_str)&numin=$(numin)&numax=$(numax)"

    # Ensure directory exists
    edition_dir = joinpath(_hitran_scratch_dir(), edition)
    mkpath(edition_dir)

    @info "Downloading HITRAN data for $name from hitran.org..." url=url

    # Download with error handling.
    # HITRAN's server uses chunked transfer encoding without Content-Length,
    # which causes Downloads.jl to report "bytes missing" RequestError even on
    # successful downloads. We use throw=false and validate the file content.
    local resp
    open(par_path, "w") do output
        resp = Downloads.request(url; output=output, throw=false)
    end

    # Check for real HTTP errors (rate limiting, server errors)
    http_status = if resp isa Downloads.Response
        resp.status
    elseif resp isa Downloads.RequestError
        resp.response isa Downloads.Response ? resp.response.status : 0
    else
        0
    end
    if http_status == 403
        rm(par_path; force=true)
        error("HITRAN API rate limit exceeded. You have exceeded the daily limit of API queries. Try again tomorrow.")
    elseif http_status >= 400
        rm(par_path; force=true)
        error("HITRAN API request failed with HTTP status $http_status")
    end

    # Verify we actually got valid HITRAN data
    if !isfile(par_path) || filesize(par_path) == 0
        rm(par_path; force=true)
        error("HITRAN download produced an empty file. Check your query parameters (numin=$numin, numax=$numax).")
    end

    # Compute SHA256 of downloaded file
    file_hash = open(par_path, "r") do io
        bytes2hex(sha256(io))
    end

    # Write metadata
    meta_path = _hitran_meta_path(name, edition)
    _write_meta_toml(meta_path;
        molecule=name,
        edition=edition,
        download_date=string(Dates.now()),
        numin=numin,
        numax=numax,
        iso_ids=ids_str,
        source_url=url,
        sha256=file_hash,
        file_size_bytes=filesize(par_path)
    )

    @info "Downloaded $(name).par ($(Base.format_bytes(filesize(par_path)))) to $edition_dir"
    return par_path
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
    for line in eachline(path)
        line = strip(line)
        (isempty(line) || startswith(line, '#')) && continue
        if occursin(" = ", line)
            key, val = split(line, " = "; limit=2)
            # Strip quotes from string values
            if startswith(val, '"') && endswith(val, '"')
                meta[strip(key)] = val[2:end-1]
            else
                # Try numeric parse
                meta[strip(key)] = tryparse(Float64, val) !== nothing ? tryparse(Float64, val) : val
            end
        end
    end
    return meta
end
