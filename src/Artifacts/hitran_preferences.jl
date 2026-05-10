#=

HITRAN edition preference system.

Controls whether artifact() routes through the legacy Julia Artifacts (HITRAN 2016)
or the new scratch-space download cache (e.g., HITRAN2024 from hitran.org).

=#

"""Global HITRAN edition setting. Default "artifact" uses legacy Artifacts.toml (HITRAN 2016)."""
const _HITRAN_EDITION = Ref{String}("artifact")

"""
    $(FUNCTIONNAME)(edition::AbstractString)

Set the active HITRAN edition. Valid values:
- `"artifact"` — use legacy HITRAN 2016 data from Julia Artifacts (default)
- Any other string (e.g., `"HITRAN2024"`) — use data downloaded from hitran.org,
  stored in a local scratch-space cache. Data is auto-downloaded on first access.

# Example
```julia
set_hitran_edition!("HITRAN2024")  # switch to latest HITRAN
set_hitran_edition!("artifact")    # back to legacy
```
"""
function set_hitran_edition!(edition::AbstractString)
    _HITRAN_EDITION[] = edition
    return edition
end

"""
    $(FUNCTIONNAME)()

Return the currently active HITRAN edition string.
"""
get_hitran_edition() = _HITRAN_EDITION[]

"""
    $(FUNCTIONNAME)()

List all locally available HITRAN editions. Always includes `"artifact"` (legacy HITRAN 2016).
Additional editions appear as they are downloaded to the scratch cache.
"""
function available_hitran_editions()
    editions = ["artifact"]
    scratch = _hitran_scratch_dir()
    if isdir(scratch)
        for entry in readdir(scratch)
            if isdir(joinpath(scratch, entry))
                push!(editions, entry)
            end
        end
    end
    return editions
end

"""
    $(FUNCTIONNAME)(molecule::AbstractString)

Return metadata for the given molecule under the current HITRAN edition.

For `"artifact"` edition, returns basic info. For downloaded editions,
returns the full metadata from the `.meta.toml` file (edition, download date,
wavenumber range, SHA256, source URL).
"""
function hitran_info(molecule::AbstractString)
    edition = get_hitran_edition()
    if edition == "artifact"
        return Dict{String,Any}(
            "edition" => "artifact",
            "description" => "HITRAN 2016 via Julia Artifacts",
            "molecule" => molecule,
            "source" => "http://web.gps.caltech.edu/~cfranken/hitran_2016/"
        )
    else
        meta_path = _hitran_meta_path(molecule, edition)
        if !isfile(meta_path)
            return Dict{String,Any}(
                "edition" => edition,
                "molecule" => molecule,
                "status" => "not downloaded"
            )
        end
        return _read_meta_toml(meta_path)
    end
end
