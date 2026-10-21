#=

This file contains helper functions related to downloading/using HITRAN artifacts.

The artifact() function dispatches between the legacy Julia Artifacts system (HITRAN 2016)
and the new scratch-space download cache based on the active HITRAN edition setting.

=#

""" Shorthand for @artifact_str — resolves a legacy artifact name to its .par file path."""
function artifact_helper(name::AbstractString)
    joinpath(ensure_artifact_installed(name, find_artifacts_toml(@__DIR__),
                                       quiet_download = true), name) * ".par"
end

"""
    $(FUNCTIONNAME)(molecule::AbstractString, database::AbstractString = "hitran")

Given a molecule name, retrieve the path to its HITRAN line-by-line .par file.

The data source depends on the active HITRAN edition (see `set_hitran_edition!`):
- `"artifact"` (default): Uses legacy HITRAN 2016 data from Julia Artifacts
- Any other edition (e.g., `"HITRAN2024"`): Uses data downloaded from hitran.org.
  The cached file is used only if it verifiably holds the full dataset
  (0–150000 cm⁻¹, all isotopologues of the molecule) and its SHA-256 matches
  its metadata. Otherwise the full dataset is downloaded; if that fails, an
  error explains why the cache was rejected. A narrower file written by
  `fetch_hitran(...; numin, numax)` under the same edition is not accepted.

# Example
```julia
# Default: HITRAN 2016 from Artifacts
path = artifact("CO2")

# After switching edition:
set_hitran_edition!("HITRAN2024")
path = artifact("CO2")  # downloads from hitran.org on first call
```
"""
function artifact(molecule::AbstractString,
                  database::AbstractString = "hitran")
    @assert (database == "hitran") "This package currently only supports the hitran database!"
    @assert molecule in Absorption.mol_names "Not a supported molecule!"

    edition = get_hitran_edition()

    if edition == "artifact"
        # Legacy path — HITRAN 2016 via Julia Artifacts (unchanged behavior)
        mol_id = Absorption.mol_number(molecule)
        return artifact_helper("hitran_molec_id_$(mol_id)_$(molecule)")
    else
        # Scratch-space cache: reuse only a verified full dataset (same rule
        # as fetch_hitran), otherwise download the full dataset.
        verified, reason = _hitran_cache_status(molecule, edition,
                                                _molecule_iso_ids(molecule),
                                                _HITRAN_FULL_NUMIN, _HITRAN_FULL_NUMAX)
        verified && return _hitran_par_path(molecule, edition)
        @info "Downloading $molecule from hitran.org (edition: $edition; cache: $reason)..."
        try
            return fetch_hitran(molecule; edition=edition, force=true)
        catch err
            err isa InterruptException && rethrow()
            error("artifact(\"$molecule\") needs the full HITRAN dataset for edition " *
                  "\"$edition\" ($(_HITRAN_FULL_NUMIN)–$(_HITRAN_FULL_NUMAX) cm⁻¹, all " *
                  "isotopologues). The cache was not usable ($reason) and the download " *
                  "failed: $(sprint(showerror, err)). Retry when hitran.org is reachable, " *
                  "or call set_hitran_edition!(\"artifact\") to use the HITRAN2016 Julia artifact.")
        end
    end
end
