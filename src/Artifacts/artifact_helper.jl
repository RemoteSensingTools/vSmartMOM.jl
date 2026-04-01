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
- Any other edition (e.g., `"HITRAN2024"`): Uses data downloaded from hitran.org,
  auto-downloading on first access if not already cached.

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
        # New path — scratch-space cache with auto-download
        par_path = _hitran_par_path(molecule, edition)
        if !isfile(par_path)
            @info "Downloading $molecule from hitran.org (edition: $edition)..."
            fetch_hitran(molecule; edition=edition)
        end
        return par_path
    end
end
