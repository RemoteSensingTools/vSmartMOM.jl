""" Shorthand for @artifact_str """
function artifact_helper(name::AbstractString) 

    # The test here is whether to include the keyword quiet_download
    # The progress bar on artifacts >1.5 is wonky (especially when generating docs)
    if (VERSION > v"1.5")
        joinpath(ensure_artifact_installed(name, find_artifacts_toml(@__DIR__), 
                                       quiet_download = true), name) * ".par"
    else
        joinpath(ensure_artifact_installed(name, find_artifacts_toml(@__DIR__)), name) * ".par"
    end
    
end

""" Given a molecule name and a database name, retrieve the transition file """
function artifact(molecule::AbstractString, 
                  database::AbstractString = "hitran")
    @assert (database == "hitran") "This package currently only supports the hitran database!"
    @assert molecule in CrossSection.mol_names "Not a supported molecule!"
    mol_id = CrossSection.mol_number(molecule)
    return artifact_helper("hitran_molec_id_$(mol_id)_$(molecule)")
end