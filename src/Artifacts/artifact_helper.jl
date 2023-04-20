#=
 
This file contains helper functions related to downloading/using artifacts
 
=#

""" Shorthand for @artifact_str """
function artifact_helper(name::AbstractString) 
    @show @__DIR__
    joinpath(ensure_artifact_installed(name, find_artifacts_toml(@__DIR__), 
                                       quiet_download = true), name) * ".par"
    
end

""" 
    $(FUNCTIONNAME)(molecule::AbstractString, database::AbstractString = "hitran")

Given a molecule name and a database name, retrieve the transition file 
"""
function artifact(molecule::AbstractString, 
                  database::AbstractString = "hitran")
    @assert (database == "hitran") "This package currently only supports the hitran database!"
    @assert molecule in Absorption.mol_names "Not a supported molecule!"
    mol_id = Absorption.mol_number(molecule)
    return artifact_helper("hitran_molec_id_$(mol_id)_$(molecule)")
end