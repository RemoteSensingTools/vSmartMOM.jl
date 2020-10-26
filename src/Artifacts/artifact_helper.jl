""" Shorthand for @artifact_str """
function artifact(path_name::AbstractString) 
    return joinpath(@artifact_str(path_name), path_name) * ".par"
end