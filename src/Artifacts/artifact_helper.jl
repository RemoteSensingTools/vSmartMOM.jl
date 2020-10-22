""" Shorthand for @artifact_str """
artifact(name) = joinpath(@artifact_str(name), name) * ".par"