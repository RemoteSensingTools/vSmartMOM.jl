""" Shorthand for @artifact_str """
artifact(name::AbstractString) = joinpath(ensure_artifact_installed(name, find_artifacts_toml(@__DIR__)), name) * ".par"
