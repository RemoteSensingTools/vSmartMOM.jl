using Pkg.Artifacts
usinh HTTP

# This is the path to the Artifacts.toml we will manipulate
artifact_toml = joinpath(@__DIR__, "Artifacts.toml")

# Query the `Artifacts.toml` file for the hash bound to the name "iris"
# (returns `nothing` if no such binding exists)
hitran_hash = artifact_hash("hitran", artifact_toml)

# If the name was not bound, or the hash it was bound to does not exist, create it!
if hitran_hash == nothing || !artifact_exists(hitran_hash)
    # create_artifact() returns the content-hash of the artifact directory once we're finished creating it
    hitran_hash = create_artifact() do artifact_dir
        # We create the artifact by simply downloading a few files into the new artifact directory
        hitran_url_base = "https://hitran.org/lbl/api?numin=0&numax=50000&iso_ids_list="
        for i=1:130
            println(i)
            try
                download(hitran_url_base*string(i), joinpath(artifact_dir, "hitran_global_id_" * string(i)*".par"))
            catch
            end
        end
    end

    # Now bind that hash within our `Artifacts.toml`.  `force = true` means that if it already exists,
    # just overwrite with the new content-hash.  Unless the source files change, we do not expect
    # the content hash to change, so this should not cause unnecessary version control churn.
    bind_artifact!(artifact_toml, "hitran", hitran_hash)
end

# Get the path of the iris dataset, either newly created or previously generated.
# this should be something like `~/.julia/artifacts/dbd04e28be047a54fbe9bf67e934be5b5e0d357a`
hitran_dataset_path = artifact_path(hitran_hash)
