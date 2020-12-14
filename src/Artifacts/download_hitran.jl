using Pkg.Artifacts
using RadiativeTransfer.CrossSection
using HTTP

deploy = true

isos = CrossSection.iso_info

hitran_url_base    = "https://hitran.org/lbl/api?numin=0&numax=150000&iso_ids_list="
local_download_url = "/net/fluo/data1/ftp/XYZT_hitran/"
target_url = "http://web.gps.caltech.edu/~cfranken/hitran_2016/"
# Loop over all current HITRAN molecules:
for molec=1:49
    if haskey(isos,(molec,1))
        locs = [k for (k,v) in isos if k[1]==molec]
        
        global_ids = ""
        for i in locs
            global_ids *= ","*string(isos[i][1])
            molec_name = isos[i][5] 
        end
        println(global_ids[2:end])
        url = hitran_url_base*global_ids[2:end]
        
        filename = "hitran_molec_id_" * string(molec)*"_" * isos[molec,1][5] * ".par"
        name = splitext(filename)[1]
        @info("Generating artifact for $(filename)")
        HTTP.open(:GET, url) do http
            open(local_download_url*filename, "w") do file
                write(file, http)
            end
        end
        hash = create_artifact() do artifact_dir
            # Copy in weights
            cp(joinpath(local_download_url, filename), joinpath(artifact_dir, filename))
        end

        # Spit tarballs to be hosted out to local temporary directory:
        if deploy
            tarball_hash = archive_artifact(hash, joinpath(local_download_url, "$(name).tar.gz"))

            # Calculate tarball url
            tarball_url = target_url*name*".tar.gz"
            # Bind this to an Artifacts.toml file
            @info("Binding $(name) in Artifacts.toml...")
            bind_artifact!(joinpath(@__DIR__, "Artifacts.toml"), name, hash; download_info=[(tarball_url, tarball_hash)], lazy=true, force=true)
        end
    end
end