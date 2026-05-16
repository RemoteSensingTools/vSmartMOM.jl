# Full-grid GCHP aerosol optical-depth diagnostic.
#
# This is intentionally outside CoreRT: it validates the microphysics → Mie
# extinction boundary without changing the production RT aerosol path.

"""
    write_gchp_aod_diagnostic(out_path, gchp_path; wavelengths_um, ...)

Compute aerosol optical depth for every requested GCHP grid cell and write a
NetCDF file with dimensions `(x, y, face, wavelength)`.

The default wavelengths are `[0.760, 1.600, 2.200]` μm. The GCHP file is
opened once and aerosol ingest defaults to `:auto`.
"""
function write_gchp_aod_diagnostic(out_path::AbstractString,
                                   gchp_path::AbstractString;
                                   wavelengths_um = [0.760, 1.600, 2.200],
                                   aerosol_scheme = :auto,
                                   faces = nothing,
                                   xrange = nothing,
                                   yrange = nothing,
                                   integration = nothing,
                                   mie_lut = nothing,
                                   ri_round_digits = nothing,
                                   FT::DataType = Float64)
    λs = FT.(wavelengths_um)
    mie_cache = Dict{Any, FT}()

    GCHPFile(gchp_path; aerosol_scheme=aerosol_scheme, FT=FT) do f
        face_iter = faces === nothing ? (1:f.nf) : faces
        x_iter    = xrange === nothing ? (1:f.nx) : xrange
        y_iter    = yrange === nothing ? (1:f.ny) : yrange

        NCDataset(String(out_path), "c") do ds
            defDim(ds, "x", f.nx)
            defDim(ds, "y", f.ny)
            defDim(ds, "face", f.nf)
            defDim(ds, "wavelength", length(λs))

            λ_var = defVar(ds, "wavelength_um", Float64, ("wavelength",))
            λ_var[:] = Float64.(λs)
            λ_var.attrib["units"] = "micrometer"
            λ_var.attrib["long_name"] = "wavelength"

            lat = defVar(ds, "lat", Float64, ("x", "y", "face"))
            lon = defVar(ds, "lon", Float64, ("x", "y", "face"))
            lat[:, :, :] = f.lats
            lon[:, :, :] = f.lons
            lat.attrib["units"] = "degree_north"
            lon.attrib["units"] = "degree_east"

            aod = defVar(ds, "aod", Float64, ("x", "y", "face", "wavelength"))
            aod[:, :, :, :] .= NaN
            aod.attrib["long_name"] = "aerosol optical depth from sectional aerosol microphysics"
            aod.attrib["units"] = "1"

            for idf in face_iter, idy in y_iter, idx in x_iter
                scene = scene_at(f, idx, idy, idf; FT=FT)
                aod_kwargs = integration === nothing ?
                    (; mie_cache=mie_cache, mie_lut=mie_lut,
                     ri_round_digits=ri_round_digits) :
                    (; integration=integration, mie_cache=mie_cache,
                     mie_lut=mie_lut, ri_round_digits=ri_round_digits)
                aod[idx, idy, idf, :] = Float64.(compute_scene_aod(
                    scene, λs; aod_kwargs...))
            end

            ds.attrib["source_file"] = String(gchp_path)
            ds.attrib["aerosol_scheme"] = string(aerosol_scheme)
            ds.attrib["time"] = string(f.time)
            ds.attrib["aod_algorithm"] = "Mie extinction at sectional bin centers with volume-weighted effective RI"
            ds.attrib["ri_round_digits"] = ri_round_digits === nothing ? "none" : string(ri_round_digits)
            ds.attrib["mie_lut"] = mie_lut === nothing ? "none" : string(typeof(mie_lut))
            ds.attrib["integration"] = integration === nothing ? "scene default" : string(integration)
        end
    end

    return String(out_path)
end
