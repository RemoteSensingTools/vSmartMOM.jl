# =============================================================================
# Benchmark NetCDF writer + scene-loop driver (gchp-io PR).
#
# Emits NetCDF-4 (HDF5 on disk) via NCDatasets — already a dep, no new
# packages. Python consumers read the files with `netCDF4`/`xarray`/`h5py`.
# =============================================================================

"""
    write_scene_result(path, scene::GCHPScene, result::SweepResult;
                       mode = :overwrite)

Write a sweep result + scene metadata to a NetCDF-4 file. The file
contains:

- Dimensions: `sza`, `view_pair`, `brdf`, `pol`, `spec`
- Variables: `R(sza, view_pair, brdf, pol, spec)`,
  `T(sza, view_pair, brdf, pol, spec)`, plus the sweep axes (`sza`, `vza`,
  `vaz`) and `lat`/`lon`/`time` scene metadata.

`mode = :overwrite` truncates the file. Pass `:append_scene` to add this
scene to an existing file as a new group `scene_<NNNN>` — used by
[`generate_benchmark`](@ref) when looping over many scenes.
"""
function write_scene_result(path::AbstractString,
                            scene::GCHPScene,
                            result::SweepResult;
                            mode::Symbol = :overwrite,
                            scene_index::Integer = 0)
    if mode === :overwrite
        NCDataset(String(path), "c") do ds
            _write_scene_group(ds, scene, result)
        end
    elseif mode === :append_scene
        NCDataset(String(path), "a") do ds
            grp = defGroup(ds, @sprintf("scene_%04d", scene_index))
            _write_scene_group(grp, scene, result)
        end
    else
        throw(ArgumentError("Unknown mode :$mode (use :overwrite or :append_scene)"))
    end
    return path
end

function _write_scene_group(ds, scene::GCHPScene, result::SweepResult)
    N_sza, N_vp, N_brdf, pol_n, nSpec = size(result.R)

    defDim(ds, "sza",        N_sza)
    defDim(ds, "view_pair",  N_vp)
    defDim(ds, "brdf",       N_brdf)
    defDim(ds, "pol",        pol_n)
    defDim(ds, "spec",       nSpec)

    sza_var = defVar(ds, "sza", Float64, ("sza",))
    sza_var[:] = Float64.(result.sweep.sza)
    sza_var.attrib["units"] = "degrees"
    sza_var.attrib["long_name"] = "solar zenith angle"

    vza_var = defVar(ds, "vza", Float64, ("view_pair",))
    vaz_var = defVar(ds, "vaz", Float64, ("view_pair",))
    vza_var[:] = Float64[first(p) for p in result.sweep.view_pairs]
    vaz_var[:] = Float64[last(p)  for p in result.sweep.view_pairs]
    vza_var.attrib["units"] = "degrees"
    vaz_var.attrib["units"] = "degrees"
    vza_var.attrib["long_name"] = "view zenith angle (paired with vaz)"
    vaz_var.attrib["long_name"] = "view azimuth angle (paired with vza)"

    R = defVar(ds, "R", Float64, ("sza", "view_pair", "brdf", "pol", "spec"))
    T = defVar(ds, "T", Float64, ("sza", "view_pair", "brdf", "pol", "spec"))
    R[:, :, :, :, :] = Float64.(result.R)
    T[:, :, :, :, :] = Float64.(result.T)
    R.attrib["long_name"] = "TOA reflectance"
    T.attrib["long_name"] = "BOA transmittance"

    ds.attrib["scene_path"]   = scene.path
    ds.attrib["scene_indices"] = string(scene.indices)
    ds.attrib["scene_lat"]    = Float64(scene.lat)
    ds.attrib["scene_lon"]    = Float64(scene.lon)
    ds.attrib["scene_time"]   = string(scene.time)
    ds.attrib["vsmartmom_writer"] = "gchp-io PR (Phase D)"
    return nothing
end

# =============================================================================

"""
    generate_benchmark(gchp_path, base_params, sweep::ScenarioSweep, outdir;
                       faces=nothing, xrange=nothing, yrange=nothing,
                       FT=Float64) -> Vector{String}

Top-level scene-loop driver. Opens the GCHP file once, iterates the
requested grid cells, builds a per-scene model + runs the sweep, and
writes one NetCDF-4 per scene under `outdir`.

The benchmark-dataset repo wraps this with its own HDF5/Zarr packaging
on top; here we keep the writer minimal and dependency-free.
"""
function generate_benchmark(gchp_path::AbstractString,
                            base_params::vSmartMOM_Parameters{FT},
                            sweep::ScenarioSweep{FT},
                            outdir::AbstractString;
                            faces = nothing,
                            xrange = nothing,
                            yrange = nothing) where {FT}
    isdir(outdir) || mkpath(outdir)
    output_paths = String[]

    GCHPFile(gchp_path; FT=FT) do f
        face_iter = faces === nothing ? (1:f.nf) : faces
        x_iter    = xrange === nothing ? (1:f.nx) : xrange
        y_iter    = yrange === nothing ? (1:f.ny) : yrange

        for idf in face_iter, idy in y_iter, idx in x_iter
            scene = scene_at(f, idx, idy, idf; FT=FT)
            params = parameters_from_scene(scene, base_params)
            model = model_from_parameters(params)
            so = scene_optics(model, params)
            result = run_sweep(so, sweep)
            out_name = @sprintf("scene_%02d_%03d_%03d.nc4", idf, idy, idx)
            out_path = joinpath(outdir, out_name)
            write_scene_result(out_path, scene, result)
            push!(output_paths, out_path)
        end
    end

    return output_paths
end
