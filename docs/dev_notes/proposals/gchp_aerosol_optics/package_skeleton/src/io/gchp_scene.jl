# GCHP column reader.
#
# PORT TARGET: gchp-io:src/IO/NetCDF/GCHPScene.jl already implements GCHPFile / GCHPScene /
# scene_at / scenes with an open-once-read-many design and TOMAS-15 ingestion. Move that
# file here and replace internal `..Aerosols` / `..Scattering` references with `vSmartMOM`.
#
# This stub only records the intended surface so the rest of the package can compile
# against it once ported.

"""
    GCHPFile

Handle to an open GCHP NetCDF file (cubed-sphere). Open once, serve many `scene_at` calls.
Carries grid metadata (`nx, ny, nf, nlev`), `has_tomas`, the gas list, and an optional
aerosol scheme. See verified file layout in STUDENT_BRIEF §2.
"""
struct GCHPFile
    # TODO (port from gchp-io): path, ds::NCDataset, nx/ny/nf/nlev, has_tomas,
    # gas_list, aer_scheme, lats, lons, time
end

"""
    GCHPScene{FT}

One extracted column: `p_half, T, q, air_mass (Met_AD), air_vol (Met_AIRVOL), vmr`, and
`aerosols::Union{Nothing,SectionalAerosolData{FT}}` (per-layer × per-bin × per-species).
"""
struct GCHPScene{FT}
    # TODO (port from gchp-io): indices, lat, lon, time, p_half, T, q,
    # air_mass, air_vol, vmr, molecules, aerosols
end

"""
    GCHPFile(path; aerosol_scheme=nothing, FT=Float64) -> GCHPFile

`aerosol_scheme`: `nothing` (skip), `:auto` (detect TOMAS via `NK01`), a `Symbol`
(`:tomas15`), or a prebuilt scheme object.
"""
function GCHPFile(path::AbstractString; aerosol_scheme=nothing, FT::Type=Float64)
    error("not implemented — port gchp-io:src/IO/NetCDF/GCHPScene.jl")
end

"""    scene_at(f::GCHPFile, ix, iy, iface; FT=Float64) -> GCHPScene{FT}"""
function scene_at end

"""    scenes(f::GCHPFile; faces, xrange, yrange, FT=Float64) -> iterator of GCHPScene"""
function scenes end

Base.close(f::GCHPFile) = nothing  # TODO: close underlying NCDataset
