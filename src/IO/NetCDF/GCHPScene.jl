# GCHP scene reader: open-once / read-many handle + per-cell scenes
#
# Lifecycle model is the inverse of the legacy `geoschem_to_dict` reader
# which opens/closes the NetCDF on every call. `GCHPFile(path; ...) do f`
# loads the lat/lon mesh once, leaves the dataset open, and serves
# `scene_at(f, idx, idy, idf)` cheaply.
#
# Behaviour changes vs. the legacy reader (intentional bug fixes):
#   * `q` is converted from g/kg to kg/kg when the file reports `g kg-1`.
#     The downstream `compute_atmos_profile_fields` treats q as kg/kg
#     (atmo_prof.jl:39); the old reader passed g/kg through unconverted.
#   * H2O is kept in `scene.vmr` but is NEVER placed in
#     `scene_to_dict`'s `absorption.molecules` list. The vSmartMOM
#     parser explicitly rejects H2O in that list (Parameters.jl:918).
#     H2O VMR is still consumed implicitly via `q`.

using ..Aerosols: AerosolScheme, AbstractAerosolScheme, AbstractAerosolBinData,
                  SectionalAerosolData, TOMASScheme
using ..Aerosols: read_aerosol_cell, compute_column_aod

# ============================================================================
# Types
# ============================================================================

"""
    GCHPFile

Open handle to a GCHP/GEOS-Chem NetCDF4 file. Keeps the `NCDataset` open
across many per-cell reads; close with `close(f)` or use the do-block
form to dispose deterministically.

# Fields
- `path::String`: Source file path
- `ds::NCDataset`: Open NetCDF handle
- `nx, ny, nf::Int`: Grid extents (Xdim, Ydim, nf)
- `nlev::Int`: Number of model layers
- `has_tomas::Bool`: Whether TOMAS `SpeciesConcVV_NK01` is present
- `gas_list::Vector{String}`: Discovered trace-gas species (intersection of
  known molecules and `SpeciesConcVV_*` variables)
- `aer_scheme::Union{Nothing, AerosolScheme}`: Aerosol scheme handle
  (optional — set via `GCHPFile(path; aerosol_scheme=...)`)
- `lats, lons::Array{Float64,3}`: `(nx, ny, nf)` location grids, loaded once
- `time::Any`: Timestamp from the file
"""
struct GCHPFile
    path::String
    ds::NCDataset
    nx::Int
    ny::Int
    nf::Int
    nlev::Int
    has_tomas::Bool
    gas_list::Vector{String}
    aer_scheme::Union{Nothing, AerosolScheme}
    lats::Array{Float64,3}
    lons::Array{Float64,3}
    time::Any
end

"""
    GCHPScene{FT}

One scene (single grid cell × single timestep) extracted from a `GCHPFile`.
All profile arrays are oriented TOA→BOA. The `aerosols` field is populated
when the parent file has TOMAS data and an aerosol scheme is configured;
otherwise it is `nothing`.
"""
struct GCHPScene{FT}
    path::String
    indices::NTuple{3,Int}
    lat::FT
    lon::FT
    time::Any
    p_half::Vector{FT}
    T::Vector{FT}
    q::Vector{FT}
    air_mass::Vector{FT}
    air_vol::Vector{FT}
    vmr::Dict{String,Vector{FT}}
    molecules::Vector{String}
    aerosols::Union{Nothing, AbstractAerosolBinData{FT}}
end

# ============================================================================
# Construction / lifecycle
# ============================================================================

const _GCHP_TRACE_GASES = ("N2O", "CH4", "C2H6", "CO2", "CO", "H2O")

"""
    GCHPFile(path; aerosol_scheme=nothing, FT=Float64) -> GCHPFile
    GCHPFile(f, path; ...)

Open a GCHP/GEOS-Chem NetCDF4 file and inspect its grid + variables. The
returned handle is safe to pass to `scene_at` / `scenes` repeatedly; the
underlying dataset stays open until `close(f)` (or end-of-do-block).

`aerosol_scheme` accepts:
- `nothing` (the default) — **no aerosol ingest**. Scenes have `aerosols = nothing`
  even when the file carries TOMAS variables. Use this for atmosphere-only
  pipelines (the legacy `geoschem_to_dict` shim passes `nothing`).
- `:auto` — autodetect from file contents. Currently maps to `:tomas15` when
  `SpeciesConcVV_NK01` is present; otherwise leaves aerosols disabled.
- A `Symbol` naming a concrete TOMAS variant: `:tomas12`, `:tomas15`,
  `:tomas30`, `:tomas40`. Built via [`TOMASScheme(variant)`](@ref).
- A pre-built `AbstractAerosolScheme` instance (full control over species
  selection, RI keys, etc.) — stored as-is.
"""
function GCHPFile(path::AbstractString; aerosol_scheme=nothing, FT::DataType=Float64)
    ds = NCDataset(String(path))
    lats = Array{Float64,3}(ds["lats"][:, :, :])
    lons = Array{Float64,3}(ds["lons"][:, :, :])
    nx, ny, nf = size(lats)
    nlev = length(ds["lev"])
    time = ds["time"][1]

    # TOMAS detection — NK01 is the canonical first bin
    has_tomas = haskey(ds, "SpeciesConcVV_NK01")

    # Discover available trace gases
    gas_list = String[]
    for g in _GCHP_TRACE_GASES
        haskey(ds, "SpeciesConcVV_$g") && push!(gas_list, g)
    end

    # Resolve scheme argument:
    #   nothing  → disabled (no aerosol ingest)
    #   :auto    → autodetect (TOMAS if NK01 present, else disabled)
    #   :tomasNN → concrete TOMAS variant
    #   <scheme> → use as-is
    scheme = if aerosol_scheme === nothing
        nothing
    elseif aerosol_scheme === :auto
        has_tomas ? TOMASScheme(:tomas15; FT=FT) : nothing
    elseif aerosol_scheme isa Symbol
        TOMASScheme(aerosol_scheme; FT=FT)
    else
        aerosol_scheme
    end

    return GCHPFile(String(path), ds, nx, ny, nf, nlev, has_tomas,
                    gas_list, scheme, lats, lons, time)
end

function GCHPFile(f::Function, path::AbstractString; kwargs...)
    handle = GCHPFile(path; kwargs...)
    try
        return f(handle)
    finally
        close(handle)
    end
end

Base.close(f::GCHPFile) = close(f.ds)

Base.show(io::Base.IO, f::GCHPFile) = print(io,
    "GCHPFile(", basename(f.path), "; nx=", f.nx, ", ny=", f.ny,
    ", nf=", f.nf, ", nlev=", f.nlev,
    ", has_tomas=", f.has_tomas,
    ", gases=", length(f.gas_list), ")")

# ============================================================================
# Per-cell read
# ============================================================================

"""
    scene_at(f::GCHPFile, idx, idy, idf; FT=Float64) -> GCHPScene{FT}

Read one cell `(idx, idy, idf)` from the open file. Profile arrays are
flipped from GCHP's BOA→TOA storage to vSmartMOM's TOA→BOA convention.

Aerosols are loaded only when both (a) the file carries TOMAS variables
(`f.has_tomas == true`) **and** (b) the file handle has a concrete
`AbstractAerosolScheme` attached. With `GCHPFile(path; aerosol_scheme=nothing)`
(the default) the `aerosols` field is `nothing`. Pass
`aerosol_scheme = :auto` (or a specific variant like `:tomas15`) to
populate it.
"""
function scene_at(f::GCHPFile, idx::Integer, idy::Integer, idf::Integer;
                  FT::DataType=Float64)
    1 ≤ idx ≤ f.nx || throw(BoundsError(f, (idx, idy, idf)))
    1 ≤ idy ≤ f.ny || throw(BoundsError(f, (idx, idy, idf)))
    1 ≤ idf ≤ f.nf || throw(BoundsError(f, (idx, idy, idf)))

    ds  = f.ds
    lat = FT(f.lats[idx, idy, idf])
    lon = FT(f.lons[idx, idy, idf])

    # Pressure half-grid (TOA→BOA), hPa. Met_DELP layer thickness in hPa;
    # Met_PS2WET surface pressure in hPa. GCHP stores BOA→TOA, so the half
    # grid built from sp .+ cumsum(-dp) goes from surface up; reverse it.
    dp = ds["Met_DELP"].var[idx, idy, idf, :, 1]
    sp = ds["Met_PS2WET"].var[idx, idy, idf, 1]
    p_half = reverse([sp; sp .+ cumsum(-dp)])

    # Temperature (K), TOA→BOA
    _require_geoschem(ds["Met_T"].attrib["units"] == "K",
                      "Met_T units expected to be 'K'")
    T = reverse(ds["Met_T"].var[idx, idy, idf, :, 1])

    # Specific humidity — convert g/kg → kg/kg when needed
    q_units = ds["Met_SPHU"].attrib["units"]
    q_raw   = reverse(ds["Met_SPHU"].var[idx, idy, idf, :, 1])
    q = if q_units == "g kg-1"
        q_raw ./ 1000
    elseif q_units == "kg kg-1"
        q_raw
    else
        _geoschem_error("Unsupported Met_SPHU units: '$(q_units)' (expected 'g kg-1' or 'kg kg-1')")
    end

    # Dry-air mass per layer (kg) and grid-box volume per layer (m³).
    # Both are consumed by the TOMAS aerosol unit conversion (see
    # `read_aerosol_cell`); kept on the scene so non-aerosol callers can
    # still introspect the column.
    air_mass = haskey(ds, "Met_AD") ?
        reverse(ds["Met_AD"].var[idx, idy, idf, :, 1]) :
        fill(NaN, f.nlev)
    air_vol = haskey(ds, "Met_AIRVOL") ?
        reverse(ds["Met_AIRVOL"].var[idx, idy, idf, :, 1]) :
        fill(NaN, f.nlev)

    # Trace-gas VMRs (mol/mol dry), TOA→BOA. H2O is kept here for callers
    # that want it explicitly; `scene_to_dict` will exclude it from the
    # absorption.molecules list (the parser rejects H2O there).
    vmr = Dict{String, Vector{FT}}()
    for g in f.gas_list
        var_name = "SpeciesConcVV_$g"
        _require_geoschem(ds[var_name].attrib["units"] == "mol mol-1 dry",
                          "$(var_name) units expected to be 'mol mol-1 dry'")
        vmr[g] = FT.(reverse(ds[var_name].var[idx, idy, idf, :, 1]))
    end

    # Molecules list for absorption — exclude H2O (auto-included via q)
    molecules = String[g for g in f.gas_list if g != "H2O"]

    # Aerosols (only when a TOMAS scheme has been configured AND the file
    # carries the expected TOMAS variables). `f.aer_scheme === nothing`
    # means the caller opted out of aerosol ingest — honour it.
    aerosols = if f.has_tomas && f.aer_scheme isa AbstractAerosolScheme
        read_aerosol_cell(f.aer_scheme, ds, idx, idy, idf)
    else
        nothing
    end

    return GCHPScene{FT}(f.path, (Int(idx), Int(idy), Int(idf)),
                         lat, lon, f.time,
                         FT.(p_half), FT.(T), FT.(q),
                         FT.(air_mass), FT.(air_vol),
                         vmr, molecules, aerosols)
end

Base.show(io::Base.IO, s::GCHPScene{FT}) where {FT} = print(io,
    "GCHPScene{", FT, "}(idx=", s.indices,
    ", lat=", round(s.lat; digits=2), "°, lon=", round(s.lon; digits=2), "°, ",
    "nlev=", length(s.T), ", gases=", length(s.vmr),
    ", aerosols=", isnothing(s.aerosols) ? "nothing" : nameof(typeof(s.aerosols)),
    ")")

"""
    compute_scene_aod(scene::GCHPScene, wavelengths_um; kwargs...)

Convenience wrapper around `Aerosols.compute_column_aod` for a GCHP scene.
The scene must carry `scene.aerosols`; open the file with
`GCHPFile(path; aerosol_scheme=:auto)` or an explicit TOMAS scheme.
"""
function compute_scene_aod(scene::GCHPScene, wavelengths_um::AbstractVector; kwargs...)
    scene.aerosols === nothing &&
        throw(ArgumentError("scene has no aerosol data; open GCHPFile with aerosol_scheme=:auto or an explicit scheme"))
    return compute_column_aod(scene.aerosols, scene.p_half, scene.T, scene.q,
                              wavelengths_um; kwargs...)
end

# ============================================================================
# Iteration
# ============================================================================

"""
    scenes(f::GCHPFile; faces=1:f.nf, xrange=1:f.nx, yrange=1:f.ny, FT=Float64)

Lazy iterator yielding one `GCHPScene` per `(idx, idy, idf)` in the
specified ranges. The dataset is read once per scene; the file handle
stays open across the whole iteration.
"""
function scenes(f::GCHPFile;
                faces = 1:f.nf,
                xrange = 1:f.nx,
                yrange = 1:f.ny,
                FT::DataType = Float64)
    return (scene_at(f, idx, idy, idf; FT=FT)
            for idf in faces, idy in yrange, idx in xrange)
end

"""
    read_gchp_scene(path, idx, idy, idf; FT=Float64, aerosol_scheme=nothing)

One-shot open/read/close helper — equivalent to:

```julia
GCHPFile(path; aerosol_scheme=aerosol_scheme, FT=FT) do f
    scene_at(f, idx, idy, idf; FT=FT)
end
```

Prefer the do-block form when reading many cells from the same file.
"""
function read_gchp_scene(path::AbstractString,
                         idx::Integer, idy::Integer, idf::Integer;
                         FT::DataType=Float64, aerosol_scheme=nothing)
    return GCHPFile(path; aerosol_scheme=aerosol_scheme, FT=FT) do f
        scene_at(f, idx, idy, idf; FT=FT)
    end
end

# ============================================================================
# Conversion / cloning APIs
# ============================================================================

"""
    scene_to_dict(scene::GCHPScene) -> Dict

Produce a configuration `Dict` compatible with `parameters_from_dict`.
This is the new-path equivalent of the legacy `geoschem_to_dict`, with
two intentional behaviour fixes:

- `q` is in kg/kg (the legacy dict left g/kg unconverted)
- H2O is excluded from `absorption.molecules` (the parser rejects it)

The H2O VMR remains accessible via `scene.vmr["H2O"]` for callers that
want it explicitly.
"""
function scene_to_dict(scene::GCHPScene{FT}) where {FT}
    cfg = Dict{String, Any}()

    cfg["atmospheric_profile"] = Dict{String, Any}(
        "T" => Float64.(scene.T),
        "p" => Float64.(scene.p_half),
        "q" => Float64.(scene.q),
        "vmr" => Dict{String, Vector{Float64}}(k => Float64.(v) for (k, v) in scene.vmr),
        "profile_reduction" => nothing,
    )

    if !isempty(scene.molecules)
        # Outer vector wraps per-band (single band here). We list ONLY the
        # molecules for which the scene provides a VMR profile —
        # `parameters_from_dict` requires a VMR for every listed species
        # (Parameters.jl:_parse_absorption). The legacy reader injected
        # `"O2"` without an O2 VMR, which caused the dict to fail validation.
        # Callers who want O2 absorption supply it via their base config
        # and merge with `parameters_from_scene`.
        all_molecules = [collect(scene.molecules)]
        # Filter VMR dict to molecules actually listed (excludes H2O even
        # if present in scene.vmr — required by the parser).
        abs_vmr = Dict{String, Vector{Float64}}(
            k => Float64.(scene.vmr[k]) for k in scene.molecules)
        cfg["absorption"] = Dict{String, Any}(
            "molecules"    => all_molecules,
            "vmr"          => abs_vmr,
            "broadening"   => "Voigt()",
            "CEF"          => "HumlicekWeidemann32SDErrorFunction()",
            "wing_cutoff"  => 25,
        )
    end

    cfg["_metadata"] = Dict{String, Any}(
        "source_file"   => scene.path,
        "source_type"   => "GEOSChem",
        "latitude"      => Float64(scene.lat),
        "longitude"     => Float64(scene.lon),
        "time"          => scene.time,
        "grid_indices"  => scene.indices,
    )

    return cfg
end

"""
    parameters_from_scene(scene, base::vSmartMOM_Parameters) -> vSmartMOM_Parameters

Clone `base` and overwrite its atmospheric profile + absorption VMRs from
`scene`. Returns a new `vSmartMOM_Parameters` that does NOT share mutable
state with `base` (verified by `test/test_parameters_from_scene_alias.jl`).

`base` provides the static configuration (geometry, quadrature, spec_bands,
surfaces, aerosol scattering setup); `scene` provides the per-scene
atmosphere. This is the canonical entry point for the scene-loop driver.
"""
function parameters_from_scene(scene::GCHPScene{FT},
                               base::vSmartMOM_Parameters) where {FT}
    params = deepcopy(base)
    PFT = params.float_type

    params.T = PFT.(scene.T)
    params.p = PFT.(scene.p_half)
    params.q = PFT.(scene.q)

    if params.absorption_params !== nothing && !isempty(scene.molecules)
        # Replace VMRs for molecules present in the scene; leave any
        # molecules `base` carried but the scene doesn't provide untouched.
        new_vmr = deepcopy(params.absorption_params.vmr)
        for mol in scene.molecules
            if haskey(scene.vmr, mol)
                new_vmr[mol] = PFT.(scene.vmr[mol])
            end
        end
        params.absorption_params.vmr = new_vmr
    end

    return params
end
