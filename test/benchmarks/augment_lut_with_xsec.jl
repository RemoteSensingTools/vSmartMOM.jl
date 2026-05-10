# =============================================================================
# augment_lut_with_xsec.jl — Merge HITRAN broadband cross-sections (.xsc)
# into an existing line-based LUT (JLD2), or produce a xsec-only LUT for
# species with no HITRAN line catalogue entry (e.g. BrO, ClO₂/OClO).
#
# Motivation: HITRAN .par line files do not cover UV-Vis photodissociation
# bands of species like NO2, SO2, H2CO, BrO. Those live in the HITRAN
# cross-section database. This script bakes the .xsc data onto the LUT's
# (ν, p, T) grid by linear interpolation in (ν, T), replicated across the
# pressure axis (xsec files tabulate σ at zero pressure only, i.e. the
# unbroadened band shape; pressure dependence for these broad electronic
# absorbers is weak).
#
# If the target species has a line-based LUT at `{species}.jld2`, its σ is
# added to the xsec σ. If not (xsec-only species), an all-zero line table
# is created on the standard (p, T, ν) grid.
#
# HITRAN XSC header (100 fixed-width cols; see
# https://hitran.org/docs/cross-sections-definitions/ ):
#   1–20   molecule name
#  21–30   ν_min    (f10.4, cm⁻¹)
#  31–40   ν_max    (f10.4, cm⁻¹)
#  41–47   N_pts    (i7)
#  48–54   T        (f7.1, K)
#  55–60   p        (f6.1, Torr)
#  61–70   σ_max    (e10.3, cm²/molec)
#  71–75   res      (f5.3)
#  76–95   common name
#  96–100  ref id
# Data: 10 σ values per line in e10.4 format.
#
# Usage:
#   julia --project=test test/benchmarks/augment_lut_with_xsec.jl SO2
# =============================================================================

using vSmartMOM
using vSmartMOM.Absorption
using Interpolations
using JLD2
using Printf

const LUT_DIR = expanduser("~/data/HITRAN_LUTs")

# Grids for xsec-only LUTs (match build_hitran_luts.jl exactly).
const XSEC_ONLY_Δν_CM = 0.01
const XSEC_ONLY_ν_LO  = 1e7 / 3500.0
const XSEC_ONLY_ν_HI  = 1e7 /  280.0
const XSEC_ONLY_ν_GRID = XSEC_ONLY_ν_LO:XSEC_ONLY_Δν_CM:XSEC_ONLY_ν_HI
const XSEC_ONLY_P_GRID = range(0.01, 1080.01; length = 15)
const XSEC_ONLY_T_GRID = range(180.0, 360.0;  length = 10)

# Sentinel mol IDs for species that are not in HITRAN's line catalogue.
# Purely informational; vSmartMOM RT does not dispatch on these.
const XSEC_ONLY_MOL_ID = Dict(
    "BrO"  => 900,  # not in HITRAN .par
    "OClO" => 901,  # HITRAN main list is ambiguous; xsec only for our range
    "ClO"  => 902,  # for when ClO xsec .files are dropped in
)

struct XscBlock
    origin::String   # "HITRAN" or "MPI-Mainz"
    T::Float64
    p_Torr::Float64
    ν::Vector{Float64}
    σ::Vector{Float64}
    source::String
end

"""
    parse_xsc_file(path) -> XscBlock

Read one HITRAN .xsc file. Reconstructs the ν grid from the header
(ν_min, ν_max, N_pts) and reads all σ values following the header.
"""
function parse_xsc_file(path::AbstractString)
    open(path) do io
        header = readline(io)
        length(header) >= 60 || error("xsc header too short in $path")

        ν_min = parse(Float64, strip(header[21:30]))
        ν_max = parse(Float64, strip(header[31:40]))
        n_pts = parse(Int,     strip(header[41:47]))
        T     = parse(Float64, strip(header[48:54]))
        # p field exists but not used here; a few newer files omit it.
        p     = try; parse(Float64, strip(header[55:60])); catch; 0.0; end

        σ = Vector{Float64}(undef, n_pts)
        k = 0
        for line in eachline(io)
            for tok in split(strip(line))
                k += 1
                k > n_pts && break
                σ[k] = parse(Float64, tok)
            end
            k >= n_pts && break
        end
        # Some HITRAN .xsc files end with an extra sentinel value; tolerate it.
        k >= n_pts || error("xsc $path: expected $n_pts values, got $k")

        ν = range(ν_min, ν_max; length = n_pts)
        return XscBlock("HITRAN", T, p, collect(ν), σ, basename(path))
    end
end

"""
    linear_interp(ν_src, σ_src, ν_query) -> Vector

Simple linear interpolator. Returns 0.0 outside the source ν range.
"""
function linear_interp(ν_src::AbstractVector, σ_src::AbstractVector,
                       ν_query::AbstractVector)
    out = zeros(eltype(σ_src), length(ν_query))
    for (k, ν) in enumerate(ν_query)
        if ν < ν_src[1] || ν > ν_src[end]
            continue
        end
        j = searchsortedfirst(ν_src, ν)
        if j == 1
            out[k] = σ_src[1]
        elseif j > length(ν_src)
            out[k] = σ_src[end]
        else
            ν1, ν2 = ν_src[j-1], ν_src[j]
            s1, s2 = σ_src[j-1], σ_src[j]
            w = (ν - ν1) / (ν2 - ν1)
            out[k] = (1 - w) * s1 + w * s2
        end
    end
    return out
end

"""
    parse_mpi_mainz_file(path) -> XscBlock

Parse an MPI-Mainz UV/Vis spectral atlas-style tab file. Assumes two
whitespace-separated columns: λ (nm) and σ (cm²/molec). Temperature is
extracted from the filename via a `_<T>K_` token (e.g. `ClO_Simon1990_300K_...`).
"""
function parse_mpi_mainz_file(path::AbstractString)
    fname = basename(path)
    m = match(r"_(\d+(?:\.\d+)?)K_", fname)
    m === nothing && error("Cannot parse temperature from filename: $fname  (expected a _<T>K_ token)")
    T = parse(Float64, m.captures[1])

    λ_nm = Float64[]
    σ    = Float64[]
    open(path) do io
        for line in eachline(io)
            s = strip(line)
            (isempty(s) || startswith(s, '#') || startswith(s, '!')) && continue
            toks = split(s)
            length(toks) < 2 && continue
            λv = tryparse(Float64, toks[1])
            σv = tryparse(Float64, toks[2])
            (λv === nothing || σv === nothing) && continue
            push!(λ_nm, λv); push!(σ, σv)
        end
    end
    # Convert λ (nm, ascending) → ν (cm⁻¹, descending → reverse to ascending).
    ν = reverse(1e7 ./ λ_nm)
    σ = reverse(σ)
    return XscBlock("MPI-Mainz", T, 0.0, ν, σ, fname)
end

"""
    build_xsec_table(xsc_dir, ν_grid, t_grid) -> Matrix(ν × T)

Reads every .xsc and .txt file in `xsc_dir`, groups by temperature, and
linearly interpolates onto `ν_grid × t_grid`. `.xsc` = HITRAN fixed-width
format; `.txt` = MPI-Mainz-style λ(nm) σ(cm²/molec) columns.
"""
function build_xsec_table(xsc_dir::AbstractString,
                          ν_grid::AbstractRange, t_grid::AbstractRange)
    xsc_paths = filter(endswith(".xsc"), readdir(xsc_dir; join = true))
    # Only pick up .txt files that match the MPI-Mainz naming convention
    # (contain a _<T>K_ token). This skips readmes, Emacs backup files, etc.
    txt_paths = filter(readdir(xsc_dir; join = true)) do p
        endswith(p, ".txt") && !endswith(p, "~") && occursin(r"_\d+(?:\.\d+)?K_", basename(p))
    end
    isempty(xsc_paths) && isempty(txt_paths) &&
        error("No .xsc or .txt xsec files found in $xsc_dir")
    blocks = vcat([parse_xsc_file(p) for p in xsc_paths],
                  [parse_mpi_mainz_file(p) for p in txt_paths])
    @info "Parsed $(length(blocks)) xsc blocks"
    Ts  = sort!(unique([b.T for b in blocks]))
    @info "Xsc temperatures available" Ts

    # σ on (ν, T_xsc) grid: accumulate per-T by stitching ν-ranges.
    σ_νT = zeros(Float32, length(ν_grid), length(Ts))
    for (jt, T_xsc) in enumerate(Ts)
        for b in blocks
            b.T == T_xsc || continue
            contrib = linear_interp(b.ν, b.σ, collect(ν_grid))
            # In overlap regions, prefer the last-parsed block (good enough for
            # our files, since HITRAN provides non-overlapping ν windows per T).
            mask = contrib .!= 0
            @views σ_νT[mask, jt] = Float32.(contrib[mask])
        end
    end

    # Interpolate (or extrapolate by nearest) onto the target t_grid.
    σ_out = zeros(Float32, length(ν_grid), length(t_grid))
    for (jt, T_out) in enumerate(t_grid)
        if T_out ≤ Ts[1]
            σ_out[:, jt] .= σ_νT[:, 1]
        elseif T_out ≥ Ts[end]
            σ_out[:, jt] .= σ_νT[:, end]
        else
            j = searchsortedfirst(Ts, T_out)
            T1, T2 = Ts[j-1], Ts[j]
            w = Float32((T_out - T1) / (T2 - T1))
            @views σ_out[:, jt] = (1 - w) .* σ_νT[:, j - 1] .+ w .* σ_νT[:, j]
        end
    end
    return σ_out, Ts
end

"""
    augment!(species_name) — main entry

1. Load `{LUT_DIR}/{species_name}.jld2`.
2. Build xsec σ(ν, T) from `{LUT_DIR}/{species_name}/*.xsc`.
3. Back up the line-only LUT to `{species_name}_lines_only.jld2`.
4. Add xsec σ to every pressure slice of the existing itp.coefs
   (xsec is p-independent at the T-column densities we care about).
5. Rebuild the BSpline(Linear()) interpolator and save over the
   existing file.
"""
function augment!(species_name::AbstractString)
    lut_path   = joinpath(LUT_DIR, species_name * ".jld2")
    xsc_dir    = joinpath(LUT_DIR, species_name)
    backup     = joinpath(LUT_DIR, species_name * "_lines_only.jld2")

    isdir(xsc_dir) || error("No xsc directory at $xsc_dir")

    # Line-LUT present → augment; absent → create xsec-only LUT.
    local m
    if isfile(lut_path)
        @info "Loading line LUT" lut_path
        m = Absorption.load_interpolation_model(lut_path)
        if !isfile(backup)
            @info "Backing up line-only LUT" backup
            cp(lut_path, backup; force = false)
        end
    else
        mol_id = get(XSEC_ONLY_MOL_ID, species_name, 999)
        @info "No line LUT → building xsec-only LUT" species_name mol_id
        zero_coefs = zeros(Float32, length(XSEC_ONLY_ν_GRID),
                                    length(XSEC_ONLY_P_GRID),
                                    length(XSEC_ONLY_T_GRID))
        itp0 = interpolate(zero_coefs, BSpline(Linear()))
        m = Absorption.InterpolationModel(itp0, mol_id, 1,
                                          XSEC_ONLY_ν_GRID,
                                          XSEC_ONLY_P_GRID,
                                          XSEC_ONLY_T_GRID)
    end
    @info "  model" mol = m.mol iso = m.iso ν_range = (first(m.ν_grid), last(m.ν_grid)) size_coefs = size(m.itp.coefs)

    σ_xsec, Ts_xsc = build_xsec_table(xsc_dir, m.ν_grid, m.t_grid)
    @info "Xsec table built" size_νT = size(σ_xsec) max_σ = maximum(σ_xsec)

    # Broadcast-add σ_xsec (ν, T) to each pressure slice of itp.coefs (ν, p, T).
    coefs = Array(m.itp.coefs)  # concrete 3-D array
    @info "  line-only coefs extrema" min = minimum(coefs) max = maximum(coefs)
    @inbounds for jt in axes(coefs, 3)
        for jp in axes(coefs, 2)
            coefs[:, jp, jt] .+= σ_xsec[:, jt]
        end
    end
    @info "  combined coefs extrema"  min = minimum(coefs) max = maximum(coefs)

    # Rebuild interpolator with BSpline(Linear()).
    itp = interpolate(coefs, BSpline(Linear()))
    new_model = Absorption.InterpolationModel(itp, m.mol, m.iso, m.ν_grid, m.p_grid, m.t_grid)
    save_interpolation_model(new_model, lut_path)
    @info "Saved augmented LUT" lut_path size_MB = round(filesize(lut_path) / 1024^2; digits = 1)
    return new_model
end

if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) < 1
        println("usage: julia --project=test test/benchmarks/augment_lut_with_xsec.jl <SPECIES>")
        exit(1)
    end
    augment!(ARGS[1])
end
