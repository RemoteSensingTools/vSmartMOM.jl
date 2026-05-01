# =============================================================================
# cia.jl — HITRAN Collision-Induced Absorption (CIA) support.
#
# Reads HITRAN .cia files (one or more (T, ν, σ) blocks per pair, e.g. O2-O2,
# N2-N2, O2-N2), pre-interpolates σ onto the model spectral grid per block
# temperature, then adds the per-layer CIA contribution to τ_abs at model
# build time. Promotion of the standalone pilot at
# `test/benchmarks/pilot_cia.jl`.
#
# HITRAN CIA Eq. (3): α(ν, T) = σ_AB(ν, T) · n_A · n_B  [cm⁻¹], with
# τ_layer = α · Δz. σ has units cm⁵·molec⁻²; n in molec/cm³; Δz in cm.
# =============================================================================

const _CIA_K_B            = 1.380649e-23   # J/K (Boltzmann)
const _CIA_VMR_O2_DEFAULT = 0.20946        # USS dry-air O2 mixing ratio
const _CIA_VMR_N2_DEFAULT = 0.78084        # USS dry-air N2 mixing ratio

struct CIABlock{FT}
    formula::String
    T::FT
    ν::Vector{FT}
    σ::Vector{FT}     # cm⁵ / molec²
end

"""
    CIATable

Pre-interpolated CIA cross-section table keyed by (ν_model_grid, T_block).
Stored in Float64 because σ values are ~1e-44 to 1e-46 cm⁵/molec² — well
below Float32's smallest normal (~1.18e-38). Built once per band at model
construction time; cheap T-interpolation per layer at runtime.
"""
struct CIATable
    pair::String                        # e.g. "O2-O2"
    species_a::String
    species_b::String
    Ts::Vector{Float64}                 # block temperatures, ascending
    σ_νT::Matrix{Float64}               # (length(ν_grid), length(Ts))
end

"""
    parse_cia_file(path) -> Vector{CIABlock{Float64}}

Read a HITRAN .cia file. Each block is a single header line followed by
n_pts lines of (ν σ) pairs. Header layout (fixed-width, 100 cols):
    cols  1–20   formula (e.g. "O2-O2")
    cols 21–30   ν_min (f10.3)
    cols 31–40   ν_max (f10.3)
    cols 41–47   n_pts (i7)
    cols 48–54   T (f7.1, K)
    cols 55–64   σ_max (e10.3)
"""
function parse_cia_file(path::AbstractString)
    blocks = CIABlock{Float64}[]
    open(path) do io
        while !eof(io)
            header = readline(io)
            (length(header) < 54 || isempty(strip(header))) && continue
            formula = String(strip(header[1:20]))
            n_pts   = parse(Int,     strip(header[41:47]))
            T_K     = parse(Float64, strip(header[48:54]))
            νs      = Vector{Float64}(undef, n_pts)
            σs      = Vector{Float64}(undef, n_pts)
            for i in 1:n_pts
                vals = split(strip(readline(io)))
                νs[i] = parse(Float64, vals[1])
                σs[i] = parse(Float64, vals[2])
            end
            push!(blocks, CIABlock{Float64}(formula, T_K, νs, σs))
        end
    end
    return blocks
end

"""
    build_cia_table(blocks, ν_grid) -> CIATable

Project all blocks onto `ν_grid`, grouping by block temperature. Outside
each block's ν-range the contribution is 0. Always Float64 internally
(see `CIATable` docstring).
"""
function build_cia_table(blocks::AbstractVector{<:CIABlock},
                         ν_grid::AbstractVector)
    isempty(blocks) && error("build_cia_table: no blocks")
    pair_str = blocks[1].formula
    a, b = _split_pair(pair_str)
    Ts = sort!(unique(Float64[blk.T for blk in blocks]))
    ν_grid_f64 = Float64.(ν_grid)
    σ_νT = zeros(Float64, length(ν_grid), length(Ts))
    for (jt, T_pick) in enumerate(Ts)
        for blk in blocks
            blk.T == T_pick || continue
            _accumulate_σ!(view(σ_νT, :, jt), ν_grid_f64, blk.ν, blk.σ)
        end
    end
    return CIATable(pair_str, a, b, Ts, σ_νT)
end

# Linear interpolate a single block onto ν_grid (in-place).
@inline function _accumulate_σ!(σ_out::AbstractVector,
                                 ν_grid::AbstractVector,
                                 ν_blk::AbstractVector,
                                 σ_blk::AbstractVector)
    νlo, νhi = ν_blk[1], ν_blk[end]
    n = length(ν_blk)
    for (k, νq) in enumerate(ν_grid)
        (νq < νlo || νq > νhi) && continue
        j = searchsortedfirst(ν_blk, νq)
        if j ≤ 1
            σ_out[k] = σ_blk[1]
        elseif j > n
            σ_out[k] = σ_blk[end]
        else
            ν1, ν2 = ν_blk[j-1], ν_blk[j]
            s1, s2 = σ_blk[j-1], σ_blk[j]
            w = (νq - ν1) / (ν2 - ν1)
            σ_out[k] = (1 - w) * s1 + w * s2
        end
    end
    return σ_out
end

"""
    cia_σ_at_T!(σ_out, table, T_layer)

Fill `σ_out[k] = σ(ν_grid[k], T_layer)` by linear interpolation in T.
Constant extrapolation outside the block-temperature range. `σ_out`
must be a Float64 vector matching the table's ν grid.
"""
function cia_σ_at_T!(σ_out::AbstractVector{Float64},
                     table::CIATable, T_layer::Real)
    Ts = table.Ts
    if T_layer ≤ Ts[1]
        @views σ_out .= table.σ_νT[:, 1]
    elseif T_layer ≥ Ts[end]
        @views σ_out .= table.σ_νT[:, end]
    else
        j = searchsortedfirst(Ts, Float64(T_layer))
        T1, T2 = Ts[j-1], Ts[j]
        w = (Float64(T_layer) - T1) / (T2 - T1)
        @views @. σ_out = (1 - w) * table.σ_νT[:, j-1] + w * table.σ_νT[:, j]
    end
    return σ_out
end

"""
    compute_τ_cia!(τ_abs, table, profile, vmr_dict)

Add the CIA contribution `σ(ν, T) · n_A · n_B · Δz` to `τ_abs[ν, layer]`
for every layer of `profile`. `vmr_dict` provides per-layer (or scalar)
mixing ratios; `O2`/`N2` fall back to USS defaults if absent.

All intermediate arithmetic is Float64 to avoid Float32 underflow on σ
(~1e-45) and overflow on `n²·Δz` (~1e41); the result is converted to
`eltype(τ_abs)` only at accumulation.
"""
function compute_τ_cia!(τ_abs::AbstractMatrix,
                         table::CIATable,
                         profile,
                         vmr_dict::AbstractDict)
    nλ, nlay = size(τ_abs)
    nλ == size(table.σ_νT, 1) ||
        error("compute_τ_cia!: ν-grid length mismatch — table $(size(table.σ_νT, 1)), τ_abs $(nλ)")
    FT_τ = eltype(τ_abs)
    σ_layer = Vector{Float64}(undef, nλ)
    @inbounds for iz in 1:nlay
        T = Float64(profile.T[iz])
        cia_σ_at_T!(σ_layer, table, T)
        # Number density (molec/cm³). p in hPa → ×100 for Pa, then ×1e-6 for m⁻³ → cm⁻³.
        n_air = Float64(profile.p_full[iz]) * 1e2 / (_CIA_K_B * T) * 1e-6
        v_a   = Float64(_layer_vmr(table.species_a, vmr_dict, iz))
        v_b   = Float64(_layer_vmr(table.species_b, vmr_dict, iz))
        Δz_cm = Float64(profile.Δz[iz]) * 100.0   # m → cm
        prod  = v_a * v_b * n_air * n_air * Δz_cm
        for k in 1:nλ
            τ_abs[k, iz] += FT_τ(σ_layer[k] * prod)
        end
    end
    return τ_abs
end

# Helper: parse "O2-O2" / "O2-N2" / "N2-N2" into species pair.
function _split_pair(formula::AbstractString)
    parts = split(strip(formula), '-')
    length(parts) ≥ 2 ||
        error("CIA pair formula \"$formula\" not recognised (expected \"A-B\")")
    return String(parts[1]), String(parts[2])
end

# Helper: layer VMR for a species, with USS fallbacks for the bulk colliders.
@inline function _layer_vmr(species::AbstractString, vmr_dict::AbstractDict, iz::Integer)
    if haskey(vmr_dict, species)
        v = vmr_dict[species]
        return v isa AbstractVector ? v[iz] : v
    elseif species == "O2"
        return _CIA_VMR_O2_DEFAULT
    elseif species == "N2"
        return _CIA_VMR_N2_DEFAULT
    else
        error("CIA: no vmr for \"$species\" and no default defined")
    end
end

"""
    load_cia_table(path, ν_grid) -> CIATable

Convenience: parse a `.cia` file and build the model-grid table in one call.
The `FT` keyword is accepted for back-compat but ignored — the table is
always Float64 (see `CIATable` docstring).
"""
function load_cia_table(path::AbstractString,
                        ν_grid::AbstractVector;
                        FT::Type = Float64)
    blocks = parse_cia_file(path)
    return build_cia_table(blocks, ν_grid)
end
