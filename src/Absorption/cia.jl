# =============================================================================
# cia.jl — HITRAN Collision-Induced Absorption (CIA) support.
#
# Reads HITRAN .cia files (one or more (T, ν, σ) blocks per pair, e.g. O2-O2,
# N2-N2, O2-N2), pre-interpolates σ onto the model spectral grid per block
# temperature, then adds the per-layer CIA contribution to τ_abs at model
# build time. Promotion of the standalone pilot at
# `test/benchmarks/pilot_cia.jl`.
#
# HITRAN CIA Eq. (3): α(ν, T) = σ_AB(ν, T) · n_A · n_B  [cm⁻¹].
# Layer integration uses the symmetric midpoint-column form
# σ · [n_A(mid)N_B + n_B(mid)N_A]/2. This is algebraically identical in the
# thin-layer limit and avoids the severe midpoint-n² × geometric-thickness
# bias in coarse pressure layers.
# =============================================================================

const _CIA_K_B            = 1.380649e-23   # J/K (Boltzmann)
const _CIA_VMR_O2_DEFAULT = 0.20946        # USS dry-air O2 mixing ratio
const _CIA_VMR_N2_DEFAULT = 0.78084        # USS dry-air N2 mixing ratio

struct CIABlock{FT}
    formula::String
    T::FT
    ν::Vector{FT}
    σ::Vector{FT}     # cm⁵ / molec²
    reference_code::String
end

CIABlock{FT}(formula::String, T::FT, ν::Vector{FT}, σ::Vector{FT}) where {FT} =
    CIABlock{FT}(formula, T, ν, σ, "")

# Before `reference_code` was added, Julia supplied this inferred four-argument
# outer constructor automatically. Keep that public call shape working while
# promoting mixed real input types to one concrete storage type.
function CIABlock(formula::AbstractString, T::Real,
                  ν::AbstractVector{<:Real}, σ::AbstractVector{<:Real})
    FT = float(promote_type(typeof(T), eltype(ν), eltype(σ)))
    return CIABlock{FT}(String(formula), FT(T), FT.(ν), FT.(σ), "")
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
    reference_codes::Vector{String}     # HITRAN CIA header reference codes
    coverage::BitVector                 # any selected source covers ν[k]
    negative_policy::Symbol
end

function CIATable(pair::String, species_a::String, species_b::String,
                  Ts::Vector{Float64}, σ_νT::Matrix{Float64})
    coverage = BitVector(vec(any(isfinite, σ_νT; dims=2)))
    return CIATable(pair, species_a, species_b, Ts, σ_νT,
                    String[], coverage, :error)
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
    cols 98–100  HITRAN CIA reference code
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
            reference_code = length(header) >= 100 ?
                String(strip(header[98:100])) : ""
            push!(blocks, CIABlock{Float64}(
                formula, T_K, νs, σs, reference_code))
        end
    end
    return blocks
end

"""
    build_cia_table(blocks, ν_grid) -> CIATable

Project all blocks onto `ν_grid`, grouping by block temperature. Missing
spectral coverage is stored as `NaN`, not as a physical zero, so subsequent
temperature interpolation can use only panels that actually tabulate a given
wavenumber. Overlapping blocks at the same temperature use the finer native
spectral sampling. Negative interpolated values are rejected by default;
`:clamp_zero` must be requested explicitly. Always Float64 internally (see
`CIATable` docstring).
"""
function build_cia_table(blocks::AbstractVector{<:CIABlock},
                         ν_grid::AbstractVector;
                         negative_policy::Symbol=:error)
    isempty(blocks) && error("build_cia_table: no blocks")
    negative_policy in (:error, :clamp_zero) || throw(ArgumentError(
        "negative_policy must be :error or :clamp_zero"))
    pair_str = blocks[1].formula
    all(block -> block.formula == pair_str, blocks) ||
        error("build_cia_table: all blocks must describe the same CIA pair")
    a, b = _split_pair(pair_str)
    Ts = sort!(unique(Float64[blk.T for blk in blocks]))
    ν_grid_f64 = Float64.(ν_grid)
    σ_νT = fill(NaN, length(ν_grid), length(Ts))
    resolution_νT = fill(Inf, length(ν_grid), length(Ts))
    for (jt, T_pick) in enumerate(Ts)
        for blk in blocks
            blk.T == T_pick || continue
            _accumulate_σ!(view(σ_νT, :, jt),
                           view(resolution_νT, :, jt),
                           ν_grid_f64, blk.ν, blk.σ)
        end
    end
    negative = isfinite.(σ_νT) .& (σ_νT .< 0)
    if any(negative)
        if negative_policy === :error
            error("build_cia_table: selected CIA source produces " *
                  "$(count(negative)) negative cross sections on the model grid")
        end
        σ_νT[negative] .= 0.0
    end
    coverage = BitVector(vec(any(isfinite, σ_νT; dims=2)))
    references = sort!(unique(filter(code -> !isempty(code),
        String[block.reference_code for block in blocks])))
    return CIATable(pair_str, a, b, Ts, σ_νT, references,
                    coverage, negative_policy)
end

# Linear interpolate a single block onto ν_grid (in-place).
@inline function _accumulate_σ!(σ_out::AbstractVector,
                                 resolution_out::AbstractVector,
                                 ν_grid::AbstractVector,
                                 ν_blk::AbstractVector,
                                 σ_blk::AbstractVector)
    νlo, νhi = ν_blk[1], ν_blk[end]
    n = length(ν_blk)
    resolution = n > 1 ? (νhi - νlo) / (n - 1) : Inf
    for (k, νq) in enumerate(ν_grid)
        (νq < νlo || νq > νhi) && continue
        (!isfinite(σ_out[k]) || resolution < resolution_out[k]) || continue
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
        resolution_out[k] = resolution
    end
    return σ_out
end

"""
    cia_σ_at_T!(σ_out, table, T_layer)

Fill `σ_out[k] = σ(ν_grid[k], T_layer)` by linear interpolation in T using
only measured panels that cover that wavenumber. Constant extrapolation uses
the nearest covering temperature. A wavenumber absent from every selected
panel is returned as zero, but this means *outside tabulated support*, not a
demonstrated physical zero; callers can inspect `table.coverage`. `σ_out` must
be a Float64 vector matching the table's ν grid.
"""
function cia_σ_at_T!(σ_out::AbstractVector{Float64},
                     table::CIATable, T_layer::Real)
    Ts = table.Ts
    T = Float64(T_layer)
    @inbounds for k in eachindex(σ_out)
        lower = 0
        upper = 0
        for j in eachindex(Ts)
            isfinite(table.σ_νT[k, j]) || continue
            if Ts[j] <= T
                lower = j
            end
            if Ts[j] >= T
                upper = j
                break
            end
        end
        if lower == 0 && upper == 0
            σ_out[k] = 0.0
        elseif lower == 0
            σ_out[k] = table.σ_νT[k, upper]
        elseif upper == 0 || lower == upper
            σ_out[k] = table.σ_νT[k, lower]
        else
            w = (T - Ts[lower]) / (Ts[upper] - Ts[lower])
            σ_out[k] = (1 - w) * table.σ_νT[k, lower] +
                       w * table.σ_νT[k, upper]
        end
    end
    return σ_out
end

"""
    compute_τ_cia!(τ_abs, table, profile, vmr_dict)

Add the symmetric midpoint-column CIA contribution
`σ(ν,T) · [n_A(mid)N_B + n_B(mid)N_A] / 2` to `τ_abs[ν,layer]`
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
        n_air = Float64(profile.p_full[iz]) * 1e2 / (_CIA_K_B * T) * 1e-6
        r_h2o = Float64(profile.vmr_h2o[iz])
        n_dry = n_air / (1.0 + r_h2o)
        N_dry = Float64(profile.vcd_dry[iz])
        # `vmr_dict` gas values are dry-air mole fractions. H2O is carried
        # separately by AtmosphericProfile. Using each collider's exact
        # hydrostatic column with the other collider's midpoint density makes
        # constant-composition CIA invariant to vertical subdivision. The
        # symmetric average also makes A-B and B-A tables numerically
        # equivalent, including mixed dry-gas/H2O pairs.
        n_a, N_a = _cia_collider_density_column(
            table.species_a, vmr_dict, profile, iz, n_air, n_dry, N_dry)
        n_b, N_b = _cia_collider_density_column(
            table.species_b, vmr_dict, profile, iz, n_air, n_dry, N_dry)
        prod = 0.5 * (n_a * N_b + n_b * N_a)
        for k in 1:nλ
            τ_abs[k, iz] += FT_τ(σ_layer[k] * prod)
        end
    end
    return τ_abs
end

@inline function _cia_collider_density_column(species::AbstractString,
                                                vmr_dict::AbstractDict,
                                                profile, iz::Integer,
                                                n_air::Float64,
                                                n_dry::Float64,
                                                N_dry::Float64)
    if species == "H2O"
        r = Float64(profile.vmr_h2o[iz])
        return r * n_dry, Float64(profile.vcd_h2o[iz])
    end
    v = Float64(_layer_vmr(species, vmr_dict, iz))
    return v * n_dry, v * N_dry
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
    load_cia_table(path, ν_grid; reference_codes=nothing,
                   negative_policy=:error) -> CIATable

Convenience: parse a `.cia` file and build the model-grid table in one call.
`reference_codes` can select one documented HITRAN source family rather than
implicitly stitching overlapping experiments. `negative_policy` is forwarded
to `build_cia_table`.
The `FT` keyword is accepted for back-compat but ignored — the table is
always Float64 (see `CIATable` docstring).
"""
function load_cia_table(path::AbstractString,
                        ν_grid::AbstractVector;
                        FT::Type = Float64,
                        reference_codes=nothing,
                        negative_policy::Symbol=:error)
    blocks = parse_cia_file(path)
    if reference_codes !== nothing
        requested = Set(_normalise_reference_codes(reference_codes))
        blocks = filter(block -> block.reference_code in requested, blocks)
        isempty(blocks) && error(
            "no CIA blocks in $path match reference code(s) " *
            join(sort!(collect(requested)), ","))
    else
        _assert_unambiguous_references(blocks, ν_grid, path)
    end
    return build_cia_table(blocks, ν_grid;
                           negative_policy=negative_policy)
end

function _normalise_reference_codes(reference_codes)
    raw = if reference_codes isa AbstractString
        [reference_codes]
    elseif reference_codes isa AbstractVector || reference_codes isa Tuple ||
           reference_codes isa AbstractSet
        collect(reference_codes)
    else
        throw(ArgumentError(
            "reference_codes must be a nonempty string or collection of strings"))
    end
    isempty(raw) && throw(ArgumentError("reference_codes must not be empty"))
    all(code -> code isa AbstractString, raw) || throw(ArgumentError(
        "reference_codes entries must be strings (for example \"54\")"))
    codes = unique(strip.(String.(raw)))
    any(isempty, codes) && throw(ArgumentError(
        "reference_codes entries must not be empty"))
    return codes
end

"""
    _assert_unambiguous_references(blocks, ν_grid, path)

Fail closed when two HITRAN reference families cover the same requested model
wavenumber and the caller did not explicitly select `reference_codes`. Blocks
from the same reference may overlap (for example, panels with different native
resolution); `build_cia_table` resolves those deterministically.
"""
function _assert_unambiguous_references(blocks, ν_grid, path)
    length(unique(block.reference_code for block in blocks)) <= 1 && return nothing

    for νq in ν_grid
        covering = Set{String}()
        for block in blocks
            νlo, νhi = minmax(first(block.ν), last(block.ν))
            νlo <= νq <= νhi && push!(covering, block.reference_code)
        end
        if length(covering) > 1
            labels = sort!([isempty(code) ? "<missing>" : code
                            for code in covering])
            throw(ArgumentError(
                "load_cia_table: $(path) has overlapping HITRAN CIA " *
                "reference families $(join(labels, ", ")) at ν=$(νq) cm⁻¹; " *
                "specify reference_codes explicitly"))
        end
    end
    return nothing
end
