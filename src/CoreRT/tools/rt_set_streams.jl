#=
Quadrature streams for the RT model. Supported types: `GaussLegQuad`,
`RadauQuad`. VZAs are appended as zero-weight output nodes. In the default
Gauss representation (`external_solar=false`), SZA is also appended as a
zero-weight operator node. In the opt-in external-solar path
(`external_solar=true`, TOA-only), scalar `μ₀` remains outside the
diffuse operator and is appended only to the phase-evaluation grid, where it
provides the exact direct-beam source column (see
`docs/src/pages/conventions.md`). Radau always embeds μ₀ as a weighted node.
Appended ordinates (VZAs and an embedded μ₀) within `4·eps(1)` of a node are
merged into it, and a node that close to μ₀ is moved onto the exact μ₀
(`_build_output_grid`); weighted nodes are not merged with each other.

Two stream counts live on the resulting `QuadPoints`:
- `Nstreams` — count of nonzero weights, the user-facing resolving-power
  knob. Public contract: `stream_l_cap = 2·Nstreams - 1`. For Gauss this
  matches `(Ltrunc + 2) ÷ 2` (Sanghavi: `+2` avoids the `Ltrunc=0` corner).
  For Radau the value is implementation-derived because the SZA-not-on-node
  branch builds `2·Nquad_subinterval` weighted nodes around μ₀.
- `Nquad` — total diffuse-operator node count, including zero-weight VZAs and,
  in the default embedded Gauss representation, the zero-weight SZA. The
  external-solar phase grid is tracked separately by `phase_qp_μ` and does not
  enlarge the square diffuse operators.
=#

# Ordinates closer than `_NODE_MERGE_ULPS · eps(1)` (absolute, on μ ∈ [0, 1])
# are treated as one direction. Distinct but numerically equal ordinates are
# harmful: kernels with 1/(μᵢ − μⱼ) or 1/(μᵢ − μ₀) cancellation lose all
# precision, so a μ₀ or VZA one rounding error from a node would give a
# visibly different answer from both the on-node and the clearly-off cases
# (e.g. `sza = acosd(node)`, whose cosine is not exactly `node`).
const _NODE_MERGE_ULPS = 4

_nodes_coincide(a, b) = abs(a - b) <=
    _NODE_MERGE_ULPS * max(eps(one(typeof(a))), eps(one(typeof(b))))

"""
    _build_output_grid(FT, nodes, μ_vza, μ₀; append_μ₀) -> Vector{FT}

Append the VZA cosines `μ_vza` (and, when `append_μ₀`, the solar cosine `μ₀`)
to the quadrature `nodes`, so that no two distinct ordinates of the result
coincide (`_nodes_coincide`) and `μ₀` is kept exact:

1. unless `μ₀` is already a node, a node coinciding with `μ₀` is moved onto
   `μ₀`, keeping its weight, so the kernels see `μᵢ == μ₀` exactly, as for an
   SZA exactly on the node, instead of a rounding-error difference that their
   `1/(μᵢ − μ₀)` terms cannot resolve (the node moves by at most
   `_NODE_MERGE_ULPS · eps(1)`);
2. each VZA coinciding with `μ₀` is represented by `μ₀` itself; a VZA
   coinciding with a node already present is skipped (served by that node
   through `nearest_point`); any other VZA is appended;
3. `μ₀` is appended if requested and not yet present.

Doing μ₀ first means a moved node cannot end up next to an already-appended
VZA. Exact duplicates are skipped as `unique` did and the order is unchanged,
so the grid is identical to `unique([nodes; μ_vza; μ₀])` unless a
near-coincidence occurs.
"""
function _build_output_grid(::Type{FT}, nodes, μ_vza, μ₀::FT;
                            append_μ₀::Bool) where {FT}
    grid = collect(FT, nodes)
    if !(μ₀ in grid)
        i = findfirst(y -> _nodes_coincide(y, μ₀), grid)
        i === nothing || (grid[i] = μ₀)
    end
    for v in μ_vza
        v = FT(v)
        μ = _nodes_coincide(v, μ₀) ? μ₀ : v
        any(y -> _nodes_coincide(y, μ), grid) || push!(grid, μ)
    end
    append_μ₀ && !(μ₀ in grid) && push!(grid, μ₀)
    return grid
end

"""
    rt_set_streams(::GaussLegQuad, Ltrunc, obs_geom, pol_type, arr_type;
                   external_solar=false)

Half-space Gauss-Legendre quadrature on `[0, 1]` for the upper hemisphere.
Builds `(Ltrunc + 2) ÷ 2` weighted nodes (Sanghavi: the `+2` avoids
collapsing to zero streams when `Ltrunc = 0`). All VZAs are appended as
zero-weight output nodes. The default (`external_solar=false`) also appends
SZA to the square operator grid — the historical representation.

With `external_solar=true`, the collimated solar direction is instead a scalar
source geometry plus the exact `Zₘ(μᵢ,μ₀)` column on `phase_qp_μ`; it is not a
diffuse stream. In this mode `iμ₀ == iμ₀Nstart == 0` are deliberate sentinels,
and source code must use `iμ₀_phase`. For example, five weighted streams plus
one distinct VZA with `Stokes_IQU` gives `Nquad=6` and an `18×18` diffuse
operator, rather than the legacy `21×21` operator that also embeds SZA.
"""
function rt_set_streams(::GaussLegQuad,
                        Ltrunc::Int,
                        obs_geom::ObsGeometry{FT},
                        pol_type,
                        arr_type;
                        external_solar::Bool = false) where {FT}
    (; sza, vza) = obs_geom
    Nquad = (Ltrunc + 2) ÷ 2

    qp_μ, wt_μ = Scattering.gauleg(Nquad, zero(FT), one(FT))
    μ₀ = cosd.(sza)

    # VZAs are diffuse-operator output nodes. In the default (embedded-μ₀)
    # representation, the solar direction is also an operator node. The
    # opt-in external-solar representation keeps μ₀ only on the phase grid
    # so the square diffuse operators do not carry its zero-weight row and
    # column.
    # Ordinates within a few ulps of each other are merged, keeping μ₀ exact
    # (`_build_output_grid`), so the kernels see the on-node geometry.
    qp_μ = _build_output_grid(FT, qp_μ, cosd.(vza), μ₀;
                              append_μ₀ = !external_solar)
    wt_μ = FT[wt_μ; zeros(FT, length(qp_μ) - length(wt_μ))]
    Nquad = length(qp_μ)
    Nstreams = count(!iszero, wt_μ)

    phase_qp_μ = external_solar ? unique(FT[qp_μ; μ₀]) : qp_μ
    iμ₀_phase = nearest_point(phase_qp_μ, μ₀)
    iμ₀Nstart_phase = pol_type.n * (iμ₀_phase - 1) + 1

    # These legacy indices address the square operator grid. A zero sentinel
    # prevents external-solar consumers from accidentally treating μ₀ as
    # an operator node; source-column consumers use the phase-grid indices.
    iμ₀ = external_solar ? 0 : iμ₀_phase
    i_start = external_solar ? 0 : iμ₀Nstart_phase

    qp_μN = arr_type(reduce(vcat, fill.(qp_μ, [pol_type.n])))
    wt_μN = arr_type(reduce(vcat, fill.(wt_μ, [pol_type.n])))
    return QuadPoints(μ₀, iμ₀, i_start, arr_type(qp_μ), arr_type(wt_μ),
                      qp_μN, wt_μN, Nquad, Nstreams, arr_type(phase_qp_μ),
                      iμ₀_phase, iμ₀Nstart_phase, external_solar)
end

"""
    _merge_coincident_nodes(qp_μ, wt_μ)

Merge exactly coincident quadrature ordinates, summing their weights; the
first-occurrence order of `qp_μ` is kept.
"""
function _merge_coincident_nodes(qp_μ::AbstractVector, wt_μ::AbstractVector)
    allunique(qp_μ) && return qp_μ, wt_μ
    nodes = unique(qp_μ)
    weights = [sum(wt_μ[qp_μ .== μ]) for μ in nodes]
    return nodes, weights
end

"""
$(FUNCTIONNAME)(::RadauQuad, 
                        Ltrunc::Int, 
                        obs_geom::ObsGeometry, 
                        pol_type,
                        arr_type;
                        external_solar=false)

Computes hemispheric quadrature points with Gauss-Radau method in two intervals between [0,SZA] and [SZA,1] to include SZA as full quadrature point for DNI
SZA included as full weighted node, VZAs included with 0 weight
`external_solar=true` is rejected because Radau's defining construction embeds
the solar direction as a weighted ordinate.
Returns computed quadrature points as [`QuadPoints`](@ref) 
"""
function rt_set_streams(::RadauQuad, 
                        Ltrunc::Int, 
                        obs_geom::ObsGeometry{FT}, 
                        pol_type,
                        arr_type;
                        external_solar::Bool = false) where {FT}
    external_solar && throw(ArgumentError(
        "external_solar=true is supported only for GaussLegQuad; RadauQuad includes the solar direction in its weighted quadrature"))
    
    (; obs_alt, sza, vza, vaz) = obs_geom
    
    # Ltrunc + 1 = number of spherical coefficients considered (+1 for l=0)
    # quadtype = 'r' for Radau (with DNI), 'g' for Gauss (no DNI)
    # sza = single solar zenith angle
    # lza = array of or single line-of-sight zenith angle
    Nquad = (Ltrunc + 1) ÷ 2
    
    tqp_μ₀, twt_μ₀ = gaussradau(Nquad)
    qp_μ₀ = -reverse(tqp_μ₀)
    wt_μ₀ =  reverse(twt_μ₀)
    μ₀ = cosd(sza) # check for degree to radian conversion

    # μ₀ is tested against the nodes MAPPED to [0,1] (the fixed Radau node
    # maps to μ = 1, i.e. SZA = 0).
    qp_μ₀_mapped = (1 .+ qp_μ₀) ./ 2

    # A μ₀ within a few ulps of a node counts as on the node; that node is
    # then moved onto the exact μ₀ (`_build_output_grid`).
    if any(y -> _nodes_coincide(y, μ₀), qp_μ₀_mapped)

        qp_μ = zeros(FT, Nquad)
        wt_μ = zeros(FT, Nquad)

        # Map [-1,1] → [0,1]: μ = (1 + x)/2, dμ = dx/2 (as in the split branch).
        for i = 1:Nquad
            qp_μ[i] = qp_μ₀_mapped[i]
            wt_μ[i] = wt_μ₀[i] / 2
        end

        Ncam = length(vza)
        μ = cosd.(vza) # check for degree to radian conversion
        
        # Screen out duplicate (and few-ulp-coincident) camera zenith angles;
        # μ₀ is already a node of this grid (or moved onto it).
        qp_μ = _build_output_grid(FT, qp_μ, cosd.(vza), μ₀; append_μ₀ = false)

        # Assign zero-weights to remaining camera zenith angles
        wt_μ = FT[wt_μ; zeros(length(qp_μ) - length(wt_μ))]
        
        Nquad = length(qp_μ)
        
    else

        qp_μ = zeros(FT, 2Nquad)
        wt_μ = zeros(FT, 2Nquad)

        for i = 1:Nquad
            qp_μ[i] = (μ₀ + μ₀ * qp_μ₀[i]) / 2
            wt_μ[i] = μ₀ * wt_μ₀[i] / 2
            qp_μ[Nquad + i] = ((1 + μ₀) + (1 - μ₀) * qp_μ₀[i]) / 2
            wt_μ[Nquad + i] = (1 - μ₀) * wt_μ₀[i] / 2
        end
        # Near-overhead sun: the sub-interval [μ₀, 1] can be narrower than the
        # FT spacing of its nodes, which then coincide (with each other, with
        # μ₀ or with 1). Merge coincident ordinates, summing their weights, so
        # each ordinate appears once and ∑w = 1 is kept.
        qp_μ, wt_μ = _merge_coincident_nodes(qp_μ, wt_μ)

        Ncam = length(vza)
        μ = cosd.(vza) # check for degree to radian conversion

        # Screen out duplicate (and few-ulp-coincident) camera zenith angles;
        # μ₀ is already a node of this grid (or moved onto it).
        qp_μ = _build_output_grid(FT, qp_μ, cosd.(vza), μ₀; append_μ₀ = false)

        # Assign zero-weights to remaining camera zenith angles
        wt_μ = FT[wt_μ; zeros(length(qp_μ) - length(wt_μ))]
        Nquad = length(qp_μ)

    end

    Nstreams = count(!iszero, wt_μ)
    iμ₀ = nearest_point(qp_μ, μ₀);
    qp_μN = arr_type(reduce(vcat, (fill.(qp_μ, [pol_type.n]))))
    wt_μN = arr_type(reduce(vcat, (fill.(wt_μ, [pol_type.n]))))
    i_start = pol_type.n*(iμ₀-1) + 1
    return QuadPoints(μ₀, iμ₀, i_start, arr_type(qp_μ), arr_type(wt_μ), qp_μN,
                      wt_μN, Nquad, Nstreams, arr_type(qp_μ), iμ₀, i_start, false)
end

# ---------------------------------------------------------------------------
# Public `nstreams` API (kwarg form — preferred from Phase D onwards)
# ---------------------------------------------------------------------------
#
# Users set `radiative_transfer.nstreams = N` (weighted streams per
# hemisphere). The builder derives the equivalent legacy `Ltrunc` per
# scheme so the existing positional-`Ltrunc` machinery keeps producing
# the same node arrays. Public contract: `stream_l_cap = 2·N - 1`.
#
# The kwarg form disambiguates the second positional `Int` from the
# existing `Ltrunc` API; both methods coexist for one release while
# legacy YAMLs (which still carry `l_trunc`) migrate.

"""
    rt_set_streams(::GaussLegQuad, obs_geom, pol_type, arr_type;
                   nstreams::Int, external_solar=false)

Build a half-space Gauss-Legendre quadrature with exactly `nstreams`
weighted nodes per hemisphere. Equivalent to the positional form with
`Ltrunc = 2·nstreams - 2` (so `(Ltrunc+2)÷2 == nstreams` matches the
inverse of Sanghavi's even-`Ltrunc` formula).

By default, the direct solar direction remains outside the diffuse operator
while retaining its exact phase-matrix source column. Pass
`external_solar=false` explicitly for the legacy embedded-μ₀ representation.
"""
function rt_set_streams(q::GaussLegQuad,
                        obs_geom::ObsGeometry,
                        pol_type,
                        arr_type;
                        nstreams::Int,
                        external_solar::Bool = false)
    nstreams >= 1 || throw(ArgumentError("nstreams must be ≥ 1; got $nstreams"))
    Ltrunc = 2 * nstreams - 2
    return rt_set_streams(q, Ltrunc, obs_geom, pol_type, arr_type;
                          external_solar = external_solar)
end

"""
    rt_set_streams(::RadauQuad, obs_geom, pol_type, arr_type;
                   nstreams::Int, external_solar=false)

Build a Gauss-Radau quadrature whose public-contract `stream_l_cap` is
`2·nstreams - 1`. The actual weighted-stream count `qp.Nstreams` may
exceed `nstreams` in the typical SZA-not-on-node branch (doubled around
μ₀); the public contract `stream_l_cap = 2·nstreams - 1` still holds.
Per Sanghavi (2026-05-07), `RadauQuad` is supported but not recommended
— `GaussLegQuad` is the default and the cheaper choice for the same
resolving power.
`external_solar=true` is unsupported; use `GaussLegQuad` for that
representation.
"""
function rt_set_streams(q::RadauQuad,
                        obs_geom::ObsGeometry,
                        pol_type,
                        arr_type;
                        nstreams::Int,
                        external_solar::Bool = false)
    nstreams >= 1 || throw(ArgumentError("nstreams must be ≥ 1; got $nstreams"))
    external_solar && throw(ArgumentError(
        "external_solar=true is supported only for GaussLegQuad; RadauQuad includes the solar direction in its weighted quadrature"))
    Ltrunc = 2 * nstreams - 1
    return rt_set_streams(q, Ltrunc, obs_geom, pol_type, arr_type)
end
