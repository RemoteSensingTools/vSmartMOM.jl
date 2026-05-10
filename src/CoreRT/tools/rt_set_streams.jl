#=
Quadrature streams for the RT model. Supported types: `GaussLegQuad`,
`RadauQuad`. SZA and VZAs are appended as zero-weight nodes for output
postprocessing (see `docs/src/pages/conventions.md`).

Two stream counts live on the resulting `QuadPoints`:
- `Nstreams` — count of nonzero weights, the user-facing resolving-power
  knob. Public contract: `stream_l_cap = 2·Nstreams - 1`. For Gauss this
  matches `(Ltrunc + 2) ÷ 2` (Sanghavi: `+2` avoids the `Ltrunc=0` corner).
  For Radau the value is implementation-derived because the SZA-not-on-node
  branch builds `2·Nquad_subinterval` weighted nodes around μ₀.
- `Nquad` — total node count including zero-weight SZA/VZA output nodes,
  used for kernel sizing.
=#

"""
    rt_set_streams(::GaussLegQuad, Ltrunc, obs_geom, pol_type, arr_type)

Half-space Gauss-Legendre quadrature on `[0, 1]` for the upper hemisphere.
Builds `(Ltrunc + 2) ÷ 2` weighted nodes (Sanghavi: the `+2` avoids
collapsing to zero streams when `Ltrunc = 0`). SZA and all VZAs are
appended as zero-weight output nodes; `Nstreams` records the count of
nonzero weights, `Nquad` records the augmented total.
"""
function rt_set_streams(::GaussLegQuad,
                        Ltrunc::Int,
                        obs_geom::ObsGeometry{FT},
                        pol_type,
                        arr_type) where {FT}
    (; sza, vza) = obs_geom
    Nquad = (Ltrunc + 2) ÷ 2

    qp_μ, wt_μ = Scattering.gauleg(Nquad, zero(FT), one(FT))
    μ₀ = cosd.(sza)

    # Append SZA and VZA cosines as zero-weight nodes
    qp_μ = unique(FT[qp_μ; cosd.(vza); μ₀])
    wt_μ = FT[wt_μ; zeros(FT, length(qp_μ) - length(wt_μ))]
    Nquad = length(qp_μ)
    Nstreams = count(!iszero, wt_μ)
    iμ₀ = nearest_point(qp_μ, μ₀)

    qp_μN = arr_type(reduce(vcat, fill.(qp_μ, [pol_type.n])))
    wt_μN = arr_type(reduce(vcat, fill.(wt_μ, [pol_type.n])))
    i_start = pol_type.n * (iμ₀ - 1) + 1
    return QuadPoints(μ₀, iμ₀, i_start, arr_type(qp_μ), arr_type(wt_μ), qp_μN, wt_μN, Nquad, Nstreams)
end

"""
$(FUNCTIONNAME)(::RadauQuad, 
                        Ltrunc::Int, 
                        obs_geom::ObsGeometry, 
                        pol_type,
                        arr_type)

Computes hemispheric quadrature points with Gauss-Radau method in two intervals between [0,SZA] and [SZA,1] to include SZA as full quadrature point for DNI
SZA included as full weighted node, VZAs included with 0 weight
Returns computed quadrature points as [`QuadPoints`](@ref) 
"""
function rt_set_streams(::RadauQuad, 
                        Ltrunc::Int, 
                        obs_geom::ObsGeometry{FT}, 
                        pol_type,
                        arr_type) where {FT}
    
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
    
    if μ₀ ∈ qp_μ₀

        qp_μ = zeros(FT, Nquad)
        wt_μ = zeros(FT, Nquad)

        for i = 1:Nquad
            qp_μ[i] = (1 + qp_μ₀[i]) / 2
            wt_μ[i] = wt_μ₀[i] 
        end

        Ncam = length(vza)
        μ = cosd.(vza) # check for degree to radian conversion
        
        # Screen out duplicate camera zenith angles
        qp_μ = unique([qp_μ; cosd.(vza)])

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

        Ncam = length(vza)
        μ = cosd.(vza) # check for degree to radian conversion

        # Screen out duplicate camera zenith angles
        qp_μ = unique([qp_μ; cosd.(vza)])

        # Assign zero-weights to remaining camera zenith angles
        wt_μ = FT[wt_μ; zeros(length(qp_μ) - length(wt_μ))]
        Nquad = length(qp_μ)

    end

    Nstreams = count(!iszero, wt_μ)
    iμ₀ = nearest_point(qp_μ, μ₀);
    qp_μN = arr_type(reduce(vcat, (fill.(qp_μ, [pol_type.n]))))
    wt_μN = arr_type(reduce(vcat, (fill.(wt_μ, [pol_type.n]))))
    i_start = pol_type.n*(iμ₀-1) + 1
    return QuadPoints(μ₀, iμ₀, i_start, arr_type(qp_μ), arr_type(wt_μ), qp_μN, wt_μN, Nquad, Nstreams)
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
    rt_set_streams(::GaussLegQuad, obs_geom, pol_type, arr_type; nstreams::Int)

Build a half-space Gauss-Legendre quadrature with exactly `nstreams`
weighted nodes per hemisphere. Equivalent to the positional form with
`Ltrunc = 2·nstreams - 2` (so `(Ltrunc+2)÷2 == nstreams` matches the
inverse of Sanghavi's even-`Ltrunc` formula).
"""
function rt_set_streams(q::GaussLegQuad,
                        obs_geom::ObsGeometry,
                        pol_type,
                        arr_type;
                        nstreams::Int)
    nstreams >= 1 || throw(ArgumentError("nstreams must be ≥ 1; got $nstreams"))
    Ltrunc = 2 * nstreams - 2
    return rt_set_streams(q, Ltrunc, obs_geom, pol_type, arr_type)
end

"""
    rt_set_streams(::RadauQuad, obs_geom, pol_type, arr_type; nstreams::Int)

Build a Gauss-Radau quadrature whose public-contract `stream_l_cap` is
`2·nstreams - 1`. The actual weighted-stream count `qp.Nstreams` may
exceed `nstreams` in the typical SZA-not-on-node branch (doubled around
μ₀); the public contract `stream_l_cap = 2·nstreams - 1` still holds.
Per Sanghavi (2026-05-07), `RadauQuad` is supported but not recommended
— `GaussLegQuad` is the default and the cheaper choice for the same
resolving power.
"""
function rt_set_streams(q::RadauQuad,
                        obs_geom::ObsGeometry,
                        pol_type,
                        arr_type;
                        nstreams::Int)
    nstreams >= 1 || throw(ArgumentError("nstreams must be ≥ 1; got $nstreams"))
    Ltrunc = 2 * nstreams - 1
    return rt_set_streams(q, Ltrunc, obs_geom, pol_type, arr_type)
end
