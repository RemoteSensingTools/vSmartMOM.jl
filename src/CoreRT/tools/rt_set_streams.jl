#=
Quadrature streams for the RT model. Supported types: `GaussLegQuad`,
`RadauQuad`. SZA and VZAs are appended as zero-weight nodes for output
postprocessing (see `docs/src/pages/conventions.md`).
=#

"""
    rt_set_streams(::GaussLegQuad, Ltrunc, obs_geom, pol_type, arr_type)

Half-space Gauss-Legendre quadrature on `[0, 1]` for the upper hemisphere.
For even `Ltrunc` use `Nstreams = (Ltrunc + 2) ÷ 2` (the `+2` form avoids
collapsing to zero streams when `Ltrunc = 0`); for odd `Ltrunc` this reduces
to `(Ltrunc + 1) ÷ 2`. SZA and all VZAs are appended as zero-weight nodes.
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
    iμ₀ = nearest_point(qp_μ, μ₀)

    qp_μN = arr_type(reduce(vcat, fill.(qp_μ, [pol_type.n])))
    wt_μN = arr_type(reduce(vcat, fill.(wt_μ, [pol_type.n])))
    i_start = pol_type.n * (iμ₀ - 1) + 1
    return QuadPoints(μ₀, iμ₀, i_start, arr_type(qp_μ), arr_type(wt_μ), qp_μN, wt_μN, Nquad)
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

    iμ₀ = nearest_point(qp_μ, μ₀);
    qp_μN = arr_type(reduce(vcat, (fill.(qp_μ, [pol_type.n]))))
    wt_μN = arr_type(reduce(vcat, (fill.(wt_μ, [pol_type.n]))))
    i_start = pol_type.n*(iμ₀-1) + 1
    return QuadPoints(μ₀, iμ₀, i_start, arr_type(qp_μ), arr_type(wt_μ), qp_μN, wt_μN, Nquad) 
end
