"""
    chain_rule_combine_dτ(dL_dτ, dτ_dp)

Contract f₂ Jacobians with parameter derivatives of layer optical depth.
`dL_dτ` has shape `(nGeom, nStokes, nSpec, nLayers)` and `dτ_dp` has shape
`(nLayers, nSpec, nParam)`. The returned array has shape
`(nGeom, nStokes, nSpec, nParam)`.
"""
function chain_rule_combine_dτ(dL_dτ::AbstractArray{FT1,4},
                               dτ_dp::AbstractArray{FT2,3}) where {FT1,FT2}
    n_geom, n_stokes, n_spec, n_layers = size(dL_dτ)
    n_layers_p, n_spec_p, n_param = size(dτ_dp)
    n_layers == n_layers_p ||
        throw(ArgumentError("layer dimension mismatch: dL_dτ has $n_layers layers, dτ_dp has $n_layers_p"))
    n_spec == n_spec_p ||
        throw(ArgumentError("spectral dimension mismatch: dL_dτ has $n_spec spectra, dτ_dp has $n_spec_p"))

    FT = promote_type(FT1, FT2)
    dL_dp = zeros(FT, n_geom, n_stokes, n_spec, n_param)

    @inbounds for ip in 1:n_param, ispec in 1:n_spec,
                  istokes in 1:n_stokes, ig in 1:n_geom
        acc = zero(FT)
        for iz in 1:n_layers
            acc += convert(FT, dL_dτ[ig, istokes, ispec, iz]) *
                   convert(FT, dτ_dp[iz, ispec, ip])
        end
        dL_dp[ig, istokes, ispec, ip] = acc
    end

    return dL_dp
end

chain_rule_combine_dτ(dL_dτ, dτ_dp) =
    throw(ArgumentError("chain_rule_combine_dτ expects a 4D dL_dτ array and a 3D dτ_dp array"))

"""
    chain_rule_combine_surface_brdf(dL_dρ, dρ_dp)

Contract f₂ Jacobians with parameter derivatives of scalar surface BRDF.
`dL_dρ` has shape `(nGeom, nStokes, nSpec)` and `dρ_dp` has shape
`(nGeom, nSpec, nParam)`. The returned array has shape
`(nGeom, nStokes, nSpec, nParam)`.
"""
function chain_rule_combine_surface_brdf(
    dL_dρ::AbstractArray{FT1,3},
    dρ_dp::AbstractArray{FT2,3},
) where {FT1,FT2}
    n_geom, n_stokes, n_spec = size(dL_dρ)
    n_geom_p, n_spec_p, n_param = size(dρ_dp)
    n_geom == n_geom_p ||
        throw(ArgumentError("geometry dimension mismatch: dL_dρ has " *
                            "$n_geom geometries, dρ_dp has $n_geom_p"))
    n_spec == n_spec_p ||
        throw(ArgumentError("spectral dimension mismatch: dL_dρ has " *
                            "$n_spec spectra, dρ_dp has $n_spec_p"))

    FT = promote_type(FT1, FT2)
    dL_dp = zeros(FT, n_geom, n_stokes, n_spec, n_param)

    @inbounds for ip in 1:n_param, ispec in 1:n_spec,
                  istokes in 1:n_stokes, ig in 1:n_geom
        dL_dp[ig, istokes, ispec, ip] =
            convert(FT, dL_dρ[ig, istokes, ispec]) *
            convert(FT, dρ_dp[ig, ispec, ip])
    end

    return dL_dp
end

chain_rule_combine_surface_brdf(dL_dρ, dρ_dp) =
    throw(ArgumentError("chain_rule_combine_surface_brdf expects a 3D " *
                        "dL_dρ array and a 3D dρ_dp array"))
