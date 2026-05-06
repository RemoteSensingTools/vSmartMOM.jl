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

chain_rule_combine_dτ(dL_dτ::AbstractArray{<:Any,4},
                      dτ_dp::AbstractArray{<:Any,3},
                      selector::SSMeasurementSelector) =
    selected_measurement_jacobian(chain_rule_combine_dτ(dL_dτ, dτ_dp),
                                  selector)

"""
    chain_rule_combine_dϖ(dL_dϖ, dϖ_dp)

Contract f₂ Jacobians with parameter derivatives of effective
single-scattering albedo. `dL_dϖ` has shape
`(nGeom, nStokes, nSpec, nLayers)` and `dϖ_dp` has shape
`(nLayers, nSpec, nParam)`. The returned array has shape
`(nGeom, nStokes, nSpec, nParam)`.
"""
function chain_rule_combine_dϖ(dL_dϖ::AbstractArray{FT1,4},
                               dϖ_dp::AbstractArray{FT2,3}) where {FT1,FT2}
    n_geom, n_stokes, n_spec, n_layers = size(dL_dϖ)
    n_layers_p, n_spec_p, n_param = size(dϖ_dp)
    n_layers == n_layers_p ||
        throw(ArgumentError("layer dimension mismatch: dL_dϖ has $n_layers layers, dϖ_dp has $n_layers_p"))
    n_spec == n_spec_p ||
        throw(ArgumentError("spectral dimension mismatch: dL_dϖ has $n_spec spectra, dϖ_dp has $n_spec_p"))

    FT = promote_type(FT1, FT2)
    dL_dp = zeros(FT, n_geom, n_stokes, n_spec, n_param)

    @inbounds for ip in 1:n_param, ispec in 1:n_spec,
                  istokes in 1:n_stokes, ig in 1:n_geom
        acc = zero(FT)
        for iz in 1:n_layers
            acc += convert(FT, dL_dϖ[ig, istokes, ispec, iz]) *
                   convert(FT, dϖ_dp[iz, ispec, ip])
        end
        dL_dp[ig, istokes, ispec, ip] = acc
    end

    return dL_dp
end

chain_rule_combine_dϖ(dL_dϖ, dϖ_dp) =
    throw(ArgumentError("chain_rule_combine_dϖ expects a 4D dL_dϖ array and a 3D dϖ_dp array"))

chain_rule_combine_dϖ(dL_dϖ::AbstractArray{<:Any,4},
                      dϖ_dp::AbstractArray{<:Any,3},
                      selector::SSMeasurementSelector) =
    selected_measurement_jacobian(chain_rule_combine_dϖ(dL_dϖ, dϖ_dp),
                                  selector)

"""
    chain_rule_combine_dP(dL_dP, dP_dp)

Contract f₂ Jacobians with parameter derivatives of geometry-dependent
effective phase values. `dL_dP` has shape
`(nGeom, nStokes, nSpec, nLayers)` and `dP_dp` has shape
`(nGeom, nLayers, nSpec, nParam)`. The returned array has shape
`(nGeom, nStokes, nSpec, nParam)`.
"""
function chain_rule_combine_dP(dL_dP::AbstractArray{FT1,4},
                               dP_dp::AbstractArray{FT2,4}) where {FT1,FT2}
    n_geom, n_stokes, n_spec, n_layers = size(dL_dP)
    n_geom_p, n_layers_p, n_spec_p, n_param = size(dP_dp)
    n_geom == n_geom_p ||
        throw(ArgumentError("geometry dimension mismatch: dL_dP has " *
                            "$n_geom geometries, dP_dp has $n_geom_p"))
    n_layers == n_layers_p ||
        throw(ArgumentError("layer dimension mismatch: dL_dP has " *
                            "$n_layers layers, dP_dp has $n_layers_p"))
    n_spec == n_spec_p ||
        throw(ArgumentError("spectral dimension mismatch: dL_dP has " *
                            "$n_spec spectra, dP_dp has $n_spec_p"))

    FT = promote_type(FT1, FT2)
    dL_dp = zeros(FT, n_geom, n_stokes, n_spec, n_param)

    @inbounds for ip in 1:n_param, ispec in 1:n_spec,
                  istokes in 1:n_stokes, ig in 1:n_geom
        acc = zero(FT)
        for iz in 1:n_layers
            acc += convert(FT, dL_dP[ig, istokes, ispec, iz]) *
                   convert(FT, dP_dp[ig, iz, ispec, ip])
        end
        dL_dp[ig, istokes, ispec, ip] = acc
    end

    return dL_dp
end

chain_rule_combine_dP(dL_dP, dP_dp) =
    throw(ArgumentError("chain_rule_combine_dP expects a 4D dL_dP array and a 4D dP_dp array"))

chain_rule_combine_dP(dL_dP::AbstractArray{<:Any,4},
                      dP_dp::AbstractArray{<:Any,4},
                      selector::SSMeasurementSelector) =
    selected_measurement_jacobian(chain_rule_combine_dP(dL_dP, dP_dp),
                                  selector)

"""
    chain_rule_combine_surface_brdf(dL_dρ, dρ_dp)

Contract f₂ Jacobians with parameter derivatives of surface BRDF. `dL_dρ`
has shape `(nGeom, nStokes, nSpec)`. For scalar BRDF derivatives, `dρ_dp`
has shape `(nGeom, nSpec, nParam)`. For vector-Stokes BRDF derivatives,
`dρ_dp` has shape `(nGeom, nStokes, nSpec, nParam)`. The returned array has
shape `(nGeom, nStokes, nSpec, nParam)`.
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

function chain_rule_combine_surface_brdf(
    dL_dρ::AbstractArray{FT1,3},
    dρ_dp::AbstractArray{FT2,4},
) where {FT1,FT2}
    n_geom, n_stokes, n_spec = size(dL_dρ)
    n_geom_p, n_stokes_p, n_spec_p, n_param = size(dρ_dp)
    n_geom == n_geom_p ||
        throw(ArgumentError("geometry dimension mismatch: dL_dρ has " *
                            "$n_geom geometries, dρ_dp has $n_geom_p"))
    n_stokes == n_stokes_p ||
        throw(ArgumentError("Stokes dimension mismatch: dL_dρ has " *
                            "$n_stokes components, dρ_dp has $n_stokes_p"))
    n_spec == n_spec_p ||
        throw(ArgumentError("spectral dimension mismatch: dL_dρ has " *
                            "$n_spec spectra, dρ_dp has $n_spec_p"))

    FT = promote_type(FT1, FT2)
    dL_dp = zeros(FT, n_geom, n_stokes, n_spec, n_param)

    @inbounds for ip in 1:n_param, ispec in 1:n_spec,
                  istokes in 1:n_stokes, ig in 1:n_geom
        dL_dp[ig, istokes, ispec, ip] =
            convert(FT, dL_dρ[ig, istokes, ispec]) *
            convert(FT, dρ_dp[ig, istokes, ispec, ip])
    end

    return dL_dp
end

chain_rule_combine_surface_brdf(dL_dρ, dρ_dp) =
    throw(ArgumentError("chain_rule_combine_surface_brdf expects a 3D " *
                        "dL_dρ array and either a 3D scalar or 4D " *
                        "vector-Stokes dρ_dp array"))

chain_rule_combine_surface_brdf(dL_dρ::AbstractArray{<:Any,3},
                                dρ_dp::AbstractArray{<:Any,3},
                                selector::SSMeasurementSelector) =
    selected_measurement_jacobian(
        chain_rule_combine_surface_brdf(dL_dρ, dρ_dp),
        selector)

chain_rule_combine_surface_brdf(dL_dρ::AbstractArray{<:Any,3},
                                dρ_dp::AbstractArray{<:Any,4},
                                selector::SSMeasurementSelector) =
    selected_measurement_jacobian(
        chain_rule_combine_surface_brdf(dL_dρ, dρ_dp),
        selector)
