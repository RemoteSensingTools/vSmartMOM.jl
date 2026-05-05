const _SUPPORTED_JACOBIAN_PATHS = (:path1, :path2, :paths_1_2)

function _validate_jacobian_paths(paths::Symbol)
    paths in _SUPPORTED_JACOBIAN_PATHS && return paths
    throw(ArgumentError("StandaloneSS f2 Jacobians currently support :path1, :path2, and :paths_1_2; got :$paths"))
end

function _validate_jacobian_config(config::ExactSSConfig)
    config.architecture isa CPU ||
        throw(ArgumentError("StandaloneSS f2 Jacobians currently support CPU() configs only"))
    return nothing
end

function _path1_jacobians(config::ExactSSConfig, optics, I0)
    FT = _config_numeric_type(config)
    n_layers, n_spec = size(optics.ϖ_eff)
    n_geom = length(config.geometry.μv)
    dτ = zeros(FT, n_geom, 1, n_spec, n_layers)
    dϖ = zeros(FT, n_geom, 1, n_spec, n_layers)
    dP = zeros(FT, n_geom, 1, n_spec, n_layers)

    @inbounds for ispec in 1:n_spec, iv in 1:n_geom
        μ₀ = convert(FT, config.geometry.μ₀)
        μᵥ = convert(FT, config.geometry.μv[iv])
        a = inv(μ₀) + inv(μᵥ)
        no_a_prefactor = I0[ispec] / (FT(4π) * μᵥ)
        prefactor = no_a_prefactor / a

        for iz in 1:n_layers
            ϖ = optics.ϖ_eff[iz, ispec]
            P = optics.P_eff[iv, iz, ispec]
            τ_top = optics.τ_cum[iz, ispec]
            τ_bot = optics.τ_cum[iz + 1, ispec]
            layer_factor = exp(-τ_top * a) - exp(-τ_bot * a)

            dϖ[iv, 1, ispec, iz] += prefactor * P * layer_factor
            dP[iv, 1, ispec, iz] += prefactor * ϖ * layer_factor

            same_layer = no_a_prefactor * ϖ * P * exp(-τ_bot * a)
            dτ[iv, 1, ispec, iz] += same_layer

            attenuation_from_above = -no_a_prefactor * ϖ * P * layer_factor
            for i_above in 1:(iz - 1)
                dτ[iv, 1, ispec, i_above] += attenuation_from_above
            end
        end
    end

    return (; τ_layer=dτ, ϖ_eff=dϖ, P_eff=dP)
end

function _path2_jacobians(config::ExactSSConfig, optics, surface_brdf, I0,
                          path2)
    FT = _config_numeric_type(config)
    n_layers = size(optics.ϖ_eff, 1)
    n_spec = size(optics.ϖ_eff, 2)
    n_geom = length(config.geometry.μv)
    dτ = zeros(FT, n_geom, 1, n_spec, n_layers)
    dρ = zeros(FT, n_geom, 1, n_spec)

    @inbounds for ispec in 1:n_spec, iv in 1:n_geom
        μ₀ = convert(FT, config.geometry.μ₀)
        μᵥ = convert(FT, config.geometry.μv[iv])
        τ_total = optics.τ_total_column[ispec]
        attenuation = exp(-τ_total / μ₀) * exp(-τ_total / μᵥ)
        dρ[iv, 1, ispec] = μ₀ * I0[ispec] * attenuation

        slope = inv(μ₀) + inv(μᵥ)
        for iz in 1:n_layers
            dτ[iv, 1, ispec, iz] = -path2[iv, 1, ispec] * slope
        end
    end

    return (; τ_layer=dτ, surface_brdf=dρ)
end

function _zero_path1_jacobians(::Type{FT}, n_geom::Int, n_spec::Int,
                               n_layers::Int) where {FT}
    z = zeros(FT, n_geom, 1, n_spec, n_layers)
    return (; τ_layer=copy(z), ϖ_eff=copy(z), P_eff=copy(z))
end

function _zero_path2_jacobians(::Type{FT}, n_geom::Int, n_spec::Int,
                               n_layers::Int) where {FT}
    return (; τ_layer=zeros(FT, n_geom, 1, n_spec, n_layers),
            surface_brdf=zeros(FT, n_geom, 1, n_spec))
end

function _combine_path_jacobians(path1_jac, path2_jac)
    return (;
        τ_layer = path1_jac.τ_layer .+ path2_jac.τ_layer,
        ϖ_eff = path1_jac.ϖ_eff,
        P_eff = path1_jac.P_eff,
        surface_brdf = path2_jac.surface_brdf)
end

"""
    run_exact_ss_with_jacobians(config; paths=:paths_1_2)

Run the standalone solver and return handcoded f₂ Jacobians for paths 1 and 2.
The Jacobians are with respect to the solver seam variables
`τ_layer`, `ϖ_eff`, `P_eff`, and `surface_brdf`, holding the other seam
variables fixed. This Phase 5a prototype is CPU-only.
"""
function run_exact_ss_with_jacobians(config::ExactSSConfig;
                                     paths::Symbol = :paths_1_2)
    _validate_jacobian_paths(paths)
    _validate_jacobian_config(config)

    result = run_exact_ss(config; paths)
    FT = _config_numeric_type(config)
    optics = _precompute_optics(config)
    n_layers, n_spec = size(optics.ϖ_eff)
    n_geom = length(config.geometry.μv)
    I0 = _vectorize_I0(config.I0, n_spec, FT)

    path1_jac = _wants_path1(paths) ?
        _path1_jacobians(config, optics, I0) :
        _zero_path1_jacobians(FT, n_geom, n_spec, n_layers)

    if _wants_path2(paths)
        surface_brdf = _precompute_surface_brdf(config.surface, config.geometry,
                                                n_spec, FT)
        path2_jac = _path2_jacobians(config, optics, surface_brdf, I0,
                                     result.path2)
    else
        path2_jac = _zero_path2_jacobians(FT, n_geom, n_spec, n_layers)
    end

    jacobians = (; total=_combine_path_jacobians(path1_jac, path2_jac),
                 path1=path1_jac,
                 path2=path2_jac)
    return merge(result, (; jacobians))
end
