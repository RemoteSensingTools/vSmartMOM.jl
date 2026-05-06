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

function _selector_path_requested(path::Symbol, paths::Symbol)
    path === :total && return true
    path === :path1 && return _wants_path1(paths)
    path === :path2 && return _wants_path2(paths)
    path === :path3 && return _wants_path3(paths)
    path === :path4 && return _wants_path4(paths)
    return false
end

function _validate_jacobian_selector(selector::SSMeasurementSelector,
                                     paths::Symbol)
    for path in selector.paths
        _selector_path_requested(path, paths) ||
            throw(ArgumentError("SSMeasurementSelector path :$path was not " *
                                "requested by run_exact_ss_with_jacobians " *
                                "paths=:$paths"))
    end
    return nothing
end

function _path1_jacobians(::Val{N}, config::ExactSSConfig, optics, I0) where {N}
    FT = _config_numeric_type(config)
    n_layers, n_spec = size(optics.ϖ_eff)
    n_geom = length(config.geometry.μv)
    dτ = zeros(FT, n_geom, N, n_spec, n_layers)
    dϖ = zeros(FT, n_geom, N, n_spec, n_layers)
    dP = zeros(FT, n_geom, N, n_spec, n_layers)

    @inbounds for ispec in 1:n_spec, istokes in 1:N, iv in 1:n_geom
        μ₀ = convert(FT, config.geometry.μ₀)
        μᵥ = convert(FT, config.geometry.μv[iv])
        a = inv(μ₀) + inv(μᵥ)
        no_a_prefactor = I0[ispec] / (FT(4π) * μᵥ)
        prefactor = no_a_prefactor / a

        for iz in 1:n_layers
            ϖ = optics.ϖ_eff[iz, ispec]
            P = _path1_jacobian_phase(Val(N), optics.P_eff, iv, istokes, iz,
                                       ispec)
            τ_top = optics.τ_cum[iz, ispec]
            τ_bot = optics.τ_cum[iz + 1, ispec]
            layer_factor = exp(-τ_top * a) - exp(-τ_bot * a)

            dϖ[iv, istokes, ispec, iz] += prefactor * P * layer_factor
            dP[iv, istokes, ispec, iz] += prefactor * ϖ * layer_factor

            same_layer = no_a_prefactor * ϖ * P * exp(-τ_bot * a)
            dτ[iv, istokes, ispec, iz] += same_layer

            attenuation_from_above = -no_a_prefactor * ϖ * P * layer_factor
            for i_above in 1:(iz - 1)
                dτ[iv, istokes, ispec, i_above] += attenuation_from_above
            end
        end
    end

    return (; τ_layer=dτ, ϖ_eff=dϖ, P_eff=dP)
end

_path1_jacobian_phase(::Val{1}, P_eff, iv, istokes, iz, ispec) =
    P_eff[iv, iz, ispec]
_path1_jacobian_phase(::Val{N}, P_eff, iv, istokes, iz, ispec) where {N} =
    P_eff[iv, istokes, iz, ispec]

function _path2_jacobians(::Val{N}, config::ExactSSConfig, optics, surface_brdf,
                          I0, path2) where {N}
    FT = _config_numeric_type(config)
    n_layers = size(optics.ϖ_eff, 1)
    n_spec = size(optics.ϖ_eff, 2)
    n_geom = length(config.geometry.μv)
    dτ = zeros(FT, n_geom, N, n_spec, n_layers)
    dρ = zeros(FT, n_geom, N, n_spec)

    @inbounds for ispec in 1:n_spec, istokes in 1:N, iv in 1:n_geom
        μ₀ = convert(FT, config.geometry.μ₀)
        μᵥ = convert(FT, config.geometry.μv[iv])
        τ_total = optics.τ_total_column[ispec]
        attenuation = exp(-τ_total / μ₀) * exp(-τ_total / μᵥ)
        dρ[iv, istokes, ispec] = μ₀ * I0[ispec] * attenuation

        slope = inv(μ₀) + inv(μᵥ)
        for iz in 1:n_layers
            dτ[iv, istokes, ispec, iz] =
                -path2[iv, istokes, ispec] * slope
        end
    end

    return (; τ_layer=dτ, surface_brdf=dρ)
end

function _zero_path1_jacobians(::Type{FT}, n_geom::Int, n_spec::Int,
                               n_layers::Int) where {FT}
    return _zero_path1_jacobians(FT, n_geom, 1, n_spec, n_layers)
end

function _zero_path1_jacobians(::Type{FT}, n_geom::Int, n_stokes::Int,
                               n_spec::Int, n_layers::Int) where {FT}
    z = zeros(FT, n_geom, n_stokes, n_spec, n_layers)
    return (; τ_layer=copy(z), ϖ_eff=copy(z), P_eff=copy(z))
end

function _zero_path2_jacobians(::Type{FT}, n_geom::Int, n_spec::Int,
                               n_layers::Int) where {FT}
    return _zero_path2_jacobians(FT, n_geom, 1, n_spec, n_layers)
end

function _zero_path2_jacobians(::Type{FT}, n_geom::Int, n_stokes::Int,
                               n_spec::Int, n_layers::Int) where {FT}
    return (; τ_layer=zeros(FT, n_geom, n_stokes, n_spec, n_layers),
            surface_brdf=zeros(FT, n_geom, n_stokes, n_spec))
end

function _combine_path_jacobians(path1_jac, path2_jac)
    return (;
        τ_layer = path1_jac.τ_layer .+ path2_jac.τ_layer,
        ϖ_eff = path1_jac.ϖ_eff,
        P_eff = path1_jac.P_eff,
        surface_brdf = path2_jac.surface_brdf)
end

"""
    run_exact_ss_with_jacobians(config; paths=:paths_1_2,
                                selector=SSMeasurementSelector())

Run the standalone solver and return handcoded f₂ Jacobians for paths 1 and 2.
The Jacobians are with respect to the solver seam variables
`τ_layer`, `ϖ_eff`, `P_eff`, and `surface_brdf`, holding the other seam
variables fixed. The returned NamedTuple also includes `measurements`, the
selected retrieval measurement vector from `selected_measurements(result,
selector)`. This Phase 5a prototype is CPU-only.
"""
function run_exact_ss_with_jacobians(config::ExactSSConfig;
                                     paths::Symbol = :paths_1_2,
                                     selector::SSMeasurementSelector =
                                         SSMeasurementSelector())
    _validate_jacobian_paths(paths)
    _validate_jacobian_config(config)
    _validate_jacobian_selector(selector, paths)

    result = run_exact_ss(config; paths)
    FT = _config_numeric_type(config)
    n_stokes = _config_n_stokes(config)
    stokes_dispatch = Val(n_stokes)
    optics = _precompute_optics(stokes_dispatch, config)
    n_layers, n_spec = size(optics.ϖ_eff)
    n_geom = length(config.geometry.μv)
    I0 = _vectorize_I0(config.I0, n_spec, FT)

    path1_jac = _wants_path1(paths) ?
        _path1_jacobians(stokes_dispatch, config, optics, I0) :
        _zero_path1_jacobians(FT, n_geom, n_stokes, n_spec, n_layers)

    if _wants_path2(paths)
        surface_brdf = _precompute_surface_brdf(
            stokes_dispatch, config.surface, config.geometry, n_spec, FT)
        path2_jac = _path2_jacobians(stokes_dispatch, config, optics,
                                     surface_brdf, I0, result.path2)
    else
        path2_jac = _zero_path2_jacobians(FT, n_geom, n_stokes, n_spec,
                                          n_layers)
    end

    jacobians = (; total=_combine_path_jacobians(path1_jac, path2_jac),
                 path1=path1_jac,
                 path2=path2_jac)
    measurements = selected_measurements(result, selector)
    return merge(result, (; jacobians, measurement_selector=selector,
                          measurements))
end
