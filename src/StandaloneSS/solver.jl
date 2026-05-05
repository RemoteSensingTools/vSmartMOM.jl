const _SUPPORTED_PATHS = (:path1, :path2, :path3, :path4, :paths_1_2,
                          :all, :all_four)

function _validate_paths(paths::Symbol)
    paths in _SUPPORTED_PATHS && return paths
    throw(ArgumentError("StandaloneSS supports paths :path1, :path2, :path3, :path4, :paths_1_2, :all, and :all_four; got :$paths"))
end

_wants_path1(paths::Symbol) = paths in (:path1, :paths_1_2, :all, :all_four)
_wants_path2(paths::Symbol) = paths in (:path2, :paths_1_2, :all, :all_four)
_wants_path3(paths::Symbol) = paths in (:path3, :all, :all_four)
_wants_path4(paths::Symbol) = paths in (:path4, :all, :all_four)

function _validate_geometry(geometry::SSGeometry)
    geometry.μ₀ > zero(geometry.μ₀) ||
        throw(ArgumentError("SSGeometry.μ₀ must be positive"))
    all(>(zero(eltype(geometry.μv))), geometry.μv) ||
        throw(ArgumentError("all SSGeometry.μv entries must be positive"))
    length(geometry.μv) == length(geometry.Δϕ) ||
        throw(ArgumentError("SSGeometry.μv and SSGeometry.Δϕ must have the same length"))
end

function _dims(contributors)
    isempty(contributors) && throw(ArgumentError("ExactSSConfig requires at least one contributor"))
    first_size = size(τ_matrix(first(contributors)))
    length(first_size) == 2 ||
        throw(ArgumentError("contributor optical depths must be matrices with shape (nLayer, nSpec)"))
    for c in contributors
        size(τ_matrix(c)) == first_size ||
            throw(ArgumentError("all contributor optical-depth matrices must have the same shape"))
    end
    return first_size
end

function _vectorize_albedo(surface::LambertianSSSurface, n_spec::Int, ::Type{FT}) where {FT}
    a = surface.albedo
    if a isa Number
        return fill(convert(FT, a), n_spec)
    end
    length(a) == n_spec ||
        throw(ArgumentError("LambertianSSSurface.albedo vector length must match nSpec"))
    return convert(Vector{FT}, collect(a))
end

function _vectorize_I0(I0, n_spec::Int, ::Type{FT}) where {FT}
    if I0 isa Number
        return fill(convert(FT, I0), n_spec)
    end
    length(I0) == n_spec ||
        throw(ArgumentError("ExactSSConfig.I0 length must match nSpec"))
    return convert(Vector{FT}, collect(I0))
end

function _scattering_angle_cosine(μ₀, μᵥ, Δϕ)
    FT = promote_type(typeof(μ₀), typeof(μᵥ), typeof(Δϕ))
    μ₀ = convert(FT, μ₀)
    μᵥ = convert(FT, μᵥ)
    Δϕ = convert(FT, Δϕ)
    sin₀ = sqrt(max(zero(FT), one(FT) - μ₀^2))
    sinᵥ = sqrt(max(zero(FT), one(FT) - μᵥ^2))
    return -μ₀ * μᵥ + sin₀ * sinᵥ * cos(Δϕ)
end

function _gauss_legendre_01(n::Int, ::Type{FT}) where {FT}
    n > 0 || throw(ArgumentError("inner_nquad must be positive"))
    x, w = gausslegendre(n)
    return convert(Vector{FT}, (x .+ 1) ./ 2),
           convert(Vector{FT}, w ./ 2)
end

function _azimuthal_average_phase(c::AbstractSSContributor, μ_a, μ_b,
                                  n_phi::Int)
    FT = promote_type(typeof(μ_a), typeof(μ_b), _contributor_numeric_type(c))
    μ_a = convert(FT, μ_a)
    μ_b = convert(FT, μ_b)
    n_phi > 0 || throw(ArgumentError("azimuth_nquad must be positive"))
    a = μ_a * μ_b
    b = sqrt(max(zero(FT), one(FT) - μ_a^2)) *
        sqrt(max(zero(FT), one(FT) - μ_b^2))
    total = zero(FT)
    for k in 0:(n_phi - 1)
        ϕ = FT(2) * FT(pi) * FT(k) / FT(n_phi)
        cosΘ = clamp(a + b * cos(ϕ), -one(FT), one(FT))
        total += exact_phase_function(c, cosΘ)
    end
    return total / FT(n_phi)
end

function _precompute_optics(config::ExactSSConfig)
    FT = _config_numeric_type(config)
    contributors = config.contributors
    n_layers, n_spec = _dims(contributors)
    n_geom = length(config.geometry.μv)

    τ_total_layer = zeros(FT, n_layers, n_spec)
    τ_scat_layer = zeros(FT, n_layers, n_spec)
    weighted_phase = zeros(FT, n_geom, n_layers, n_spec)

    for c in contributors
        τ = τ_matrix(c)
        ϖ = convert(FT, single_scattering_albedo(c))
        @inbounds for ispec in 1:n_spec, iz in 1:n_layers
            τᵢ = convert(FT, τ[iz, ispec])
            τ_total_layer[iz, ispec] += τᵢ
            scat_weight = τᵢ * ϖ
            τ_scat_layer[iz, ispec] += scat_weight
            if scat_weight != zero(FT)
                for iv in 1:n_geom
                    cosΘ = _scattering_angle_cosine(config.geometry.μ₀,
                                                    config.geometry.μv[iv],
                                                    config.geometry.Δϕ[iv])
                    weighted_phase[iv, iz, ispec] +=
                        scat_weight * exact_phase_function(c, cosΘ)
                end
            end
        end
    end

    ϖ_eff = zeros(FT, n_layers, n_spec)
    P_eff = zeros(FT, n_geom, n_layers, n_spec)
    τ_cum = zeros(FT, n_layers + 1, n_spec)

    @inbounds for ispec in 1:n_spec
        for iz in 1:n_layers
            τ = τ_total_layer[iz, ispec]
            τ_scat = τ_scat_layer[iz, ispec]
            ϖ_eff[iz, ispec] = τ == zero(FT) ? zero(FT) : τ_scat / τ
            τ_cum[iz + 1, ispec] = τ_cum[iz, ispec] + τ
            if τ_scat != zero(FT)
                for iv in 1:n_geom
                    P_eff[iv, iz, ispec] = weighted_phase[iv, iz, ispec] / τ_scat
                end
            end
        end
    end

    τ_total_column = vec(τ_cum[end, :])
    return (; τ_total_layer, τ_scat_layer, τ_cum, τ_total_column, ϖ_eff, P_eff)
end

function _precompute_azimuthal_phase(config::ExactSSConfig, τ_scat_layer, μ_nodes,
                                     reference_μ::AbstractVector)
    FT = _config_numeric_type(config)
    contributors = config.contributors
    n_layers, n_spec = _dims(contributors)
    n_geom = length(reference_μ)
    n_quad = length(μ_nodes)
    P̄ = zeros(FT, n_geom, n_layers, n_spec, n_quad)

    for c in contributors
        τ = τ_matrix(c)
        ϖ = convert(FT, single_scattering_albedo(c))
        ϖ == zero(FT) && continue
        @inbounds for ispec in 1:n_spec, iz in 1:n_layers
            scat_weight = convert(FT, τ[iz, ispec]) * ϖ
            scat_weight == zero(FT) && continue
            for iv in 1:n_geom, k in 1:n_quad
                P̄[iv, iz, ispec, k] += scat_weight *
                    _azimuthal_average_phase(c, convert(FT, μ_nodes[k]),
                                             convert(FT, reference_μ[iv]),
                                             config.azimuth_nquad)
            end
        end
    end

    return _normalize_azimuthal_phase!(P̄, τ_scat_layer, FT)
end

function _normalize_azimuthal_phase!(P̄, τ_scat_layer, ::Type{FT}) where {FT}
    n_geom, n_layers, n_spec, n_quad = size(P̄)
    @inbounds for ispec in 1:n_spec, iz in 1:n_layers
        τ_scat = τ_scat_layer[iz, ispec]
        τ_scat == zero(FT) && continue
        for iv in 1:n_geom, k in 1:n_quad
            P̄[iv, iz, ispec, k] /= τ_scat
        end
    end
    return P̄
end

function _precompute_azimuthal_phase_pair(config::ExactSSConfig, τ_scat_layer,
                                          μ_nodes, reference_a::AbstractVector,
                                          reference_b::AbstractVector)
    length(reference_a) == length(reference_b) ||
        throw(ArgumentError("paired azimuthal phase references must have the same length"))

    FT = _config_numeric_type(config)
    contributors = config.contributors
    n_layers, n_spec = _dims(contributors)
    n_geom = length(reference_a)
    n_quad = length(μ_nodes)
    P̄a = zeros(FT, n_geom, n_layers, n_spec, n_quad)
    P̄b = zeros(FT, n_geom, n_layers, n_spec, n_quad)

    for c in contributors
        τ = τ_matrix(c)
        ϖ = convert(FT, single_scattering_albedo(c))
        ϖ == zero(FT) && continue
        @inbounds for ispec in 1:n_spec, iz in 1:n_layers
            scat_weight = convert(FT, τ[iz, ispec]) * ϖ
            scat_weight == zero(FT) && continue
            for iv in 1:n_geom, k in 1:n_quad
                μ_node = convert(FT, μ_nodes[k])
                P̄a[iv, iz, ispec, k] += scat_weight *
                    _azimuthal_average_phase(c, μ_node,
                                             convert(FT, reference_a[iv]),
                                             config.azimuth_nquad)
                P̄b[iv, iz, ispec, k] += scat_weight *
                    _azimuthal_average_phase(c, μ_node,
                                             convert(FT, reference_b[iv]),
                                             config.azimuth_nquad)
            end
        end
    end

    _normalize_azimuthal_phase!(P̄a, τ_scat_layer, FT)
    _normalize_azimuthal_phase!(P̄b, τ_scat_layer, FT)
    return P̄a, P̄b
end

function _validate_config(config::ExactSSConfig, paths::Symbol)
    _validate_paths(paths)
    config.n_stokes == 1 ||
        throw(ArgumentError("StandaloneSS Phase 1 supports only n_stokes == 1"))
    if (_wants_path3(paths) || _wants_path4(paths)) &&
       !(config.surface isa LambertianSSSurface)
        throw(ArgumentError("StandaloneSS paths 3 and 4 currently support LambertianSSSurface only"))
    end
    _validate_geometry(config.geometry)
    config.inner_nquad > 0 ||
        throw(ArgumentError("ExactSSConfig.inner_nquad must be positive"))
    config.azimuth_nquad > 0 ||
        throw(ArgumentError("ExactSSConfig.azimuth_nquad must be positive"))
    n_layers, n_spec = _dims(config.contributors)
    n_spec > 0 && n_layers > 0 ||
        throw(ArgumentError("ExactSSConfig requires at least one layer and one spectral point"))
    return n_layers, n_spec
end

"""
    run_exact_ss(config::ExactSSConfig; paths=:paths_1_2)

Run standalone exact single-scattering for selected scalar paths. Path 2
supports Lambertian and Cox-Munk surfaces; paths 3 and 4 currently support
Lambertian surfaces. Radiance arrays have shape `(nGeom, 1, nSpec)`.
"""
function run_exact_ss(config::ExactSSConfig; paths::Symbol = :paths_1_2)
    FT = _config_numeric_type(config)
    _validate_config(config, paths)
    optics = _precompute_optics(config)

    n_geom = length(config.geometry.μv)
    n_spec = size(optics.τ_cum, 2)
    path1 = zeros(FT, n_geom, 1, n_spec)
    path2 = zeros(FT, n_geom, 1, n_spec)
    path3 = zeros(FT, n_geom, 1, n_spec)
    path4 = zeros(FT, n_geom, 1, n_spec)

    I0 = _vectorize_I0(config.I0, n_spec, FT)
    if _wants_path1(paths)
        _run_path1_kernel!(path1, optics.τ_cum, optics.ϖ_eff, optics.P_eff,
                           config.geometry, I0)
    end
    if _wants_path2(paths)
        surface_brdf = _precompute_surface_brdf(config.surface, config.geometry,
                                                n_spec, FT)
        _run_path2_kernel!(path2, optics.τ_total_column, config.geometry,
                           surface_brdf, I0)
    end
    if _wants_path3(paths) || _wants_path4(paths)
        albedo = _vectorize_albedo(config.surface, n_spec, FT)
        μ_nodes, μ_weights = _gauss_legendre_01(config.inner_nquad, FT)
        if _wants_path3(paths) && _wants_path4(paths)
            reference_μ₀ = fill(convert(FT, config.geometry.μ₀), n_geom)
            P̄3, P̄4 = _precompute_azimuthal_phase_pair(
                config, optics.τ_scat_layer, μ_nodes, reference_μ₀,
                config.geometry.μv)
            _run_path3_kernel!(path3, optics.τ_cum, optics.ϖ_eff, P̄3,
                               config.geometry, albedo, I0, μ_nodes, μ_weights)
            _run_path4_kernel!(path4, optics.τ_cum, optics.ϖ_eff, P̄4,
                               config.geometry, albedo, I0, μ_nodes, μ_weights)
        elseif _wants_path3(paths)
            reference_μ₀ = fill(convert(FT, config.geometry.μ₀), n_geom)
            P̄3 = _precompute_azimuthal_phase(config, optics.τ_scat_layer,
                                              μ_nodes, reference_μ₀)
            _run_path3_kernel!(path3, optics.τ_cum, optics.ϖ_eff, P̄3,
                               config.geometry, albedo, I0, μ_nodes, μ_weights)
        else
            P̄4 = _precompute_azimuthal_phase(config, optics.τ_scat_layer,
                                              μ_nodes, config.geometry.μv)
            _run_path4_kernel!(path4, optics.τ_cum, optics.ϖ_eff, P̄4,
                               config.geometry, albedo, I0, μ_nodes, μ_weights)
        end
    end

    total = path1 .+ path2 .+ path3 .+ path4
    quadrature_info = (; paths, kernel_backend = :KernelAbstractionsCPU,
                       inner_quadrature = config.inner_nquad,
                       azimuth_quadrature = config.azimuth_nquad)
    metadata = (; n_layers = size(optics.ϖ_eff, 1), n_spec, n_geom,
                n_stokes = config.n_stokes)
    return (; total, path1, path2, path3, path4, quadrature_info, metadata)
end
