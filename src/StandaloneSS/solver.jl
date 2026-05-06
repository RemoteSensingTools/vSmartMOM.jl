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

_vector_stokes_paths_supported(paths::Symbol) =
    paths in (:path1, :path2, :paths_1_2)
_optics_dispatch(paths::Symbol, stokes_dispatch) =
    _wants_path1(paths) ? stokes_dispatch : Val(1)

_architecture_backend(architecture::AbstractArchitecture) = devi(architecture)
_architecture_array_type(architecture::AbstractArchitecture) = array_type(architecture)

_to_arch(::CPU, x::AbstractArray) = x
_to_arch(architecture::AbstractArchitecture, x::AbstractArray) =
    _architecture_array_type(architecture)(x)

_zeros_arch(::CPU, ::Type{FT}, dims::Int...) where {FT} = zeros(FT, dims...)
_zeros_arch(architecture::AbstractArchitecture, ::Type{FT}, dims::Int...) where {FT} =
    _architecture_array_type(architecture)(zeros(FT, dims...))

_path_zeros(architecture::AbstractArchitecture, ::Type{FT}, n_geom::Int,
            n_spec::Int) where {FT} =
    _path_zeros(architecture, FT, n_geom, 1, n_spec)

_path_zeros(architecture::AbstractArchitecture, ::Type{FT}, n_geom::Int,
            n_stokes::Int, n_spec::Int) where {FT} =
    _zeros_arch(architecture, FT, n_geom, n_stokes, n_spec)

_path_or_zeros(path, architecture::AbstractArchitecture, ::Type{FT},
               n_geom::Int, n_spec::Int) where {FT} =
    _path_or_zeros(path, architecture, FT, n_geom, 1, n_spec)

_path_or_zeros(path, architecture::AbstractArchitecture, ::Type{FT},
               n_geom::Int, n_stokes::Int, n_spec::Int) where {FT} =
    path === nothing ? _path_zeros(architecture, FT, n_geom, n_stokes, n_spec) : path

function _requested_total(paths::Symbol, path1, path2, path3, path4)
    if _wants_path1(paths)
        total = copy(path1)
        _wants_path2(paths) && (total .+= path2)
        _wants_path3(paths) && (total .+= path3)
        _wants_path4(paths) && (total .+= path4)
        return total
    elseif _wants_path2(paths)
        total = copy(path2)
        _wants_path3(paths) && (total .+= path3)
        _wants_path4(paths) && (total .+= path4)
        return total
    elseif _wants_path3(paths)
        total = copy(path3)
        _wants_path4(paths) && (total .+= path4)
        return total
    else
        return copy(path4)
    end
end

_kernel_backend_name(::CPU) = :KernelAbstractionsCPU
_kernel_backend_name(::GPU) = :CUDA
_kernel_backend_name(::MetalGPU) = :Metal
_kernel_backend_name(architecture::AbstractArchitecture) =
    Symbol(nameof(typeof(architecture)))

_validate_architecture(::CPU) = nothing
function _validate_architecture(architecture::AbstractArchitecture)
    try
        _architecture_backend(architecture)
        _architecture_array_type(architecture)
    catch err
        throw(ArgumentError("StandaloneSS architecture $(typeof(architecture)) " *
                            "is not available. Load a functional backend " *
                            "package, e.g. `using CUDA`, before using GPU()."))
    end
    return nothing
end

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

function _azimuthal_average_phase(c::RayleighSSContributor, μ_a, μ_b,
                                  n_phi::Int)
    n_phi > 0 || throw(ArgumentError("azimuth_nquad must be positive"))
    FT = promote_type(typeof(μ_a), typeof(μ_b), typeof(c.depol))
    return _rayleigh_azimuthal_average(
        convert(FT, μ_a), convert(FT, μ_b), convert(FT, c.depol))
end

_precompute_optics(config::ExactSSConfig) =
    _precompute_optics(Val(1), config)

function _unsupported_vector_path1_contributor(c)
    throw(ArgumentError("StandaloneSS vector path1 currently supports RayleighSSContributor, GreekCoefsSSContributor, and AbsorptionSSContributor only; got $(typeof(c))"))
end

_validate_vector_path1_contributor(::RayleighSSContributor) = nothing
_validate_vector_path1_contributor(::GreekCoefsSSContributor) = nothing
_validate_vector_path1_contributor(::AbsorptionSSContributor) = nothing
_validate_vector_path1_contributor(c::AbstractSSContributor) =
    _unsupported_vector_path1_contributor(c)

function _validate_vector_path1_contributors(contributors)
    foreach(_validate_vector_path1_contributor, contributors)
    return nothing
end

function _phase_first_column(c::RayleighSSContributor, pol::Val{N}, μ₀::FT,
                             μᵥ::FT, Δϕ::FT) where {N,FT}
    greek = Scattering.get_greek_rayleigh(convert(FT, c.depol))
    return Scattering.phase_matrix_first_column(greek, μ₀, μᵥ, Δϕ, pol)
end

function _phase_first_column(c::GreekCoefsSSContributor, pol::Val{N}, μ₀::FT,
                             μᵥ::FT, Δϕ::FT) where {N,FT}
    return Scattering.phase_matrix_first_column(c.greek_coefs, μ₀, μᵥ, Δϕ,
                                                pol)
end

function _precompute_optics(::Val{1}, config::ExactSSConfig)
    FT = _config_numeric_type(config)
    contributors = config.contributors
    n_layers, n_spec = _dims(contributors)
    n_geom = length(config.geometry.μv)

    τ_total_layer = zeros(FT, n_layers, n_spec)
    τ_scat_layer = zeros(FT, n_layers, n_spec)
    weighted_phase = zeros(FT, n_geom, n_layers, n_spec)
    cosΘ = Vector{FT}(undef, n_geom)
    @inbounds for iv in 1:n_geom
        cosΘ[iv] = _scattering_angle_cosine(config.geometry.μ₀,
                                            config.geometry.μv[iv],
                                            config.geometry.Δϕ[iv])
    end

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
                    weighted_phase[iv, iz, ispec] +=
                        scat_weight * exact_phase_function(c, cosΘ[iv])
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

function _precompute_optics(pol::Val{N}, config::ExactSSConfig) where {N}
    N > 1 || return _precompute_optics(Val(1), config)
    _validate_vector_path1_contributors(config.contributors)

    FT = _config_numeric_type(config)
    contributors = config.contributors
    n_layers, n_spec = _dims(contributors)
    n_geom = length(config.geometry.μv)

    τ_total_layer = zeros(FT, n_layers, n_spec)
    τ_scat_layer = zeros(FT, n_layers, n_spec)
    weighted_phase = zeros(FT, n_geom, N, n_layers, n_spec)

    μ₀ = convert(FT, config.geometry.μ₀)
    phase_cols = Vector{NTuple{N,FT}}(undef, n_geom)

    for c in contributors
        τ = τ_matrix(c)
        ϖ = convert(FT, single_scattering_albedo(c))
        if c isa Union{RayleighSSContributor, GreekCoefsSSContributor}
            @inbounds for iv in 1:n_geom
                phase_cols[iv] = _phase_first_column(
                    c, pol, μ₀, convert(FT, config.geometry.μv[iv]),
                    convert(FT, config.geometry.Δϕ[iv]))
            end
        end

        @inbounds for ispec in 1:n_spec, iz in 1:n_layers
            τᵢ = convert(FT, τ[iz, ispec])
            τ_total_layer[iz, ispec] += τᵢ
            scat_weight = τᵢ * ϖ
            τ_scat_layer[iz, ispec] += scat_weight
            if scat_weight != zero(FT) &&
               c isa Union{RayleighSSContributor, GreekCoefsSSContributor}
                for iv in 1:n_geom, istokes in 1:N
                    weighted_phase[iv, istokes, iz, ispec] +=
                        scat_weight * phase_cols[iv][istokes]
                end
            end
        end
    end

    ϖ_eff = zeros(FT, n_layers, n_spec)
    P_eff = zeros(FT, n_geom, N, n_layers, n_spec)
    τ_cum = zeros(FT, n_layers + 1, n_spec)

    @inbounds for ispec in 1:n_spec
        for iz in 1:n_layers
            τ = τ_total_layer[iz, ispec]
            τ_scat = τ_scat_layer[iz, ispec]
            ϖ_eff[iz, ispec] = τ == zero(FT) ? zero(FT) : τ_scat / τ
            τ_cum[iz + 1, ispec] = τ_cum[iz, ispec] + τ
            if τ_scat != zero(FT)
                for iv in 1:n_geom, istokes in 1:N
                    P_eff[iv, istokes, iz, ispec] =
                        weighted_phase[iv, istokes, iz, ispec] / τ_scat
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

_contributor_kind(::RayleighSSContributor) = Int32(1)
_contributor_kind(::HGAerosolSSContributor) = Int32(2)
_contributor_kind(::AbsorptionSSContributor) = Int32(0)

_contributor_g(c::RayleighSSContributor, ::Type{FT}) where {FT} =
    convert(FT, c.depol)
_contributor_g(c::HGAerosolSSContributor, ::Type{FT}) where {FT} =
    convert(FT, c.g)
_contributor_g(::AbsorptionSSContributor, ::Type{FT}) where {FT} = zero(FT)

function _pack_contributors(contributors, n_layers::Int, n_spec::Int,
                            ::Type{FT}) where {FT}
    n_contrib = length(contributors)
    τ_contrib = zeros(FT, n_contrib, n_layers, n_spec)
    ϖ_contrib = zeros(FT, n_contrib)
    g_contrib = zeros(FT, n_contrib)
    kind_contrib = zeros(Int32, n_contrib)

    for (ic, c) in enumerate(contributors)
        τ = τ_matrix(c)
        @inbounds for ispec in 1:n_spec, iz in 1:n_layers
            τ_contrib[ic, iz, ispec] = convert(FT, τ[iz, ispec])
        end
        ϖ_contrib[ic] = convert(FT, single_scattering_albedo(c))
        g_contrib[ic] = _contributor_g(c, FT)
        kind_contrib[ic] = _contributor_kind(c)
    end

    return (; τ_contrib, ϖ_contrib, g_contrib, kind_contrib)
end

function _precompute_azimuthal_phase_pair_arch(config::ExactSSConfig,
                                               τ_scat_layer, μ_nodes,
                                               reference_a::AbstractVector,
                                               reference_b::AbstractVector,
                                               architecture::CPU, backend)
    return _precompute_azimuthal_phase_pair(config, τ_scat_layer, μ_nodes,
                                            reference_a, reference_b)
end

function _precompute_azimuthal_phase_pair_arch(config::ExactSSConfig,
                                               τ_scat_layer, μ_nodes,
                                               reference_a::AbstractVector,
                                               reference_b::AbstractVector,
                                               architecture::AbstractArchitecture,
                                               backend)
    length(reference_a) == length(reference_b) ||
        throw(ArgumentError("paired azimuthal phase references must have the same length"))

    FT = _config_numeric_type(config)
    n_layers, n_spec = _dims(config.contributors)
    n_geom = length(reference_a)
    n_quad = length(μ_nodes)
    P̄a = _zeros_arch(architecture, FT, n_geom, n_layers, n_spec, n_quad)
    P̄b = _zeros_arch(architecture, FT, n_geom, n_layers, n_spec, n_quad)
    packed = _pack_contributors(config.contributors, n_layers, n_spec, FT)

    _run_azimuthal_phase_pair_kernel!(
        P̄a, P̄b,
        _to_arch(architecture, τ_scat_layer),
        _to_arch(architecture, packed.τ_contrib),
        _to_arch(architecture, packed.ϖ_contrib),
        _to_arch(architecture, packed.g_contrib),
        _to_arch(architecture, packed.kind_contrib),
        _to_arch(architecture, μ_nodes),
        _to_arch(architecture, reference_a),
        _to_arch(architecture, reference_b),
        config.azimuth_nquad,
        backend)

    return P̄a, P̄b
end

function _validate_config(config::ExactSSConfig, paths::Symbol)
    _validate_paths(paths)
    n_stokes = _config_n_stokes(config)
    if n_stokes != 1 && !_vector_stokes_paths_supported(paths)
        throw(ArgumentError("StandaloneSS vector Stokes kernels currently support only `paths=:path1`, `paths=:path2`, and `paths=:paths_1_2`; received n_stokes=$n_stokes and paths=:$paths"))
    end
    n_stokes == 1 || !_wants_path1(paths) ||
        _validate_vector_path1_contributors(config.contributors)
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
    _validate_architecture(config.architecture)
    return n_layers, n_spec
end

"""
    run_exact_ss(config::ExactSSConfig; paths=:paths_1_2)

Run standalone exact single-scattering for selected paths. Stokes-vector
path 1 currently supports Rayleigh scattering; path 2 supports Lambertian and
Cox-Munk direct-beam reflection. Paths 3 and 4 currently support scalar
Lambertian surfaces. Radiance arrays have shape `(nGeom, nStokes, nSpec)`.
"""
function run_exact_ss(config::ExactSSConfig; paths::Symbol = :paths_1_2)
    FT = _config_numeric_type(config)
    _validate_config(config, paths)
    architecture = config.architecture
    backend = _architecture_backend(architecture)

    n_geom = length(config.geometry.μv)
    n_stokes = _config_n_stokes(config)
    stokes_dispatch = Val(n_stokes)
    optics = _precompute_optics(_optics_dispatch(paths, stokes_dispatch), config)
    n_spec = size(optics.τ_cum, 2)
    path1 = nothing
    path2 = nothing
    path3 = nothing
    path4 = nothing

    μv = _to_arch(architecture, config.geometry.μv)
    τ_cum = _to_arch(architecture, optics.τ_cum)
    τ_total_column = _to_arch(architecture, optics.τ_total_column)
    ϖ_eff = _to_arch(architecture, optics.ϖ_eff)
    P_eff = _to_arch(architecture, optics.P_eff)

    I0 = _to_arch(architecture, _vectorize_I0(config.I0, n_spec, FT))
    if _wants_path1(paths)
        path1 = _path_zeros(architecture, FT, n_geom, n_stokes, n_spec)
        _run_path1_kernel!(stokes_dispatch, path1, τ_cum, ϖ_eff, P_eff,
                           config.geometry.μ₀, μv, I0, backend)
    end
    if _wants_path2(paths)
        path2 = _path_zeros(architecture, FT, n_geom, n_stokes, n_spec)
        surface_brdf = _precompute_surface_brdf(
            stokes_dispatch, config.surface, config.geometry, n_spec, FT)
        surface_brdf = _to_arch(architecture, surface_brdf)
        _run_path2_kernel!(stokes_dispatch, path2, τ_total_column,
                           config.geometry.μ₀, μv, surface_brdf, I0, backend)
    end
    if _wants_path3(paths) || _wants_path4(paths)
        albedo = _to_arch(architecture, _vectorize_albedo(config.surface, n_spec, FT))
        μ_nodes_host, μ_weights_host = _gauss_legendre_01(config.inner_nquad, FT)
        μ_nodes = _to_arch(architecture, μ_nodes_host)
        μ_weights = _to_arch(architecture, μ_weights_host)
        if _wants_path3(paths) && _wants_path4(paths)
            path3 = _path_zeros(architecture, FT, n_geom, n_stokes, n_spec)
            path4 = _path_zeros(architecture, FT, n_geom, n_stokes, n_spec)
            reference_μ₀ = fill(convert(FT, config.geometry.μ₀), n_geom)
            P̄3, P̄4 = _precompute_azimuthal_phase_pair_arch(
                config, optics.τ_scat_layer, μ_nodes_host, reference_μ₀,
                config.geometry.μv, architecture, backend)
            _run_path34_kernel!(stokes_dispatch, path3, path4, τ_cum, ϖ_eff,
                                P̄3, P̄4, config.geometry.μ₀, μv, albedo,
                                I0, μ_nodes, μ_weights, backend)
        elseif _wants_path3(paths)
            path3 = _path_zeros(architecture, FT, n_geom, n_stokes, n_spec)
            reference_μ₀ = fill(convert(FT, config.geometry.μ₀), n_geom)
            P̄3 = _precompute_azimuthal_phase(config, optics.τ_scat_layer,
                                              μ_nodes_host, reference_μ₀)
            P̄3 = _to_arch(architecture, P̄3)
            _run_path3_kernel!(stokes_dispatch, path3, τ_cum, ϖ_eff, P̄3,
                               config.geometry.μ₀, μv, albedo, I0, μ_nodes,
                               μ_weights, backend)
        else
            path4 = _path_zeros(architecture, FT, n_geom, n_stokes, n_spec)
            P̄4 = _precompute_azimuthal_phase(config, optics.τ_scat_layer,
                                              μ_nodes_host, config.geometry.μv)
            P̄4 = _to_arch(architecture, P̄4)
            _run_path4_kernel!(stokes_dispatch, path4, τ_cum, ϖ_eff, P̄4,
                               config.geometry.μ₀, μv, albedo, I0, μ_nodes,
                               μ_weights, backend)
        end
    end

    total = _requested_total(paths, path1, path2, path3, path4)
    path1 = _path_or_zeros(path1, architecture, FT, n_geom, n_stokes, n_spec)
    path2 = _path_or_zeros(path2, architecture, FT, n_geom, n_stokes, n_spec)
    path3 = _path_or_zeros(path3, architecture, FT, n_geom, n_stokes, n_spec)
    path4 = _path_or_zeros(path4, architecture, FT, n_geom, n_stokes, n_spec)
    quadrature_info = (; paths, kernel_backend = _kernel_backend_name(architecture),
                       inner_quadrature = config.inner_nquad,
                       azimuth_quadrature = config.azimuth_nquad)
    metadata = (; n_layers = size(optics.ϖ_eff, 1), n_spec, n_geom,
                n_stokes)
    return (; total, path1, path2, path3, path4, quadrature_info, metadata)
end
