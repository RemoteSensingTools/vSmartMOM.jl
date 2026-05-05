@kernel function _path1_kernel!(path1, @Const(τ_cum), @Const(ϖ_eff),
                               @Const(P_eff), μ₀, @Const(μv), @Const(I0))
    iv, ispec = @index(Global, NTuple)
    FT = eltype(path1)
    μᵥ = μv[iv]
    a = inv(μ₀) + inv(μᵥ)
    prefactor = I0[ispec] / (FT(4) * FT(pi) * μᵥ * a)

    value = zero(FT)
    n_layers = size(ϖ_eff, 1)
    for iz in 1:n_layers
        layer_factor = exp(-τ_cum[iz, ispec] * a) -
                       exp(-τ_cum[iz + 1, ispec] * a)
        value += prefactor * ϖ_eff[iz, ispec] *
                 P_eff[iv, iz, ispec] * layer_factor
    end
    path1[iv, 1, ispec] = value
end

@kernel function _path2_kernel!(path2, @Const(τ_total), μ₀, @Const(μv),
                               @Const(surface_brdf), @Const(I0))
    iv, ispec = @index(Global, NTuple)
    FT = eltype(path2)
    μᵥ = μv[iv]
    τ = τ_total[ispec]
    path2[iv, 1, ispec] = (μ₀ * I0[ispec] * surface_brdf[iv, ispec]) *
                          exp(-τ / μ₀) * exp(-τ / μᵥ)
end

@inline function _τ_integral(τ_top, τ_bot, τ_total, μ_first, μ_second)
    FT = typeof(τ_top)
    b = inv(μ_first) - inv(μ_second)
    f_top = exp(-τ_top / μ_first - (τ_total - τ_top) / μ_second)
    f_bot = exp(-τ_bot / μ_first - (τ_total - τ_bot) / μ_second)
    if abs(b) < FT(1e-10)
        return FT(0.5) * (f_top + f_bot) * (τ_bot - τ_top)
    end
    return (f_top - f_bot) / b
end

@inline function _inner_scatter_sum(τ_cum, ϖ_eff, P̄, τ_total, iv, ispec,
                                    iz, μ_first, μ_nodes, μ_weights)
    FT = typeof(τ_total)
    ϖ = ϖ_eff[iz, ispec]
    ϖ == zero(FT) && return zero(FT)

    τ_top = τ_cum[iz, ispec]
    τ_bot = τ_cum[iz + 1, ispec]
    layer_sum = zero(FT)
    n_quad = length(μ_nodes)
    for k in 1:n_quad
        μ_inner = μ_nodes[k]
        τ_int = _τ_integral(τ_top, τ_bot, τ_total, μ_first, μ_inner)
        layer_sum += μ_weights[k] * P̄[iv, iz, ispec, k] *
                     τ_int * FT(0.5)
    end
    return ϖ * layer_sum
end

@inline function _rayleigh_azimuthal_average(μ_a::FT, μ_b::FT) where {FT}
    a = μ_a * μ_b
    b = sqrt(max(zero(FT), one(FT) - μ_a^2)) *
        sqrt(max(zero(FT), one(FT) - μ_b^2))
    return FT(0.75) * (one(FT) + a^2 + FT(0.5) * b^2)
end

@inline function _phase_azimuth_average_kind(kind::Int32, g::FT, μ_a::FT,
                                             μ_b::FT, n_phi::Int) where {FT}
    if kind == Int32(1)
        return _rayleigh_azimuthal_average(μ_a, μ_b)
    elseif kind == Int32(2)
        a = μ_a * μ_b
        b = sqrt(max(zero(FT), one(FT) - μ_a^2)) *
            sqrt(max(zero(FT), one(FT) - μ_b^2))
        total = zero(FT)
        for j in 0:(n_phi - 1)
            ϕ = FT(2) * FT(pi) * FT(j) / FT(n_phi)
            cosΘ = clamp(a + b * cos(ϕ), -one(FT), one(FT))
            total += (one(g) - g^2) /
                     (one(g) + g^2 - FT(2) * g * cosΘ)^(FT(1.5))
        end
        return total / FT(n_phi)
    end
    return zero(FT)
end

@kernel function _azimuthal_phase_pair_kernel!(
    P̄a, P̄b, @Const(τ_scat_layer), @Const(τ_contrib),
    @Const(ϖ_contrib), @Const(g_contrib), @Const(kind_contrib),
    @Const(μ_nodes), @Const(reference_a), @Const(reference_b), n_phi)

    iv, iz, ispec, k = @index(Global, NTuple)
    FT = eltype(P̄a)
    τ_scat = τ_scat_layer[iz, ispec]
    if τ_scat == zero(FT)
        P̄a[iv, iz, ispec, k] = zero(FT)
        P̄b[iv, iz, ispec, k] = zero(FT)
    else
        μ_node = μ_nodes[k]
        weighted_a = zero(FT)
        weighted_b = zero(FT)
        n_contrib = size(τ_contrib, 1)
        for ic in 1:n_contrib
            scat_weight = τ_contrib[ic, iz, ispec] * ϖ_contrib[ic]
            if scat_weight != zero(FT)
                kind = kind_contrib[ic]
                g = g_contrib[ic]
                weighted_a += scat_weight *
                              _phase_azimuth_average_kind(kind, g, μ_node,
                                                          reference_a[iv], n_phi)
                weighted_b += scat_weight *
                              _phase_azimuth_average_kind(kind, g, μ_node,
                                                          reference_b[iv], n_phi)
            end
        end

        P̄a[iv, iz, ispec, k] = weighted_a / τ_scat
        P̄b[iv, iz, ispec, k] = weighted_b / τ_scat
    end
end

@kernel function _path3_kernel!(path3, @Const(τ_cum), @Const(ϖ_eff),
                               @Const(P̄), μ₀, @Const(μv), @Const(albedo),
                               @Const(I0), @Const(μ_nodes), @Const(μ_weights))
    iv, ispec = @index(Global, NTuple)
    FT = eltype(path3)
    μᵥ = μv[iv]
    τ_total = τ_cum[size(τ_cum, 1), ispec]

    F_surface = zero(FT)
    n_layers = size(ϖ_eff, 1)
    for iz in 1:n_layers
        F_surface += I0[ispec] *
                     _inner_scatter_sum(τ_cum, ϖ_eff, P̄, τ_total, iv,
                                        ispec, iz, μ₀, μ_nodes, μ_weights)
    end

    path3[iv, 1, ispec] =
        (albedo[ispec] / FT(pi)) * F_surface * exp(-τ_total / μᵥ)
end

@kernel function _path4_kernel!(path4, @Const(τ_cum), @Const(ϖ_eff),
                               @Const(P̄), μ₀, @Const(μv), @Const(albedo),
                               @Const(I0), @Const(μ_nodes), @Const(μ_weights))
    iv, ispec = @index(Global, NTuple)
    FT = eltype(path4)
    μᵥ = μv[iv]
    τ_total = τ_cum[size(τ_cum, 1), ispec]
    L_surface = (albedo[ispec] / FT(pi)) * μ₀ * I0[ispec] *
                exp(-τ_total / μ₀)

    value = zero(FT)
    n_layers = size(ϖ_eff, 1)
    for iz in 1:n_layers
        value += L_surface *
                 _inner_scatter_sum(τ_cum, ϖ_eff, P̄, τ_total, iv, ispec,
                                    iz, μᵥ, μ_nodes, μ_weights) / μᵥ
    end

    path4[iv, 1, ispec] = value
end

@kernel function _path34_kernel!(path3, path4, @Const(τ_cum), @Const(ϖ_eff),
                                @Const(P̄3), @Const(P̄4), μ₀, @Const(μv),
                                @Const(albedo), @Const(I0), @Const(μ_nodes),
                                @Const(μ_weights))
    iv, ispec = @index(Global, NTuple)
    FT = eltype(path3)
    μᵥ = μv[iv]
    τ_total = τ_cum[size(τ_cum, 1), ispec]
    ρ = albedo[ispec] / FT(pi)

    if ρ == zero(FT)
        path3[iv, 1, ispec] = zero(FT)
        path4[iv, 1, ispec] = zero(FT)
    else
        F_surface = zero(FT)
        path4_sum = zero(FT)
        n_layers = size(ϖ_eff, 1)
        for iz in 1:n_layers
            F_surface += I0[ispec] *
                         _inner_scatter_sum(τ_cum, ϖ_eff, P̄3, τ_total,
                                            iv, ispec, iz, μ₀, μ_nodes,
                                            μ_weights)
            path4_sum += _inner_scatter_sum(τ_cum, ϖ_eff, P̄4, τ_total,
                                            iv, ispec, iz, μᵥ, μ_nodes,
                                            μ_weights)
        end

        path3[iv, 1, ispec] =
            ρ * F_surface * exp(-τ_total / μᵥ)
        path4[iv, 1, ispec] =
            ρ * μ₀ * I0[ispec] * exp(-τ_total / μ₀) * path4_sum / μᵥ
    end
end

@inline function _synchronize_kernel(event, backend)
    event === nothing || wait(event)
    KernelAbstractions.synchronize(backend)
    return nothing
end

function _run_path1_kernel!(path1, τ_cum, ϖ_eff, P_eff, μ₀, μv, I0, backend)
    kernel! = _path1_kernel!(backend)
    event = kernel!(path1, τ_cum, ϖ_eff, P_eff, μ₀, μv, I0;
                    ndrange=(size(path1, 1), size(path1, 3)))
    _synchronize_kernel(event, backend)
    return path1
end

function _run_path2_kernel!(path2, τ_total, μ₀, μv, surface_brdf, I0, backend)
    kernel! = _path2_kernel!(backend)
    event = kernel!(path2, τ_total, μ₀, μv, surface_brdf, I0;
                    ndrange=(size(path2, 1), size(path2, 3)))
    _synchronize_kernel(event, backend)
    return path2
end

function _run_path3_kernel!(path3, τ_cum, ϖ_eff, P̄, μ₀, μv, albedo, I0,
                            μ_nodes, μ_weights, backend)
    kernel! = _path3_kernel!(backend)
    event = kernel!(path3, τ_cum, ϖ_eff, P̄, μ₀, μv, albedo, I0,
                    μ_nodes, μ_weights;
                    ndrange=(size(path3, 1), size(path3, 3)))
    _synchronize_kernel(event, backend)
    return path3
end

function _run_path4_kernel!(path4, τ_cum, ϖ_eff, P̄, μ₀, μv, albedo, I0,
                            μ_nodes, μ_weights, backend)
    kernel! = _path4_kernel!(backend)
    event = kernel!(path4, τ_cum, ϖ_eff, P̄, μ₀, μv, albedo, I0,
                    μ_nodes, μ_weights;
                    ndrange=(size(path4, 1), size(path4, 3)))
    _synchronize_kernel(event, backend)
    return path4
end

function _run_path34_kernel!(path3, path4, τ_cum, ϖ_eff, P̄3, P̄4, μ₀, μv,
                             albedo, I0, μ_nodes, μ_weights, backend)
    kernel! = _path34_kernel!(backend)
    event = kernel!(path3, path4, τ_cum, ϖ_eff, P̄3, P̄4, μ₀, μv,
                    albedo, I0, μ_nodes, μ_weights;
                    ndrange=(size(path3, 1), size(path3, 3)))
    _synchronize_kernel(event, backend)
    return path3, path4
end

function _run_azimuthal_phase_pair_kernel!(P̄a, P̄b, τ_scat_layer, τ_contrib,
                                           ϖ_contrib, g_contrib, kind_contrib,
                                           μ_nodes, reference_a, reference_b,
                                           n_phi::Int, backend)
    kernel! = _azimuthal_phase_pair_kernel!(backend)
    event = kernel!(P̄a, P̄b, τ_scat_layer, τ_contrib, ϖ_contrib,
                    g_contrib, kind_contrib, μ_nodes, reference_a, reference_b,
                    n_phi; ndrange=size(P̄a))
    _synchronize_kernel(event, backend)
    return P̄a, P̄b
end
