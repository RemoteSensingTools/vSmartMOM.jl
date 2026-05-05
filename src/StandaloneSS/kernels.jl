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
                               @Const(albedo), @Const(I0))
    iv, ispec = @index(Global, NTuple)
    FT = eltype(path2)
    μᵥ = μv[iv]
    τ = τ_total[ispec]
    path2[iv, 1, ispec] = (μ₀ * I0[ispec] * albedo[ispec] / FT(pi)) *
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

@kernel function _path3_kernel!(path3, @Const(τ_cum), @Const(ϖ_eff),
                               @Const(P̄), μ₀, @Const(μv), @Const(albedo),
                               @Const(I0), @Const(μ_nodes), @Const(μ_weights))
    iv, ispec = @index(Global, NTuple)
    FT = eltype(path3)
    μᵥ = μv[iv]
    τ_total = τ_cum[size(τ_cum, 1), ispec]

    F_surface = zero(FT)
    n_layers = size(ϖ_eff, 1)
    n_quad = length(μ_nodes)
    for iz in 1:n_layers
        ϖ = ϖ_eff[iz, ispec]
        if ϖ != zero(FT)
            τ_top = τ_cum[iz, ispec]
            τ_bot = τ_cum[iz + 1, ispec]
            layer_sum = zero(FT)
            for k in 1:n_quad
                μ_d = μ_nodes[k]
                τ_int = _τ_integral(τ_top, τ_bot, τ_total, μ₀, μ_d)
                layer_sum += μ_weights[k] * P̄[iv, iz, ispec, k] *
                             τ_int * FT(0.5)
            end
            F_surface += ϖ * I0[ispec] * layer_sum
        end
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
    n_quad = length(μ_nodes)
    for iz in 1:n_layers
        ϖ = ϖ_eff[iz, ispec]
        if ϖ != zero(FT)
            τ_top = τ_cum[iz, ispec]
            τ_bot = τ_cum[iz + 1, ispec]
            layer_sum = zero(FT)
            for k in 1:n_quad
                μ_u = μ_nodes[k]
                τ_int = _τ_integral(τ_top, τ_bot, τ_total, μᵥ, μ_u)
                layer_sum += μ_weights[k] * P̄[iv, iz, ispec, k] *
                             τ_int * FT(0.5)
            end
            value += L_surface * ϖ * layer_sum / μᵥ
        end
    end

    path4[iv, 1, ispec] = value
end

function _run_path1_kernel!(path1, τ_cum, ϖ_eff, P_eff, geometry, I0)
    backend = KernelAbstractions.CPU()
    kernel! = _path1_kernel!(backend)
    event = kernel!(path1, τ_cum, ϖ_eff, P_eff, geometry.μ₀, geometry.μv, I0;
                    ndrange=(size(path1, 1), size(path1, 3)))
    event === nothing || wait(event)
    return path1
end

function _run_path2_kernel!(path2, τ_total, geometry, albedo, I0)
    backend = KernelAbstractions.CPU()
    kernel! = _path2_kernel!(backend)
    event = kernel!(path2, τ_total, geometry.μ₀, geometry.μv, albedo, I0;
                    ndrange=(size(path2, 1), size(path2, 3)))
    event === nothing || wait(event)
    return path2
end

function _run_path3_kernel!(path3, τ_cum, ϖ_eff, P̄, geometry, albedo, I0,
                            μ_nodes, μ_weights)
    backend = KernelAbstractions.CPU()
    kernel! = _path3_kernel!(backend)
    event = kernel!(path3, τ_cum, ϖ_eff, P̄, geometry.μ₀, geometry.μv,
                    albedo, I0, μ_nodes, μ_weights;
                    ndrange=(size(path3, 1), size(path3, 3)))
    event === nothing || wait(event)
    return path3
end

function _run_path4_kernel!(path4, τ_cum, ϖ_eff, P̄, geometry, albedo, I0,
                            μ_nodes, μ_weights)
    backend = KernelAbstractions.CPU()
    kernel! = _path4_kernel!(backend)
    event = kernel!(path4, τ_cum, ϖ_eff, P̄, geometry.μ₀, geometry.μv,
                    albedo, I0, μ_nodes, μ_weights;
                    ndrange=(size(path4, 1), size(path4, 3)))
    event === nothing || wait(event)
    return path4
end
