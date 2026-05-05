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
