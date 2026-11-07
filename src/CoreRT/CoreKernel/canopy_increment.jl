"""
    canopy_increment_doubling!(pol_type, SFI, dτ, properties, m, ndoubl,
                               quad_points, added_layer, I_static, architecture)

Double a directional canopy layer using `Δ = T - I` and the direct-beam
increment `β = exp(-δτ G₀/μ₀) - 1`. The diagonal increments are constructed
from optical properties with `expm1`, before near-unit transmission has been
rounded. This permits thin elemental layers in Float32 without accumulating
their lost diagonal increments into a spurious absorption or gain.

Call after directional `elemental!`, which supplies the off-diagonal
transmission, reflection and solar sources in the usual D-transformed basis.
On return all four R/T operators and the solar sources use the ordinary
production representation. No flux, phase matrix or absorption is rescaled.

This canopy-specific path currently carries the solar source only. A future
distributed canopy source must supply its own attenuation increment explicitly.
"""
function canopy_increment_doubling!(pol_type, SFI, dτ,
        properties::CoreDirectionalScatteringOpticalProperties, m::Int, ndoubl::Int,
        quad_points, added_layer, I_static, architecture)
    ndoubl == 0 && return nothing
    isempty(added_layer.j₀_by_src) || throw(ArgumentError(
        "Canopy increment doubling requires explicit increments for non-solar sources"))
    (; r⁻⁺, r⁺⁻, t⁺⁺, t⁻⁻, j₀⁺, j₀⁻, temp1, temp2,
       dbl_gp_refl, dbl_j₁⁺, dbl_j₁⁻, dbl_v1, dbl_v2) = added_layer
    FT = eltype(t⁺⁺)
    (; ϖ, G, Z⁺⁺) = properties
    (; qp_μN, wt_μN, iμ₀Nstart, μ₀) = quad_points
    arr_type = array_type(architecture)
    weights = arr_type(scaled_weights(m, wt_μN))
    kernel! = canopy_transmission_increment!(devi(architecture))
    kernel!(t⁺⁺, ϖ, dτ, G, Z⁺⁺, arr_type(qp_μN), weights;
        ndrange=(size(t⁺⁺, 1), size(t⁺⁺, 3)))
    synchronize_if_gpu()

    # The reverse-direction operators are reconstructed after doubling, so
    # their buffers can hold H=(I-RR)^(-1)RR and C=T(I-RR)^(-1)-I meanwhile.
    h, c = r⁺⁻, t⁻⁻
    rr = temp1 === nothing ? similar(t⁺⁺) : temp1
    system = temp2 === nothing ? similar(t⁺⁺) : temp2
    product = dbl_gp_refl === nothing ? similar(t⁺⁺) : dbl_gp_refl
    jd1 = dbl_j₁⁺ === nothing ? similar(j₀⁺) : dbl_j₁⁺
    ju1 = dbl_j₁⁻ === nothing ? similar(j₀⁻) : dbl_j₁⁻
    v1 = dbl_v1 === nothing ? similar(j₀⁺) : dbl_v1
    v2 = dbl_v2 === nothing ? similar(j₀⁺) : dbl_v2
    beta = expm1.(-dτ .* view(G, iμ₀Nstart:iμ₀Nstart) ./ μ₀)
    beta_row = reshape(beta, 1, 1, :)

    for _ in 1:ndoubl
        _bmm!(rr, r⁻⁺, r⁻⁺)
        system .= I_static .- rr
        batch_solve!(h, system, rr)
        _bmm!(product, t⁺⁺, h)
        c .= t⁺⁺ .+ h .+ product

        if SFI
            jd1 .= j₀⁺ .+ beta_row .* j₀⁺
            ju1 .= j₀⁻ .+ beta_row .* j₀⁻
            _bmm!(v1, r⁻⁺, j₀⁺)
            v1 .+= ju1
            _bmm!(v2, c, v1)
            j₀⁻ .+= v1 .+ v2
            _bmm!(v1, r⁻⁺, ju1)
            v1 .+= j₀⁺
            _bmm!(v2, c, v1)
            j₀⁺ .= jd1 .+ v1 .+ v2
        end

        # Q=R+RΔ; R'=R+Q+CQ; Δ'=C+Δ+CΔ. Form no near-identity T.
        _bmm!(rr, r⁻⁺, t⁺⁺)
        rr .+= r⁻⁺
        _bmm!(product, c, rr)
        r⁻⁺ .+= rr .+ product
        _bmm!(product, c, t⁺⁺)
        t⁺⁺ .= c .+ t⁺⁺ .+ product
        beta .= beta .* (FT(2) .+ beta)
    end
    t⁺⁺ .+= I_static
    apply_D_matrix!(pol_type.n, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻)
    SFI && apply_D_matrix_SFI!(pol_type.n, j₀⁻)
    synchronize_if_gpu()
    return nothing
end

"Construct the diagonal of T-I directly; off-diagonal entries already equal T."
@kernel function canopy_transmission_increment!(delta, @Const(omega), @Const(dt),
        @Const(G), @Const(Zpp), @Const(mu), @Const(weights))
    i, n = @index(Global, NTuple)
    FT = eltype(delta)
    spectral = size(Zpp, 3) > 1 ? n : 1
    x = dt[n] * G[i] / mu[i]
    scattering = weights[i] > rt_weight_tol(eltype(weights)) ?
        omega[n] * G[i] * Zpp[i,i,spectral] * (dt[n]/mu[i]) * weights[i] : zero(FT)
    delta[i,i,n] = expm1(-x) + exp(-x)*scattering
end
