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
distributed canopy source needs its own sublayer source construction, rather
than inheriting the solar attenuation recurrence.
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
    backend = KernelAbstractions.get_backend(h)
    fused_solve = _use_fused_gp(h)
    # The existing 3N² solve gate is conservative for both update kernels.
    # Honor the ordinary doubling switch as well for independent comparisons.
    fused_updates = fused_solve && _FUSED_DOUBLING_ENABLED[]

    for _ in 1:ndoubl
        _bmm!(rr, r⁻⁺, r⁻⁺)
        if fused_solve
            # RR commutes with its resolvent: RR*(I-RR)^(-1) equals
            # (I-RR)^(-1)*RR. Keep RR as the small RHS; never recover H by
            # subtracting identity from a near-identity rounded inverse.
            ka_fused_gp_solve!(h,r⁻⁺,rr,backend)
        else
            system .= I_static .- rr
            canopy_increment_solve!(h, system, rr)
        end
        _bmm!(product, t⁺⁺, h)
        c .= t⁺⁺ .+ h .+ product

        if SFI
            if fused_updates
                ka_canopy_increment_sources!(j₀⁺,j₀⁻,r⁻⁺,c,beta,backend)
            else
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
        end

        # Q=R+RΔ; R'=R+Q+CQ; Δ'=C+Δ+CΔ. Form no near-identity T.
        if fused_updates
            ka_canopy_increment_operators!(r⁻⁺,t⁺⁺,c,beta,backend)
        else
            _bmm!(rr, r⁻⁺, t⁺⁺)
            rr .+= r⁻⁺
            _bmm!(product, c, rr)
            r⁻⁺ .+= rr .+ product
            _bmm!(product, c, t⁺⁺)
            t⁺⁺ .= c .+ t⁺⁺ .+ product
            beta .= beta .* (FT(2) .+ beta)
        end
    end
    t⁺⁺ .+= I_static
    apply_D_matrix!(pol_type.n, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻)
    SFI && apply_D_matrix_SFI!(pol_type.n, j₀⁻)
    synchronize_if_gpu()
    return nothing
end

"Backend fallback for the canopy geometric-series solve."
canopy_increment_solve!(X,A,B) = batch_solve!(X,A,B)

"Use pivoted LU for square CPU systems; preserve the general backend solve API."
function canopy_increment_solve!(X::Array{FT,3},A::Array{FT,3},B::Array{FT,3}) where
        {FT<:LinearAlgebra.BLAS.BlasFloat}
    Threads.@threads for k in axes(A,3)
        @views begin
            copyto!(X[:,:,k],B[:,:,k])
            ldiv!(lu!(A[:,:,k]),X[:,:,k])
        end
    end
    X
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
