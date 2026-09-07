"""
Finite-thickness solar-source factors from SF2023-II (11) and their thickness
partials. With j± = wₘ ϖ (Z± F₀) f±, the other partials follow by removing
one linear factor: ∂j±/∂ϖ = wₘ (Z± F₀) f± and ∂j±/∂(Z±F₀) = wₘ ϖ f±.
Evaluate these products directly: dividing j by a vanishing factor would
incorrectly erase a nonzero derivative. The equal-μ limit is regular at τ=0.
"""
@inline function _single_scatter_source_factors(τ, μᵢ, μ₀)
    if μᵢ == μ₀
        x = τ/μ₀
        f⁺ = x * exp(-x)
        df⁺ = exp(-x) * (one(x)-x) / μ₀
    else
        c = μ₀/(μᵢ-μ₀)
        f⁺ = c * expdiff_neg(τ/μᵢ,τ/μ₀)
        df⁺ = c * (-exp(-τ/μᵢ)/μᵢ + exp(-τ/μ₀)/μ₀)
    end
    q = inv(μᵢ)+inv(μ₀)
    c = μ₀/(μᵢ+μ₀)
    f⁻ = c * (-expm1(-τ*q))
    df⁻ = c * exp(-τ*q) * q
    return f⁺,f⁻,df⁺,df⁻
end

# Finite-thickness elemental evaluation and supplied-direction chain rule.
# ============================================================================
# Fused kernels: combine elemental RT + chain rule in a single pass.
# These eliminate the separate lin_added_layer_all_params! call by computing
# ap_ṙ⁻⁺, ap_ṫ⁺⁺ (and optionally ap_ṙ⁺⁻, ap_ṫ⁻⁻) directly from local
# 3-core scalar intermediates, avoiding ~12 full-array reads.
# ============================================================================

"""
    get_elem_rt_fused!(...)

Fused elemental R/T kernel: computes forward r⁻⁺, t⁺⁺ and their per-parameter
derivatives ap_ṙ⁻⁺, ap_ṫ⁺⁺ (and ap_ṙ⁺⁻, ap_ṫ⁻⁻ for ndoubl < 1) in a single pass.

The 3-core derivatives (ṙ⁻⁺[1:3], ṫ⁺⁺[1:3]) are kept as local scalars and used
directly for the chain rule, then also written to their arrays for backward
compatibility with the 3-core doubling path.
"""
@kernel function get_elem_rt_fused!(r⁻⁺, t⁺⁺,
                        ṙ⁻⁺, ṫ⁺⁺,
                        ap_ṙ⁻⁺, ap_ṫ⁺⁺, ap_ṙ⁺⁻, ap_ṫ⁻⁻,
                        @Const(ϖ_λ), @Const(dτ_λ),
                        @Const(Z⁻⁺), @Const(Z⁺⁺),
                        @Const(dτ̇), @Const(ϖ̇),
                        @Const(Ż⁻⁺), @Const(Ż⁺⁺_lin),
                        @Const(qp_μN), @Const(wct),
                        nparams, ndoubl, pol_n)
    FT = eltype(r⁻⁺)
    i, j, n = @index(Global, NTuple)
    n2 = 1
    if size(Z⁻⁺, 3) > 1
        n2 = n
    end
    n2_lin = 1
    if size(Ż⁻⁺, 3) > 1
        n2_lin = n
    end

    # D-matrix Stokes signs
    i_stokes = mod(i, pol_n)
    j_stokes = mod(j, pol_n)
    i12 = (pol_n == 1) | ((1 <= i_stokes) & (i_stokes <= 2))
    j12 = (pol_n == 1) | ((1 <= j_stokes) & (j_stokes <= 2))
    same_block = (i12 & j12) | (!i12 & !j12)
    d_sign = ifelse(same_block, one(FT), -one(FT))
    di = ifelse(i12, one(FT), -one(FT))
    dj = ifelse(j12, one(FT), -one(FT))

    # R⁻⁺ row-sign correction for elemental D-matrix (ndoubl >= 1 negates Stokes 3,4 rows)
    sign_r = ifelse((ndoubl >= 1) & !i12, -one(FT), one(FT))

    # Local 3-core derivative scalars
    ṙ_tau = FT(0); ṙ_w = FT(0); ṙ_Z = FT(0)
    ṫ_tau = FT(0); ṫ_w = FT(0); ṫ_Z = FT(0)

    if (wct[j] > eps(FT))
        # SF2023-II (10), with Fourier quadrature weight wct[j]:
        # rᵢⱼ = ϖ Zᵢⱼ cᵢⱼ (1-exp(-aᵢⱼ dτ)),
        # cᵢⱼ = μⱼ wct[j]/(μᵢ+μⱼ), aᵢⱼ = 1/μᵢ+1/μⱼ.
        # ∂r/∂dτ = ϖ Zᵢⱼ wct[j] exp(-aᵢⱼ dτ)/μᵢ;
        # ∂r/∂ϖ and ∂r/∂Z remove the corresponding linear factor.
        # These Z partials are element-local ONLY before doubling.
        r⁻⁺[i,j,n] =
            ϖ_λ[n] * Z⁻⁺[i,j,n2] *
            (qp_μN[j] / (qp_μN[i] + qp_μN[j])) * wct[j] *
            -expm1(-dτ_λ[n] * ((1 / qp_μN[i]) + (1 / qp_μN[j])))

        ṙ_tau = ϖ_λ[n] * Z⁻⁺[i,j,n2] *
            (1/qp_μN[i]) * wct[j] *
            exp(-dτ_λ[n] * ((1 / qp_μN[i]) + (1 / qp_μN[j])))
        ṙ_w = Z⁻⁺[i,j,n2] * (qp_μN[j] / (qp_μN[i] + qp_μN[j])) * wct[j] *
            -expm1(-dτ_λ[n] * ((1 / qp_μN[i]) + (1 / qp_μN[j])))
        ṙ_Z = ϖ_λ[n] *
            (qp_μN[j] / (qp_μN[i] + qp_μN[j])) * wct[j] *
            -expm1(-dτ_λ[n] * ((1 / qp_μN[i]) + (1 / qp_μN[j])))

        # Write 3-core (backward compat)
        ṙ⁻⁺[i,j,n,1] = ṙ_tau
        ṙ⁻⁺[i,j,n,2] = ṙ_w
        ṙ⁻⁺[i,j,n,3] = ṙ_Z

        # SF2023-II (10): transmission adds δᵢⱼ exp(-dτ/μᵢ) to
        # ϖ Zᵢⱼ wct[j] μⱼ/(μᵢ-μⱼ) [exp(-dτ/μᵢ)-exp(-dτ/μⱼ)].
        # For equal μ, use its finite limit (dτ/μ) exp(-dτ/μ).
        # Different Stokes rows can share μ without sharing δᵢⱼ.
        if (qp_μN[i] == qp_μN[j])
            if i == j
                t⁺⁺[i,j,n] =
                    exp(-dτ_λ[n] / qp_μN[i]) *
                    (1 + ϖ_λ[n] * Z⁺⁺[i,i,n2] * (dτ_λ[n] / qp_μN[i]) * wct[i])
                ṫ_tau =
                    exp(-dτ_λ[n] / qp_μN[i]) * (1 / qp_μN[i]) *
                    (-1 + ϖ_λ[n] * Z⁺⁺[i,i,n2] * wct[i] * (1 - dτ_λ[n] / qp_μN[i]))
                ṫ_w =
                    exp(-dτ_λ[n] / qp_μN[i]) *
                    Z⁺⁺[i,i,n2] * (dτ_λ[n] / qp_μN[i]) * wct[i]
                ṫ_Z =
                    exp(-dτ_λ[n] / qp_μN[i]) *
                    ϖ_λ[n] * (dτ_λ[n] / qp_μN[i]) * wct[i]
            else
                t⁺⁺[i,j,n] = exp(-dτ_λ[n] / qp_μN[j]) *
                    (ϖ_λ[n] * Z⁺⁺[i,j,n2] * (dτ_λ[n] / qp_μN[i]) * wct[j])
                ṫ_tau = (exp(-dτ_λ[n] / qp_μN[j]) *
                        ϖ_λ[n] * Z⁺⁺[i,j,n2] / qp_μN[i]) *
                        (1 - dτ_λ[n] / qp_μN[j]) * wct[j]
                ṫ_w = exp(-dτ_λ[n] / qp_μN[j]) *
                    Z⁺⁺[i,j,n2] * (dτ_λ[n] / qp_μN[i]) * wct[j]
                ṫ_Z = exp(-dτ_λ[n] / qp_μN[j]) *
                    ϖ_λ[n] * (dτ_λ[n] / qp_μN[i]) * wct[j]
            end
        else
            t⁺⁺[i,j,n] =
                ϖ_λ[n] * Z⁺⁺[i,j,n2] *
                (qp_μN[j] / (qp_μN[i] - qp_μN[j])) * wct[j] *
                expdiff_neg(dτ_λ[n] / qp_μN[i], dτ_λ[n] / qp_μN[j])
            ṫ_tau = -ϖ_λ[n] * Z⁺⁺[i,j,n2] *
                (qp_μN[j] / (qp_μN[i] - qp_μN[j])) * wct[j] *
                (exp(-dτ_λ[n] / qp_μN[i])/ qp_μN[i] -
                exp(-dτ_λ[n] / qp_μN[j])/ qp_μN[j])
            ṫ_w = Z⁺⁺[i,j,n2] * (qp_μN[j] / (qp_μN[i] - qp_μN[j])) * wct[j] *
                expdiff_neg(dτ_λ[n] / qp_μN[i], dτ_λ[n] / qp_μN[j])
            ṫ_Z = ϖ_λ[n] *
                (qp_μN[j] / (qp_μN[i] - qp_μN[j])) * wct[j] *
                expdiff_neg(dτ_λ[n] / qp_μN[i], dτ_λ[n] / qp_μN[j])
        end

        # Write 3-core (backward compat)
        ṫ⁺⁺[i,j,n,1] = ṫ_tau
        ṫ⁺⁺[i,j,n,2] = ṫ_w
        ṫ⁺⁺[i,j,n,3] = ṫ_Z

        # S2014 (C.25)–(C.26), applied before matrix products mix indices:
        # drᵢⱼ/dp = (∂rᵢⱼ/∂dτ) dτ/dp + (∂rᵢⱼ/∂ϖ) dϖ/dp
        #           + (∂rᵢⱼ/∂Zᵢⱼ) dZᵢⱼ/dp; likewise for t.
        # This contraction also accepts arbitrary supplied core-optics
        # directions; p need not denote an aerosol or gas parameter.
        for iparam = 1:nparams
            val_r = ṙ_tau * dτ̇[n,iparam] + ṙ_w * ϖ̇[n,iparam] + ṙ_Z * Ż⁻⁺[i,j,n2_lin,iparam]
            val_t = ṫ_tau * dτ̇[n,iparam] + ṫ_w * ϖ̇[n,iparam] + ṫ_Z * Ż⁺⁺_lin[i,j,n2_lin,iparam]

            ap_ṙ⁻⁺[i,j,n,iparam] = sign_r * val_r
            ap_ṫ⁺⁺[i,j,n,iparam] = val_t

            if ndoubl < 1
                # For ndoubl < 1 (no doubling): compute ⁺⁻ and ⁻⁻ via D-matrix
                # ṙ⁺⁻ = d_sign * (ṙ_tau*dτ̇ + ṙ_w*ϖ̇ + ṙ_Z * D·Ż⁻⁺·D)
                # where D·Ż·D at (i,j) = di*dj*Ż
                ap_ṙ⁺⁻[i,j,n,iparam] = d_sign * (ṙ_tau * dτ̇[n,iparam] + ṙ_w * ϖ̇[n,iparam] + ṙ_Z * di * dj * Ż⁻⁺[i,j,n2_lin,iparam])
                ap_ṫ⁻⁻[i,j,n,iparam] = d_sign * (ṫ_tau * dτ̇[n,iparam] + ṫ_w * ϖ̇[n,iparam] + ṫ_Z * di * dj * Ż⁺⁺_lin[i,j,n2_lin,iparam])
            end
        end
    else
        # No scattering weight: only diagonal transmission
        if i == j
            t⁺⁺[i,j,n] = exp(-dτ_λ[n] / qp_μN[i])
            ṫ_tau = -exp(-dτ_λ[n] / qp_μN[i]) / qp_μN[i]
            ṫ⁺⁺[i,j,n,1] = ṫ_tau
            for iparam = 1:nparams
                val_t = ṫ_tau * dτ̇[n,iparam]
                ap_ṫ⁺⁺[i,j,n,iparam] = val_t
                if ndoubl < 1
                    # diagonal: same_block=true so d_sign=1
                    ap_ṫ⁻⁻[i,j,n,iparam] = val_t
                end
            end
        end
    end
    nothing
end

"""
    get_elem_rt_SFI_fused!(...)

Fused SFI source kernel: computes J₀⁺, J₀⁻ and their per-parameter derivatives
ap_J̇₀⁺, ap_J̇₀⁻ in a single pass, including above-layer beam attenuation.

Eliminates the separate chain-rule pass for SFI terms and the per-parameter
τ̇_sum correction loop in rt_kernel!.
"""
@kernel function get_elem_rt_SFI_fused!(J₀⁺, J₀⁻,
                J̇₀⁺, J̇₀⁻,
                ap_J̇₀⁺, ap_J̇₀⁻,
                @Const(ϖ_λ), @Const(dτ_λ),
                @Const(τ_sum), @Const(τ̇_sum),
                @Const(Z⁻⁺), @Const(Z⁺⁺), @Const(F₀),
                @Const(dτ̇), @Const(ϖ̇),
                @Const(Ż⁻⁺), @Const(Ż⁺⁺_lin),
                @Const(qp_μN), ndoubl, wct02, nStokes,
                @Const(I₀), μ0, iμ0, @Const(D), nparams)
    i_start  = nStokes*(iμ0-1) + 1
    i_end    = nStokes*iμ0

    i, _, n = @index(Global, NTuple)
    FT = eltype(I₀)

    n2 = 1
    if size(Z⁻⁺, 3) > 1
        n2 = n
    end
    n2_lin = 1
    if size(Ż⁻⁺, 3) > 1
        n2_lin = n
    end

    # Forward Z·I₀ products
    Z⁺⁺_I₀ = zero(FT)
    Z⁻⁺_I₀ = zero(FT)
    for ii = i_start:i_end
        # `n2` follows the Z spectral axis, which may be length 1 on flat-Z
        # fast paths.  F₀ remains band-resolved and must be sampled at `n`.
        Z⁺⁺_I₀ += Z⁺⁺[i,ii,n2] * F₀[ii-i_start+1,n]
        Z⁻⁺_I₀ += Z⁻⁺[i,ii,n2] * F₀[ii-i_start+1,n]
    end

    # SF2023-II (11), direct solar term at the top of this layer:
    # j± = wₘ ϖ (Z± F₀) f±(dτ), where
    # f⁺ = μ₀/(μᵢ-μ₀) [exp(-dτ/μᵢ)-exp(-dτ/μ₀)],
    # f⁻ = μ₀/(μᵢ+μ₀) [1-exp(-dτ(1/μᵢ+1/μ₀))].
    # At μᵢ=μ₀, f⁺ = (dτ/μ₀) exp(-dτ/μ₀). The three local
    # partials differentiate dτ, ϖ, and the projected phase column Z F₀.
    f⁺,f⁻,df⁺,df⁻ = _single_scatter_source_factors(dτ_λ[n],qp_μN[i],μ0)
    J₀⁺[i,1,n] = wct02 * ϖ_λ[n] * Z⁺⁺_I₀ * f⁺
    J₀⁻[i,1,n] = wct02 * ϖ_λ[n] * Z⁻⁺_I₀ * f⁻
    J̇⁺_tau = wct02 * ϖ_λ[n] * Z⁺⁺_I₀ * df⁺
    J̇⁻_tau = wct02 * ϖ_λ[n] * Z⁻⁺_I₀ * df⁻
    J̇⁺_w = wct02 * Z⁺⁺_I₀ * f⁺
    J̇⁻_w = wct02 * Z⁻⁺_I₀ * f⁻
    J̇⁺_Z = wct02 * ϖ_λ[n] * f⁺
    J̇⁻_Z = wct02 * ϖ_λ[n] * f⁻

    # ---- Apply beam attenuation exp(-τ_sum/μ₀) ----
    beam_atten = exp(-τ_sum[n]/μ0)
    J₀⁺[i, 1, n] *= beam_atten
    J₀⁻[i, 1, n] *= beam_atten
    J̇⁺_tau *= beam_atten
    J̇⁺_w   *= beam_atten
    J̇⁺_Z   *= beam_atten
    J̇⁻_tau *= beam_atten
    J̇⁻_w   *= beam_atten
    J̇⁻_Z   *= beam_atten

    # Write 3-core arrays (backward compat)
    J̇₀⁺[i, 1, n, 1] = J̇⁺_tau
    J̇₀⁺[i, 1, n, 2] = J̇⁺_w
    J̇₀⁺[i, 1, n, 3] = J̇⁺_Z
    J̇₀⁻[i, 1, n, 1] = J̇⁻_tau
    J̇₀⁻[i, 1, n, 2] = J̇⁻_w
    J̇₀⁻[i, 1, n, 3] = J̇⁻_Z

    # ---- D-matrix for J₀⁻ (ndoubl >= 1) ----
    if ndoubl >= 1
        J₀⁻[i, 1, n] = D[i,i]*J₀⁻[i, 1, n]
        J̇⁻_tau = D[i,i]*J̇⁻_tau
        J̇⁻_w   = D[i,i]*J̇⁻_w
        J̇⁻_Z   = D[i,i]*J̇⁻_Z
        J̇₀⁻[i, 1, n, 1] = J̇⁻_tau
        J̇₀⁻[i, 1, n, 2] = J̇⁻_w
        J̇₀⁻[i, 1, n, 3] = J̇⁻_Z
    end

    # Complete solar tangent, by S2014 (C.6) applied to SF2023-II (11):
    # dj±/dp = e_top [j±_τ dτ/dp + j±_ϖ dϖ/dp + j±_Z (dZ±/dp) F₀]
    #          - j±_attenuated (dτ_sum/dp)/μ₀.
    # Above, e_top has already been applied to each local partial. The
    # last term differentiates OVERLYING attenuation, distinct from the
    # current layer's thickness derivative; both are required for SFI.
    for iparam = 1:nparams
        # Compute Ż·I₀ dot products for this parameter
        Ż⁺⁺_I₀_p = FT(0)
        Ż⁻⁺_I₀_p = FT(0)
        for ii = i_start:i_end
            Ż⁺⁺_I₀_p += Ż⁺⁺_lin[i, ii, n2_lin, iparam] * F₀[ii-i_start+1, n]
            Ż⁻⁺_I₀_p += Ż⁻⁺[i, ii, n2_lin, iparam] * F₀[ii-i_start+1, n]
        end

        # Chain rule: ap_J̇ = J̇_tau*dτ̇ + J̇_w*ϖ̇ + J̇_Z*Ż_I₀
        ap_J̇₀⁺[i, 1, n, iparam] = J̇⁺_tau * dτ̇[n,iparam] + J̇⁺_w * ϖ̇[n,iparam] + J̇⁺_Z * Ż⁺⁺_I₀_p
        ap_J̇₀⁻[i, 1, n, iparam] = J̇⁻_tau * dτ̇[n,iparam] + J̇⁻_w * ϖ̇[n,iparam] + J̇⁻_Z * Ż⁻⁺_I₀_p

        # Per-direction above-layer beam attenuation derivative
        # d(exp(-τ_sum/μ₀))/dp_j * J₀ = -τ̇_sum[j]/μ₀ * J₀
        ap_J̇₀⁺[i, 1, n, iparam] += J₀⁺[i, 1, n] * (-τ̇_sum[n, iparam] / μ0)
        ap_J̇₀⁻[i, 1, n, iparam] += J₀⁻[i, 1, n] * (-τ̇_sum[n, iparam] / μ0)
    end

    nothing
end

"""
Tangent-linear counterpart of `get_elem_rt_solar_columns!`.

The output layout is `(NquadN, nStokes, nSpec, nParams)`.  Derivatives include
the elemental optical-depth, single-scatter-albedo, and phase-column terms,
but intentionally exclude attenuation above the layer; that final chain-rule
term is applied when the columns are contracted into `J₀±`.
"""
@kernel function get_elem_rt_solar_columns_lin!(Ṙ₀⁻⁺, Ṙ₀⁺⁻, Ṫ₀⁺⁺, Ṫ₀⁻⁻,
                                                 @Const(ϖ), @Const(dτ),
                                                 @Const(Z₀⁻), @Const(Z₀⁺),
                                                 @Const(dτ̇), @Const(ϖ̇),
                                                 @Const(Ż₀⁻), @Const(Ż₀⁺),
                                                 @Const(μ), μ₀, wct02,
                                                 nStokes, @Const(Dpol))
    i, s, n, p = @index(Global, NTuple)
    iz = size(Z₀⁺,3) == 1 ? 1 : n
    izd = size(Ż₀⁺,3) == 1 ? 1 : n
    μᵢ = μ[i]
    f_t,f_r,df_t,df_r = _single_scatter_source_factors(dτ[n],μᵢ,μ₀)
    common_t = ϖ̇[n,p]*Z₀⁺[i,s,iz] + ϖ[n]*Ż₀⁺[i,s,izd,p]
    common_r = ϖ̇[n,p]*Z₀⁻[i,s,iz] + ϖ[n]*Ż₀⁻[i,s,izd,p]
    td = wct02 * (common_t*f_t + ϖ[n]*Z₀⁺[i,s,iz]*df_t*dτ̇[n,p])
    rd = wct02 * (common_r*f_r + ϖ[n]*Z₀⁻[i,s,iz]*df_r*dτ̇[n,p])
    Ṫ₀⁺⁺[i,s,n,p] = td
    Ṙ₀⁻⁺[i,s,n,p] = rd
    parity = Dpol[mod1(i,nStokes)] * Dpol[s]
    Ṫ₀⁻⁻[i,s,n,p] = parity * td
    Ṙ₀⁺⁻[i,s,n,p] = parity * rd
    nothing
end

"""Contract forward/tangent solar columns and add the above-layer extinction tangent."""
@kernel function apply_solar_columns_lin!(J₀⁺, J₀⁻, J̇₀⁺, J̇₀⁻,
                                          @Const(R₀⁻⁺), @Const(T₀⁺⁺),
                                          @Const(Ṙ₀⁻⁺), @Const(Ṫ₀⁺⁺),
                                          @Const(F₀), @Const(τ_above), @Const(τ̇_above),
                                          μ₀, ndoubl, nStokes, @Const(Dpol))
    i, _, n, p = @index(Global, NTuple)
    FT = eltype(J̇₀⁺)
    down = zero(FT); up = zero(FT)
    down_dot = zero(FT); up_dot = zero(FT)
    for s in 1:size(R₀⁻⁺,2)
        f = F₀[s,n]
        down += T₀⁺⁺[i,s,n] * f
        up += R₀⁻⁺[i,s,n] * f
        down_dot += Ṫ₀⁺⁺[i,s,n,p] * f
        up_dot += Ṙ₀⁻⁺[i,s,n,p] * f
    end
    beam = exp(-τ_above[n]/μ₀)
    jplus = down * beam
    jminus = up * beam
    jdplus = beam * (down_dot - down*τ̇_above[n,p]/μ₀)
    jdminus = beam * (up_dot - up*τ̇_above[n,p]/μ₀)
    if ndoubl >= 1
        parity = Dpol[mod1(i,nStokes)]
        jminus *= parity
        jdminus *= parity
    end
    # Avoid a same-address write race across the parameter workitems on GPU.
    if p == 1
        J₀⁺[i,1,n] = jplus
        J₀⁻[i,1,n] = jminus
    end
    J̇₀⁺[i,1,n,p] = jdplus
    J̇₀⁻[i,1,n,p] = jdminus
    nothing
end
