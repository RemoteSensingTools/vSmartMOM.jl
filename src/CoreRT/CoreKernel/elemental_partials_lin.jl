# Elemental core partials for the legacy/reference API.

"""
    get_elem_rt!(r⁻⁺, t⁺⁺, ṙ⁻⁺, ṫ⁺⁺, ϖ_λ, dτ_λ, Z⁻⁺, Z⁺⁺, qp_μN, wct)

KernelAbstractions tangent-linear elemental R/T kernel. Each workitem owns one
matrix/spectral element `(i, j, n)`, writes the exact finite-δ forward
reflection/transmission entries, and stores the three local core derivatives
with respect to `(dτ, ϖ, Z)` in the fourth dimension of `ṙ⁻⁺` and `ṫ⁺⁺`.
Zero-weight quadrature columns receive Beer-law diagonal transmission only.
"""
@kernel function get_elem_rt!(r⁻⁺, t⁺⁺,
                        ṙ⁻⁺, ṫ⁺⁺, 
                        @Const(ϖ_λ), @Const(dτ_λ),
                        @Const(Z⁻⁺), @Const(Z⁺⁺),
                        @Const(qp_μN), @Const(wct))
    FT = eltype(r⁻⁺)
    n2 = 1
    i, j, n = @index(Global, NTuple) 
    if size(Z⁻⁺,3)>1
        n2 = n
    end
    
    if (wct[j] > eps(FT)) 
        # 𝐑⁻⁺(μᵢ, μⱼ) = ϖ ̇𝐙⁻⁺(μᵢ, μⱼ) ̇(μⱼ/(μᵢ+μⱼ)) ̇(1 - exp{-τ ̇(1/μᵢ + 1/μⱼ)}) ̇𝑤ⱼ
        # d𝐑⁻⁺(μᵢ, μⱼ)/dτ = ϖ ̇𝐙⁻⁺(μᵢ, μⱼ) ̇(1/μᵢ) ̇exp{-τ ̇(1/μᵢ + 1/μⱼ)}  ̇𝑤ⱼ
        # d𝐑⁻⁺(μᵢ, μⱼ)/dϖ = 𝐙⁻⁺(μᵢ, μⱼ) ̇(μⱼ/(μᵢ+μⱼ)) ̇(1 - exp{-τ ̇(1/μᵢ + 1/μⱼ)}) ̇𝑤ⱼ
        # d𝐑⁻⁺(μᵢ, μⱼ)/dZ = ϖ ̇(μⱼ/(μᵢ+μⱼ)) ̇(1 - exp{-τ ̇(1/μᵢ + 1/μⱼ)}) ̇𝑤ⱼ
        r⁻⁺[i,j,n] = 
            ϖ_λ[n] * Z⁻⁺[i,j,n2] * 
            #Z⁻⁺[i,j] * 
            (qp_μN[j] / (qp_μN[i] + qp_μN[j])) * wct[j] * 
            -expm1(-dτ_λ[n] * ((1 / qp_μN[i]) + (1 / qp_μN[j])))
        # derivative wrt τ_λ
        ṙ⁻⁺[i,j,n,1] = 
            ϖ_λ[n] * Z⁻⁺[i,j,n2] * 
            (1/qp_μN[i]) * wct[j] * 
            exp(-dτ_λ[n] * ((1 / qp_μN[i]) + (1 / qp_μN[j]))) 
        # derivative wrt ϖ
        ṙ⁻⁺[i,j,n,2] = Z⁻⁺[i,j,n2] * (qp_μN[j]/(qp_μN[i]+qp_μN[j])) * wct[j] *
            -expm1(-dτ_λ[n]*(inv(qp_μN[i])+inv(qp_μN[j])))
        # derivative wrt Z
        # derivative wrt Z: direct formula avoids 0/0 when Z=0
        ṙ⁻⁺[i,j,n,3] = ϖ_λ[n] * 
            (qp_μN[j] / (qp_μN[i] + qp_μN[j])) * wct[j] * 
            -expm1(-dτ_λ[n] * ((1 / qp_μN[i]) + (1 / qp_μN[j])))
                    
        if (qp_μN[i] == qp_μN[j])
            # 𝐓⁺⁺(μᵢ, μᵢ) = (exp{-τ/μᵢ}(1 + ϖ ̇𝐙⁺⁺(μᵢ, μᵢ) ̇(τ/μᵢ))) ̇𝑤ᵢ
            # d𝐓⁺⁺(μᵢ, μᵢ)/dτ_λ = (exp{-τ/μⱼ}/μᵢ)⋅(ϖ ̇𝐙⁺⁺(μᵢ, μᵢ)⋅(1-τ/μⱼ)-1) ̇𝑤ⱼ  
            # d𝐓⁺⁺(μᵢ, μᵢ)/dϖ_λ = 𝐙⁺⁺(μᵢ, μᵢ)⋅(τ/μᵢ) ̇exp{-τ/μᵢ} ̇𝑤ᵢ
            # d𝐓⁺⁺(μᵢ, μᵢ)/dZ   = ϖ ̇(τ/μᵢ) ̇exp{-τ/μᵢ} ̇𝑤ᵢ
            if i == j
                t⁺⁺[i,j,n] = 
                    exp(-dτ_λ[n] / qp_μN[i]) *
                    (1 + ϖ_λ[n] * Z⁺⁺[i,i,n2] * (dτ_λ[n] / qp_μN[i]) * wct[i])
                # derivative wrt τ_λ
                ṫ⁺⁺[i,j,n,1] = 
                    exp(-dτ_λ[n] / qp_μN[i]) * (1 / qp_μN[i]) *
                    (-1 + ϖ_λ[n] * Z⁺⁺[i,i,n2] * wct[i] * (1 - dτ_λ[n] / qp_μN[i]))
                # derivative wrt ϖ_λ
                ṫ⁺⁺[i,j,n,2] = 
                    exp(-dτ_λ[n] / qp_μN[i]) *
                    Z⁺⁺[i,i,n2] * (dτ_λ[n] / qp_μN[i]) * wct[i]    
                # derivative wrt Z
                ṫ⁺⁺[i,j,n,3] = 
                    exp(-dτ_λ[n] / qp_μN[i]) *
                    ϖ_λ[n] * (dτ_λ[n] / qp_μN[i]) * wct[i]
            else
                # 𝐓⁺⁺(μᵢ, μⱼ) = (exp{-τ/μⱼ}(ϖ ̇𝐙⁺⁺(μᵢ, μⱼ) ̇(τ/μᵢ))) ̇𝑤ⱼ        
                # d𝐓⁺⁺(μᵢ, μⱼ)/dτ_λ = (exp{-τ/μⱼ}⋅ϖ ̇𝐙⁺⁺(μᵢ, μᵢ)/μᵢ)⋅(1 - τ/μⱼ) ̇𝑤ⱼ
                # d𝐓⁺⁺(μᵢ, μᵢ)/dϖ_λ = 𝐙⁺⁺(μᵢ, μᵢ)⋅(τ/μᵢ) ̇exp{-τ/μᵢ} ̇𝑤ᵢ
                # d𝐓⁺⁺(μᵢ, μᵢ)/dZ   = ϖ ̇(τ/μᵢ) ̇exp{-τ/μᵢ} ̇𝑤ᵢ
                t⁺⁺[i,j,n] = exp(-dτ_λ[n] / qp_μN[j]) *
                    (ϖ_λ[n] * Z⁺⁺[i,j,n2] * (dτ_λ[n] / qp_μN[i]) * wct[j])
                # derivative wrt τ_λ
                ṫ⁺⁺[i,j,n,1] = (exp(-dτ_λ[n] / qp_μN[j]) *
                        ϖ_λ[n] * Z⁺⁺[i,j,n2] / qp_μN[i]) * 
                        (1 - dτ_λ[n] / qp_μN[j]) * wct[j]
                # derivative wrt ϖ_λ
                ṫ⁺⁺[i,j,n,2] = exp(-dτ_λ[n]/qp_μN[j]) *
                    Z⁺⁺[i,j,n2] * (dτ_λ[n]/qp_μN[i]) * wct[j]
                # derivative wrt Z
                # derivative wrt Z: direct formula avoids 0/0
                ṫ⁺⁺[i,j,n,3] = exp(-dτ_λ[n] / qp_μN[j]) *
                    ϖ_λ[n] * (dτ_λ[n] / qp_μN[i]) * wct[j]
            end
        else
    
            # 𝐓⁺⁺(μᵢ, μⱼ) = ϖ ̇𝐙⁺⁺(μᵢ, μⱼ) ̇(μⱼ/(μᵢ-μⱼ)) ̇(exp{-τ/μᵢ} - exp{-τ/μⱼ}) ̇𝑤ⱼ
            # d𝐓⁺⁺(μᵢ, μⱼ)/dτ_λ = -ϖ ̇𝐙⁺⁺(μᵢ, μⱼ) ̇(μⱼ/(μᵢ-μⱼ)) ̇(exp{-τ/μᵢ}/μᵢ - exp{-τ/μⱼ}/μⱼ) ̇𝑤ⱼ
            # (𝑖 ≠ 𝑗)
            t⁺⁺[i,j,n] = 
                ϖ_λ[n] * Z⁺⁺[i,j,n2] * 
                #Z⁺⁺[i,j] * 
                (qp_μN[j] / (qp_μN[i] - qp_μN[j])) * wct[j] * 
                expdiff_neg(dτ_λ[n] / qp_μN[i], dτ_λ[n] / qp_μN[j])
            # derivative wrt τ_λ
            ṫ⁺⁺[i,j,n,1] = -ϖ_λ[n] * Z⁺⁺[i,j,n2] * 
                (qp_μN[j] / (qp_μN[i] - qp_μN[j])) * wct[j] * 
                (exp(-dτ_λ[n] / qp_μN[i])/ qp_μN[i] - 
                exp(-dτ_λ[n] / qp_μN[j])/ qp_μN[j]) 
            # derivative wrt ϖ_λ
            ṫ⁺⁺[i,j,n,2] = Z⁺⁺[i,j,n2] * (qp_μN[j]/(qp_μN[i]-qp_μN[j])) * wct[j] *
                expdiff_neg(dτ_λ[n]/qp_μN[i],dτ_λ[n]/qp_μN[j])
            # derivative wrt Z
            # derivative wrt Z: direct formula avoids 0/0
            ṫ⁺⁺[i,j,n,3] = ϖ_λ[n] * 
                (qp_μN[j] / (qp_μN[i] - qp_μN[j])) * wct[j] * 
                expdiff_neg(dτ_λ[n] / qp_μN[i], dτ_λ[n] / qp_μN[j])
        end
    else
        #r⁻⁺[i,j,n] = 0.0
        #ṙ⁻⁺[i,j,n,:] = 0.0
        if i==j
            t⁺⁺[i,j,n] = exp(-dτ_λ[n] / qp_μN[i]) #Suniti
            # derivative wrt τ_λ
            ṫ⁺⁺[i,j,n,1] = -exp(-dτ_λ[n] / qp_μN[i]) / qp_μN[i]
        #else
        #    t⁺⁺[i,j,n] = 0.0
            # derivative wrt τ_λ
        #    ṫ⁺⁺[i,j,n,1] = 0.0
        end
        # derivative wrt ϖ_λ
        #ṫ⁺⁺[i,j,n,2] = 0.0
        # derivative wrt Z
        #ṫ⁺⁺[i,j,n,3] = 0.0
    end
    nothing
end

"""
    get_elem_rt_SFI!(J₀⁺, J₀⁻, J̇₀⁺, J̇₀⁻, ϖ_λ, dτ_λ, τ_sum, τ̇_sum,
                     Z⁻⁺, Z⁺⁺, F₀, qp_μN, ndoubl, wct02, nStokes, I₀, iμ0, D)

KernelAbstractions tangent-linear elemental source-function kernel. Each
workitem computes the direct-beam source vectors for one stream/spectral point
and stores the three core derivatives of those source terms with respect to
`(dτ, ϖ, Z)`. Beam attenuation from the optical depth above the layer is
included in both the forward and derivative outputs.
"""
@kernel function get_elem_rt_SFI!(J₀⁺, J₀⁻, 
                J̇₀⁺, J̇₀⁻, 
                @Const(ϖ_λ), @Const(dτ_λ),
                @Const(τ_sum), @Const(τ̇_sum),
                @Const(Z⁻⁺), @Const(Z⁺⁺), @Const(F₀),
                @Const(qp_μN), ndoubl, wct02, nStokes,
                @Const(I₀), iμ0, @Const(D))
    i_start  = nStokes*(iμ0-1) + 1 
    i_end    = nStokes*iμ0
    
    i, _, n = @index(Global, NTuple) ##Suniti: What are Global and Ntuple?
    FT = eltype(I₀)
    #J₀⁺[i, 1, n]=0
    #J₀⁻[i, 1, n]=0
    #J̇₀⁺[i, 1, n, 1:3]=0
    #J̇₀⁻[i, 1, n, 1:3]=0
    n2=1
    if size(Z⁻⁺,3)>1
        n2 = n
    end
    
    Z⁺⁺_I₀ = zero(FT);
    Z⁻⁺_I₀ = zero(FT);
    
    for ii = i_start:i_end
        # `n2` follows the Z spectral axis, which may be length 1 on flat-Z
        # fast paths.  F₀ remains band-resolved and must be sampled at `n`.
        Z⁺⁺_I₀ += Z⁺⁺[i,ii,n2] * F₀[ii-i_start+1,n] #I₀[ii-i_start+1]
        Z⁻⁺_I₀ += Z⁻⁺[i,ii,n2] * F₀[ii-i_start+1,n] #I₀[ii-i_start+1] 
    end

    # Direct products preserve the τ=0, ϖ=0 and projected-Z=0 limits.
    f⁺,f⁻,df⁺,df⁻ = _single_scatter_source_factors(dτ_λ[n],qp_μN[i],qp_μN[i_start])
    J₀⁺[i,1,n] = wct02 * ϖ_λ[n] * Z⁺⁺_I₀ * f⁺
    J₀⁻[i,1,n] = wct02 * ϖ_λ[n] * Z⁻⁺_I₀ * f⁻
    J̇₀⁺[i,1,n,1] = wct02 * ϖ_λ[n] * Z⁺⁺_I₀ * df⁺
    J̇₀⁻[i,1,n,1] = wct02 * ϖ_λ[n] * Z⁻⁺_I₀ * df⁻
    J̇₀⁺[i,1,n,2] = wct02 * Z⁺⁺_I₀ * f⁺
    J̇₀⁻[i,1,n,2] = wct02 * Z⁻⁺_I₀ * f⁻
    J̇₀⁺[i,1,n,3] = wct02 * ϖ_λ[n] * f⁺
    J̇₀⁻[i,1,n,3] = wct02 * ϖ_λ[n] * f⁻

    # TODO: Move this out until after doubling (it is not necessary to consider this here already if Raman scattering is not involved)
    J₀⁺[i, 1, n] *= exp(-τ_sum[n]/qp_μN[i_start])
    J₀⁻[i, 1, n] *= exp(-τ_sum[n]/qp_μN[i_start])

    # Bug 22 fix: Remove τ̇_sum[1,n] contribution from core derivative.
    # The τ̇_sum beam attenuation derivative is per-physical-parameter and must be
    # added AFTER the chain rule (in rt_kernel!), not here in the 3-core framework.
    # Old code used τ̇_sum[1,n] which only captured parameter 1's contribution.
    J̇₀⁺[i, 1, n, 1] = J̇₀⁺[i, 1, n, 1]*exp(-τ_sum[n]/qp_μN[i_start])
    J̇₀⁻[i, 1, n, 1] = J̇₀⁻[i, 1, n, 1]*exp(-τ_sum[n]/qp_μN[i_start])
    J̇₀⁺[i, 1, n, 2] = J̇₀⁺[i, 1, n, 2]*exp(-τ_sum[n]/qp_μN[i_start]) #+
                        #J₀⁺[i, 1, n] * (-τ̇_sum[1,n]/qp_μN[i_start])
    J̇₀⁻[i, 1, n, 2] = J̇₀⁻[i, 1, n, 2]*exp(-τ_sum[n]/qp_μN[i_start]) #+
                        #J₀⁻[i, 1, n] * (-τ̇_sum[1,n]/qp_μN[i_start])
    J̇₀⁺[i, 1, n, 3] = J̇₀⁺[i, 1, n, 3]*exp(-τ_sum[n]/qp_μN[i_start]) #+
                        #J₀⁺[i, 1, n] * (-τ̇_sum[1,n]/qp_μN[i_start])
    J̇₀⁻[i, 1, n, 3] = J̇₀⁻[i, 1, n, 3]*exp(-τ_sum[n]/qp_μN[i_start]) #+
                        #J₀⁻[i, 1, n] * (-τ̇_sum[1,n]/qp_μN[i_start])


    if ndoubl >= 1
        J₀⁻[i, 1, n] = D[i,i]*J₀⁻[i, 1, n] #D = Diagonal{1,1,-1,-1,...Nquad times}
        J̇₀⁻[i, 1, n, 1] = D[i,i]*J̇₀⁻[i, 1, n, 1]
        J̇₀⁻[i, 1, n, 2] = D[i,i]*J̇₀⁻[i, 1, n, 2]
        J̇₀⁻[i, 1, n, 3] = D[i,i]*J̇₀⁻[i, 1, n, 3]
    end  
    #if (n==840||n==850)
    #    @show i, n, J₀⁺[i, 1, n], J₀⁻[i, 1, n]
    #end
    nothing
end

