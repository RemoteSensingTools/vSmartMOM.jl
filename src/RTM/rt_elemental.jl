"Elemental single-scattering layer"
function rt_elemental_helper!(pol_type, dτ, ϖ, Z⁺⁺, Z⁻⁺, m, 
                              ndoubl, scatter, qp_μ, wt_μ, 
                              r⁻⁺::AbstractArray{FT,3}, 
                              t⁺⁺::AbstractArray{FT,3}, 
                              r⁺⁻::AbstractArray{FT,3}, 
                              t⁻⁻::AbstractArray{FT,3}, 
                              D::AbstractArray{FT,3},
                              I_static::AbstractArray) where {FT}

    # ToDo: Main output is r⁺⁻, r⁻⁺, t⁻⁻, t⁺⁺ (can be renamed to t⁺⁺, etc)
    # Need to check with paper nomenclature. This is basically eqs. 19-20 in vSmartMOM

    # dτ: optical depth of elemental layer
    # ϖ: single scattering albedo of elemental layer
    # bb: thermal source function at the upper boundary of the elemental layer
    # m: fourier moment
    # n: layer of which this is an elemental
    # ndoubl: number of doubling computations needed to progress from the elemental layer 
    #         to the full homogeneous layer n
    # scatter: flag indicating scattering

    nSpec = size(r⁻⁺, 3)

    # Z⁺⁺_ = repeat(Z⁺⁺, 1, 1, nSpec)
    # Z⁻⁺_ = repeat(Z⁻⁺, 1, 1, nSpec)

    if scatter

        #TODO: import vector containing quadrature cosines qp_μ of length Nquad4
        #TODO: import vector containing quadrature weights wt_μ of length Nquad4
        #TODO: construct composite, post-truncation dτ=τ/2^{ndoubl} , ϖ, p⁺⁺, p⁻⁺ matrices and import them here

        qp_μ4 = reduce(vcat, (fill.(qp_μ,[pol_type.n])))
        wt_μ4 = reduce(vcat, (fill.(wt_μ,[pol_type.n])))

        NquadN = length(qp_μ4)

        # wct = m==0 ? 0.50 * ϖ .* wt_μ4  : 0.25 .* ϖ .* wt_μ4
        wct = m==0 ? 0.50 * 1 .* wt_μ4  : 0.25 .* 1 .* wt_μ4

        # Get the diagonal matrices first
        d_qp = Array(Diagonal(1 ./ qp_μ4)) 
        d_wct = Array(Diagonal(wct))

        # Calculate r⁻⁺ and t⁺⁺
        
        # Version 1: no absorption in batch mode (like before), need to separate these modes
        # if maximum(dτ) < 0.0001 
        #     r⁻⁺[:] = d_qp ⊠ Z⁻⁺ ⊠ (d_wct * dτ)
        #     t⁺⁺[:] = I_static .- (d_qp ⊠ ((I_static .- Z⁺⁺ ⊠ d_wct) * dτ))
        
        # else    
        # Version 2: with absorption in batch mode, low tau_scatt but higher tau_total, needs different equations
        # This is not yet GPU ready as it has element wise operations (should work for CPU)

            for i = 1:NquadN, j=1:NquadN

                # 𝐑⁻⁺(μᵢ, μⱼ) = ϖ ̇𝐙⁻⁺(μᵢ, μⱼ) ̇(μⱼ/(μᵢ+μⱼ)) ̇(1 - exp{-τ ̇(1/μᵢ + 1/μⱼ)}) ̇𝑤ⱼ
                r⁻⁺[i,j,:] = ϖ .* Z⁻⁺[i,j] .* (qp_μ4[j]/(qp_μ4[i]+qp_μ4[j])) .* (1 .- exp.(-dτ .* ((1/qp_μ4[i])+(1/qp_μ4[j])))) .* (wct[j]) 

                if (i==j)

                    # 𝐓⁺⁺(μᵢ, μᵢ) = ϖ ̇𝐙⁺⁺(μᵢ, μᵢ) ̇(τ/μᵢ) ̇exp{-τ/μᵢ} ̇𝑤ᵢ
                    t⁺⁺[i,j,:] = ϖ .* Z⁺⁺[i,i] .* (dτ ./ qp_μ4[i]) .* exp.(-dτ./qp_μ4[i]) .* wct[i]

                else

                    # 𝐓⁺⁺(μᵢ, μⱼ) = ϖ ̇𝐙⁺⁺(μᵢ, μⱼ) ̇(μⱼ/(μᵢ-μⱼ)) ̇(exp{-τ/μᵢ} - exp{-τ/μⱼ}) ̇𝑤ⱼ
                    # (𝑖 ≠ 𝑗)
                    t⁺⁺[i,j,:] = ϖ .* Z⁺⁺[i,j] .* (qp_μ4[j]/(qp_μ4[i]-qp_μ4[j])) .* (exp.(-dτ ./qp_μ4[i]) - exp.(-dτ./qp_μ4[j])) .* wct[j]
                end
            end

        # end

        if ndoubl<1
            r⁺⁻[:] = D ⊠ r⁻⁺ ⊠ D
            t⁻⁻[:] = D ⊠ t⁺⁺ ⊠ D
        else
            r⁻⁺[:] = D ⊠ r⁻⁺
        end
    else 
        # Note: τ is not defined here
        t⁺⁺[:] = Diagonal{exp(-τ./qp_μ4)}
        t⁻⁻[:] = Diagonal{exp(-τ./qp_μ4)}
    end    

end

function rt_elemental!(pol_type, dτ, ϖ, Z⁺⁺, Z⁻⁺, m, 
                              ndoubl, scatter, qp_μ, wt_μ, 
                              r⁻⁺::AbstractArray{FT,3}, 
                              t⁺⁺::AbstractArray{FT,3}, 
                              r⁺⁻::AbstractArray{FT,3}, 
                              t⁻⁻::AbstractArray{FT,3}, 
                              D::AbstractArray{FT,3},
                              I_static::AbstractArray) where {FT}

    rt_elemental_helper!(pol_type, dτ, ϖ, Z⁺⁺, Z⁻⁺, m, ndoubl, scatter, qp_μ, wt_μ, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻, D, I_static)
    synchronize()
end