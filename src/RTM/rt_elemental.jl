"Elemental single-scattering layer"
function rt_elemental!(pol_type, dτ, ϖ, Z⁺⁺, Z⁻⁺, m, ndoubl, scatter, qp_μ, wt_μ, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻, D)
    # ToDo: Main output is r⁺⁻, r⁻⁺, t⁻⁻, t⁺⁺ (can be renamed to t⁺⁺, etc)
    # Need to check with paper nomenclature. This is basically eqs. 19-20 in vSmartMOM

    #dτ: optical depth of elemental layer
    #ϖ: single scattering albedo of elemental layer
    #bb: thermal source function at the upper boundary of the elemental layer
    #m: fourier moment
    #n: layer of which this is an elemental
    #ndoubl: number of doubling computations needed to progress from the elemental layer to the full homogeneous layer n
    #scatter: flag indicating scattering

    # Create temporary matrices
    I_static = one(similar(Z⁻⁺))
    aux1 = similar(Z⁻⁺)

    if scatter
        #TODO: import vector containing quadrature cosines qp_μ of length Nquad4
        #TODO: import vector containing quadrature weights wt_μ of length Nquad4
        #TODO: construct composite, post-truncation dτ=τ/2^{ndoubl} , ϖ, p⁺⁺, p⁻⁺ matrices and import them here
        
        qp_μ4 = reduce(vcat, (fill.(qp_μ,[pol_type.n])))
        wt_μ4 = reduce(vcat, (fill.(wt_μ,[pol_type.n])))
        Nquadn = length(qp_μ4)

        wct = m==0 ? 0.50 * ϖ * wt_μ4  : 0.25 * ϖ * wt_μ4


        # The following section performs: 
        # r⁻⁺[:] = Diagonal(1 ./ qp_μ4) * Z⁻⁺ * Diagonal(wct) * dτ
        # t⁺⁺[:] = I - (Diagonal(1 ./ qp_μ4) * (I - Z⁺⁺ * Diagonal(wct)) * dτ)

        # Get the diagonal matrices first
        d_qp = Diagonal(1 ./ qp_μ4)
        d_wct = Diagonal(wct)

        # Calculate r⁻⁺
        mul!(aux1, d_qp, Z⁻⁺)        # Diagonal(1 ./ qp_μ4) * Z⁻⁺
        mul!(r⁻⁺, aux1, d_wct * dτ)  # r⁻⁺ = (Diagonal(1 ./ qp_μ4) * Z⁻⁺) * Diagonal(wct) * dτ

        # Calculate t⁺⁺
        mul!(aux1, Z⁺⁺, d_wct)      # Z⁺⁺ * Diagonal(wct)
        @. aux1 = I_static - aux1   # (I - Z⁺⁺ * Diagonal(wct))
        rmul!(aux1, dτ)             # (I - Z⁺⁺ * Diagonal(wct)) * dτ
        mul!(aux1, d_qp, aux1)      # (Diagonal(1 ./ qp_μ4) * (I - Z⁺⁺ * Diagonal(wct)) * dτ)
        @. t⁺⁺ = I_static - aux1    # t⁺⁺ = I - ⤴

        #test = I - Diagonal(1 ./ qp_μ4) * dτ
        #@show t⁺⁺[1],  test[1], 1 ./ qp_μ4[1]
        
        if ndoubl<1
            mul!(aux1, D, r⁻⁺)
            mul!(r⁺⁻, aux1, D)

            mul!(aux1, D, t⁺⁺)
            mul!(t⁻⁻, aux1, D)
            #=for iμ = 1:Nquad4, jμ = 1:Nquad4
                # That "4" and Nquad4 needs to be dynamic, coming from the PolType struct.
                i=mod(iμ-1,4)
                j=mod(jμ-1,4)
                if ((i<=1)&(j<=1)) | ((i>=2)&(j>=2))
                    r⁺⁻[iμ,jμ] = r⁻⁺[iμ,jμ]
                    t⁻⁻[iμ,jμ] = t⁺⁺[iμ,jμ]
                else
                    r⁺⁻[iμ,jμ] = - r⁻⁺[iμ,jμ]
                    t⁻⁻[iμ,jμ] = - t⁺⁺[iμ,jμ]
                end
            end =#  
        else
            #For doubling, transform R->DR, where D = Diagonal{1,1,-1,-1}
            mul!(aux1, D, r⁻⁺)
            @. r⁻⁺ = aux1
            #= for iμ = 1:Nquad4;
                i=mod(iμ-1,4)    
                if (i>=2)
                    r⁻⁺[iμ,:] = - r⁻⁺[iμ,:]
                end
            end =#  
        end
    else
        t⁺⁺[:] = Diagonal{exp(-τ./qp_μ4)}
        t⁻⁻[:] = Diagonal{exp(-τ./qp_μ4)}
    end 
    return nothing 
end