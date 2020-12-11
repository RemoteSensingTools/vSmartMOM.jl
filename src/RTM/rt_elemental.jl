"Elemental single-scattering layer"
function rt_elemental!(dτ, ϖ, Z⁺⁺, Z⁻⁺, m, ndoubl, scatter, qp_μ, wt_μ, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻)
    # ToDo: Main output is r⁺⁻, r⁻⁺, t⁻⁻, t⁺⁺ (can be renamed to t⁺⁺, etc)
    # Need to check with paper nomenclature. This is basically eqs. 19-20 in vSmartMOM

    #dτ: optical depth of elemental layer
    #ϖ: single scattering albedo of elemental layer
    #bb: thermal source function at the upper boundary of the elemental layer
    #m: fourier moment
    #n: layer of which this is an elemental
    #ndoubl: number of doubling computations needed to progress from the elemental layer to the full homogeneous layer n
    #scatter: flag indicating scattering
    # dims = size(Z⁺⁺)
    # r⁻⁺ = zeros(dims)
    # t⁺⁺ = zeros(dims)
    # r⁺⁻ = zeros(dims)
    # t⁻⁻ = zeros(dims)
    if scatter
        #TODO: import vector containing quadrature cosines qp_μ of length Nquad4
        #TODO: import vector containing quadrature weights wt_μ of length Nquad4
        #TODO: construct composite, post-truncation dτ=τ/2^{ndoubl} , ϖ, p⁺⁺, p⁻⁺ matrices and import them here
        
        qp_μ4 = reduce(vcat, (fill.(qp_μ,[4])))
        wt_μ4 = reduce(vcat, (fill.(wt_μ,[4])))
        Nquad4 = length(qp_μ4)

        if m==0
            wct=0.50 * ϖ * wt_μ4 
        else    
            wct=0.25 * ϖ * wt_μ4
        end
        #@show size(Diagonal(1 ./ qp_μ4)), size(Z⁻⁺), size(Diagonal(wct) * dτ)
        r⁻⁺[:] = Diagonal(1 ./ qp_μ4) * Z⁻⁺ * Diagonal(wct) * dτ
        t⁺⁺[:] = I - (Diagonal(1 ./ qp_μ4) * (I - Z⁺⁺ * Diagonal(wct)) * dτ)
        #test = I - Diagonal(1 ./ qp_μ4) * dτ
        #@show t⁺⁺[1],  test[1], 1 ./ qp_μ4[1]
        if ndoubl<1
            for iμ = 1:Nquad4, jμ = 1:Nquad4
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
            end 
        else
            #For doubling, transform R->DR, where D = Diagonal{1,1,-1,-1}
            for iμ = 1:Nquad4;
                i=mod(iμ-1,4)    
                if (i>=2)
                    r⁻⁺[iμ,:] = - r⁻⁺[iμ,:]
                end
            end 
        end
    else
        t⁺⁺[:] = Diagonal{exp(-τ./qp_μ4)}
        t⁻⁻[:] = Diagonal{exp(-τ./qp_μ4)}
    end 
    return nothing # r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻
end