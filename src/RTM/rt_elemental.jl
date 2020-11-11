"Elemental single-scattering layer"
function rt_elemental!(mo::Stokes_IQUV, dτ, ϖ, bb, m, n, ndoubl, scatter)
    # ToDo: Main output is r_elt_pm, r_elt_mp, t_elt_mm, t_elt_pp (can be renamed to t⁺⁺, etc)
    # Need to check with paper nomenclature. This is basically eqs. 19-20 in vSmartMOM

    #dτ: optical depth of elemental layer
    #ϖ: single scattering albedo of elemental layer
    #bb: thermal source function at the upper boundary of the elemental layer
    #m: fourier moment
    #n: layer of which this is an elemental
    #ndoubl: number of doubling computations needed to progress from the elemental layer to the full homogeneous layer n
    #scatter: flag indicating scattering

    r_elt_mp = zeros(Nquad4,Nquad4)
    t_elt_pp = zeros(Nquad4,Nquad4)
    if scatter
        #TODO: import vector containing quadrature cosines qp_μ of length Nquad4
        #TODO: import vector containing quadrature weights wt_μ of length Nquad4
        #TODO: construct composite, post-truncation dτ=τ/2^{ndoubl} , ϖ, p⁺⁺, p⁻⁺ matrices and import them here
        if m==0
            wct=0.50 * ϖ * wt_μ 
        else    
            wct=0.25 * ϖ * wt_μ
        end
        r_elt_mp = diag(1./qp_μ) * p⁻⁺ * diag(wct) * dτ
        t_elt_pp = I - (diag(1./qp_μ) * (I - p⁺⁺ * diag(wct)) * dτ)
        if ndoubl<1
            for iμ in 1:Nquad4; jμ in 1:Nquad4
                # That "4" and Nquad4 needs to be dynamic, coming from the PolType struct.
                i=mod(iμ-1,4)
                j=mod(jμ-1,4)
                if ((i<=1)&(j<=1)) | ((i>=2)&(j>=2))
                    r_elt_pm[iμ,jμ] = r_elt_mp[iμ,jμ]
                    t_elt_mm[iμ,jμ] = t_elt_pp[iμ,jμ]
                else
                    r_elt_pm[iμ,jμ] = - r_elt_mp[iμ,jμ]
                    t_elt_mm[iμ,jμ] = - t_elt_pp[iμ,jμ]
                end
            end 
        else
            #For doubling, transform R->DR, where D = diag{1,1,-1,-1}
            for iμ in 1:Nquad4;
                i=mod(iμ-1,4)    
                if (i>=2)
                    r_elt_pm[iμ,:] = - r_elt_mp[iμ,:]
                end
            end 
        end
    else
        t_elt_pp = diag{exp(τ./qp_μ)}
        t_elt_mm = diag{exp(τ./qp_μ)}
    end 
    return nothing
end