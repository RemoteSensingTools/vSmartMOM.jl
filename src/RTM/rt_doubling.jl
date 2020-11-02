"Prototype doubling methods, compute homogenous layer matrices from its elemental layer in `ndoubl` doubling steps"
function rt_doubling!(mo::Stokes_IQUV, dτ, τ_total, m, ndoubl)
    # # ToDo: Important output doubling applied to elemental layer, using same variables r_elt_pm, r_elt_mp, t_elt_mm, t_elt_pp (can be renamed to t⁺⁺, etc)
    # Need to check with paper nomenclature. This is basically eqs. 23-28 in vSmartMOM but using simplifications in eq. 29-32)
    τ_total=dτ
    if (ndoubl==0)
        return
    end

    for n in 1:ndoubl
        M1=int(I-r_elt_mp*r_elt_mp)
        tr_elt_mp = r_elt_mp + t_elt_pp*r_elt_mp*M1*t_elt_pp
        tt_elt_pp = t_elt_pp*M1*t_elt_pp

        r_elt_mp = tr_elt_mp
        t_elt_pp = tt_elt_pp

        τ_total = 2 * τ_total
    end
    for iμ in 1:Nquad4; jμ in 1:Nquad4
        # That "4" and Nquad4 needs to be dynamic, coming from the PolType struct.
        i=mod(iμ-1,4)
        j=mod(jμ-1,4)
        if (i<=2)
            r_elt_mp[iμ,jμ] = - r_elt_mp[iμ, jμ]
        end
        if ((i<=1)&(j<=1)) | ((i>=2)&(j>=2))
            r_elt_pm[iμ,jμ] = r_elt_mp[iμ,jμ]
            t_elt_mm[iμ,jμ] = t_elt_pp[iμ,jμ]
        else
            r_elt_pm[iμ,jμ] = - r_elt_mp[iμ,jμ]
            t_elt_mm[iμ,jμ] = - t_elt_pp[iμ,jμ]
        end
    end 
    return nothing
end

"minimum number of doublings needed to reach an optical depth τ_end, starting with an optical depth τ_min.
#The starting optical depth τ_min is also determined from its maximum possible value, τ"
function doubling_number!(mo::Stokes_IQUV, τ, τ_end, τ_min, ndoubl)
    #minimum number of doublings needed to reach an optical depth τ_end, starting with an optical depth τ_min.
    #The starting optical depth τ_min is also determined from its maximum possible value, τ
    if τ_end<=τ
        τ_min = τ_end
        ndoubl = 0
        return nothing
    else
        q1 = log10(2.0)
        q2 = log10(τ)
        q3 = log10(τ_end)
        tlimit = (q3 - q2)/q1
        nlimit = floor(tlimit)
        diff = tlimit - nlimit
        if diff < 1.e-16
            τ_min = τ
            ndoubl = nlimit
        else
            ndoubl = nlimit+1       
            x = q3 - q1 * ndoubl
            τ_min = 10.0^x
        end 
        return nothing
    end
end