
"Prototype doubling methods, compute homogenous layer matrices from its elemental layer in `ndoubl` doubling steps"
function rt_doubling!(mo::Stokes_IQUV, dτ, τ_total, m, ndoubl, r_elt_mp, t_elt_pp, r_elt_pm, t_elt_mm)
    # # ToDo: Important output doubling applied to elemental layer, using same variables r_elt_pm, r_elt_mp, t_elt_mm, t_elt_pp (can be renamed to t⁺⁺, etc)
    # Need to check with paper nomenclature. This is basically eqs. 23-28 in vSmartMOM but using simplifications in eq. 29-32)
    τ_total=dτ
    if (ndoubl==0)
        return
    end

    for n = 1:ndoubl
        M1=int(I-r_elt_mp*r_elt_mp)
        tr_elt_mp = r_elt_mp + t_elt_pp*r_elt_mp*M1*t_elt_pp
        tt_elt_pp = t_elt_pp*M1*t_elt_pp

        r_elt_mp = tr_elt_mp
        t_elt_pp = tt_elt_pp

        τ_total = 2 * τ_total
    end
    for iμ = 1:Nquad4; jμ = 1:Nquad4
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

"minimum number of doublings needed to reach an optical depth τ_end, starting with an optical depth dτ.
#The starting optical depth dτ is also determined from its maximum possible value, τ"
function doubling_number(dτ_max, τ_end) #check if τ_end can be replaced by τ_end*ϖ for absorbing atmospheres
    #minimum number of doublings needed to reach an optical depth τ_end, starting with an optical depth dτ.
    #The starting optical depth dτ is also determined from its maximum possible value, dτ_max
    if τ_end<=dτ_max
        dτ = τ_end
        ndoubl = 0
        return dτ, ndoubl
    else
        q1 = log10(2.0)
        q2 = log10(dτ_max)
        q3 = log10(τ_end)
        tlimit = (q3 - q2)/q1
        nlimit = floor(Int,tlimit)
        diff = tlimit - nlimit
        if diff < 1.e-16
            dτ = dτ_max
            ndoubl = nlimit
        else
            ndoubl = nlimit+1       
            x = q3 - q1 * ndoubl
            dτ = 10.0^x
        end 
        return dτ, ndoubl
    end
end