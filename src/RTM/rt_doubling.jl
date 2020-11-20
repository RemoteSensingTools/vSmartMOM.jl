
"Prototype doubling methods, compute homogenous layer matrices from its elemental layer in `ndoubl` doubling steps"
function rt_doubling(dτ, τ_total, ndoubl, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻)
    # # ToDo: Important output doubling applied to elemental layer, using same variables r⁺⁻, r⁻⁺, t⁻⁻, t⁺⁺ (can be renamed to t⁺⁺, etc)
    # Need to check with paper nomenclature. This is basically eqs. 23-28 in vSmartMOM but using simplifications in eq. 29-32)
    Nquad4 = size(r⁻⁺, 1)
    if (ndoubl==0)
        @assert (τ_total==dτ*2^ndoubl)
        return r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻
    end
    τ_total=dτ
    for n = 1:ndoubl
        M1=inv(I - r⁻⁺ * r⁻⁺)
        tr⁻⁺ = r⁻⁺ + t⁺⁺ * r⁻⁺ * M1 * t⁺⁺
        tt⁺⁺ = t⁺⁺ * M1 * t⁺⁺

        r⁻⁺ = tr⁻⁺
        t⁺⁺ = tt⁺⁺

        τ_total = 2 * τ_total
    end
    for iμ = 1:Nquad4, jμ = 1:Nquad4
        # That "4" and Nquad4 needs to be dynamic, coming from the PolType struct.
        i=mod(iμ-1,4)
        j=mod(jμ-1,4)
        #@show i,j
        if (i>=2)
            r⁻⁺[iμ,jμ] = - r⁻⁺[iμ, jμ]
        end
        if ((i<=1)&(j<=1)) | ((i>=2)&(j>=2))
            r⁺⁻[iμ,jμ] = r⁻⁺[iμ,jμ]
            t⁻⁻[iμ,jμ] = t⁺⁺[iμ,jμ]
        else
            r⁺⁻[iμ,jμ] = - r⁻⁺[iμ,jμ]
            t⁻⁻[iμ,jμ] = - t⁺⁺[iμ,jμ]
        end
    end 
    @assert (τ_total==dτ*2^ndoubl)
    return r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻
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