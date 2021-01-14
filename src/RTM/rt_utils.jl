function get_kn(kn, scatter, iz)

    # First layer (TOA)
    if (iz == 1)

        # If scattering, 4. If non-scattering, 1. 
        kn = scatter ? 4 : 1
    
    # Not the first layer (not TOA)
    elseif (kn >= 1)

        # If kn was 1, then toggle between 0/0 and 0/1 
        # Else, toggle between 1/0 and 1/1
        kn = (kn == 1) ? (!scatter ? 1 : 2) : (!scatter ? 3 : 4)
    end

    return kn
end

"minimum number of doublings needed to reach an optical depth τ_end, starting with an optical depth dτ.
#The starting optical depth dτ is also determined from its maximum possible value, τ"
function doubling_number(dτ_max, τ_end) # check if τ_end can be replaced by τ_end*ϖ for absorbing atmospheres
    FT = eltype(dτ_max)
    # @show FT, eltype(τ_end), dτ_max, τ_end
    # minimum number of doublings needed to reach an optical depth τ_end, starting with an optical depth dτ.
    # The starting optical depth dτ is also determined from its maximum possible value, dτ_max
    if τ_end <= dτ_max
        dτ = τ_end
        ndoubl = 0
        return dτ, ndoubl
    else
        q1 = log10(2.0)
        q2 = log10(dτ_max)
        q3 = log10(τ_end)
        tlimit = (q3 - q2) / q1
        nlimit = floor(Int, tlimit)
        diff = tlimit - nlimit
        if diff < eps(FT)
            dτ = dτ_max
            ndoubl = nlimit
        else
            ndoubl = nlimit + 1       
            x = q3 - q1 * ndoubl
            dτ = 10.0^x
        end 
        return dτ, ndoubl
    end
end

# Finds index i of f_array (i) which is nearest point to f
function nearest_point(f_array, f)
    d0 = 999.9
    index = 0
    for i in 1:length(f_array)
        d = abs(f_array[i] - f)
        if d < d0
            d0 = d
            index = i
        end     
    end
    return index
end