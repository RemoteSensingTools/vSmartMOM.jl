function recursive_wigner!(m, n, l, wig_values)
    
    if ((m>=abs(n-l))&&(l>=abs(n-m))&&(n>=abs(m-l)))

        if (wig_values[m+1,n+1] === nothing) 

            if m == (l+n+1)
                wig_values[m+1,n+1] = 0.0
                return 0.0
        
            elseif (m == n+l)
                
                if (n == 1)
                    result = (-1)^l * sqrt(((l+1) * (l+2))/((2l+1) * (2l+2) * (2l+3)))
                    wig_values[m+1,n+1] = result
                    return result
                else 
                    result = -1 * recursive_wigner!(l+n-1, n-1, l, wig_values) * sqrt((n * (2n-1) * ((n+l)^2-1))/ ((n+l) * (2(n+l)+1) * (n^2-1)))
                    wig_values[m+1,n+1] = result
                    return result
                end
        
            else 
        
                k = m + 1
                D_k_nl = (1/k) * sqrt((k^2 - 1) * (k^2 - (l-n)^2) * ((n+l+1)^2 - k^2)) # Checked, correct
                M_k = 1 - (n * (n+1) - l * (l+1))/(k * (k+1)) # Checked, correct
                kP1 = k+1
                D_kP1_nl = (1/kP1) * sqrt((kP1^2 - 1) * (kP1^2 - (l-n)^2) * ((n+l+1)^2 - kP1^2))

                # The correct term is (2m + 3)
                # Not sqrt(2m+3)   
                # !!! 
                result = (1/D_k_nl) * (M_k * (2*m+3) * recursive_wigner!(m+1, n, l, wig_values) - D_kP1_nl*recursive_wigner!(m+2, n, l, wig_values))

                wig_values[m+1,n+1] = result
                return result
            end

        else
            return wig_values[m+1,n+1]
        end
    else 
        return 0.0
    end
end

m = 5
n = 3
l = 3
arr_to_fill = Array{Any}(nothing, l+n+1, n+1);

recursive_wigner!(m, n, l, arr_to_fill)

using WignerSymbols
Float64(wigner3j(m, n, l, -1, 1, 0))

