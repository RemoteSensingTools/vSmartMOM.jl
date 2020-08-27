"""
    $(FUNCTIONNAME)(m, n, l, wig_values)

Compute the Wigner 3j values for all n to nmax and all m from mstart to (n+l+1)

"""
function recursive_wigner_m110!(m::Integer, n::Integer, l::Integer, wig_values)
    
    # Only compute the value if this condition is met
    if ((m>=abs(n-l))&&(l>=abs(n-m))&&(n>=abs(m-l)))

        # Only perform this computation if it has not been calculated yet. 
        # In these cases, always store the result before returning
        # This technique is called ** memoization ** (forgot the word til just now)
        if (wig_values[m+1,n+1] === nothing) 

            # Base Case 1: m is l+n+1, in which case the result is 0
            if m == (l+n+1)

                wig_values[m+1,n+1] = 0.0
                return 0.0
        
            # Base Case 2: m is l+n
            elseif (m == n+l)
                
                # If n is 1, return fixed expression (Eqn. 28)
                if (n == 1)
                    result = (-1)^l * sqrt(((l+1) * (l+2))/((2l+1) * (2l+2) * (2l+3)))
                    wig_values[m+1,n+1] = result
                    return result

                # If n is not 1, return recursive expression (recursing on n, Eqn. 27)
                else 
                    result = -recursive_wigner_m110!(l+n-1, n-1, l, wig_values) * 
                              sqrt((n * (2n-1) * ((n+l)^2-1))/ ((n+l) * (2(n+l)+1) * (n^2-1)))
                    wig_values[m+1,n+1] = result
                    return result
                end
        
            # Main Recursive Case
            else 
        
                # Define variables (Eqn. 26)
                k = m + 1
                D_k_nl = (1/k) * sqrt((k^2 - 1) * (k^2 - (l-n)^2) * ((n+l+1)^2 - k^2)) 
                M_k = 1 - (n * (n+1) - l * (l+1))/(k * (k+1)) 
                kP1 = k+1
                D_kP1_nl = (1/kP1) * sqrt((kP1^2 - 1) * (kP1^2 - (l-n)^2) * ((n+l+1)^2 - kP1^2))

                # Perform recursion (Eqn. 25)
                result = (1/D_k_nl) * 
                         (M_k * (2*m+3) * recursive_wigner_m110!(m+1, n, l, wig_values) 
                             - D_kP1_nl * recursive_wigner_m110!(m+2, n, l, wig_values))

                wig_values[m+1,n+1] = result
                return result
            end
        
        # This value has already been calculated, so just return it
        else
            return wig_values[m+1,n+1]
        end

    # Condition is not met: return 0
    else 
        return 0.0
    end
end

