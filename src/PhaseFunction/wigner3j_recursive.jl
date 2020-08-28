###
### Compute Wigner 3-j symbols recursively, as specified in paper: 
###     "Revisiting the Fourier expansion of Mie scattering matrices 
###      in generalized spherical functions", Suniti Sanghavi, 2014
###     
### Note: Minor typos from the paper are corrected in this implementation: 
###     1. Eqn 25:
###       -   sqrt(2*m+3) 
###       +   (2*m+3)
### 
###     2. Eqn 31: 
###       -   sqrt((l-1) * l * (l+1) * (l+2))
###       +   ((l-1) * l * (l+1) * (l+2)) ^ (-0.5)
###

# Compute wigner 3j symbol where m₁ = -1, m₂ == 1, and m₃ == 0
function wigner_m110!(m::Integer, n::Integer, l::Integer, wig_values)
    
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
                    result = -wigner_m110!(l+n-1, n-1, l, wig_values) * 
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
                # NOTE: Typo found in paper's equation. Corrected here. 
                # Incorrect: sqrt(2*m+3)
                # Correct: 2*m+3
                result = (1/D_k_nl) * 
                         (M_k * (2*m+3) * wigner_m110!(m+1, n, l, wig_values) 
                             - D_kP1_nl * wigner_m110!(m+2, n, l, wig_values))

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

# Compute wigner 3j symbol where m₁ = 0, m₂ == 0, and m₃ == 0
function wigner_000!(m::Integer, n::Integer, l::Integer, wig_values)

    # If m = n+l+1, return 0 (does not satisfy condition)
    if (m == n + l + 1)
        return 0.0
    # If m = n+l, return the result specified (Eqn. 30)
    elseif (m == n + l)
        return wigner_m110!(n+l, n, l, wig_values) * 2 * 
               sqrt((n+l) * (n+l+1) * n * (n+1))/(l*(l+1) - (l+n)*(l+n+1) - n*(n+1))
    # If m < n+l,
    else 

        # If m+n+l is even, return the result specified (Eqn. 29)
        if (iseven(m + n + l == 0))
            return -wigner_000!(m+2, n, l, wig_values) * 
                    sqrt(((m+2)^2 - (n-l)^2 ) / ((m+1)^2 - (n-l)^2 )) * 
                    sqrt((1-(1/(n+l-m))) * (1+(1/(m+n+l+2))))
        # If m+n+l is odd, return 0 
        else
            return 0.0
        end

    end
    
end

# Compute wigner 3j symbol where m₁ = -1, m₂ == -1, and m₃ == 2
function wigner_m1m12!(m::Integer, n::Integer, l::Integer, wig_values)

    # Simply return the result specified in Eqn. 31 
    return (-1)^(m+n+l) * (((l-1) * l * (l+1) * (l+2))^(-0.5)) * 
           ( (m*(m+1) + (-1)^(m+n+l)*n*(n+1)) * 
             wigner_m110!(m, n, l, wig_values) +
             2 * sqrt(m*(m+1)*n*(n+1)) * 
             wigner_000!(m, n, l, wig_values) )
end

"""
    $(FUNCTIONNAME)(m, n, l, wig_values)

Compute the Wigner 3j values for all n to nmax and all m from mstart to (n+l+1)

"""
function wigner!(j₁::Integer, j₂::Integer, j₃::Integer, 
                 m₁::Integer, m₂::Integer, m₃::Integer)

    # Switch between the three different configurations and call appropriate function
    # (m₁, m₂, m₃): (-1, 1, 0), (-1, -1, 2), (0, 0, 0)
    if (m₁ == -1 && m₂ == 1 && m₃ == 0)

        m110_values = Array{Any}(nothing, j₂+j₃+1, j₂+1);
        return wigner_m110!(j₁, j₂, j₃, m110_values)

    elseif (m₁ == -1 && m₂ == -1 && m₃ == 2)

        m110_values = Array{Any}(nothing, j₂+j₃+1, j₂+1);
        m1m12_values = Array{Any}(nothing, j₂+j₃+1, j₂+1);
        return wigner_m1m12!(j₁, j₂, j₃, m110_values)

    elseif (m₁ == 0 && m₂ == 0 && m₃ == 0)

        m110_values = Array{Any}(nothing, j₂+j₃+1, j₂+1);
        return wigner_000!(j₁, j₂, j₃, m110_values)

    else
        throw(DomainError((m₁, m₂, m₃,), 
            "Unsupported configuration. Supported: (-1,1,0), (-1,-1,2), (0,0,0)"))
    end

end