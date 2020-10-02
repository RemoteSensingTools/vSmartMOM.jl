function compute_Sl(l, ν₁, ν₂, ν₂_positive_flag, k, N_max, an, bn, mat_anam, mat_bnbm, mat_anbm, mat_bnam, wigner_A, wigner_B,an_m_bn,an_p_bn, w)

    @assert ((ν₁, ν₂) == (2, -2) || (ν₁, ν₂) == (2, 2) || (ν₁, ν₂) == (0, 0) || (ν₁, ν₂) == (0, 2)) "Not a valid input to calculate Sl"
    @assert ((ν₂ == 2 && ν₂_positive_flag == true) || (ν₂ == -2 && ν₂_positive_flag == false) || ν₂ == 0) "Sign flag doesn't match ν₂"

    ll = l-1
    coef = (2ll+1)*π / k^2

    first_term = 0.0
    second_term = 0.0

    if (ν₁, ν₂) == (0, 2)

        @inbounds for n = 1:N_max
            m_star = max(ll-n,n+1)
            m_max = min(ll+n,N_max)

            @inbounds for m = m_star:m_max
                avg = (exp_m1(ll + n + m)) * (mat_anam[m,n] + mat_bnam[m,n] - mat_anbm[m,n] - mat_bnbm[m,n]) + (mat_anam[m,n]' - mat_bnam[m,n]' + mat_anbm[m,n]' - mat_bnbm[m,n]')
                first_term += avg * 2 * (2n+1) * (2m+1) * wigner_A[m, n, l] * wigner_B[m, n, l]
            end

            anan,anbn,bnan,bnbn = (mat_anam[n,n], mat_anbm[n,n], mat_bnam[n,n], mat_bnbm[n,n])
            avg = anan - anbn + bnan - bnbn
            second_term += avg * 2 * (2n+1) * (2n+1) * wigner_A[n, n, l] * wigner_B[n, n, l]
        end

    else
        wig_lnm = ν₁ == 0 ? wigner_A : wigner_B

        @inbounds for n = 1:N_max
            m_star = max(ll-n,n+1)
            m_max = min(ll+n,N_max)
            
            @inbounds for m = m_star:m_max
                anam,bnbm,anbm,bnam = (mat_anam[m,n], mat_bnbm[m,n], mat_anbm[m,n], mat_bnam[m,n])
                real_avg = ν₂_positive_flag ? real(anam+anbm+bnam+bnbm) : real(anam-anbm-bnam+bnbm);
                coef_lnm2 = 2 * (2m + 1)  * (2n + 1) * (ν₂_positive_flag ? 1 : (exp_m1(ll + n + m)))
                first_term += (real_avg * (wig_lnm[m, n, l])^2 * coef_lnm2)
            end

            coef_lnn2 = (2n + 1)^2 * (ν₂_positive_flag ? 1 : (exp_m1(ll)))
            
            if (ν₂_positive_flag)
                second_term += an_p_bn[n] * (wig_lnm[n, n, l])^2 * coef_lnn2
            else 
                second_term += an_m_bn[n] * (wig_lnm[n, n, l])^2 * coef_lnn2
            end

        end
    end

    return coef * (first_term + second_term)
end