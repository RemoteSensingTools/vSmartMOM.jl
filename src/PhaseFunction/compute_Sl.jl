exp_m1(x) = iseven(x) ? 1 : -1

function get_abnabm(an,bn,n,m,w)
    FT = eltype(an)
    anam = bnbm = anbm = bnam = FT(0);
    @inbounds for i=1:size(an)[2]
        anam += w[i] * an[n,i]' * an[m,i]
        bnbm += w[i] * bn[n,i]' * bn[m,i]
        anbm += w[i] * an[n,i]' * bn[m,i]
        bnam += w[i] * bn[n,i]' * an[m,i]
    end
    # println(size(w))
    # println(size(an[n,:]'))
    # println(size(an[m,:]))
    # @inbounds anam = (w .* an[n,:])' * an[m,:]
    # @inbounds bnbm = (w .* bn[n,:])' * bn[m,:]
    # @inbounds anbm = (w .* an[n,:])' * bn[m,:]
    # @inbounds bnam = (w .* bn[n,:])' * an[m,:]
    
    return anam,bnbm,anbm,bnam
end

function compute_Sl(l, ν₁, ν₂, ν₂_positive_flag, k, N_max, an, bn, wigner_A, wigner_B, w)

    @assert ((ν₁, ν₂) == (2, -2) || (ν₁, ν₂) == (2, 2) || (ν₁, ν₂) == (0, 0) || (ν₁, ν₂) == (0, 2)) "Not a valid input to calculate Sl"
    @assert ((ν₂ == 2 && ν₂_positive_flag == true) || (ν₂ == -2 && ν₂_positive_flag == false) || ν₂ == 0) "Sign flag doesn't match ν₂"

    ll = l-1
    coef = (2ll+1)*π / k^2

    first_term = 0
    second_term = 0

    if (ν₁, ν₂) == (0, 2)

        for n = 1:N_max
            m_star = max(ll-n,n+1)
            m_max = min(ll+n,N_max)

            for m = m_star:m_max

                anam,bnbm,anbm,bnam = get_abnabm(an,bn,n,m,w)
                aman,bmbn,ambn,bman = get_abnabm(an,bn,m,n,w)

                avg = (exp_m1(ll + n + m)) * (anam + bnam - anbm - bnbm) + (aman - ambn + bman - bmbn)
                first_term += avg * (2n+1) * (2m+1) * wigner_A[m, n, l] * wigner_B[m, n, l]
            end

            anan,anbn,bnan,bnbn = get_abnabm(an,bn,n,n,w)
            avg = anan - anbn + bnan - bnbn
            second_term += avg * (2n+1) * (2n+1) * wigner_A[n, n, l] * wigner_B[n, n, l]
        end

    else
        wig_lnm = ν₁ == 0 ? wigner_A : wigner_B

        for n = 1:N_max
            m_star = max(ll-n,n+1)
            m_max = min(ll+n,N_max)
            
            for m = m_star:m_max
                anam,bnbm,anbm,bnam = get_abnabm(an,bn,n,m,w)
                real_avg = ν₂_positive_flag ? real(anam+anbm+bnam+bnbm) : real(anam-anbm-bnam+bnbm);
                coef_lnm2 = 2 * (2m + 1)  * (2n + 1) * (ν₂_positive_flag ? 1 : (exp_m1(ll + n + m)))
                first_term += (real_avg * (wig_lnm[m, n, l])^2 * coef_lnm2)
            end

            coef_lnn2 = (2n + 1)^2 * (ν₂_positive_flag ? 1 : (exp_m1(ll)))

            # for i = 1:size(an,2)
            #     if (ν₂_positive_flag)
            #         second_term += w[i]' * (abs2(an[n,i] + bn[n,i])) * (wig_lnm[n, n, l])^2 * coef_lnn2
            #     else 
            #         second_term += w[i]' * (abs2(an[n,i] - bn[n,i])) * (wig_lnm[n, n, l])^2 * coef_lnn2
            #     end
            # end

            if (ν₂_positive_flag)
                second_term += w' *(abs2.(an[n,:] .+ bn[n,:])) * (wig_lnm[n, n, l])^2 * coef_lnn2
            else 
                second_term += w' *(abs2.(an[n,:] .- bn[n,:])) * (wig_lnm[n, n, l])^2 * coef_lnn2
            end
        end
    end

    return coef * (first_term + second_term)
end