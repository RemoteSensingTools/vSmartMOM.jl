exp_m1(x) = iseven(x) ? 1 : -1

function get_abnabm(an,bn,n,m,w)
    FT = eltype(an)
    anam = bnbm = anbm = bnam = FT(0);
    @inbounds for i=1:size(an)[2]
        anam += w[i] * an[i,n]' * an[i,m]
        bnbm += w[i] * bn[i,n]' * bn[i,m]
        anbm += w[i] * an[i,n]' * bn[i,m]
        bnam += w[i] * bn[i,n]' * an[i,m]
    end
    return anam,bnbm,anbm,bnam
end

# This can add another flag/number to avoid multiplications by 0 (e.g. where an,bn is 0)
@kernel function avg_anbn!(@Const(an), @Const(bn), mat_anam, mat_bnbm, mat_anbm, mat_bnam, @Const(w), @Const(nMax))
    FT = eltype(an)
    # Indices over n and m
    m, n, i = @index(Global, NTuple)
    if m>=n && m<nMax[i] && n < nMax[i]
        @inbounds mat_anam[m,n] += (w[i] * (an[i,n]' * an[i,m]));
        @inbounds mat_bnbm[m,n] += (w[i] * (bn[i,n]' * bn[i,m]));
        @inbounds mat_anbm[m,n] += (w[i] * (an[i,n]' * bn[i,m]));
        @inbounds mat_bnam[m,n] += (w[i] * (bn[i,n]' * an[i,m]));
 
    end
end

function fill_avg_anbns!(an, bn, mat_anam, mat_bnbm, mat_anbm, mat_bnam, wₓ, N_max, N_max_, architecture)

    # Fill all matrices with 0
    [fill!(mat,0) for mat in [mat_anam, mat_bnbm, mat_anbm, mat_bnam]]

    # Set the kernel device
    kernel! = avg_anbn!(architecture)

    # Let it run
    event = kernel!(an, bn, mat_anam, mat_bnbm, mat_anbm, mat_bnam, wₓ, N_max_, ndrange = (N_max, N_max, length(wₓ))); 
    wait(CPU(), event) 

    return nothing  
end

function goCUDA!(mat_anamC, mat_bnbmC,mat_anbmC,mat_bnamC)
    fill!(mat_anamC,0);
    fill!(mat_bnbmC,0);
    fill!(mat_bnamC,0);
    fill!(mat_anbmC,0);
    kernel! = avg_anbn!(CUDADevice())
    event = kernel!(anC,bnC,mat_anamC, mat_bnbmC,mat_anbmC,mat_bnamC,wₓC,N_max_C, ndrange=(Nmax,Nmax,length(wₓ))); 
    wait(CUDADevice(), event)    ;
    #@show  Array(mat_anamC)[12,10]  
    return nothing
end

function compute_Sl(l, ν₁, ν₂, ν₂_positive_flag, k, N_max, an, bn, mat_anam, mat_bnbm, mat_anbm, mat_bnam, wigner_A, wigner_B, w)

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

                anam,bnbm,anbm,bnam = (mat_anam[m,n], mat_bnbm[m,n], mat_anbm[m,n], mat_bnam[m,n])
                aman,bmbn,ambn,bman = get_abnabm(an,bn,m,n,w)

                avg = (exp_m1(ll + n + m)) * (anam + bnam - anbm - bnbm) + (aman - ambn + bman - bmbn)
                first_term += avg * (2n+1) * (2m+1) * wigner_A[m, n, l] * wigner_B[m, n, l]
            end

            anan,anbn,bnan,bnbn = (mat_anam[n,n], mat_bnbm[n,n], mat_anbm[n,n], mat_bnam[n,n])
            avg = anan - anbn + bnan - bnbn
            second_term += avg * (2n+1) * (2n+1) * wigner_A[n, n, l] * wigner_B[n, n, l]
        end

    else
        wig_lnm = ν₁ == 0 ? wigner_A : wigner_B

        for n = 1:N_max
            m_star = max(ll-n,n+1)
            m_max = min(ll+n,N_max)
            
            for m = m_star:m_max
                anam,bnbm,anbm,bnam = (mat_anam[m,n], mat_bnbm[m,n], mat_anbm[m,n], mat_bnam[m,n]) # get_abnabm(an,bn,n,m,w)
                real_avg = ν₂_positive_flag ? real(anam+anbm+bnam+bnbm) : real(anam-anbm-bnam+bnbm);
                coef_lnm2 = 2 * (2m + 1)  * (2n + 1) * (ν₂_positive_flag ? 1 : (exp_m1(ll + n + m)))
                first_term += (real_avg * (wig_lnm[m, n, l])^2 * coef_lnm2)
            end

            coef_lnn2 = (2n + 1)^2 * (ν₂_positive_flag ? 1 : (exp_m1(ll)))

            for i = 1:size(an,1)
                if (ν₂_positive_flag)
                    second_term += w[i]' * (abs2(an[i,n] + bn[i,n])) * (wig_lnm[n, n, l])^2 * coef_lnn2
                else 
                    second_term += w[i]' * (abs2(an[i,n] - bn[i,n])) * (wig_lnm[n, n, l])^2 * coef_lnn2
                end
            end

            # if (ν₂_positive_flag)
            #     second_term += w' *(abs2.(an[n,:] .+ bn[n,:])) * (wig_lnm[n, n, l])^2 * coef_lnn2
            # else 
            #     second_term += w' *(abs2.(an[n,:] .- bn[n,:])) * (wig_lnm[n, n, l])^2 * coef_lnn2
            # end
        end
    end

    return coef * (first_term + second_term)
end