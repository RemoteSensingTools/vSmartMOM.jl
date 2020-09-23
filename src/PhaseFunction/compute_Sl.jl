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

function compute_avg_anbn!(an, bn,  mat_anam, mat_bnbm, mat_anbm, mat_bnam,w, Nmax, N_max_)
    FT2 = eltype(an)
    fill!(mat_anam,0)
    fill!(mat_bnbm,0)
    fill!(mat_anbm,0)
    fill!(mat_bnam,0)
    @inbounds for n=1:Nmax
        @inbounds for m=n:Nmax
            anam = FT2(0);
            bnbm = FT2(0);
            anbm = FT2(0);
            bnam = FT2(0);
            @inbounds for i = 1:size(an,1)
                if m < N_max_[i] && n < N_max_[i]
                    anam += w[i] * an[i,n]' * an[i,m]
                    bnbm += w[i] * bn[i,n]' * bn[i,m]
                    anbm += w[i] * an[i,n]' * bn[i,m]
                    bnam += w[i] * bn[i,n]' * an[i,m]
                end
            end 
            @inbounds mat_anam[m,n] = anam;
            @inbounds mat_bnbm[m,n] = bnbm;
            @inbounds mat_anbm[m,n] = anbm;
            @inbounds mat_bnam[m,n] = bnam;
        end
    end
    return nothing
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
    return nothing
end

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
                # first_term += avg * (2n+1) * (2m+1) * get(wigner_A, (m, n, l), 0.0) * get(wigner_B, (m, n, l), 0.0)
            end

            # anan,anbn,bnan,bnbn = (mat_anam[n,n], mat_bnbm[n,n], mat_anbm[n,n], mat_bnam[n,n])
            anan,anbn,bnan,bnbn = (mat_anam[n,n], mat_anbm[n,n], mat_bnam[n,n], mat_bnbm[n,n])
            avg = anan - anbn + bnan - bnbn
            second_term += avg * 2 * (2n+1) * (2n+1) * wigner_A[n, n, l] * wigner_B[n, n, l]
            # second_term += avg * (2n+1) * (2n+1) * get(wigner_A, (n, n, l), 0.0) * get(wigner_B, (n, n, l), 0.0)
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
                # first_term += (real_avg * (get(wig_lnm, (m, n, l), 0.0))^2 * coef_lnm2)
            end

            coef_lnn2 = (2n + 1)^2 * (ν₂_positive_flag ? 1 : (exp_m1(ll)))
            
            if (ν₂_positive_flag)
                second_term += an_p_bn[n] * (wig_lnm[n, n, l])^2 * coef_lnn2
                # second_term += an_p_bn[n] * (get(wig_lnm, (n, n, l), 0.0))^2 * coef_lnn2
            else 
                second_term += an_m_bn[n] * (wig_lnm[n, n, l])^2 * coef_lnn2
                # second_term += an_m_bn[n] * (get(wig_lnm, (n, n, l), 0.0))^2 * coef_lnn2
            end

        end
    end

    return coef * (first_term + second_term)
end