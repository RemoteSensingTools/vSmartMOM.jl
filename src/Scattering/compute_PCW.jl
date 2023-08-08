#=
 
This file specifies how to compute aerosol optical properties using the Domke-PCW method
 
=#

"""
    $(FUNCTIONNAME)(model::MieModel{FDT}) where FDT<:PCW

Reference: Suniti Sanghavi 2014, https://doi.org/10.1016/j.jqsrt.2013.12.015

Compute the aerosol optical properties using the Domke-PCW method 
Input: MieModel, holding all computation and aerosol properties 
Output: AerosolOptics, holding all Greek coefficients and Cross-Sectional information
"""
function compute_aerosol_optical_properties(model::MieModel{FDT}, FT2::Type=Float64) where FDT<:PCW

    # Unpack the model
    @unpack computation_type, aerosol, r_max, nquad_radius, λ, polarization_type, truncation_type, wigner_A, wigner_B = model

    # Extract variables from aerosol struct:
    # @unpack aerosol, r_max, nquad_radius = mie_aerosol
    @unpack size_distribution, nᵣ, nᵢ = aerosol

    # Get the refractive index's real part type
    FT = eltype(nᵣ);

    # Compute radii and weights
    # start,stop = quantile(size_distribution,[0.0025,0.9975])
    r, wᵣ = gauleg(nquad_radius, 0.0, r_max ; norm=true) 
    # r, wᵣ = gauleg(nquad_radius, start, min(stop,r_max) ; norm=true) 
    #r, wᵣ = gauleg(nquad_radius, 0.0, r_max ; norm=true)
    wₓ = compute_wₓ(size_distribution, wᵣ, r, r_max) 

    # Find overall N_max from the maximum radius
    N_max = Scattering.get_n_max(2 * π * r_max/ λ)

    # Compute an, bn values
    an, bn = compute_anbn(model, λ, r)

    # Compute average cross-sectional scattering and extinction
    k = 2π / λ
    avg_C_scatt, avg_C_ext = compute_avg_C_scatt_ext(k, an, bn, wₓ)

    # l range to go over
    ls = 1:(2 * N_max - 1)

    # Where to store the Greek Coefficients
    greek_coefs = zeros(FT, 6, size(ls, 1))

    # Pre-compute anbn averages
    # That is, pre-compute (an✶)am, (an✶)bm, (bn✶)am, (bn✶)bm 
    # So that you can quickly compute (an✶ + bn✶) × (am + bm)
    ComplexFT = Complex{FT}
    mat_anam = LowerTriangular(zeros(ComplexFT, N_max, N_max));
    mat_anbm = LowerTriangular(zeros(ComplexFT, N_max, N_max));
    mat_bnam = LowerTriangular(zeros(ComplexFT, N_max, N_max));
    mat_bnbm = LowerTriangular(zeros(ComplexFT, N_max, N_max));
    ab_pairs = (mat_anam, mat_anbm, mat_bnam, mat_bnbm)

    # Get the N_max value corresponding to each radius 
    N_max_ = Scattering.get_n_max.(2π * r/ λ)

    # Then pre-compute anbn averages
    Scattering.compute_avg_anbns!(an, bn, ab_pairs, wₓ, N_max, N_max_)

    # Pre-compute |an ± bn|² * w
    an_m_bn = transpose(abs2.(an-bn)) * wₓ
    an_p_bn = transpose(abs2.(an+bn)) * wₓ

    # For each l
    @showprogress 1 "Computing S functions ..." for l in ls

        # Compute all necessary Sl_νν values 
        # Eq 22 in Sanghavi 2014

        Sl_00  = compute_Sl(l, 0, 0,  true,  k, N_max, ab_pairs, an_m_bn, an_p_bn, wigner_A, wigner_B)
        Sl_0m0 = compute_Sl(l, 0, 0,  false, k, N_max, ab_pairs, an_m_bn, an_p_bn, wigner_A, wigner_B)
        Sl_22  = compute_Sl(l, 2, 2,  true,  k, N_max, ab_pairs, an_m_bn, an_p_bn, wigner_A, wigner_B)
        Sl_2m2 = compute_Sl(l, 2, -2, false, k, N_max, ab_pairs, an_m_bn, an_p_bn, wigner_A, wigner_B)
        Sl_02  = compute_Sl(l, 0, 2,  true,  k, N_max, ab_pairs, an_m_bn, an_p_bn, wigner_A, wigner_B)

        # Use Sl_νν values to compute Greek coefficients 
        # Eq 24 in Sanghavi 2014

        @inbounds greek_coefs[1,l] = (1/avg_C_scatt) * (Sl_00 + Sl_0m0)
        @inbounds greek_coefs[2,l] = (1/avg_C_scatt) * (Sl_00 - Sl_0m0)
        @inbounds greek_coefs[3,l] = (1/avg_C_scatt) * (Sl_22 + Sl_2m2)
        @inbounds greek_coefs[4,l] = (1/avg_C_scatt) * (Sl_22 - Sl_2m2)
        @inbounds greek_coefs[5,l] = (1/avg_C_scatt) * real(Sl_02)
        @inbounds greek_coefs[6,l] = (1/avg_C_scatt) * imag(Sl_02)
    end

    # Create GreekCoefs object with α, β, γ, δ, ϵ, ζ (in that order)
    greek_coefs = GreekCoefs(convert.(FT2, greek_coefs[3,:]), 
                             convert.(FT2, greek_coefs[1,:]), 
                             convert.(FT2, greek_coefs[5,:]), 
                             convert.(FT2, greek_coefs[2,:]), 
                             convert.(FT2, greek_coefs[6,:]), 
                             convert.(FT2, greek_coefs[4,:]))

    # Return the packaged AerosolOptics object
    return AerosolOptics(greek_coefs=greek_coefs, ω̃=FT2(avg_C_scatt/avg_C_ext), k=FT2(avg_C_ext), fᵗ=FT(1)) 
end


"""
    $(FUNCTIONNAME)(l, ν₁, ν₂, ν₂_positive_flag, k, N_max, ab_pairs, an_m_bn, an_p_bn, wigner_A, wigner_B)

Reference: Suniti Sanghavi 2014, https://doi.org/10.1016/j.jqsrt.2013.12.015

Compute Sl_νν, given by Eq 22
Input: l, ν₁, ν₂, ν₂_positive_flag (to differentiate b/w 0 and -0), 
k, N_max, ab_pairs (precomputed), an_m_bn, an_p_bn, wigner_A, wigner_B, 
Output: Complex{Float64}
"""
function compute_Sl(l::Integer, ν₁::Integer, ν₂::Integer, ν₂_positive_flag::Bool, 
                    k, N_max, ab_pairs, an_m_bn, an_p_bn, wigner_A, wigner_B)

    @assert ((ν₁, ν₂) == (2, -2) || (ν₁, ν₂) == (2, 2) || (ν₁, ν₂) == (0, 0) || (ν₁, ν₂) == (0, 2)) "Not a valid input to calculate Sl"
    @assert ((ν₂ == 2 && ν₂_positive_flag == true) || (ν₂ == -2 && ν₂_positive_flag == false) || ν₂ == 0) "Sign flag doesn't match ν₂"

    # Extract conjugate tranpose pairs
    mat_anam, mat_anbm, mat_bnam, mat_bnbm = ab_pairs

    # Calculate the coefficient beforehand
    ll = l-1
    coef = (2ll+1)*π / k^2

    # Accumulate the terms of the Sl expression in these variables
    first_term = 0.0    # First summation over n
    second_term = 0.0   # Second summation over n

    # Use wigner_A or wigner_B, depending on value of ν
    wig_lnm = ν₁ == 0 ? wigner_A : wigner_B

    # Summation over n from 1 to N_max
    @inbounds for n = 1:N_max

        # Obtain appropriate m range to loop over
        m_star, m_max = (max(ll-n,n+1), min(ll+n,N_max))
        
        # Compute and accumulate the first/second terms if (ν₁, ν₂) == (0, 2)
        if (ν₁, ν₂) == (0, 2)

            # Loop over m
            @inbounds for m = m_star:m_max

                # First term
                # [(-1)^(l+m+n) * (an⋆ + bn⋆)(am - bm) + (am⋆ + bm⋆)(an - bn)] × A_lnm × B_lnm
                avg = (exp_m1(ll + n + m)) * (mat_anam[m,n] + mat_bnam[m,n] - mat_anbm[m,n] - mat_bnbm[m,n]) + (mat_anam[m,n]' - mat_bnam[m,n]' + mat_anbm[m,n]' - mat_bnbm[m,n]')
                first_term += avg * 2 * (2n+1) * (2m+1) * wigner_A[m, n, l] * wigner_B[m, n, l]
            end

            # Second term
            # (an⋆ + bn⋆)(an - bn) × A_lnm × B_lnm
            anan,anbn,bnan,bnbn = (mat_anam[n,n], mat_anbm[n,n], 
                                   mat_bnam[n,n], mat_bnbm[n,n])
            avg = anan - anbn + bnan - bnbn
            second_term += avg * 2 * (2n+1) * (2n+1) * wigner_A[n, n, l] * wigner_B[n, n, l]
        
        # Compute and accumulate the first/second terms if (ν₁, ν₂) ≠ (0, 2)
        else

            # Small change in following equations, depending on ν₂ parity (See eq. 34)

            # Loop over m
            @inbounds for m = m_star:m_max

                # First term
                # [(-1)^(l+m+n)? * Re[(an⋆ ± bn⋆)(am ± bm)] * A/B_lnm]
                anam,bnbm,anbm,bnam = (mat_anam[m,n], mat_bnbm[m,n], 
                                        mat_anbm[m,n], mat_bnam[m,n])
                real_avg = ν₂_positive_flag ? real(anam+anbm+bnam+bnbm) : real(anam-anbm-bnam+bnbm);
                coef_lnm2 = 2 * (2m + 1)  * (2n + 1) * (ν₂_positive_flag ? 1 : (exp_m1(ll + n + m)))
                first_term += (real_avg * (wig_lnm[m, n, l])^2 * coef_lnm2)
            end

            # Second term
            # (-1)^l * |an-bn|² * B_lnm

            coef_lnn2 = (2n + 1)^2 * (ν₂_positive_flag ? 1 : (exp_m1(ll)))
            if (ν₂_positive_flag)
                second_term += an_p_bn[n] * (wig_lnm[n, n, l])^2 * coef_lnn2
            else 
                second_term += an_m_bn[n] * (wig_lnm[n, n, l])^2 * coef_lnn2
            end
        end
    end

    # Combine the result, add coefficients, and return
    return coef * (first_term + second_term)
end