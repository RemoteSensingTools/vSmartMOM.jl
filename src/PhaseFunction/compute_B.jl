


function compute_B(aerosol::UnivariateAerosol, wigner_A, wigner_B, wl, r, w)

    # Find overall N_max from the maximum radius
    N_max = PhaseFunction.get_n_max(2 * π * aerosol.r_max/ wl)

    # Compute an, bn values
    an, bn = compute_anbn(aerosol::UnivariateAerosol, wl, r)

    # Compute the average cross-sectional scattering
    k = 2 * π / wl

    # Compute average cross-sectional scattering and extinction
    avg_C_scatt, avg_C_ext = compute_avg_C_scatt_ext(k, an, bn, w)

    # Only do these l's for now
    ls = 1:(2 * N_max - 1)

    # Where to store the values
    greek_coefs = zeros(6, size(ls, 1))

    # Pre-compute anbn averages
    FT2 = Complex{Float64}
    mat_anam = LowerTriangular(zeros(FT2, N_max, N_max));
    mat_bnbm = LowerTriangular(zeros(FT2, N_max, N_max));
    mat_anbm = LowerTriangular(zeros(FT2, N_max, N_max));
    mat_bnam = LowerTriangular(zeros(FT2, N_max, N_max));
    
    N_max_ = PhaseFunction.get_n_max.(2π * r/ wl)
    PhaseFunction.compute_avg_anbn!(an, bn, mat_anam, mat_bnbm, mat_anbm, mat_bnam, w, N_max, N_max_)
    an_m_bn = transpose(abs2.(an-bn)) * w
    an_p_bn = transpose(abs2.(an+bn)) * w

    elapsed_Sl = zeros(5, (2 * N_max - 1))

    # For each l
    @showprogress 1 "Computing S functions ..." for l in ls

        # Compute all necessary Sl_νν values 
        # Eq 22 in Sanghavi 2013

        Sl_00  = compute_Sl(l, 0, 0,  true,  k, N_max, an,bn, mat_anam, mat_bnbm, mat_anbm, mat_bnam, wigner_A, wigner_B, an_m_bn, an_p_bn, w)
        Sl_0m0 = compute_Sl(l, 0, 0,  false, k, N_max, an,bn, mat_anam, mat_bnbm, mat_anbm, mat_bnam, wigner_A, wigner_B, an_m_bn, an_p_bn, w)
        Sl_22  = compute_Sl(l, 2, 2,  true,  k, N_max, an,bn, mat_anam, mat_bnbm, mat_anbm, mat_bnam, wigner_A, wigner_B, an_m_bn, an_p_bn, w)
        Sl_2m2 = compute_Sl(l, 2, -2, false, k, N_max, an,bn, mat_anam, mat_bnbm, mat_anbm, mat_bnam, wigner_A, wigner_B, an_m_bn, an_p_bn, w)
        Sl_02  = compute_Sl(l, 0, 2,  true,  k, N_max, an,bn, mat_anam, mat_bnbm, mat_anbm, mat_bnam, wigner_A, wigner_B, an_m_bn, an_p_bn, w)

        # Use Sl_νν values to compute Greek coefficients 
        # Eq 24 in Sanghavi 2013

        @inbounds greek_coefs[1,l] = (1/avg_C_scatt) * (Sl_00 + Sl_0m0)
        @inbounds greek_coefs[2,l] = (1/avg_C_scatt) * (Sl_00 - Sl_0m0)
        @inbounds greek_coefs[3,l] = (1/avg_C_scatt) * (Sl_22 + Sl_2m2)
        @inbounds greek_coefs[4,l] = (1/avg_C_scatt) * (Sl_22 - Sl_2m2)
        @inbounds greek_coefs[5,l] = (1/avg_C_scatt) * real(Sl_02)
        @inbounds greek_coefs[6,l] = (1/avg_C_scatt) * imag(Sl_02)
    end

    # Create GreekCoefs object with α, β, γ, δ, ϵ, ζ (in that order)
    greek_coefs = GreekCoefs(greek_coefs[3,:], greek_coefs[1,:], greek_coefs[5,:], 
                  greek_coefs[2,:], greek_coefs[6,:], greek_coefs[4,:])

    # Return the packaged AerosolOptics object
    return AerosolOptics(greek_coefs, avg_C_scatt, avg_C_ext) 

end