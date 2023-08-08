"Construct constants for N₂:"
function getMolecularConstants(::N₂, vmr::FT) where {FT}
    @info "Constructing N₂ constants, using VMR (unitless) of " vmr
    @assert 0 ≤ vmr ≤ 1 "VMR has to be between [0,1]"
    @assert FT <: Float64 "Requires 64 bit floating type"
    p = PolarizationTensor{FT}(
            α̅₀₀     = 1.7406e-24, #[cm^3]
            α₀₀_prime = 1.86e-24, #[cm^3] #to be multiplied by √(Bₑ/ωₑ)
            ω₀      = 2.6049e16, 
            α_b     = 1.8e-6, 
            α_c     = 0.0,
            γ̅₀₀       = 0.71e-24,  #[cm^3] #to be multiplied by √(Bₑ/ωₑ)
            γ₀₀_prime = 2.23e-24)   #[cm^3]
    
    Y = zeros(FT, 5,5)
    # Fill in all non-zero elements:
    Y[1,2] =  1.99824  # rotational constant in equilibrium position (cm^{-1})
    Y[1,3] = -5.76E-6  # centrifugal distortion constant (cm^{-1})
    Y[2,1] =  2358.57  # vibrational constant - first term (cm^{-1}) 
    Y[2,2] = -0.017318 # rotational constant - first term (cm^{-1})
    Y[3,1] = -14.324
    Y[4,1] = -2.26e-3  #vibrational constant - third term (cm^{-1})

    gₛ = [3, 6] #{odd J_initial, even J_initial}
    MolecularConstants(vmr=vmr,
                        PolTensor=p, 
                        Y=Y,
                        gₛ = gₛ, 
                        effCoeff=EffectiveCoefficients(
                                T = FT(273), 
                                α̅ = FT(0), 
                                γ̅ = FT(0),
                                α_prime = FT(0),
                                γ_prime = FT(0),     
                                ϵ = FT(0), 
                                ϵ_prime = FT(0), 
                                γ_C_Rayl = FT(0), 
                                γ_C_RotRaman = FT(0), 
                                γ_C_VibRaman = FT(0), 
                                γ_C_RoVibRaman = FT(0), 
                                rho_depol_Rayl = FT(0), 
                                rho_depol_RotRaman = FT(0), 
                                rho_depol_VibRaman = FT(0), 
                                rho_depol_RoVibRaman = FT(0),
                                σ_Rayl_coeff = FT(0), 
                                σ_Rayl_coeff_hires = FT[1.0], 
                                Δν̃_Rayl_coeff_hires = FT[1.0], 
                                σ_VibRaman_coeff_0to1 = FT(0), 
                                σ_VibRaman_coeff_1to0 = FT(0), 
                                σ_VibRaman_coeff_0to1_hires = FT[1.0], 
                                σ_VibRaman_coeff_1to0_hires = FT[1.0],
                                Δν̃_VibRaman_coeff_0to1 = FT(0), 
                                Δν̃_VibRaman_coeff_1to0 = FT(0), 
                                Δν̃_VibRaman_coeff_0to1_hires = FT[1.0], 
                                Δν̃_VibRaman_coeff_1to0_hires = FT[1.0],
                                σ_RoVibRaman_coeff_0to1_JtoJm2 = FT[1.0],
                                σ_RoVibRaman_coeff_0to1_JtoJp2 = FT[1.0],
                                σ_RoVibRaman_coeff_1to0_JtoJm2 = FT[1.0],
                                σ_RoVibRaman_coeff_1to0_JtoJp2 = FT[1.0],
                                Δν̃_RoVibRaman_coeff_0to1_JtoJm2 = FT[1.0],
                                Δν̃_RoVibRaman_coeff_0to1_JtoJp2 = FT[1.0],
                                Δν̃_RoVibRaman_coeff_1to0_JtoJm2 = FT[1.0],
                                Δν̃_RoVibRaman_coeff_1to0_JtoJp2 = FT[1.0],
                                σ_RoRaman_coeff_JtoJm2 = FT[1.0],
                                σ_RoRaman_coeff_JtoJp2 = FT[1.0],
                                Δν̃_RoRaman_coeff_JtoJm2 = FT[1.0],
                                Δν̃_RoRaman_coeff_JtoJp2 = FT[1.0],
                                E_vJ = FT[1.0] )
     )   # Fill with dummies for start
end

"Construct constants for O₂:"
function getMolecularConstants(::O₂, vmr::FT) where {FT}
    @info "Constructing O₂ constants, using VMR (unitless) of " vmr
    @assert 0 ≤ vmr ≤ 1 "VMR has to be between [0,1]"
    @assert FT <: Float64 "Requires 64 bit floating type"
    p = PolarizationTensor{FT}(
            α̅₀₀     = 1.5658e-24, #[cm^3]
            α₀₀_prime = 1.76e-24, #[cm^3]
            ω₀      = 2.1801e16, 
            α_b     = -2.369e-6, 
            α_c     = 8.687e-9,
            γ̅₀₀       = 1.080e-24,  #[cm^3] Ref: Asawaroengchai & Rosenblatt (1980), J. Chem. Phys. (gamma is written as beta_e in the paper) );
            γ₀₀_prime = 3.19e-24)   #[cm^3]
    
    Y = zeros(FT, 5,5)
    # Fill in all non-zero elements:
    Y[1,2] =  1.4376766  # rotational constant in equilibrium position, Bₑ (cm^{-1})
    Y[1,3] = -4.839E-6  # centrifugal distortion constant, -Dₑ (cm^{-1})
    Y[2,1] =  1580.19  # vibrational constant - first term, ωₑ (cm^{-1}) 
    Y[2,2] = -0.01590  # rotational constant - first term, -αₑ (cm^{-1})
    Y[3,1] = -11.98    # -ωₑxₑ
    Y[4,1] =  0.0      # vibrational constant - third term, γₑ (cm^{-1})

    gₛ = [1, 0] #{odd J_initial, even J_initial}
    MolecularConstants(vmr=vmr,
                        PolTensor=p, 
                        Y=Y,
                        gₛ = gₛ, 
                        effCoeff=EffectiveCoefficients(
                                T = FT(273), 
                                α̅ = FT(0),    
                                γ̅ = FT(0),
                                α_prime = FT(0),
                                γ_prime = FT(0),  
                                ϵ = FT(0), 
                                ϵ_prime = FT(0), 
                                γ_C_Rayl = FT(0), 
                                γ_C_RotRaman = FT(0), 
                                γ_C_VibRaman = FT(0), 
                                γ_C_RoVibRaman = FT(0), 
                                rho_depol_Rayl = FT(0), 
                                rho_depol_RotRaman = FT(0), 
                                rho_depol_VibRaman = FT(0), 
                                rho_depol_RoVibRaman = FT(0),
                                σ_Rayl_coeff = FT(0), 
                                σ_Rayl_coeff_hires = FT[1.0], 
                                Δν̃_Rayl_coeff_hires = FT[1.0], 
                                σ_VibRaman_coeff_0to1 = FT(0), 
                                σ_VibRaman_coeff_1to0 = FT(0), 
                                σ_VibRaman_coeff_0to1_hires = FT[1.0], 
                                σ_VibRaman_coeff_1to0_hires = FT[1.0],
                                Δν̃_VibRaman_coeff_0to1 = FT(0), 
                                Δν̃_VibRaman_coeff_1to0 = FT(0), 
                                Δν̃_VibRaman_coeff_0to1_hires = FT[1.0], 
                                Δν̃_VibRaman_coeff_1to0_hires = FT[1.0],
                                σ_RoVibRaman_coeff_0to1_JtoJm2 = FT[1.0],
                                σ_RoVibRaman_coeff_0to1_JtoJp2 = FT[1.0],
                                σ_RoVibRaman_coeff_1to0_JtoJm2 = FT[1.0],
                                σ_RoVibRaman_coeff_1to0_JtoJp2 = FT[1.0],
                                Δν̃_RoVibRaman_coeff_0to1_JtoJm2 = FT[1.0],
                                Δν̃_RoVibRaman_coeff_0to1_JtoJp2 = FT[1.0],
                                Δν̃_RoVibRaman_coeff_1to0_JtoJm2 = FT[1.0],
                                Δν̃_RoVibRaman_coeff_1to0_JtoJp2 = FT[1.0],
                                σ_RoRaman_coeff_JtoJm2 = FT[1.0],
                                σ_RoRaman_coeff_JtoJp2 = FT[1.0],
                                Δν̃_RoRaman_coeff_JtoJm2 = FT[1.0],
                                Δν̃_RoRaman_coeff_JtoJp2 = FT[1.0],
                                E_vJ = FT[1.0] )
                        )
end

"Construct constants for H₂:"
function getMolecularConstants(::H₂, vmr::FT) where {FT}
    @info "Constructing H₂ constants, using VMR (unitless) of " vmr
    @assert 0 ≤ vmr ≤ 1 "VMR has to be between [0,1]"
    @assert FT <: Float64 "Requires 64 bit floating type"
    p = PolarizationTensor{FT}(
            α̅₀₀     = 0.8032e-24, #[cm^3]
            α₀₀_prime = 0.90e-24, #[cm^3]
            ω₀      = 2.1399e16, 
            α_b     = 5.870e-6, 
            α_c     = 7.544e-9,
            γ̅₀₀       = 0.288e-24,  #[cm^3] Ref: Asawaroengchai & Rosenblatt (1980), J. Chem. Phys. (gamma is written as beta_e in the paper) );
            γ₀₀_prime = 1.02e-24)   #[cm^3]
    
    Y = zeros(FT, 5,5)
    # Fill in all non-zero elements:
    Y[1,2] =  60.853    # rotational constant in equilibrium position, Bₑ (cm^{-1})
    Y[1,3] = -0.0471    # centrifugal distortion constant, -Dₑ (cm^{-1})
    Y[2,1] =  4401.21   # vibrational constant - first term, ωₑ (cm^{-1}) 
    Y[2,2] = -3.062     # rotational constant - first term, -αₑ (cm^{-1})
    Y[3,1] = -121.33    # -ωₑxₑ
    Y[4,1] =  0.0      # vibrational constant - third term, γₑ (cm^{-1})

    gₛ = [3, 1] #{odd J_initial, even J_initial}
    MolecularConstants(vmr=vmr,
                        PolTensor=p, 
                        Y=Y,
                        gₛ = gₛ, 
                        effCoeff=EffectiveCoefficients(
                                T = FT(273), 
                                α̅ = FT(0),    
                                γ̅ = FT(0),
                                α_prime = FT(0),
                                γ_prime = FT(0),  
                                ϵ = FT(0), 
                                ϵ_prime = FT(0), 
                                γ_C_Rayl = FT(0), 
                                γ_C_RotRaman = FT(0), 
                                γ_C_VibRaman = FT(0), 
                                γ_C_RoVibRaman = FT(0), 
                                rho_depol_Rayl = FT(0), 
                                rho_depol_RotRaman = FT(0), 
                                rho_depol_VibRaman = FT(0), 
                                rho_depol_RoVibRaman = FT(0),
                                σ_Rayl_coeff = FT(0), 
                                σ_Rayl_coeff_hires = FT[1.0], 
                                Δν̃_Rayl_coeff_hires = FT[1.0], 
                                σ_VibRaman_coeff_0to1 = FT(0), 
                                σ_VibRaman_coeff_1to0 = FT(0), 
                                σ_VibRaman_coeff_0to1_hires = FT[1.0], 
                                σ_VibRaman_coeff_1to0_hires = FT[1.0],
                                Δν̃_VibRaman_coeff_0to1 = FT(0), 
                                Δν̃_VibRaman_coeff_1to0 = FT(0), 
                                Δν̃_VibRaman_coeff_0to1_hires = FT[1.0], 
                                Δν̃_VibRaman_coeff_1to0_hires = FT[1.0],
                                σ_RoVibRaman_coeff_0to1_JtoJm2 = FT[1.0],
                                σ_RoVibRaman_coeff_0to1_JtoJp2 = FT[1.0],
                                σ_RoVibRaman_coeff_1to0_JtoJm2 = FT[1.0],
                                σ_RoVibRaman_coeff_1to0_JtoJp2 = FT[1.0],
                                Δν̃_RoVibRaman_coeff_0to1_JtoJm2 = FT[1.0],
                                Δν̃_RoVibRaman_coeff_0to1_JtoJp2 = FT[1.0],
                                Δν̃_RoVibRaman_coeff_1to0_JtoJm2 = FT[1.0],
                                Δν̃_RoVibRaman_coeff_1to0_JtoJp2 = FT[1.0],
                                σ_RoRaman_coeff_JtoJm2 = FT[1.0],
                                σ_RoRaman_coeff_JtoJp2 = FT[1.0],
                                Δν̃_RoRaman_coeff_JtoJm2 = FT[1.0],
                                Δν̃_RoRaman_coeff_JtoJp2 = FT[1.0],
                                E_vJ = FT[1.0] )
                        )
end
