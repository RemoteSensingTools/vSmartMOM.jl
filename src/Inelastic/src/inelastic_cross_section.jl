function compute_effective_coefficents!(ν_eff, T, mol::MolecularConstants{FT}) where {FT}#molecules::Array{MolecularConstants{FT}}) where {FT}   
    #@unpack Y = mol
    @unpack α̅,  γ̅, α_prime, γ_prime, ϵ, ϵ_prime = mol.effCoeff
    @unpack α̅₀₀, γ̅₀₀, ω₀, α_b, α_c, α₀₀_prime, γ₀₀_prime = mol.PolTensor
    #Computing α̅
    α̅ = α̅₀₀*(1 + α_b*T + α_c*T^2)/(1-(2π*c*ν_eff/ω₀)^2)
    γ̅ = γ̅₀₀
    ϵ = α̅/γ̅ 
    α_prime = α₀₀_prime * sqrt(mol.Y[1,2]/mol.Y[2,1]) # see Eqs (36a-39b) of Buldakov et al. 1996 (Spectrochimica Acta Part A) 
    γ_prime = γ₀₀_prime * sqrt(mol.Y[1,2]/mol.Y[2,1]) # see Eqs (36a-39b) of Buldakov et al. 1996 (Spectrochimica Acta Part A) 
    γ_C_Rayl = 3/(45(ϵ)^2+4)
    γ_C_RotRaman = 3/4
    ϵ_prime = α_prime/γ_prime
    γ_C_VibRaman =  3/(45(ϵ_prime)^2+4)
    γ_C_RoVibRaman = 3/4
    rho_depol_Rayl = 2γ_C_Rayl/(1 + γ_C_Rayl) # depolarization factor, as defined in Table 1 (Spurr 2006)
    rho_depol_RotRaman = 2γ_C_RotRaman/(1 + γ_C_RotRaman) # depolarization factor, as defined in Table 1 (Spurr 2006)
    rho_depol_VibRaman = 2γ_C_VibRaman/(1 + γ_C_VibRaman) # depolarization factor, as defined in Table 1 (Spurr 2006)
    rho_depol_RoVibRaman = 2γ_C_RoVibRaman/(1 + γ_C_RoVibRaman) # depolarization factor, as defined in Table 1 (Spurr 2006)

    @pack! mol.effCoeff = rho_depol_Rayl, rho_depol_RotRaman, rho_depol_VibRaman, rho_depol_RoVibRaman,
        α̅, γ̅, α_prime, γ_prime, ϵ, ϵ_prime, γ_C_Rayl, γ_C_RotRaman, γ_C_VibRaman, γ_C_RoVibRaman
end

#Compute elastic scattering cross-section (Cabannes line)
function compute_σ_Rayl_coeff!(mol::MolecularConstants{FT}) where {FT}#ν, molecules::Array{MolecularConstants{FT}}) where {FT}
    @unpack α̅, γ_C_Rayl, σ_Rayl_coeff = mol.effCoeff
    σ_Rayl_coeff = 128π^5 * α̅^2 * (1+2*γ_C_Rayl)/(3-4*γ_C_Rayl)# * ν^4
    #@show σ_Rayl_coeff,  α̅^2 
    @pack! mol.effCoeff = σ_Rayl_coeff,α̅,γ_C_Rayl
end

#Compute elastic scattering cross-section (Cabannes line)
function compute_σ_Rayl_VibRaman_coeff_hires!(T, mol::MolecularConstants{FT}; Jmax=30) where {FT}#ν, molecules::Array{MolecularConstants{FT}}) where {FT}
    @unpack α̅, γ̅, α_prime, γ_prime, E_vJ, σ_Rayl_coeff_hires, σ_VibRaman_coeff_0to1_hires, σ_VibRaman_coeff_1to0_hires, Δν̃_Rayl_coeff_hires, Δν̃_VibRaman_coeff_0to1_hires, Δν̃_VibRaman_coeff_1to0_hires = mol.effCoeff
    @unpack gₛ = mol
    σ_Rayl_coeff_hires = OffsetArray(zeros(FT, Jmax+1), 0:Jmax);
    Δν̃_Rayl_coeff_hires = OffsetArray(zeros(FT, Jmax+1), 0:Jmax);
    σ_VibRaman_coeff_0to1_hires = OffsetArray(zeros(FT, Jmax+1), 0:Jmax);
    σ_VibRaman_coeff_1to0_hires = OffsetArray(zeros(FT, Jmax+1), 0:Jmax);
    Δν̃_VibRaman_coeff_0to1_hires = OffsetArray(zeros(FT, Jmax+1), 0:Jmax);
    Δν̃_VibRaman_coeff_1to0_hires = OffsetArray(zeros(FT, Jmax+1), 0:Jmax);

    kᵥ = 128*π^5
    #σ_Rayl_coeff_hires = 128π^5 * α̅^2 * (1+2*γ_C_Rayl)/(3-4*γ_C_Rayl)# * ν^4
    Z_pf = 0
    for Ji in 0:Jmax   
        b_JJ   = Ji*(Ji+1)/((2Ji-1)*(2Ji+3))
        g_N = gₛ[2-Ji%2]    

        # Rayleigh
        vi = 0
        vf = 0 
        Ni_by_N = exp(-h*c*E_vJ[vi,Ji]/(k_B*T))  
        γ_C = 3/(4+45(α̅/(b_JJ*γ̅))^2)
        Z_pf += g_N * (2Ji+1) * (exp(-h*c*E_vJ[0,Ji]/(k_B*T)) + exp(-h*c*E_vJ[1,Ji]/(k_B*T))) #accounting for both the ground- and the first excited vibrational state - depending on the ambient temperature (e.g. on hot exoplanets), more states may have to be accounted for   
        # Molecular energy level
        Δν =  E_vJ[vi, Ji]
        # Wavenumber change of scattered light
        Δνₛ = Δν 
    
        σ_Rayl_coeff_hires[Ji] = kᵥ * g_N * (2Ji+1) * Ni_by_N * α̅^2 * (1+2*γ_C)/(3-4*γ_C)# * νₛ[vi, Ji, vf, Jf]^4
        Δν̃_Rayl_coeff_hires[Ji] = Δνₛ

        # Vib. Raman 0->1
        vi = 0
        vf = 1 
        Ni_by_N = exp(-h*c*E_vJ[vi,Ji]/(k_B*T))  
        γ_C = 3/(4+45(α_prime/(b_JJ*γ_prime))^2)
        # Molecular energy level
        Δν =  E_vJ[vf, Ji] - E_vJ[vi, Ji]
        # Wavenumber change of scattered light
        Δνₛ = -Δν 
    
        σ_VibRaman_coeff_0to1_hires[Ji] = kᵥ * g_N * (2Ji+1) * Ni_by_N * α_prime^2 * (1+2*γ_C)/(3-4*γ_C)# * νₛ[vi, Ji, vf, Jf]^4
        Δν̃_VibRaman_coeff_0to1_hires[Ji] = Δνₛ
        #@show Ji, α_prime^2, (1+2*γ_C)/(3-4*γ_C), α_prime^2 * (1+2*γ_C)/(3-4*γ_C)
        # Vib. Raman 1->0
        vi = 1
        vf = 0 
        Ni_by_N = exp(-h*c*E_vJ[vi,Ji]/(k_B*T))  
        
        # Molecular energy level
        Δν =  E_vJ[vf, Ji] - E_vJ[vi, Ji]
        # Wavenumber change of scattered light
        Δνₛ = -Δν 
    
        σ_VibRaman_coeff_1to0_hires[Ji] = kᵥ * g_N * (2Ji+1) * Ni_by_N * α_prime^2 * (1+2*γ_C)/(3-4*γ_C)# * νₛ[vi, Ji, vf, Jf]^4
        Δν̃_VibRaman_coeff_1to0_hires[Ji] = Δνₛ
    end
    σ_Rayl_coeff_hires /= Z_pf
    σ_VibRaman_coeff_0to1_hires /= Z_pf
    σ_VibRaman_coeff_1to0_hires /= Z_pf
    #@show Z_pf
    @pack! mol.effCoeff = σ_Rayl_coeff_hires, σ_VibRaman_coeff_0to1_hires, σ_VibRaman_coeff_1to0_hires, Δν̃_Rayl_coeff_hires, Δν̃_VibRaman_coeff_0to1_hires, Δν̃_VibRaman_coeff_1to0_hires
end

#Compute energy levels [in wavenumbers [cm^{-1}]] corresponding to v={0, 1, 2} and J={0, 1, 2,..., 10}
function compute_energy_levels!(mol::MolecularConstants{FT}; vmax=2, Jmax=30) where {FT}#molecules::Array{MolecularConstants{FT}}) where {FT}
    #for mol in molecules
        @unpack E_vJ = mol.effCoeff
        @unpack Y = mol
        E_vJ = OffsetArray(zeros(FT, vmax+1,Jmax+1), 0:vmax, 0:Jmax);
        for v in 0:vmax, J in 0:Jmax
            E₁ = J*(J+1)
            for k in 1:5
                E₂ = (v+0.5)^(k-1)
                for l in 1:5
                    E₃ = E₁^(l-1) * E₂ * Y[k,l]
                    E_vJ[v,J] += E₃
                end
            end
        end
        @pack! mol.effCoeff = E_vJ
    #end
end

#Compute vibrational Raman scattering coefficient (for Δν=±1)
function compute_σ_VibRaman_coeff!(T, mol::MolecularConstants{FT}; vmax=2, Jmax=30) where {FT}#ν, molecules::Array{MolecularConstants{FT}}) where {FT}
    @unpack α_prime, γ_C_VibRaman, E_vJ,σ_VibRaman_coeff_0to1,σ_VibRaman_coeff_1to0, Δν̃_VibRaman_coeff_0to1, Δν̃_VibRaman_coeff_1to0 = mol.effCoeff
    
    #σ_Rayl_coeff = (128/3)π^5 * α̅^2 * F_King# * ν^4
    #Stokes
    vi = 0
    Δν̃ = E_vJ[1,0]-E_vJ[0,0]
    Nvib = (1 - exp(-h*c*Δν̃/(k_B*T)))^(-1)
    #sq_bₖ = h/(8π²c*E_νJ[vi,0])#bₖ = √(h/8π²cνₖ)
    σ_VibRaman_coeff_0to1 = 128π^5 * α_prime^2 * Nvib * (1+2*γ_C_VibRaman)/(3-4*γ_C_VibRaman) 
    Δν̃_VibRaman_coeff_0to1 = -Δν̃
    #Anti-Stokes
    vi = 1

    Nvib = (exp(h*c*Δν̃/(k_B*T))-1)^(-1)
    #sq_bₖ = h/(8π²c*E_νJ[vi,0])#bₖ = √(h/8π²cνₖ)
    σ_VibRaman_coeff_1to0 = 128π^5 * α_prime^2 * Nvib * (1+2*γ_C_VibRaman)/(3-4*γ_C_VibRaman) 
    Δν̃_VibRaman_coeff_1to0 = Δν̃

    #@show σ_VibRaman_coeff_0to1,σ_VibRaman_coeff_1to0  
    @pack! mol.effCoeff = σ_VibRaman_coeff_0to1, σ_VibRaman_coeff_1to0, Δν̃_VibRaman_coeff_0to1, Δν̃_VibRaman_coeff_1to0
end

function compute_σ_RoVibRaman_coeff!(T, mol::MolecularConstants{FT}; vmax=2, Jmax=30) where {FT}#molecules::Array{MolecularConstants{FT}}) where {FT}
    kᵥ = (256/27)*π^5
    #for mol in molecules
    @unpack γ̅, γ_prime,
        E_vJ, σ_RoRaman_coeff_JtoJm2,σ_RoRaman_coeff_JtoJp2,
        σ_RoVibRaman_coeff_0to1_JtoJm2,σ_RoVibRaman_coeff_0to1_JtoJp2, 
        σ_RoVibRaman_coeff_1to0_JtoJm2,σ_RoVibRaman_coeff_1to0_JtoJp2,
        Δν̃_RoRaman_coeff_JtoJm2, Δν̃_RoRaman_coeff_JtoJp2,
        Δν̃_RoVibRaman_coeff_0to1_JtoJm2, Δν̃_RoVibRaman_coeff_0to1_JtoJp2,
        Δν̃_RoVibRaman_coeff_1to0_JtoJm2, Δν̃_RoVibRaman_coeff_1to0_JtoJp2 = mol.effCoeff
    @unpack gₛ = mol 
    
    #σ_VibRaman_coeff = OffsetArray(zeros(FT, vmax+1,vmax+1), 0:vmax, 0:vmax); 
    σ_RoVibRaman_coeff_0to1_JtoJm2 = OffsetArray(zeros(FT, Jmax+1), 0:Jmax); #zeros(3, 11, 3, 11) #(vi, ji, vf, jf)
    σ_RoVibRaman_coeff_0to1_JtoJp2 = OffsetArray(zeros(FT, Jmax+1), 0:Jmax);
    σ_RoVibRaman_coeff_1to0_JtoJm2 = OffsetArray(zeros(FT, Jmax+1), 0:Jmax); #zeros(3, 11, 3, 11) #(vi, ji, vf, jf)
    σ_RoVibRaman_coeff_1to0_JtoJp2 = OffsetArray(zeros(FT, Jmax+1), 0:Jmax);
    σ_RoRaman_coeff_JtoJm2 = OffsetArray(zeros(FT, Jmax+1), 0:Jmax);
    σ_RoRaman_coeff_JtoJp2 = OffsetArray(zeros(FT, Jmax+1), 0:Jmax);

    Δν̃_RoVibRaman_coeff_0to1_JtoJm2 = OffsetArray(zeros(FT, Jmax+1), 0:Jmax); #zeros(3, 11, 3, 11) #(vi, ji, vf, jf)
    Δν̃_RoVibRaman_coeff_0to1_JtoJp2 = OffsetArray(zeros(FT, Jmax+1), 0:Jmax);
    Δν̃_RoVibRaman_coeff_1to0_JtoJm2 = OffsetArray(zeros(FT, Jmax+1), 0:Jmax); #zeros(3, 11, 3, 11) #(vi, ji, vf, jf)
    Δν̃_RoVibRaman_coeff_1to0_JtoJp2 = OffsetArray(zeros(FT, Jmax+1), 0:Jmax);
    Δν̃_RoRaman_coeff_JtoJm2 = OffsetArray(zeros(FT, Jmax+1), 0:Jmax);
    Δν̃_RoRaman_coeff_JtoJp2 = OffsetArray(zeros(FT, Jmax+1), 0:Jmax);
    #Δνₛ = OffsetArray(zeros(FT, vmax+1,Jmax+1,vmax+1,Jmax+1), 0:vmax, 0:Jmax, 0:vmax, 0:Jmax); #zeros(3, 11, 3, 11)
        
    Z_pf = 0
    for Ji in 0:Jmax         
        b_JJm2 = 3Ji*(Ji-1)/(2*(2Ji+1)*(2Ji-1))
        #b_JJ   = Ji*(Ji+1)/((2Ji-1)*(2Ji+3))
        b_JJp2 = 3*(Ji+1)*(Ji+2)/(2*(2Ji+1)*(2Ji+3))
        g_N = gₛ[2-Ji%2]    
        
        #Rotational Raman scattering
        vi = 0
        vf = 0 
        Ni_by_N = exp(-h*c*E_vJ[vi,Ji]/(k_B*T))  
        #@show vi, Ji, Ni_by_N
        Z_pf += g_N * (2Ji+1) * (exp(-h*c*E_vJ[0,Ji]/(k_B*T)) + exp(-h*c*E_vJ[1,Ji]/(k_B*T))) #accounting for both the ground- and the first excited vibrational state - depending on the ambient temperature (e.g. on hot exoplanets), more states may have to be accounted for   
        
        Jf=Ji-2
        if (Ji-2)>=0 #anti-Stokes
            #Change in molecular energy level
            Δν = E_vJ[vf, Jf] - E_vJ[vi, Ji]
            #Wavenumber change of scattered light
            Δνₛ = - Δν 
        
            σ_RoRaman_coeff_JtoJm2[Ji] = kᵥ * g_N * (2Ji+1) * b_JJm2 * Ni_by_N * γ̅^2# * νₛ[vi, Ji, vf, Jf]^4
            Δν̃_RoRaman_coeff_JtoJm2[Ji] = Δνₛ
        else
            σ_RoRaman_coeff_JtoJm2[Ji] = 0
            Δν̃_RoRaman_coeff_JtoJm2[Ji] = 0
        end        
        Jf=Ji+2
        if (Ji+2)<=Jmax #Stokes
            #Change in molecular energy level
            Δν = E_vJ[vf, Jf] - E_vJ[vi, Ji]
            #Wavenumber change of scattered light
            Δνₛ = - Δν 

            σ_RoRaman_coeff_JtoJp2[Ji] = kᵥ * g_N * (2Ji+1) * b_JJp2 * Ni_by_N * γ̅^2# * νₛ[vi, Ji, vf, Jf]^4 
            Δν̃_RoRaman_coeff_JtoJp2[Ji] = Δνₛ
        else
            σ_RoRaman_coeff_JtoJp2[Ji] = 0
            Δν̃_RoRaman_coeff_JtoJp2[Ji] = 0
        end
    
        #Rovibrational Raman scattering: Stokes
        vi = 0
        vf = 1
            
        Ni_by_N = exp(-h*c*E_vJ[vi,Ji]/(k_B*T))   
    
        Jf=Ji-2
        #b_v = (1/2π)*(h/(2c*E_vJ[vf,Jf]))^0.5 
        if (Ji-2)>=0 #anti-Stokes
            #Change in molecular energy level
            Δν = E_vJ[vf, Jf] - E_vJ[vi, Ji]
            #Wavenumber change of scattered light
            Δνₛ = - Δν 
        
            σ_RoVibRaman_coeff_0to1_JtoJm2[Ji] = kᵥ * g_N * (2Ji+1) * b_JJm2 * Ni_by_N * (γ_prime)^2# * νₛ[vi, Ji, vf, Jf]^4 
            Δν̃_RoVibRaman_coeff_0to1_JtoJm2[Ji] = Δνₛ
        else
            σ_RoVibRaman_coeff_0to1_JtoJm2[Ji] = 0
            Δν̃_RoVibRaman_coeff_0to1_JtoJm2[Ji] = 0
        end        
        Jf=Ji+2
        #b_v = (1/2π)*(h/(2c*E_vJ[vf,Jf]))^0.5 
        if (Ji+2)<=Jmax #Stokes
            #Change in molecular energy level
            Δν = E_vJ[vf, Jf] - E_vJ[vi, Ji]
            #Wavenumber change of scattered light
            Δνₛ = - Δν 
        
            σ_RoVibRaman_coeff_0to1_JtoJp2[Ji] = kᵥ * g_N * (2Ji+1) * b_JJp2 * Ni_by_N * (γ_prime)^2# * νₛ[vi, Ji, vf, Jf]^4 
            Δν̃_RoVibRaman_coeff_0to1_JtoJp2[Ji] = Δνₛ
        else
            σ_RoVibRaman_coeff_0to1_JtoJp2[Ji] = 0
            Δν̃_RoVibRaman_coeff_0to1_JtoJp2[Ji] = 0
        end
    
        #Rovibrational Raman scattering: anti-Stokes
        vi = 1
        vf = 0

        Ni_by_N = exp(-h*c*E_vJ[vi,Ji]/(k_B*T))   
    
        Jf=Ji-2
        if (Ji-2)>=0 #anti-Stokes
            #Change in molecular energy level
            Δν = E_vJ[vf, Jf] - E_vJ[vi, Ji]
            #Wavenumber change of scattered light
            Δνₛ = - Δν 
    
            σ_RoVibRaman_coeff_1to0_JtoJm2[Ji] = kᵥ * g_N * (2Ji+1) * b_JJm2 * Ni_by_N * (γ_prime)^2# * νₛ[vi, Ji, vf, Jf]^4 
            Δν̃_RoVibRaman_coeff_1to0_JtoJm2[Ji] = Δνₛ
        else
            σ_RoVibRaman_coeff_1to0_JtoJm2[Ji] = 0
            Δν̃_RoVibRaman_coeff_1to0_JtoJm2[Ji] = 0
        end        
        Jf=Ji+2
        if (Ji+2)<=Jmax #Stokes
            #Change in molecular energy level
            Δν = E_vJ[vf, Jf] - E_vJ[vi, Ji]
            #Wavenumber change of scattered light
            Δνₛ = - Δν 
    
            σ_RoVibRaman_coeff_1to0_JtoJp2[Ji] = kᵥ * g_N * (2Ji+1) * b_JJp2 * Ni_by_N * (γ_prime)^2# * νₛ[vi, Ji, vf, Jf]^4 
            Δν̃_RoVibRaman_coeff_1to0_JtoJp2[Ji] = Δνₛ
        else
            σ_RoVibRaman_coeff_1to0_JtoJp2[Ji] = 0
            Δν̃_RoVibRaman_coeff_1to0_JtoJp2[Ji] = 0
        end
    end
    σ_RoRaman_coeff_JtoJm2 /= Z_pf
    σ_RoRaman_coeff_JtoJp2 /= Z_pf
    σ_RoVibRaman_coeff_0to1_JtoJm2 /= Z_pf
    σ_RoVibRaman_coeff_0to1_JtoJp2 /= Z_pf
    σ_RoVibRaman_coeff_1to0_JtoJm2 /= Z_pf
    σ_RoVibRaman_coeff_1to0_JtoJp2 /= Z_pf

    #@show Z_pf
    @pack! mol.effCoeff = σ_RoRaman_coeff_JtoJm2, σ_RoRaman_coeff_JtoJp2,
        σ_RoVibRaman_coeff_0to1_JtoJm2, σ_RoVibRaman_coeff_0to1_JtoJp2, 
        σ_RoVibRaman_coeff_1to0_JtoJm2, σ_RoVibRaman_coeff_1to0_JtoJp2,
        Δν̃_RoRaman_coeff_JtoJm2, Δν̃_RoRaman_coeff_JtoJp2,
        Δν̃_RoVibRaman_coeff_0to1_JtoJm2, Δν̃_RoVibRaman_coeff_0to1_JtoJp2, 
        Δν̃_RoVibRaman_coeff_1to0_JtoJm2, Δν̃_RoVibRaman_coeff_1to0_JtoJp2  
end
