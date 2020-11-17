# atmospheric RTM
function RTM(polarization_type,τRayl,ϖRayl, τAer, ϖAer, fᵗ, qp_μ, Ltrunc, Nz, aerosol_optics, GreekRayleigh)
    FT = eltype(τRayl)
    #get vertical grid
    #get solar+viewing geometry, compute streams
    #compute Aersol SSP
    #compute Rayleigh SSP

    for m=0:Ltrunc
        #compute Zmp_Aer, Zpp_Aer, Zmp_Rayl, Zpp_Rayl
        # For m>=3, Rayleigh matrices will be 0, can catch with if statement if wanted 
        Rayl𝐙⁺⁺, Rayl𝐙⁻⁺ = PhaseFunction.compute_Z_moments(polarization_type, qp_μ, GreekRayleigh, m);
        dims = size(Rayl𝐙⁺⁺)
        nAer = length(aerosol_optics)

        Aer𝐙⁺⁺ = [zeros(FT,dims) for i in 1:nAer]
        Aer𝐙⁻⁺ = similar(Aer𝐙⁺⁺)

        for i = 1:nAer
            Aer𝐙⁺⁺[i], Aer𝐙⁻⁺[i] = PhaseFunction.compute_Z_moments(polarization_type, qp_μ, aerosol_optics[i].greek_coefs, m)
        end
        
        # Homogenous R and T matrices
        r⁻⁺ = zeros(dims)
        t⁺⁺ = zeros(dims)
        r⁺⁻ = zeros(dims)
        t⁻⁻ = zeros(dims)
        # Composite layer R and T matrices
        R⁻⁺ = zeros(dims)
        R⁺⁻ = zeros(dims)
        T⁺⁺ = zeros(dims)
        T⁻⁻ = zeros(dims)

        # loop over vertical layers:
        for iz=1:Nz  #Count from TOA to BOA
            τ, ϖ, Z⁺⁺, Z⁻⁺ = construct_atm_layer(m, iz, Nquad4, Naer, τRayl[iz], τAer[iz], ϖRayl[iz], ϖAer, fᵗ, RaylZ⁺⁺, RaylZ⁻⁺, AerZ⁺⁺, AerZ⁻⁺)
            dτ_max = min(τ, min(qp_μ)/5)
            doubling_number!(Stokes_IQUV, dτ_max, τ, dτ, ndoubl)
            scatter=0
            if (sum(τAer)>1.e-8)
                scatter=1
            elseif (τRayl>1.e-8) & (m<3)
                scatter=1
            end        
            if (scatter)
                r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻ = rt_elemental(polarization_type, dτ, ϖ, bb=0, m, iz, ndoubl, scatter)
                rt_doubling!(polarization_type, dτ, τ, m, ndoubl, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻)
            else
                r⁻⁺ = 0
                r⁺⁻ = 0
                for i in 1:Nquad4
                    ii=1+floor((i-1)/4)
                t⁺⁺[i,i] = exp(-τ/qp_μ[ii])
                t⁻⁻[i,i] = exp(-τ/qp_μ[ii])
            end
            if (iz==1)
                if (scatter==1)
                    kn=4
                else
                    kn=1
                end
            elseif 
                if (kn==1) & (scatter==0)
                    kn = 1
                elseif (kn==1) & (scatter==1)
                    kn = 2
                elseif (kn>1) & (scatter==0)
                    kn = 3
                elseif (kn>1) & (scatter==1)
                    kn = 4
                end 
            end
            if (iz==1)
                T⁺⁺ = t⁺⁺
                T⁻⁻ = t⁻⁻
                R⁻⁺ = r⁻⁺
                R⁺⁻ = r⁺⁻
            elseif
                rt_interaction!(Stokes_IQUV, kn, iz, m, T⁺⁺, T⁻⁻, R⁻⁺, R⁺⁻, t⁺⁺, t⁻⁻, r⁻⁺, r⁺⁻)
            end
        end #z
    end #m
end