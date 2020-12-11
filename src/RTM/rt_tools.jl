# atmospheric RTM
function run_RTM(polarization_type, sza, vza, vaz, τRayl,ϖRayl, τAer, ϖAer, fᵗ, qp_μ, wt_μ, Ltrunc, aerosol_optics, GreekRayleigh)
    FT = eltype(τRayl)

    #Output variables: Reflected and transmitted solar irradiation at TOA and BOA respectively
    R = zeros(length(vza),4)
    T = zeros(length(vza),4)    
    μ0 = cosd(sza)
    iμ0 = nearest_point(qp_μ, μ0) # input μ0 = cos(SZA)
    I0 = [1, 0, 0, 0] #assuming completely unpolarized incident stellar radiation
    #get vertical grid
    #get solar+viewing geometry, compute streams
    #compute Aersol SSP
    #compute Rayleigh SSP
    Nz = length(τRayl)
    Naer = length(aerosol_optics)
    for m=0:Ltrunc-1
        @show m
        if (m==0)
            weight=0.5
        else
            weight=1.0
        end
        #compute Zmp_Aer, Zpp_Aer, Zmp_Rayl, Zpp_Rayl
        # For m>=3, Rayleigh matrices will be 0, can catch with if statement if wanted 
        Rayl𝐙⁺⁺, Rayl𝐙⁻⁺ = PhaseFunction.compute_Z_moments(polarization_type, qp_μ, GreekRayleigh, m);
        dims = size(Rayl𝐙⁺⁺)
        nAer = length(aerosol_optics)
        Nquad4 = dims[1]
        Aer𝐙⁺⁺ = [zeros(FT,dims) for i in 1:nAer]
        Aer𝐙⁻⁺ = similar(Aer𝐙⁺⁺)

        for i = 1:nAer
            Aer𝐙⁺⁺[i], Aer𝐙⁻⁺[i] = PhaseFunction.compute_Z_moments(polarization_type, qp_μ, aerosol_optics[i].greek_coefs, m)
        end
        
        # Homogenous R and T matrices
        r⁻⁺ = zeros(FT, dims)
        t⁺⁺ = zeros(FT, dims)
        r⁺⁻ = zeros(FT, dims)
        t⁻⁻ = zeros(FT, dims)
        # Composite layer R and T matrices
        R⁻⁺ = zeros(dims)
        R⁺⁻ = zeros(dims)
        T⁺⁺ = zeros(dims)
        T⁻⁻ = zeros(dims)

        kn=0
        # loop over vertical layers:
        for iz=1:Nz  #Count from TOA to BOA
            τ, ϖ, Z⁺⁺, Z⁻⁺ = construct_atm_layer(τRayl[iz], τAer[iz,:], ϖRayl[iz], ϖAer, fᵗ, Rayl𝐙⁺⁺, Rayl𝐙⁻⁺, Aer𝐙⁺⁺, Aer𝐙⁻⁺)
            dτ_max = minimum([τ, 0.2*minimum(qp_μ)])
            dτ, ndoubl = doubling_number(dτ_max, τ)
            scatter=false
            if (sum(τAer)>1.e-8)
                scatter=true
            elseif (τRayl[iz]>1.e-8) & (m<3)
                scatter=true
            end        
            if (scatter)
                @timeit "elemental" rt_elemental!(dτ, ϖ, Z⁺⁺, Z⁻⁺, m, ndoubl, scatter,qp_μ, wt_μ, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻)
                @timeit "doubling" rt_doubling!(dτ, τ, ndoubl, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻)
            else
                r⁻⁺ = 0
                r⁺⁻ = 0
                for i = 1:Nquad4
                    ii=1+floor(Int,(i-1)/4)
                    t⁺⁺[i,i] = exp(-τ/qp_μ[ii])
                    t⁻⁻[i,i] = exp(-τ/qp_μ[ii])
                end
            end
            kn = get_kn(kn, scatter, iz)
            
            if (iz==1)
                T⁺⁺ = t⁺⁺
                T⁻⁻ = t⁻⁻
                R⁻⁺ = r⁻⁺
                R⁺⁻ = r⁺⁻
            else
                @timeit "interaction" R⁻⁺, T⁺⁺, R⁺⁻, T⁻⁻ = rt_interaction(kn, R⁻⁺, T⁺⁺, R⁺⁻, T⁻⁻, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻)
            end
        end #z
        

        # include surface function
        # TBD
        
        for i = 1:length(vza)
            iμ = nearest_point(qp_μ, cosd(vza[i])) #input vaz, vza as arrays
            #@show i, vza[i], cosd(vza[i]), iμ, qp_μ[iμ]
            # compute bigCS
            cos_m_phi = cosd(m * vaz[i])
            sin_m_phi = sind(m * vaz[i])
            bigCS = Diagonal([cos_m_phi, cos_m_phi, sin_m_phi, sin_m_phi])
            # accumulate Fourier moments after azimuthal weighting
            #Measurement at the TOA
            st_iμ  = (iμ-1)*4+1
            st_iμ0 = (iμ0-1)*4+1
            #@show st_iμ:st_iμ+3, iμ0,st_iμ0:st_iμ0+3
            #@show size(R⁻⁺)
            R[i,:] += weight * bigCS * (R⁻⁺[st_iμ:st_iμ+3, st_iμ0:st_iμ0+3]/wt_μ[iμ0]) * I0
            #Measurement at the BOA
            T[i,:] += weight * bigCS * (T⁺⁺[st_iμ:st_iμ+3, st_iμ0:st_iμ0+3]/wt_μ[iμ0]) * I0     
            #if m==0
            #    @show bigCS
            #    @show m, i, iμ, bigCS[1,1], weight*R⁻⁺[(iμ-1)*4+1, (iμ0-1)*4+1]/wt_μ[iμ0]   
            #end
        end
    end  #m

    
    print_timer()
    reset_timer!()

    return R, T  
end

function get_kn(kn, scatter, iz)
    if (iz==1)
        if (scatter)
            kn=4
        else
            kn=1
        end
    else 
        if (kn==1) & (!scatter)
            kn = 1
        elseif (kn==1) & (scatter)
            kn = 2
        elseif (kn>1) & (!scatter)
            kn = 3
        elseif (kn>1) & (scatter)
            kn = 4
        end 
    end

    return kn
end