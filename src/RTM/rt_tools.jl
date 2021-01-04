using ..Architectures: devi, default_architecture

# atmospheric RTM
function run_RTM(pol_type, sza, vza, vaz, τRayl, ϖRayl, τAer, ϖAer, fᵗ, qp_μ, wt_μ, Ltrunc, aerosol_optics, GreekRayleigh)
    FT = eltype(τRayl)
    @show FT
    # Output variables: Reflected and transmitted solar irradiation at TOA and BOA respectively
    R = zeros(length(vza), pol_type.n)
    T = zeros(length(vza), pol_type.n)    
    μ0 = cosd(sza)
    # @show(μ0)
    iμ0 = nearest_point(qp_μ, μ0) # input μ0 = cos(SZA)
    # @show(iμ0)
    dims = size(qp_μ)
    # @show dims
    Nquadn = pol_type.n * dims[1]
    # I0 = [1, 0, 0, 0] #assuming completely unpolarized incident stellar radiation
    D = Diagonal(repeat(pol_type.D, size(qp_μ)[1]))
    # @show D
    # get vertical grid
    # get solar+viewing geometry, compute streams
    # compute Aersol SSP
    # compute Rayleigh SSP
    Nz = length(τRayl)
    Naer = length(aerosol_optics)
    qp_μ4 = reduce(vcat, (fill.(qp_μ, [pol_type.n])))
    for m = 0:Ltrunc - 1
        # @show m

        weight = m == 0 ? 0.5 : 1.0

        # compute Zmp_Aer, Zpp_Aer, Zmp_Rayl, Zpp_Rayl
        # For m>=3, Rayleigh matrices will be 0, can catch with if statement if wanted 
        Rayl𝐙⁺⁺, Rayl𝐙⁻⁺ = PhaseFunction.compute_Z_moments(pol_type, qp_μ, GreekRayleigh, m);
        @show size(Rayl𝐙⁺⁺)
        nAer = length(aerosol_optics)
        dims = size(Rayl𝐙⁺⁺)
        
        Aer𝐙⁺⁺ = [zeros(FT, dims) for i in 1:nAer]
        Aer𝐙⁻⁺ = similar(Aer𝐙⁺⁺)

        @timeit "Aerosol Z" for i = 1:nAer
            Aer𝐙⁺⁺[i], Aer𝐙⁻⁺[i] = PhaseFunction.compute_Z_moments(pol_type, qp_μ, aerosol_optics[i].greek_coefs, m)
        end
        
        # Homogenous R and T matrices
        r⁻⁺ = zeros(FT, dims)
        t⁺⁺ = zeros(FT, dims)
        r⁺⁻ = zeros(FT, dims)
        t⁻⁻ = zeros(FT, dims)

        # Composite layer R and T matrices
        R⁻⁺ = zeros(FT, dims)
        R⁺⁻ = zeros(FT, dims)
        T⁺⁺ = zeros(FT, dims)
        T⁻⁻ = zeros(FT, dims)

        kn = 0
        # loop over vertical layers:
        for iz = 1:Nz  # Count from TOA to BOA
            @timeit "Constructing" τ, ϖ, Z⁺⁺, Z⁻⁺ = construct_atm_layer(τRayl[iz], τAer[iz,:], ϖRayl[iz], ϖAer, fᵗ, Rayl𝐙⁺⁺, Rayl𝐙⁻⁺, Aer𝐙⁺⁺, Aer𝐙⁻⁺)
            dτ_max = minimum([τ, 0.2 * minimum(qp_μ)])
            # @show dτ_max, τ, 0.2 * minimum(qp_μ)
            dτ, ndoubl = doubling_number(dτ_max, τ)
            scatter = false
            if (sum(τAer) > 1.e-8)
                scatter = true
            elseif (τRayl[iz] > 1.e-8) & (m < 3)
                scatter = true
            end        
            if (scatter)
                @timeit "elemental" rt_elemental!(pol_type, dτ, ϖ, Z⁺⁺, Z⁻⁺, m, ndoubl, scatter, qp_μ, wt_μ, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻, D)
                @timeit "doubling" rt_doubling!(dτ, τ, ndoubl, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻, D)
            else
                r⁻⁺ = 0
                r⁺⁻ = 0
                t⁺⁺ = Diagonal(exp(-τ / qp_μ4))
                t⁻⁻ = Diagonal(exp(-τ / qp_μ4))
                #= for i = 1:Nquadn
                    ii=1+floor(Int,(i-1)/pol_type.n)
                    t⁺⁺[i,i] = exp(-τ/qp_μ[ii])
                    t⁻⁻[i,i] = exp(-τ/qp_μ[ii])
                end =#
            end
            kn = get_kn(kn, scatter, iz)
            
            if (iz == 1)
                T⁺⁺[:] = t⁺⁺
                T⁻⁻[:] = t⁻⁻
                R⁻⁺[:] = r⁻⁺
                R⁺⁻[:] = r⁺⁻
            else
                @timeit "interaction" rt_interaction!(kn, R⁻⁺, T⁺⁺, R⁺⁻, T⁻⁻, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻)
            end
        end # z

        # include surface function
        # TBD
        st_iμ0 = (iμ0 - 1) * pol_type.n
        istart0 = st_iμ0 + 1
        iend0   = st_iμ0 + pol_type.n
        for i = 1:length(vza)
            iμ = nearest_point(qp_μ, cosd(vza[i])) # input vaz, vza as arrays
            # @show i, vza[i], cosd(vza[i]), iμ, qp_μ[iμ]
            # compute bigCS
            cos_m_phi = cosd(m * vaz[i])
            sin_m_phi = sind(m * vaz[i])
            if pol_type.n == 4
                bigCS = Diagonal([cos_m_phi, cos_m_phi, sin_m_phi, sin_m_phi])
            elseif pol_type.n == 3    
                bigCS = Diagonal([cos_m_phi, cos_m_phi, sin_m_phi])
            elseif pol_type.n == 1
                bigCS = Diagonal([cos_m_phi])
            end
            # accumulate Fourier moments after azimuthal weighting
            # Measurement at the TOA
            st_iμ  = (iμ - 1) * pol_type.n
            istart = st_iμ + 1
            iend   = st_iμ + pol_type.n
            # @show st_iμ+1:st_iμ+pol_type.n, iμ0,st_iμ0+1:st_iμ0+pol_type.n
            # @show size(R⁻⁺)
            
            Δ = weight * bigCS * (R⁻⁺[istart:iend, istart0:iend0] / wt_μ[iμ0]) * pol_type.I0
            # @show m, mean(abs.((Δ / R[i,:] * 100)))
            
            R[i,:] += Δ
            # @show wt_μ[iμ0]
            # Measurement at the BOA
            T[i,:] += weight * bigCS * (T⁺⁺[istart:iend, istart0:iend0] / wt_μ[iμ0]) * pol_type.I0
            # Needs something like this but working :-)
            # if mean(abs.((Δ / R[i,:] * 100))) < 0.1 # if smaller than 0.1%
            #    println("Breaking m loop at ", m, "; Max diff is now ",  mean(abs.((Δ / R[i,:] * 100))), "%")
            #    m = Ltrunc
                
            # end     
            # if m==0
            #    @show bigCS
            #    @show m, i, iμ, bigCS[1,1], weight*R⁻⁺[(iμ-1)*4+1, (iμ0-1)*4+1]/wt_μ[iμ0]   
            # end
        end
    end  # m

    
    print_timer()
    reset_timer!()

    return R, T  
end

function get_kn(kn, scatter, iz)
    if (iz == 1)
        kn = scatter ? 4 : 1
    elseif (kn >= 1)
        kn = (kn == 1) ? (!scatter ? 1 : 2) : (!scatter ? 3 : 4)
    end

    return kn
end