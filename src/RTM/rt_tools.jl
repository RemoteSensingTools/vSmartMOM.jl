using ..Architectures: devi, default_architecture


function run_RTM(pol_type,          # Polarization type (IQUV)
                 sza, vza, vaz,     # Solar Zenith, Viewing Zenith, Viewing Azimuthal 
                 τRayl, ϖRayl,      # Rayleigh optical depth and single-scattering albedo
                 τAer, ϖAer,        # Aerosol optical depth and single-scattering albedo
                 fᵗ,                # Truncation factor
                 qp_μ, wt_μ,        # Quadrature points and weights
                 Ltrunc,            # Trunction length for legendre terms
                 aerosol_optics,    # AerosolOptics (greek_coefs, ω̃, k, fᵗ)
                 GreekRayleigh,     # Greek coefficients of Rayleigh Phase Function
                 τ_abs)             # nSpec x Nz matrix of absorption

    # Get the float-type to use
    FT = eltype(τRayl)

    nSpec = size(τ_abs, 1)

    # Output variables: Reflected and transmitted solar irradiation at TOA and BOA respectively
    R = zeros(length(vza), pol_type.n, nSpec)
    T = zeros(length(vza), pol_type.n, nSpec)    

    # μ0 defined as cos(θ); θ = sza
    μ0 = cosd(sza)

    # Find the closest point to μ0 in qp_μ
    iμ0 = nearest_point(qp_μ, μ0)

    # Dimensions of quadrature points array
    dims = size(qp_μ)

    # Number of quadrature points (qp_μ array size * Stokes Vector size)
    Nquadn = pol_type.n * dims[1]

    

    # I0 = [1, 0, 0, 0] 
    # assuming completely unpolarized incident stellar radiation
    # This should depend on pol_type right? 
    D = Diagonal(repeat(pol_type.D, size(qp_μ)[1]))

    # Number of vertical slices
    Nz = length(τRayl)

    # Copy qp_μ "pol_type.n" times
    qp_μ4 = reduce(vcat, (fill.(qp_μ, [pol_type.n])))

    # get vertical grid
    # get solar+viewing geometry, compute streams
    # compute Aersol SSP
    # compute Rayleigh SSP

    # Loop over number of truncation terms
    for m = 0:Ltrunc - 1

        # Azimuthal weighting
        weight = m == 0 ? 0.5 : 1.0

        # Compute Z-moments of the Rayleigh phase matrix 
        # For m>=3, Rayleigh matrices will be 0, can catch with if statement if wanted 
        Rayl𝐙⁺⁺, Rayl𝐙⁻⁺ = PhaseFunction.compute_Z_moments(pol_type, qp_μ, GreekRayleigh, m);

        # Number of aerosols
        nAer = length(aerosol_optics)
        dims = size(Rayl𝐙⁺⁺)
        
        # Compute aerosol Z-matrices
        Aer𝐙⁺⁺ = [zeros(FT, dims) for i in 1:nAer]
        Aer𝐙⁻⁺ = similar(Aer𝐙⁺⁺)

        @timeit "Aerosol Z" for i = 1:nAer
            Aer𝐙⁺⁺[i], Aer𝐙⁻⁺[i] = PhaseFunction.compute_Z_moments(pol_type, qp_μ, aerosol_optics[i].greek_coefs, m)
        end
        
        # Note: The following are n x n x 1 now, but need to be n x n x nSpec

        # Homogenous R and T matrices
        r⁻⁺ = zeros(FT, tuple(dims[1], dims[2], nSpec))
        t⁺⁺ = zeros(FT, tuple(dims[1], dims[2], nSpec))
        r⁺⁻ = zeros(FT, tuple(dims[1], dims[2], nSpec))
        t⁻⁻ = zeros(FT, tuple(dims[1], dims[2], nSpec))

        # Composite layer R and T matrices
        R⁻⁺ = zeros(FT, tuple(dims[1], dims[2], nSpec))
        R⁺⁻ = zeros(FT, tuple(dims[1], dims[2], nSpec))
        T⁺⁺ = zeros(FT, tuple(dims[1], dims[2], nSpec))
        T⁻⁻ = zeros(FT, tuple(dims[1], dims[2], nSpec))

        I_static = Diagonal{FT}(ones(dims[1]))
        I_static_ = repeat(I_static, 1, 1, 1)

        kn = 0

        # Loop over vertical layers:
        for iz = 1:Nz  # Count from TOA to BOA

            @show iz

            # Construct the atmospheric layer
            # From Rayleigh and aerosol τ, ϖ, compute overall layer τ, ϖ
            @timeit "Constructing" τ_nSpec, ϖ_nSpec, τ, ϖ, Z⁺⁺, Z⁻⁺ = construct_atm_layer(τRayl[iz], τAer[iz,:], ϖRayl[iz], ϖAer, fᵗ, Rayl𝐙⁺⁺, Rayl𝐙⁻⁺, Aer𝐙⁺⁺, Aer𝐙⁻⁺, τ_abs[:,iz])

            # τ * ϖ should remain constant even though they individually change over wavelength
            @assert all(i->(i==τ*ϖ), τ_nSpec .* ϖ_nSpec)

            dτ_max = minimum([τ * ϖ, 0.2 * minimum(qp_μ)])
            dτ_tmp, ndoubl = doubling_number(dτ_max, τ*ϖ)

            dτ = τ_nSpec ./ (2^ndoubl)
            
            scatter = false
            if (sum(τAer) > 1.e-8)
                scatter = true

            elseif (τRayl[iz] > 1.e-8) && (m < 3)
                scatter = true
            end      

            if (scatter)
                # @timeit "elemental" rt_elemental!(pol_type, dτ, ϖ, Z⁺⁺, Z⁻⁺, m, ndoubl, scatter, qp_μ, wt_μ, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻, D)
                
                @timeit "elemental" rt_elemental!(pol_type, dτ, dτ_max, ϖ_nSpec, ϖ, Z⁺⁺, Z⁻⁺, m, ndoubl, scatter, qp_μ, wt_μ, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻, Array{Float64,3}(repeat(D, 1, 1, nSpec)), I_static)

                @timeit "doubling" rt_doubling!(ndoubl, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻, Array{Float64,3}(repeat(D, 1, 1, nSpec)), I_static_)
                # @timeit "doubling" rt_doubling!(dτ, τ, ndoubl, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻, D)
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
                
                @timeit "interaction" rt_interaction!(kn, R⁻⁺, T⁺⁺, R⁺⁻, T⁻⁻, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻, I_static_)
                
                # @timeit "interaction" rt_interaction!(kn, R⁻⁺, T⁺⁺, R⁺⁻, T⁻⁻, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻)
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
            
            for s = 1:nSpec
                Δ = weight * bigCS * (R⁻⁺[istart:iend, istart0:iend0, s] / wt_μ[iμ0]) * pol_type.I0
                R[i,:,s] += Δ
                T[i,:,s] += weight * bigCS * (T⁺⁺[istart:iend, istart0:iend0, s] / wt_μ[iμ0]) * pol_type.I0
            end

            

            
            # @show m, mean(abs.((Δ / R[i,:] * 100)))
            
            
            # @show wt_μ[iμ0]
            # Measurement at the BOA
            
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