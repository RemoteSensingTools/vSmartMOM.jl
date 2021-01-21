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


    #= 
    Define types, variables, and static quantities =#
    
    FT = eltype(τRayl)                  # Get the float-type to use
    Nz = length(τRayl)                  # Number of vertical slices
    nSpec = size(τ_abs, 1)              # Number of spectral points
    Nquadn = pol_type.n * size(qp_μ)[1] # Number of quadrature points 
                                        # (qp_μ array size * Stokes Vector size)
    μ0 = cosd(sza)                      # μ0 defined as cos(θ); θ = sza
    iμ0 = nearest_point(qp_μ, μ0)       # Find the closest point to μ0 in qp_μ

    # Output variables: Reflected and transmitted solar irradiation at TOA and BOA respectively
    R = zeros(length(vza), pol_type.n, nSpec)
    T = zeros(length(vza), pol_type.n, nSpec)    

    # Assuming completely unpolarized incident stellar radiation
    # This should depend on pol_type right? 
    D = Diagonal(repeat(pol_type.D, size(qp_μ)[1]))

    # Copy qp_μ "pol_type.n" times
    qp_μ4 = reduce(vcat, (fill.(qp_μ, [pol_type.n])))

    #= 
    Loop over number of truncation terms =#

    for m = 0:Ltrunc - 1

        @show m

        # Azimuthal weighting
        weight = m == 0 ? 0.5 : 1.0

        # Compute Z-moments of the Rayleigh phase matrix 
        # For m>=3, Rayleigh matrices will be 0, can catch with if statement if wanted 
        Rayl𝐙⁺⁺, Rayl𝐙⁻⁺ = PhaseFunction.compute_Z_moments(pol_type, qp_μ, GreekRayleigh, m);

        # Number of aerosols
        nAer = length(aerosol_optics)
        dims = size(Rayl𝐙⁺⁺)
        
        # Compute aerosol Z-matrices for all aerosols
        Aer𝐙⁺⁺ = [zeros(FT, dims) for i in 1:nAer]
        Aer𝐙⁻⁺ = similar(Aer𝐙⁺⁺)

        @timeit "Aerosol Z" for i = 1:nAer
            Aer𝐙⁺⁺[i], Aer𝐙⁻⁺[i] = PhaseFunction.compute_Z_moments(pol_type, qp_μ, aerosol_optics[i].greek_coefs, m)
        end

        # Create R and T matrices for this m

        # Homogenous R and T matrices

        default_matrix = zeros(FT, tuple(dims[1], dims[2], nSpec))

        added_layer = AddedLayer(copy(default_matrix), copy(default_matrix), 
                                 copy(default_matrix), copy(default_matrix))

        composite_layer = CompositeLayer(copy(default_matrix), copy(default_matrix), 
                                         copy(default_matrix), copy(default_matrix))

        I_static = Diagonal{FT}(ones(dims[1]))
        I_static_ = repeat(I_static, 1, 1, 1)

        kn = 0

        # Loop over vertical layers:
        @showprogress 1 "Looping over layers ..." for iz = 1:Nz  # Count from TOA to BOA

            # Construct the atmospheric layer
            # From Rayleigh and aerosol τ, ϖ, compute overall layer τ, ϖ
            @timeit "Constructing" τ_nSpec, ϖ_nSpec, τ, ϖ, Z⁺⁺, Z⁻⁺ = construct_atm_layer(τRayl[iz], τAer[iz,:], ϖRayl[iz], ϖAer, fᵗ, Rayl𝐙⁺⁺, Rayl𝐙⁻⁺, Aer𝐙⁺⁺, Aer𝐙⁻⁺, τ_abs[:,iz])

            # τ * ϖ should remain constant even though they individually change over wavelength
            @assert all(i -> (i ≈ τ * ϖ), τ_nSpec .* ϖ_nSpec)

            # Compute doubling number
            dτ_max = minimum([τ * ϖ, 0.02 * minimum(qp_μ)])
            dτ_tmp, ndoubl = doubling_number(dτ_max, τ * ϖ)

            dτ = τ_nSpec ./ (2^ndoubl)
            
            # Determine whether there is scattering
            scatter = (  sum(τAer) > 1.e-8 || 
                      (( τRayl[iz] > 1.e-8 ) && (m < 3))) ? 
                      true : false

            # If there is scattering, perform the elemental and doubling steps
            if (scatter)
                
                @timeit "elemental" rt_elemental!(pol_type, dτ, dτ_max, ϖ_nSpec, ϖ, Z⁺⁺, Z⁻⁺, m, ndoubl, scatter, qp_μ, wt_μ, added_layer, Array{Float64,3}(repeat(D, 1, 1, nSpec)), I_static)

                @timeit "doubling" rt_doubling!(ndoubl, added_layer, Array{Float64,3}(repeat(D, 1, 1, nSpec)), I_static_)
            else
                added_layer.r⁻⁺ = 0
                added_layer.r⁺⁻ = 0
                added_layer.t⁺⁺ = Diagonal(exp(-τ / qp_μ4))
                added_layer.t⁻⁻ = Diagonal(exp(-τ / qp_μ4))
            end

            # kn is an index that tells whether there is scattering in the 
            # added layer, composite layer, neither or both
            kn = get_kn(kn, scatter, iz)

            @assert !any(isnan.(added_layer.t⁺⁺))
            
            # If this TOA, just copy the added layer into the composite layer
            if (iz == 1)

                composite_layer.T⁺⁺[:] = added_layer.t⁺⁺
                composite_layer.T⁻⁻[:] = added_layer.t⁻⁻
                composite_layer.R⁻⁺[:] = added_layer.r⁻⁺
                composite_layer.R⁺⁻[:] = added_layer.r⁺⁻
            
            # If this is not the TOA, perform the interaction step
            else
                @timeit "interaction" rt_interaction!(kn, composite_layer, added_layer, I_static_)
            end
        end # z

        # include surface function

        # idx of μ0 = cos(sza)
        st_iμ0 = (iμ0 - 1) * pol_type.n
        istart0 = st_iμ0 + 1
        iend0   = st_iμ0 + pol_type.n

        # Loop over all viewing zenith angles
        for i = 1:length(vza)

            # Find the nearest quadrature point idx
            iμ = nearest_point(qp_μ, cosd(vza[i])) # input vaz, vza as arrays
            
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

            # Accumulate Fourier moments after azimuthal weighting

            st_iμ  = (iμ - 1) * pol_type.n
            istart = st_iμ + 1
            iend   = st_iμ + pol_type.n
            
            for s = 1:nSpec
                Δ = weight * bigCS * (composite_layer.R⁻⁺[istart:iend, istart0:iend0, s] / wt_μ[iμ0]) * pol_type.I0
                R[i,:,s] += Δ
                T[i,:,s] += weight * bigCS * (composite_layer.T⁺⁺[istart:iend, istart0:iend0, s] / wt_μ[iμ0]) * pol_type.I0
            end
            
        end
    end

    print_timer()
    reset_timer!()

    return R, T  
end