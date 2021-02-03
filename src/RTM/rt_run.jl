using ..Architectures: devi, default_architecture, AbstractArchitecture


function rt_run(pol_type,              # Polarization type (IQUV)
                obs_geom::ObsGeometry, # Solar Zenith, Viewing Zenith, Viewing Azimuthal 
                τRayl, ϖRayl,          # Rayleigh optical depth and single-scattering albedo
                τAer, ϖAer,            # Aerosol optical depth and single-scattering albedo
                fᵗ,                    # Truncation factor
                qp_μ, wt_μ,            # Quadrature points and weights
                Ltrunc,                # Trunction length for legendre terms
                aerosol_optics,        # AerosolOptics (greek_coefs, ω̃, k, fᵗ)
                GreekRayleigh,         # Greek coefficients of Rayleigh Phase Function
                τ_abs,                 # nSpec x Nz matrix of absorption
                architecture::AbstractArchitecture) # Whether to use CPU / GPU

    println("Processing on: ", architecture)

    #= 
    Define types, variables, and static quantities =#
    
    @unpack obs_alt, sza, vza, vaz = obs_geom   # Observational geometry properties
    FT = eltype(sza)                  # Get the float-type to use
    # FT = Float32
    Nz = length(τRayl)                  # Number of vertical slices
    nSpec = size(τ_abs, 1)              # Number of spectral points
    Nquadn = pol_type.n * size(qp_μ)[1] # Number of quadrature points 
                                        # (qp_μ array size * Stokes Vector size)
    μ0 = cosd(sza)                      # μ0 defined as cos(θ); θ = sza
    iμ0 = nearest_point(qp_μ, μ0)       # Find the closest point to μ0 in qp_μ

    arr_type = array_type(architecture)
    @show FT
    # Output variables: Reflected and transmitted solar irradiation at TOA and BOA respectively
    R = (zeros(FT, length(vza), pol_type.n, nSpec))
    T = (zeros(FT, length(vza), pol_type.n, nSpec))

    # Assuming completely unpolarized incident stellar radiation
    # This should depend on pol_type right? 
    # D = arr_type(Diagonal(repeat(pol_type.D, size(qp_μ)[1])))

    # Copy qp_μ "pol_type.n" times
    qp_μN = arr_type(repeat(qp_μ, pol_type.n)) # reduce(vcat, (fill.(arr_type(qp_μ), [pol_type.n])))

    #= 
    Loop over number of truncation terms =#

    for m = 0:Ltrunc - 1

        println("Fourier Moment: ", m)

        # Azimuthal weighting
        weight = m == 0 ? FT(0.5) : FT(1.0)

        # Compute Z-moments of the Rayleigh phase matrix 
        # For m>=3, Rayleigh matrices will be 0, can catch with if statement if wanted 
        Rayl𝐙⁺⁺, Rayl𝐙⁻⁺ = PhaseFunction.compute_Z_moments(pol_type, qp_μ, GreekRayleigh, m);
        Rayl𝐙⁺⁺, Rayl𝐙⁻⁺ = (arr_type(Rayl𝐙⁺⁺), arr_type(Rayl𝐙⁻⁺))

        # Number of aerosols
        nAer = length(aerosol_optics)
        dims = size(Rayl𝐙⁺⁺)
        
        # Compute aerosol Z-matrices for all aerosols
        # Aer𝐙⁺⁺ = [zeros(FT, dims) for i in 1:nAer]
        Aer𝐙⁺⁺ = arr_type(zeros(FT, (dims[1], dims[2], nAer)))
        Aer𝐙⁻⁺ = similar(Aer𝐙⁺⁺)

        @timeit "Aerosol Z" for i = 1:nAer
            Aer𝐙⁺⁺_curr, Aer𝐙⁻⁺_curr = PhaseFunction.compute_Z_moments(pol_type, qp_μ, aerosol_optics[i].greek_coefs, m)
            Aer𝐙⁺⁺[:,:,i], Aer𝐙⁻⁺[:,:,i] = (arr_type(Aer𝐙⁺⁺_curr), arr_type(Aer𝐙⁻⁺_curr))
        end

        # Create R and T matrices for this m

        # Homogenous R and T matrices
        # @show FT
        default_matrix = arr_type(zeros(FT, tuple(dims[1], dims[2], nSpec)))

        added_layer = AddedLayer(deepcopy(default_matrix), deepcopy(default_matrix), 
        deepcopy(default_matrix), deepcopy(default_matrix))

        composite_layer = CompositeLayer(deepcopy(default_matrix), deepcopy(default_matrix), 
        deepcopy(default_matrix), deepcopy(default_matrix))

        I_static  = Diagonal{FT}(ones(dims[1]))
        I_static_ = Diagonal(arr_type(I_static));
        # I_static_ = arr_type(repeat(I_static, 1, 1))

        scattering_interface = ScatteringInterface_00()

        # Loop over vertical layers:
        @showprogress 1 "Looping over layers ..." for iz = 1:Nz  # Count from TOA to BOA

            # Construct the atmospheric layer
            # From Rayleigh and aerosol τ, ϖ, compute overall layer τ, ϖ
            @timeit "Constructing" τ_λ, ϖ_λ, τ, ϖ, Z⁺⁺, Z⁻⁺ = construct_atm_layer(τRayl[iz], τAer[iz,:], ϖRayl[iz], ϖAer, fᵗ, Rayl𝐙⁺⁺, Rayl𝐙⁻⁺, Aer𝐙⁺⁺, Aer𝐙⁻⁺, τ_abs[:,iz], arr_type)

            # τ * ϖ should remain constant even though they individually change over wavelength
            # @assert all(i -> (i ≈ τ * ϖ), τ_λ .* ϖ_λ)

            # Compute doubling number
            dτ_max = minimum([τ * ϖ, FT(0.1) * minimum(qp_μ)])
            dτ_tmp, ndoubl = doubling_number(dτ_max, τ * ϖ)
            # @show ndoubl
            # Compute dτ vector
            # Assert that dτ .* ϖ_λ are the same
            dτ = arr_type(τ_λ ./ (FT(2)^ndoubl))
            
            # Determine whether there is scattering
            scatter = (  sum(τAer) > 1.e-8 || 
                      (( τRayl[iz] > 1.e-8 ) && (m < 3))) ? 
                      true : false

            # If there is scattering, perform the elemental and doubling steps
            if (scatter)
                
                @timeit "elemental" rt_elemental!(pol_type, dτ, dτ_max, ϖ_λ, ϖ, Z⁺⁺, Z⁻⁺, m, ndoubl, scatter, qp_μ, wt_μ, added_layer,  I_static, arr_type, architecture)
                
                @timeit "doubling" rt_doubling!(pol_type, ndoubl, added_layer, I_static_, architecture)
            else
                added_layer.r⁻⁺ = 0
                added_layer.r⁺⁻ = 0
                added_layer.t⁺⁺ = Diagonal(exp(-τ / qp_μN))
                added_layer.t⁻⁻ = Diagonal(exp(-τ / qp_μN))
            end

            # kn is an index that tells whether there is scattering in the 
            # added layer, composite layer, neither or both
            scattering_interface = get_scattering_interface(scattering_interface, scatter, iz)

            @assert !any(isnan.(added_layer.t⁺⁺))
            
            # If this TOA, just copy the added layer into the composite layer
            if (iz == 1)

                composite_layer.T⁺⁺[:] = added_layer.t⁺⁺
                composite_layer.T⁻⁻[:] = added_layer.t⁻⁻
                composite_layer.R⁻⁺[:] = added_layer.r⁻⁺
                composite_layer.R⁺⁻[:] = added_layer.r⁺⁻
            
            # If this is not the TOA, perform the interaction step
            else
                @timeit "interaction" rt_interaction!(scattering_interface, composite_layer, added_layer, I_static_)
            end
        end # z

        # include surface function

        # idx of μ0 = cos(sza)
        st_iμ0 = (iμ0 - 1) * pol_type.n
        istart0 = st_iμ0 + 1
        iend0   = st_iμ0 + pol_type.n

        R⁻⁺ = Array(composite_layer.R⁻⁺)
        T⁺⁺ = Array(composite_layer.T⁺⁺)
        # Loop over all viewing zenith angles
        for i = 1:length(vza)

            # Find the nearest quadrature point idx
            iμ = nearest_point(qp_μ, cosd(vza[i])) # input vaz, vza as arrays
            
            # TODO: Write as function, make type stable:
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

            # TODO: Write as function, make type stable:
            # Accumulate Fourier moments after azimuthal weighting
            st_iμ  = (iμ - 1) * pol_type.n
            istart = st_iμ + 1
            iend   = st_iμ + pol_type.n
            
            for s = 1:nSpec
                Δ = weight * bigCS
                Δ *= (R⁻⁺[istart:iend, istart0:iend0, s] / wt_μ[iμ0])
                Δ *= pol_type.I0
                R[i,:,s] += Δ
                T[i,:,s] += weight * bigCS * (T⁺⁺[istart:iend, istart0:iend0, s] / wt_μ[iμ0]) * pol_type.I0
            end
            
        end
    end

    print_timer()
    reset_timer!()

    return R, T  
end