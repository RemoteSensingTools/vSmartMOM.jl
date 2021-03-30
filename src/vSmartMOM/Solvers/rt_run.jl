using ..Architectures: devi, default_architecture, AbstractArchitecture
using Plots

function rt_run(pol_type,              # Polarization type (IQUV)
                obs_geom::ObsGeometry, # Solar Zenith, Viewing Zenith, Viewing Azimuthal 
                τRayl,                 # Rayleigh optical depth 
                τAer,                  # Aerosol optical depth and single-scattering albedo
                quadPoints,            # Quadrature points and weights
                max_m,                 # Max Fourier terms
                aerosol_optics,        # AerosolOptics (greek_coefs, ω̃, k, fᵗ)
                GreekRayleigh,         # Greek coefficients of Rayleigh Phase Function
                τ_abs,                 # nSpec x Nz matrix of absorption
                brdf,                  # BRDF surface type
                architecture::AbstractArchitecture) # Whether to use CPU / GPU

    @unpack obs_alt, sza, vza, vaz = obs_geom   # Observational geometry properties
    @unpack qp_μ, wt_μ, qp_μN, wt_μN, iμ₀Nstart,μ₀, iμ₀,Nquad = quadPoints
    FT = eltype(sza)                    # Get the float-type to use
    Nz = length(τRayl)                  # Number of vertical slices
    nSpec = size(τ_abs, 1)              # Number of spectral points
    
    arr_type = array_type(architecture)

    # Output variables: Reflected and transmitted solar irradiation at TOA and BOA respectively
    R = zeros(FT, length(vza), pol_type.n, nSpec)
    T = zeros(FT, length(vza), pol_type.n, nSpec)
    R_SFI = zeros(FT, length(vza), pol_type.n, nSpec)
    T_SFI = zeros(FT, length(vza), pol_type.n, nSpec)

    println("Processing on: ", architecture)
    println("With FT: ", FT)

    #= 
    Loop over number of truncation terms =#
    SFI = true # true #true
    @show SFI
    print("Creating arrays")
    NquadN = Nquad * pol_type.n
    dims = (NquadN,NquadN)
    @timeit "Creating layers" added_layer         = make_added_layer(FT, arr_type, dims, nSpec)
    # For surface:
    @timeit "Creating layers" added_layer_surface = make_added_layer(FT, arr_type, dims, nSpec)
    # For atmosphere+surface:
    @timeit "Creating layers" composite_layer     = make_composite_layer(FT, arr_type, dims, nSpec)

    # Compute aerosol Z-matrices for all aerosols
    nAer  = length(aerosol_optics)
    @timeit "Creating arrays" Aer𝐙⁺⁺ = arr_type(zeros(FT, (dims[1], dims[2], nAer)))
    @timeit "Creating arrays" Aer𝐙⁻⁺ = similar(Aer𝐙⁺⁺)

    @timeit "Creating arrays" I_static = Diagonal(arr_type(Diagonal{FT}(ones(dims[1]))));
    println("...done")
    for m = 0:max_m - 1

        println("Fourier Moment: ", m)

        # Azimuthal weighting
        weight = m == 0 ? FT(0.5) : FT(1.0)

        # Compute Z-moments of the Rayleigh phase matrix 
        # For m>=3, Rayleigh matrices will be 0, can catch with if statement if wanted 
        @timeit "Z moments" Rayl𝐙⁺⁺, Rayl𝐙⁻⁺ = Scattering.compute_Z_moments(pol_type, Array(qp_μ), GreekRayleigh, m, arr_type = arr_type);

        # Just for now (will change this later):
        # iBand = 1

        #nAer, nBand = size(aerosol_optics)
        #@show nAer#, nBand
        
        @show size(Rayl𝐙⁺⁺)

        # Need to make sure arrays are 0:
        # TBD here
        
        for i = 1:nAer
            #@show aerosol_optics[i,1]
            @timeit "Z moments"  Aer𝐙⁺⁺[:,:,i], Aer𝐙⁻⁺[:,:,i] = Scattering.compute_Z_moments(pol_type, Array(qp_μ), aerosol_optics[i].greek_coefs, m, arr_type = arr_type)
        end

        scattering_interface = ScatteringInterface_00()

        @timeit "Creating arrays" τ_sum = arr_type(zeros(FT,nSpec)) #Suniti: declaring τ_sum to be of length nSpec
        @timeit "Creating arrays" τ_λ   = arr_type(zeros(FT,nSpec))

        # Loop over vertical layers:
        @showprogress 1 "Looping over layers ..." for iz = 1:Nz  # Count from TOA to BOA
            # Suniti: compute sum of optical thicknesses of all layers above the current layer
            # Suniti: Remember to always place the following if-else statements before the calling construct_atm_layer for the current layer!!
            if iz==1
                τ_sum  = τ_λ
            else
                τ_sum += τ_λ     
            end
            #@show(iz, Nz)
            # Construct the atmospheric layer
            # From Rayleigh and aerosol τ, ϖ, compute overall layer τ, ϖ
            #@timeit "Constructing" 
            @timeit "Constructing" τ_λ, ϖ_λ, τ, ϖ, Z⁺⁺, Z⁻⁺ = construct_atm_layer(τRayl[iz], τAer[:,iz], aerosol_optics, Rayl𝐙⁺⁺, Rayl𝐙⁻⁺, Aer𝐙⁺⁺, Aer𝐙⁻⁺, τ_abs[:,iz], arr_type)
         
            # τ * ϖ should remain constant even though they individually change over wavelength
            # @assert all(i -> (i ≈ τ * ϖ), τ_λ .* ϖ_λ)

            # Compute doubling number
            dτ_max = minimum([τ * ϖ, FT(0.01) * minimum(qp_μ)])
            dτ, ndoubl = doubling_number(dτ_max, τ * ϖ) #Suniti
            #@show(ndoubl, dτ_max, τ)
            # Compute dτ vector
            dτ_λ = arr_type(τ_λ ./ (FT(2)^ndoubl))
            expk = exp.(-dτ_λ /μ₀) #Suniti
            
            # Determine whether there is scattering
            scatter = (  sum(τAer[:,iz]) > 1.e-8 || 
                      (( τRayl[iz] > 1.e-8 ) && (m < 3))) ? 
                      true : false
            #@show(iz, scatter)
            # If there is scattering, perform the elemental and doubling steps
            if scatter
                @timeit "elemental" elemental!(pol_type, SFI, τ_sum, dτ_λ, dτ, ϖ_λ, ϖ, Z⁺⁺, Z⁻⁺, m, ndoubl, scatter, quadPoints,  added_layer,  I_static, architecture)
                @timeit "doubling"   doubling!(pol_type, SFI, expk, ndoubl, added_layer, I_static, architecture)
            else # This might not work yet on GPU!
                # If not, there is no reflectance. Assign r/t appropriately
                added_layer.r⁻⁺[:] .= 0;
                added_layer.r⁺⁻[:] .= 0;
                added_layer.J₀⁻[:] .= 0;
                temp = Array(exp.(-τ_λ./qp_μN'))
                #added_layer.t⁺⁺, added_layer.t⁻⁻ = (Diagonal(exp(-τ_λ / qp_μN)), Diagonal(exp(-τ_λ / qp_μN)))   
                for iλ = 1:length(τ_λ)
                    added_layer.t⁺⁺[:,:,iλ] = Diagonal(temp[iλ,:]);
                    added_layer.t⁻⁻[:,:,iλ] = Diagonal(temp[iλ,:]);
                end
            end

            # Whether there is scattering in the added layer, composite layer, neither or both
            scattering_interface = get_scattering_interface(scattering_interface, scatter, iz)

            # @assert !any(isnan.(added_layer.t⁺⁺))
            
            # If this TOA, just copy the added layer into the composite layer
            if (iz == 1)
                composite_layer.T⁺⁺[:], composite_layer.T⁻⁻[:] = (added_layer.t⁺⁺, added_layer.t⁻⁻)
                composite_layer.R⁻⁺[:], composite_layer.R⁺⁻[:] = (added_layer.r⁻⁺, added_layer.r⁺⁻)
                composite_layer.J₀⁺[:], composite_layer.J₀⁻[:] = (added_layer.J₀⁺, added_layer.J₀⁻ )
            
            # If this is not the TOA, perform the interaction step
            else
                @timeit "interaction" interaction!(scattering_interface, SFI, composite_layer, added_layer, I_static)
            end
            # At the bottom of the atmosphere, we have to compute total τ_sum (bottom of lowest layer), for the surface interaction later
            if iz==Nz
                τ_sum = τ_sum + τ_λ     
            end
        end 

        #surf = LambertianSurfaceScalar(0.5)
        create_surface_layer!(brdf, added_layer, SFI, m, pol_type, quadPoints, τ_sum, architecture);
        #@show added_layer.J₀⁻[:,:,1]
        @timeit "interaction" interaction!(scattering_interface, SFI, composite_layer, added_layer, I_static)

        # All of this is "postprocessing" now, can move this into a separate function:

        # idx of μ0 = cos(sza)
        st_iμ0, istart0, iend0 = get_indices(iμ₀, pol_type);

        # Convert these to Arrays (if CuArrays), so they can be accessed by index
        R⁻⁺ = Array(composite_layer.R⁻⁺)
        T⁺⁺ = Array(composite_layer.T⁺⁺)
        J₀⁺ = Array(composite_layer.J₀⁺)
        J₀⁻ = Array(composite_layer.J₀⁻)
        # Loop over all viewing zenith angles
        for i = 1:length(vza)

            # Find the nearest quadrature point idx
            iμ = nearest_point(qp_μ, cosd(vza[i]));
            st_iμ, istart, iend = get_indices(iμ, pol_type);
            
            # Compute bigCS
            cos_m_phi, sin_m_phi = (cosd(m * vaz[i]), sind(m * vaz[i]));
            bigCS = weight * Diagonal([cos_m_phi, cos_m_phi, sin_m_phi, sin_m_phi][1:pol_type.n]);

            # Accumulate Fourier moments after azimuthal weighting
            
            for s = 1:nSpec
                R[i,:,s] += bigCS * (R⁻⁺[istart:iend, istart0:iend0, s] / μ₀) * pol_type.I₀;
                T[i,:,s] += bigCS * (T⁺⁺[istart:iend, istart0:iend0, s] / μ₀) * pol_type.I₀;
                if SFI
                    R_SFI[i,:,s] += bigCS * J₀⁻[istart:iend,1, s];
                    T_SFI[i,:,s] += bigCS * J₀⁺[istart:iend,1, s];
                end
                #@show(m,R[i,1,s], R_SFI[i,1,s])
            end
            
        end
    end

    print_timer()
    reset_timer!()

    return R, T, R_SFI, T_SFI  
end

function rt_run(model::vSmartMOM_Model)

    return rt_run(model.params.polarization_type,
                  model.obs_geom::ObsGeometry,
                  model.τRayl, 
                  model.τAer, 
                  model.quadPoints,
                  model.params.max_m,
                  model.aerosol_optics,
                  model.greek_rayleigh,
                  model.τ_abs,
                  model.brdf,
                  model.params.architecture)
end
