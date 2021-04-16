using ..Architectures: devi, default_architecture, AbstractArchitecture
using Plots

"""
    $(FUNCTIONNAME)(pol_type, obs_geom::ObsGeometry, τRayl, τAer, quadPoints, max_m, aerosol_optics, GreekRayleigh, τ_abs, brdf, architecture::AbstractArchitecture)

<< Rupesh >>

"""
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
    # Need to check this a bit better in the future!
    FT_dual = typeof(τAer[1])

    # Output variables: Reflected and transmitted solar irradiation at TOA and BOA respectively # Might need Dual later!!
    R = zeros(FT_dual, length(vza), pol_type.n, nSpec)
    T = zeros(FT_dual, length(vza), pol_type.n, nSpec)
    R_SFI = zeros(FT_dual, length(vza), pol_type.n, nSpec)
    T_SFI = zeros(FT_dual, length(vza), pol_type.n, nSpec)

    println("Processing on: ", architecture)
    println("With FT: ", FT)

    #= 
    Loop over number of truncation terms =#
    SFI = true # true #true
    @show SFI
    print("Creating arrays")
    NquadN = Nquad * pol_type.n
    dims = (NquadN,NquadN)
    
    println("Dimensions(n,n,nSpec) = (", NquadN,"," ,NquadN,"," ,nSpec,")")
    @timeit "Creating layers" added_layer         = make_added_layer(FT_dual, arr_type, dims, nSpec)
    # For surface:
    @timeit "Creating layers" added_layer_surface = make_added_layer(FT_dual, arr_type, dims, nSpec)
    # For atmosphere+surface:
    @timeit "Creating layers" composite_layer     = make_composite_layer(FT_dual, arr_type, dims, nSpec)

    # Compute aerosol Z-matrices for all aerosols
    nAer  = length(aerosol_optics)
    @timeit "Creating arrays" Aer𝐙⁺⁺ = arr_type(zeros(FT_dual, (dims[1], dims[2], nAer)))
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

        @timeit "Creating arrays" τ_sum = arr_type(zeros(FT, nSpec)) #Suniti: declaring τ_sum to be of length nSpec
        @timeit "Creating arrays" τ_λ   = arr_type(zeros(FT, nSpec))

        computed_atmosphere_properties = construct_all_atm_layers(FT, nSpec, Nz, NquadN, τRayl, τAer, aerosol_optics, Rayl𝐙⁺⁺, Rayl𝐙⁻⁺, Aer𝐙⁺⁺, Aer𝐙⁻⁺, τ_abs, arr_type, qp_μ, μ₀, m)

        # Loop over vertical layers:
        @showprogress 1 "Looping over layers ..." for iz = 1:Nz  # Count from TOA to BOA
            # Suniti: compute sum of optical thicknesses of all layers above the current layer
            # Suniti: Remember to always place the following if-else statements before the calling construct_atm_layer for the current layer!!
            if iz==1
                τ_sum  = τ_λ
            else
                τ_sum += τ_λ     
            end

            # Construct the atmospheric layer
            # From Rayleigh and aerosol τ, ϖ, compute overall layer τ, ϖ

            computed_layer_properties = get_layer_properties(computed_atmosphere_properties, iz, arr_type)

            # Whether there is scattering in the added layer, composite layer, neither or both
            scattering_interface = get_scattering_interface(scattering_interface, computed_layer_properties.scatter, iz)

            # τ * ϖ should remain constant even though they individually change over wavelength
            # @assert all(i -> (i ≈ τ * ϖ), τ_λ .* ϖ_λ)

            τ_λ = computed_layer_properties.τ_λ
            
            rt_kernel!(pol_type, SFI, added_layer, composite_layer, τ_sum, computed_layer_properties, m, quadPoints, I_static, architecture, qp_μN, scattering_interface, iz) 

            # At the bottom of the atmosphere, we have to compute total τ_sum (bottom of lowest layer), for the surface interaction later
            if iz==Nz
                τ_sum = τ_sum + τ_λ     
            end
        end 

        # Create surface matrices:
        create_surface_layer!(brdf, added_layer, SFI, m, pol_type, quadPoints, τ_sum, architecture);

        # One last interaction with surface:
        @timeit "interaction" interaction!(scattering_interface, SFI, composite_layer, added_layer, I_static)

        # Postprocess and weight according to vza
        postprocessing_vza!(iμ₀, pol_type, composite_layer, vza, qp_μ, m, vaz, μ₀, weight, nSpec, SFI, R, R_SFI, T, T_SFI)
    end

    print_timer()
    reset_timer!()

    #return R, T, R_SFI, T_SFI;
    if SFI
        return R_SFI;  
    else
        return R;
    end

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
