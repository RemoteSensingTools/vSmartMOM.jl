"""
    $(FUNCTIONNAME)(pol_type, obs_geom::ObsGeometry, τRayl, τAer, quadPoints::QuadPoints, max_m, aerosol_optics, GreekRayleigh, τ_abs, brdf, architecture::AbstractArchitecture)

Perform Radiative Transfer calculations using given parameters

"""
function rt_run(pol_type::AbstractPolarizationType,   # Polarization type (IQUV)
                obs_geom::ObsGeometry,          # Solar Zenith, Viewing Zenith, Viewing Azimuthal 
                τRayl,                          # Rayleigh optical depth 
                τAer,                           # Aerosol optical depth and single-scattering albedo
                quadPoints::QuadPoints,         # Quadrature points and weights
                max_m,                          # Max Fourier terms
                aerosol_optics,                 # AerosolOptics (greek_coefs, ω̃, k, fᵗ)
                GreekRayleigh::GreekCoefs,      # Greek coefficients of Rayleigh Phase Function
                τ_abs,                          # nSpec x Nz matrix of absorption
                brdf,                           # BRDF surface type
                architecture::AbstractArchitecture) # Whether to use CPU / GPU

    @unpack obs_alt, sza, vza, vaz = obs_geom   # Observational geometry properties
    @unpack qp_μ, wt_μ, qp_μN, wt_μN, iμ₀Nstart,μ₀, iμ₀,Nquad = quadPoints
    FT = eltype(sza)                    # Get the float-type to use
    Nz = length(τRayl)                  # Number of vertical slices
    nSpec = size(τ_abs, 1)              # Number of spectral points
    arr_type = array_type(architecture) # Type of array to use
    SFI = true                          # SFI flag
    NquadN = Nquad * pol_type.n         # Nquad (multiplied by Stokes n)
    dims = (NquadN,NquadN)              # nxn dims
    nAer  = length(aerosol_optics)      # Number of aerosols

    # Need to check this a bit better in the future!
    FT_dual = typeof(τAer[1])

    # Output variables: Reflected and transmitted solar irradiation at TOA and BOA respectively # Might need Dual later!!
    R = zeros(FT_dual, length(vza), pol_type.n, nSpec)
    T = zeros(FT_dual, length(vza), pol_type.n, nSpec)
    R_SFI = zeros(FT_dual, length(vza), pol_type.n, nSpec)
    T_SFI = zeros(FT_dual, length(vza), pol_type.n, nSpec)

    # Notify user of processing parameters
    println("Processing on:\t", architecture)
    println("With FT:\t", FT)
    @show SFI
    println("Dimensions(n,n,nSpec) = (", NquadN,"," ,NquadN,"," ,nSpec,")")
    println("Creating arrays...")

    # Create arrays
    @timeit "Creating layers" added_layer         = make_added_layer(FT_dual, arr_type, dims, nSpec)
    @timeit "Creating layers" added_layer_surface = make_added_layer(FT_dual, arr_type, dims, nSpec)
    @timeit "Creating layers" composite_layer     = make_composite_layer(FT_dual, arr_type, dims, nSpec)
    @timeit "Creating arrays" Aer𝐙⁺⁺ = arr_type(zeros(FT_dual, (dims[1], dims[2], nAer)))
    @timeit "Creating arrays" Aer𝐙⁻⁺ = similar(Aer𝐙⁺⁺)
    @timeit "Creating arrays" I_static = Diagonal(arr_type(Diagonal{FT}(ones(dims[1]))));
    
    println("Finished initializing arrays")

    # Loop over fourier moments
    for m = 0:max_m - 1

        println("Fourier Moment: ", m)

        # Azimuthal weighting
        weight = m == 0 ? FT(0.5) : FT(1.0)

        # Compute Z-moments of the Rayleigh phase matrix 
        # For m>=3, Rayleigh matrices will be 0, can catch with if statement if wanted 
        @timeit "Z moments" Rayl𝐙⁺⁺, Rayl𝐙⁻⁺ = Scattering.compute_Z_moments(pol_type, Array(qp_μ), GreekRayleigh, m, arr_type = arr_type);

        # Just for now (will change this later):
        # iBand = 1

        # Need to make sure arrays are 0:
        # TBD here
        
        # Compute aerosol Z-matrices for all aerosols
        for i = 1:nAer
            @timeit "Z moments"  Aer𝐙⁺⁺[:,:,i], Aer𝐙⁻⁺[:,:,i] = Scattering.compute_Z_moments(pol_type, Array(qp_μ), aerosol_optics[i].greek_coefs, m, arr_type = arr_type)
        end

        @timeit "Creating arrays" τ_sum_old = arr_type(zeros(FT, nSpec)) #Suniti: declaring τ_sum to be of length nSpec
        # @timeit "Creating arrays" τ_λ   = arr_type(zeros(FT, nSpec))

        # Loop over all layers and pre-compute all properties before performing core RT
        @timeit "Computing Layer Properties" computed_atmosphere_properties = construct_all_atm_layers(FT, nSpec, Nz, NquadN, τRayl, τAer, aerosol_optics, Rayl𝐙⁺⁺, Rayl𝐙⁻⁺, Aer𝐙⁺⁺, Aer𝐙⁻⁺, τ_abs, arr_type, qp_μ, μ₀, m)

        # Loop over vertical layers:
        @showprogress 1 "Looping over layers ..." for iz = 1:Nz  # Count from TOA to BOA

            # Construct the atmospheric layer
            # From Rayleigh and aerosol τ, ϖ, compute overall layer τ, ϖ
            computed_layer_properties = get_layer_properties(computed_atmosphere_properties, iz, arr_type)

            # Perform Core RT (doubling/elemental/interaction)
            rt_kernel!(pol_type, SFI, added_layer, composite_layer, computed_layer_properties, m, quadPoints, I_static, architecture, qp_μN, iz) 
        end 

        # Create surface matrices:
        # arr_type(computed_atmosphere_properties.τ_sum_all[:,end])
        create_surface_layer!(brdf, added_layer, SFI, m, pol_type, quadPoints, arr_type(computed_atmosphere_properties.τ_sum_all[:,end]), architecture);

        # One last interaction with surface:
        @timeit "interaction" interaction!(computed_atmosphere_properties.scattering_interfaces_all[end], SFI, composite_layer, added_layer, I_static)

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

"""
    $(FUNCTIONNAME)(model::vSmartMOM_Model)

Perform Radiative Transfer calculations using parameters passed in through the 
vSmartMOM_Model struct

"""
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
