"""
    $(FUNCTIONNAME)(RS_type, pol_type, obs_geom::ObsGeometry, τ_rayl, τ_aer, quad_points::QuadPoints, max_m, aerosol_optics, greek_rayleigh, τ_abs, brdf, architecture::AbstractArchitecture)

Perform Radiative Transfer calculations using given parameters 
for multiple sensors at different atmospheric levels: e.g. TOA (for satellite 
observations), BOA (for ground-based observations), and (a user-defined number of) 
intermediate vertical levels for airborne (airplane, balloon, tower, hilltop) 
measurements.
Sensors are labeled from TOA->BOA, and each share the same set of viewing angles.
Upward looking viewing angles are denoted by -90<VZA<0 and downward viewing angles
are denoted by 0<VZA<90.
"""
#=
function rt_run_ms(RS_type::AbstractRamanType, #Default - no Raman scattering (noRS)
                pol_type::AbstractPolarizationType,   # Polarization type (IQUV)
                obs_geom::ObsGeometry,                # Solar Zenith, Viewing Zenith, Viewing Azimuthal 
                sensor_levels::Vector{Int64},#::ObsGeometry,           # [0] (R for satellite, T for groundbased, default), 
                                                      # [Nz] (only groundbased),
                                                      # [0, (1:Nz-1)] (general case)
                τ_rayl,                               # Rayleigh optical depth 
                τ_aer,                                # Aerosol optical depth and single-scattering albedo
                quad_points::QuadPoints,              # Quadrature points and weights
                max_m,                                # Max Fourier terms
                aerosol_optics,                       # AerosolOptics (greek_coefs, ω̃, k, fᵗ)
                greek_rayleigh::GreekCoefs,           # Greek coefficients of Rayleigh Phase Function
                τ_abs,                                # nSpec x Nz matrix of absorption
                brdf,                                 # BRDF surface type
                architecture::AbstractArchitecture)   # Whether to use CPU / GPU

    @unpack obs_alt, sza, vza, vaz = obs_geom   # Observational geometry properties
    @unpack qp_μ, wt_μ, qp_μN, wt_μN, iμ₀Nstart, μ₀, iμ₀, Nquad = quad_points # All quadrature points
    @unpack ϖ_Cabannes = RS_type
    FT = eltype(sza)                    # Get the float-type to use
    Nz = length(τ_rayl)                 # Number of vertical slices
    nSpec = size(τ_abs, 1)              # Number of spectral points
    arr_type = array_type(architecture) # Type of array to use
    SFI = true                          # SFI flag
    NquadN = Nquad * pol_type.n         # Nquad (multiplied by Stokes n)
    dims = (NquadN,NquadN)              # nxn dims
    nAer  = length(aerosol_optics)      # Number of aerosols
 
    # Need to check this a bit better in the future!
    FT_dual = length(τ_aer) > 0 ? typeof(τ_aer[1]) : FT
    @show FT_dual

    # Output variables: Reflected and transmitted solar irradiation at TOA and BOA respectively # Might need Dual later!!
    #R = zeros(FT_dual, length(sensor_levels), length(vza), pol_type.n, nSpec)
    #T = zeros(FT_dual, length(sensor_levels), length(vza), pol_type.n, nSpec)
    #R_SFI 
    uwJ = zeros(FT_dual, length(sensor_levels), length(vza), pol_type.n, nSpec)
    #T_SFI 
    dwJ = zeros(FT_dual, length(sensor_levels), length(vza), pol_type.n, nSpec)
    #ieR_SFI
    uwieJ = zeros(FT_dual, length(sensor_levels), length(vza), pol_type.n, nSpec)
    #ieT_SFI
    dwieJ = zeros(FT_dual, length(sensor_levels), length(vza), pol_type.n, nSpec)
    # Notify user of processing parameters
    msg = 
    """
    Processing on: $(architecture)
    With FT: $(FT)
    Source Function Integration: $(SFI)
    Dimensions: $((NquadN, NquadN, nSpec))
    """
    @info msg

    # Create arrays
    @timeit "Creating layers" added_layer         = make_added_layer(RS_type,FT_dual, arr_type, dims, nSpec)
    # Just for now, only use noRS here
    @timeit "Creating layers" added_layer_surface = make_added_layer(RS_type,FT_dual, arr_type, dims, nSpec)
    @timeit "Creating layers" composite_layer     = 
            make_composite_layer(RS_type, FT_dual, arr_type, length(sensor_levels), dims, nSpec)
    @timeit "Creating arrays" Aer𝐙⁺⁺ = arr_type(zeros(FT_dual, (dims[1], dims[2], nAer)))
    @timeit "Creating arrays" Aer𝐙⁻⁺ = similar(Aer𝐙⁺⁺)
    @timeit "Creating arrays" I_static = Diagonal(arr_type(Diagonal{FT}(ones(dims[1]))));
    #TODO: if RS_type!=noRS, create ϖ_λ₁λ₀, i_λ₁λ₀, fscattRayl, Z⁺⁺_λ₁λ₀, Z⁻⁺_λ₁λ₀ (for input), and ieJ₀⁺, ieJ₀⁻, ieR⁺⁻, ieR⁻⁺, ieT⁻⁻, ieT⁺⁺, ier⁺⁻, ier⁻⁺, iet⁻⁻, iet⁺⁺ (for output)
    #getRamanSSProp(RS_type, λ, grid_in)
    
    println("Finished initializing arrays")

    # Loop over fourier moments
    for m = 0:max_m - 1

        println("Fourier Moment: ", m, "/", max_m-1)

        # Azimuthal weighting
        weight = m == 0 ? FT(0.5) : FT(1.0)
        # Compute Z-moments of the Rayleigh phase matrix 
        # For m>=3, Rayleigh matrices will be 0, can catch with if statement if wanted 
        @timeit "Z moments" Rayl𝐙⁺⁺, Rayl𝐙⁻⁺ = Scattering.compute_Z_moments(pol_type, Array(qp_μ), greek_rayleigh, m, arr_type = arr_type);
        if !(typeof(RS_type) <: noRS)
            @timeit "Z moments" RS_type.Z⁺⁺_λ₁λ₀, RS_type.Z⁻⁺_λ₁λ₀ = Scattering.compute_Z_moments(pol_type, Array(qp_μ), RS_type.greek_raman, m, arr_type = arr_type);
            #@show size(RS_type.Z⁺⁺_λ₁λ₀), size(RS_type.Z⁻⁺_λ₁λ₀)
        end
        # Need to make sure arrays are 0:
        # TBD here
        
        # Compute aerosol Z-matrices for all aerosols
        for i = 1:nAer
            @timeit "Z moments"  Aer𝐙⁺⁺[:,:,i], Aer𝐙⁻⁺[:,:,i] = Scattering.compute_Z_moments(pol_type, Array(qp_μ), aerosol_optics[i].greek_coefs, m, arr_type = arr_type)
        end

        @show RS_type.ϖ_Cabannes, ϖ_Cabannes
        # Loop over all layers and pre-compute all properties before performing core RT
        @timeit "Computing Layer Properties" computed_atmosphere_properties = 
                construct_all_atm_layers(FT, nSpec, Nz, NquadN, 
                                        τ_rayl, τ_aer, aerosol_optics, 
                                        Rayl𝐙⁺⁺, Rayl𝐙⁻⁺, Aer𝐙⁺⁺, Aer𝐙⁻⁺, 
                                        τ_abs, ϖ_Cabannes,
                                        arr_type, qp_μ, μ₀, m)

        # Loop over vertical layers:
        @showprogress 1 "Looping over layers ..." for iz = 1:Nz  # Count from TOA to BOA

            # Construct the atmospheric layer
            # From Rayleigh and aerosol τ, ϖ, compute overall layer τ, ϖ
            # Suniti: modified to return fscattRayl as the last element of  computed_atmosphere_properties
            # Computing Rayleigh scattering fraction, fscattRayl = τRayl*ϖRayl/τ
            computed_layer_properties = get_layer_properties(computed_atmosphere_properties, iz, arr_type)
            #@show computed_layer_properties.fscattRayl
            @show RS_type.fscattRayl
            if !(typeof(RS_type) <: noRS)
                RS_type.fscattRayl = [computed_layer_properties.fscattRayl]
            end
            #@show RS_type.fscattRayl, RS_type.ϖ_Cabannes
            # Perform Core RT (doubling/elemental/interaction)
            rt_kernel_multisensor!(RS_type, 
                sensor_levels, 
                pol_type, 
                SFI, added_layer, composite_layer, 
                computed_layer_properties, m, quad_points, I_static, architecture, qp_μN, iz) 
        end 

        # Create surface matrices:
        create_surface_layer!(brdf, added_layer_surface, SFI, m, pol_type, quad_points, arr_type(computed_atmosphere_properties.τ_sum_all[:,end]), architecture);

        # One last interaction with surface:
        for ims=1:length(sensor_levels)
            if sensor_levels[ims]==0 #include ims==Nz with ims==0
                @timeit "interaction_multisensor" interaction_top!(
                    ims, RS_type, 
                    scattering_interface, SFI, 
                    composite_layer, added_layer_surface, I_static)                
            else
                @timeit "interaction_multisensor" interaction_bot!(
                    ims, RS_type, 
                    scattering_interface, SFI, 
                    composite_layer, added_layer_surface, I_static)
            end
        end
        #@timeit "interaction" interaction_ms!(RS_type,
        #    computed_atmosphere_properties.scattering_interfaces_all[end], 
        #    SFI, composite_layer, added_layer_surface, I_static)
        
            #interaction_inelastic!(RS_type,computed_atmosphere_properties.scattering_interfaces_all[end], 
        #    SFI, composite_layer, added_layer_surface, I_static)
        # Postprocess and weight according to vza
        postprocessing_vza_ms!(RS_type, 
            iμ₀, pol_type, 
            composite_layer, 
            vza, qp_μ, m, vaz, μ₀, 
            weight, nSpec, SFI, 
            uwJ, dwJ, uwieJ, dwieJ,
            I_static)
    end

    # Show timing statistics
    print_timer()
    reset_timer!()

    # Return R_SFI or R, depending on the flag
    #return SFI ? (R_SFI, T_SFI, ieR_SFI, ieT_SFI) : (R, T)
    return uwJ, dwJ, uwieJ, dwieJ
end
=#

function rt_run_test_ms(RS_type::AbstractRamanType, 
                        sensor_levels::Vector{Int64},
                        model::vSmartMOM_Model, iBand)
    @unpack obs_alt, sza, vza, vaz = model.obs_geom   # Observational geometry properties
    @unpack qp_μ, wt_μ, qp_μN, wt_μN, iμ₀Nstart, μ₀, iμ₀, Nquad = model.quad_points # All quadrature points
    pol_type = model.params.polarization_type
    @unpack max_m = model.params
    @unpack quad_points = model

    # Also to be changed!!
    brdf = model.params.brdf[iBand[1]]
    @unpack ϖ_Cabannes = RS_type


    FT = eltype(sza)                    # Get the float-type to use
    Nz = length(model.profile.p_full)   # Number of vertical slices
    # CFRANKEN NEEDS to be changed for concatenated arrays!!


    RS_type.bandSpecLim = [] # (1:τ_abs[iB])#zeros(Int64, iBand, 2) #Suniti: how to do this?
    #Suniti: make bandSpecLim a part of RS_type (including noRS) so that it can be passed into rt_kernel and elemental/doubling/interaction and postprocessing_vza without major syntax changes
    #put this code in model_from_parameters
    nSpec = 0;
    for iB in iBand
        nSpec0 = nSpec+1;
        nSpec += size(model.τ_abs[iB], 1); # Number of spectral points
        push!(RS_type.bandSpecLim,nSpec0:nSpec);                
    end

    arr_type = array_type(model.params.architecture) # Type of array to use
    SFI = true                          # SFI flag
    NquadN = Nquad * pol_type.n         # Nquad (multiplied by Stokes n)
    dims = (NquadN,NquadN)              # nxn dims

    # Need to check this a bit better in the future!
    #FT_dual = length(model.τ_aer[1][1]) > 0 ? typeof(model.τ_aer[1][1]) : FT
    FT_dual = FT

    # Output variables: Reflected and transmitted solar irradiation at TOA and BOA respectively # Might need Dual later!!
    #Suniti: consider adding a new dimension (iBand) to these arrays. The assignment of simulated spectra to their specific bands will take place after batch operations, thereby leaving the computational time unaffected 
    uwJ   = [zeros(FT, (length(vza), pol_type.n, nSpec)) for i=1:length(sensor_levels)]
    #zeros(FT_dual, length(sensor_levels), length(vza), pol_type.n, nSpec)
    dwJ   = [zeros(FT, (length(vza), pol_type.n, nSpec)) for i=1:length(sensor_levels)]
    #zeros(FT_dual, length(sensor_levels), length(vza), pol_type.n, nSpec)
    #R_SFI   = zeros(FT_dual, length(vza), pol_type.n, nSpec)
    #T_SFI   = zeros(FT_dual, length(vza), pol_type.n, nSpec)
    uwieJ = [zeros(FT, (length(vza), pol_type.n, nSpec)) for i=1:length(sensor_levels)]
    #zeros(FT_dual, length(sensor_levels), length(vza), pol_type.n, nSpec)
    dwieJ = [zeros(FT, (length(vza), pol_type.n, nSpec)) for i=1:length(sensor_levels)]
    #zeros(FT_dual, length(sensor_levels), length(vza), pol_type.n, nSpec)
    # Notify user of processing parameters
    msg = 
    """
    Processing on: $(architecture)
    With FT: $(FT)
    Source Function Integration: $(SFI)
    Dimensions: $((length(sensor_levels), NquadN, NquadN, nSpec))
    """
    @info msg

    # Create arrays
    @timeit "Creating layers" added_layer         = 
        make_added_layer(RS_type, FT_dual, arr_type, dims, nSpec)
    # Just for now, only use noRS here
    @timeit "Creating layers" added_layer_surface = 
        make_added_layer(RS_type, FT_dual, arr_type, dims, nSpec)
    @timeit "Creating layers" composite_layer     = 
        make_composite_layer(RS_type, FT_dual, arr_type, length(sensor_levels), dims, nSpec)
    #@timeit "Creating layers" composite_layer     = 
    #make_composite_layer(RS_type, FT_dual, arr_type, dims, nSpec)
    @timeit "Creating arrays" I_static = 
        Diagonal(arr_type(Diagonal{FT}(ones(dims[1]))));
    #TODO: if RS_type!=noRS, create ϖ_λ₁λ₀, i_λ₁λ₀, fscattRayl, Z⁺⁺_λ₁λ₀, Z⁻⁺_λ₁λ₀ (for input), and ieJ₀⁺, ieJ₀⁻, ieR⁺⁻, ieR⁻⁺, ieT⁻⁻, ieT⁺⁺, ier⁺⁻, ier⁻⁺, iet⁻⁻, iet⁺⁺ (for output)
    #getRamanSSProp(RS_type, λ, grid_in)

    println("Finished initializing arrays")

    # Loop over fourier moments
    for m = 0:max_m - 1

        println("Fourier Moment: ", m, "/", max_m-1)

        # Azimuthal weighting
        weight = m == 0 ? FT(0.5) : FT(1.0)
        # Set the Zλᵢλₒ interaction parameters for Raman (or nothing for noRS)
        InelasticScattering.computeRamanZλ!(RS_type, pol_type, Array(qp_μ), m, arr_type)
        # Compute the core layer optical properties:
        @show iBand
        layer_opt_props, fScattRayleigh   = 
            constructCoreOpticalProperties(RS_type,iBand,m,model);
        # Determine the scattering interface definitions:
        scattering_interfaces_all, τ_sum_all = 
            extractEffectiveProps(layer_opt_props);

        # Loop over vertical layers: 
        @showprogress 1 "Looping over layers ..." for iz = 1:Nz  # Count from TOA to BOA

            # Construct the atmospheric layer
            # From Rayleigh and aerosol τ, ϖ, compute overall layer τ, ϖ
            # computed_layer_properties = get_layer_properties(computed_atmosphere_properties, iz, arr_type)
            
            # Suniti: modified to return fscattRayl as the last element of  computed_atmosphere_properties
            if !(typeof(RS_type) <: noRS)
                #@show fScattRayleigh[iz]
                RS_type.fscattRayl = expandBandScalars(RS_type, fScattRayleigh[iz]) ; #Array(fScattRayleigh[iz])
            end

            # Expand all layer optical properties to their full dimension:
            layer_opt = 
                expandOpticalProperties(layer_opt_props[iz], arr_type) 
                
            # Perform Core RT (doubling/elemental/interaction)
            rt_kernel_multisensor!(RS_type,
                    sensor_levels, 
                    pol_type, SFI, 
                    #bandSpecLim, 
                    added_layer, composite_layer, 
                    layer_opt,
                    scattering_interfaces_all[iz], 
                    τ_sum_all[:,iz], 
                    m, quad_points, 
                    I_static, 
                    model.params.architecture, 
                    qp_μN, iz, arr_type) 
        end 

        # Create surface matrices:
        create_surface_layer!(brdf, 
                    added_layer_surface, 
                    SFI, m, 
                    pol_type, 
                    quad_points, 
                    arr_type(τ_sum_all[:,end]), 
                    model.params.architecture);

        # One last interaction with surface:
        for ims=1:length(sensor_levels)
            #if sensor_levels[ims]==0 #include ims==Nz with ims==0
            #    @timeit "interaction_multisensor" interaction_top!(
            #        ims, RS_type, 
            #        scattering_interfaces_all[end], 
            #        SFI, 
            #        composite_layer, added_layer_surface, I_static, arr_type)                
            #else
            @timeit "interaction_multisensor" interaction_bot!(
                    ims, RS_type, 
                    scattering_interfaces_all[end], 
                    SFI, 
                    composite_layer, added_layer_surface, I_static, arr_type)
            #end
        end
        # @timeit "interaction" interaction_ms!(RS_type,
        #                    #bandSpecLim,
        #                    scattering_interfaces_all[end], 
        #                    SFI, 
        #                    composite_layer, 
        #                    added_layer_surface, 
        #                    I_static)

        # Postprocess and weight according to vza
        postprocessing_vza_ms!(RS_type, 
                    sensor_levels,
                    iμ₀, pol_type, 
                    composite_layer, 
                    vza, qp_μ, m, vaz, μ₀, 
                    weight, nSpec, SFI, 
                    uwJ, dwJ, uwieJ, dwieJ,
                    I_static, arr_type)
    end

    # Show timing statistics
    print_timer()
    reset_timer!()

    # Return R_SFI or R, depending on the flag
    return uwJ, dwJ, uwieJ, dwieJ #SFI ? (R_SFI, T_SFI, ieR_SFI, ieT_SFI) : (R, T)
end