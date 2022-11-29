#=

This file contains the entry point for running the RT simulation, rt_run. 

There are two implementations: one that accepts the raw parameters, and one that accepts
the model. The latter should generally be used by users. 

=#

"""
    $(FUNCTIONNAME)(RS_type, pol_type, obs_geom::ObsGeometry, τ_rayl, τ_aer, quad_points::QuadPoints, max_m, aerosol_optics, greek_rayleigh, τ_abs, brdf, architecture::AbstractArchitecture)

Perform Radiative Transfer calculations using given parameters

"""
function rt_run_bck(RS_type::AbstractRamanType, #Default - no Raman scattering (noRS)
                pol_type::AbstractPolarizationType,   # Polarization type (IQUV)
                obs_geom::ObsGeometry,                # Solar Zenith, Viewing Zenith, Viewing Azimuthal 
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
    @unpack ϖ_Cabannes, F₀ = RS_type
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
    #@show FT_dual

    # Output variables: Reflected and transmitted solar irradiation at TOA and BOA respectively # Might need Dual later!!
    R = zeros(FT_dual, length(vza), pol_type.n, nSpec)
    T = zeros(FT_dual, length(vza), pol_type.n, nSpec)
    R_SFI = zeros(FT_dual, length(vza), pol_type.n, nSpec)
    T_SFI = zeros(FT_dual, length(vza), pol_type.n, nSpec)
    ieR_SFI = zeros(FT_dual, length(vza), pol_type.n, nSpec)
    ieT_SFI = zeros(FT_dual, length(vza), pol_type.n, nSpec)
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
    @timeit "Creating layers" composite_layer     = make_composite_layer(RS_type,FT_dual, arr_type, dims, nSpec)
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
        weight = m == 0 ? FT(0.5/π) : FT(1.0/π)
        # Compute Z-moments of the Rayleigh phase matrix 
        # For m>=3, Rayleigh matrices will be 0, can catch with if statement if wanted 
        @timeit "Z moments" Rayl𝐙⁺⁺, Rayl𝐙⁻⁺ = Scattering.compute_Z_moments(pol_type, Array(qp_μ), greek_rayleigh[1], m, arr_type = arr_type);
        if !(typeof(RS_type) <: noRS)
            @timeit "Z moments" RS_type.Z⁺⁺_λ₁λ₀, RS_type.Z⁻⁺_λ₁λ₀ = Scattering.compute_Z_moments(pol_type, Array(qp_μ), RS_type.greek_raman, m, arr_type = arr_type);
            #@show size(RS_type.Z⁺⁺_λ₁λ₀), size(RS_type.Z⁻⁺_λ₁λ₀)
        end
        #=
        aa = RS_type.ϖ_Cabannes*Rayl𝐙⁺⁺+sum(RS_type.ϖ_λ₁λ₀)*RS_type.Z⁺⁺_λ₁λ₀
        bb = RS_type.ϖ_Cabannes*Rayl𝐙⁻⁺+sum(RS_type.ϖ_λ₁λ₀)*RS_type.Z⁻⁺_λ₁λ₀
        for ia=1:NquadN
            for ib=1:NquadN
                @show ia, ib, aa[ia, ib]
            end
        end
        for ia=1:NquadN
            for ib=1:NquadN
                @show ia, ib, bb[ia, ib]
            end
        end
        blabla
        =#
        # Need to make sure arrays are 0:
        # TBD here
        
        # Compute aerosol Z-matrices for all aerosols
        for i = 1:nAer
            @timeit "Z moments"  Aer𝐙⁺⁺[:,:,i], Aer𝐙⁻⁺[:,:,i] = Scattering.compute_Z_moments(pol_type, Array(qp_μ), aerosol_optics[i].greek_coefs, m, arr_type = arr_type)
        end

        #@show RS_type.ϖ_Cabannes, ϖ_Cabannes
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
            #@show RS_type.fscattRayl
            if !(typeof(RS_type) <: noRS)
                RS_type.fscattRayl = [computed_layer_properties.fscattRayl]
            end
            #@show RS_type.fscattRayl, RS_type.ϖ_Cabannes
            # Perform Core RT (doubling/elemental/interaction)
            rt_kernel!(RS_type, pol_type, SFI, added_layer, composite_layer, computed_layer_properties, m, quad_points, I_static, architecture, qp_μN, iz) 
        end 

        # Create surface matrices:
        create_surface_layer!(RS_type, brdf, added_layer_surface, 
                    SFI, m, pol_type, quad_points, 
                    arr_type(computed_atmosphere_properties.τ_sum_all[:,end]), 
                    arr_type(F₀), architecture);

        # One last interaction with surface:
        @timeit "interaction" interaction!(RS_type,
            computed_atmosphere_properties.scattering_interfaces_all[end], 
            SFI, composite_layer, added_layer_surface, I_static)
        
            #interaction_inelastic!(RS_type,computed_atmosphere_properties.scattering_interfaces_all[end], 
        #    SFI, composite_layer, added_layer_surface, I_static)
        # Postprocess and weight according to vza
        postprocessing_vza!(RS_type, iμ₀, pol_type, composite_layer, vza, qp_μ, m, vaz, μ₀, weight, nSpec, SFI, R, R_SFI, T, T_SFI, ieR_SFI, ieT_SFI)
    end

    # Show timing statistics
    print_timer()
    reset_timer!()

    # Return R_SFI or R, depending on the flag
    return SFI ? (R_SFI, T_SFI, ieR_SFI, ieT_SFI) : (R, T)
end


"""
    $(FUNCTIONNAME)(model::vSmartMOM_Model, i_band::Integer = -1)

Perform Radiative Transfer calculations using parameters passed in through the 
vSmartMOM_Model struct

"""
function rt_run_bck(model::vSmartMOM_Model; i_band::Integer = -1)

    # Number of bands total
    n_bands = length(model.params.spec_bands)

    # Check that i_band is valid
    @assert (i_band == -1 || i_band in collect(1:n_bands)) "i_band is $(i_band) but there are only $(n_bands) bands"

    # User wants a specific band
    if i_band != -1
        return rt_run_bck(noRS(),model.params.polarization_type,
                      model.obs_geom,
                      model.τ_rayl[i_band], 
                      model.τ_aer[i_band], 
                      model.quad_points,
                      model.params.max_m,
                      model.aerosol_optics[i_band],
                      model.greek_rayleigh[i_band],
                      model.τ_abs[i_band],
                      model.params.brdf[i_band],
                      model.params.architecture)

    # User doesn't specify band, but there's only one 
    elseif n_bands == 1

        return rt_run_bck(noRS(),
                      model.params.polarization_type,
                      model.obs_geom,
                      model.τ_rayl[1], 
                      model.τ_aer[1], 
                      model.quad_points,
                      model.params.max_m,
                      model.aerosol_optics[1],
                      model.greek_rayleigh[1],
                      model.τ_abs[1],
                      model.params.brdf[1],
                      model.params.architecture)

    # User wants all bands
    else

        Rs = []

        for i in 1:n_bands

            println("------------------------------")
            println("Computing R for band #$(i)")
            println("------------------------------")

            R = rt_run_bck(noRS(),
                    model.params.polarization_type,
                       model.obs_geom,
                       model.τ_rayl[i], 
                       model.τ_aer[i], 
                       model.quad_points,
                       model.params.max_m,
                       model.aerosol_optics[i],
                       model.greek_rayleigh[i],
                       model.τ_abs[i],
                       model.params.brdf[i],
                       model.params.architecture)
            push!(Rs, R);
        end

        return Rs
    end

    
end

# Mockup if no Raman type is chosen:
function rt_run(model::vSmartMOM_Model; i_band::Integer = 1)
    rt_run(noRS(), model, i_band)
end

# Just to make sure we still have it:
function rt_run_test(RS_type::AbstractRamanType, 
        model::vSmartMOM_Model, 
        iBand)
    rt_run(RS_type,model,iBand)
end

# Mockup if no Raman type is chosen:
function rt_run_ss(model::vSmartMOM_Model; i_band::Integer = 1)
    rt_run_ss(noRS(), model, i_band)
end

# Just to make sure we still have it:
function rt_run_test_ss(RS_type::AbstractRamanType, 
        model::vSmartMOM_Model, 
        iBand)
    rt_run_ss(RS_type,model,iBand)
end

# Full multiple scattering
function rt_run(RS_type::AbstractRamanType, 
                    model::vSmartMOM_Model, iBand)
    @unpack obs_alt, sza, vza, vaz = model.obs_geom   # Observational geometry properties
    @unpack qp_μ, wt_μ, qp_μN, wt_μN, iμ₀Nstart, μ₀, iμ₀, Nquad = model.quad_points # All quadrature points
    pol_type = model.params.polarization_type
    @unpack max_m = model.params
    @unpack quad_points = model

    if obs_alt != 0
        @info "Run ms as height !=0"
        return rt_run_test_ms(RS_type, model, iBand)
    end
    
    # Also to be changed!!
    brdf = model.params.brdf[iBand[1]]
    @unpack F₀ = RS_type
    if (typeof(RS_type)<:Union{RRS,RRS_plus})
        RS_type.ϖ_λ₁λ₀ .*=  (1. - model.ϖ_Cabannes[iBand[1]])/sum(RS_type.ϖ_λ₁λ₀) 
    end   
    
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
    R       = zeros(FT_dual, length(vza), pol_type.n, nSpec)
    T       = zeros(FT_dual, length(vza), pol_type.n, nSpec)
    R_SFI   = zeros(FT_dual, length(vza), pol_type.n, nSpec)
    T_SFI   = zeros(FT_dual, length(vza), pol_type.n, nSpec)
    ieR_SFI = zeros(FT_dual, length(vza), pol_type.n, nSpec)
    ieT_SFI = zeros(FT_dual, length(vza), pol_type.n, nSpec)
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
    @timeit "Creating layers" added_layer         = 
        make_added_layer(RS_type, FT_dual, arr_type, dims, nSpec)
    # Just for now, only use noRS here
    @timeit "Creating layers" added_layer_surface = 
        make_added_layer(RS_type, FT_dual, arr_type, dims, nSpec)
    @timeit "Creating layers" composite_layer     = 
        make_composite_layer(RS_type, FT_dual, arr_type, dims, nSpec)
    @timeit "Creating arrays" I_static = 
        Diagonal(arr_type(Diagonal{FT}(ones(dims[1]))));
    #TODO: if RS_type!=noRS, create ϖ_λ₁λ₀, i_λ₁λ₀, fscattRayl, Z⁺⁺_λ₁λ₀, Z⁻⁺_λ₁λ₀ (for input), and ieJ₀⁺, ieJ₀⁻, ieR⁺⁻, ieR⁻⁺, ieT⁻⁻, ieT⁺⁺, ier⁺⁻, ier⁻⁺, iet⁻⁻, iet⁺⁺ (for output)
    #getRamanSSProp(RS_type, λ, grid_in)
    
    println("Finished initializing arrays")

    # Loop over fourier moments
    for m = 0:max_m - 1

        println("Fourier Moment: ", m, "/", max_m-1)

        # Azimuthal weighting
        weight = m == 0 ? FT(0.5/π) : FT(1.0/π)
        # Set the Zλᵢλₒ interaction parameters for Raman (or nothing for noRS)
        InelasticScattering.computeRamanZλ!(RS_type, pol_type,Array(qp_μ), m, arr_type)
        # Compute the core layer optical properties:
        @timeit "OpticalProps" layer_opt_props, fScattRayleigh   = 
            constructCoreOpticalProperties(RS_type,iBand,m,model);
        #@show size(fScattRayleigh)
        #@show size(fScattRayleigh[1])
            # Determine the scattering interface definitions:
        scattering_interfaces_all, τ_sum_all = 
            extractEffectiveProps(layer_opt_props);
        
        
        # Loop over vertical layers: 
        @showprogress 1 "Looping over layers ..." for iz = 1:Nz  # Count from TOA to BOA
            
            # Construct the atmospheric layer
            # From Rayleigh and aerosol τ, ϖ, compute overall layer τ, ϖ
            # Suniti: modified to return fscattRayl as the last element of  computed_atmosphere_properties
            if !(typeof(RS_type) <: noRS)
                RS_type.fscattRayl = expandBandScalars(RS_type, fScattRayleigh[iz]) 
            end
            
            # Expand all layer optical properties to their full dimension:
            @timeit "OpticalProps" layer_opt = 
                expandOpticalProperties(layer_opt_props[iz], arr_type)
            #@show size(layer_opt.Z⁺⁺[:,:,1]), size(RS_type.Z⁺⁺_λ₁λ₀)
            #@show typeof(layer_opt.Z⁺⁺[:,:,1]), typeof(RS_type.Z⁺⁺_λ₁λ₀)
            #aa = Array(layer_opt.Z⁺⁺[:,:,1]) #Array(RS_type.ϖ_Cabannes[1]*layer_opt.Z⁺⁺[:,:,1]) .+ (sum(RS_type.ϖ_λ₁λ₀)*RS_type.Z⁺⁺_λ₁λ₀)
            #bb = Array(layer_opt.Z⁻⁺[:,:,1]) #Array(RS_type.ϖ_Cabannes[1]*layer_opt.Z⁻⁺[:,:,1]) .+ (sum(RS_type.ϖ_λ₁λ₀)*RS_type.Z⁻⁺_λ₁λ₀)
            
            #aa = Array((RS_type.ϖ_Cabannes[1]*layer_opt.Z⁺⁺[:,:,1]) .+ (sum(RS_type.ϖ_λ₁λ₀)*RS_type.Z⁺⁺_λ₁λ₀))[1,:]
            #=
            for ia=1:NquadN
                for ib=1:NquadN
                    @show ia, ib, aa[ia, ib]
                end
            end
            
            for ia=1:NquadN
                for ib=1:NquadN
                    @show ia, ib, bb[ia, ib]
                end
            end
            bbb
            =#

            # Perform Core RT (doubling/elemental/interaction)
            rt_kernel!(RS_type, pol_type, SFI, 
                        #bandSpecLim, 
                        added_layer, composite_layer, 
                        layer_opt,
                        scattering_interfaces_all[iz], 
                        τ_sum_all[:,iz], 
                        m, quad_points, 
                        I_static, 
                        model.params.architecture, 
                        qp_μN, iz) 
        end 

        # Create surface matrices:
        create_surface_layer!(RS_type, brdf, 
                            added_layer_surface, 
                            SFI, m, 
                            pol_type, 
                            quad_points, 
                            arr_type(τ_sum_all[:,end]), 
                            arr_type(F₀),
                            model.params.architecture);
        @show F₀[:,1]
        @show scattering_interfaces_all[end]
                            #@show scattering_interfaces_all[end]
        #blapl
        # One last interaction with surface:
        
        @timeit "interaction" interaction!(RS_type,
                                    #bandSpecLim,
                                    scattering_interfaces_all[end], 
                                    SFI, 
                                    composite_layer, 
                                    added_layer_surface, 
                                    I_static)
                                    
        
        # Postprocess and weight according to vza
        postprocessing_vza!(RS_type, 
                            iμ₀, pol_type, 
                            composite_layer, 
                            vza, qp_μ, m, vaz, μ₀, 
                            weight, nSpec, 
                            SFI, 
                            R, R_SFI, 
                            T, T_SFI,
                            ieR_SFI, ieT_SFI)
    end

    # Show timing statistics
    print_timer()
    reset_timer!()

    # Return R_SFI or R, depending on the flag
    return SFI ? (R_SFI, T_SFI, ieR_SFI, ieT_SFI) : (R, T)
    #return Array(added_layer.ieJ₀⁻), Array(composite_layer.ieJ₀⁻)#
end

# Single scattering only
function rt_run_ss(RS_type::AbstractRamanType, 
    model::vSmartMOM_Model, iBand)
    @unpack obs_alt, sza, vza, vaz = model.obs_geom   # Observational geometry properties
    @unpack qp_μ, wt_μ, qp_μN, wt_μN, iμ₀Nstart, μ₀, iμ₀, Nquad = model.quad_points # All quadrature points
    pol_type = model.params.polarization_type
    @unpack max_m = model.params
    @unpack quad_points = model

    if obs_alt != 0
        return rt_run_test_ms_ss(RS_type, model, iBand)
    end

    # Also to be changed!!
    brdf = model.params.brdf[iBand[1]]
    @unpack F₀ = RS_type
    if (typeof(RS_type)<:Union{RRS,RRS_plus})
        RS_type.ϖ_λ₁λ₀ *=  (1. - model.ϖ_Cabannes[iBand[1]])/sum(RS_type.ϖ_λ₁λ₀) 
    end   

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
    R       = zeros(FT_dual, length(vza), pol_type.n, nSpec)
    T       = zeros(FT_dual, length(vza), pol_type.n, nSpec)
    R_SFI   = zeros(FT_dual, length(vza), pol_type.n, nSpec)
    T_SFI   = zeros(FT_dual, length(vza), pol_type.n, nSpec)
    ieR_SFI = zeros(FT_dual, length(vza), pol_type.n, nSpec)
    ieT_SFI = zeros(FT_dual, length(vza), pol_type.n, nSpec)
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
    @timeit "Creating layers" added_layer         = 
        make_added_layer(RS_type, FT_dual, arr_type, dims, nSpec)
    # Just for now, only use noRS here
    @timeit "Creating layers" added_layer_surface = 
        make_added_layer(RS_type, FT_dual, arr_type, dims, nSpec)
    @timeit "Creating layers" composite_layer     = 
        make_composite_layer(RS_type, FT_dual, arr_type, dims, nSpec)
    @timeit "Creating arrays" I_static = 
        Diagonal(arr_type(Diagonal{FT}(ones(dims[1]))));
    #TODO: if RS_type!=noRS, create ϖ_λ₁λ₀, i_λ₁λ₀, fscattRayl, Z⁺⁺_λ₁λ₀, Z⁻⁺_λ₁λ₀ (for input), and ieJ₀⁺, ieJ₀⁻, ieR⁺⁻, ieR⁻⁺, ieT⁻⁻, ieT⁺⁺, ier⁺⁻, ier⁻⁺, iet⁻⁻, iet⁺⁺ (for output)
    #getRamanSSProp(RS_type, λ, grid_in)

    println("Finished initializing arrays")

    # Loop over fourier moments
    for m = 0:max_m - 1

        println("Fourier Moment: ", m, "/", max_m-1)

        # Azimuthal weighting
        weight = m == 0 ? FT(0.5/π) : FT(1.0/π)
        # Set the Zλᵢλₒ interaction parameters for Raman (or nothing for noRS)
        InelasticScattering.computeRamanZλ!(RS_type, pol_type, Array(qp_μ), m, arr_type)
        # Compute the core layer optical properties:
        @timeit "OpticalProps" layer_opt_props, fScattRayleigh   = 
            constructCoreOpticalProperties(RS_type,iBand,m,model);
        # Determine the scattering interface definitions:
        scattering_interfaces_all, τ_sum_all = 
            extractEffectiveProps(layer_opt_props);

        # Loop over vertical layers: 
        @showprogress 1 "Looping over layers ..." for iz = 1:Nz  # Count from TOA to BOA

            # Construct the atmospheric layer
            # From Rayleigh and aerosol τ, ϖ, compute overall layer τ, ϖ
            # Suniti: modified to return fscattRayl as the last element of  computed_atmosphere_properties
            if !(typeof(RS_type) <: noRS)
                RS_type.fscattRayl = expandBandScalars(RS_type, fScattRayleigh[iz]) 
            end

            # Expand all layer optical properties to their full dimension:
            @timeit "OpticalProps" layer_opt = 
            expandOpticalProperties(layer_opt_props[iz], arr_type)

            # Perform Core RT (doubling/elemental/interaction)
            rt_kernel_ss!(RS_type, pol_type, SFI, 
                    #bandSpecLim, 
                    added_layer, composite_layer, 
                    layer_opt,
                    scattering_interfaces_all[iz], 
                    τ_sum_all[:,iz], 
                    m, quad_points, 
                    I_static, 
                    model.params.architecture, 
                    qp_μN, iz) 
        end 

        # Create surface matrices:
        create_surface_layer!(RS_type, brdf, 
                    added_layer_surface, 
                    SFI, m, 
                    pol_type, 
                    quad_points, 
                    arr_type(τ_sum_all[:,end]), 
                    arr_type(F₀),
                    model.params.architecture);

        # One last interaction with surface:
        #@timeit "interaction" interaction!(RS_type,
        #                    #bandSpecLim,
        #                    scattering_interfaces_all[end], 
        #                    SFI, 
        #                    composite_layer, 
        #                    added_layer_surface, 
        #                    I_static)
        
        τsurf = zeros(FT,length(τ_sum_all[:,Nz+1]))
        interaction_ss!(SFI,
                    composite_layer, 
                    added_layer_surface, 
                    τ_sum_all[:,Nz+1],
                    τsurf,
                    quad_points,
                    model.params.architecture)
        #=if !(typeof(RS_type) <: noRS)
            interaction_inelastic_ss!(RS_type,
                        SFI,
                        composite_layer, 
                        added_layer_surface, 
                        τ_sum_all[:,Nz+1],
                        τsurf,
                        quad_points,
                        model.params.architecture)
        end=#
        # Postprocess and weight according to vza
        postprocessing_vza!(RS_type, 
                    iμ₀, pol_type, 
                    composite_layer, 
                    vza, qp_μ, m, vaz, μ₀, 
                    weight, nSpec, 
                    SFI, 
                    R, R_SFI, 
                    T, T_SFI,
                    ieR_SFI, ieT_SFI)
    end

    # Show timing statistics
    print_timer()
    reset_timer!()

    # Return R_SFI or R, depending on the flag
    return SFI ? (R_SFI, T_SFI, ieR_SFI, ieT_SFI) : (R, T)
end