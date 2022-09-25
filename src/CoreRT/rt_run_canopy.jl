#=

This file contains the entry point for running the RT simulation, rt_run. 

There are two implementations: one that accepts the raw parameters, and one that accepts
the model. The latter should generally be used by users. 

=#



# Mockup if no Raman type is chosen:
function rt_run_canopy(model::vSmartMOM_Model; i_band::Integer = 1)
    rt_run_canopy(noRS(), model, i_band)
end

function rt_run_canopy(RS_type::AbstractRamanType, 
                    model::vSmartMOM_Model, iBand)
    @unpack obs_alt, sza, vza, vaz = model.obs_geom   # Observational geometry properties
    @unpack qp_μ, wt_μ, qp_μN, wt_μN, iμ₀Nstart, μ₀, iμ₀, Nquad = model.quad_points # All quadrature points
    pol_type = model.params.polarization_type
    @unpack max_m = model.params
    @unpack quad_points = model
    FT = model.params.float_type

    n_aer = isnothing(model.params.scattering_params) ? 0 : length(model.params.scattering_params.rt_aerosols)
    
    # Also to be changed if more than 1 band is used!!
    # CFRANKEN NEEDS to be changed for concatenated arrays!!
    brdf = model.params.brdf[iBand[1]]

    @unpack ϖ_Cabannes = RS_type

    #FT = eltype(sza)                   # Get the float-type to use
    Nz = length(model.profile.p_full)   # Number of vertical slices
    
    
    
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
    NquadN =  Nquad * pol_type.n         # Nquad (multiplied by Stokes n)
    dims   = (NquadN,NquadN)              # nxn dims
    
    # Need to check this a bit better in the future!
    FT_dual = n_aer > 0 ? typeof(model.τ_aer[1][1]) : FT
    #FT_dual = FT
    @show FT_dual
    # Output variables: Reflected and transmitted solar irradiation at TOA and BOA respectively # Might need Dual later!!
    #Suniti: consider adding a new dimension (iBand) to these arrays. The assignment of simulated spectra to their specific bands will take place after batch operations, thereby leaving the computational time unaffected 
    R       = zeros(FT_dual, length(vza), pol_type.n, nSpec)
    T       = zeros(FT_dual, length(vza), pol_type.n, nSpec)
    R_SFI   = zeros(FT_dual, length(vza), pol_type.n, nSpec)
    T_SFI   = zeros(FT_dual, length(vza), pol_type.n, nSpec)
    ieR_SFI = zeros(FT_dual, length(vza), pol_type.n, nSpec)
    ieT_SFI = zeros(FT_dual, length(vza), pol_type.n, nSpec)
    hdr     = zeros(FT_dual, length(vza), pol_type.n, nSpec) # for RAMI
    bhr_dw     = zeros(FT_dual, pol_type.n, nSpec) # for RAMI
    bhr_uw     = zeros(FT_dual, pol_type.n, nSpec) # for RAMI
    hdr_J₀⁻    = zeros(FT_dual, length(vza), pol_type.n, nSpec) # for RAMI
    #  bhr[i] = bhr_uw[i,:]./bhr_dw[1,:]   
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
        weight = m == 0 ? FT(0.5) : FT(1.0)
        # Set the Zλᵢλₒ interaction parameters for Raman (or nothing for noRS)
        InelasticScattering.computeRamanZλ!(RS_type, pol_type,Array(qp_μ), m, arr_type)
        # Compute the core layer optical properties:
        @timeit "OpticalProps" layer_opt_props, fScattRayleigh   = 
            constructCoreOpticalProperties(RS_type,iBand,m,model);

        if(m==0)
            𝐙⁺⁺, 𝐙⁻⁺ = CanopyOptics.compute_Z_matrices(CanopyMod, μ, LD, m)
            
            # This basically multiplies with G again, needs to be fixed later (or removed from compute_Z_matrices)
            G1 = arr_type(CanopyOptics.G(qp_μN))
    
            canopyCore = CoreRT.CoreDirectionalScatteringOpticalProperties(LAI, ϖ_canopy, arr_type(𝐙⁺⁺), arr_type(𝐙⁻⁺), G1)
            @show canopyCore.ϖ
            # Add Canopy at the bottom here:
            layer_opt_props =  [layer_opt_props; canopyCore]
        end
        
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
            rt_kernel_canopy!(RS_type, pol_type, SFI, 
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
        create_surface_layer!(brdf, 
                            added_layer_surface, 
                            SFI, m, 
                            pol_type, 
                            quad_points, 
                            arr_type(τ_sum_all[:,end]), 
                            model.params.architecture);
        
        #@show composite_layer.J₀⁺[iμ₀,1,1:3]
        # One last interaction with surface:
        @timeit "interaction" interaction!(RS_type,
                                    #bandSpecLim,
                                    scattering_interfaces_all[end], 
                                    SFI, 
                                    composite_layer, 
                                    added_layer_surface, 
                                    I_static)
       #@show composite_layer.J₀⁺[iμ₀,1,1:3]
        hdr_J₀⁻ = similar(composite_layer.J₀⁻)
        # One last interaction with surface:
        @timeit "interaction_HDRF" interaction_hdrf!(#RS_type,
                                    #bandSpecLim,
                                    #scattering_interfaces_all[end], 
                                    SFI, 
                                    composite_layer, 
                                    added_layer_surface, 
                                    m, pol_type, quad_points,
                                    hdr_J₀⁻, bhr_uw, bhr_dw)
        
        # Postprocess and weight according to vza
        @timeit "Postprocessing" postprocessing_vza!(RS_type, 
                            iμ₀, pol_type, 
                            composite_layer, 
                            vza, qp_μ, m, vaz, μ₀, 
                            weight, nSpec, 
                            SFI, 
                            R, R_SFI, 
                            T, T_SFI,
                            ieR_SFI, ieT_SFI)

        @timeit "Postprocessing" postprocessing_vza_hdrf!(RS_type, 
            iμ₀, pol_type, 
            hdr_J₀⁻, 
            vza, qp_μ, m, vaz, μ₀, 
            weight, nSpec, 
            hdr)
            
    end

    # Show timing statistics
    print_timer()
    reset_timer!()

    # Return R_SFI or R, depending on the flag
    #if RAMI
    #@show size(hdr), size(bhr_dw)
    #hdr = hdr[:,1,:] ./ bhr_dw[1,:]
    return SFI ? (R_SFI, T_SFI, ieR_SFI, ieT_SFI, hdr, bhr_uw[1,:], bhr_dw[1,:]) : (R, T)
    #else
    #return SFI ? (R_SFI, T_SFI, ieR_SFI, ieT_SFI) : (R, T)
    #end
end