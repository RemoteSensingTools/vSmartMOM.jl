function rt_run_test(RS_type::CoreRT.AbstractRamanType, 
    model, CanopyMod, LAI, LD, ϖ_canopy, iBand)
    @unpack obs_alt, sza, vza, vaz = model.obs_geom   # Observational geometry properties
    @unpack qp_μ, wt_μ, qp_μN, wt_μN, iμ₀Nstart, μ₀, iμ₀, Nquad = model.quad_points # All quadrature points
    pol_type = model.params.polarization_type
    @unpack max_m = model.params
    @unpack quad_points = model
    arr_type = array_type(model.params.architecture)
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
    FT_dual = length(model.τ_aer[1][1]) > 0 ? typeof(model.τ_aer[1][1]) : FT

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
    CoreRT.make_added_layer(RS_type, FT_dual, arr_type, dims, nSpec)
    # Just for now, only use noRS here
    @timeit "Creating layers" added_layer_surface = 
    CoreRT.make_added_layer(RS_type, FT_dual, arr_type, dims, nSpec)
    @timeit "Creating layers" composite_layer     = 
    CoreRT.make_composite_layer(RS_type, FT_dual, arr_type, dims, nSpec)
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
        vSmartMOM.InelasticScattering.computeRamanZλ!(RS_type, pol_type,Array(qp_μ), m, arr_type)
        # Compute the core layer optical properties:
        layer_opt_props, fScattRayleigh   = CoreRT.constructCoreOpticalProperties(RS_type,iBand,m,model);
        μ = Array(qp_μ)
        𝐙⁺⁺, 𝐙⁻⁺ = CanopyOptics.compute_Z_matrices(CanopyMod, μ, LD, m)
        #@show sum(𝐙⁺⁺ ), sum(layer_opt_props[1].Z⁺⁺)
        #@show sum(𝐙⁻⁺ ), sum(layer_opt_props[1].Z⁻⁺)
        # This basically multiplies with G again, needs to be fixed later (or removed from compute_Z_matrices)
        G1 = 0.5
        #𝐙⁻⁺ *= G1
        #𝐙⁺⁺ *= G1

        canopyCore = CoreRT.CoreScatteringOpticalProperties(G1*LAI, ϖ_canopy, arr_type(𝐙⁺⁺), arr_type(𝐙⁻⁺))
        @show canopyCore.ϖ
        # Add Canopy at the bottom here:
        layer_opt_props[end] += canopyCore

        # Determine the scattering interface definitions:
        scattering_interfaces_all, τ_sum_all = CoreRT.extractEffectiveProps(layer_opt_props);

    # Loop over vertical layers: 
    for iz = 1:Nz  # Count from TOA to BOA

    # Construct the atmospheric layer
    # From Rayleigh and aerosol τ, ϖ, compute overall layer τ, ϖ
    # Suniti: modified to return fscattRayl as the last element of  computed_atmosphere_properties
    if !(typeof(RS_type) <: CoreRT.noRS)
    RS_type.fscattRayl = fScattRayleigh[iz]
    end

    # Expand all layer optical properties to their full dimension:
    layer_opt = CoreRT.expandOpticalProperties(layer_opt_props[iz], arr_type)
    
    # Perform Core RT (doubling/elemental/interaction)
    CoreRT.rt_kernel!(RS_type, pol_type, SFI, 
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
    CoreRT.create_surface_layer!(brdf, 
                added_layer_surface, 
                SFI, m, 
                pol_type, 
                quad_points, 
                arr_type(τ_sum_all[:,end]), 
                model.params.architecture);

    # One last interaction with surface:
    @timeit "interaction" CoreRT.interaction!(RS_type,
                        #bandSpecLim,
                        scattering_interfaces_all[end], 
                        SFI, 
                        composite_layer, 
                        added_layer_surface, 
                        I_static)

    # Postprocess and weight according to vza
    CoreRT.postprocessing_vza!(RS_type, 
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

function convOCO(RR, ν, ils)
    earth_out = reverse(RR[:])
    λ_grid = reverse(1e4 ./ ν)
    @time interp_I = LinearInterpolation(λ_grid, earth_out);
    res = 0.001e-3;
    off = 0.5e-3
    wl = ils.ν_out[1]-off:res:ils.ν_out[end]+off;
    #@show wl[1],wl[end], λ_grid[1],λ_grid[end]
    I_wl = interp_I(wl);
    I_conv = InstrumentOperator.conv_spectra(ils, wl, I_wl)
end