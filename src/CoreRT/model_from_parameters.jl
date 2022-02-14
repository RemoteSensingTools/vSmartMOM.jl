#=

This file contains the `model_from_parameters` function, which computes all derived information
like optical thicknesses, from the input parameters. Produces a vSmartMOM_Model object. 

=#

"Generate default set of parameters for Radiative Transfer calculations (from ModelParameters/)"
default_parameters() = parameters_from_yaml(joinpath(dirname(pathof(vSmartMOM)), "CoreRT", "DefaultParameters.yaml"))

"Take the parameters specified in the vSmartMOM_Parameters struct, and calculate derived attributes into a vSmartMOM_Model" 
function model_from_parameters(params::vSmartMOM_Parameters)

    # Number of total bands and aerosols (for convenience)
    n_bands = length(params.spec_bands)
    n_aer = isnothing(params.scattering_params) ? 0 : length(params.scattering_params.rt_aerosols)

    # Create observation geometry
    obs_geom = ObsGeometry(params.sza, params.vza, params.vaz, params.obs_alt)

    # Create truncation type
    truncation_type = Scattering.δBGE{params.float_type}(params.l_trunc, params.Δ_angle)

    # Set quadrature points for streams
    quad_points = rt_set_streams(params.quadrature_type, params.l_trunc, obs_geom, params.polarization_type, array_type(params.architecture))

    # Get AtmosphericProfile from parameters
    vmr = isnothing(params.absorption_params) ? Dict() : params.absorption_params.vmr
    p_full, p_half, vmr_h2o, vcd_dry, vcd_h2o, new_vmr = compute_atmos_profile_fields(params.T, params.p, params.q, vmr)

    profile = AtmosphericProfile(params.T, p_full, params.q, p_half, vmr_h2o, vcd_dry, vcd_h2o, new_vmr)
    
    # Reduce the profile to the number of target layers (if specified)
    if params.profile_reduction_n != -1
        profile = reduce_profile(params.profile_reduction_n, profile);
    end

    # Rayleigh optical properties calculation
    greek_rayleigh = Scattering.get_greek_rayleigh(params.depol)
    τ_rayl = [zeros(params.float_type, length(params.T)) for i=1:n_bands];
    
    # This is a kludge for now, tau_abs sometimes needs to be a dual. Suniti & us need to rethink this all!!
    # i.e. code the rt core with fixed amount of derivatives as in her paper, then compute chain rule for dtau/dVMr, etc...
    FT2 = isnothing(params.absorption_params) || !haskey(params.absorption_params.vmr,"CO2") ? params.float_type : eltype(params.absorption_params.vmr["CO2"])
    τ_abs     = [zeros(FT2, length(params.spec_bands[i]), length(profile.p_full)) for i in 1:n_bands]
    
    # Loop over all bands:
    for i_band=1:n_bands

        # i'th spectral band (convert from cm⁻¹ to μm)
        curr_band_λ = 1e4 ./ params.spec_bands[i_band]

        # Compute Rayleigh properties per layer for `i_band` band center 
        τ_rayl[i_band]   = getRayleighLayerOptProp(profile.p_half[end], 
                                (maximum(curr_band_λ) + minimum(curr_band_λ)/2), 
                                params.depol, profile.vcd_dry);

        # If no absorption, continue to next band
        isnothing(params.absorption_params) && continue
        
        # Loop over all molecules in this band, obtain profile for each, and add them up
        for molec_i in 1:length(params.absorption_params.molecules[i_band])
            # This can be precomputed as well later in my mind, providing an absorption_model or an interpolation_model!
            if isempty(params.absorption_params.luts)
                # Obtain hitran data for this molecule
                @timeit "Read HITRAN"  hitran_data = read_hitran(artifact(params.absorption_params.molecules[i_band][molec_i]), iso=1)

                println("Computing profile for $(params.absorption_params.molecules[i_band][molec_i]) with vmr $(profile.vmr[params.absorption_params.molecules[i_band][molec_i]]) for band #$(i_band)")
                # Create absorption model with parameters beforehand now:
                absorption_model = make_hitran_model(hitran_data, 
                    params.absorption_params.broadening_function, 
                    wing_cutoff = params.absorption_params.wing_cutoff, 
                    CEF = params.absorption_params.CEF, 
                    architecture = params.architecture, 
                    vmr = 0);#mean(profile.vmr[params.absorption_params.molecules[i_band][molec_i]]))
                # Calculate absorption profile
                @timeit "Absorption Coeff"  compute_absorption_profile!(τ_abs[i_band], absorption_model, params.spec_bands[i_band],profile.vmr[params.absorption_params.molecules[i_band][molec_i]], profile);
            # Use LUT directly
            else
                compute_absorption_profile!(τ_abs[i_band], params.absorption_params.luts[i_band][molec_i], params.spec_bands[i_band],profile.vmr[params.absorption_params.molecules[i_band][molec_i]], profile);
            end
        end
    end

    # aerosol_optics[iBand][iAer]
    aerosol_optics = [Array{AerosolOptics}(undef, (n_aer)) for i=1:n_bands];
        
    FT2 = isnothing(params.scattering_params) ? params.float_type : typeof(params.scattering_params.rt_aerosols[1].τ_ref)
    #FT2 =  params.float_type 

    # τ_aer[iBand][iAer,iZ]
    τ_aer = [zeros(FT2, n_aer, length(profile.p_full)) for i=1:n_bands];

    # Loop over aerosol type
    for i_aer=1:n_aer

        # Get curr_aerosol
        c_aero = params.scattering_params.rt_aerosols[i_aer]
        curr_aerosol = c_aero.aerosol
        
        # Create Aerosol size distribution for each aerosol species
        size_distribution = curr_aerosol.size_distribution

        # Create a univariate aerosol distribution
        mie_aerosol = Aerosol(size_distribution, curr_aerosol.nᵣ, curr_aerosol.nᵢ)
        @show typeof(curr_aerosol.nᵣ)
        #mie_aerosol = make_mie_aerosol(size_distribution, curr_aerosol.nᵣ, curr_aerosol.nᵢ, params.scattering_params.r_max, params.scattering_params.nquad_radius) #Suniti: why is the refractive index needed here?

        # Create the aerosol extinction cross-section at the reference wavelength:
        mie_model      = make_mie_model(params.scattering_params.decomp_type, 
                                        mie_aerosol, 
                                        params.scattering_params.λ_ref, 
                                        params.polarization_type, 
                                        truncation_type, 
                                        params.scattering_params.r_max, 
                                        params.scattering_params.nquad_radius)       
        k_ref          = compute_ref_aerosol_extinction(mie_model, params.float_type)

        #params.scattering_params.rt_aerosols[i_aer].p₀, params.scattering_params.rt_aerosols[i_aer].σp
        # Loop over bands
        for i_band=1:n_bands
            
            # i'th spectral band (convert from cm⁻¹ to μm)
            curr_band_λ = 1e4 ./ params.spec_bands[i_band]

            # Create the aerosols:
            mie_model      = make_mie_model(params.scattering_params.decomp_type, 
                                            mie_aerosol, 
                                            (maximum(curr_band_λ)+minimum(curr_band_λ))/2, 
                                            params.polarization_type, 
                                            truncation_type, 
                                            params.scattering_params.r_max, 
                                            params.scattering_params.nquad_radius)

            # Compute raw (not truncated) aerosol optical properties (not needed in RT eventually) 
            # @show FT2
            @timeit "Mie calc"  aerosol_optics_raw = compute_aerosol_optical_properties(mie_model, FT2);

            # Compute truncated aerosol optical properties (phase function and fᵗ), consistent with Ltrunc:
            #@show i_aer, i_band
            aerosol_optics[i_band][i_aer] = Scattering.truncate_phase(truncation_type, 
                                                    aerosol_optics_raw; reportFit=false)

            # Compute nAer aerosol optical thickness profiles
            τ_aer[i_band][i_aer,:] = 
                params.scattering_params.rt_aerosols[i_aer].τ_ref[] * 
                (aerosol_optics[i_band][i_aer].k/k_ref) * 
                getAerosolLayerOptProp(1.0, c_aero.p₀, c_aero.σp, profile.p_full)
        end 
    end

    # Check the floating-type output matches specified FT

    # Return the model 
    return vSmartMOM_Model(params, 
                        aerosol_optics,  
                        greek_rayleigh, 
                        quad_points, 
                        τ_abs, 
                        τ_rayl, 
                        τ_aer, 
                        obs_geom, 
                        profile)
end


#=

Modified version for vibrational Ramnan scattering

=#


"Take the parameters specified in the vSmartMOM_Parameters struct, and calculate derived attributes into a vSmartMOM_Model" 
function model_from_parameters(RS_type::Union{VS_0to1_plus, VS_1to0_plus}, 
                    λ₀,
                    params::vSmartMOM_Parameters)
    @show params.absorption_params.molecules
    # Number of total bands and aerosols (for convenience)
    n_bands = 3 #length(params.spec_bands)
    n_aer = isnothing(params.scattering_params) ? 0 : length(params.scattering_params.rt_aerosols)

    # Create observation geometry
    obs_geom = ObsGeometry(params.sza, params.vza, params.vaz, params.obs_alt)

    # Create truncation type
    truncation_type = Scattering.δBGE{params.float_type}(params.l_trunc, params.Δ_angle)

    # Set quadrature points for streams
    quad_points = rt_set_streams(params.quadrature_type, params.l_trunc, obs_geom, params.polarization_type, array_type(params.architecture))

    # Get AtmosphericProfile from parameters
    vmr = isnothing(params.absorption_params) ? Dict() : params.absorption_params.vmr
    p_full, p_half, vmr_h2o, vcd_dry, vcd_h2o, new_vmr = compute_atmos_profile_fields(params.T, params.p, params.q, vmr)
    profile = AtmosphericProfile(params.T, p_full, params.q, p_half, vmr_h2o, vcd_dry, vcd_h2o, new_vmr)
    
    # Reduce the profile to the number of target layers (if specified)
    if params.profile_reduction_n != -1
        profile = reduce_profile(params.profile_reduction_n, profile);
    end

    effT = (profile.vcd_dry' * profile.T) / sum(profile.vcd_dry);
    # Define RS type
    # Compute N2 and O2
    RS_type.n2, RS_type.o2 = 
        InelasticScattering.getRamanAtmoConstants(1.e7/λ₀,effT);
    println("here 0")
    InelasticScattering.getRamanSSProp!(RS_type, λ₀);
    println("here 1")
    n_bands = length(RS_type.iBand)
    @show RS_type.grid_in
    params.spec_bands = RS_type.grid_in
    @show params.spec_bands

    # Rayleigh optical properties calculation
    greek_rayleigh = Scattering.get_greek_rayleigh(params.depol)
    τ_rayl = [zeros(params.float_type, length(params.T)) for i=1:n_bands];

    # τ_abs[iBand][iSpec,iZ]
    τ_abs     = [zeros(params.float_type, length(params.spec_bands[i]), length(profile.p_full)) for i in 1:n_bands]
    @show params.absorption_params.molecules[2]
    # Loop over all bands:
    for i_band=1:n_bands
        @show params.spec_bands[i_band]
        # i'th spectral band (convert from cm⁻¹ to μm)
        curr_band_λ = 1e4 ./ params.spec_bands[i_band]

        # Compute Rayleigh properties per layer for `i_band` band center 
        τ_rayl[i_band]   = getRayleighLayerOptProp(profile.p_half[end], 
                                (maximum(curr_band_λ) + minimum(curr_band_λ)/2), 
                                params.depol, profile.vcd_dry);

        # If no absorption, continue to next band
        isnothing(params.absorption_params) && continue
        @show i_band, params.absorption_params.molecules[i_band]
        # Loop over all molecules in this band, obtain profile for each, and add them up
        for molec_i in 1:length(params.absorption_params.molecules[i_band])
            # This can be precomputed as well later in my mind, providing an absorption_model or an interpolation_model!
            if isempty(params.absorption_params.luts)
                # Obtain hitran data for this molecule
                @timeit "Read HITRAN"  hitran_data = read_hitran(artifact(params.absorption_params.molecules[i_band][molec_i]), iso=1)

                println("Computing profile for $(params.absorption_params.molecules[i_band][molec_i]) with vmr $(profile.vmr[params.absorption_params.molecules[i_band][molec_i]]) for band #$(i_band)")
                # Create absorption model with parameters beforehand now:
                absorption_model = make_hitran_model(hitran_data, 
                    params.absorption_params.broadening_function, 
                    wing_cutoff = params.absorption_params.wing_cutoff, 
                    CEF = params.absorption_params.CEF, 
                    architecture = params.architecture, 
                    vmr = 0);#mean(profile.vmr[params.absorption_params.molecules[i_band][molec_i]]))
                # Calculate absorption profile
                @timeit "Absorption Coeff"  compute_absorption_profile!(τ_abs[i_band], absorption_model, params.spec_bands[i_band],profile.vmr[params.absorption_params.molecules[i_band][molec_i]], profile);
            # Use LUT directly
            else
                compute_absorption_profile!(τ_abs[i_band], params.absorption_params.luts[i_band][molec_i], params.spec_bands[i_band],profile.vmr[params.absorption_params.molecules[i_band][molec_i]], profile);
            end
        end
    end

    # aerosol_optics[iBand][iAer]
    aerosol_optics = [Array{AerosolOptics}(undef, (n_aer)) for i=1:n_bands];
        
    FT2 = isnothing(params.scattering_params) ? params.float_type : typeof(params.scattering_params.rt_aerosols[1].τ_ref)
    FT2 =  params.float_type 

    # τ_aer[iBand][iAer,iZ]
    τ_aer = [zeros(FT2, n_aer, length(profile.p_full)) for i=1:n_bands];

    # Loop over aerosol type
    for i_aer=1:n_aer

        # Get curr_aerosol
        c_aero = params.scattering_params.rt_aerosols[i_aer]
        curr_aerosol = c_aero.aerosol
        
        # Create Aerosol size distribution for each aerosol species
        size_distribution = curr_aerosol.size_distribution

        # Create a univariate aerosol distribution
        mie_aerosol = Aerosol(size_distribution, curr_aerosol.nᵣ, curr_aerosol.nᵢ)
        @show typeof(curr_aerosol.nᵣ)
        #mie_aerosol = make_mie_aerosol(size_distribution, curr_aerosol.nᵣ, curr_aerosol.nᵢ, params.scattering_params.r_max, params.scattering_params.nquad_radius) #Suniti: why is the refractive index needed here?

        # Create the aerosol extinction cross-section at the reference wavelength:
        mie_model      = make_mie_model(params.scattering_params.decomp_type, 
                                        mie_aerosol, 
                                        params.scattering_params.λ_ref, 
                                        params.polarization_type, 
                                        truncation_type, 
                                        params.scattering_params.r_max, 
                                        params.scattering_params.nquad_radius)       
        k_ref          = compute_ref_aerosol_extinction(mie_model, params.float_type)

        #params.scattering_params.rt_aerosols[i_aer].p₀, params.scattering_params.rt_aerosols[i_aer].σp
        # Loop over bands
        for i_band=1:n_bands
            
            # i'th spectral band (convert from cm⁻¹ to μm)
            curr_band_λ = 1e4 ./ params.spec_bands[i_band]

            # Create the aerosols:
            mie_model      = make_mie_model(params.scattering_params.decomp_type, 
                                            mie_aerosol, 
                                            (maximum(curr_band_λ)+minimum(curr_band_λ))/2, 
                                            params.polarization_type, 
                                            truncation_type, 
                                            params.scattering_params.r_max, 
                                            params.scattering_params.nquad_radius)

            # Compute raw (not truncated) aerosol optical properties (not needed in RT eventually) 
            # @show FT2
            @timeit "Mie calc"  aerosol_optics_raw = compute_aerosol_optical_properties(mie_model, FT2);

            # Compute truncated aerosol optical properties (phase function and fᵗ), consistent with Ltrunc:
            @show i_aer, i_band
            aerosol_optics[i_band][i_aer] = Scattering.truncate_phase(truncation_type, 
                                                    aerosol_optics_raw; reportFit=false)

            # Compute nAer aerosol optical thickness profiles
            τ_aer[i_band][i_aer,:] = 
                params.scattering_params.rt_aerosols[i_aer].τ_ref[] * 
                (aerosol_optics[i_band][i_aer].k/k_ref) * 
                CoreRT.getAerosolLayerOptProp(1.0, c_aero.p₀, c_aero.σp, profile.p_full)
        end 
    end

    # Check the floating-type output matches specified FT

    # Return the model 
    return vSmartMOM_Model(params, aerosol_optics,  greek_rayleigh, quad_points, τ_abs, τ_rayl, τ_aer, obs_geom, profile)

end

function loadAbsco(file; scale=(1.0))
    absco = Dataset(file)
    mol = absco["Gas_Index"][1]
    
    cs_name = "Gas_"* mol * "_Absorption"
    # Loading cross sections:
    σ = Float32(scale)*absco[cs_name][:]
    # Temperature
    T = absco["Temperature"][:]
    p = absco["Pressure"][:]/100
    ν = absco["Wavenumber"][:]
    return Absorption.AbscoTable(parse(Int,mol), -1, ν, σ, p, T )
end
