#=

This file contains the `model_from_parameters` function, which computes all derived information
like optical thicknesses, from the input parameters. Produces an RTModel object.

=#

"Generate default set of parameters for Radiative Transfer calculations (from ModelParameters/)"
default_parameters() = vSmartMOM.IO.parameters_from_yaml(joinpath(dirname(pathof(vSmartMOM)), "CoreRT", "DefaultParameters.yaml"))

"""
    model_from_parameters(params::vSmartMOM_Parameters) -> RTModel

Construct an [`RTModel`](@ref) from user-supplied [`vSmartMOM_Parameters`](@ref).

Computes all derived quantities needed by the RT solver:
- Observation geometry and quadrature points
- Rayleigh scattering coefficients (Greek/Cabannes)
- Aerosol optics via Mie theory with δ-M truncation
- Gas absorption cross-sections and optical depth profiles
- Per-band Fourier truncation limits

# Arguments
- `params::vSmartMOM_Parameters`: Configuration struct (typically from `parameters_from_yaml`).

# Returns
- `model::RTModel{ARCH, FT}`: Hierarchical model ready for `rt_run(model)`.

# See also
- `model_from_parameters(LinMode(), params)` for the linearized (Jacobian) variant.
- `parameters_from_yaml(path)` to load parameters from a YAML file.
"""
function model_from_parameters(params::vSmartMOM_Parameters)
    FT = params.float_type
    #@show FT
    # Number of total bands and aerosols (for convenience)
    n_bands = length(params.spec_bands)
    n_aer = isnothing(params.scattering_params) ? 0 : length(params.scattering_params.rt_aerosols)

    # Create observation geometry
    obs_geom = ObsGeometry{FT}(params.sza, params.vza, params.vaz, params.obs_alt)

    # Create truncation type
    truncation_type = Scattering.δBGE{FT}(params.l_trunc, params.Δ_angle)
    #@show truncation_type
    # Set quadrature points for streams
    quad_points = rt_set_streams(params.quadrature_type, params.l_trunc, obs_geom, params.polarization_type, array_type(params.architecture))

    # Get AtmosphericProfile from parameters
    vmr = isnothing(params.absorption_params) ? Dict() : params.absorption_params.vmr
    p_full, p_half, vmr_h2o, vcd_dry, vcd_h2o, new_vmr, Δz = compute_atmos_profile_fields(params.T, params.p, params.q, vmr)

    profile = AtmosphericProfile(params.T, p_full, params.q, p_half, vmr_h2o, vcd_dry, vcd_h2o, new_vmr,Δz)
    
    # Reduce the profile to the number of target layers (if specified)
    if params.profile_reduction_n != -1
        profile = reduce_profile(params.profile_reduction_n, profile);
    end

    # Rayleigh optical properties calculation
    greek_rayleigh = Scattering.get_greek_rayleigh(FT(params.depol))
    # Rayleigh optical depth per spectral point per layer (uses reduced profile size)
    τ_rayl = [zeros(FT,length(params.spec_bands[i]), length(profile.p_full)) for i=1:n_bands];
    
    # Per-band Cabannes / Rayleigh depolarization (for inelastic scattering support)
    greek_cabannes = typeof(greek_rayleigh)[]
    ϖ_Cabannes = zeros(FT, n_bands)
    
    τ_abs     = [zeros(FT, length(params.spec_bands[i]), length(profile.p_full)) for i in 1:n_bands]
    
    # Track per-band l_max from aerosol greek coef lengths
    l_max_aer = zeros(Int, max(n_aer, 1), n_bands)
    
    # Loop over all bands:
    for i_band=1:n_bands

        # i'th spectral band (convert from cm⁻¹ to μm)
        curr_band_λ = FT.(1e4 ./ params.spec_bands[i_band])
        
        # Compute per-band Cabannes properties for inelastic scattering support
        νₘ = FT(0.5) * (params.spec_bands[i_band][1] + params.spec_bands[i_band][end])
        λₘ = FT(1.0e7) / νₘ
        ϖ_Cab, γ_air_Cab, _ = InelasticScattering.compute_γ_air_Rayleigh!(λₘ)
        ϖ_Cabannes[i_band] = FT(ϖ_Cab)
        depol_air_Cab = 2γ_air_Cab / (1 + γ_air_Cab)
        push!(greek_cabannes, Scattering.get_greek_rayleigh(FT(depol_air_Cab)))
        
        # Compute Rayleigh properties per layer for `i_band` band center  
        τ_rayl[i_band]   .= getRayleighLayerOptProp(profile.p_half[end], 
                                curr_band_λ, #(mean(curr_band_λ)), 
                                params.depol, profile.vcd_dry);
        #@show τ_rayl[i_band]
        # If no absorption, continue to next band
        isnothing(params.absorption_params) && continue
        
        # Loop over all molecules in this band, obtain profile for each, and add them up
        for molec_i in 1:length(params.absorption_params.molecules[i_band])
            #@show params.absorption_params.molecules[i_band][molec_i]
            # This can be precomputed as well later in my mind, providing an absorption_model or an interpolation_model!
            if isempty(params.absorption_params.luts)
                # Obtain hitran data for this molecule
                @timeit "Read HITRAN"  hitran_data = read_hitran(artifact(params.absorption_params.molecules[i_band][molec_i]), iso=-1)

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
        
    # τ_aer[iBand][iAer,iZ]
    τ_aer = [zeros(FT, n_aer, length(profile.p_full)) for i=1:n_bands];

    # Loop over aerosol type
    for i_aer=1:n_aer

        # Get curr_aerosol
        c_aero = params.scattering_params.rt_aerosols[i_aer]
        curr_aerosol = c_aero.aerosol
        
        # Create Aerosol size distribution for each aerosol species
        size_distribution = curr_aerosol.size_distribution

        # Create a univariate aerosol distribution
        mie_aerosol = Aerosol(size_distribution, curr_aerosol.nᵣ, curr_aerosol.nᵢ)
        #@show typeof(curr_aerosol.nᵣ)
        #mie_aerosol = make_mie_aerosol(size_distribution, curr_aerosol.nᵣ, curr_aerosol.nᵢ, params.scattering_params.r_max, params.scattering_params.nquad_radius) #Suniti: why is the refractive index needed here?

        # Create the aerosol extinction cross-section at the reference wavelength:
        mie_model      = make_mie_model(params.scattering_params.decomp_type, 
                                        mie_aerosol, 
                                        params.scattering_params.λ_ref, 
                                        params.polarization_type, 
                                        truncation_type, 
                                        params.scattering_params.r_max, 
                                        params.scattering_params.nquad_radius)   
        mie_model.aerosol.nᵣ = real(params.scattering_params.n_ref)
        mie_model.aerosol.nᵢ = -imag(params.scattering_params.n_ref)
        # k for reference wavelength
        k_ref          = compute_ref_aerosol_extinction(mie_model, params.float_type)
        
        # Loop over bands
        for i_band=1:n_bands
            
            # i'th spectral band (convert from cm⁻¹ to μm)
            curr_band_λ = FT.(1e4 ./ params.spec_bands[i_band])

            # Create the aerosols:
            mie_model      = make_mie_model(params.scattering_params.decomp_type, 
                                            mie_aerosol, 
                                            (maximum(curr_band_λ)+minimum(curr_band_λ))/2, 
                                            params.polarization_type, 
                                            truncation_type, 
                                            params.scattering_params.r_max, 
                                            params.scattering_params.nquad_radius)
            n_ref = params.scattering_params.n_ref
            k = compute_ref_aerosol_extinction(mie_model,  params.float_type)
            
            #@show k
            # Compute raw (not truncated) aerosol optical properties (not needed in RT eventually) 
            @timeit "Mie calc"  aerosol_optics_raw = compute_aerosol_optical_properties(mie_model, FT);
            # Compute truncated aerosol optical properties (phase function and fᵗ), consistent with Ltrunc:
            #@show i_aer, i_band
            aerosol_optics[i_band][i_aer] = Scattering.truncate_phase(truncation_type, 
                                                    aerosol_optics_raw; reportFit=false)
            #aerosol_optics[i_band][i_aer] =  aerosol_optics_raw
            
            # Track greek coef length for l_max computation
            l_max_aer[i_aer, i_band] = min(length(aerosol_optics[i_band][i_aer].greek_coefs.β), 
                                            truncation_type.l_max)

            #@show aerosol_optics[i_band][i_aer].fᵗ
            # Compute nAer aerosol optical thickness profiles
            τ_aer[i_band][i_aer,:] = 
                params.scattering_params.rt_aerosols[i_aer].τ_ref * 
                (aerosol_optics[i_band][i_aer].k/k_ref) * 
                getAerosolLayerOptProp(1, c_aero.profile, profile)
            @info "AOD at band $i_band : $(sum(τ_aer[i_band][i_aer,:])), truncation factor = $(aerosol_optics[i_band][i_aer].fᵗ)"
        end 
    end

    # Compute per-band l_max and max_m from aerosol greek coefficient lengths
    l_max = zeros(Int, n_bands)
    max_m_bands = zeros(Int, n_bands)
    for i_band = 1:n_bands
        if n_aer > 0
            l_max[i_band] = maximum(l_max_aer[:, i_band])
        else
            l_max[i_band] = params.l_trunc
        end
        max_m_bands[i_band] = Int(ceil((l_max[i_band] + 1) / 2))
        # Clamp to user-specified max_m if it's lower (hard cutoff)
        max_m_bands[i_band] = min(max_m_bands[i_band], params.max_m)
    end

    # Build the hierarchical RTModel
    FT_ = FT
    solver = SolverConfig{FT_, typeof(params.polarization_type), typeof(params.quadrature_type)}(
        params.polarization_type,
        params.quadrature_type,
        params.max_m,
        max_m_bands,
        l_max,
        params.l_trunc,
        FT_(params.Δ_angle),
        FT_(params.depol),
    )
    atm = Atmosphere(profile, params.spec_bands)
    rayleigh = RayleighScattering(greek_rayleigh, greek_cabannes, FT_.(ϖ_Cabannes))
    aerosols_s = AerosolState(aerosol_optics, τ_aer)
    optics = Optics(rayleigh, aerosols_s, τ_abs, τ_rayl)
    return RTModel(params.architecture, solver, obs_geom, quad_points, atm, optics, params.brdf)
end


#=

Modified version for vibrational Raman scattering

=#


"Take the parameters specified in the vSmartMOM_Parameters struct, and calculate derived attributes into an RTModel"
function model_from_parameters(RS_type::Union{VS_0to1_plus, VS_1to0_plus},
                    λ₀,
                    params::vSmartMOM_Parameters)
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
    p_full, p_half, vmr_h2o, vcd_dry, vcd_h2o, new_vmr, Δz = compute_atmos_profile_fields(params.T, params.p, params.q, vmr)
    profile = AtmosphericProfile(params.T, p_full, params.q, p_half, vmr_h2o, vcd_dry, vcd_h2o, new_vmr,Δz)
    
    # Reduce the profile to the number of target layers (if specified)
    if params.profile_reduction_n != -1
        profile = reduce_profile(params.profile_reduction_n, profile);
    end

    effT = (profile.vcd_dry' * profile.T) / sum(profile.vcd_dry);
    # Define RS type
    # Compute N2 and O2
    RS_type.n2, RS_type.o2 = 
        InelasticScattering.getRamanAtmoConstants(1.e7/λ₀,effT);
    InelasticScattering.getRamanSSProp!(RS_type, λ₀);
    n_bands = length(RS_type.iBand)
    params.spec_bands = RS_type.grid_in

    # Rayleigh optical properties calculation
    greek_rayleigh = Scattering.get_greek_rayleigh(params.depol)
    τ_rayl = [zeros(params.float_type,1, length(profile.p_full)) for i=1:n_bands];

    # Per-band Cabannes / Rayleigh depolarization (for inelastic scattering support)
    FT_vrs = params.float_type
    greek_cabannes = typeof(greek_rayleigh)[]
    ϖ_Cabannes = zeros(FT_vrs, n_bands)
    l_max_aer = zeros(Int, max(n_aer, 1), n_bands)

    # Pre-allocated absorption arrays
    τ_abs     = [zeros(params.float_type, length(params.spec_bands[i]), length(profile.p_full)) for i in 1:n_bands]
    # Loop over all bands:
    for i_band=1:n_bands
        # i'th spectral band (convert from cm⁻¹ to μm)
        curr_band_λ = 1e4 ./ params.spec_bands[i_band]

        # Compute per-band Cabannes properties
        νₘ = 0.5 * (params.spec_bands[i_band][1] + params.spec_bands[i_band][end])
        λₘ = 1.0e7 / νₘ
        ϖ_Cab, γ_air_Cab, _ = InelasticScattering.compute_γ_air_Rayleigh!(λₘ)
        ϖ_Cabannes[i_band] = FT_vrs(ϖ_Cab)
        depol_air_Cab = 2γ_air_Cab / (1 + γ_air_Cab)
        push!(greek_cabannes, Scattering.get_greek_rayleigh(FT_vrs(depol_air_Cab)))

        # Compute Rayleigh properties per layer for `i_band` band center
        τ_rayl[i_band]   .= getRayleighLayerOptProp(profile.p_half[end], 
                            (maximum(curr_band_λ) + minimum(curr_band_λ))/2, 
                            params.depol, profile.vcd_dry);

        # If no absorption, continue to next band
        isnothing(params.absorption_params) && continue
        # Loop over all molecules in this band, obtain profile for each, and add them up
        for molec_i in 1:length(params.absorption_params.molecules[i_band])
            # This can be precomputed as well later in my mind, providing an absorption_model or an interpolation_model!
            if isempty(params.absorption_params.luts)
                # Obtain hitran data for this molecule
                @timeit "Read HITRAN"  hitran_data = read_hitran(artifact(params.absorption_params.molecules[i_band][molec_i]), iso=-1)

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
        
    # τ_aer[iBand][iAer,iZ]
    τ_aer = [zeros(FT, n_aer, length(profile.p_full)) for i=1:n_bands];

    # Loop over aerosol type
    for i_aer=1:n_aer

        # Get curr_aerosol
        c_aero = params.scattering_params.rt_aerosols[i_aer]
        curr_aerosol = c_aero.aerosol
        
        # Create Aerosol size distribution for each aerosol species
        size_distribution = curr_aerosol.size_distribution

        # Create a univariate aerosol distribution
        mie_aerosol = Aerosol(size_distribution, curr_aerosol.nᵣ, curr_aerosol.nᵢ)
        #mie_aerosol = make_mie_aerosol(size_distribution, curr_aerosol.nᵣ, curr_aerosol.nᵢ, params.scattering_params.r_max, params.scattering_params.nquad_radius) #Suniti: why is the refractive index needed here?

        # Create the aerosol extinction cross-section at the reference wavelength:
        mie_model      = make_mie_model(params.scattering_params.decomp_type, 
                                        mie_aerosol, 
                                        params.scattering_params.λ_ref, 
                                        params.polarization_type, 
                                        truncation_type, 
                                        params.scattering_params.r_max, 
                                        params.scattering_params.nquad_radius)
        mie_model.aerosol.nᵣ = real(params.scattering_params.n_ref)
        mie_model.aerosol.nᵢ = -imag(params.scattering_params.n_ref)
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

            @timeit "Mie calc"  aerosol_optics_raw = compute_aerosol_optical_properties(mie_model, FT2);

            # Compute truncated aerosol optical properties (phase function and fᵗ), consistent with Ltrunc:
            aerosol_optics[i_band][i_aer] = Scattering.truncate_phase(truncation_type, 
                                                    aerosol_optics_raw; reportFit=false)

            # Compute nAer aerosol optical thickness profiles
            τ_aer[i_band][i_aer,:] = 
                params.scattering_params.rt_aerosols[i_aer].τ_ref * 
                (aerosol_optics[i_band][i_aer].k/k_ref) * 
                CoreRT.getAerosolLayerOptProp(1.0, c_aero.p₀, c_aero.σp, profile)
        end 
    end

    # Compute per-band l_max and max_m
    l_max = zeros(Int, n_bands)
    max_m_bands = zeros(Int, n_bands)
    for i_band = 1:n_bands
        if n_aer > 0
            l_max[i_band] = maximum(l_max_aer[:, i_band])
        else
            l_max[i_band] = params.l_trunc
        end
        max_m_bands[i_band] = Int(ceil((l_max[i_band] + 1) / 2))
        max_m_bands[i_band] = min(max_m_bands[i_band], params.max_m)
    end

    # Build the hierarchical RTModel
    FT_vrs2 = params.float_type
    solver = SolverConfig{FT_vrs2, typeof(params.polarization_type), typeof(params.quadrature_type)}(
        params.polarization_type,
        params.quadrature_type,
        params.max_m,
        max_m_bands,
        l_max,
        params.l_trunc,
        FT_vrs2(params.Δ_angle),
        FT_vrs2(params.depol),
    )
    atm = Atmosphere(profile, params.spec_bands)
    rayleigh_s = RayleighScattering(greek_rayleigh, greek_cabannes, FT_vrs2.(ϖ_Cabannes))
    aerosols_s = AerosolState(aerosol_optics, τ_aer)
    optics = Optics(rayleigh_s, aerosols_s, τ_abs, τ_rayl)
    return RTModel(params.architecture, solver, obs_geom, quad_points, atm, optics, params.brdf)
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
