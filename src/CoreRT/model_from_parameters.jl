#=

This file contains the `model_from_parameters` function, which computes all derived information
like optical thicknesses, from the input parameters. Produces a vSmartMOM_Model object. 

=#

"Generate default set of parameters for Radiative Transfer calculations (from ModelParameters/)"
default_parameters() = parameters_from_yaml(joinpath(dirname(pathof(vSmartMOM)), "CoreRT", "DefaultParameters.yaml"))

"""
Given observational altitudes [km], pressure [hPa] and temperature [K] profiles, 
return new profiles with interpolated layers added and new layer indices for sensors
"""
function resize_layers_for_ms(obs_alt, p, T, q)
    
    if obs_alt==0.0
        return obs_alt, p, T, q;
    else
        p₀ = 1013.25 # hPa
        H₀ = 8       # km

        # Convert from km to hPa
        obs_p = p₀ * exp.(-obs_alt/H₀)

        p_fixed = copy(0.5*(p[2:end].+p[1:(length(p)-1)]))

        new_p = copy(0.5*(p[2:end].+p[1:(length(p)-1)]))
        new_T = copy(T)
        new_q = copy(q)

        sensor_idx = []

        filter_out_idx = (obs_alt .<= 500 .&& obs_alt .!= 0)
        obs_p = obs_p[filter_out_idx]

        if !isempty(filter(x-> x > 500 || x ≈ 0, obs_alt))
            push!(sensor_idx, 0)
        end
        obs_T_int = LinearInterpolation(p_fixed, T)
        obs_q_int = LinearInterpolation(p_fixed, q)
        # obs_T = zeros(length(obs_p))
        for i=1:length(obs_p)
            obs_T = obs_T_int(obs_p[i]);
            obs_q = obs_q_int(obs_p[i]);
            sort!(push!(new_p, obs_p[i]))

            curr_idx = indexin(obs_p[i], new_p)[1]
            insert!(new_T, curr_idx, obs_T)
            insert!(new_q, curr_idx, obs_q)

            push!(sensor_idx, curr_idx)
        end 

        push!(new_p, new_p[end])
        @show "NOTE: T and q are currently half-levels, but they should be at the boundary"

        return convert(Vector{Int64}, sensor_idx), new_p, new_T, new_q
    end
end

"Set a common l_max in a given band for each aerosol type"
function set_uniform_lmax!(lmax::Vector{Int}, aerosol_optics::Vector{Vector{vSmartMOM.Scattering.AerosolOptics}})
    n_bands = length(lmax)
    n_aer = length(aerosol_optics[1])
    for i_band = 1:n_bands
        pref_lmax = lmax[i_band]
        for i_aer = 1:n_aer
            aer_lmax = length(aerosol_optics[i_band][i_aer].greek_coefs.β)
            if aer_lmax < pref_lmax
                α = aerosol_optics[i_band][i_aer].greek_coefs.α
                β = aerosol_optics[i_band][i_aer].greek_coefs.β
                γ = aerosol_optics[i_band][i_aer].greek_coefs.γ
                δ = aerosol_optics[i_band][i_aer].greek_coefs.δ
                ϵ = aerosol_optics[i_band][i_aer].greek_coefs.ϵ
                ζ = aerosol_optics[i_band][i_aer].greek_coefs.ζ

                aerosol_optics[i_band][i_aer].greek_coefs.α = zeros(eltype(α), pref_lmax)
                aerosol_optics[i_band][i_aer].greek_coefs.α[1:aer_lmax] .= α
                
                aerosol_optics[i_band][i_aer].greek_coefs.β = zeros(eltype(β), pref_lmax)
                aerosol_optics[i_band][i_aer].greek_coefs.β[1:aer_lmax] .= β

                aerosol_optics[i_band][i_aer].greek_coefs.γ = zeros(eltype(γ), pref_lmax)
                aerosol_optics[i_band][i_aer].greek_coefs.γ[1:aer_lmax] .= γ

                aerosol_optics[i_band][i_aer].greek_coefs.δ = zeros(eltype(δ), pref_lmax)
                aerosol_optics[i_band][i_aer].greek_coefs.δ[1:aer_lmax] .= δ

                aerosol_optics[i_band][i_aer].greek_coefs.ϵ = zeros(eltype(ϵ), pref_lmax)
                aerosol_optics[i_band][i_aer].greek_coefs.ϵ[1:aer_lmax] .= ϵ

                aerosol_optics[i_band][i_aer].greek_coefs.ζ = zeros(eltype(ζ), pref_lmax)
                aerosol_optics[i_band][i_aer].greek_coefs.ζ[1:aer_lmax] .= ζ
            end
        end
    end
end

"Take the parameters specified in the vSmartMOM_Parameters struct, and calculate derived attributes into a vSmartMOM_Model" 
function model_from_parameters(params::vSmartMOM_Parameters)

    # Number of total bands and aerosols (for convenience)
    n_bands = length(params.spec_bands)
    n_aer = isnothing(params.scattering_params) ? 0 : length(params.scattering_params.rt_aerosols)

    
    # Get new p/T profiles using obs_alt
    #Suniti TODO
    #p, T, q = resize_layers(params.obs_alt, params.p, params.T, params.q)
    # sensor_levels, p, T, q = resize_layers_for_ms(params.obs_alt, params.p, params.T, params.q)


    # return sensor_levels
    # Create observation geometry
    # Suniti TODO
    obs_geom = ObsGeometry(params.sza, params.vza, params.vaz, params.obs_alt)   
    # obs_geom = ObsGeometry(params.sza, params.vza, params.vaz, params.obs_alt, Array(sensor_levels))

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
    greek_cabannes = Vector{vSmartMOM.Scattering.GreekCoefs{Float64}}()
    greek_rayleigh = Vector{vSmartMOM.Scattering.GreekCoefs{Float64}}()
    ϖ_Cabannes = zeros(n_bands)
    # Remove rayleigh for testing:
    τ_rayl = [zeros(params.float_type,length(params.spec_bands[i]), length(profile.T)) for i=1:n_bands];
    #τ_rayl = [zeros(params.float_type,1,length(profile.T)) for i=1:n_bands];
    
    # This is a kludge for now, tau_abs sometimes needs to be a dual. Suniti & us need to rethink this all!!
    # i.e. code the rt core with fixed amount of derivatives as in her paper, then compute chain rule for dtau/dVMr, etc...
    FT2 = isnothing(params.absorption_params) || !haskey(params.absorption_params.vmr,"CO2") ? params.float_type : eltype(params.absorption_params.vmr["CO2"])
    τ_abs     = [zeros(FT2, length(params.spec_bands[i]), length(profile.p_full)) for i in 1:n_bands]
    max_m = zeros(Int, n_bands)
    l_max = zeros(Int, n_bands)
    l_max_aer = zeros(Int, n_aer, n_bands)
    # Loop over all bands:
    for i_band=1:n_bands

        # i'th spectral band (convert from cm⁻¹ to μm)
        curr_band_λ = 1e4 ./ params.spec_bands[i_band]
        # @show profile.vcd_dry, size(τ_rayl[i_band])
        # Compute Rayleigh properties per layer for `i_band` band center  
        #Suniti: the following two lines are temporary. this Cabannes depolarization applies only to the Earth's atmosphere
        # It has been added to make sure that the code sees the same Rayleigh cross section, regardless of elastic or inelastic RT 
        # Needs better (more general) formulation 
        νₘ = 0.5*(params.spec_bands[i_band][1]+params.spec_bands[i_band][end])
        λₘ = 1.e7/νₘ
        #@show i_band
        ϖ_Cabannes[i_band], γ_air_Cabannes, γ_air_Rayleigh = 
            InelasticScattering.compute_γ_air_Rayleigh!(λₘ)
        #@show ϖ_Cabannes[i_band]
        depol_air_Cabannes = 2γ_air_Cabannes/(1+γ_air_Cabannes)
        depol_air_Rayleigh = 2γ_air_Rayleigh/(1+γ_air_Rayleigh)
        
        a = Scattering.get_greek_rayleigh(depol_air_Cabannes)
        push!(greek_cabannes, a)
        a = Scattering.get_greek_rayleigh(depol_air_Rayleigh)
        push!(greek_rayleigh, a)

        τ_rayl[i_band]   .= getRayleighLayerOptProp(profile.p_half[end], 
            curr_band_λ, 
            depol_air_Rayleigh, profile.vcd_dry);
        #τ_rayl[i_band] *= (RS_type.n2.vmr + RS_type.o2.vmr) #this factor has been added to account for computations of isolated o2 and n2 contributions to VS
        
        #=
        # Original form
        τ_rayl[i_band]   .= getRayleighLayerOptProp(profile.p_half[end], 
                                (mean(curr_band_λ)), 
                                params.depol, profile.vcd_dry);
        =#
        
                                #@show τ_rayl[i_band]

        # If no absorption, continue to next band
        (isnothing(params.absorption_params) && isnothing(params.q)) && continue

        if !isnothing(params.q)
            # Obtain hitran data for this H₂O
            @timeit "Read HITRAN" hitran_data = read_hitran(artifact("H2O"), iso=1)

            println("Computing profile for water vapor in band #$(i_band)")
                # Create absorption model with parameters beforehand now:
                absorption_model = make_hitran_model(hitran_data, 
                    vSmartMOM.Absorption.Voigt(), 
                    wing_cutoff = 150, 
                    CEF = vSmartMOM.Absorption.HumlicekWeidemann32SDErrorFunction(), 
                    architecture = params.architecture, 
                    vmr = 0);#mean(profile.vmr[params.absorption_params.molecules[i_band][molec_i]]))
                # Calculate absorption profile
                
                @timeit "Absorption Coeff"  compute_absorption_profile!(τ_abs[i_band], 
                    absorption_model, 
                    params.spec_bands[i_band],
                    profile.vmr_h2o, 
                    profile);
        end
        #@show "hello"
        #@show params.absorption_params.molecules[i_band]
        #@show length(params.absorption_params.molecules[i_band])
        # Loop over all molecules in this band, obtain profile for each, and add them up
        if !isnothing(params.absorption_params)
            for molec_i in 1:length(params.absorption_params.molecules[i_band])
                
                # This can be precomputed as well later in my mind, providing an absorption_model or an interpolation_model!
                if isempty(params.absorption_params.luts)
                    # Obtain hitran data for this molecule
                    @timeit "Read HITRAN" hitran_data = read_hitran(artifact(params.absorption_params.molecules[i_band][molec_i]), iso=1)

                    println("Computing profile for $(params.absorption_params.molecules[i_band][molec_i]) with vmr $(profile.vmr[params.absorption_params.molecules[i_band][molec_i]]) for band #$(i_band)")
                    # Create absorption model with parameters beforehand now:
                    absorption_model = make_hitran_model(hitran_data, 
                        params.absorption_params.broadening_function, 
                        wing_cutoff = params.absorption_params.wing_cutoff, 
                        CEF = params.absorption_params.CEF, 
                        architecture = params.architecture, 
                        vmr = 0);#mean(profile.vmr[params.absorption_params.molecules[i_band][molec_i]]))
                    # Calculate absorption profile
                    
                    @timeit "Absorption Coeff"  compute_absorption_profile!(τ_abs[i_band], 
                        absorption_model, 
                            params.spec_bands[i_band],
                                profile.vmr[params.absorption_params.molecules[i_band][molec_i]], 
                                profile);
                # Use LUT directly
                else
                    compute_absorption_profile!(τ_abs[i_band], params.absorption_params.luts[i_band][molec_i], params.spec_bands[i_band],profile.vmr[params.absorption_params.molecules[i_band][molec_i]], profile);
                end
            end
        end
    end

    # aerosol_optics[iBand][iAer]
    aerosol_optics = [Array{AerosolOptics}(undef, (n_aer)) for i=1:n_bands];
        
    FT2 = isnothing(params.scattering_params) ? params.float_type : typeof(params.scattering_params.rt_aerosols[1].τ_ref)
    FT2 =  params.float_type 

    # τ_aer[iBand][iAer,iZ]        
    τ_aer = [zeros(FT2, n_aer, length(params.spec_bands[i]), length(profile.p_full)) for i=1:n_bands];

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
        k_ref          = compute_ref_aerosol_extinction(mie_model, params.float_type)

        #params.scattering_params.rt_aerosols[i_aer].p₀, params.scattering_params.rt_aerosols[i_aer].σp
        # Loop over bands
        for i_band=1:n_bands
            
            # i'th spectral band (convert from cm⁻¹ to μm)
            curr_band_λ = 1e4 ./ params.spec_bands[i_band]

            # Create the aerosols:
            if length(curr_band_λ)==1
                mie_model      = make_mie_model(params.scattering_params.decomp_type, 
                                            mie_aerosol, 
                                            (maximum(curr_band_λ)+minimum(curr_band_λ))/2, 
                                            params.polarization_type, 
                                            truncation_type, 
                                            params.scattering_params.r_max, 
                                            params.scattering_params.nquad_radius)
                # Compute raw (not truncated) aerosol optical properties (not needed in RT eventually) 
                @timeit "Mie calc"  aerosol_optics_raw = compute_aerosol_optical_properties(mie_model, FT2);

            else
                # NOTE: in the following, the spectrum is ordered by increasing wavenumber. 
                # As a result, curr_band_λ[0] > curr_band_λ[end]. To avoid confusion, interpolation 
                # is carried out on a wavenumber grid
                mie_model_0    = make_mie_model(params.scattering_params.decomp_type, 
                                            mie_aerosol, 
                                            curr_band_λ[1], 
                                            params.polarization_type, 
                                            truncation_type, 
                                            params.scattering_params.r_max, 
                                            params.scattering_params.nquad_radius)
                # Compute raw (not truncated) aerosol optical properties (not needed in RT eventually) 
                @timeit "Mie calc"  aerosol_optics_raw_0 = compute_aerosol_optical_properties(mie_model_0, FT2);

                mie_model_1    = make_mie_model(params.scattering_params.decomp_type, 
                                            mie_aerosol, 
                                            curr_band_λ[end], 
                                            params.polarization_type, 
                                            truncation_type, 
                                            params.scattering_params.r_max, 
                                            params.scattering_params.nquad_radius)
                # Compute raw (not truncated) aerosol optical properties (not needed in RT eventually) 
                @timeit "Mie calc"  aerosol_optics_raw_1 = compute_aerosol_optical_properties(mie_model_1, FT2);
                Nl = length(aerosol_optics_raw_1.greek_coefs.α)
                Nl_ = length(aerosol_optics_raw_0.greek_coefs.α)
                α = zeros(Nl)
                β = zeros(Nl)
                γ = zeros(Nl)
                δ = zeros(Nl)
                ϵ = zeros(Nl)
                ζ = zeros(Nl)
                
                α[1:Nl_] = 0.5*(aerosol_optics_raw_0.greek_coefs.α[1:Nl_] + aerosol_optics_raw_1.greek_coefs.α[1:Nl_])
                β[1:Nl_] = 0.5*(aerosol_optics_raw_0.greek_coefs.β[1:Nl_] + aerosol_optics_raw_1.greek_coefs.β[1:Nl_])
                γ[1:Nl_] = 0.5*(aerosol_optics_raw_0.greek_coefs.γ[1:Nl_] + aerosol_optics_raw_1.greek_coefs.γ[1:Nl_])
                δ[1:Nl_] = 0.5*(aerosol_optics_raw_0.greek_coefs.δ[1:Nl_] + aerosol_optics_raw_1.greek_coefs.δ[1:Nl_])
                ϵ[1:Nl_] = 0.5*(aerosol_optics_raw_0.greek_coefs.ϵ[1:Nl_] + aerosol_optics_raw_1.greek_coefs.ϵ[1:Nl_])
                ζ[1:Nl_] = 0.5*(aerosol_optics_raw_0.greek_coefs.ζ[1:Nl_] + aerosol_optics_raw_1.greek_coefs.ζ[1:Nl_])
            
                α[1+Nl_:Nl] = 0.5*(aerosol_optics_raw_1.greek_coefs.α[1+Nl_:Nl])
                β[1+Nl_:Nl] = 0.5*(aerosol_optics_raw_1.greek_coefs.β[1+Nl_:Nl])
                γ[1+Nl_:Nl] = 0.5*(aerosol_optics_raw_1.greek_coefs.γ[1+Nl_:Nl])
                δ[1+Nl_:Nl] = 0.5*(aerosol_optics_raw_1.greek_coefs.δ[1+Nl_:Nl])
                ϵ[1+Nl_:Nl] = 0.5*(aerosol_optics_raw_1.greek_coefs.ϵ[1+Nl_:Nl])
                ζ[1+Nl_:Nl] = 0.5*(aerosol_optics_raw_1.greek_coefs.ζ[1+Nl_:Nl])
                
                greek_coeffs = GreekCoefs(α, β, γ, δ, ϵ, ζ)
                ν_grid = [1e4/curr_band_λ[1], 1e4/curr_band_λ[end]] 
                #@show 1e4/curr_band_λ[1], 1e4/curr_band_λ[end]
                #@show aerosol_optics_raw_0.k, aerosol_optics_raw_1.k
                #@show aerosol_optics_raw_0.k*aerosol_optics_raw_0.ω̃, aerosol_optics_raw_1.k*aerosol_optics_raw_1.ω̃
                kext_grid = [aerosol_optics_raw_0.k, aerosol_optics_raw_1.k]
                ksca_grid = [aerosol_optics_raw_0.k*aerosol_optics_raw_0.ω̃, aerosol_optics_raw_1.k*aerosol_optics_raw_1.ω̃] 
                interp_linear_kext = linear_interpolation(ν_grid, kext_grid)
                interp_linear_ksca = linear_interpolation(ν_grid, ksca_grid)
                k = zeros(length(curr_band_λ))
                ω̃ = zeros(length(curr_band_λ))
                for i = 1:length(curr_band_λ)
                    k[i] = interp_linear_kext(1e4/curr_band_λ[i])
                    ω̃[i] = interp_linear_ksca(1e4/curr_band_λ[i])/k[i]
                end
                aerosol_optics_raw = AerosolOptics(greek_coefs=greek_coeffs, ω̃=ω̃, k=k, fᵗ=0.0)
            end
            
            # Compute truncated aerosol optical properties (phase function and fᵗ), consistent with Ltrunc:
            #@show i_aer, i_band
            #@show length(aerosol_optics_raw.greek_coefs.β), truncation_type.l_max
            
            if length(aerosol_optics_raw.greek_coefs.β) > truncation_type.l_max
                aerosol_optics[i_band][i_aer] = Scattering.truncate_phase(truncation_type, 
                                                    aerosol_optics_raw; reportFit=false)
                l_max_aer[i_aer, i_band] = truncation_type.l_max
                #max_m[i_band] = Int(ceil(l_max[i_band] + 1)/2)                                 
            else
                #@show truncation_type.l_max
                aerosol_optics[i_band][i_aer] = aerosol_optics_raw
                l_max_aer[i_aer, i_band] = length(aerosol_optics_raw.greek_coefs.β)
                #max_m[i_band] = Int(ceil(l_max[i_band] + 1)/2)
                #@show max_m[i_band]
            end
            #@show size(getAerosolLayerOptProp(1, c_aero.z₀, c_aero.σ₀, profile.p_half, profile.T))
            #@show size(aerosol_optics[i_band][i_aer].k), k_ref
            #@show size(τ_aer[i_band][i_aer,:,:])
            # Compute nAer aerosol optical thickness profiles
            τ_aer[i_band][i_aer,:,:] = 
                (params.scattering_params.rt_aerosols[i_aer].τ_ref/k_ref) * 
                aerosol_optics[i_band][i_aer].k *
                (getAerosolLayerOptProp(1, c_aero.z₀, c_aero.σ₀, profile.p_half, profile.T))' 
                             
            #@show size(τ_aer[i_band][i_aer,:,:])
            #@show params.scattering_params.rt_aerosols[i_aer].τ_ref
        end                  
    end
    for i_band = 1:n_bands
        l_max[i_band] = maximum(l_max_aer[:,i_band]) 
        max_m[i_band] = Int(ceil(l_max[i_band] + 1)/2)  
    end
    set_uniform_lmax!(l_max, aerosol_optics)
#@show typeof(τ_aer)
#@show typeof(τ_rayl)

    # Check the floating-type output matches specified FT
#@show size(ϖ_Cabannes), typeof(ϖ_Cabannes)
    # Return the model 
    return vSmartMOM_Model(
                        max_m,
                        l_max,
                        params, 
                        aerosol_optics, 
                        ϖ_Cabannes, 
                        greek_cabannes,  
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


"Take the parameters specified in the vSmartMOM_Parameters struct, 
and calculate derived attributes into a vSmartMOM_Model" 
function model_from_parameters(RS_type::Union{VS_0to1_plus, VS_1to0_plus}, 
                    λ₀,
                    params::vSmartMOM_Parameters)
    #@show params.absorption_params.molecules
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
    #RS_type.n2, RS_type.o2 = 
    #    InelasticScattering.getRamanAtmoConstants(1.e7/λ₀,effT);
    ϖ_Cabannes = zeros(n_bands)
    ϖ_Cabannes[1], γ_air_Cabannes, γ_air_Rayleigh = 
        InelasticScattering.compute_γ_air_Rayleigh!(λ₀)
    depol_air_Cabannes = 2γ_air_Cabannes/(1+γ_air_Cabannes)
    depol_air_Rayleigh = 2γ_air_Rayleigh/(1+γ_air_Rayleigh)
    #println("here 0")
    InelasticScattering.getRamanSSProp!(RS_type, depol_air_Rayleigh, λ₀);
    #println("here 1")
    n_bands = length(RS_type.iBand)
    #@show RS_type.grid_in
    params.spec_bands = RS_type.grid_in
    #@show params.spec_bands

    # Rayleigh optical properties calculation
    #greek_rayleigh = Scattering.get_greek_rayleigh(params.depol)
    # Rayleigh optical properties calculation
    
    greek_cabannes = Vector{vSmartMOM.Scattering.GreekCoefs{Float64}}()
    greek_rayleigh = Vector{vSmartMOM.Scattering.GreekCoefs{Float64}}()
    τ_rayl = [zeros(params.float_type,length(params.spec_bands[i]), length(profile.T)) for i=1:n_bands];
    #τ_rayl = [zeros(params.float_type,1, length(profile.p_full)) for i=1:n_bands];

    # τ_abs[iBand][iSpec,iZ]
    τ_abs     = [zeros(params.float_type, length(params.spec_bands[i]), length(profile.p_full)) for i in 1:n_bands]
    #@show params.absorption_params.molecules[2]
    max_m = zeros(Int, n_bands)
    l_max = zeros(Int, n_bands)
    l_max_aer = zeros(Int, n_aer, n_bands)
    # Loop over all bands:
    for i_band=1:n_bands
        #@show params.spec_bands[i_band]
        # i'th spectral band (convert from cm⁻¹ to μm)
        curr_band_λ = 1e4 ./ params.spec_bands[i_band]
        νₘ = 0.5*(params.spec_bands[i_band][1]+params.spec_bands[i_band][end])
        λₘ = 1.e7/νₘ
        ϖ_Cabannes[i_band], γ_air_Cabannes, γ_air_Rayleigh = 
            InelasticScattering.compute_γ_air_Rayleigh!(λₘ)
        if(i_band==1) 
            ϖ_Cabannes[i_band]=1.
        end
        depol_air_Cabannes = 2γ_air_Cabannes/(1+γ_air_Cabannes)
        depol_air_Rayleigh = 2γ_air_Rayleigh/(1+γ_air_Rayleigh)

        #ϖ_Cabannes_VS[i_band], γ_air_Cabannes_VS, γ_air_Rayleigh_VS = 
        #InelasticScattering.compute_γ_air_Rayleigh_VS!(λₘ)

        #depol_air_Cabannes_VS = 2γ_air_Cabannes_VS/(1+γ_air_Cabannes_VS)
        #depol_air_Rayleigh_VS = 2γ_air_Rayleigh_VS/(1+γ_air_Rayleigh_VS)

        a = Scattering.get_greek_rayleigh(depol_air_Cabannes)
        push!(greek_cabannes, a)
        a = Scattering.get_greek_rayleigh(depol_air_Rayleigh)
        push!(greek_rayleigh, a)
        # Compute Rayleigh properties per layer for `i_band` band center
        # (Using depol_air_Rayleigh here to obtain the same Rayleigh optical thickness as in the case of RS_type::noRS or RS_type::RRS) 
        τ_rayl[i_band]   .= getRayleighLayerOptProp(profile.p_half[end], 
            curr_band_λ, #λₘ/1000., 
            depol_air_Rayleigh, profile.vcd_dry);
        τ_rayl[i_band] *= (RS_type.n2.vmr + RS_type.o2.vmr) #this factor has been added to account for computations of isolated o2 and n2 contributions to VS
        #=
        # Compute Rayleigh properties per layer for `i_band` band center
        τ_rayl[i_band]   .= getRayleighLayerOptProp(profile.p_half[end], 
                            (maximum(curr_band_λ) + minimum(curr_band_λ))/2, 
                            params.depol, profile.vcd_dry);
        =#
        # If no absorption, continue to next band
        isnothing(params.absorption_params) && continue
        #@show i_band, params.absorption_params.molecules[i_band]
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
        k_ref          = compute_ref_aerosol_extinction(mie_model, params.float_type)

        #params.scattering_params.rt_aerosols[i_aer].p₀, params.scattering_params.rt_aerosols[i_aer].σp
        # Loop over bands
        for i_band=1:n_bands
            
            # i'th spectral band (convert from cm⁻¹ to μm)
            curr_band_λ = 1e4 ./ params.spec_bands[i_band]

            # Create the aerosols:
            if length(curr_band_λ)==1
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
            else
                # NOTE: in the following, the spectrum is ordered by increasing wavenumber. 
                # As a result, curr_band_λ[0] > curr_band_λ[end]. To avoid confusion, interpolation 
                # is carried out on a wavenumber grid
                mie_model_0    = make_mie_model(params.scattering_params.decomp_type, 
                                            mie_aerosol, 
                                            curr_band_λ[1], 
                                            params.polarization_type, 
                                            truncation_type, 
                                            params.scattering_params.r_max, 
                                            params.scattering_params.nquad_radius)
                # Compute raw (not truncated) aerosol optical properties (not needed in RT eventually) 
                @timeit "Mie calc"  aerosol_optics_raw_0 = compute_aerosol_optical_properties(mie_model_0, FT2);

                mie_model_1    = make_mie_model(params.scattering_params.decomp_type, 
                                            mie_aerosol, 
                                            curr_band_λ[end], 
                                            params.polarization_type, 
                                            truncation_type, 
                                            params.scattering_params.r_max, 
                                            params.scattering_params.nquad_radius)
                # Compute raw (not truncated) aerosol optical properties (not needed in RT eventually) 
                @timeit "Mie calc"  aerosol_optics_raw_1 = compute_aerosol_optical_properties(mie_model_1, FT2);
                Nl = length(aerosol_optics_raw_1.greek_coefs.α)
                Nl_ = length(aerosol_optics_raw_0.greek_coefs.α)
                α = zeros(Nl)
                β = zeros(Nl)
                γ = zeros(Nl)
                δ = zeros(Nl)
                ϵ = zeros(Nl)
                ζ = zeros(Nl)
                
                α[1:Nl_] = 0.5*(aerosol_optics_raw_0.greek_coefs.α[1:Nl_] + aerosol_optics_raw_1.greek_coefs.α[1:Nl_])
                β[1:Nl_] = 0.5*(aerosol_optics_raw_0.greek_coefs.β[1:Nl_] + aerosol_optics_raw_1.greek_coefs.β[1:Nl_])
                γ[1:Nl_] = 0.5*(aerosol_optics_raw_0.greek_coefs.γ[1:Nl_] + aerosol_optics_raw_1.greek_coefs.γ[1:Nl_])
                δ[1:Nl_] = 0.5*(aerosol_optics_raw_0.greek_coefs.δ[1:Nl_] + aerosol_optics_raw_1.greek_coefs.δ[1:Nl_])
                ϵ[1:Nl_] = 0.5*(aerosol_optics_raw_0.greek_coefs.ϵ[1:Nl_] + aerosol_optics_raw_1.greek_coefs.ϵ[1:Nl_])
                ζ[1:Nl_] = 0.5*(aerosol_optics_raw_0.greek_coefs.ζ[1:Nl_] + aerosol_optics_raw_1.greek_coefs.ζ[1:Nl_])
            
                α[1+Nl_:Nl] = 0.5*(aerosol_optics_raw_1.greek_coefs.α[1+Nl_:Nl])
                β[1+Nl_:Nl] = 0.5*(aerosol_optics_raw_1.greek_coefs.β[1+Nl_:Nl])
                γ[1+Nl_:Nl] = 0.5*(aerosol_optics_raw_1.greek_coefs.γ[1+Nl_:Nl])
                δ[1+Nl_:Nl] = 0.5*(aerosol_optics_raw_1.greek_coefs.δ[1+Nl_:Nl])
                ϵ[1+Nl_:Nl] = 0.5*(aerosol_optics_raw_1.greek_coefs.ϵ[1+Nl_:Nl])
                ζ[1+Nl_:Nl] = 0.5*(aerosol_optics_raw_1.greek_coefs.ζ[1+Nl_:Nl])
                
                greek_coeffs = GreekCoefs(α, β, γ, δ, ϵ, ζ)
                ν_grid = [1e4/curr_band_λ[1], 1e4/curr_band_λ[end]] 
                #@show 1e4/curr_band_λ[1], 1e4/curr_band_λ[end]
                #@show aerosol_optics_raw_0.k, aerosol_optics_raw_1.k
                #@show aerosol_optics_raw_0.k*aerosol_optics_raw_0.ω̃, aerosol_optics_raw_1.k*aerosol_optics_raw_1.ω̃
                kext_grid = [aerosol_optics_raw_0.k, aerosol_optics_raw_1.k]
                ksca_grid = [aerosol_optics_raw_0.k*aerosol_optics_raw_0.ω̃, aerosol_optics_raw_1.k*aerosol_optics_raw_1.ω̃] 
                interp_linear_kext = linear_interpolation(ν_grid, kext_grid)
                interp_linear_ksca = linear_interpolation(ν_grid, ksca_grid)
                k = zeros(length(curr_band_λ))
                ω̃ = zeros(length(curr_band_λ))
                for i = 1:length(curr_band_λ)
                    k[i] = interp_linear_kext(1e4/curr_band_λ[i])
                    ω̃[i] = interp_linear_ksca(1e4/curr_band_λ[i])/k[i]
                end
                aerosol_optics_raw = AerosolOptics(greek_coefs=greek_coeffs, ω̃=ω̃, k=k, fᵗ=0.0)
            end
            # Compute truncated aerosol optical properties (phase function and fᵗ), consistent with Ltrunc:
            #@show i_aer, i_band
            if length(aerosol_optics_raw.greek_coefs.β) > truncation_type.l_max
                aerosol_optics[i_band][i_aer] = Scattering.truncate_phase(truncation_type, 
                                                    aerosol_optics_raw; reportFit=false)
                l_max_aer[i_aer, i_band] = truncation_type.l_max
                #max_m[i_band] = Int(ceil(l_max[i_band] + 1)/2)                                 
            else
                #@show truncation_type.l_max
                aerosol_optics[i_band][i_aer] = aerosol_optics_raw
                l_max_aer[i_aer, i_band] = length(aerosol_optics_raw.greek_coefs.β)
                #max_m[i_band] = Int(ceil(l_max[i_band] + 1)/2)
                #@show max_m[i_band]
            end
            #aerosol_optics[i_band][i_aer] = Scattering.truncate_phase(truncation_type, 
            #                                        aerosol_optics_raw; reportFit=false)

            # Compute nAer aerosol optical thickness profiles
            τ_aer[i_band][i_aer,:,:] = 
                (params.scattering_params.rt_aerosols[i_aer].τ_ref/k_ref) * 
                aerosol_optics[i_band][i_aer].k *
                (getAerosolLayerOptProp(1, c_aero.z₀, c_aero.σ₀, profile.p_half, profile.T))' 
            #τ_aer[i_band][i_aer,:] = 
            #    params.scattering_params.rt_aerosols[i_aer].τ_ref * 
            #    (aerosol_optics[i_band][i_aer].k/k_ref) * 
            #    CoreRT.getAerosolLayerOptProp(1.0, c_aero.p₀, c_aero.σp, profile.p_half)
        end 
    end
    for i_band = 1:n_bands
        l_max[i_band] = maximum(l_max_aer[:,i_band]) 
        max_m[i_band] = Int(ceil(l_max[i_band] + 1)/2)  
    end
    set_uniform_lmax!(l_max, aerosol_optics)
    # Check the floating-type output matches specified FT

    # Return the model 
    return vSmartMOM_Model(
                        max_m,
                        l_max,
                        params, 
                        aerosol_optics, 
                        ϖ_Cabannes, 
                        greek_cabannes, 
                        greek_rayleigh, 
                        quad_points, 
                        τ_abs, 
                        τ_rayl, 
                        τ_aer, 
                        obs_geom, 
                        profile)
end

"Take the parameters specified in the vSmartMOM_Parameters struct, 
and calculate derived attributes into a vSmartMOM_Model for VRS and 
RVRS to a specified target (spectral) grid" 
function model_from_parameters(RS_type::Union{VS_0to1_plus, VS_1to0_plus}, 
                    λ₀,
                    params::vSmartMOM_Parameters, 
                    target_grid::AbstractArray{FT})where {FT<:Union{AbstractFloat, ForwardDiff.Dual}}
    
    #@show params.absorption_params.molecules
    # Number of total bands and aerosols (for convenience)
    n_bands = 2 #length(params.spec_bands)
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
    #RS_type.n2, RS_type.o2 = 
    #    InelasticScattering.getRamanAtmoConstants(1.e7/λ₀,effT);
    #println("here 0")
    InelasticScattering.getRamanSSProp!(RS_type, λ₀, target_grid);
    #println("here 1")
    n_bands = length(RS_type.iBand)
    #@show RS_type.grid_in
    params.spec_bands = RS_type.grid_in
    #@show params.spec_bands

    # Rayleigh optical properties calculation
    ϖ_Cabannes = zeros(n_bands)
    greek_cabannes = Vector{vSmartMOM.Scattering.GreekCoefs{Float64}}()
    greek_rayleigh = Vector{vSmartMOM.Scattering.GreekCoefs{Float64}}()
    #greek_rayleigh = Scattering.get_greek_rayleigh(params.depol)
    τ_rayl = [zeros(params.float_type,length(params.spec_bands[i]), length(profile.T)) for i=1:n_bands];
    #τ_rayl = [zeros(params.float_type,1, length(profile.p_full)) for i=1:n_bands];

    # τ_abs[iBand][iSpec,iZ]
    τ_abs     = [zeros(params.float_type, length(params.spec_bands[i]), length(profile.p_full)) for i in 1:n_bands]
    #@show params.absorption_params.molecules[2]
    max_m = zeros(Int, n_bands)
    l_max = zeros(Int, n_bands)
    l_max_aer = zeros(Int, n_aer, n_bands)
    # Loop over all bands:
    for i_band=1:n_bands
        #@show params.spec_bands[i_band]
        # i'th spectral band (convert from cm⁻¹ to μm)
        curr_band_λ = 1e4 ./ params.spec_bands[i_band]
        νₘ = 0.5*(params.spec_bands[i_band][1]+params.spec_bands[i_band][end])
        λₘ = 1.e7/νₘ
        ϖ_Cabannes[i_band], γ_air_Cabannes, γ_air_Rayleigh = 
            InelasticScattering.compute_γ_air_Rayleigh!(λₘ)
        depol_air_Cabannes = 2γ_air_Cabannes/(1+γ_air_Cabannes)
        depol_air_Rayleigh = 2γ_air_Rayleigh/(1+γ_air_Rayleigh)
        a = Scattering.get_greek_rayleigh(depol_air_Cabannes)
        push!(greek_cabannes, a)
        a = Scattering.get_greek_rayleigh(depol_air_Rayleigh)
        push!(greek_rayleigh, a)
        # Compute Rayleigh properties per layer for `i_band` band center
        τ_rayl[i_band]   .= getRayleighLayerOptProp(profile.p_half[end], 
            curr_band_λ, #λₘ/1000., 
            depol_air_Rayleigh, profile.vcd_dry);
        τ_rayl[i_band] *= (RS_type.n2.vmr + RS_type.o2.vmr) #this factor has been added to account for computations of isolated o2 and n2 contributions to VS
        
#=
        # Compute Rayleigh properties per layer for `i_band` band center
        τ_rayl[i_band]   .= getRayleighLayerOptProp(profile.p_half[end], 
                            (minimum(curr_band_λ) + maximum(curr_band_λ))/2, 
                            params.depol, profile.vcd_dry);
=#
        # If no absorption, continue to next band
        isnothing(params.absorption_params) && continue
        #@show i_band, params.absorption_params.molecules[i_band]
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
        k_ref          = compute_ref_aerosol_extinction(mie_model, params.float_type)

        #params.scattering_params.rt_aerosols[i_aer].p₀, params.scattering_params.rt_aerosols[i_aer].σp
        # Loop over bands
        for i_band=1:n_bands
            
            # i'th spectral band (convert from cm⁻¹ to μm)
            curr_band_λ = 1e4 ./ params.spec_bands[i_band]

            # Create the aerosols:
            if length(curr_band_λ)==1
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
            else
                # NOTE: in the following, the spectrum is ordered by increasing wavenumber. 
                # As a result, curr_band_λ[0] > curr_band_λ[end]. To avoid confusion, interpolation 
                # is carried out on a wavenumber grid
                mie_model_0    = make_mie_model(params.scattering_params.decomp_type, 
                                            mie_aerosol, 
                                            curr_band_λ[1], 
                                            params.polarization_type, 
                                            truncation_type, 
                                            params.scattering_params.r_max, 
                                            params.scattering_params.nquad_radius)
                # Compute raw (not truncated) aerosol optical properties (not needed in RT eventually) 
                @timeit "Mie calc"  aerosol_optics_raw_0 = compute_aerosol_optical_properties(mie_model_0, FT2);

                mie_model_1    = make_mie_model(params.scattering_params.decomp_type, 
                                            mie_aerosol, 
                                            curr_band_λ[end], 
                                            params.polarization_type, 
                                            truncation_type, 
                                            params.scattering_params.r_max, 
                                            params.scattering_params.nquad_radius)
                # Compute raw (not truncated) aerosol optical properties (not needed in RT eventually) 
                @timeit "Mie calc"  aerosol_optics_raw_1 = compute_aerosol_optical_properties(mie_model_1, FT2);
                Nl = length(aerosol_optics_raw_1.greek_coefs.α)
                Nl_ = length(aerosol_optics_raw_0.greek_coefs.α)
                α = zeros(Nl)
                β = zeros(Nl)
                γ = zeros(Nl)
                δ = zeros(Nl)
                ϵ = zeros(Nl)
                ζ = zeros(Nl)
                
                α[1:Nl_] = 0.5*(aerosol_optics_raw_0.greek_coefs.α[1:Nl_] + aerosol_optics_raw_1.greek_coefs.α[1:Nl_])
                β[1:Nl_] = 0.5*(aerosol_optics_raw_0.greek_coefs.β[1:Nl_] + aerosol_optics_raw_1.greek_coefs.β[1:Nl_])
                γ[1:Nl_] = 0.5*(aerosol_optics_raw_0.greek_coefs.γ[1:Nl_] + aerosol_optics_raw_1.greek_coefs.γ[1:Nl_])
                δ[1:Nl_] = 0.5*(aerosol_optics_raw_0.greek_coefs.δ[1:Nl_] + aerosol_optics_raw_1.greek_coefs.δ[1:Nl_])
                ϵ[1:Nl_] = 0.5*(aerosol_optics_raw_0.greek_coefs.ϵ[1:Nl_] + aerosol_optics_raw_1.greek_coefs.ϵ[1:Nl_])
                ζ[1:Nl_] = 0.5*(aerosol_optics_raw_0.greek_coefs.ζ[1:Nl_] + aerosol_optics_raw_1.greek_coefs.ζ[1:Nl_])
            
                α[1+Nl_:Nl] = 0.5*(aerosol_optics_raw_1.greek_coefs.α[1+Nl_:Nl])
                β[1+Nl_:Nl] = 0.5*(aerosol_optics_raw_1.greek_coefs.β[1+Nl_:Nl])
                γ[1+Nl_:Nl] = 0.5*(aerosol_optics_raw_1.greek_coefs.γ[1+Nl_:Nl])
                δ[1+Nl_:Nl] = 0.5*(aerosol_optics_raw_1.greek_coefs.δ[1+Nl_:Nl])
                ϵ[1+Nl_:Nl] = 0.5*(aerosol_optics_raw_1.greek_coefs.ϵ[1+Nl_:Nl])
                ζ[1+Nl_:Nl] = 0.5*(aerosol_optics_raw_1.greek_coefs.ζ[1+Nl_:Nl])
                
                greek_coeffs = GreekCoefs(α, β, γ, δ, ϵ, ζ)
                ν_grid = [1e4/curr_band_λ[1], 1e4/curr_band_λ[end]] 
                #@show 1e4/curr_band_λ[1], 1e4/curr_band_λ[end]
                #@show aerosol_optics_raw_0.k, aerosol_optics_raw_1.k
                #@show aerosol_optics_raw_0.k*aerosol_optics_raw_0.ω̃, aerosol_optics_raw_1.k*aerosol_optics_raw_1.ω̃
                kext_grid = [aerosol_optics_raw_0.k, aerosol_optics_raw_1.k]
                ksca_grid = [aerosol_optics_raw_0.k*aerosol_optics_raw_0.ω̃, aerosol_optics_raw_1.k*aerosol_optics_raw_1.ω̃] 
                interp_linear_kext = linear_interpolation(ν_grid, kext_grid)
                interp_linear_ksca = linear_interpolation(ν_grid, ksca_grid)
                k = zeros(length(curr_band_λ))
                ω̃ = zeros(length(curr_band_λ))
                for i = 1:length(curr_band_λ)
                    k[i] = interp_linear_kext(1e4/curr_band_λ[i])
                    ω̃[i] = interp_linear_ksca(1e4/curr_band_λ[i])/k[i]
                end
                aerosol_optics_raw = AerosolOptics(greek_coefs=greek_coeffs, ω̃=ω̃, k=k, fᵗ=0.0)
            end
            # Compute truncated aerosol optical properties (phase function and fᵗ), consistent with Ltrunc:
            #@show i_aer, i_band
            if length(aerosol_optics_raw.greek_coefs.β) > truncation_type.l_max
                aerosol_optics[i_band][i_aer] = Scattering.truncate_phase(truncation_type, 
                                                    aerosol_optics_raw; reportFit=false)
                l_max_aer[i_aer, i_band] = truncation_type.l_max
                #max_m[i_band] = Int(ceil(l_max[i_band] + 1)/2)                                 
            else
                #@show truncation_type.l_max
                aerosol_optics[i_band][i_aer] = aerosol_optics_raw
                l_max_aer[i_aer, i_band] = length(aerosol_optics_raw.greek_coefs.β)
                #max_m[i_band] = Int(ceil(l_max[i_band] + 1)/2)
                #@show max_m[i_band]
            end
            #aerosol_optics[i_band][i_aer] = Scattering.truncate_phase(truncation_type, 
            #                                        aerosol_optics_raw; reportFit=false)

            # Compute nAer aerosol optical thickness profiles
            τ_aer[i_band][i_aer,:,:] = 
                (params.scattering_params.rt_aerosols[i_aer].τ_ref/k_ref) * 
                aerosol_optics[i_band][i_aer].k *
                (getAerosolLayerOptProp(1, c_aero.z₀, c_aero.σ₀, profile.p_half, profile.T))' 
            #τ_aer[i_band][i_aer,:] = 
            #    params.scattering_params.rt_aerosols[i_aer].τ_ref * 
            #    (aerosol_optics[i_band][i_aer].k/k_ref) * 
            #    CoreRT.getAerosolLayerOptProp(1.0, c_aero.p₀, c_aero.σp, profile.p_half)
        end 
    end
    for i_band = 1:n_bands
        l_max[i_band] = maximum(l_max_aer[:,i_band]) 
        max_m[i_band] = Int(ceil(l_max[i_band] + 1)/2)  
    end
    set_uniform_lmax!(l_max, aerosol_optics)
    # Check the floating-type output matches specified FT

    # Return the model 
    return vSmartMOM_Model(
                        max_m,
                        l_max,
                        params, 
                        aerosol_optics, 
                        ϖ_Cabannes, 
                        greek_cabannes, 
                        greek_rayleigh, 
                        quad_points, 
                        τ_abs, 
                        τ_rayl, 
                        τ_aer, 
                        obs_geom, 
                        profile)

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
