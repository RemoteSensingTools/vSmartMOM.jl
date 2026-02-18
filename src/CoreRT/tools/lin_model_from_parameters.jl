"Set a common l_max in a given band for each aerosol type"
function set_uniform_lmax!(lmax::Vector{Int}, aerosol_optics)
    n_bands = length(lmax)
    n_aer = length(aerosol_optics[1])
    for i_band = 1:n_bands
        pref_lmax = lmax[i_band]
        for i_aer = 1:n_aer
            aer_lmax = length(aerosol_optics[i_band][i_aer].greek_coefs.β)
            if aer_lmax < pref_lmax
                for fname in (:α, :β, :γ, :δ, :ϵ, :ζ)
                    old = getfield(aerosol_optics[i_band][i_aer].greek_coefs, fname)
                    new_arr = zeros(eltype(old), pref_lmax)
                    new_arr[1:aer_lmax] .= old
                    setfield!(aerosol_optics[i_band][i_aer].greek_coefs, fname, new_arr)
                end
            end
        end
    end
end

"""
    model_from_parameters(::LinMode, params::vSmartMOM_Parameters)

Construct both the forward `vSmartMOM_Model` and the linearized `vSmartMOM_Lin` objects
from the input parameters, for use in linearized (Jacobian) RT computations.

This is the **linearized** counterpart of `model_from_parameters(params)`. It computes:

1. **Forward optical properties** (same as the forward-only model):
   - Rayleigh optical depth per layer and spectral point.
   - Aerosol optical properties (via Mie theory) with δ-M truncation.
   - Trace gas absorption cross-sections and optical depths.

2. **Linearized optical properties** (derivatives of the above):
   - `τ̇_aer[iB][iaer, 7, nSpec, nLayers]`: Derivatives of aerosol τ w.r.t. 7 sub-parameters
     `[τ_ref, nᵣ, nᵢ, rₘ, σᵣ, p₀, σp]` per aerosol type.
   - `τ̇_abs[iB][NGas, nSpec, nLayers]`: Derivatives of gas absorption τ w.r.t. VMR.
   - `lin_aerosol_optics[iB][iaer]`: Derivatives of Mie properties (ω̃, fᵗ, greek coefficients)
     w.r.t. Mie parameters `[nᵣ, nᵢ, rₘ, σᵣ]`.

# Returns
- `model::vSmartMOM_Model`: Forward model (optical properties, geometry, quadrature).
- `lin_model::vSmartMOM_Lin`: Linearized model (all derivative arrays).

# Notes
- The atmospheric profile may be truncated to the observer altitude for tower/airborne sensors.
- Aerosol Mie calculations use `ForwardDiff.Dual` numbers to simultaneously obtain
  derivatives of the extinction cross-section, single-scattering albedo, truncation
  factor, and greek coefficients with respect to `[nᵣ, nᵢ, rₘ, σᵣ]`.
"""
function model_from_parameters(lin::LinMode, 
    params::vSmartMOM_Parameters)
    FT =  params.float_type 
    # Number of total bands and aerosols (for convenience)
    n_bands = length(params.spec_bands)
    n_aer = isnothing(params.scattering_params) ? 0 : length(params.scattering_params.rt_aerosols)

    
    # Get new p/T profiles using obs_alt
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
    p_full, p_half, vmr_h2o, vcd_dry, vcd_h2o, new_vmr, Δz = compute_atmos_profile_fields(params.T, params.p, params.q, vmr)

    profile = AtmosphericProfile(params.T, p_full, params.q, p_half, vmr_h2o, vcd_dry, vcd_h2o, new_vmr, Δz)
    
    # Reduce the profile to the number of target layers (if specified)
    if params.profile_reduction_n != -1
        profile = reduce_profile(params.profile_reduction_n, profile);
    end

    # Rayleigh optical properties calculation
    greek_cabannes = Vector{vSmartMOM.Scattering.GreekCoefs{FT}}()
    greek_rayleigh = Vector{vSmartMOM.Scattering.GreekCoefs{FT}}()
    ϖ_Cabannes = zeros(n_bands)
    # Remove rayleigh for testing:
    τ_rayl = [zeros(params.float_type,length(params.spec_bands[i]), length(profile.p_full)) for i=1:n_bands];
    #τ_rayl = [zeros(params.float_type,1,length(profile.T)) for i=1:n_bands];
    
    # This is a kludge for now, tau_abs sometimes needs to be a dual. Suniti & us need to rethink this all!!
    # i.e. code the rt core with fixed amount of derivatives as in her paper, then compute chain rule for dtau/dVMr, etc...
    FT2 = isnothing(params.absorption_params) || !haskey(params.absorption_params.vmr,"CO2") ? params.float_type : eltype(params.absorption_params.vmr["CO2"])
    τ_abs     = [zeros(FT2, length(params.spec_bands[i]), length(profile.p_full)) for i in 1:n_bands]
    # Define N_fix_gas as the number of fixed abundance gases (like O2, N2, etc.) - these will not be included in the computation of the Jacobian matrix
    N_fix_gas = length(unique(Iterators.flatten(params.absorption_params.fixed_molecules)))
    # Define N_var_gas as the number of variable gases whose abundance is to be determined - these will be included in the Jacobian computation
    N_var_gas = length(unique(Iterators.flatten(params.absorption_params.variable_molecules)))
    # dertivatives with respect to gas concentration per layer
    τ̇_abs     = [zeros(FT2, 1+N_var_gas, length(params.spec_bands[i]), length(profile.p_full)) for i in 1:n_bands] # reserve the first derivative for H2O. When q is not provided, this derivative can be set to zero 
    max_m = zeros(Int, n_bands)
    l_max = zeros(Int, n_bands)
    l_max_aer = zeros(Int, n_aer, n_bands)
    # Loop over all bands:
    for i_band=1:n_bands

        # i'th spectral band (convert from cm⁻¹ to μm)
        curr_band_λ = params.float_type(1e4) ./ params.spec_bands[i_band]
        # Compute Rayleigh properties per layer for `i_band` band center  
        # It has been added to make sure that the code sees the same Rayleigh cross section, regardless of elastic or inelastic RT 
        # Needs better (more general) formulation 
        νₘ = FT(0.5)*(params.spec_bands[i_band][1]+params.spec_bands[i_band][end])
        λₘ = FT(1.e7)/νₘ  # ← FIX: Use params.float_type
        ϖ_Cabannes[i_band], γ_air_Cabannes, γ_air_Rayleigh = 
            InelasticScattering.compute_γ_air_Rayleigh!(λₘ)
        depol_air_Cabannes = 2γ_air_Cabannes/(1+γ_air_Cabannes)
        depol_air_Rayleigh = 2γ_air_Rayleigh/(1+γ_air_Rayleigh)
        
        a = Scattering.get_greek_rayleigh(depol_air_Cabannes)
        push!(greek_cabannes, a)
        a = Scattering.get_greek_rayleigh(depol_air_Rayleigh)
        push!(greek_rayleigh, a)

        # Use params.depol for consistency with forward model_from_parameters
        τ_rayl[i_band]   .= getRayleighLayerOptProp(profile.p_half[end], 
                                curr_band_λ, 
                                params.depol, profile.vcd_dry);
        

        # If no absorption, continue to next band
        (isnothing(params.absorption_params) && isnothing(params.q)) && continue

        if !isnothing(params.q)
            # Obtain hitran data for H₂O
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
            # TODO: write the following linearized function
            jac_idx = 1 # first index reserved for H2O
            @timeit "Absorption Coeff"  compute_absorption_profile!(
                τ_abs[i_band], 
                τ̇_abs[i_band], 
                jac_idx,
                absorption_model, 
                params.spec_bands[i_band],
                profile.vmr_h2o, 
                profile);
        end
        # Loop over all molecules in this band, obtain profile for each, and add them up
        if !isnothing(params.absorption_params) 
            if !isempty(params.absorption_params.fixed_molecules[i_band])
                # fixed abundance molecules
                for molec_i in 1:length(params.absorption_params.fixed_molecules[i_band])
                    
                    # This can be precomputed as well later in my mind, providing an absorption_model or an interpolation_model!
                    if isempty(params.absorption_params.luts)
                        # Obtain hitran data for this molecule
                        @timeit "Read HITRAN" hitran_data = 
                            read_hitran(artifact(params.absorption_params.fixed_molecules[i_band][molec_i]), iso=1)

                        println("Computing profile for $(params.absorption_params.fixed_molecules[i_band][molec_i]) with vmr $(profile.vmr[params.absorption_params.fixed_molecules[i_band][molec_i]]) for band #$(i_band)")
                        # Create absorption model with parameters beforehand now:
                        absorption_model = make_hitran_model(hitran_data, 
                            params.absorption_params.broadening_function, 
                            wing_cutoff = params.absorption_params.wing_cutoff, 
                            CEF = params.absorption_params.CEF, 
                            architecture = params.architecture, 
                            vmr = 0);#mean(profile.vmr[params.absorption_params.molecules[i_band][molec_i]]))
                        # Calculate absorption profile
                        # TODO: write the following linearized function
                        @timeit "Absorption Coeff"  compute_absorption_profile!(
                                τ_abs[i_band], 
                                absorption_model, 
                                params.spec_bands[i_band],
                                profile.vmr[params.absorption_params.fixed_molecules[i_band][molec_i]], 
                                profile);
                    # Use LUT directly
                    else
                        # TODO: write the following linearized function
                        compute_absorption_profile!(τ_abs[i_band], 
                            params.absorption_params.luts[i_band][molec_i], 
                            params.spec_bands[i_band],
                            profile.vmr[params.absorption_params.fixed_molecules[i_band][molec_i]], 
                            profile);
                    end
                end
            end
            # Variable abundance molecules
            if !isempty(params.absorption_params.variable_molecules[i_band])
                for molec_i in 1:length(params.absorption_params.variable_molecules[i_band])
                    jac_idx = molec_i + 1 # jacobian index for variable abundance molecules
                    # This can be precomputed as well later in my mind, providing an absorption_model or an interpolation_model!
                    if isempty(params.absorption_params.luts)
                        # Obtain hitran data for this molecule
                        @timeit "Read HITRAN" hitran_data = 
                            read_hitran(artifact(params.absorption_params.variable_molecules[i_band][molec_i]), iso=1)
                        println("Computing profile for $(params.absorption_params.variable_molecules[i_band][molec_i]) with vmr $(profile.vmr[params.absorption_params.variable_molecules[i_band][molec_i]]) for band #$(i_band)")
        
                        # Create absorption model with parameters beforehand now:
                        absorption_model = make_hitran_model(hitran_data, 
                            params.absorption_params.broadening_function, 
                            wing_cutoff = params.absorption_params.wing_cutoff, 
                            CEF = params.absorption_params.CEF, 
                            architecture = params.architecture, 
                            vmr = 0);#mean(profile.vmr[params.absorption_params.molecules[i_band][molec_i]]))
                        # Calculate absorption profile
                        # TODO: write the following linearized function
                        @timeit "Absorption Coeff"  compute_absorption_profile!(
                                τ_abs[i_band], 
                                τ̇_abs[i_band], 
                                jac_idx,
                                absorption_model, 
                                params.spec_bands[i_band],
                                profile.vmr[params.absorption_params.variable_molecules[i_band][molec_i]], 
                                profile);
                    # Use LUT directly
                    else
                        # TODO: write the following linearized function
                        # The following form of the function does not seem to exist. May still need to write it (also the forward model version)
                        compute_absorption_profile!(
                            τ_abs[i_band], 
                            τ̇_abs[i_band], 
                            jac_idx,
                            params.absorption_params.luts[i_band][molec_i], 
                            params.spec_bands[i_band],
                            profile.vmr[params.absorption_params.variable_molecules[i_band][molec_i]], 
                            profile);
                    end
                end
            end
        end
    end

    # aerosol_optics[iBand][iAer]
    aerosol_optics = [Array{AerosolOptics}(undef, (n_aer)) for i=1:n_bands];
    lin_aerosol_optics = [Array{linAerosolOptics}(undef, (n_aer)) for i=1:n_bands];

    FT2 = isnothing(params.scattering_params) ? params.float_type : typeof(params.scattering_params.rt_aerosols[1].τ_ref)
    FT2 =  params.float_type 

    # τ_aer[iBand][iAer,iZ]        
    τ_aer = [zeros(FT2, n_aer, length(params.spec_bands[i]), length(profile.p_full)) for i=1:n_bands];
    τ̇_aer = [zeros(FT2, n_aer, 7, length(params.spec_bands[i]), length(profile.p_full)) for i=1:n_bands];
    # the 7 aerosol derivative parameters: τ_ref, nᵣ, nᵢ, rₚ, σₚ, p₀, σp
    # Loop over aerosol type
    for i_aer=1:n_aer

        # Get curr_aerosol
        c_aero = params.scattering_params.rt_aerosols[i_aer]
        curr_aerosol = c_aero.aerosol
        
        # Create Aerosol size distribution for each aerosol species
        size_distribution = curr_aerosol.size_distribution

        # Create a univariate aerosol distribution
        mie_aerosol = Aerosol(size_distribution, curr_aerosol.nᵣ, curr_aerosol.nᵢ)

        # Create the aerosol extinction cross-section at the reference wavelength:
        mie_model      = make_mie_model(params.scattering_params.decomp_type, 
                                        mie_aerosol, 
                                        params.scattering_params.λ_ref, 
                                        params.polarization_type, 
                                        truncation_type, 
                                        params.scattering_params.r_max, 
                                        params.scattering_params.nquad_radius)       
        k_ref, k̇_ref   = compute_ref_aerosol_extinction(lin, mie_model, params.float_type)

        #params.scattering_params.rt_aerosols[i_aer].p₀, params.scattering_params.rt_aerosols[i_aer].σp
        # Loop over bands
        for i_band=1:n_bands
            
            # i'th spectral band (convert from cm⁻¹ to μm)
            curr_band_λ = params.float_type(1e4) ./ params.spec_bands[i_band]

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
                @timeit "Mie calc"  aerosol_optics_raw, lin_aerosol_optics_raw = 
                                compute_aerosol_optical_properties(lin, mie_model, FT2);

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
                @timeit "Mie calc"  aerosol_optics_raw_0, lin_aerosol_optics_raw_0 = 
                                compute_aerosol_optical_properties(lin, mie_model_0, FT2);

                mie_model_1    = make_mie_model(params.scattering_params.decomp_type, 
                                            mie_aerosol, 
                                            curr_band_λ[end], 
                                            params.polarization_type, 
                                            truncation_type, 
                                            params.scattering_params.r_max, 
                                            params.scattering_params.nquad_radius)
                # Compute raw (not truncated) aerosol optical properties (not needed in RT eventually) 
                @timeit "Mie calc"  aerosol_optics_raw_1, lin_aerosol_optics_raw_1 = 
                                compute_aerosol_optical_properties(lin, mie_model_1, FT2);
                Nl = length(aerosol_optics_raw_1.greek_coefs.α)
                Nl_ = length(aerosol_optics_raw_0.greek_coefs.α)
                α = zeros(Nl)
                β = zeros(Nl)
                γ = zeros(Nl)
                δ = zeros(Nl)
                ϵ = zeros(Nl)
                ζ = zeros(Nl)

                α̇ = zeros(4,Nl)
                β̇ = zeros(4,Nl)
                γ̇ = zeros(4,Nl)
                δ̇ = zeros(4,Nl)
                ϵ̇ = zeros(4,Nl)
                ζ̇ = zeros(4,Nl)
                
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

                α̇[:,1:Nl_] = 0.5*(lin_aerosol_optics_raw_0.lin_greek_coefs.α̇[:,1:Nl_] + lin_aerosol_optics_raw_1.lin_greek_coefs.α̇[:,1:Nl_])
                β̇[:,1:Nl_] = 0.5*(lin_aerosol_optics_raw_0.lin_greek_coefs.β̇[:,1:Nl_] + lin_aerosol_optics_raw_1.lin_greek_coefs.β̇[:,1:Nl_])
                γ̇[:,1:Nl_] = 0.5*(lin_aerosol_optics_raw_0.lin_greek_coefs.γ̇[:,1:Nl_] + lin_aerosol_optics_raw_1.lin_greek_coefs.γ̇[:,1:Nl_])
                δ̇[:,1:Nl_] = 0.5*(lin_aerosol_optics_raw_0.lin_greek_coefs.δ̇[:,1:Nl_] + lin_aerosol_optics_raw_1.lin_greek_coefs.δ̇[:,1:Nl_])
                ϵ̇[:,1:Nl_] = 0.5*(lin_aerosol_optics_raw_0.lin_greek_coefs.ϵ̇[:,1:Nl_] + lin_aerosol_optics_raw_1.lin_greek_coefs.ϵ̇[:,1:Nl_])
                ζ̇[:,1:Nl_] = 0.5*(lin_aerosol_optics_raw_0.lin_greek_coefs.ζ̇[:,1:Nl_] + lin_aerosol_optics_raw_1.lin_greek_coefs.ζ̇[:,1:Nl_])
            
                α̇[:,1+Nl_:Nl] = 0.5*(lin_aerosol_optics_raw_1.lin_greek_coefs.α̇[:,1+Nl_:Nl])
                β̇[:,1+Nl_:Nl] = 0.5*(lin_aerosol_optics_raw_1.lin_greek_coefs.β̇[:,1+Nl_:Nl])
                γ̇[:,1+Nl_:Nl] = 0.5*(lin_aerosol_optics_raw_1.lin_greek_coefs.γ̇[:,1+Nl_:Nl])
                δ̇[:,1+Nl_:Nl] = 0.5*(lin_aerosol_optics_raw_1.lin_greek_coefs.δ̇[:,1+Nl_:Nl])
                ϵ̇[:,1+Nl_:Nl] = 0.5*(lin_aerosol_optics_raw_1.lin_greek_coefs.ϵ̇[:,1+Nl_:Nl])
                ζ̇[:,1+Nl_:Nl] = 0.5*(lin_aerosol_optics_raw_1.lin_greek_coefs.ζ̇[:,1+Nl_:Nl])
                
                greek_coeffs = GreekCoefs(α, β, γ, δ, ϵ, ζ)
                lin_greek_coeffs = linGreekCoefs(α̇, β̇, γ̇, δ̇, ϵ̇, ζ̇)
                ν_grid = [1e4/curr_band_λ[1], 1e4/curr_band_λ[end]] 
                kext_grid = [aerosol_optics_raw_0.k, aerosol_optics_raw_1.k]
                ksca_grid = [aerosol_optics_raw_0.k*aerosol_optics_raw_0.ω̃, aerosol_optics_raw_1.k*aerosol_optics_raw_1.ω̃] 
                interp_linear_kext = LinearInterpolation(ν_grid, kext_grid)
                interp_linear_ksca = LinearInterpolation(ν_grid, ksca_grid)
                k = zeros(length(curr_band_λ))
                ω̃ = zeros(length(curr_band_λ))
                fᵗ = zeros(length(curr_band_λ))
                for i = 1:length(curr_band_λ)
                    k[i] = interp_linear_kext(1e4/curr_band_λ[i])
                    ω̃[i] = interp_linear_ksca(1e4/curr_band_λ[i])/k[i]
                end
                k̇ = zeros(4,length(curr_band_λ))
                ω̃̇ = zeros(4,length(curr_band_λ))
                ḟᵗ= zeros(4,length(curr_band_λ))
                for ctr=1:4
                    k̇ext_grid = [lin_aerosol_optics_raw_0.k̇[ctr], lin_aerosol_optics_raw_1.k̇[ctr]]
                    # Bug 20 fix: use product rule d(k*ω̃)/dp = dk/dp*ω̃ + k*dω̃/dp
                    # (was: k̇*ω̃̇, i.e. product of two derivatives — incorrect)
                    k̇sca_grid = [lin_aerosol_optics_raw_0.k̇[ctr] * aerosol_optics_raw_0.ω̃ +
                                  aerosol_optics_raw_0.k * lin_aerosol_optics_raw_0.ω̃̇[ctr], 
                                  lin_aerosol_optics_raw_1.k̇[ctr] * aerosol_optics_raw_1.ω̃ +
                                  aerosol_optics_raw_1.k * lin_aerosol_optics_raw_1.ω̃̇[ctr]] 
                    interp_linear_k̇ext = LinearInterpolation(ν_grid, k̇ext_grid)
                    interp_linear_k̇sca = LinearInterpolation(ν_grid, k̇sca_grid)
                    
                    for i = 1:length(curr_band_λ)
                        k̇[ctr,i] = interp_linear_k̇ext(1e4/curr_band_λ[i])
                        # Bug 20 fix: use quotient rule d(ksca/kext)/dp = (dk_sca/dp - ω̃*dk_ext/dp)/k_ext
                        # (was: k̇sca/k̇ext — divides by derivative, blows up when k̇≈0)
                        ω̃̇[ctr,i] = (interp_linear_k̇sca(1e4/curr_band_λ[i]) - ω̃[i] * k̇[ctr,i]) / k[i]
                    end
                end

                aerosol_optics_raw = AerosolOptics(greek_coefs=greek_coeffs, ω̃=ω̃, k=k, fᵗ=fᵗ)
                lin_aerosol_optics_raw = linAerosolOptics(lin_greek_coefs=lin_greek_coeffs, ω̃̇=ω̃̇, k̇=k̇, ḟᵗ=ḟᵗ)
            end
            
            # Compute truncated aerosol optical properties (phase function and fᵗ), consistent with Ltrunc:
            
            if length(aerosol_optics_raw.greek_coefs.β) > truncation_type.l_max
                aerosol_optics[i_band][i_aer], lin_aerosol_optics[i_band][i_aer] = 
                    Scattering.truncate_phase(truncation_type, 
                                aerosol_optics_raw, lin_aerosol_optics_raw; reportFit=false)
                l_max_aer[i_aer, i_band] = truncation_type.l_max
                
                #max_m[i_band] = Int(ceil(l_max[i_band] + 1)/2)                                 
            else
                aerosol_optics[i_band][i_aer] = aerosol_optics_raw
                lin_aerosol_optics[i_band][i_aer] = lin_aerosol_optics_raw
                l_max_aer[i_aer, i_band] = length(aerosol_optics_raw.greek_coefs.β)
                #max_m[i_band] = Int(ceil(l_max[i_band] + 1)/2)
            end
            # Extract p₀ and σp from the aerosol's vertical distribution (Normal(p₀, σp))
            aer_p₀ = mean(c_aero.profile)
            aer_σp = std(c_aero.profile)
            # Compute aerosol optical thickness profile and derivatives w.r.t. p₀ and σp
            τₚ, dτₚdp₀, dτₚdσp = getAerosolLayerOptProp(lin, 1, aer_p₀, aer_σp, profile.p_half)

            # Compute nAer aerosol optical thickness profiles
            τ_aer[i_band][i_aer,:,:] = 
                (params.scattering_params.rt_aerosols[i_aer].τ_ref/k_ref) * 
                aerosol_optics[i_band][i_aer].k *
                τₚ' 

            # Derivative w.r.t. τ_ref (parameter 1)
            τ̇_aer[i_band][i_aer,1,:,:] .= 
                (aerosol_optics[i_band][i_aer].k/k_ref) *
                τₚ'  
                
            # Derivatives w.r.t. nᵣ, nᵢ, rₚ, σₚ (parameters 2-5, via linearized Mie)
            # τ_aer = (τ_ref/k_ref) * k(λ) * τₚ  ⟹  dτ_aer/dp = τ_ref * d(k/k_ref)/dp * τₚ
            #       = (τ_ref/k_ref * dk/dp - τ_ref * k * dk_ref/dp / k_ref²) * τₚ
            for ctr=1:4
                τ̇_aer[i_band][i_aer,ctr+1,:,:] = 
                    ((params.scattering_params.rt_aerosols[i_aer].τ_ref/k_ref) * 
                    lin_aerosol_optics[i_band][i_aer].k̇[ctr,:] .- 
                    # Bug 21 fix: was missing k_band factor (aerosol_optics.k) in the k_ref term
                    (params.scattering_params.rt_aerosols[i_aer].τ_ref/k_ref^2) * 
                    k̇_ref[ctr] .* aerosol_optics[i_band][i_aer].k)*τₚ'
            end
            # Derivatives w.r.t. p₀ and σp (parameters 6-7, via linearized vertical profile)
            for ctr=5:6
                τ̇_aer[i_band][i_aer,ctr+1,:,:] = 
                    (params.scattering_params.rt_aerosols[i_aer].τ_ref/k_ref) * 
                    aerosol_optics[i_band][i_aer].k *
                    (ctr==5 ? dτₚdp₀' : dτₚdσp')
            end
                             
        end                  
    end
    for i_band = 1:n_bands
        l_max[i_band] = maximum(l_max_aer[:,i_band]) 
        max_m[i_band] = Int(ceil(l_max[i_band] + 1)/2)  
    end
    set_uniform_lmax!(l_max, aerosol_optics)

    # Check the floating-type output matches specified FT
    # Return the model 
    return  vSmartMOM_Model(
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
                        profile), 
            vSmartMOM_Lin(τ̇_abs, τ̇_aer, lin_aerosol_optics)
end

