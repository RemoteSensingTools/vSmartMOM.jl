"Generate default set of parameters for Radiative Transfer calculations (from ModelParameters/)"
default_parameters() = parameters_from_yaml(joinpath(dirname(pathof(RadiativeTransfer)), "vSmartMOM", "ModelParameters", "DefaultParameters.yaml"))

# Check that a field exists in yaml file
function check_yaml_field(dict::Dict{Any, Any}, full_keys::Array{String}, curr_keys::Array{String}, final_type::Type, valid_options::Array{String})

    # Should have at least one key
    @assert length(curr_keys) >= 1

    # Only one key (that should exist)
    if length(curr_keys) == 1
        @assert curr_keys[1] in keys(dict) "Missing key in parameters yaml: $(join(full_keys, '/'))"
        @assert dict[curr_keys[1]] isa final_type "Improper type for $(join(full_keys, '/')); must be a $(final_type)" 
        if length(valid_options) > 0
            @assert dict[curr_keys[1]] in valid_options "Field $(join(full_keys, '/')) must be one of $(valid_options)"
        end
        return true
    
    # Recursive case, check sub-dictionary
    else
        @assert curr_keys[1] in keys(dict) "Missing key in parameters yaml: $(join(full_keys, '/'))"
        return check_yaml_field(dict[curr_keys[1]], full_keys, curr_keys[2:end], final_type, valid_options)
    end
end

# Given a list of aerosol parameter dictionaries, validate all of them
function validate_aerosols(aerosols::Union{Array{Dict{Any, Any}}, Vector{Any}})

    @assert !isempty(aerosols) "Aerosols list shouldn't be empty if scattering block is included"

    fields = [(["τ_ref"], Real),
              (["μ"], Real),
              (["σ"], Real),
              (["nᵣ"], Real),
              (["nᵢ"], Real),
              (["p₀"], Real),
              (["p₀"], Real)]
              
    for aerosol in aerosols
        for i in 1:length(fields)
            key, elty = fields[i][1:2]
            @assert check_yaml_field(aerosol, key, key, elty, Array{String}([])) 
        end
    end
end

function aerosol_params_to_obj(aerosols::Union{Array{Dict{Any, Any}}, Vector{Any}}, FT)
    # return Aerosol[]
    aerosol_obj_list = Aerosol[]
    for aerosol in aerosols
        new_aerosol_obj = Aerosol(FT(aerosol["τ_ref"]),
                                  FT(aerosol["μ"]), 
                                  FT(aerosol["σ"]),
                                  FT(aerosol["nᵣ"]),
                                  FT(aerosol["nᵢ"]),
                                  FT(aerosol["p₀"]),
                                  FT(aerosol["σp"]))
        push!(aerosol_obj_list, new_aerosol_obj)
    end

    return aerosol_obj_list
end

"Check that the vmr's in the atmospheric profile match the molecules in the parameters" 
function validate_vmrs(molecules::Array, vmr::Dict)
    for molec in unique(vcat(molecules...))
        @assert molec in keys(vmr) "$(molec) listed as molecule in parameters yaml, but no vmr given in atmospheric profile"
        @assert vmr[molec] isa Real || vmr[molec] isa Vector "The vmr for $(molec) in the atmospheric profile must either be a real-valued number, or an array of nodal points from surface to 0hPa (TOA)"
    end
end

"Given a parameter dictionary from a YAML file, validate the dictionary"
function validate_yaml_parameters(params)

    # All the fields that *should* exist and their types
    fields = [

              # radiative_transfer group
              (["radiative_transfer", "spec_bands"], Array{String}),  
              (["radiative_transfer", "surface"], Array{String}),   
              (["radiative_transfer", "quadrature_type"], String), 
              (["radiative_transfer", "polarization_type"], String),   
              (["radiative_transfer", "max_m"], Integer), 
              (["radiative_transfer", "Δ_angle"], Real),
              (["radiative_transfer", "l_trunc"], Integer),
              (["radiative_transfer", "depol"], Real),
              (["radiative_transfer", "float_type"], String),
              (["radiative_transfer", "architecture"], String),

              # geometry group
              (["geometry", "sza"], Real),
              (["geometry", "vza"], Array{<:Real}),
              (["geometry", "vaz"], Array{<:Real}),
              (["geometry", "obs_alt"], Real),

              # atmospheric_profile group
              (["atmospheric_profile", "T"], Array{<:Real}),
              (["atmospheric_profile", "p"], Array{<:Real}),
              (["atmospheric_profile", "profile_reduction"], Union{Integer, Nothing}),

              # absorption group
              (["absorption", "molecules"], Array),
              (["absorption", "vmr"], Dict),
              (["absorption", "broadening"], String, ["Voigt()", "Doppler()", "Lorentz()"]), 
              (["absorption", "CEF"], String, ["HumlicekWeidemann32SDErrorFunction()"]),
              (["absorption", "wing_cutoff"], Integer),
              
              # scattering group
              (["scattering", "aerosols"], Array),
              (["scattering", "r_max"], Real),
              (["scattering", "nquad_radius"], Real),
              (["scattering", "λ_ref"], Real),
              (["scattering", "decomp_type"], String, ["NAI2()", "PCW()"]),
              ]

    # Check that every field exists and is correctly typed
    for i in 1:length(fields)
        key, elty = fields[i][1:2]
        valid_options = length(fields[i]) == 3 ? fields[i][3] : Array{String}([])

        # Check if the section exists in the YAML file and is an optional section
        if ("absorption" in keys(params) && key[1] == "absorption") || 
           ("scattering" in keys(params) && key[1] == "scattering") || 
           (key[1] != "absorption" && key[1] != "scattering")
            @assert check_yaml_field(params, key, key, elty, valid_options) 
        end
    end

    # Check the aerosols separately if specified
    if "scattering" in keys(params)
        validate_aerosols(params["scattering"]["aerosols"])
    end
end

"Given a path to a YAML parameters file, load all the parameters into a vSmartMOM_Parameters struct" 
function parameters_from_yaml(file_path)

    # Load the YAML file
    params_dict = YAML.load_file(file_path)

    # Validate to make sure it has all necessary fields
    validate_yaml_parameters(params_dict)

    # #########################################################
    # Evaluate type-specifying fields to get their proper types
    # #########################################################

    # radiative_transfer group
    spec_bands = convert.(Array, map(x -> collect(eval(Meta.parse(x))), params_dict["radiative_transfer"]["spec_bands"]))
    BRDF_per_band = map(x -> eval(Meta.parse(x)), params_dict["radiative_transfer"]["surface"]) 
    quadrature_type = eval(Meta.parse(params_dict["radiative_transfer"]["quadrature_type"]))
    polarization_type = eval(Meta.parse(params_dict["radiative_transfer"]["polarization_type"]))
    FT = eval(Meta.parse(params_dict["radiative_transfer"]["float_type"]))
    architecture = eval(Meta.parse(params_dict["radiative_transfer"]["architecture"]))

    # atmospheric_profile group
    T = convert.(FT, params_dict["atmospheric_profile"]["T"]) # Level
    p = convert.(FT, params_dict["atmospheric_profile"]["p"]) # Boundaries
    q = "q" in keys(params_dict["atmospheric_profile"]) ? # Specific humidity, if it's specified. 
            convert.(Float64, params_dict["atmospheric_profile"]["q"]) : zeros(length(T)) # Otherwise 0

    # absorption group
    if "absorption" in keys(params_dict)
        molecules = Array(params_dict["absorption"]["molecules"])
        vmr = convert(Dict{String, Union{Real, Vector}}, params_dict["absorption"]["vmr"])
        validate_vmrs(params_dict["absorption"]["molecules"], vmr)
        broadening_function = eval(Meta.parse(params_dict["absorption"]["broadening"]))
        CEF = eval(Meta.parse(params_dict["absorption"]["CEF"]))
        wing_cutoff = params_dict["absorption"]["wing_cutoff"]

        absorption_params = AbsorptionParameters(molecules, vmr, broadening_function, 
                                                 CEF, wing_cutoff)
    else 
        absorption_params = nothing
    end
    
    # scattering group
    if "scattering" in keys(params_dict)
        
        # Get aerosols from types
        aerosols = aerosol_params_to_obj(params_dict["scattering"]["aerosols"], FT)
        r_max = FT(params_dict["scattering"]["r_max"])
        nquad_radius = params_dict["scattering"]["nquad_radius"]
        λ_ref = FT(params_dict["scattering"]["λ_ref"])
        decomp_type = eval(Meta.parse(params_dict["scattering"]["decomp_type"]))

        scattering_params = ScatteringParameters(aerosols, r_max, nquad_radius, 
                                                 λ_ref, decomp_type)
    else
        scattering_params = nothing
    end

    return vSmartMOM_Parameters(
                                # radiative_transfer group
                                spec_bands, 
                                BRDF_per_band,
                                quadrature_type,
                                polarization_type,
                                params_dict["radiative_transfer"]["max_m"],
                                FT(params_dict["radiative_transfer"]["Δ_angle"]),
                                params_dict["radiative_transfer"]["l_trunc"],
                                FT(params_dict["radiative_transfer"]["depol"]),
                                FT, 
                                architecture,
                                
                                # geometry group
                                FT(params_dict["geometry"]["sza"]),
                                convert.(FT, params_dict["geometry"]["vza"]),
                                convert.(FT, params_dict["geometry"]["vaz"]),
                                FT(params_dict["geometry"]["obs_alt"]),

                                # atmospheric_profile group
                                T, p, q, 
                                params_dict["atmospheric_profile"]["profile_reduction"],

                                # absorption group
                                absorption_params,

                                # scattering group
                                scattering_params
                                );
end

"Take the parameters specified in the vSmartMOM_Parameters struct, and calculate derived attributes into a vSmartMOM_Model" 
function model_from_parameters(params::vSmartMOM_Parameters)

    # Number of total bands and aerosols (for convenience)
    n_bands = length(params.spec_bands)
    n_aer = isnothing(params.scattering_params) ? 0 : length(params.scattering_params.aerosols)

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
        profile = vSmartMOM.reduce_profile(params.profile_reduction_n, profile);
    end

    # Rayleigh optical properties calculation
    greek_rayleigh = Scattering.get_greek_rayleigh(params.depol)
    τ_rayl = [zeros(params.float_type, length(params.T)) for i=1:n_bands];

    # τ_abs[iBand][iSpec,iZ]
    τ_abs     = [zeros(params.float_type, length(params.spec_bands[i]), length(profile.p_full)) for i in 1:n_bands]
    
    # Loop over all bands:
    for i_band=1:n_bands

        # i'th spectral band (convert from cm⁻¹ to μm)
        curr_band_λ = 1e4 ./ params.spec_bands[i_band]

        # Compute Rayleigh properties per layer for `i_band` band center 
        τ_rayl[i_band]   = getRayleighLayerOptProp(profile.p_half[end], (maximum(curr_band_λ) + minimum(curr_band_λ)/2), params.depol, profile.vcd_dry);

        # If no absorption, continue to next band
        isnothing(params.absorption_params) && continue
        
        # Loop over all molecules in this band, obtain profile for each, and add them up
        for molec_i in 1:length(params.absorption_params.molecules[i_band])

            # Obtain hitran data for this molecule
            @timeit "Read HITRAN"  hitran_data = read_hitran(artifact(params.absorption_params.molecules[i_band][molec_i]), iso=1)

            println("Computing profile for $(params.absorption_params.molecules[i_band][molec_i]) with vmr $(profile.vmr[params.absorption_params.molecules[i_band][molec_i]]) for band #$(i_band)")

            # Calculate absorption profile
            @timeit "Absorption Coeff"  compute_absorption_profile!(τ_abs[i_band], hitran_data, params.absorption_params.broadening_function, params.absorption_params.wing_cutoff, params.absorption_params.CEF, params.architecture, profile.vmr[params.absorption_params.molecules[i_band][molec_i]], params.spec_bands[i_band], profile);
        end
    end

    # aerosol_optics[iBand][iAer]
    aerosol_optics = [Array{AerosolOptics}(undef, (n_aer)) for i=1:length(0)];
        
    FT2 = isnothing(params.scattering_params) ? params.float_type : typeof(params.scattering_params.aerosols[1].τ_ref)

    # τ_aer[iBand][iAer,iZ]
    τ_aer = [zeros(FT2, n_aer, length(profile.p_full)) for i=1:n_bands];

    # Loop over aerosol type
    for i_aer=1:n_aer

        # Get curr_aerosol
        curr_aerosol = params.scattering_params.aerosols[i_aer]

        # Create Aerosol size distribution for each aerosol species
        size_distribution = LogNormal(log(curr_aerosol.μ), log(curr_aerosol.σ))

        # Create a univariate aerosol distribution
        aerosol = make_univariate_aerosol(size_distribution, params.scattering_params.r_max, params.scattering_params.nquad_radius, curr_aerosol.nᵣ, curr_aerosol.nᵢ) #Suniti: why is the refractive index needed here?

        # Create the aerosol extinction cross-section at the reference wavelength:
        mie_model      = make_mie_model(params.scattering_params.decomp_type, aerosol, params.scattering_params.λ_ref, params.polarization_type, truncation_type)       
        k_ref          = compute_ref_aerosol_extinction(mie_model, params.float_type)

        # Loop over bands
        for i_band=1:n_bands

            # i'th spectral band (convert from cm⁻¹ to μm)
            curr_band_λ = 1e4 ./ params.spec_bands[i_band]

            # Create the aerosols:
            mie_model      = make_mie_model(params.scattering_params.decomp_type, aerosol, (maximum(curr_band_λ)+minimum(curr_band_λ))/2, params.polarization_type, truncation_type)

            # Compute raw (not truncated) aerosol optical properties (not needed in RT eventually) 
            # @show FT2
            @timeit "Mie calc"  aerosol_optics_raw = compute_aerosol_optical_properties(mie_model, Float64);

            # Compute truncated aerosol optical properties (phase function and fᵗ), consistent with Ltrunc:
            @show i_aer, i_band
            aerosol_optics[i_band][i_aer] = Scattering.truncate_phase(truncation_type, aerosol_optics_raw; reportFit=false)

            # Compute nAer aerosol optical thickness profiles
            τ_aer[i_band][i_aer,:] = params.scattering_params.aerosols[i_aer].τ_ref[] * (aerosol_optics[i_band][i_aer].k/k_ref) * vSmartMOM.getAerosolLayerOptProp(1.0, params.scattering_params.aerosols[i_aer].p₀, params.scattering_params.aerosols[i_aer].σp, profile.p_full)
            
        end 
    end

    # Check the floating-type output matches specified FT

    # Return the model 
    return vSmartMOM_Model(params, aerosol_optics,  greek_rayleigh, quad_points, τ_abs, τ_rayl, τ_aer, obs_geom, profile)

end