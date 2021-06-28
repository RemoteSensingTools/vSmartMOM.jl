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
function validate_aerosols(aerosols::Array{Dict{Any, Any}})
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

"Check that the vmr's in the atmospheric profile match the molecules in the parameters" 
function validate_vmrs(params::vSmartMOM_Parameters, profile::AtmosphericProfile)
    for molec in unique(vcat(params.molecules...))
        @assert molec in keys(profile.vmr) "$(molec) listed as molecule in parameters yaml, but no vmr given in atmospheric profile"
        @assert profile.vmr[molec] isa Real || profile.vmr[molec] isa Vector "The vmr for $(molec) in the atmospheric profile must either be a real-valued number, or an array of nodal points from surface to 0hPa (TOA)"
    end
end

"Given a parameter dictionary from a YAML file, validate the dictionary"
function validate_yaml_parameters(params)

    # All the fields that *should* exist and their types
    fields = [(["absorption", "broadening"], String, ["Voigt()", "Doppler()", "Lorentz()"]), 
              (["absorption", "CEF"], String, ["HumlicekWeidemann32SDErrorFunction()"]),
              (["absorption", "wing_cutoff"], Integer),
              (["absorption", "spec_bands"], Array{String}),
              (["absorption", "molecules"], Array{Array{String, 1}, 1}),
              
              (["scattering", "aerosols"], Array{Dict{Any,Any}}),
              (["scattering", "r_max"], Real),
              (["scattering", "nquad_radius"], Real),
              (["scattering", "λ"], Array{<:Real}),
              (["scattering", "λ_ref"], Real),
              (["scattering", "depol"], Real),
              (["scattering", "polarization_type"], String, ["Stokes_I()", "Stokes_IQ()", "Stokes_IQU()", "Stokes_IQUV()"]),
              (["scattering", "decomp_type"], String, ["NAI2()", "PCW()"]),

              (["geometry", "vza"], Array{<:Real}),
              (["geometry", "vaz"], Array{<:Real}),
              (["geometry", "sza"], Real),
              (["geometry", "obs_alt"], Real),

              (["atmospheric_profile", "file"], String),
              (["atmospheric_profile", "profile_reduction"], Integer),

              (["surface", "BRDFs"], Array{String}),

              (["truncation_type", "l_trunc"], Integer),
              (["truncation_type", "Δ_angle"], Real),

              (["vSmartMOM", "quadrature_type"], String),
              (["vSmartMOM", "max_m"], Integer),
              (["vSmartMOM", "architecture"], String),
              (["vSmartMOM", "float_type"], String)
              ]

    # Check that every field exists and is correctly typed
    for i in 1:length(fields)
        key, elty = fields[i][1:2]
        valid_options = length(fields[i]) == 3 ? fields[i][3] : Array{String}([])
        @assert check_yaml_field(params, key, key, elty, valid_options) 
    end

    # Check the aerosols separately
    validate_aerosols(params["scattering"]["aerosols"])
end

"Given a path to a YAML parameters file, load all the parameters into a vSmartMOM_Parameters struct" 
function parameters_from_yaml(file_path)

    # Load the YAML file
    params_dict = YAML.load_file(file_path)

    # Validate to make sure it has all necessary fields
    validate_yaml_parameters(params_dict)

    # Evaluate certain fields to get their proper types
    broadening_function = eval(Meta.parse(params_dict["absorption"]["broadening"]))
    CEF = eval(Meta.parse(params_dict["absorption"]["CEF"]))
    polarization_type = eval(Meta.parse(params_dict["scattering"]["polarization_type"]))
    decomp_type = eval(Meta.parse(params_dict["scattering"]["decomp_type"]))
    quadrature_type = eval(Meta.parse(params_dict["vSmartMOM"]["quadrature_type"]))
    architecture = eval(Meta.parse(params_dict["vSmartMOM"]["architecture"]))
    FT = eval(Meta.parse(params_dict["vSmartMOM"]["float_type"]))
    BRDF_per_band = map(x -> eval(Meta.parse(x)), params_dict["surface"]["BRDFs"]) 
    lengths = convert.(Integer, map(x -> length(collect(eval(Meta.parse(x)))), params_dict["absorption"]["spec_bands"]))

    # Evaluate profile path if it's not a plain string (need to joinpath)
    if occursin("joinpath", params_dict["atmospheric_profile"]["file"])
        profile_path = eval(Meta.parse(params_dict["atmospheric_profile"]["file"]))
    else
        profile_path = params_dict["atmospheric_profile"]["file"]
    end
    
    # If the specified atmospheric profile is YAML, then set time_index/lat/lon to nothing
    if endswith(profile_path, ".yaml")
        time_index = nothing
        lat = nothing
        lon = nothing

    # Otherwise, store them
    else
        time_index = params_dict["atmospheric_profile"]["time_index"]
        lat = params_dict["atmospheric_profile"]["lat"]
        lon = params_dict["atmospheric_profile"]["lon"]
    end

    return vSmartMOM_Parameters(FT,
                                convert.(FT, params_dict["scattering"]["λ"]),
                                FT(params_dict["scattering"]["λ_ref"]),
                                FT(params_dict["scattering"]["depol"]),
                                params_dict["truncation_type"]["l_trunc"],
                                FT(params_dict["truncation_type"]["Δ_angle"]),
                                params_dict["vSmartMOM"]["max_m"],
                                polarization_type,
                                FT(params_dict["geometry"]["obs_alt"]),
                                FT(params_dict["geometry"]["sza"]),
                                convert.(FT, params_dict["geometry"]["vza"]),
                                convert.(FT, params_dict["geometry"]["vaz"]),
                                length(params_dict["scattering"]["aerosols"]), 
                                convert.(FT, map(x -> x["τ_ref"], params_dict["scattering"]["aerosols"])),
                                convert.(FT, map(x -> x["μ"], params_dict["scattering"]["aerosols"])),
                                convert.(FT, map(x -> x["σ"], params_dict["scattering"]["aerosols"])),
                                FT(params_dict["scattering"]["r_max"]),
                                params_dict["scattering"]["nquad_radius"],
                                convert.(FT, map(x -> x["nᵣ"], params_dict["scattering"]["aerosols"])),
                                convert.(FT, map(x -> x["nᵢ"], params_dict["scattering"]["aerosols"])),
                                convert.(FT, map(x -> x["p₀"], params_dict["scattering"]["aerosols"])),
                                convert.(FT, map(x -> x["σp"], params_dict["scattering"]["aerosols"])),
                                profile_path,
                                time_index,
                                lat,
                                lon,
                                params_dict["atmospheric_profile"]["profile_reduction"],
                                convert.(FT, map(x -> first(collect(eval(Meta.parse(x)))), params_dict["absorption"]["spec_bands"])),
                                convert.(FT, map(x -> last(collect(eval(Meta.parse(x)))), params_dict["absorption"]["spec_bands"])),
                                convert(Array{Integer}, convert.(Integer, map(x -> length(collect(eval(Meta.parse(x)))), params_dict["absorption"]["spec_bands"]))),
                                params_dict["absorption"]["molecules"],
                                broadening_function,
                                params_dict["absorption"]["wing_cutoff"],
                                CEF,
                                BRDF_per_band,
                                decomp_type,
                                quadrature_type,
                                architecture)
end

"Take the parameters specified in the vSmartMOM_Parameters struct, and calculate derived attributes into a vSmartMOM_Model" 
function model_from_parameters(params::vSmartMOM_Parameters)

    @unpack λ_band = params
    nBands = length(λ_band)

    # Create observation geometry
    obs_geom = ObsGeometry(params.obs_alt, params.sza, params.vza, params.vaz)

    # Create truncation type
    truncation_type = Scattering.δBGE{params.float_type}(params.l_trunc, params.Δ_angle)

    # Set quadrature points
    quad_points = rt_set_streams(params.quadrature_type, params.l_trunc, obs_geom, params.polarization_type, array_type(params.architecture))

    # Read profile (and generate dry/wet VCDs per layer)
    if isnothing(params.timeIndex)
        profile_hr = vSmartMOM.read_atmos_profile(params.file);
    else
        profile_hr = vSmartMOM.read_atmos_profile(params.file, params.lat, params.lon, params.timeIndex);
    end

    # Validate that the vmr's in the atmospheric profile, match those in the parameters
    validate_vmrs(params, profile_hr)
    
    # Reduce the profile to the number of target layers
    profile = vSmartMOM.reduce_profile(params.profile_reduction_n, profile_hr);

    #Rayleigh
    greek_rayleigh = Scattering.get_greek_rayleigh(params.depol)
    τRayl = [zeros(params.float_type, length(profile.p)) for i=1:nBands];

    # Create empty array (can be changed later)
    # spec_grid[iBand][iSpec]
    spec_grid = [zeros(params.float_type,n) for n in params.spec_grid_n];

    # τ_abs[iBand][iSpec,iZ]
    τ_abs     = [zeros(params.float_type, n, length(profile.p)) for n in params.spec_grid_n]
    
    # Loop over all bands:
    for ib=1:length(λ_band)

        # Compute Rayleigh properties per layer for `ib` band center 
        τRayl[ib]   = getRayleighLayerOptProp(profile.psurf / 100, λ_band[ib], params.depol, profile.vcd_dry);

        # Define spectral grid per band
        spec_grid[ib] = range(params.spec_grid_start[ib], params.spec_grid_end[ib], length=params.spec_grid_n[ib]);
        
        # Loop over all molecules in this band, obtain profile for each, and add them up
        for molec_i in 1:length(params.molecules[ib])

            # Obtain hitran data for this molecule
            hitran_data = read_hitran(artifact(params.molecules[ib][molec_i]), iso=1)

            println("Computing profile for $(params.molecules[ib][molec_i]) with vmr $(profile_hr.vmr[params.molecules[ib][molec_i]]) for band #$(ib)")

            # Calculate absorption profile
            compute_absorption_profile!(τ_abs[ib], hitran_data, params.broadening_function, params.wing_cutoff, params.CEF, params.architecture, profile.vmr[params.molecules[ib][molec_i]], spec_grid[ib], profile);
        end
    end

    # aerosol_optics[iBand][iAer]
    aerosol_optics = [Array{AerosolOptics}(undef, (params.nAer)) for i=1:nBands];

    FT2 = typeof(params.nᵢ[1])

    # τAer[iBand][iAer,iZ]
    τAer = [zeros(FT2, params.nAer, length(profile.p)) for i=1:nBands];

    # Loop over aerosol type
    for iaer=1:params.nAer

        # Create Aerosol size distribution for each aerosol species
        size_distribution = LogNormal(log(params.μ[iaer]), log(params.σ[iaer]))

        # Create a univariate aerosol distribution
        aerosol = make_univariate_aerosol(size_distribution, params.r_max, params.nquad_radius, params.nᵣ[iaer], params.nᵢ[iaer]) #Suniti: why is the refractive index needed here?

        # Create the aerosol extinction cross-section at the reference wavelength:
        mie_model      = make_mie_model(params.decomp_type, aerosol, params.λ_ref, params.polarization_type, truncation_type)       
        k_ref          = compute_ref_aerosol_extinction(mie_model, params.float_type)

        # Loop over bands
        for ib=1:nBands

            # Create the aerosols:
            mie_model      = make_mie_model(params.decomp_type, aerosol, params.λ_band[ib], params.polarization_type, truncation_type)

            # Compute raw (not truncated) aerosol optical properties (not needed in RT eventually) 
            aerosol_optics_raw = compute_aerosol_optical_properties(mie_model, FT2);

            # Compute truncated aerosol optical properties (phase function and fᵗ), consistent with Ltrunc:
            @show iaer,ib
            aerosol_optics[ib][iaer] = Scattering.truncate_phase(truncation_type, aerosol_optics_raw; reportFit=false)

            # Compute nAer aerosol optical thickness profiles
            τAer[ib][iaer,:] = params.τAer_ref[iaer] * (aerosol_optics[ib][iaer].k/k_ref) * vSmartMOM.getAerosolLayerOptProp(1.0, params.p₀[iaer], params.σp[iaer], profile.p)
            
        end 
    end
    
    # Return the model 
    return vSmartMOM_Model(params, aerosol_optics,  greek_rayleigh, quad_points, τ_abs, τRayl, τAer, obs_geom, profile)

end