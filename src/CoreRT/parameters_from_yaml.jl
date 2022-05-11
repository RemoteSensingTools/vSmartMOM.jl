#=

This file contains functions that help convert a parameters YAML file into a vSmartMOM_Parameters
object. This process includes converting data types, specifying default values, validating
the input, and producing the output object. 

=#

"Check that a field exists in yaml file"
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

"Given a list of aerosol parameter dictionaries, validate all of them"
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

"Convert the input dictionary of aerosols into a list of RT_aerosols"
function aerosol_params_to_obj(aerosols::Union{Array{Dict{Any, Any}}, Vector{Any}}, FT)

    rt_aerosol_obj_list = RT_Aerosol[]

    for aerosol in aerosols
        @assert aerosol["σ"] ≥ 1 "Geometric standard deviation has to be ≥ 1"    
        size_distribution = LogNormal(log(FT(aerosol["μ"])), log(FT(aerosol["σ"])))
        #@show size_distribution
        new_aerosol_obj = Aerosol(size_distribution,
                                  FT(aerosol["nᵣ"]),
                                  FT(aerosol["nᵢ"]))
        
        new_rt_aerosol_obj = RT_Aerosol(new_aerosol_obj, FT(aerosol["τ_ref"]), FT(aerosol["p₀"]), FT(aerosol["σp"]))

        push!(rt_aerosol_obj_list, new_rt_aerosol_obj)
    end

    return rt_aerosol_obj_list
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
    FT = eval(Meta.parse(params_dict["radiative_transfer"]["float_type"]))
    # @show params_dict["radiative_transfer"]["float_type"]
    # Each spec band can have units in wavenumber/wavelength, or not. Regardless, convert to cm⁻¹
    spec_bands = []
    for spec_band in params_dict["radiative_transfer"]["spec_bands"]
        parsed_band = eval(Meta.parse(spec_band))
        # If no units, assume cm⁻¹
        if (all(x-> x == NoUnits, unit.(parsed_band)))
            wn_band = collect(parsed_band)
        else
            # If units but cm⁻¹ already, don't uconvert (issue trying to convert to already existing type)
            if (all(x-> x == u"cm^-1", unit.(parsed_band)))
                wn_band = sort(ustrip(collect(parsed_band)))
            # If units but not cm⁻¹, uconvert to cm⁻¹
            else
                wn_band = sort(ustrip(uconvert.( u"cm^-1", collect(parsed_band), Spectral())))
            end
        end
        final_band = convert(Array{FT}, wn_band)
        push!(spec_bands, final_band)
    end
    BRDF_per_band     = map(x -> eval(Meta.parse(x)), params_dict["radiative_transfer"]["surface"]) 
    quadrature_type   = eval(Meta.parse(params_dict["radiative_transfer"]["quadrature_type"]))
    
    # Make type stable, pol type has to be in the right FT:
    pol_type = replace(params_dict["radiative_transfer"]["polarization_type"], "()" => "{$FT}()")
    #@show pol_type
    polarization_type = eval(Meta.parse(pol_type))
    
    architecture      = eval(Meta.parse(params_dict["radiative_transfer"]["architecture"]))
    #@show polarization_type, quadrature_type 
    # atmospheric_profile group
    T = convert.(FT, params_dict["atmospheric_profile"]["T"]) # Level
    p = convert.(FT, params_dict["atmospheric_profile"]["p"]) # Boundaries
    q = "q" in keys(params_dict["atmospheric_profile"]) ? # Specific humidity, if it's specified. 
            convert.(FT, params_dict["atmospheric_profile"]["q"]) : zeros(length(T)) # Otherwise 0

    # absorption group
    if "absorption" in keys(params_dict)
        molecules = Array(params_dict["absorption"]["molecules"])
        #@show params_dict["absorption"]["vmr"]
        vmr = convert(Dict{String, Union{Real, Vector}}, params_dict["absorption"]["vmr"])

        validate_vmrs(params_dict["absorption"]["molecules"], vmr)
        broadening_function = eval(Meta.parse(params_dict["absorption"]["broadening"]))
        CEF = eval(Meta.parse(params_dict["absorption"]["CEF"]))
        wing_cutoff = FT(params_dict["absorption"]["wing_cutoff"])

        
        # Option to load lookup tables!
        luts = []
        if "LUTfiles" in keys(params_dict["absorption"])
            files_lut = Array(params_dict["absorption"]["LUTfiles"])
            @assert size(files_lut) == size(molecules) "Size of LUTfiles has to match molecules"
            @show size(files_lut)
            
            for i in eachindex(files_lut)
                #@show i, files_lut[i]
                #@show typeof(load_interpolation_model(files_lut[1]))
                push!(luts,[load_interpolation_model(file) for file in files_lut[i]])
            end
        end
        absorption_params = AbsorptionParameters(molecules, vmr, broadening_function, 
                                                 CEF, wing_cutoff, luts)



    else 
        absorption_params = nothing
    end
    
    # scattering group
    if "scattering" in keys(params_dict)
        #Force aerosols to be Float64!!
        FTa = Float64
        # Get aerosols from types
        aerosols = aerosol_params_to_obj(params_dict["scattering"]["aerosols"], FTa)
        r_max = FTa(params_dict["scattering"]["r_max"])
        nquad_radius = params_dict["scattering"]["nquad_radius"]
        λ_ref = FTa(params_dict["scattering"]["λ_ref"])
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
                                convert.(FT,T), convert.(FT,p), convert.(FT,q), 
                                params_dict["atmospheric_profile"]["profile_reduction"],

                                # absorption group
                                absorption_params,

                                # scattering group
                                scattering_params
                                );
end