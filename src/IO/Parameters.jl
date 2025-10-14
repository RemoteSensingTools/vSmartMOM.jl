# Parameter IO: YAML and Dict converters

# Bring needed names into scope
using YAML
using Unitful
using UnitfulEquivalences
using Distributions
using Parameters: @unpack
using ..CoreRT: vSmartMOM_Parameters, AbsorptionParameters, ScatteringParameters, RT_Aerosol, AtmosphericProfile
using ..Absorption: AbstractBroadeningFunction, AbstractComplexErrorFunction, load_interpolation_model
using ..Scattering: AbstractFourierDecompositionType, Aerosol
using ..Architectures

"Check that a field exists in yaml file"
function check_yaml_field(dict::AbstractDict, full_keys::Vector{String}, curr_keys::Vector{String}, final_type, valid_options::Vector{String})
    @assert length(curr_keys) >= 1
    if length(curr_keys) == 1
        @assert curr_keys[1] in keys(dict) "Missing key in parameters yaml: $(join(full_keys, '/'))"
        @assert dict[curr_keys[1]] isa final_type "Improper type for $(join(full_keys, '/')); must be a $(final_type)" 
        if !isempty(valid_options)
            @assert dict[curr_keys[1]] in valid_options "Field $(join(full_keys, '/')) must be one of $(valid_options)"
        end
        return true
    else
        @assert curr_keys[1] in keys(dict) "Missing key in parameters yaml: $(join(full_keys, '/'))"
        return check_yaml_field(dict[curr_keys[1]], full_keys, curr_keys[2:end], final_type, valid_options)
    end
end

"Given a list of aerosol parameter dictionaries, validate all of them"
function validate_aerosols(aerosols)
    @assert !isempty(aerosols) "Aerosols list shouldn't be empty if scattering block is included"
    fields = [(["τ_ref"], Real), (["μ"], Real), (["σ"], Real), (["nᵣ"], Real), (["nᵢ"], Real), (["p₀"], Real), (["p₀"], Real)]
    for aerosol in aerosols
        for f in fields
            key, elty = f[1:2]
            @assert check_yaml_field(aerosol, key, key, elty, String[]) 
        end
    end
end

"Convert the input dictionary of aerosols into a list of RT_aerosols"
function aerosol_params_to_obj(aerosols, FT)
    rt_aerosol_obj_list = RT_Aerosol[]
    for aerosol in aerosols
        @assert aerosol["σ"] ≥ 1 "Geometric standard deviation has to be ≥ 1"    
        size_distribution = LogNormal(log(FT(aerosol["μ"])), log(FT(aerosol["σ"])))
        new_aerosol_obj = Aerosol(size_distribution, FT(aerosol["nᵣ"]), FT(aerosol["nᵢ"]))
        new_rt_aerosol_obj = RT_Aerosol(new_aerosol_obj, FT(aerosol["τ_ref"]), Normal(FT(aerosol["p₀"]), FT(aerosol["σp"])))
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
    fields = [
        (["radiative_transfer", "spec_bands"], Array{String}),  
        (["radiative_transfer", "surface"], Array{String}),   
        (["radiative_transfer", "quadrature_type"], String), 
        (["radiative_transfer", "polarization_type"], String),   
        (["radiative_transfer", "max_m"], Integer), 
        (["radiative_transfer", "Δ_angle"], Real),
        (["radiative_transfer", "l_trunc"], Integer),
        (["radiative_transfer", "depol"], Real),
        (["radiative_transfer", "float_type"], String),
        (["radiative_transfer", "architecture"], String, ["default_architecture", "Architectures.GPU()", "Architectures.CPU()", "GPU()", "CPU()"]),
        (["geometry", "sza"], Real),
        (["geometry", "vza"], Array{<:Real}),
        (["geometry", "vaz"], Array{<:Real}),
        (["geometry", "obs_alt"], Real),
        (["atmospheric_profile", "T"], Array{<:Real}),
        (["atmospheric_profile", "p"], Array{<:Real}),
        (["atmospheric_profile", "profile_reduction"], Union{Integer, Nothing}),
        (["absorption", "molecules"], Array),
        (["absorption", "vmr"], Dict),
        (["absorption", "broadening"], String, ["Voigt()", "Doppler()", "Lorentz()"]), 
        (["absorption", "CEF"], String, ["HumlicekWeidemann32SDErrorFunction()"]),
        (["absorption", "wing_cutoff"], Integer),
        (["scattering", "aerosols"], Array),
        (["scattering", "r_max"], Real),
        (["scattering", "nquad_radius"], Real),
        (["scattering", "λ_ref"], Real),
        (["scattering", "decomp_type"], String, ["NAI2()", "PCW()"]),
    ]
    for f in fields
        key, elty = f[1:2]
        valid_options = length(f) == 3 ? f[3] : String[]
        top = key[1]
        if top == "absorption"
            if haskey(params, "absorption")
                @assert check_yaml_field(params, key, key, elty, valid_options)
            end
        elseif top == "scattering"
            if haskey(params, "scattering")
                @assert check_yaml_field(params, key, key, elty, valid_options)
            end
        else
            @assert check_yaml_field(params, key, key, elty, valid_options)
        end
    end
    if "scattering" in keys(params)
        validate_aerosols(params["scattering"]["aerosols"])
    end
end

"Build parameters from a Dict (e.g., parsed YAML)"
function parameters_from_dict(params_dict::Dict)
    validate_yaml_parameters(params_dict)
    FT = eval(Meta.parse(params_dict["radiative_transfer"]["float_type"]))
    spec_bands = []
    for spec_band in params_dict["radiative_transfer"]["spec_bands"]
        parsed_band = eval(Meta.parse(spec_band))
        if (all(x-> x == NoUnits, unit.(parsed_band)))
            wn_band = collect(parsed_band)
        else
            if (all(x-> x == u"cm^-1", unit.(parsed_band)))
                wn_band = sort(ustrip(collect(parsed_band)))
            else
                wn_band = sort(ustrip(uconvert.( u"cm^-1", collect(parsed_band), Spectral())))
            end
        end
        final_band = convert(Array{FT}, wn_band)
        push!(spec_bands, final_band)
    end
    BRDF_per_band   = map(x -> eval(Meta.parse(x)), params_dict["radiative_transfer"]["surface"]) 
    quadrature_type = eval(Meta.parse(params_dict["radiative_transfer"]["quadrature_type"]))
    pol_type_str    = replace(params_dict["radiative_transfer"]["polarization_type"], "()" => "{$FT}()")
    polarization_type = eval(Meta.parse(pol_type_str))
    arch_string = params_dict["radiative_transfer"]["architecture"]
    arch_string = !startswith(arch_string, "Architectures.") ? "Architectures." * arch_string : arch_string
    architecture = eval(Meta.parse(arch_string))

    T = convert.(FT, params_dict["atmospheric_profile"]["T"]) # Level
    p = convert.(FT, params_dict["atmospheric_profile"]["p"]) # Boundaries
    q = "q" in keys(params_dict["atmospheric_profile"]) ? convert.(FT, params_dict["atmospheric_profile"]["q"]) : zeros(length(T))

    if "absorption" in keys(params_dict)
        molecules = Array(params_dict["absorption"]["molecules"])
        vmr = convert(Dict{String, Union{Real, Vector}}, params_dict["absorption"]["vmr"])
        validate_vmrs(params_dict["absorption"]["molecules"], vmr)
        broadening_function = eval(Meta.parse(params_dict["absorption"]["broadening"]))
        CEF = eval(Meta.parse(params_dict["absorption"]["CEF"]))
        wing_cutoff = FT(params_dict["absorption"]["wing_cutoff"])
        luts = []
        if "LUTfiles" in keys(params_dict["absorption"])
            files_lut = Array(params_dict["absorption"]["LUTfiles"])
            @assert size(files_lut) == size(molecules) "Size of LUTfiles has to match molecules"
            for i in eachindex(files_lut)
                push!(luts,[load_interpolation_model(file) for file in files_lut[i]])
            end
        end
        absorption_params = AbsorptionParameters(molecules, vmr, broadening_function, CEF, wing_cutoff, luts)
    else 
        absorption_params = nothing
    end

    if "scattering" in keys(params_dict)
        FTa = Float64
        aerosols = aerosol_params_to_obj(params_dict["scattering"]["aerosols"], FTa)
        r_max = FTa(params_dict["scattering"]["r_max"])
        nquad_radius = params_dict["scattering"]["nquad_radius"]
        λ_ref = FTa(params_dict["scattering"]["λ_ref"])
        decomp_type = eval(Meta.parse(params_dict["scattering"]["decomp_type"]))
        n_ref = if !haskey(params_dict["scattering"],"n_ref")
            aerosols[1].aerosol.nᵣ - im*aerosols[1].aerosol.nᵢ
        else
            parse(Complex{FTa},params_dict["scattering"]["n_ref"])
        end
        scattering_params = ScatteringParameters(aerosols, r_max, nquad_radius, λ_ref, n_ref, decomp_type)
    else
        scattering_params = nothing
    end

    return vSmartMOM_Parameters(
        spec_bands, BRDF_per_band, quadrature_type, polarization_type,
        params_dict["radiative_transfer"]["max_m"], FT(params_dict["radiative_transfer"]["Δ_angle"]),
        params_dict["radiative_transfer"]["l_trunc"], FT(params_dict["radiative_transfer"]["depol"]), FT,
        architecture,
        FT(params_dict["geometry"]["sza"]), convert.(FT, params_dict["geometry"]["vza"]),
        convert.(FT, params_dict["geometry"]["vaz"]), FT(params_dict["geometry"]["obs_alt"]),
        convert.(FT,T), convert.(FT,p), convert.(FT,q), params_dict["atmospheric_profile"]["profile_reduction"],
        absorption_params, scattering_params
    )
end

"Load parameters from YAML file path"
function parameters_from_yaml(file_path::AbstractString)
    params_dict = YAML.load_file(file_path)
    return parameters_from_dict(params_dict)
end

# Convenience wrappers
read_parameters(x::String) = parameters_from_yaml(x)
read_parameters(x::Dict)   = parameters_from_dict(x)
read_parameters(p::vSmartMOM_Parameters) = p
