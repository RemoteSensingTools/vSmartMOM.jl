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

"""
FLOAT_MAP

Maps configuration strings to floating-point types.
Supported: "Float32", "Float64".
"""
const FLOAT_MAP = Dict(
    "Float32" => Float32,
    "Float64" => Float64,
)

"""
POLARIZATION_MAP

String → constructor mapping for polarization types parameterized by FT.
Keys: "Stokes_I", "Stokes_IQ", "Stokes_IQU", "Stokes_IQUV".
"""
const POLARIZATION_MAP = Dict(
    "Stokes_I"   => (FT)->Scattering.Stokes_I{FT}(),
    "Stokes_IQ"  => (FT)->Scattering.Stokes_IQ{FT}(),
    "Stokes_IQU" => (FT)->Scattering.Stokes_IQU{FT}(),
    "Stokes_IQUV"=> (FT)->Scattering.Stokes_IQUV{FT}()
)

"""
ARCH_MAP

String → constructor/value mapping for compute architecture.
Keys: "CPU", "GPU", "default_architecture".
"""
const ARCH_MAP = Dict(
    "CPU" => ()->Architectures.CPU(),
    "GPU" => ()->Architectures.GPU(),
    # default_architecture is a value (GPU() or CPU()), not a function
    "default_architecture" => ()->Architectures.default_architecture()
)

"""
QUAD_MAP

String → quadrature type mapping for stream calculations.
Keys: "RadauQuad", "GaussQuadHemisphere", "GaussQuadFullSphere".
"""
const QUAD_MAP = Dict(
    "RadauQuad"           => ()->CoreRT.RadauQuad(),
    "GaussQuadHemisphere" => ()->CoreRT.GaussQuadHemisphere(),
    "GaussQuadFullSphere" => ()->CoreRT.GaussQuadFullSphere(),
)

"""
DECOMP_MAP

String → Fourier decomposition method mapping for scattering.
Keys: "NAI2", "PCW".
"""
const DECOMP_MAP = Dict(
    "NAI2" => ()->Scattering.NAI2(),
    "PCW"  => ()->Scattering.PCW()
)

"""
BROADENING_MAP

String → absorption broadening function mapping.
Keys: "Voigt", "Lorentz", "Doppler".
"""
const BROADENING_MAP = Dict(
    "Voigt"   => ()->Absorption.Voigt(),
    "Lorentz" => ()->Absorption.Lorentz(),
    "Doppler" => ()->Absorption.Doppler()
)

"""
CEF_MAP

String → complex error function mapping for Voigt computations.
Keys: "HumlicekWeidemann32SDErrorFunction".
"""
const CEF_MAP = Dict(
    "HumlicekWeidemann32SDErrorFunction" => ()->Absorption.HumlicekWeidemann32SDErrorFunction()
)

# --- BRDF (surface) mapping ---
"""
BRDF_MAP

String → factory mapping for surface BRDF constructors.
Supported:
- LambertianSurfaceScalar(albedo)
- LambertianSurfaceSpectrum([albedo_per_band])
- LambertianSurfaceLegendre([coefficients])
- rpvSurfaceScalar(ρ₀, ρ_c, k, Θ)
- RossLiSurfaceScalar(fvol, fgeo, fiso)
"""
const BRDF_MAP = Dict{String, Function}(
    "LambertianSurfaceScalar" => (FT, args)-> begin
        @assert length(args) == 1 "LambertianSurfaceScalar expects 1 argument (albedo)."
        CoreRT.LambertianSurfaceScalar(FT(args[1]))
    end,
    "LambertianSurfaceSpectrum" => (FT, args)-> begin
        @assert length(args) == 1 && args[1] isa AbstractVector "LambertianSurfaceSpectrum expects 1 vector argument."
        CoreRT.LambertianSurfaceSpectrum(convert.(FT, collect(args[1])))
    end,
    "LambertianSurfaceLegendre" => (FT, args)-> begin
        @assert length(args) == 1 && args[1] isa AbstractVector "LambertianSurfaceLegendre expects 1 vector argument."
        CoreRT.LambertianSurfaceLegendre(convert.(FT, collect(args[1])))
    end,
    "rpvSurfaceScalar" => (FT, args)-> begin
        @assert length(args) == 4 "rpvSurfaceScalar expects 4 arguments (ρ₀, ρ_c, k, Θ)."
        CoreRT.rpvSurfaceScalar(FT(args[1]), FT(args[2]), FT(args[3]), FT(args[4]))
    end,
    "RossLiSurfaceScalar" => (FT, args)-> begin
        @assert length(args) == 3 "RossLiSurfaceScalar expects 3 arguments (fvol, fgeo, fiso)."
        CoreRT.RossLiSurfaceScalar(FT(args[1]), FT(args[2]), FT(args[3]))
    end,
)

"Split a comma-separated argument string into top-level arguments while respecting brackets"
function _split_args(args_str::AbstractString)
    args = String[]
    buf = IOBuffer()
    depth = 0
    for c in collect(args_str)
        if c == '['
            depth += 1
            print(buf, c)
        elseif c == ']'
            depth -= 1
            print(buf, c)
        elseif c == ',' && depth == 0
            push!(args, String(take!(buf)) |> x->strip(x))
        else
            print(buf, c)
        end
    end
    s = String(take!(buf)) |> x->strip(x)
    if !isempty(s)
        push!(args, s)
    end
    return args
end

"Parse a single argument into a number or a vector of numbers (no eval)"
function _parse_arg(arg::AbstractString, FT)
    a = strip(arg)
    if isempty(a)
        return nothing
    end
    if startswith(a, "[") && endswith(a, "]")
        inner = strip(a[2:end-1])
        isempty(inner) && return FT[]
        # split by comma or whitespace
        tokens = split(inner, r"[,\s]+"; keepempty=false)
        return [FT(parse(Float64, t)) for t in tokens]
    else
        return FT(parse(Float64, a))
    end
end

"""
Parse a surface spec like "LambertianSurfaceScalar(0.1)" into a CoreRT surface instance, without eval.
"""
function parse_surface_str(s::AbstractString, FT)
    # Match "Name(args)" or "Name{Type}(args)" -- strip optional type parameter
    name_args = match(r"^(\w+)(?:\{[^}]*\})?\((.*)\)$", strip(String(s)))
    @assert name_args !== nothing "Invalid surface specification: '$(s)'"
    name = name_args.captures[1]
    args_raw = name_args.captures[2]
    args_list = isempty(strip(args_raw)) ? Any[] : [_parse_arg(x, FT) for x in _split_args(args_raw)]
    @assert haskey(BRDF_MAP, name) "Unknown surface type $(name)."
    return BRDF_MAP[name](FT, args_list)
end

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
function _parse_float_type(params_dict::Dict)
    key = String(params_dict["radiative_transfer"]["float_type"])
    @assert haskey(FLOAT_MAP, key) "Unknown float_type $(key)."
    return FLOAT_MAP[key]
end

function _parse_spec_bands(params_dict::Dict, FT)
    spec_bands = Vector{Vector{FT}}()
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
        push!(spec_bands, convert.(FT, vec(wn_band)))
    end
    return spec_bands
end

function _parse_surfaces(params_dict::Dict, FT)
    return [parse_surface_str(String(s), FT) for s in params_dict["radiative_transfer"]["surface"]]
end

function _parse_quadrature(params_dict::Dict)
    qkey = replace(String(params_dict["radiative_transfer"]["quadrature_type"]),"()"=>"")
    @assert haskey(QUAD_MAP, qkey) "Unknown quadrature_type $(qkey)."
    return QUAD_MAP[qkey]()
end

function _parse_polarization(params_dict::Dict, FT)
    pol_key = replace(String(params_dict["radiative_transfer"]["polarization_type"]),"()"=>"")
    @assert haskey(POLARIZATION_MAP, pol_key) "Unknown polarization_type $(pol_key)."
    return POLARIZATION_MAP[pol_key](FT)
end

function _parse_architecture(params_dict::Dict)
    arch_key = String(params_dict["radiative_transfer"]["architecture"]) |> x->replace(x, "Architectures."=>"") |> x->replace(x, "()"=>"")
    @assert haskey(ARCH_MAP, arch_key) "Unknown architecture $(arch_key)."
    return ARCH_MAP[arch_key]()
end

function _parse_atmosphere(params_dict::Dict, FT)
    T = convert.(FT, params_dict["atmospheric_profile"]["T"]) # Level
    p = convert.(FT, params_dict["atmospheric_profile"]["p"]) # Boundaries
    q = haskey(params_dict["atmospheric_profile"], "q") ? convert.(FT, params_dict["atmospheric_profile"]["q"]) : zeros(FT, length(T))
    prof_red = params_dict["atmospheric_profile"]["profile_reduction"]
    return T, p, q, prof_red
end

function _parse_absorption(params_dict::Dict, FT)
    if !haskey(params_dict, "absorption")
        return nothing
    end
    molecules = Array(params_dict["absorption"]["molecules"])
    vmr = convert(Dict{String, Union{Real, Vector}}, params_dict["absorption"]["vmr"])
    validate_vmrs(params_dict["absorption"]["molecules"], vmr)
    bkey = replace(String(params_dict["absorption"]["broadening"]),"()"=>"")
    @assert haskey(BROADENING_MAP, bkey) "Unknown broadening $(bkey)."
    broadening_function = BROADENING_MAP[bkey]()
    ckey = replace(String(params_dict["absorption"]["CEF"]),"()"=>"")
    @assert haskey(CEF_MAP, ckey) "Unknown CEF $(ckey)."
    CEF = CEF_MAP[ckey]()
    wing_cutoff = FT(params_dict["absorption"]["wing_cutoff"])
    luts = []
    if haskey(params_dict["absorption"], "LUTfiles")
        files_lut = Array(params_dict["absorption"]["LUTfiles"])
        @assert size(files_lut) == size(molecules) "Size of LUTfiles has to match molecules"
        for i in eachindex(files_lut)
            push!(luts, [load_interpolation_model(file) for file in files_lut[i]])
        end
    end
    # Partition molecules into fixed (no Jacobian) and variable (Jacobian computed).
    # Default: all molecules are fixed unless explicitly listed as variable in YAML.
    if haskey(params_dict["absorption"], "variable_molecules")
        variable_molecules = Array(params_dict["absorption"]["variable_molecules"])
        # Fixed = molecules minus variable per band
        fixed_molecules = [setdiff(molecules[i], 
            i <= length(variable_molecules) ? variable_molecules[i] : String[]) 
            for i in 1:length(molecules)]
        # Pad variable_molecules to match n_bands if shorter
        while length(variable_molecules) < length(molecules)
            push!(variable_molecules, String[])
        end
    else
        # No linearization requested: all molecules are fixed, none variable
        fixed_molecules = deepcopy(molecules)
        variable_molecules = [String[] for _ in 1:length(molecules)]
    end
    return AbsorptionParameters(molecules, fixed_molecules, variable_molecules, vmr, broadening_function, CEF, wing_cutoff, luts)
end

function _parse_scattering(params_dict::Dict, FT::Type{<:AbstractFloat}=Float64)
    if !haskey(params_dict, "scattering")
        return nothing
    end
    FTa = FT
    aerosols = aerosol_params_to_obj(params_dict["scattering"]["aerosols"], FTa)
    r_max = FTa(params_dict["scattering"]["r_max"])
    nquad_radius = params_dict["scattering"]["nquad_radius"]
    λ_ref = FTa(params_dict["scattering"]["λ_ref"])
    dkey = replace(String(params_dict["scattering"]["decomp_type"]),"()"=>"")
    @assert haskey(DECOMP_MAP, dkey) "Unknown decomp_type $(dkey)."
    decomp_type = DECOMP_MAP[dkey]()
    n_ref = if !haskey(params_dict["scattering"],"n_ref")
        aerosols[1].aerosol.nᵣ - im*aerosols[1].aerosol.nᵢ
    else
        parse(Complex{FTa}, params_dict["scattering"]["n_ref"])
    end
    return ScatteringParameters(aerosols, r_max, nquad_radius, λ_ref, n_ref, decomp_type)
end

"""
parameters_from_dict(params_dict::Dict) -> vSmartMOM_Parameters

Convert a configuration Dict (e.g., loaded from YAML) into vSmartMOM_Parameters.
Validates schema and parses radiative_transfer, geometry, atmospheric_profile,
and optional absorption/scattering sections.
"""
function parameters_from_dict(params_dict::Dict)
    validate_yaml_parameters(params_dict)
    FT = _parse_float_type(params_dict)
    spec_bands = _parse_spec_bands(params_dict, FT)
    BRDF_per_band = _parse_surfaces(params_dict, FT)
    quadrature_type = _parse_quadrature(params_dict)
    polarization_type = _parse_polarization(params_dict, FT)
    architecture = _parse_architecture(params_dict)
    T, p, q, profile_reduction = _parse_atmosphere(params_dict, FT)
    absorption_params = _parse_absorption(params_dict, FT)
    scattering_params = _parse_scattering(params_dict, FT)

    return vSmartMOM_Parameters(
        spec_bands, BRDF_per_band, quadrature_type, polarization_type,
        params_dict["radiative_transfer"]["max_m"], FT(params_dict["radiative_transfer"]["Δ_angle"]),
        params_dict["radiative_transfer"]["l_trunc"], FT(params_dict["radiative_transfer"]["depol"]), FT,
        architecture,
        FT(params_dict["geometry"]["sza"]), convert.(FT, params_dict["geometry"]["vza"]),
        convert.(FT, params_dict["geometry"]["vaz"]), FT(params_dict["geometry"]["obs_alt"]),
        T, p, q, profile_reduction,
        absorption_params, scattering_params
    )
end

# Convenience wrappers
"""
read_parameters(x)

Multi-dispatch convenience wrapper:
- String (file path) → parameters_from_yaml
- Dict → parameters_from_dict
- vSmartMOM_Parameters → pass-through
"""
read_parameters(x::String) = parameters_from_yaml(x)
read_parameters(x::Dict)   = parameters_from_dict(x)
read_parameters(p::vSmartMOM_Parameters) = p
