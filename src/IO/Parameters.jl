# Parameter IO: YAML and Dict converters

# Bring needed names into scope
using YAML
using Distributions
using Unitful
using UnitfulEquivalences
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
- CoxMunkSurface(wind_speed=U) or CoxMunkSurface(U)
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
    "CoxMunkSurface" => (FT, args)-> begin
        # Called with kwargs dict as first positional arg (from parse_surface_str)
        # or with a single positional wind_speed
        if length(args) == 1 && args[1] isa Dict
            kw = args[1]
            ws = FT(kw[:wind_speed])
            n_raw = get(kw, :n_water, nothing)
            n_water = n_raw === nothing ? nothing : Complex{FT}(n_raw)
            wc_alb = FT(get(kw, :whitecap_albedo, FT(0.22)))
            inc_wc = Bool(get(kw, :include_whitecaps, true))
            shad   = Bool(get(kw, :shadowing, true))
            CoreRT.CoxMunkSurface(wind_speed=ws, n_water=n_water,
                                  whitecap_albedo=wc_alb,
                                  include_whitecaps=inc_wc,
                                  shadowing=shad)
        else
            @assert length(args) == 1 "CoxMunkSurface expects 1 argument (wind_speed) or keyword arguments."
            CoreRT.CoxMunkSurface(wind_speed=FT(args[1]))
        end
    end,
)

"""
    _parse_canopy_section(canopy_dict, FT, brdf_list)

Parse the optional `canopy:` YAML section and wrap each band's existing BRDF
in a `CanopySurface`, replacing the `brdf_list` entries in-place.

Expected YAML keys:
- `LAI`: total leaf area index (required)
- `n_layers`: number of canopy sub-layers (default 1)
- `leaf_reflectance`: scalar or list (default 0.4)
- `leaf_transmittance`: scalar or list (default 0.05)
- `leaf_optics_grid`: wavelength/wavenumber grid for spectral leaf R/T (optional)
- `grid_unit`: "nm" or "cm_inv" (default "nm")
- `include_atm`: bool (default false)
- `canopy_dp`: pressure thickness of canopy air column in hPa (optional)
- `soil`: a BRDF string *or* "from_surface" to reuse the band's existing BRDF
"""
function _parse_canopy_section(canopy_dict::Dict, FT, brdf_list::Vector)
    LAI       = FT(get(canopy_dict, "LAI", 3.0))
    n_layers  = get(canopy_dict, "n_layers", 1)
    leaf_R    = get(canopy_dict, "leaf_reflectance", 0.4)
    leaf_T    = get(canopy_dict, "leaf_transmittance", 0.05)
    incl_atm  = get(canopy_dict, "include_atm", false)
    dp_raw    = get(canopy_dict, "canopy_dp", nothing)
    dp        = dp_raw === nothing ? nothing : FT(dp_raw)

    lr = leaf_R isa AbstractVector ? convert(Vector{FT}, leaf_R) : FT(leaf_R)
    lt = leaf_T isa AbstractVector ? convert(Vector{FT}, leaf_T) : FT(leaf_T)

    lg_raw = get(canopy_dict, "leaf_optics_grid", nothing)
    lg = lg_raw === nothing ? nothing : convert(Vector{FT}, lg_raw)
    gu = Symbol(get(canopy_dict, "grid_unit", "nm"))

    soil_spec = get(canopy_dict, "soil", "from_surface")

    for (i, existing_brdf) in enumerate(brdf_list)
        if soil_spec == "from_surface"
            soil = existing_brdf
        else
            soil = _parse_brdf_string(soil_spec, FT)
        end
        brdf_list[i] = CoreRT.CanopySurface(;
            soil=soil, LAI=LAI, n_layers=n_layers,
            leaf_reflectance=lr, leaf_transmittance=lt,
            leaf_optics_grid=lg, grid_unit=gu,
            include_atm=incl_atm, canopy_dp=dp)
    end
    return brdf_list
end

_parse_brdf_string(s::AbstractString, FT) = parse_surface_str(s, FT)

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

"Parse a single argument into a number, boolean, or a vector of numbers (no eval)"
function _parse_arg(arg::AbstractString, FT)
    a = strip(arg)
    if isempty(a)
        return nothing
    end
    # Boolean literals
    if lowercase(a) == "true"
        return true
    elseif lowercase(a) == "false"
        return false
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
Parse a surface spec like "LambertianSurfaceScalar(0.1)" or
"CoxMunkSurface(wind_speed=5.0)" into a CoreRT surface instance, without eval.

Supports both positional arguments and keyword arguments (key=value syntax).
When keyword arguments are detected, a `Dict{Symbol,Any}` is passed as the
sole element of the args list to the BRDF_MAP constructor.
"""
function parse_surface_str(s::AbstractString, FT)
    # Match "Name(args)" or "Name{Type}(args)" -- strip optional type parameter
    name_args = match(r"^(\w+)(?:\{[^}]*\})?\((.*)\)$", strip(String(s)))
    @assert name_args !== nothing "Invalid surface specification: '$(s)'"
    name = name_args.captures[1]
    args_raw = name_args.captures[2]
    @assert haskey(BRDF_MAP, name) "Unknown surface type $(name)."

    if isempty(strip(args_raw))
        return BRDF_MAP[name](FT, Any[])
    end

    raw_parts = _split_args(args_raw)

    # Detect keyword arguments (key=value, excluding vector args like [1=...])
    has_kwargs = any(p -> occursin('=', p) && !startswith(strip(p), "["), raw_parts)

    if has_kwargs
        # Pack keyword arguments into a Dict passed as args[1]
        kwargs = Dict{Symbol,Any}()
        for part in raw_parts
            p = strip(part)
            if occursin('=', p) && !startswith(p, "[")
                kv = split(p, '='; limit=2)
                kwargs[Symbol(strip(kv[1]))] = _parse_arg(strip(kv[2]), FT)
            end
        end
        return BRDF_MAP[name](FT, Any[kwargs])
    else
        args_list = [_parse_arg(x, FT) for x in raw_parts]
        return BRDF_MAP[name](FT, args_list)
    end
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

"Given a list of aerosol parameter dictionaries, validate all of them.

Vertical distribution may be specified either as (z₀, σ₀) — altitude-form log-normal,
preferred going forward to match sanghavi's YAML convention — or as (p₀, σp) —
pressure-form normal, the legacy unified convention. Exactly one form must be
present per aerosol."
function validate_aerosols(aerosols)
    @assert !isempty(aerosols) "Aerosols list shouldn't be empty if scattering block is included"
    base_fields = [(["τ_ref"], Real), (["μ"], Real), (["σ"], Real),
                   (["nᵣ"], Real), (["nᵢ"], Real)]
    for aerosol in aerosols
        for f in base_fields
            key, elty = f[1:2]
            @assert check_yaml_field(aerosol, key, key, elty, String[])
        end
        has_alt  = haskey(aerosol, "z₀") && haskey(aerosol, "σ₀")
        has_pres = haskey(aerosol, "p₀") && haskey(aerosol, "σp")
        @assert has_alt || has_pres "Aerosol must specify vertical distribution as either (z₀, σ₀) [altitude-form, preferred] or (p₀, σp) [pressure-form, legacy]"
        @assert !(has_alt && has_pres) "Aerosol must specify exactly one of (z₀, σ₀) or (p₀, σp) — got both"
    end
end

"""Convert the input dictionary of aerosols into a list of RT_aerosols.

When (z₀, σ₀) is provided (altitude-form, preferred), the profile is stored as
`LogNormal(log(z₀), σ₀)` and consumed by the altitude-space `getAerosolLayerOptProp`
signature. When (p₀, σp) is provided (pressure-form, legacy), the profile is
stored as `Normal(p₀, σp)` and consumed by the pressure-space signature.

NOTE (Phase 1b): the altitude-form → pressure-grid integration path in
`getAerosolLayerOptProp(total_τ, dist::Distribution, profile::AtmosphericProfile)`
still interprets `dist` in pressure space. When τ_ref = 0 (e.g. the Phase 1b
regression gate), this is a no-op. Proper altitude-form dispatch will land in a
follow-up alongside the aerosol-module wire-in (Phase 1d) or the workspace
landing (Phase 4), whichever proves more natural."""
function aerosol_params_to_obj(aerosols, FT)
    rt_aerosol_obj_list = RT_Aerosol{FT}[]
    for aerosol in aerosols
        @assert aerosol["σ"] ≥ 1 "Geometric standard deviation has to be ≥ 1"
        size_distribution = LogNormal(log(FT(aerosol["μ"])), log(FT(aerosol["σ"])))
        new_aerosol_obj = Aerosol(size_distribution, FT(aerosol["nᵣ"]), FT(aerosol["nᵢ"]))
        profile = if haskey(aerosol, "z₀")
            LogNormal(log(FT(aerosol["z₀"])), FT(aerosol["σ₀"]))
        else
            Normal(FT(aerosol["p₀"]), FT(aerosol["σp"]))
        end
        new_rt_aerosol_obj = RT_Aerosol(new_aerosol_obj, FT(aerosol["τ_ref"]), profile)
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

"""
    _safe_parse_number(s) -> Float64

Parse a numeric expression that may contain division, e.g. "1e7/775".
Supports integers, floats, scientific notation, and a single `/` operator.
"""
function _safe_parse_number(s::AbstractString)
    s = strip(s)
    if occursin('/', s)
        parts = split(s, '/'; limit=2)
        return parse(Float64, strip(parts[1])) / parse(Float64, strip(parts[2]))
    else
        return parse(Float64, s)
    end
end

"""
    _safe_parse_spec_band(s) -> Vector{Float64}

Parse a spectral band string without `eval`. Supported formats:
- Range: `"(1e7/775):0.05:(1e7/755)"` or `"start:step:stop"`
- Explicit list: `"[19417.0 19418.0]"` or `"[19417.0, 19418.0]"`
"""
function _safe_parse_spec_band(s::AbstractString)
    s = strip(s)
    # Strip outer quotes if present
    if startswith(s, '"') && endswith(s, '"')
        s = s[2:end-1]
    end

    # Handle collect(...)u"nm" or collect(...)u"cm_inv" wrapper
    unit = :cm_inv
    m = match(r"^collect\((.+)\)\s*u\"(nm|cm_inv)\"$", s)
    if m !== nothing
        s = m.captures[1]
        unit = Symbol(m.captures[2])
    elseif startswith(s, "collect(") && endswith(s, ")")
        s = s[9:end-1]  # strip collect(...)
    end

    # Bracket list: [val1 val2 ...] or [val1, val2, ...]
    if startswith(s, '[') && endswith(s, ']')
        inner = strip(s[2:end-1])
        tokens = split(inner, r"[,\s]+"; keepempty=false)
        vals = [parse(Float64, strip(t)) for t in tokens]
        if unit == :nm
            vals = sort(1e7 ./ vals)
        end
        return vals
    end

    # Range: expr:step:expr or (expr):step:(expr)
    # Split on ':' that is NOT inside parentheses
    range_parts = String[]
    buf = IOBuffer()
    depth = 0
    for c in s
        if c == '('
            depth += 1; print(buf, c)
        elseif c == ')'
            depth -= 1; print(buf, c)
        elseif c == ':' && depth == 0
            push!(range_parts, String(take!(buf)))
        else
            print(buf, c)
        end
    end
    push!(range_parts, String(take!(buf)))

    local vals::Vector{Float64}
    if length(range_parts) == 3
        start = _safe_parse_number(strip(range_parts[1], ['(', ')']))
        step  = _safe_parse_number(strip(range_parts[2], ['(', ')']))
        stop  = _safe_parse_number(strip(range_parts[3], ['(', ')']))
        vals = collect(start:step:stop)
    elseif length(range_parts) == 2
        start = _safe_parse_number(strip(range_parts[1], ['(', ')']))
        stop  = _safe_parse_number(strip(range_parts[2], ['(', ')']))
        vals = collect(start:stop)
    else
        vals = [_safe_parse_number(s)]
    end

    if unit == :nm
        vals = sort(1e7 ./ vals)
    end
    return vals
end

function _parse_spec_bands(params_dict::Dict, FT)
    spec_bands = Vector{Vector{FT}}()
    for spec_band in params_dict["radiative_transfer"]["spec_bands"]
        wn_band = _safe_parse_spec_band(spec_band)
        push!(spec_bands, convert.(FT, vec(wn_band)))
    end
    return spec_bands
end

function _parse_surfaces(params_dict::Dict, FT)
    surfs = CoreRT.AbstractSurfaceType[parse_surface_str(String(s), FT) for s in params_dict["radiative_transfer"]["surface"]]
    return surfs
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

    if haskey(params_dict, "canopy")
        _parse_canopy_section(params_dict["canopy"], FT, BRDF_per_band)
    end

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
