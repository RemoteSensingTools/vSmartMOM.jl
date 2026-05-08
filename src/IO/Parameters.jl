# Parameter IO: YAML and Dict converters

# Bring needed names into scope
using YAML
using Distributions
using Unitful
using UnitfulEquivalences
using CanopyOptics
using ..CoreRT: vSmartMOM_Parameters, RTNumericalParameters, AbsorptionParameters, ScatteringParameters, RT_Aerosol, AtmosphericProfile
using ..Absorption: AbstractBroadeningFunction, AbstractComplexErrorFunction, load_interpolation_model
using ..Scattering
using ..Scattering: AbstractFourierDecompositionType, Aerosol
using ..Architectures

"""
    _config_error(msg)

Raise a stable `ArgumentError` for invalid user configuration input.
"""
@inline _config_error(msg) = throw(ArgumentError(msg))

"""
    _require_config(cond, msg)

Validate a user configuration condition and raise `ArgumentError` when it
fails. This keeps parameter parsing independent of Julia's assert settings.
"""
@inline _require_config(cond, msg) = cond ? nothing : _config_error(msg)

const ENV_PATH_PATTERN = r"^\$\{ENV:([A-Za-z_][A-Za-z0-9_]*)\}(.*)$"

_strip_leading_path_separators(path::AbstractString) = replace(String(path), r"^[\\/]+" => "")

"""
    _expand_env_path(path)

Resolve an optional leading `\${ENV:NAME}` marker in a configuration path.
This keeps checked-in parameter files portable while still allowing large
local assets, such as absorption LUTs, to live outside the repository.
"""
function _expand_env_path(path::AbstractString)
    text = String(path)
    match_env = match(ENV_PATH_PATTERN, text)
    match_env === nothing && return text

    env_name, suffix = match_env.captures
    _require_config(haskey(ENV, env_name), "Environment variable $(env_name) must be set to resolve $(text)")
    _require_config(!isempty(ENV[env_name]), "Environment variable $(env_name) must not be empty to resolve $(text)")
    isempty(suffix) && return normpath(ENV[env_name])
    return normpath(joinpath(ENV[env_name], _strip_leading_path_separators(suffix)))
end

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
Keys: "CPU", "GPU", "MetalGPU", "default_architecture".
"""
const ARCH_MAP = Dict(
    "CPU" => ()->Architectures.CPU(),
    "GPU" => ()->Architectures.GPU(),
    "MetalGPU" => ()->Architectures.MetalGPU(),
    # default_architecture is a value (GPU() or CPU()), not a function
    "default_architecture" => ()->Architectures.default_architecture()
)

"""
QUAD_MAP

String → quadrature type mapping. Keys: "RadauQuad", "GaussLegQuad".
"""
const QUAD_MAP = Dict(
    "RadauQuad"           => ()->CoreRT.RadauQuad(),
    "GaussLegQuad" => ()->CoreRT.GaussLegQuad(),
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
    "LambertianSurfaceScalar" => (FT, args)-> _construct_surface(Val(:LambertianSurfaceScalar), FT, args),
    "LambertianSurfaceSpectrum" => (FT, args)-> _construct_surface(Val(:LambertianSurfaceSpectrum), FT, args),
    "LambertianSurfaceLegendre" => (FT, args)-> _construct_surface(Val(:LambertianSurfaceLegendre), FT, args),
    "rpvSurfaceScalar" => (FT, args)-> _construct_surface(Val(:rpvSurfaceScalar), FT, args),
    "RossLiSurfaceScalar" => (FT, args)-> _construct_surface(Val(:RossLiSurfaceScalar), FT, args),
    "CoxMunkSurface" => (FT, args)-> _construct_surface(Val(:CoxMunkSurface), FT, args),
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
- `clumping`: optional scalar constant Ω, "none", or a dict with
  `type: constant|chen_leblanc`
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
    n_leaf    = get(canopy_dict, "n_leaf_quadrature", nothing)
    n_azimuth = get(canopy_dict, "n_azimuth_quadrature", nothing)

    lr = leaf_R isa AbstractVector ? convert(Vector{FT}, leaf_R) : FT(leaf_R)
    lt = leaf_T isa AbstractVector ? convert(Vector{FT}, leaf_T) : FT(leaf_T)

    lg_raw = get(canopy_dict, "leaf_optics_grid", nothing)
    lg = lg_raw === nothing ? nothing : convert(Vector{FT}, lg_raw)
    gu = Symbol(get(canopy_dict, "grid_unit", "nm"))
    canopy_clumping = _parse_canopy_clumping(get(canopy_dict, "clumping", nothing), FT)

    soil_spec = get(canopy_dict, "soil", "from_surface")
    default_canopy_quadrature = CanopyOptics.CanopyQuadrature()
    canopy_quadrature = if n_leaf === nothing && n_azimuth === nothing
        default_canopy_quadrature
    else
        CanopyOptics.CanopyQuadrature(
            n_leaf = n_leaf === nothing ? default_canopy_quadrature.n_leaf : n_leaf,
            n_azimuth = n_azimuth === nothing ? default_canopy_quadrature.n_azimuth : n_azimuth)
    end

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
            canopy_quadrature=canopy_quadrature,
            canopy_clumping=canopy_clumping,
            include_atm=incl_atm, canopy_dp=dp)
    end
    return brdf_list
end

function _parse_canopy_clumping(spec, FT)
    spec === nothing && return CanopyOptics.NoClumping{FT}()

    if spec isa Number
        return CanopyOptics.ConstantClumping{FT}(Ω = spec)
    elseif spec isa AbstractString
        key = lowercase(replace(String(spec), "-" => "_"))
        key in ("none", "no", "no_clumping", "random") &&
            return CanopyOptics.NoClumping{FT}()
        throw(ArgumentError("Unknown canopy clumping string: $spec"))
    elseif spec isa Dict
        clumping_type = lowercase(replace(String(get(spec, "type", "constant")),
                                          "-" => "_"))
        if clumping_type in ("none", "no", "no_clumping", "random")
            return CanopyOptics.NoClumping{FT}()
        elseif clumping_type in ("constant", "constant_clumping")
            Ω = get(spec, "Ω", get(spec, "Omega", get(spec, "Omega0",
                get(spec, "Ω₀", get(spec, "value", 1)))))
            return CanopyOptics.ConstantClumping{FT}(Ω = Ω)
        elseif clumping_type in ("chen_leblanc", "chenleblanc", "chen_leb")
            Ω₀ = get(spec, "Ω₀", get(spec, "Ω", get(spec, "Omega0",
                get(spec, "Omega", get(spec, "value", 0.7)))))
            c = get(spec, "c", 2.0)
            e = get(spec, "e", 2.0)
            return CanopyOptics.ChenLeblancClumping{FT}(Ω₀ = Ω₀, c = c, e = e)
        end
        throw(ArgumentError("Unknown canopy clumping type: $clumping_type"))
    end

    throw(ArgumentError("Unsupported canopy clumping specification: $(typeof(spec))"))
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
    _require_nargs(name, args, n, detail)

Validate the positional arity for an inline surface specification.
"""
_require_nargs(name, args, n, detail) =
    _require_config(length(args) == n, "$(name) expects $(n) argument$(n == 1 ? "" : "s") ($(detail)).")

"""
    _construct_surface(::Val{name}, FT, args)

Construct a surface from an inline YAML surface tag. New surface parsers should
add a method for their tag instead of extending a branch inside
`parse_surface_str`.
"""
function _construct_surface(::Val{:LambertianSurfaceScalar}, FT, args)
    _require_nargs("LambertianSurfaceScalar", args, 1, "albedo")
    return CoreRT.LambertianSurfaceScalar(FT(args[1]))
end

function _construct_surface(::Val{:LambertianSurfaceSpectrum}, FT, args)
    _require_nargs("LambertianSurfaceSpectrum", args, 1, "albedo vector")
    _require_config(args[1] isa AbstractVector, "LambertianSurfaceSpectrum expects 1 vector argument.")
    return CoreRT.LambertianSurfaceSpectrum(convert.(FT, collect(args[1])))
end

function _construct_surface(::Val{:LambertianSurfaceLegendre}, FT, args)
    _require_nargs("LambertianSurfaceLegendre", args, 1, "coefficient vector")
    _require_config(args[1] isa AbstractVector, "LambertianSurfaceLegendre expects 1 vector argument.")
    return CoreRT.LambertianSurfaceLegendre(convert.(FT, collect(args[1])))
end

function _construct_surface(::Val{:rpvSurfaceScalar}, FT, args)
    _require_nargs("rpvSurfaceScalar", args, 4, "ρ₀, ρ_c, k, Θ")
    return CoreRT.rpvSurfaceScalar(FT(args[1]), FT(args[2]), FT(args[3]), FT(args[4]))
end

function _construct_surface(::Val{:RossLiSurfaceScalar}, FT, args)
    _require_nargs("RossLiSurfaceScalar", args, 3, "fvol, fgeo, fiso")
    return CoreRT.RossLiSurfaceScalar(FT(args[1]), FT(args[2]), FT(args[3]))
end

function _construct_surface(::Val{:CoxMunkSurface}, FT, args)
    if length(args) == 1 && args[1] isa AbstractDict
        kw = args[1]
        _require_config(haskey(kw, :wind_speed), "CoxMunkSurface keyword arguments require wind_speed.")
        ws = FT(kw[:wind_speed])
        n_raw = get(kw, :n_water, nothing)
        n_water = n_raw === nothing ? nothing : Complex{FT}(n_raw)
        wc_alb = FT(get(kw, :whitecap_albedo, FT(0.22)))
        inc_wc = Bool(get(kw, :include_whitecaps, true))
        shad   = Bool(get(kw, :shadowing, true))
        return CoreRT.CoxMunkSurface(wind_speed=ws, n_water=n_water,
                                     whitecap_albedo=wc_alb,
                                     include_whitecaps=inc_wc,
                                     shadowing=shad)
    else
        _require_nargs("CoxMunkSurface", args, 1, "wind_speed")
        return CoreRT.CoxMunkSurface(wind_speed=FT(args[1]))
    end
end

_construct_surface(::Val{name}, FT, args) where {name} =
    _config_error("Unknown surface type $(name).")

"""
Parse a surface spec like "LambertianSurfaceScalar(0.1)" or
"CoxMunkSurface(wind_speed=5.0)" into a CoreRT surface instance, without eval.

Supports both positional arguments and keyword arguments (key=value syntax).
When keyword arguments are detected, a `Dict{Symbol,Any}` is passed as the
sole element of the args list to the dispatching surface constructor.
"""
function parse_surface_str(s::AbstractString, FT)
    # Match "Name(args)" or "Name{Type}(args)" -- strip optional type parameter
    name_args = match(r"^(\w+)(?:\{[^}]*\})?\((.*)\)$", strip(String(s)))
    _require_config(name_args !== nothing, "Invalid surface specification: '$(s)'")
    name = name_args.captures[1]
    args_raw = name_args.captures[2]

    if isempty(strip(args_raw))
        return _construct_surface(Val(Symbol(name)), FT, Any[])
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
        return _construct_surface(Val(Symbol(name)), FT, Any[kwargs])
    else
        args_list = [_parse_arg(x, FT) for x in raw_parts]
        return _construct_surface(Val(Symbol(name)), FT, args_list)
    end
end

function _construct_phase_function(::Val{:HenyeyGreensteinPhaseFunction}, FT, args)
    if length(args) == 1 && args[1] isa AbstractDict
        kw = args[1]
        _require_config(haskey(kw, :g), "HenyeyGreensteinPhaseFunction keyword arguments require g.")
        return Scattering.HenyeyGreensteinPhaseFunction(g=FT(kw[:g]))
    end
    _require_nargs("HenyeyGreensteinPhaseFunction", args, 1, "g")
    return Scattering.HenyeyGreensteinPhaseFunction(g=FT(args[1]))
end

function _construct_phase_function(::Val{:SyntheticPolarizedHenyeyGreensteinPhaseFunction},
                                   FT, args)
    if length(args) == 1 && args[1] isa AbstractDict
        kw = args[1]
        _require_config(haskey(kw, :g),
                        "SyntheticPolarizedHenyeyGreensteinPhaseFunction keyword arguments require g.")
        p = get(kw, :polarization_fraction, get(kw, :p, nothing))
        _require_config(p !== nothing,
                        "SyntheticPolarizedHenyeyGreensteinPhaseFunction keyword arguments require polarization_fraction.")
        return Scattering.SyntheticPolarizedHenyeyGreensteinPhaseFunction(
            g=FT(kw[:g]), polarization_fraction=FT(p))
    end
    _require_nargs("SyntheticPolarizedHenyeyGreensteinPhaseFunction", args, 2,
                   "g, polarization_fraction")
    return Scattering.SyntheticPolarizedHenyeyGreensteinPhaseFunction(
        g=FT(args[1]), polarization_fraction=FT(args[2]))
end

_construct_phase_function(::Val{name}, FT, args) where {name} =
    _config_error("Unknown analytic phase function $(name).")

function parse_phase_function_spec(spec, FT)
    if spec isa AbstractDict
        name = get(spec, "type", get(spec, :type, nothing))
        _require_config(name !== nothing, "Analytic phase function dictionaries require a `type` field.")
        kwargs = Dict{Symbol,Any}()
        for (k, v) in spec
            String(k) == "type" && continue
            kwargs[Symbol(k)] = v isa Real ? FT(v) : v
        end
        return _construct_phase_function(Val(Symbol(String(name))), FT,
                                         Any[kwargs])
    elseif spec isa AbstractString
        name_args = match(r"^(\w+)(?:\{[^}]*\})?\((.*)\)$", strip(String(spec)))
        _require_config(name_args !== nothing,
                        "Invalid analytic phase function specification: '$(spec)'")
        name = name_args.captures[1]
        args_raw = name_args.captures[2]
        raw_parts = isempty(strip(args_raw)) ? String[] : _split_args(args_raw)
        has_kwargs = any(p -> occursin('=', p) && !startswith(strip(p), "["),
                         raw_parts)
        if has_kwargs
            kwargs = Dict{Symbol,Any}()
            for part in raw_parts
                p = strip(part)
                if occursin('=', p) && !startswith(p, "[")
                    kv = split(p, '='; limit=2)
                    kwargs[Symbol(strip(kv[1]))] = _parse_arg(strip(kv[2]), FT)
                end
            end
            return _construct_phase_function(Val(Symbol(name)), FT, Any[kwargs])
        end
        return _construct_phase_function(Val(Symbol(name)), FT,
                                         [_parse_arg(x, FT) for x in raw_parts])
    end
    _config_error("Unsupported analytic phase function specification: $(typeof(spec)).")
end

"Check that a field exists in yaml file"
function check_yaml_field(dict::AbstractDict, full_keys::Vector{String}, curr_keys::Vector{String}, final_type, valid_options::Vector{String})
    _require_config(length(curr_keys) >= 1, "Internal error: empty YAML key path.")
    if length(curr_keys) == 1
        _require_config(curr_keys[1] in keys(dict), "Missing key in parameters yaml: $(join(full_keys, '/'))")
        _require_config(dict[curr_keys[1]] isa final_type, "Improper type for $(join(full_keys, '/')); must be a $(final_type)")
        if !isempty(valid_options)
            _require_config(dict[curr_keys[1]] in valid_options, "Field $(join(full_keys, '/')) must be one of $(valid_options)")
        end
        return true
    else
        _require_config(curr_keys[1] in keys(dict), "Missing key in parameters yaml: $(join(full_keys, '/'))")
        return check_yaml_field(dict[curr_keys[1]], full_keys, curr_keys[2:end], final_type, valid_options)
    end
end

"Given a list of aerosol parameter dictionaries, validate all of them.

Vertical distribution may be specified either as (z₀, σ₀) — altitude-form log-normal,
preferred going forward to match sanghavi's YAML convention — or as (p₀, σp) —
pressure-form normal, the legacy unified convention. Exactly one form must be
present per aerosol."
function validate_aerosols(aerosols)
    _require_config(!isempty(aerosols), "Aerosols list shouldn't be empty if scattering block is included")
    for aerosol in aerosols
        check_yaml_field(aerosol, ["τ_ref"], ["τ_ref"], Real, String[])
        has_phase = haskey(aerosol, "phase_function")
        microphysical_keys = ("μ", "σ", "nᵣ", "nᵢ")
        if has_phase
            present = count(k -> haskey(aerosol, k), microphysical_keys)
            _require_config(present == 0 || present == length(microphysical_keys),
                            "Analytic phase-function aerosols must specify either all or none of μ, σ, nᵣ, nᵢ.")
            if haskey(aerosol, "ssa")
                check_yaml_field(aerosol, ["ssa"], ["ssa"], Real, String[])
            elseif haskey(aerosol, "ϖ")
                check_yaml_field(aerosol, ["ϖ"], ["ϖ"], Real, String[])
            end
        else
            base_fields = [(["μ"], Real), (["σ"], Real),
                           (["nᵣ"], Real), (["nᵢ"], Real)]
            for f in base_fields
                key, elty = f[1:2]
                check_yaml_field(aerosol, key, key, elty, String[])
            end
        end
        has_alt  = haskey(aerosol, "z₀") && haskey(aerosol, "σ₀")
        has_pres = haskey(aerosol, "p₀") && haskey(aerosol, "σp")
        _require_config(has_alt || has_pres, "Aerosol must specify vertical distribution as either (z₀, σ₀) [altitude-form, preferred] or (p₀, σp) [pressure-form, legacy]")
        _require_config(!(has_alt && has_pres), "Aerosol must specify exactly one of (z₀, σ₀) or (p₀, σp) — got both")
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
        has_microphysics = all(k -> haskey(aerosol, k), ("μ", "σ", "nᵣ", "nᵢ"))
        if has_microphysics
            _require_config(aerosol["σ"] ≥ 1, "Geometric standard deviation has to be ≥ 1")
        end
        μ = has_microphysics ? FT(aerosol["μ"]) : FT(0.1)
        σ = has_microphysics ? FT(aerosol["σ"]) : FT(1.1)
        size_distribution = LogNormal(log(μ), log(σ))
        nᵣ = has_microphysics ? FT(aerosol["nᵣ"]) : one(FT)
        nᵢ = has_microphysics ? FT(aerosol["nᵢ"]) : zero(FT)
        new_aerosol_obj = Aerosol(size_distribution, nᵣ, nᵢ)
        profile = if haskey(aerosol, "z₀")
            LogNormal(log(FT(aerosol["z₀"])), FT(aerosol["σ₀"]))
        else
            Normal(FT(aerosol["p₀"]), FT(aerosol["σp"]))
        end
        phase_function = haskey(aerosol, "phase_function") ?
                         parse_phase_function_spec(aerosol["phase_function"], FT) :
                         nothing
        ssa = FT(get(aerosol, "ssa", get(aerosol, "ϖ", one(FT))))
        new_rt_aerosol_obj = phase_function === nothing ?
            RT_Aerosol(new_aerosol_obj, FT(aerosol["τ_ref"]), profile) :
            RT_Aerosol(new_aerosol_obj, FT(aerosol["τ_ref"]), profile,
                       phase_function; ϖ=ssa)
        push!(rt_aerosol_obj_list, new_rt_aerosol_obj)
    end
    return rt_aerosol_obj_list
end

"Check that the vmr's in the atmospheric profile match the molecules in the parameters" 
function validate_vmrs(molecules::Array, vmr::Dict)
    for molec in unique(vcat(molecules...))
        _require_config(molec in keys(vmr), "$(molec) listed as molecule in parameters yaml, but no vmr given in atmospheric profile")
        _require_config(vmr[molec] isa Real || vmr[molec] isa Vector, "The vmr for $(molec) in the atmospheric profile must either be a real-valued number, or an array of nodal points from surface to 0hPa (TOA)")
    end
end

"Given a parameter dictionary from a YAML file, validate the dictionary"
function validate_yaml_parameters(params)
    # Required fields — must be present and well-typed.
    required_fields = [
        (["radiative_transfer", "spec_bands"], Array{String}),
        (["radiative_transfer", "surface"], Array{String}),
        (["radiative_transfer", "polarization_type"], String),
        (["radiative_transfer", "Δ_angle"], Real),
        (["radiative_transfer", "depol"], Real),
        (["radiative_transfer", "float_type"], String),
        (["radiative_transfer", "architecture"], String, ["default_architecture", "Architectures.GPU()", "Architectures.MetalGPU()", "Architectures.CPU()", "GPU()", "MetalGPU()", "CPU()"]),
        (["geometry", "sza"], Real),
        (["geometry", "vza"], Array{<:Real}),
        (["geometry", "vaz"], Array{<:Real}),
        (["geometry", "obs_alt"], Real),
        (["atmospheric_profile", "T"], Array{<:Real}),
        (["atmospheric_profile", "p"], Array{<:Real}),
        (["atmospheric_profile", "profile_reduction"], Union{Integer, Nothing}),
    ]
    # Phase D — optional fields. When absent the parser supplies a sensible
    # default (`GaussLegQuad()` for quadrature, `nstreams = 13` for the new
    # schema, etc.). When present, type/value must match. `quadrature_type`
    # also accepts `null` (YAML) → defaults to `GaussLegQuad()`. (Codex
    # Phase D1 P3 finding.)
    optional_fields = [
        (["radiative_transfer", "quadrature_type"], Union{String, Nothing}),  # Phase D: default GaussLegQuad()
        (["radiative_transfer", "max_m"], Integer),              # Legacy
        (["radiative_transfer", "l_trunc"], Integer),            # Legacy
        (["radiative_transfer", "nstreams"], Integer),           # Phase D primary knob
        (["radiative_transfer", "m_max"], Union{Integer, Nothing}),  # Phase D explicit override
    ]
    section_fields = [
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

    for f in required_fields
        key, elty = f[1:2]
        valid_options = length(f) == 3 ? f[3] : String[]
        check_yaml_field(params, key, key, elty, valid_options)
    end

    for f in optional_fields
        key, elty = f[1:2]
        valid_options = length(f) == 3 ? f[3] : String[]
        # Only validate if present.
        if haskey(params, key[1]) && haskey(params[key[1]], key[2])
            check_yaml_field(params, key, key, elty, valid_options)
        end
    end

    for f in section_fields
        key, elty = f[1:2]
        valid_options = length(f) == 3 ? f[3] : String[]
        top = key[1]
        if top == "absorption" && haskey(params, "absorption")
            check_yaml_field(params, key, key, elty, valid_options)
        elseif top == "scattering" && haskey(params, "scattering")
            check_yaml_field(params, key, key, elty, valid_options)
        end
    end

    # New-schema sanity: forbid the user-confusing combination of both
    # `nstreams` and a competing `l_trunc` (they imply different
    # `stream_l_cap` values). `m_max` plus `nstreams` is fine because
    # `m_max` is a tighter override on top of `stream_l_cap`.
    rt = params["radiative_transfer"]
    if haskey(rt, "nstreams") && haskey(rt, "l_trunc")
        @warn "Both `nstreams` and `l_trunc` set in radiative_transfer; `l_trunc` is legacy and will be ignored. Use only `nstreams` for new configs."
    end

    if "scattering" in keys(params)
        validate_aerosols(params["scattering"]["aerosols"])
    end
    if haskey(params, "absorption")
        ab = params["absorption"]
        @assert haskey(ab, "molecules") || haskey(ab, "fixed_molecules") || haskey(ab, "variable_molecules") "absorption section must define `fixed_molecules` and/or `variable_molecules` (or legacy `molecules`)"
    end
end

"Build parameters from a Dict (e.g., parsed YAML)"
function _parse_float_type(params_dict::Dict)
    key = String(params_dict["radiative_transfer"]["float_type"])
    _require_config(haskey(FLOAT_MAP, key), "Unknown float_type $(key).")
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

"""
    _parse_quadrature(params_dict)

Parse `radiative_transfer.quadrature_type` from YAML. Phase D makes the
field optional; missing or `null` defaults to `GaussLegQuad()` (the
recommended scheme — see `docs/dev_notes/fourier_stream_resolution_plan.md`
and the natraj-benchmark observation that Radau is 5–50× less accurate per
stream on Rayleigh-only setups).
"""
function _parse_quadrature(params_dict::Dict)
    rt = params_dict["radiative_transfer"]
    if !haskey(rt, "quadrature_type") || rt["quadrature_type"] === nothing
        return CoreRT.GaussLegQuad()
    end
    qkey = replace(String(rt["quadrature_type"]),"()"=>"")
    _require_config(haskey(QUAD_MAP, qkey), "Unknown quadrature_type $(qkey).")
    return QUAD_MAP[qkey]()
end

function _parse_polarization(params_dict::Dict, FT)
    pol_key = replace(String(params_dict["radiative_transfer"]["polarization_type"]),"()"=>"")
    _require_config(haskey(POLARIZATION_MAP, pol_key), "Unknown polarization_type $(pol_key).")
    return POLARIZATION_MAP[pol_key](FT)
end

function _parse_architecture(params_dict::Dict)
    arch_key = String(params_dict["radiative_transfer"]["architecture"]) |> x->replace(x, "Architectures."=>"") |> x->replace(x, "()"=>"")
    _require_config(haskey(ARCH_MAP, arch_key), "Unknown architecture $(arch_key).")
    return ARCH_MAP[arch_key]()
end

function _parse_atmosphere(params_dict::Dict, FT)
    T = convert.(FT, params_dict["atmospheric_profile"]["T"]) # Level
    p = convert.(FT, params_dict["atmospheric_profile"]["p"]) # Boundaries
    q = haskey(params_dict["atmospheric_profile"], "q") ? convert.(FT, params_dict["atmospheric_profile"]["q"]) : zeros(FT, length(T))
    prof_red = params_dict["atmospheric_profile"]["profile_reduction"]
    return T, p, q, prof_red
end

function _parse_absorption(params_dict::Dict, FT, q=nothing)
    if !haskey(params_dict, "absorption")
        return nothing
    end
    abs_dict = params_dict["absorption"]

    # ── Resolve fixed_molecules / variable_molecules. ─────────────────────
    # Three accepted YAML forms:
    #   (a) `fixed_molecules:` and/or `variable_molecules:` — preferred.
    #   (b) `molecules:` only (legacy) — treated as all-fixed.
    #   (c) `molecules:` + `variable_molecules:` (legacy) — fixed = molecules \ variable.
    # H2O must NOT appear in any of these lists; when q is provided, H2O is
    # auto-included as variable (linearised path) or as fixed (forward path).
    has_fixed = haskey(abs_dict, "fixed_molecules")
    has_var   = haskey(abs_dict, "variable_molecules")
    has_legacy_molecules = haskey(abs_dict, "molecules")
    if has_fixed
        fixed_molecules = Array(abs_dict["fixed_molecules"])
    elseif has_legacy_molecules && has_var
        legacy_mol = Array(abs_dict["molecules"])
        var_legacy = Array(abs_dict["variable_molecules"])
        fixed_molecules = [setdiff(legacy_mol[i],
            i <= length(var_legacy) ? var_legacy[i] : String[])
            for i in 1:length(legacy_mol)]
    elseif has_legacy_molecules
        fixed_molecules = deepcopy(Array(abs_dict["molecules"]))
    else
        fixed_molecules = Vector{Vector{String}}()
    end
    if has_var
        variable_molecules = Array(abs_dict["variable_molecules"])
    else
        variable_molecules = [String[] for _ in eachindex(fixed_molecules)]
    end
    while length(variable_molecules) < length(fixed_molecules)
        push!(variable_molecules, String[])
    end
    while length(fixed_molecules) < length(variable_molecules)
        push!(fixed_molecules, String[])
    end
    n_bands = length(fixed_molecules)
    for ib in 1:n_bands
        _require_config(!any(==("H2O"), fixed_molecules[ib]),
            "H2O must not appear in fixed_molecules (band $ib). H2O is driven by atmospheric_profile.q.")
        _require_config(!any(==("H2O"), variable_molecules[ib]),
            "H2O must not appear in variable_molecules (band $ib). H2O is driven by atmospheric_profile.q.")
    end

    # ── VMR dict: H2O entry (if any) is ignored — we derive from q. ──────
    vmr = convert(Dict{String, Union{Real, Vector}}, abs_dict["vmr"])
    delete!(vmr, "H2O")
    # Validate every fixed/variable species (excluding H2O) has a VMR.
    all_species = unique(vcat(fixed_molecules..., variable_molecules...))
    for sp in all_species
        _require_config(haskey(vmr, sp),
            "$(sp) listed as molecule in YAML but no vmr given in absorption.vmr")
        _require_config(vmr[sp] isa Real || vmr[sp] isa Vector,
            "vmr for $(sp) must be a Real or a Vector")
    end

    # ── Broadening / CEF / wing_cutoff. ─────────────────────────────────────
    bkey = replace(String(abs_dict["broadening"]),"()"=>"")
    _require_config(haskey(BROADENING_MAP, bkey), "Unknown broadening $(bkey).")
    broadening_function = BROADENING_MAP[bkey]()
    ckey = replace(String(abs_dict["CEF"]),"()"=>"")
    _require_config(haskey(CEF_MAP, ckey), "Unknown CEF $(ckey).")
    CEF = CEF_MAP[ckey]()
    wing_cutoff = FT(abs_dict["wing_cutoff"])

    # ── LUTs. Per band, load all LUTfiles. Detect H2O LUT (mol == 1) and
    #    pop it into h2o_lut[band]. The remaining LUTs must be parallel to
    #    `vcat(fixed_molecules[band], variable_molecules[band])`. ───────────
    luts = []
    h2o_lut = Vector{Any}(nothing, n_bands)
    if haskey(abs_dict, "LUTfiles")
        files_lut = Array(abs_dict["LUTfiles"])
        _require_config(length(files_lut) == n_bands,
            "LUTfiles must have one entry per band ($(n_bands)); got $(length(files_lut))")
        for ib in 1:n_bands
            band_models = [load_interpolation_model(_expand_env_path(f)) for f in files_lut[ib]]
            h2o_idx = findfirst(m -> getfield(m, :mol) == 1, band_models)
            if h2o_idx !== nothing
                h2o_lut[ib] = band_models[h2o_idx]
                deleteat!(band_models, h2o_idx)
            end
            expected = length(fixed_molecules[ib]) + length(variable_molecules[ib])
            _require_config(length(band_models) == expected,
                "LUTfiles for band $ib has $(length(band_models)) non-H2O entries but fixed+variable_molecules sum is $(expected)")
            push!(luts, band_models)
        end
    end

    cia_files = haskey(abs_dict, "cia_files") ?
                String.(Array(abs_dict["cia_files"])) :
                String[]
    mtckd_file = haskey(abs_dict, "mtckd_file") ?
                 String(abs_dict["mtckd_file"]) :
                 ""

    return AbsorptionParameters(fixed_molecules, variable_molecules, vmr,
                                broadening_function, CEF, wing_cutoff,
                                luts, h2o_lut, cia_files, mtckd_file)
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
    _require_config(haskey(DECOMP_MAP, dkey), "Unknown decomp_type $(dkey).")
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
    absorption_params = _parse_absorption(params_dict, FT, q)
    scattering_params = _parse_scattering(params_dict, FT)

    if haskey(params_dict, "canopy")
        _parse_canopy_section(params_dict["canopy"], FT, BRDF_per_band)
    end

    # Phase D — resolve the new-vs-legacy schema knobs. See
    # `_resolve_resolution_knobs` for the exact precedence rules.
    res = _resolve_resolution_knobs(params_dict)

    Δ_angle = FT(params_dict["radiative_transfer"]["Δ_angle"])
    truncation = _parse_truncation(params_dict, res.l_trunc, Δ_angle, FT)
    numerics = _parse_numerics(params_dict, FT)

    return vSmartMOM_Parameters(
        spec_bands, BRDF_per_band, quadrature_type, polarization_type,
        res.max_m, Δ_angle,
        res.l_trunc, truncation, FT(params_dict["radiative_transfer"]["depol"]),
        numerics, FT,
        architecture,
        FT(params_dict["geometry"]["sza"]), convert.(FT, params_dict["geometry"]["vza"]),
        convert.(FT, params_dict["geometry"]["vaz"]), FT(params_dict["geometry"]["obs_alt"]),
        T, p, q, profile_reduction,
        absorption_params, scattering_params,
        res.nstreams, res.m_max_override, res.stream_l_cap, res.legacy_l_cap_override,
    )
end

"""
    _resolve_resolution_knobs(params_dict) -> NamedTuple

Phase D — normalize the resolution knobs from YAML into the fields the
struct expects. Distinguishes legacy schema (`max_m`/`l_trunc` present)
from new schema (`nstreams` present). Returns a NamedTuple with:

- `nstreams::Union{Int, Nothing}`        — verbatim if user set it, else nothing
- `m_max_override::Union{Int, Nothing}`  — verbatim `m_max` (order) if set
- `stream_l_cap::Int`                    — derived: `2·nstreams - 1` (new) or
                                            `l_trunc` (legacy)
- `legacy_l_cap_override::Union{Int, Nothing}` — `l_trunc` if user set it
- `max_m::Int`                           — count, populated for back-compat;
                                            derived from `nstreams + 1` or
                                            `m_max + 1` if the user gave one
- `l_trunc::Int`                         — populated for back-compat; derived
                                            from `2·nstreams - 1` if absent

Defaults when the user supplies neither legacy nor new fields:
`nstreams = 13` (new-schema default), `max_m = 14`, `l_trunc = 25`.
"""
function _resolve_resolution_knobs(params_dict::Dict)
    rt = params_dict["radiative_transfer"]

    nstreams_raw   = haskey(rt, "nstreams") ? Int(rt["nstreams"]) : nothing
    m_max_raw      = (haskey(rt, "m_max") && rt["m_max"] !== nothing) ?
                       Int(rt["m_max"]) : nothing
    legacy_max_m   = haskey(rt, "max_m")   ? Int(rt["max_m"])    : nothing
    legacy_l_trunc = haskey(rt, "l_trunc") ? Int(rt["l_trunc"])  : nothing

    # Decide effective `nstreams`.
    nstreams = if nstreams_raw !== nothing
        nstreams_raw
    elseif legacy_max_m === nothing && legacy_l_trunc === nothing
        13   # Phase D default for minimal new-schema configs
    else
        nothing   # legacy schema: leave `nstreams` undefined
    end

    if nstreams !== nothing
        # Solar/scattering minimum: Rayleigh declares Fourier support
        # through m=2, so nstreams ≥ 3 ⇒ stream_l_cap ≥ 5 ≥ 2. Anything
        # smaller silently drops the Rayleigh m=2 term. (Codex Phase D1
        # P2 finding.) Thermal-only / source-only scenes that genuinely
        # only need m=0 will get a separate escape valve once Phase 7
        # ThermalEmission lands.
        nstreams >= 3 || throw(ArgumentError(
            "radiative_transfer.nstreams = $nstreams; must be ≥ 3 for " *
            "solar/scattering scenes (Rayleigh contributes through m=2). " *
            "Streams here are weighted streams per hemisphere — classical " *
            "two-stream uses one stream per hemisphere, vSmartMOM's solar " *
            "scattering minimum is three."))
    end

    # `stream_l_cap` and the legacy `l_trunc` mirror.
    stream_l_cap = if nstreams !== nothing
        2 * nstreams - 1
    else
        # Legacy schema: l_trunc is the projection cap.
        legacy_l_trunc !== nothing ? legacy_l_trunc :
            (legacy_max_m !== nothing ? legacy_max_m : 25)
    end

    # When `nstreams` is the primary knob, override any legacy `l_trunc`
    # the user kept around — passing the larger `l_trunc` to
    # `rt_set_streams` would build more streams than `nstreams` requests.
    # (Codex Phase D1 P2 finding.)
    l_trunc = if nstreams !== nothing
        stream_l_cap
    elseif legacy_l_trunc !== nothing
        legacy_l_trunc
    else
        stream_l_cap
    end

    # `max_m` (legacy count). When `nstreams` drives the schema, the
    # back-compat count must reflect `stream_l_cap + 1` so the trait
    # aggregator's `legacy_order_cap = max_m - 1` matches `stream_l_cap`.
    # Otherwise `nstreams: 8` would silently truncate Cox-Munk at
    # `max_m - 1 = nstreams = 8` instead of the documented
    # `stream_l_cap = 15`. (Codex Phase D1 P1 finding.)
    max_m = if legacy_max_m !== nothing
        legacy_max_m
    elseif m_max_raw !== nothing
        m_max_raw + 1
    elseif nstreams !== nothing
        stream_l_cap + 1
    else
        14
    end

    return (; nstreams,
            m_max_override = m_max_raw,
            stream_l_cap,
            legacy_l_cap_override = legacy_l_trunc,
            max_m, l_trunc)
end

"""
    _parse_numerics(params_dict, FT)

Build the `RTNumericalParameters` from the optional `radiative_transfer.numerics`
YAML subsection. Missing keys fall back to type defaults (defined on the
struct in `src/CoreRT/types.jl`). Existing configs without a `numerics`
section keep their legacy values.
"""
function _parse_numerics(params_dict, FT)
    rt = params_dict["radiative_transfer"]
    if !haskey(rt, "numerics") || rt["numerics"] === nothing
        return RTNumericalParameters{FT}()
    end
    n = rt["numerics"]
    kwargs = Dict{Symbol, Any}()
    if haskey(n, "dτ_max_threshold") || haskey(n, "dtau_max_threshold")
        v = haskey(n, "dτ_max_threshold") ? n["dτ_max_threshold"] : n["dtau_max_threshold"]
        kwargs[:dτ_max_threshold] = FT(v)
    end
    if haskey(n, "dτ_min_floor") || haskey(n, "dtau_min_floor")
        v = haskey(n, "dτ_min_floor") ? n["dτ_min_floor"] : n["dtau_min_floor"]
        kwargs[:dτ_min_floor] = FT(v)
    end
    if haskey(n, "blas_threads")
        v = n["blas_threads"]
        # Allow `null` / `~` in YAML to mean "leave BLAS alone"; otherwise
        # parse to Int and pass through as `RTNumericalParameters.blas_threads`.
        kwargs[:blas_threads] = v === nothing ? nothing : Int(v)
    end
    return RTNumericalParameters{FT}(; kwargs...)
end

"""
    _parse_truncation(params_dict, l_trunc, Δ_angle, FT)

Resolve the truncation method from the YAML/dict config. The
new-vs-legacy distinguisher is the presence of `nstreams` —
configs that opt into the v0.7 schema get the documented `auto`
default and treat explicit `null` as `NoTruncation()` (per
docs/src/pages/IO/Schema/radiative_transfer.md). Legacy configs
(`max_m`/`l_trunc` set) keep the historical `δBGE(l_trunc, Δ_angle)`
default for backward compatibility.

| YAML form                    | new schema (nstreams set) | legacy schema |
|------------------------------|---------------------------|---------------|
| `truncation:` *(omitted)*    | `AutoTruncation()`        | `δBGE(...)`   |
| `truncation: null`           | `NoTruncation()`          | `δBGE(...)`   |
| `truncation: auto`           | `AutoTruncation()`        | same          |
| `truncation: NoTruncation()` | `NoTruncation()`          | same          |
| `truncation: δBGE(L, Δ)`     | `δBGE(L, Δ)`              | same          |

YAML/TOML string values are matched against a fixed allow-list of
constructors — `Meta.parse` + `eval` is not used, so untrusted configs
can't execute arbitrary Julia code from this field.
"""
function _parse_truncation(params_dict, l_trunc, Δ_angle, FT)
    rt = params_dict["radiative_transfer"]
    is_new_schema = haskey(rt, "nstreams")

    if haskey(rt, "truncation")
        spec = rt["truncation"]
        if spec === nothing
            # Explicit `null` — under new schema this means NoTruncation()
            # (the schema doc promises "exactly no transform; errors if
            # coefs exceed cap"). Under legacy schema we preserve the
            # historical δBGE fallback.
            return is_new_schema ? Scattering.NoTruncation() :
                                    Scattering.δBGE{FT}(l_trunc, Δ_angle)
        elseif spec isa Scattering.AbstractTruncationType
            return spec
        elseif spec isa AbstractString
            return _truncation_from_string(spec, FT)
        else
            throw(ArgumentError("radiative_transfer.truncation must be a string " *
                                "or AbstractTruncationType, got $(typeof(spec))"))
        end
    end
    # Field omitted entirely: new schema → AutoTruncation (deferred
    # decision); legacy → δBGE(l_trunc, Δ_angle) for back-compat.
    return is_new_schema ? Scattering.AutoTruncation() :
                            Scattering.δBGE{FT}(l_trunc, Δ_angle)
end

"""
    _truncation_from_string(spec, FT) -> AbstractTruncationType

Parse a YAML truncation spec string against a fixed allow-list of
constructor shapes:

* `"NoTruncation()"`
* `"NoTruncation(l_max)"` or `"NoTruncation(l_max=N)"`
* `"δBGE(l_max, Δ_angle)"` or `"δBGE{Float64}(l_max, Δ_angle)"`

Anything else throws `ArgumentError`. No `eval` — string matching only,
so untrusted configs can't execute arbitrary code.
"""
function _truncation_from_string(spec::AbstractString, ::Type{FT}) where {FT}
    s = strip(spec)
    # Phase D — `auto` is the deferred-decision marker; resolved at
    # `model_from_parameters` time per the per-band phase_lmax check.
    if s == "auto" || s == "AutoTruncation()"
        return Scattering.AutoTruncation()
    end
    if s == "NoTruncation()"
        return Scattering.NoTruncation()
    end
    # NoTruncation(N) or NoTruncation(l_max=N)
    m = match(r"^NoTruncation\(\s*(?:l_max\s*=\s*)?(\d+)\s*\)$", s)
    if m !== nothing
        return Scattering.NoTruncation(l_max = parse(Int, m.captures[1]))
    end
    # δBGE(l_max, Δ_angle) or δBGE{TYPE}(l_max, Δ_angle)
    m = match(r"^δBGE(?:\{[^}]+\})?\(\s*(\d+)\s*,\s*([+-]?\d+(?:\.\d*)?(?:[eE][+-]?\d+)?)\s*\)$", s)
    if m !== nothing
        return Scattering.δBGE{FT}(parse(Int, m.captures[1]),
                                   parse(FT, m.captures[2]))
    end
    throw(ArgumentError(
        "radiative_transfer.truncation = $(repr(spec)) does not match a " *
        "supported constructor. Allowed: \"auto\", \"NoTruncation()\", " *
        "\"NoTruncation(l_max=N)\", \"δBGE(l_max, Δ_angle)\", " *
        "\"δBGE{T}(l_max, Δ_angle)\"."))
end

# Convenience wrappers
"""
read_parameters(x)

Multi-dispatch convenience wrapper:
- String (file path) → parameters_from_file
- IOSource → parameters_from_source
- Dict → parameters_from_dict
- vSmartMOM_Parameters → pass-through
"""
read_parameters(x::AbstractString) = parameters_from_file(x)
read_parameters(x::Formats.IOSource) = parameters_from_source(x)
read_parameters(x::Dict) = parameters_from_dict(x)
read_parameters(p::vSmartMOM_Parameters) = p
