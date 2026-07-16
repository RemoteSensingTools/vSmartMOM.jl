# Atmospheric profile IO

using YAML
using Interpolations
using Parameters: @unpack
using ..CoreRT: AtmosphericProfile, _layer_centered_input, _layer_centered_vmr
import ..CoreRT: compute_atmos_profile_fields

"""
    _atmos_profile_error(msg)

Raise a stable `ArgumentError` for invalid atmospheric-profile input.
"""
@inline _atmos_profile_error(msg) = throw(ArgumentError(msg))

"""
    _require_atmos_profile(cond, msg)

Validate atmospheric-profile input and raise `ArgumentError` when invalid.
"""
@inline _require_atmos_profile(cond, msg) = cond ? nothing : _atmos_profile_error(msg)

"Normalize and validate the trace-gas mapping in a standalone profile."
function _read_profile_vmr(raw)
    _require_atmos_profile(raw isa AbstractDict,
                           "vmr must be a mapping from gas names to numbers or vectors")
    vmr = Dict{String, Union{Real, Vector}}()
    for (name, value) in raw
        if value isa Real && !(value isa Bool)
            value_ft = Float64(value)
            _require_atmos_profile(isfinite(value_ft), "vmr/$name must be finite")
            vmr[string(name)] = value_ft
        elseif value isa AbstractVector && all(x -> x isa Real && !(x isa Bool), value)
            value_ft = Float64.(value)
            _require_atmos_profile(!isempty(value_ft) && all(isfinite, value_ft),
                                   "vmr/$name must be a nonempty vector of finite numbers")
            vmr[string(name)] = value_ft
        else
            _atmos_profile_error("vmr/$name must be a number or a vector of numbers")
        end
    end
    return vmr
end

"Read atmospheric profile from a parameters dictionary."
function read_atmos_profile_dict(params_dict::AbstractDict)
    _require_atmos_profile(haskey(params_dict, "T"), "Atmospheric profile requires T")
    T = Float64.(params_dict["T"])
    _require_atmos_profile(!isempty(T) && all(isfinite, T),
                           "T must be a nonempty vector of finite values")

    if haskey(params_dict, "ak") || haskey(params_dict, "bk")
        _require_atmos_profile(haskey(params_dict, "ak") && haskey(params_dict, "bk") &&
                               haskey(params_dict, "p_surf"),
                               "Hybrid profiles require ak, bk, and p_surf")
        ak = Float64.(params_dict["ak"])
        bk = Float64.(params_dict["bk"])
        _require_atmos_profile(length(ak) == length(bk),
                               "ak and bk must have the same length")
        p_half = ak .+ bk .* Float64(params_dict["p_surf"])
    else
        _require_atmos_profile(haskey(params_dict, "p_half"),
                               "Atmospheric profile requires p_half")
        p_half = Float64.(params_dict["p_half"])
    end

    q_input = haskey(params_dict, "q") ? Float64.(params_dict["q"]) : zeros(length(T))
    _require_atmos_profile(length(p_half) == length(T) + 1,
                           "p_half must contain one more boundary than T has layers")
    _require_atmos_profile(all(isfinite, p_half) && all(>(0), p_half) &&
                           all(>(0), diff(p_half)),
                           "p_half must be finite, positive, and ordered from TOA to BOA")
    q = _layer_centered_input("specific humidity", q_input, p_half)
    _require_atmos_profile(all(x -> isfinite(x) && 0 <= x < 1, q),
                           "q must contain finite specific humidities in [0, 1)")

    vmr = _layer_centered_vmr(
        _read_profile_vmr(get(params_dict, "vmr", Dict{String, Any}())), p_half)
    p_full, p_half, vmr_h2o, vcd_dry, vcd_h2o, vmr, Δz =
        compute_atmos_profile_fields(T, p_half, q, vmr)
    return AtmosphericProfile(T, p_full, q, p_half, vmr_h2o,
                              vcd_dry, vcd_h2o, vmr, Δz)
end

"Load atmospheric profile from YAML file path"
function read_atmos_profile(file_path::AbstractString)
    _require_atmos_profile(endswith(file_path, ".yaml"), "File must be yaml")
    params_dict = YAML.load_file(file_path)
    return read_atmos_profile_dict(params_dict)
end
