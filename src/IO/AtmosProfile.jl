# Atmospheric profile IO

using YAML
using Interpolations
using Parameters: @unpack
using ..CoreRT: AtmosphericProfile
import ..CoreRT: compute_atmos_profile_fields

"Read atmospheric profile from a parameters Dict"
function read_atmos_profile_dict(params_dict::Dict)
    T = convert.(Float64, params_dict["T"])
    if ("ak" in keys(params_dict))
        psurf = convert(Float64, params_dict["p_surf"])
        q     = convert.(Float64, params_dict["q"])
        ak    = convert.(Float64, params_dict["ak"])
        bk    = convert.(Float64, params_dict["bk"])
        p_half = (ak + bk * psurf)
        p_full, p_half, vmr_h2o, vcd_dry, vcd_h2o, Δz = compute_atmos_profile_fields(T, p_half, q, Dict())
    elseif ("q" in keys(params_dict))
        p_half = convert.(Float64, params_dict["p_half"])
        psurf = p_half[end]
        q      = convert.(Float64, params_dict["q"])
        p_full, p_half, vmr_h2o, vcd_dry, vcd_h2o, Δz = compute_atmos_profile_fields(T, p_half, q, Dict())
    else
        p_half = convert.(Float64, params_dict["p_half"])
        psurf = p_half[end]
        q = zeros(length(T))
        p_full, p_half, vmr_h2o, vcd_dry, vcd_h2o, Δz = compute_atmos_profile_fields(T, p_half, q, Dict())
    end
    vmr = convert(Dict{String, Union{Real, Vector}}, params_dict["vmr"])
    return AtmosphericProfile(T, q, p_full, p_half, vmr_h2o, vcd_dry, vcd_h2o, vmr, Δz)
end

"Load atmospheric profile from YAML file path"
function read_atmos_profile(file_path::AbstractString)
    @assert endswith(file_path, ".yaml") "File must be yaml"
    params_dict = YAML.load_file(file_path)
    return read_atmos_profile_dict(params_dict)
end
