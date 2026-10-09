#!/usr/bin/env julia

"""
Export the canonical materialized 16-layer atmosphere used by the RRS-XCO2
truth and retrieval forward models.

The plotting code needs the actual post-materialization pressure, temperature,
and humidity arrays to convert retrieved CO2 VMRs into molecular columns and
layer-mean number concentrations.  Reconstructing the profile independently
from the raw YAML would duplicate the model's interpolation and reduction
logic, so this script deliberately calls the same shared constructor as the
retrieval.

Usage
-----

    julia --project=. \
        RRS_XCO2/inversion/retrieval_setup/export_retrieval_atmosphere.jl

An optional first argument overrides the output path.  The default is
`retrieval_atmosphere_16layer.nc` beside this script.
"""

using Dates
using NCDatasets
using vSmartMOM
using vSmartMOM.CoreRT

include(joinpath(@__DIR__, "..", "..", "scripts", "common.jl"))
using .RRSXCO2Common

const DEFAULT_OUTPUT = joinpath(@__DIR__, "retrieval_atmosphere_16layer.nc")
const REFERENCE_PSURF_HPA = 1000.0
const NLAYERS = 16

function canonical_profile()
    params = RRSXCO2Common.load_parameters(;
        float_type=Float32, architecture=:CPU, nstreams=9)
    RRSXCO2Common.prepare_shared_profile!(
        params; psurf=REFERENCE_PSURF_HPA, nlayers=NLAYERS)

    # Model construction reframes this already-materialized p/T/q profile once
    # more with profile_reduction=-1.  Calling the same entry point here makes
    # this artifact match the actual retrieval atmosphere, including Float32
    # interpolation and hydrostatic arithmetic.
    profile, _ = CoreRT.prepare_observer_profile(
        Vector{Float32}(params.T), Vector{Float32}(params.p),
        Vector{Float32}(params.q),
        Dict("CO2" => fill(Float32(400e-6), NLAYERS)),
        Float32[0], -1)
    length(profile.T) == NLAYERS || error(
        "expected $NLAYERS layers, got $(length(profile.T))")
    return profile
end

function write_profile(path::AbstractString)
    profile = canonical_profile()
    z_half = CoreRT.half_level_altitudes(profile)
    NCDataset(path, "c") do output
        defDim(output, "layer", NLAYERS)
        defDim(output, "interface", NLAYERS + 1)

        function layer_variable(name, values, units, long_name)
            variable = defVar(output, name, Float64, ("layer",))
            variable.attrib["units"] = units
            variable.attrib["long_name"] = long_name
            variable[:] = Float64.(values)
            return variable
        end
        function interface_variable(name, values, units, long_name)
            variable = defVar(output, name, Float64, ("interface",))
            variable.attrib["units"] = units
            variable.attrib["long_name"] = long_name
            variable[:] = Float64.(values)
            return variable
        end

        interface_variable(
            "pressure_interface", profile.p_half, "hPa",
            "reference pressure interfaces, TOA to BOA")
        layer_variable(
            "pressure_center", profile.p_full, "hPa",
            "arithmetic layer-center pressure, TOA to BOA")
        layer_variable(
            "temperature", profile.T, "K",
            "materialized layer temperature, TOA to BOA")
        layer_variable(
            "specific_humidity", profile.q, "kg kg-1",
            "materialized specific humidity, TOA to BOA")
        layer_variable(
            "h2o_dry_molar_ratio", profile.vmr_h2o, "mol mol-1 dry air",
            "water-vapor to dry-air molar ratio")
        layer_variable(
            "dry_air_vertical_column", profile.vcd_dry, "molecules cm-2",
            "reference dry-air molecular column per layer")
        layer_variable(
            "layer_thickness", profile.Δz, "m",
            "reference geometric layer thickness")
        interface_variable(
            "altitude_interface", z_half, "km above BOA",
            "reference geometric interfaces, TOA to BOA")

        output.attrib["reference_surface_pressure_hpa"] =
            REFERENCE_PSURF_HPA
        output.attrib["layer_order"] = "TOA to BOA"
        output.attrib["profile_configuration"] =
            abspath(RRSXCO2Common.CONFIG)
        output.attrib["profile_construction"] =
            "RRSXCO2Common.prepare_shared_profile!(psurf=1000,nlayers=16), " *
            "then CoreRT.prepare_observer_profile(profile_reduction=-1)"
        output.attrib["retrieval_psurf_mapping"] =
            "replace only pressure_interface[end]; retain temperature and " *
            "specific_humidity; recompute VCD, thickness, and altitude"
        output.attrib["float_type"] = "Float32 (exported as Float64 values)"
        output.attrib["created_utc"] = string(now(UTC))
        output.attrib["complete"] = Int8(1)
    end
    return path
end

output = isempty(ARGS) ? DEFAULT_OUTPUT : abspath(ARGS[1])
println(write_profile(output))
