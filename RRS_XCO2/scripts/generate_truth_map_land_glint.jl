#!/usr/bin/env julia

"""
Generate one no-SIF land-glint subset without modifying the nadir truth map.

Environment controls:

- `GLINT_STREAM_INDEX`: one of the nine-stream Gauss-Legendre indices 5:9
  (default 5). These are the native directions from 60 down to 10.237 degrees;
  the solar and viewing zenith angles are both set to the selected direction.
- `GLINT_ANGLE_DEG=60` remains accepted for the already-running 60-degree
  production job. New jobs should select a native stream index.
- `GLINT_PHASE`: `aerosol` (default) or `no_aerosol`.
- `GLINT_OUT_ROOT`: default `RRS_XCO2/truth_map/land_glint`.

The relative azimuth is fixed to 0 degrees in the vSmartMOM convention used by
this experiment. Aerosol scenes use the resumable Raman/chunked production
driver; aerosol-free scenes use the full-band driver. Both retain the original
truth-map state indices and select only the SIF-off rows. Aerosol calculations
are fixed to nine streams and 256 retained O2 points (11 spectral chunks).
"""

using Printf

const LAND_GLINT_PHASE = lowercase(get(ENV, "GLINT_PHASE", "aerosol"))
const LAND_GLINT_NATIVE_MU = Dict(
    5 => Float32(0.500000000000000),
    6 => Float32(0.662126711701905),
    7 => Float32(0.806685716350295),
    8 => Float32(0.918015553663318),
    9 => Float32(0.984080119753813),
)

"Return a Float32 angle whose `cosd` is exactly the selected Float32 node."
function exact_native_angle(mu::Float32)
    theta0 = acosd(mu)
    for offset in 0:16
        candidates = offset == 0 ? (theta0,) :
            (prevfloat(theta0, offset), nextfloat(theta0, offset))
        for theta in candidates
            cosd(theta) == mu && return theta
        end
    end
    error("could not round-trip native quadrature node mu=$mu through VZA")
end

const LAND_GLINT_STREAM_INDEX = if haskey(ENV, "GLINT_STREAM_INDEX")
    parse(Int, ENV["GLINT_STREAM_INDEX"])
elseif get(ENV, "GLINT_ANGLE_DEG", "60") == "60"
    # Compatibility for the active 60-degree job and its queued second phase.
    5
else
    error("use GLINT_STREAM_INDEX=5,6,7,8,or9 for native land-glint angles")
end

haskey(LAND_GLINT_NATIVE_MU, LAND_GLINT_STREAM_INDEX) || error(
    "GLINT_STREAM_INDEX must be one of 5, 6, 7, 8, or 9")
LAND_GLINT_PHASE in ("aerosol", "no_aerosol") || error(
    "GLINT_PHASE must be aerosol or no_aerosol")

const LAND_GLINT_MU = LAND_GLINT_NATIVE_MU[LAND_GLINT_STREAM_INDEX]
const LAND_GLINT_ANGLE = exact_native_angle(LAND_GLINT_MU)
angle_tag = LAND_GLINT_STREAM_INDEX == 5 ? "60" :
    replace(@sprintf("%.6f", Float64(LAND_GLINT_ANGLE)), "." => "p")
geometry_tag = "sza$(angle_tag)_vza$(angle_tag)_relaz00"
root = normpath(joinpath(@__DIR__, ".."))
out_root = get(ENV, "GLINT_OUT_ROOT", joinpath(root, "truth_map", "land_glint"))
phase_dir = LAND_GLINT_PHASE == "aerosol" ? "aerosol_chunked" : "no_aerosol"

ENV["TRUTH_SZA_DEG"] = string(LAND_GLINT_ANGLE)
ENV["TRUTH_VZA_DEG"] = string(LAND_GLINT_ANGLE)
ENV["TRUTH_RELATIVE_AZIMUTH_DEG"] = "0"
ENV["SIF_CASE_FILTER"] = "off"
ENV["AEROSOL_CASE_FILTER"] = LAND_GLINT_PHASE == "aerosol" ? "aerosol" : "none"
ENV["TRUTH_OUT"] = joinpath(out_root, geometry_tag, phase_dir)
ENV["AEROSOL_NSTREAMS"] = "9"
ENV["O2_CHUNK_POINTS"] = "256"

if LAND_GLINT_PHASE == "aerosol"
    include(joinpath(@__DIR__, "generate_truth_map_aerosol_chunked.jl"))
    main_aerosol_chunked()
else
    include(joinpath(@__DIR__, "generate_truth_map.jl"))
    main()
end
