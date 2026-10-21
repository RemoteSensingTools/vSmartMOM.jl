using CUDA
using DelimitedFiles
using vSmartMOM

const SOURCE_CONFIG = "/home/cfranken/.codex/attachments/863294a0-bb39-460e-950d-063f515085bf/pasted-text.txt"
const OUTPUT = joinpath(@__DIR__, "so4_fixed_vza0to80_alt1p67_spectra.csv")
const SO4_INDEX = 7
const TARGET_AOD = 0.03
const OBSERVER_ALTITUDE_KM = 1.67
const VIEWING_ZENITH_ANGLES = collect(0.0:5.0:80.0)

function load_student_config()
    source = read(SOURCE_CONFIG, String)
    expression, _ = Meta.parse(source, 1; greedy=true, raise=true)
    expression.head == :(=) || error("Expected the first expression to assign cfg")
    return Core.eval(Main, expression)
end

function prepare_config()
    cfg = deepcopy(load_student_config())
    rt = cfg["radiative_transfer"]
    aerosols = cfg["scattering"]["aerosols"]

    for aerosol in aerosols
        aerosol["τ_ref"] = 0.0
    end
    aerosols[SO4_INDEX]["τ_ref"] = TARGET_AOD

    # Fixed optical setup, retaining the original eight-stream quadrature.
    rt["nstreams"] = 8
    rt["truncation"] = "auto"
    rt["Δ_angle"] = 2.0
    delete!(rt, "max_m")
    delete!(rt, "l_trunc")
    cfg["scattering"]["aerosols"] = [aerosols[SO4_INDEX]]

    cfg["geometry"]["vza"] = VIEWING_ZENITH_ANGLES
    cfg["geometry"]["vaz"] = zeros(length(VIEWING_ZENITH_ANGLES))
    # A vector containing zero requests the endpoint outputs in addition to
    # the strict-interior 1.67 km interface.
    cfg["geometry"]["obs_alt"] = [0.0, OBSERVER_ALTITUDE_KM]
    return cfg
end

function extract_upwelling(result)
    result.toa === nothing && error("TOA output was not returned")
    length(result.levels) == 1 || error("Expected one interior observer level")
    toa = Array(result.toa[:, 1, :])
    altitude = Array(only(result.levels).upwelling[:, 1, :])
    return toa, altitude
end

CUDA.device!(0)
CUDA.allowscalar(false)

cfg = prepare_config()
println("AEROSOL=SO4; CASE=fixed; H2O=profile q")
println("NSTREAMS=8")
println("VZA=", VIEWING_ZENITH_ANGLES)
println("observer outputs=[TOA, $(OBSERVER_ALTITUDE_KM) km above BOA]")

build = @timed begin
    params = read_parameters(cfg)
    model = model_from_parameters(params)
    (params, model)
end
params, model = build.value
println("model build seconds: ", build.time)
println("solver: ", model.solver)
println("resolved interior altitude [km above BOA]: ",
        only(model.obs_geom.sensor_altitudes))
println("resolved pressure [hPa]: ",
        model.atmosphere.profile.p_half[only(model.obs_geom.sensor_levels) + 1])

τ03 = copy(model.optics.aerosols.τ_aer[1])
fill!(model.optics.aerosols.τ_aer[1], 0)
zero_run = @timed rt_run(model)
toa0, altitude0 = extract_upwelling(zero_run.value)
println("AOD=0 RT seconds: ", zero_run.time)

model.optics.aerosols.τ_aer[1] .= τ03
aerosol_run = @timed rt_run(model)
toa03, altitude03 = extract_upwelling(aerosol_run.value)
println("AOD=0.03 RT seconds: ", aerosol_run.time)

ν = collect(params.spec_bands[1])
size(toa0) == (length(VIEWING_ZENITH_ANGLES), length(ν)) ||
    error("Unexpected TOA output shape $(size(toa0))")
size(altitude0) == size(toa0) || error("TOA/interior shape mismatch")

open(OUTPUT, "w") do io
    println(io, "vza_deg,wavenumber_cm-1,toa_aod0,toa_aod003,alt1p67_aod0,alt1p67_aod003")
    for ivza in eachindex(VIEWING_ZENITH_ANGLES), iν in eachindex(ν)
        writedlm(io, reshape([
            VIEWING_ZENITH_ANGLES[ivza], ν[iν],
            toa0[ivza, iν], toa03[ivza, iν],
            altitude0[ivza, iν], altitude03[ivza, iν],
        ], 1, :), ',')
    end
end

println("TOA reverse-ratio extrema by VZA:")
for (ivza, vza) in enumerate(VIEWING_ZENITH_ANGLES)
    println("  ", vza, "°: ", extrema(toa0[ivza, :] ./ toa03[ivza, :]))
end
println("1.67 km reverse-ratio extrema by VZA:")
for (ivza, vza) in enumerate(VIEWING_ZENITH_ANGLES)
    println("  ", vza, "°: ", extrema(altitude0[ivza, :] ./ altitude03[ivza, :]))
end
println("wrote ", OUTPUT)
