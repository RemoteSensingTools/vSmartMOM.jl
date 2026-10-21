using CUDA
using DelimitedFiles
using Statistics
using vSmartMOM

const SOURCE_CONFIG = "/home/cfranken/.codex/attachments/863294a0-bb39-460e-950d-063f515085bf/pasted-text.txt"
const OUTPUT_DIR = @__DIR__
const SO4_INDEX = 7
const TARGET_AOD = 0.03

function load_student_config()
    source = read(SOURCE_CONFIG, String)
    expression, _ = Meta.parse(source, 1; greedy=true, raise=true)
    expression.head == :(=) || error("Expected the first expression to assign cfg")
    return Core.eval(Main, expression)
end

function prepare_config(kind::String, observer_altitude, dry_atmosphere::Bool,
                        vza_zero::Bool, nstreams_30::Bool)
    cfg = deepcopy(load_student_config())
    rt = cfg["radiative_transfer"]
    aerosols = cfg["scattering"]["aerosols"]
    length(aerosols) >= SO4_INDEX || error("Expected sulfate at aerosol index 7")

    # Preserve the student's sulfate microphysics but make it the only loaded
    # aerosol: all other aerosol optical depths are exactly zero.
    for aerosol in aerosols
        aerosol["τ_ref"] = 0.0
    end
    aerosols[SO4_INDEX]["τ_ref"] = TARGET_AOD

    if dry_atmosphere
        fill!(cfg["atmospheric_profile"]["q"], 0.0)
    end

    if observer_altitude !== nothing
        # suniti_multi_sensor interprets this as km above model BOA.
        cfg["geometry"]["obs_alt"] = [observer_altitude]
    end

    if vza_zero
        cfg["geometry"]["vza"] = [0.0]
    end

    if nstreams_30
        rt["nstreams"] = 30
    end

    if kind == "fixed"
        rt["truncation"] = "auto"
        rt["Δ_angle"] = 2.0
        delete!(rt, "max_m")
        delete!(rt, "l_trunc")
        cfg["scattering"]["aerosols"] = [aerosols[SO4_INDEX]]
    elseif kind != "original"
        error("case must be original or fixed")
    end

    # For the 30-stream convergence test, keep the fixed case at the same
    # Fourier cap used by its 8-stream baseline. Otherwise the automatic
    # 2N-1 cap would change both quadrature and azimuthal convergence at once.
    if nstreams_30 && kind == "fixed"
        rt["max_m"] = 15
    end
    return cfg
end

function requested_spectra(result, observer_altitude)
    if observer_altitude === nothing
        return vec(Array(result.toa[1, 1, :])), vec(Array(result.boa[1, 1, :]))
    end
    level = only(result.levels)
    return vec(Array(level.upwelling[1, 1, :])),
           vec(Array(total_downwelling(level)[1, 1, :]))
end

function finite_extrema(x)
    values = filter(isfinite, x)
    return extrema(values)
end

function run_case(kind::String, observer_altitude, dry_atmosphere::Bool,
                  vza_zero::Bool, nstreams_30::Bool)
    CUDA.device!(0)
    CUDA.allowscalar(false)
    cfg = prepare_config(
        kind, observer_altitude, dry_atmosphere, vza_zero, nstreams_30
    )
    so4 = kind == "fixed" ? only(cfg["scattering"]["aerosols"]) :
                            cfg["scattering"]["aerosols"][SO4_INDEX]

    println("AEROSOL=SO4")
    println("H2O=", dry_atmosphere ? "zero" : "profile q")
    println("VZA=", vza_zero ? "0 degrees" : "student configuration")
    println("NSTREAMS=", cfg["radiative_transfer"]["nstreams"])
    println("CASE=$kind")
    println("CUDA device: ", CUDA.name(CUDA.device()))
    println("SO4 properties: ", so4)
    println("radiative_transfer: ", cfg["radiative_transfer"])
    println("configured aerosol count: ", length(cfg["scattering"]["aerosols"]))

    build = @timed begin
        params = read_parameters(cfg)
        model = model_from_parameters(params)
        (params, model)
    end
    params, model = build.value
    println("model build seconds: ", build.time)
    println("solver: ", model.solver)
    if observer_altitude !== nothing
        boundary = only(model.obs_geom.sensor_levels)
        println("resolved observer altitude [km above BOA]: ",
                only(model.obs_geom.sensor_altitudes))
        println("resolved observer boundary index: ", boundary)
        println("resolved observer pressure [hPa]: ",
                model.atmosphere.profile.p_half[boundary + 1])
    end

    τ03 = copy(model.optics.aerosols.τ_aer[1])
    fill!(model.optics.aerosols.τ_aer[1], 0)
    zero_run = @timed rt_run(model)
    up0, down0 = requested_spectra(zero_run.value, observer_altitude)
    println("AOD=0 RT seconds: ", zero_run.time)

    model.optics.aerosols.τ_aer[1] .= τ03
    aerosol_run = @timed rt_run(model)
    up03, down03 = requested_spectra(aerosol_run.value, observer_altitude)
    println("AOD=0.03 RT seconds: ", aerosol_run.time)

    ν = collect(params.spec_bands[1])
    length(ν) == length(up0) || error("Spectral grid/output length mismatch")
    table = hcat(ν, up0, up03, down0, down03)
    suffix = observer_altitude === nothing ? "" :
             "_alt" * replace(string(observer_altitude), "." => "p")
    dry_suffix = dry_atmosphere ? "_dry" : ""
    vza_suffix = vza_zero ? "_vza0" : ""
    streams_suffix = nstreams_30 ? "_nstreams30" : ""
    output = joinpath(
        OUTPUT_DIR,
        "so4_$(kind)$(suffix)$(dry_suffix)$(vza_suffix)$(streams_suffix)_spectra.csv",
    )
    open(output, "w") do io
        println(io, "wavenumber_cm-1,up_aod0,up_aod003,down_aod0,down_aod003")
        writedlm(io, table, ',')
    end

    up_ratio = up03 ./ up0
    down_ratio = down03 ./ down0
    println("upwelling ratio extrema: ", finite_extrema(up_ratio))
    println("downwelling ratio extrema: ", finite_extrema(down_ratio))
    println("upwelling ratio mean: ", mean(filter(isfinite, up_ratio)))
    println("wrote ", output)
end

1 <= length(ARGS) <= 5 ||
    error("Usage: julia run_so4_case.jl original|fixed [observer_altitude_km] [dry] [vza0] [nstreams30]")
dry_atmosphere = any(==("dry"), ARGS[2:end])
vza_zero = any(==("vza0"), ARGS[2:end])
nstreams_30 = any(==("nstreams30"), ARGS[2:end])
altitude_args = filter(
    arg -> arg != "dry" && arg != "vza0" && arg != "nstreams30",
    ARGS[2:end],
)
length(altitude_args) <= 1 || error("At most one observer altitude may be supplied")
observer_altitude = isempty(altitude_args) ? nothing : parse(Float64, only(altitude_args))
run_case(ARGS[1], observer_altitude, dry_atmosphere, vza_zero, nstreams_30)
