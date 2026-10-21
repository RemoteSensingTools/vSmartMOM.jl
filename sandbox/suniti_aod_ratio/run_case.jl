using CUDA
using DelimitedFiles
using Statistics
using vSmartMOM

const SOURCE_CONFIG = "/home/cfranken/.codex/attachments/863294a0-bb39-460e-950d-063f515085bf/pasted-text.txt"
const OUTPUT_DIR = @__DIR__

function load_student_config()
    source = read(SOURCE_CONFIG, String)
    expression, _ = Meta.parse(source, 1; greedy=true, raise=true)
    expression.head == :(=) || error("Expected the first expression to assign cfg")
    return Core.eval(Main, expression)
end

function prepare_config(kind::String, observer_altitude)
    cfg = deepcopy(load_student_config())
    rt = cfg["radiative_transfer"]
    aerosols = cfg["scattering"]["aerosols"]

    # Build at AOD 0.03, then zero/restore the model optical-depth array so
    # both RT calls use precisely the same gas, Mie, and quadrature state.
    aerosols[1]["τ_ref"] = 0.03
    if observer_altitude !== nothing
        cfg["geometry"]["obs_alt"] = [observer_altitude]
    end

    if kind == "fixed"
        # Documented v0.7 schema. The original duplicate Dict key has already
        # resolved to NoTruncation() by the time Julia constructs the Dict.
        rt["truncation"] = "auto"
        rt["Δ_angle"] = 2.0
        delete!(rt, "max_m")
        delete!(rt, "l_trunc")

        # The other six entries have exactly zero loading. Removing them is
        # numerically neutral and avoids six unused Mie calculations.
        cfg["scattering"]["aerosols"] = [aerosols[1]]
    elseif kind != "original"
        error("case must be original or fixed")
    end
    return cfg
end

function requested_spectra(result, observer_altitude)
    if observer_altitude === nothing
        up = vec(Array(result.toa[1, 1, :]))
        down = vec(Array(result.boa[1, 1, :]))
        return up, down
    end
    level = only(result.levels)
    up = vec(Array(level.upwelling[1, 1, :]))
    down = vec(Array(total_downwelling(level)[1, 1, :]))
    return up, down
end

function finite_extrema(x)
    values = filter(isfinite, x)
    return extrema(values)
end

function run_case(kind::String, observer_altitude)
    CUDA.allowscalar(false)
    cfg = prepare_config(kind, observer_altitude)

    println("CASE=$kind")
    println("CUDA device: ", CUDA.name(CUDA.device()))
    println("radiative_transfer: ", cfg["radiative_transfer"])
    println("aerosol count: ", length(cfg["scattering"]["aerosols"]))

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
    output = joinpath(OUTPUT_DIR, "$(kind)$(suffix)_spectra.csv")
    open(output, "w") do io
        println(io, "wavenumber_cm-1,up_aod0,up_aod003,down_aod0,down_aod003")
        writedlm(io, table, ',')
    end

    up_ratio = up03 ./ up0
    down_ratio = down03 ./ down0
    println("upwelling ratio extrema: ", finite_extrema(up_ratio))
    println("downwelling ratio extrema: ", finite_extrema(down_ratio))
    println("upwelling ratio mean: ", mean(filter(isfinite, up_ratio)))
    println("downwelling ratio mean: ", mean(filter(isfinite, down_ratio)))
    println("wrote ", output)
end

1 <= length(ARGS) <= 2 ||
    error("Usage: julia run_case.jl original|fixed [observer_altitude_km]")
observer_altitude = length(ARGS) == 2 ? parse(Float64, ARGS[2]) : nothing
run_case(ARGS[1], observer_altitude)
