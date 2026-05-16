#!/usr/bin/env julia

# Development driver for fast AOD-only Mie extinction lookup tables.
# Run from the repository root with:
#
#   julia --project=. src/Aerosols/LUT/generate_mie_extinction_lut.jl out.jld2
#
# Environment controls:
#   VSMARTMOM_AER_EXT_LUT_X        start:stop:count or log:start:stop:count
#   VSMARTMOM_AER_EXT_LUT_N_REAL   start:stop:count
#   VSMARTMOM_AER_EXT_LUT_N_IMAG   start:stop:count
#   VSMARTMOM_AER_EXT_LUT_FLOAT    Float32 or Float64

include("mie_extinction_lut.jl")

function _range_env(name, default; allow_log = false)
    raw = get(ENV, name, default)
    parts = split(raw, ":")
    if allow_log && length(parts) == 4 && parts[1] == "log"
        lo = parse(Float64, parts[2])
        hi = parse(Float64, parts[3])
        n = parse(Int, parts[4])
        lo > 0 || throw(ArgumentError("$name logarithmic start must be positive"))
        n >= 2 || throw(ArgumentError("$name count must be at least 2"))
        return exp.(range(log(lo), log(hi); length = n))
    elseif length(parts) == 3
        lo = parse(Float64, parts[1])
        hi = parse(Float64, parts[2])
        n = parse(Int, parts[3])
        n >= 2 || throw(ArgumentError("$name count must be at least 2"))
        return collect(range(lo, hi; length = n))
    end
    throw(ArgumentError("$name must be start:stop:count" *
                        (allow_log ? " or log:start:stop:count" : "")))
end

function main(args = ARGS)
    out = isempty(args) ? "mie_extinction_lut.jld2" : args[1]

    x_spec = get(ENV, "VSMARTMOM_AER_EXT_LUT_X", "log:0.001:140.0:900")
    nr_spec = get(ENV, "VSMARTMOM_AER_EXT_LUT_N_REAL", "1.20:2.20:41")
    ni_spec = get(ENV, "VSMARTMOM_AER_EXT_LUT_N_IMAG", "0.0:1.0:41")
    x_grid = _range_env("VSMARTMOM_AER_EXT_LUT_X", x_spec; allow_log = true)
    n_real_grid = _range_env("VSMARTMOM_AER_EXT_LUT_N_REAL", nr_spec)
    n_imag_grid = _range_env("VSMARTMOM_AER_EXT_LUT_N_IMAG", ni_spec)
    storage = get(ENV, "VSMARTMOM_AER_EXT_LUT_FLOAT", "Float32") == "Float64" ? Float64 : Float32

    metadata = Dict{String, Any}(
        "x_grid_env" => x_spec,
        "n_real_grid_env" => nr_spec,
        "n_imag_grid_env" => ni_spec,
        "float_type" => string(storage),
    )

    @info "Building aerosol Mie extinction LUT" out storage nx=length(x_grid) nr=length(n_real_grid) ni=length(n_imag_grid)
    lut = AerosolMieExtinctionLUT.build_mie_extinction_lut(
        x_grid, n_real_grid, n_imag_grid; FT = storage, metadata = metadata)
    AerosolMieExtinctionLUT.save_mie_extinction_lut(out, lut)
    @info "Wrote aerosol Mie extinction LUT" out
    return out
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
