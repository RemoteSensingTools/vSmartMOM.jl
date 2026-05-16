#!/usr/bin/env julia

# Development driver for homogeneous-sphere aerosol Mie/Greek lookup tables.
# Run from the repository root with:
#
#   julia --project=. src/Aerosols/LUT/generate_mie_greek_lut.jl out.jld2
#
# The grids below are intentionally controlled by environment variables so the
# first committed version stays a reusable development harness, not a frozen
# production LUT recipe.

include("mie_greek_lut.jl")

function _range_env(name, default)
    raw = get(ENV, name, default)
    parts = split(raw, ":")
    length(parts) == 3 ||
        throw(ArgumentError("$name must be formatted as start:stop:count"))
    lo = parse(Float64, parts[1])
    hi = parse(Float64, parts[2])
    n = parse(Int, parts[3])
    n >= 2 || throw(ArgumentError("$name count must be at least 2"))
    return collect(range(lo, hi; length = n))
end

function main(args = ARGS)
    out = isempty(args) ? "mie_greek_lut.jld2" : args[1]

    x_grid = _range_env("VSMARTMOM_AER_LUT_X", "0.001:80.0:160")
    n_real_grid = _range_env("VSMARTMOM_AER_LUT_N_REAL", "1.25:1.95:36")
    n_imag_grid = _range_env("VSMARTMOM_AER_LUT_N_IMAG", "0.0:0.8:41")
    lmax = parse(Int, get(ENV, "VSMARTMOM_AER_LUT_LMAX", "180"))
    storage = get(ENV, "VSMARTMOM_AER_LUT_FLOAT", "Float32") == "Float64" ? Float64 : Float32

    metadata = Dict{String, Any}(
        "x_grid_env" => get(ENV, "VSMARTMOM_AER_LUT_X", "0.001:80.0:160"),
        "n_real_grid_env" => get(ENV, "VSMARTMOM_AER_LUT_N_REAL", "1.25:1.95:36"),
        "n_imag_grid_env" => get(ENV, "VSMARTMOM_AER_LUT_N_IMAG", "0.0:0.8:41"),
        "lmax_env" => string(lmax),
        "float_type" => string(storage),
    )

    @info "Building aerosol Mie/Greek LUT" out lmax storage nx=length(x_grid) nr=length(n_real_grid) ni=length(n_imag_grid)
    lut = AerosolMieGreekLUT.build_mie_greek_lut(x_grid, n_real_grid, n_imag_grid;
                                                 lmax = lmax, FT = storage,
                                                 metadata = metadata)
    AerosolMieGreekLUT.save_mie_greek_lut(out, lut)
    @info "Wrote aerosol Mie/Greek LUT" out
    return out
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
