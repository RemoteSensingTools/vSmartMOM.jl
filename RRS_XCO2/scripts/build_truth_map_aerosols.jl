#!/usr/bin/env julia

"""Create reproducible altitude-lognormal aerosol profile tables for the truth map."""

using DelimitedFiles
using Printf

const ROOT = normpath(joinpath(@__DIR__, ".."))
const INPUT = joinpath(ROOT, "truth_map", "aerosol_vertical_grids.dat")
const OUTDIR = joinpath(ROOT, "truth_map_aerosols")
const OUTPUT = joinpath(OUTDIR, "aerosol_vertical_profiles.dat")

# z_mode is the height of maximum extinction. The LogNormal constructor takes
# log(median), not log(mode), so median = mode * exp(sigma_ln^2).
const AEROSOLS = (
    ("sulfate",        0.1935471100, 1.2 * exp(0.49^2), 0.49, 1.2),
    ("organic_carbon", 0.0807084200, 1.8 * exp(0.40^2), 0.40, 1.8),
    ("utls_sulfate",   0.0057444777, 12.0 * exp(0.10^2), 0.10, 12.0),
)

function read_geometry()
    rows = readdlm(INPUT, comments=true)
    return Float64.(rows[:, 1:6])
end

normal_cdf(x) = 0.5 * (1 + ccall(:erf, Cdouble, (Cdouble,), x / sqrt(2.0)))

function lognormal_cdf(z, median, σln)
    z <= 0 && return 0.0
    return normal_cdf(log(z / median) / σln)
end

function main()
    mkpath(OUTDIR)
    geometry = read_geometry()
    model_top = maximum(geometry[:, 3])
    open(OUTPUT, "w") do io
        println(io, "# Altitude-lognormal truth-map aerosol profiles")
        println(io, "# Layer AOD is computed by exact CDF differences at the model boundaries.")
        println(io, "# nlayer layer z_top_km z_bottom_km z_center_km dz_km species tau760 z_median_km sigma_ln_z z_mode_km layer_fraction layer_aod760")
        for row in eachrow(geometry)
            nlayer, layer = Int.(round.(row[1:2]))
            ztop, zbottom, zcenter, dz = row[3:6]
            for (name, τ760, zmedian, σln, zmode) in AEROSOLS
                column_norm = lognormal_cdf(model_top, zmedian, σln)
                fraction = (lognormal_cdf(ztop, zmedian, σln) -
                            lognormal_cdf(zbottom, zmedian, σln)) / column_norm
                @printf(io, "%d %d %.9f %.9f %.9f %.9f %-16s %.10f %.9f %.6f %.9f %.12e %.12e\n",
                        nlayer, layer, ztop, zbottom, zcenter, dz, name,
                        τ760, zmedian, σln, zmode, fraction, τ760 * fraction)
            end
        end
    end
    println(OUTPUT)
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
