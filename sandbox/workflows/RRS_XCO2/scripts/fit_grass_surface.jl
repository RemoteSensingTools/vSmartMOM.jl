#!/usr/bin/env julia

# Fit P0, P1, P2 independently in each OCO band. These representative green-
# vegetation anchors describe a bright NIR plateau followed by decreasing SWIR
# reflectance. Replace TARGETS with a measured spectrum when one is available.
using DelimitedFiles
using LinearAlgebra

const ROOT = normpath(joinpath(@__DIR__, ".."))
const BAND_NAMES = ("o2a", "weak_co2", "strong_co2")
const WAVELENGTHS_NM = ([757.0, 765.0, 773.0],
                        [1589.0, 1605.5, 1622.0],
                        [2042.0, 2063.0, 2084.0])
const TARGETS = ([0.48, 0.46, 0.43],
                 [0.34, 0.36, 0.38],
                 [0.17, 0.20, 0.24])

# vSmartMOM evaluates Legendre surfaces on ascending wavenumber index. Thus
# x=-1 corresponds to the longest wavelength and x=+1 to the shortest.
function fit_band(λ_nm, reflectance)
    order = sortperm(1.0e7 ./ λ_nm)
    y = reflectance[order]
    x = collect(range(-1.0, 1.0, length=length(y)))
    basis = hcat(ones(length(x)), x, (3 .* x.^2 .- 1) ./ 2)
    return basis \ y
end

mkpath(joinpath(ROOT, "output"))
rows = Matrix{Any}(undef, length(BAND_NAMES) + 1, 4)
rows[1, :] = ["band", "P0", "P1", "P2"]
for (i, name) in enumerate(BAND_NAMES)
    coeff = fit_band(WAVELENGTHS_NM[i], TARGETS[i])
    rows[i + 1, :] = [name, coeff...]
    println(name, ": LambertianSurfaceLegendre(", collect(coeff), ")")
end
writedlm(joinpath(ROOT, "output", "grass_legendre_coefficients.csv"), rows, ',')

