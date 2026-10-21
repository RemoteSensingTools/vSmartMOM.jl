#!/usr/bin/env julia

using DelimitedFiles
using Printf
using vSmartMOM

include(joinpath(@__DIR__, "common.jl"))
using .RRSXCO2Common

const ROOT = normpath(joinpath(@__DIR__, ".."))
const OUT = get(ENV, "TRUTH_OUT", joinpath(ROOT, "truth_map"))
const SURFACES = ("urban", "rural", "desert", "forest")
const XCO2_PPM = (380, 400, 420, 440)
const AEROSOL_CASES = ("none", "aod760_0p28")
const SIF_CASES = ("off", RRSXCO2Common.SIF_CASE_ON)
const AOD550 = Dict("none" => (0.0, 0.0, 0.0),
                    # Mie-normalized so total AOD(760 nm)=0.28 while the
                    # species AOD ratio at 550 nm remains 8:1.8:0.2.
                    "aod760_0p28" => (0.366505284486836,
                                      0.082463689009538,
                                      0.009162632112170))

function load_surface_coefficients()
    raw = readdlm(joinpath(ROOT, "surface_albedos", "lambertian_legendre_inputs.dat"),
                  comments=true)
    result = Dict{Tuple{String,String},NTuple{3,Float64}}()
    for row in eachrow(raw)
        result[(String(row[1]), String(row[2]))] =
            (Float64(row[3]), Float64(row[4]), Float64(row[5]))
    end
    return result
end

function main()
    mkpath(OUT)
    coeff = load_surface_coefficients()
    sif = RRSXCO2Common.campaign_sif_state()
    path = joinpath(OUT, "true_states.dat")
    temporary = path * ".tmp.$(getpid())"
    open(temporary, "w") do io
        println(io, "# Exact truth-map state table; state order: surface -> aerosol -> SIF -> XCO2")
        println(io, "# AOD550 ratio sulfate:organic_carbon:stratospheric_sulfate = 8:1.8:0.2; total AOD760 = 0.28")
        println(io, "# SIF definition v$(RRSXCO2Common.SIF_DEFINITION_VERSION): 2pi*L_lambda(760 nm)=0.5 mW m^-2 nm^-1; every isotropic BOA stream has L_lambda=0.5/(2pi) per sr.")
        println(io, "# sif_angular_integral760 is the stated 2pi*L_lambda value; SIF760 and mSIF are the matching per-wavenumber radiance coordinates.")
        println(io, "# index surface_index surface aerosol_index aerosol_case sif_index sif_case xco2_index xco2_ppm psurf_hpa sza_deg vza_deg sulfate_aod550 organic_aod550 stratospheric_aod550 sif_angular_integral760 SIF760 mSIF o2a_P0 o2a_P1 o2a_P2 weak_P0 weak_P1 weak_P2 strong_P0 strong_P1 strong_P2")
        index = 0
        for (isurf, surface) in enumerate(SURFACES),
            (iaer, aerosol) in enumerate(AEROSOL_CASES),
            (isif, sif_case) in enumerate(SIF_CASES),
            (ixco2, xco2) in enumerate(XCO2_PPM)
            index += 1
            aod = AOD550[aerosol]
            sif_on = sif_case == RRSXCO2Common.SIF_CASE_ON
            sc = sif_on ? (
                RRSXCO2Common.SIF_ANGULAR_INTEGRAL_760,
                sif.SIF760,
                sif.mSIF,
            ) : (0.0, 0.0, 0.0)
            c1, c2, c3 = coeff[(surface, "o2a")], coeff[(surface, "weak_co2")], coeff[(surface, "strong_co2")]
            @printf(io, "%03d %d %s %d %s %d %s %d %d 1000.0 30.0 0.0 %.8f %.8f %.8f %.8e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e\n",
                    index, isurf, surface, iaer, aerosol, isif, sif_case, ixco2,
                    xco2, aod..., sc..., c1..., c2..., c3...)
        end
        @assert index == 64
    end
    mv(temporary, path; force=true)
    println(path)
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
