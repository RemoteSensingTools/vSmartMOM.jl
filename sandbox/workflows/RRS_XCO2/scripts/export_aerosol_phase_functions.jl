#!/usr/bin/env julia

"""Export full-Mie phase functions for the three truth-map aerosols."""

include(joinpath(@__DIR__, "benchmark_aerosol_untruncated.jl"))
using Printf
using vSmartMOM.Scattering

const PHASE_OUT = joinpath(ROOT, "truth_map_aerosols", "aerosol_phase_functions_757nm.dat")
const AEROSOL_NAMES = ("sulfate", "organic_carbon", "utls_sulfate")

function main_phase()
    state = read_states()[9]
    ν = FT[maximum(output_grids()[1])]
    model = model_for_case(state, 1, ν, nothing, FT(0))
    optics = model.optics.aerosols.aerosol_optics[1]
    τ_aer = model.optics.aerosols.τ_aer[1]

    # Dense in scattering angle; reconstruct_phase takes cos(Theta).
    angle = collect(FT(0):FT(0.05):FT(180))
    μ = cosd.(angle)
    matrices = [Scattering.reconstruct_phase(op.greek_coefs, μ) for op in optics]
    rayleigh = Scattering.reconstruct_phase(
        model.optics.rayleigh.greek_rayleigh[1], μ)
    weights = [sum(Array(τ_aer[i, :, :])) *
               (op.ω̃ isa AbstractArray ? first(op.ω̃) : op.ω̃)
               for (i, op) in enumerate(optics)]
    weights ./= sum(weights)

    f11_mix = sum(weights[i] .* matrices[i].f₁₁ for i in eachindex(matrices))
    f12_mix = sum(weights[i] .* matrices[i].f₁₂ for i in eachindex(matrices))

    open(PHASE_OUT, "w") do io
        println(io, "# scattering_angle_deg species scattering_weight f11 f12 minus_f12_over_f11")
        for (name, weight, matrix) in zip(AEROSOL_NAMES, weights, matrices)
            for j in eachindex(angle)
                @printf(io, "%8.3f %-16s %.9e %.9e %.9e %.9e\n",
                        angle[j], name, weight, matrix.f₁₁[j], matrix.f₁₂[j],
                        -matrix.f₁₂[j] / matrix.f₁₁[j])
            end
        end
        for j in eachindex(angle)
            @printf(io, "%8.3f %-16s %.9e %.9e %.9e %.9e\n",
                    angle[j], "mixture", one(FT), f11_mix[j], f12_mix[j],
                    -f12_mix[j] / f11_mix[j])
        end
        for j in eachindex(angle)
            @printf(io, "%8.3f %-16s %.9e %.9e %.9e %.9e\n",
                    angle[j], "rayleigh", zero(FT), rayleigh.f₁₁[j],
                    rayleigh.f₁₂[j], -rayleigh.f₁₂[j] / rayleigh.f₁₁[j])
        end
    end
    println(PHASE_OUT)
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main_phase()
