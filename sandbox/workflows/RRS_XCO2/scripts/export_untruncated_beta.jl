#!/usr/bin/env julia

"""Export the full untruncated aerosol Greek beta_l series at 757 nm."""

include(joinpath(@__DIR__, "benchmark_aerosol_untruncated.jl"))
using Printf

const BETA_OUT = joinpath(ROOT, "truth_map_aerosols", "untruncated_beta_757nm.dat")
const BETA_NAMES = ("sulfate", "organic_carbon", "utls_sulfate")

function beta_vector(op)
    β = op.greek_coefs.β
    ndims(β) == 1 ? Array(β) : vec(Array(@view β[:, 1]))
end

function main_beta()
    state = read_states()[9]
    ν = FT[maximum(output_grids()[1])]
    model = model_for_case(state, 1, ν, nothing, FT(0))
    optics = model.optics.aerosols.aerosol_optics[1]
    τ_aer = model.optics.aerosols.τ_aer[1]
    β = [beta_vector(op) for op in optics]

    # Phase-matrix mixing uses scattering optical depth τ·ω, not extinction
    # optical depth alone. This is the effective aerosol-only beta series.
    weights = [sum(Array(τ_aer[i, :, :])) *
               (op.ω̃ isa AbstractArray ? first(Array(op.ω̃)) : op.ω̃)
               for (i, op) in enumerate(optics)]
    weights ./= sum(weights)
    βmix = sum(weights[i] .* β[i] for i in eachindex(β))

    open(BETA_OUT, "w") do io
        println(io, "# l species scattering_weight beta_l")
        for (name, weight, series) in zip(BETA_NAMES, weights, β)
            for il in eachindex(series)
                @printf(io, "%d %-16s %.9e %.9e\n", il - 1, name,
                        weight, series[il])
            end
        end
        for il in eachindex(βmix)
            @printf(io, "%d %-16s %.9e %.9e\n", il - 1, "mixture",
                    one(FT), βmix[il])
        end
    end
    println(BETA_OUT)
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main_beta()
