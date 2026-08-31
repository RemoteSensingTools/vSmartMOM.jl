#!/usr/bin/env julia

"""Export the saved testing-only truncated m=0 JLD2 result as ASCII."""

using JLD2
using DelimitedFiles

const ROOT = normpath(joinpath(@__DIR__, "..", "truth_map_aerosols"))
const INPUT = joinpath(ROOT, "truncated_m0_stokes.jld2")
const OUTPUT = joinpath(ROOT, "truncated_m0_stokes.dat")

function main()
    d = load(INPUT)
    caps = Int.(d["caps"])
    vza = d["vza"]
    vaz = d["vaz"]
    m0 = d["m0_values"]
    rows = Matrix{Float64}(undef, length(caps) * length(vza) * length(vaz), 6)
    row = 1
    for cap in caps, iv in eachindex(vza), iaz in eachindex(vaz)
        rows[row, :] .= (cap, vza[iv], vaz[iaz], m0[cap][iv, iaz, 1],
                         m0[cap][iv, iaz, 2], m0[cap][iv, iaz, 3])
        row += 1
    end
    open(OUTPUT, "w") do io
        println(io, "# l_trunc vza relaz I_m0 Q_m0 U_m0")
        writedlm(io, rows)
    end
    println(OUTPUT)
end

main()
