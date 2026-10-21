#!/usr/bin/env julia

"""Record the untruncated TOA source contribution at every stream and Fourier order."""

include(joinpath(@__DIR__, "benchmark_aerosol_untruncated.jl"))
using JLD2
using DelimitedFiles

const FOURIER_OUT = get(ENV, "FOURIER_OUT", joinpath(
    ROOT, "truth_map_aerosols", "untruncated_fourier_stream_contributions.jld2"))

function main_fourier_streams()
    state = read_states()[9]
    ν = sparse_grid(1)
    # The callback API currently requires the legacy embedded-solar layout.
    # Its TOA solution is numerically equivalent to external-solar SFI; using
    # it here exposes the settled per-m source columns without changing RT.
    model = model_for_case(state, 1, ν, nothing, FT(0); external_solar=false)
    _, sources = source_set(ν, false, solar_interpolator())

    m_values = Int[]
    weights = FT[]
    raw = Matrix{FT}[]
    qp_μ = FT[]

    callback = function (step)
        isempty(qp_μ) && append!(qp_μ, FT.(Array(step.qp_μ)))
        J = Array(step.composite_layer.J₀⁻)
        for slot in values(step.composite_layer.J₀_by_src)
            J .+= Array(slot.J₀⁻)
        end
        # Stokes components are contiguous within each angular stream.
        push!(raw, copy(reshape(@view(J[:, 1, 1]), step.pol_type.n, :)[1:3, :]))
        push!(m_values, step.m)
        push!(weights, FT(step.weight))
        return nothing
    end

    FTtype = CoreRT.float_type(model)
    elapsed = @elapsed CoreRT._rt_run_column(
        vSmartMOM.InelasticScattering.noRS{FTtype}(), model, [1];
        sources, streams_callback=callback, toa_only=true)

    # Layout: (m, stream, Stokes).  `weighted` is the coefficient entering
    # azimuthal reconstruction.  I/Q multiply cos(mΔφ); U multiplies
    # sin(mΔφ), so U is a sine amplitude rather than its value at Δφ=0.
    raw_coefficients = permutedims(cat(raw...; dims=3), (3, 2, 1))
    weighted_contributions = raw_coefficients .* reshape(weights, :, 1, 1)
    vza_quadrature = acosd.(clamp.(qp_μ, -one(FT), one(FT)))

    jldsave(FOURIER_OUT; m=m_values, weights, qp_μ, vza_quadrature,
            raw_coefficients, weighted_contributions, elapsed,
            nstreams=model.quad_points.Nstreams,
            nquad=model.quad_points.Nquad,
            wavelength_nm=1e7 / ν[1], wavenumber_cm⁻¹=ν[1],
            external_solar=model.quad_points.external_solar,
            state_index=state.index)
    table = Matrix{Float64}(undef, length(m_values) * length(qp_μ), 6)
    row = 1
    for im in eachindex(m_values), iq in eachindex(qp_μ)
        table[row, :] .= (m_values[im], qp_μ[iq],
                          weighted_contributions[im, iq, 1],
                          weighted_contributions[im, iq, 2],
                          weighted_contributions[im, iq, 3], elapsed)
        row += 1
    end
    open(replace(FOURIER_OUT, r"\.jld2$" => ".dat"), "w") do io
        println(io, "# m mu weighted_I weighted_Q weighted_U elapsed_s")
        writedlm(io, table)
    end
    println(FOURIER_OUT)
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main_fourier_streams()
