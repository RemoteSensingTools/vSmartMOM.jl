#!/usr/bin/env julia

"""
Test a VLIDORT-style successive-Fourier-moment convergence criterion on the
production aerosol configuration without the cost or spectral complications
of a band calculation.

The experiment deliberately uses:

* one wavelength: the shortest wavelength in the O2 A band;
* the production three-aerosol mixture, 16-layer atmosphere, and nine streams;
* no gaseous absorption (`tau_abs` is explicitly zeroed after model build);
* the nadir/SZA=30 degree geometry of the nadir truth map; and
* the complete solver-selected Fourier series as the reference.

The production RT kernel is first run with its existing per-moment callback;
callbacks force the complete series, so convergence can be evaluated offline
without altering that reference. The ordinary noRS forward entry point is then
run with the workflow's configured early exit and checked against the complete
per-moment reconstruction. A moment passes when

    abs(delta_I_m) <= tolerance * abs(I_accumulated)

for every Stokes component and `n_consecutive` passing moments select the
early-exit point. The saved summary reports the actual I/Q/U error relative to
the complete Fourier sum.

Environment controls:

* `FOURIER_CONVERGENCE_TOLERANCES` (default `1e-2,1e-3,1e-4,1e-5,1e-6`)
* `FOURIER_CONVERGENCE_NCONSECUTIVE` (default `2`)
* `FOURIER_CONVERGENCE_MIN_M` (default and minimum `3`)
* `FOURIER_CONVERGENCE_OUT` (default below `truth_map_aerosols/`)
"""

include(joinpath(@__DIR__, "generate_truth_map.jl"))

using DelimitedFiles
using JLD2
using Printf

const CONVERGENCE_OUT = get(ENV, "FOURIER_CONVERGENCE_OUT", joinpath(
    ROOT, "truth_map_aerosols", "fourier_m_convergence_one_wavelength.jld2"))
const CONVERGENCE_TOLERANCES = FT.(parse.(Float64, split(get(
    ENV, "FOURIER_CONVERGENCE_TOLERANCES", "1e-2,1e-3,1e-4,1e-5,1e-6"), ',')))
const CONVERGENCE_NCONSECUTIVE = parse(Int, get(
    ENV, "FOURIER_CONVERGENCE_NCONSECUTIVE", "2"))
const CONVERGENCE_MIN_M = parse(Int, get(
    ENV, "FOURIER_CONVERGENCE_MIN_M", "3"))

CONVERGENCE_NCONSECUTIVE >= 1 || error(
    "FOURIER_CONVERGENCE_NCONSECUTIVE must be at least one")
CONVERGENCE_MIN_M >= 3 || error(
    "FOURIER_CONVERGENCE_MIN_M must be at least three")
all(>(zero(FT)), CONVERGENCE_TOLERANCES) || error(
    "all Fourier convergence tolerances must be positive")

"Build the production aerosol model at one wavelength and remove all gas opacity."
function monochromatic_absorption_free_model()
    state = read_states()[9] # urban, aerosol-on, SIF-off; CO2 is irrelevant here
    full_band = output_grids()[1]
    nu = FT[maximum(full_band)] # shortest wavelength in the O2 A band

    params = load_parameters(; float_type=FT, architecture=:CPU,
                             nstreams=AEROSOL_NSTREAMS)
    set_common!(params, state)
    params.architecture = CPU()
    params.vza = FT[VZA_DEG]
    params.vaz = FT[RELATIVE_AZIMUTH_DEG]

    surface_anchor = RRSXCO2Common.surface_basis_grids(FT)[1]
    surface = band_surface(state.coeff[1], surface_anchor, nu)
    select_band!(params, 1, nu, surface)

    # The current stream callback is restricted to the legacy embedded-solar
    # layout. Elastic SFI radiances are numerically equivalent to the default
    # external-solar layout; this test uses the callback only to expose every
    # individual m contribution without changing production plumbing.
    model = model_from_parameters(params; external_solar=false,
                                  aerosol_anchor_bands=[full_band])
    fill!(model.τ_abs[1], zero(FT))
    @assert all(iszero, Array(model.τ_abs[1]))
    return state, nu, model
end

"Recover the directional TOA Stokes contribution of every Fourier moment."
function directional_contributions(streams, model)
    nstokes = streams.pol_n
    nstokes >= 3 || error("the convergence experiment requires Stokes IQU")
    imu = CoreRT.nearest_point(streams.qp_μ, cosd(VZA_DEG))
    _, istart, iend = CoreRT.get_indices(imu, CoreRT.polarization_type(model))
    nmoments = length(streams.J⁻_per_m)
    contribution = zeros(FT, 3, nmoments)

    for im in 1:nmoments
        m = im - 1
        raw = @view streams.J⁻_per_m[im][istart:iend, 1, 1]
        cos_m_phi = cosd(m * RELATIVE_AZIMUTH_DEG)
        sin_m_phi = sind(m * RELATIVE_AZIMUTH_DEG)
        contribution[:, im] .= streams.weight[im] .* raw[1:3] .*
            FT[cos_m_phi, cos_m_phi, sin_m_phi]
    end
    return contribution
end

"Apply the candidate criterion to already-computed moment contributions."
function convergence_result(contribution, tolerance, n_consecutive)
    accumulated = zeros(FT, size(contribution, 1))
    passing = 0
    m_used = size(contribution, 2) - 1
    trace = zeros(FT, 8, size(contribution, 2))

    for im in axes(contribution, 2)
        delta = @view contribution[:, im]
        accumulated .+= delta
        denominators = abs.(accumulated)
        relative_increments = map(abs.(delta), denominators) do numerator, denominator
            denominator == 0 ? (iszero(numerator) ? zero(FT) : FT(Inf)) :
            numerator / denominator
        end
        relative_increment = maximum(relative_increments)
        pass = all(abs.(delta) .<= tolerance .* denominators)
        m = im - 1
        passing = m < CONVERGENCE_MIN_M ? 0 : (pass ? passing + 1 : 0)
        trace[:, im] .= (im - 1, delta[1], delta[2], delta[3],
                         accumulated[1], accumulated[2], accumulated[3],
                         relative_increment)
        if passing >= n_consecutive
            m_used = im - 1
            break
        end
    end
    return (; m_used, accumulated=copy(accumulated), trace)
end

function main_fourier_convergence_test()
    mkpath(dirname(CONVERGENCE_OUT))
    state, nu, model = monochromatic_absorption_free_model()
    @assert model.solver.m_max_bands == [15]
    @assert model.solver.n_fourier_moments_bands == [16]

    elapsed_streams = @elapsed streams = CoreRT.rt_run_streams(model)
    contribution = directional_contributions(streams, model)
    full_stokes = vec(sum(contribution; dims=2))

    # Independently pin the complete per-m reconstruction to the ordinary
    # production entry point, which now applies the configured early exit.
    elapsed_reference = @elapsed reference = CoreRT.rt_run(model)
    integrated_m_used = CoreRT._LAST_FOURIER_M_USED[]
    reference_toa, _ = reference
    reference_stokes = FT.(vec(Array(reference_toa)[1, 1:3, 1]))
    expected_integrated_m = min(CONVERGENCE_MIN_M - 1,
                                model.solver.m_max_bands[1]) +
                            CONVERGENCE_NCONSECUTIVE
    @assert integrated_m_used == expected_integrated_m
    @assert isapprox(full_stokes, reference_stokes;
                     rtol=FT(2e-5), atol=FT(2e-7))

    ntol = length(CONVERGENCE_TOLERANCES)
    summary = zeros(Float64, ntol, 11)
    results = Dict{String,Any}()
    for (it, tolerance) in enumerate(CONVERGENCE_TOLERANCES)
        result = convergence_result(
            contribution, tolerance, CONVERGENCE_NCONSECUTIVE)
        error = result.accumulated .- full_stokes
        relative_I_error = abs(error[1]) / max(abs(full_stokes[1]), eps(FT))
        moments_used = result.m_used + 1
        fraction_saved = 1 - moments_used / size(contribution, 2)
        summary[it, :] .= (tolerance, result.m_used, moments_used,
                           fraction_saved, relative_I_error,
                           error[2], error[3], result.accumulated[1],
                           result.accumulated[2], result.accumulated[3],
                           size(contribution, 2) - 1)
        results[@sprintf("%.1e", tolerance)] = result
    end

    jldsave(CONVERGENCE_OUT;
        summary, results, contribution, full_stokes, reference_stokes,
        tolerances=CONVERGENCE_TOLERANCES,
        n_consecutive=CONVERGENCE_NCONSECUTIVE,
        min_m=CONVERGENCE_MIN_M,
        m_max=model.solver.m_max_bands[1],
        nstreams=model.quad_points.Nstreams,
        nquad=model.quad_points.Nquad,
        tau_abs=Array(model.τ_abs[1]),
        wavelength_nm=1e7 / nu[1], wavenumber_cm_inverse=nu[1],
        sza_deg=SZA_DEG, vza_deg=VZA_DEG,
        relative_azimuth_deg=RELATIVE_AZIMUTH_DEG,
        state_index=state.index, elapsed_streams, elapsed_reference,
        integrated_m_used)

    dat_path = replace(CONVERGENCE_OUT, r"\.jld2$" => ".dat")
    open(dat_path, "w") do io
        println(io, "# One-wavelength, gas-free aerosol Fourier convergence test")
        println(io, "# wavelength_nm = ", 1e7 / nu[1])
        println(io, "# nstreams = ", model.quad_points.Nstreams,
                "; m_max = ", model.solver.m_max_bands[1],
                "; min_m = ", CONVERGENCE_MIN_M,
                "; n_consecutive = ", CONVERGENCE_NCONSECUTIVE)
        println(io, "# tolerance m_used moments_used fraction_saved rel_I_error abs_Q_error abs_U_error I_partial Q_partial U_partial m_max")
        writedlm(io, summary)
    end

    moments_path = replace(CONVERGENCE_OUT, r"\.jld2$" => "_moments.dat")
    cumulative = cumsum(contribution; dims=2)
    moment_table = hcat(collect(0:size(contribution, 2)-1),
                        permutedims(contribution), permutedims(cumulative))
    open(moments_path, "w") do io
        println(io, "# m delta_I delta_Q delta_U cumulative_I cumulative_Q cumulative_U")
        writedlm(io, moment_table)
    end

    println("Full Stokes = ", full_stokes)
    println("Per-m/full reconstruction difference = ",
            full_stokes .- reference_stokes)
    println("\n tolerance  m_used  saved    rel_I_error       delta_Q          delta_U")
    for row in eachrow(summary)
        @printf(" %9.1e  %6d  %6.1f%%  %13.5e  %13.5e  %13.5e\n",
                row[1], Int(row[2]), 100 * row[4], row[5], row[6], row[7])
    end
    println("\n", CONVERGENCE_OUT)
    println(dat_path)
    println(moments_path)
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main_fourier_convergence_test()
