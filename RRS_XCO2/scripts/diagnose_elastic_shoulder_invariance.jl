#!/usr/bin/env julia

"""
Diagnose why an elastic O2 core changes when carried inside the Raman solve
shoulders. This script does not modify truth products.

It builds the same centered 256-point core twice: once alone and once inside
the production ±234 cm⁻¹ solve grid. Aerosol and surface interpolation remain
anchored to the canonical full A-band in both models. Precomputed, effective
per-Fourier optical properties and final noRS Stokes spectra are compared.
"""

ENV["RRS_XCO2_FLOAT_TYPE"] = "Float32"
ENV["NLAYERS"] = "16"
ENV["AEROSOL_NSTREAMS"] = "9"
ENV["TRUTH_SZA_DEG"] = "30"
ENV["TRUTH_VZA_DEG"] = "0"
ENV["TRUTH_RELATIVE_AZIMUTH_DEG"] = "0"

using CUDA
using LinearAlgebra
using Printf

include(joinpath(@__DIR__, "generate_truth_map.jl"))

function metrics(reference, candidate)
    a, b = Float64.(reference), Float64.(candidate)
    size(a) == size(b) || error("shape mismatch $(size(a)) versus $(size(b))")
    d = b .- a
    return (maxabs=maximum(abs, d),
            rel_l2=norm(d) / max(norm(a), eps(Float64)))
end

function core_slice(array, keep)
    host = Array(array)
    ndims(host) == 1 && return host[keep]
    return host[ntuple(_ -> Colon(), ndims(host) - 1)..., keep]
end

function print_metric(label, metric)
    @printf("%-38s maxabs=%12.5e relL2=%12.5e\n",
            label, metric.maxabs, metric.rel_l2)
end

function compare_greek(label, a, b)
    for name in (:α, :β, :γ, :δ, :ϵ, :ζ)
        print_metric("$label.$name", metrics(getfield(a, name), getfield(b, name)))
    end
end

function main()
    CUDA.functional() || error("CUDA is not functional")
    CUDA.device!(parse(Int, get(ENV, "CUDA_DEVICE", "1")))

    state = read_states()[9]
    full = output_grids()[1]
    width = min(256, length(full))
    first_index = fld(length(full) - width, 2) + 1
    full_range = first_index:first_index + width - 1
    core = full[full_range]
    solve, keep = o2_solve_grid(core)

    core_model = build_o2(state, full, core)
    shoulder_model = build_o2(state, full, solve)
    println("core=$(length(core)) solve=$(length(solve)) keep=$(length(keep))")
    println("m_max core=$(core_model.solver.m_max_bands) " *
            "shoulder=$(shoulder_model.solver.m_max_bands)")
    println("l_max core=$(core_model.solver.l_max) " *
            "shoulder=$(shoulder_model.solver.l_max)")

    print_metric("wavenumber", metrics(core, solve[keep]))
    print_metric("surface albedo", metrics(
        core_model.surfaces[1].albedo,
        shoulder_model.surfaces[1].albedo[keep]))
    print_metric("tau_abs", metrics(
        core_model.τ_abs[1], shoulder_model.τ_abs[1][keep, :]))
    print_metric("tau_rayleigh", metrics(
        core_model.τ_rayl[1], shoulder_model.τ_rayl[1][keep, :]))
    print_metric("tau_aerosol", metrics(
        core_model.τ_aer[1], shoulder_model.τ_aer[1][:, keep, :]))
    print_metric("Cabannes fraction", metrics(
        core_model.ϖ_Cabannes, shoulder_model.ϖ_Cabannes))
    compare_greek("rayleigh greek", core_model.greek_rayleigh[1],
                  shoulder_model.greek_rayleigh[1])

    for iaerosol in 1:3
        a = core_model.aerosol_optics[1][iaerosol]
        b = shoulder_model.aerosol_optics[1][iaerosol]
        print_metric("aerosol $iaerosol k", metrics(a.k, b.k[keep]))
        print_metric("aerosol $iaerosol omega", metrics(a.ω̃, b.ω̃[keep]))
        print_metric("aerosol $iaerosol ftrunc", metrics(a.fᵗ, b.fᵗ[keep]))
        compare_greek("aerosol $iaerosol node1", a.phase_greek[1], b.phase_greek[1])
        compare_greek("aerosol $iaerosol node2", a.phase_greek[end], b.phase_greek[end])
    end

    solar_T = solar_interpolator()
    core_F0, core_sources = source_set(core, false, solar_T)
    shoulder_F0, shoulder_sources = source_set(solve, false, solar_T)
    print_metric("solar source", metrics(core_F0, shoulder_F0[:, keep]))

    core_rs = make_nors(core_F0)
    shoulder_rs = make_nors(shoulder_F0)
    core_cache = CoreRT.build_m_invariant_cache(core_rs, [1], core_model)
    shoulder_cache = CoreRT.build_m_invariant_cache(
        shoulder_rs, [1], shoulder_model)
    worst = Dict(name => (maxabs=0.0, rel_l2=0.0)
                 for name in (:tau, :omega, :Zpp, :Zmp, :Z0p, :Z0m))
    ndoubl_core = Int[]
    ndoubl_shoulder = Int[]
    # m=0 is sufficient to locate whether the divergence enters the effective
    # optical-property mixer. The endpoint Greek coefficients were compared
    # above and determine every higher Fourier order; downloading all large
    # shoulder Z arrays for every m would dominate this diagnostic.
    for m in (0,)
        core_layers, _ = CoreRT.constructCoreOpticalProperties(
            core_rs, [1], m, core_model, core_cache)
        shoulder_layers, _ = CoreRT.constructCoreOpticalProperties(
            shoulder_rs, [1], m, shoulder_model, shoulder_cache)
        for iz in eachindex(core_layers)
            a, b = core_layers[iz], shoulder_layers[iz]
            comparisons = (
                tau=metrics(Array(a.τ), core_slice(b.τ, keep)),
                omega=metrics(Array(a.ϖ), core_slice(b.ϖ, keep)),
                Zpp=metrics(Array(a.Z⁺⁺), core_slice(b.Z⁺⁺, keep)),
                Zmp=metrics(Array(a.Z⁻⁺), core_slice(b.Z⁻⁺, keep)),
                Z0p=metrics(Array(a.Z₀⁺), core_slice(b.Z₀⁺, keep)),
                Z0m=metrics(Array(a.Z₀⁻), core_slice(b.Z₀⁻, keep)),
            )
            for (name, metric) in pairs(comparisons)
                old = worst[name]
                worst[name] = (maxabs=max(old.maxabs, metric.maxabs),
                               rel_l2=max(old.rel_l2, metric.rel_l2))
            end
            if m == 0
                push!(ndoubl_core,
                    CoreRT.get_dtau_ndoubl(a, core_model.quad_points)[2])
                push!(ndoubl_shoulder,
                    CoreRT.get_dtau_ndoubl(b, shoulder_model.quad_points)[2])
            end
        end
    end
    for name in (:tau, :omega, :Zpp, :Zmp, :Z0p, :Z0m)
        print_metric("effective $name (worst m=0,z)", worst[name])
    end
    println("ndoubl core=$ndoubl_core")
    println("ndoubl shoulder=$ndoubl_shoulder")

    core_toa = toa3(CoreRT.rt_run_toa(
        core_model; i_band=1, sources=core_sources))
    shoulder_toa = toa3(CoreRT.rt_run_toa(
        shoulder_model; i_band=1, sources=shoulder_sources), keep)
    for istokes in 1:3
        print_metric("TOA Stokes $istokes", metrics(
            core_toa[istokes, :], shoulder_toa[istokes, :]))
    end
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
