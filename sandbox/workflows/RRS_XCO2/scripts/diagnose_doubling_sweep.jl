#!/usr/bin/env julia

"""Locate the Float32 failure stage by sweeping the elemental-layer thickness.

`get_dtau_ndoubl` picks `dτ_max = max(dτ_min_floor, min(τϖ, threshold·μ_min))`
and then doubles `ceil(log2(τϖ/dτ_max))` times. Two OPPOSING error sources act
on that choice in Float32:

- `dτ_max` too SMALL → the elemental r/t entries scale as `dτ/μ`, so the seed
  layer carries few significant Float32 digits, and that noise compounds over
  every doubling.
- `dτ_max` too LARGE → `dτ/μ_min` approaches 1 and the exact finite-δ elemental
  seed is no longer a thin layer for the near-horizon streams.

The default `threshold = 0.001` scales `dτ_max` with `μ_min`, and `μ_min`
shrinks like 1/N². Raising `nstreams` therefore drives the solver into the
FIRST regime. This sweep varies `threshold` at fixed physics (δBGE(32)) and
reports I against the Float64 answer, isolating which regime is responsible.
"""

# `diagnose_untruncated_precision.jl` already includes
# `benchmark_aerosol_untruncated.jl` (and through it `generate_truth_map.jl`
# and `common.jl`). Including both would define `RRSXCO2Common` twice and make
# `load_parameters` ambiguous, so include only the outer file and qualify the
# helper explicitly.
include(joinpath(@__DIR__, "diagnose_untruncated_precision.jl"))
using Printf

const SWEEP_OUT = joinpath(ROOT, "truth_map_aerosols", "doubling_sweep.dat")
const SWEEP_VZA = 20.0

# (nstreams, threshold, floor value or nothing, float type)
const SWEEP = [
    (111, 1e-3, nothing,   Float64),   # converged reference
    (111, 1e-3, nothing,   Float32),   # as-run: floor inactive (2.3e-13)
    (111, 1e-3, 1.2207e-4, Float32),   # the INTENDED Float32 floor, 1024·eps
    (111, 1e-2, nothing,   Float32),   # thicker elemental, ~3 fewer doublings
    (111, 1e-1, nothing,   Float32),   # thicker still
    (111, 1e-4, nothing,   Float32),   # thinner, ~3 more doublings
    ( 50, 1e-3, nothing,   Float64),
    ( 50, 1e-3, nothing,   Float32),
    ( 50, 1e-2, nothing,   Float32),
    ( 17, 1e-3, nothing,   Float64),
    ( 17, 1e-3, nothing,   Float32),
]

function sweep_model(state, ν, vza; nstr, thr, floorval, W)
    p = RRSXCO2Common.load_parameters()
    set_common_prec!(p, state, W)
    p.numerics = CoreRT.RTNumericalParameters{W}(;
        dτ_max_threshold = W(thr),
        dτ_min_floor = floorval === nothing ? W(0) : W(floorval),
        blas_threads = p.numerics.blas_threads,
        verbose = p.numerics.verbose)
    p.vza = W[vza, vza]; p.vaz = W[0, 180]
    p.l_trunc = 2 * nstr - 2
    p.nstreams = nstr
    p.stream_l_cap = 2 * nstr - 1
    p.legacy_l_cap_override = nothing
    p.max_m = 33
    p.truncation = vSmartMOM.Scattering.δBGE{W}(32, zero(W))
    p.spec_bands = [W.(ν)]
    p.brdf = [CoreRT.LambertianSurfaceScalar(W(sum(state.coeff[1])))]
    ap = p.absorption_params
    ap.fixed_molecules = [ap.fixed_molecules[1]]
    ap.variable_molecules = [ap.variable_molecules[1]]
    ap.luts = [ap.luts[1]]
    if length(ap.cia_files) > 1
        ap.cia_files = [ap.cia_files[1]]
        ap.cia_reference_codes = [ap.cia_reference_codes[1]]
        ap.cia_negative_policies = [ap.cia_negative_policies[1]]
    end
    return model_from_parameters(p; external_solar=true)
end

function main_sweep()
    state = read_states()[9]
    ν = sparse_grid(1)
    solar_T = solar_interpolator()
    rows = []
    for (nstr, thr, floorval, W) in SWEEP
        @info "doubling sweep" nstr thr floorval W
        model = sweep_model(state, ν, SWEEP_VZA; nstr, thr, floorval, W)
        qpm = Array(model.quad_points.qp_μ); wtm = Array(model.quad_points.wt_μ)
        μmin = minimum(qpm[wtm .> 0])
        src = source_prec(ν, solar_T, W)
        t = @elapsed res = CoreRT.rt_run_toa(model; i_band=1, sources=src)
        I = Float64(Array(res)[1, 1, 1])
        push!(rows, (; nstr, thr, floorval, W=string(W), μmin=Float64(μmin), I, t))
        @info "  I" I t
        open(SWEEP_OUT, "w") do io
            println(io, "# nstreams threshold floor float mu_min dtau_max dtau_over_mumin I elapsed_s")
            for r in rows
                fl = r.floorval === nothing ? 0.0 : r.floorval
                dmax = max(fl, r.thr * r.μmin)
                @printf(io, "%4d %9.1e %10.3e %8s %11.4e %11.4e %11.4e %.9e %8.1f\n",
                        r.nstr, r.thr, fl, r.W, r.μmin, dmax, dmax / r.μmin, r.I, r.t)
            end
        end
        GC.gc()
    end
    println(SWEEP_OUT)
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main_sweep()
