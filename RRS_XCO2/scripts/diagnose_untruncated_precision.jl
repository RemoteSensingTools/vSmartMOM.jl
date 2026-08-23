#!/usr/bin/env julia

"""Separate the three candidate causes of the untruncated angular artifact.

The δBGE-truncated references (l_trunc = 8…32, nstreams = 5…17) agree with one
another to ~0.1%, yet the "untruncated" runs (nstreams = 50 and 111) disagree
with them, and with each other, by up to 5%.  Three things change together in
those runs:

1. `nstreams`  — the size of the diffuse operator (Nquad·n_stokes).
2. `l_trunc`   — the retained Legendre support, coupled as `2·nstreams − 1`.
3. nothing else … but everything is computed in **Float32**.

This script varies them independently at a single wavelength and two VZAs so
the responsible factor can be identified.  Each case reports TOA I/Q/U.

Environment:
- `PREC_VZAS`  — comma-separated VZAs (default "20,45").
- `PREC_OUT`   — output ASCII table.
"""

include(joinpath(@__DIR__, "benchmark_aerosol_untruncated.jl"))
using Printf

const PREC_VZAS = parse.(Float64, split(get(ENV, "PREC_VZAS", "20,45"), ','))
const PREC_OUT = get(ENV, "PREC_OUT",
    joinpath(ROOT, "truth_map_aerosols", "untruncated_precision_probe.dat"))

# ---------------------------------------------------------------------------
# Parameter construction at an arbitrary working precision.
# ---------------------------------------------------------------------------

"Mirror of `set_common!`/`truncate_profile!` at working precision `W`."
function set_common_prec!(p, state, ::Type{W}) where {W}
    p0, T0, q0 = Float64.(p.p), Float64.(p.T), Float64.(p.q)
    psurf = 900.0
    k = searchsortedlast(p0, psurf - eps(psurf))
    pcent = (p0[1:end-1] .+ p0[2:end]) ./ 2
    function at_pressure(v, pr)
        hi = searchsortedfirst(pcent, pr)
        hi <= 1 && return v[1]
        hi > length(v) && return v[end]
        lo = hi - 1
        t = (log(pr) - log(pcent[lo])) / (log(pcent[hi]) - log(pcent[lo]))
        return v[lo] + t * (v[hi] - v[lo])
    end
    p.float_type = W
    p.architecture = CPU()
    p.sza = W(30)
    p.nstreams = 5
    p.profile_reduction_n = NLAYERS
    p.p = W.(vcat(p0[1:k], psurf))
    p.T = W.(vcat(T0[1:k-1], at_pressure(T0, psurf)))
    p.q = W.(vcat(q0[1:k-1], at_pressure(q0, psurf)))
    p.absorption_params.vmr["CO2"] = W(state.xco2_ppm * 1e-6)
    for (aer, τ) in zip(p.scattering_params.rt_aerosols, state.aod550)
        aer.τ_ref = W(τ)
    end
    return p
end

"""Build a model for one probe case.

`nstr`   — weighted streams per hemisphere.
`lcap`   — Legendre/Fourier support.  `nothing` couples it to `2·nstr − 1`.
`trunc`  — `:none` for NoTruncation, or an `Int` δBGE cap.
`W`      — working float type.
"""
function probe_model(state, ν, vza; nstr, lcap=nothing, trunc=:none, W=Float32,
                     floor_mult=nothing)
    p = load_parameters()
    set_common_prec!(p, state, W)
    if floor_mult !== nothing
        # NOTE (measured): this override is a NO-OP for these runs, and that
        # is itself the finding. `load_parameters()` parses a Float64 YAML, so
        # `params.numerics.dτ_min_floor` is already 1024·eps(Float64)=2.3e-13;
        # `_convert_numerics` only CASTS it when the scripts later set
        # `float_type = Float32`, never recomputing the Float32 default
        # (1024·eps(F32)=1.2e-4). The floor is therefore inactive throughout,
        # `dτ_max = 0.001·μ_min` is honoured, and the Float32 error at high
        # nstreams comes from elemental layers too THIN for Float32 plus the
        # long doubling ladder — not from the floor. See
        # `diagnose_doubling_sweep.jl`.
        p.numerics = CoreRT.RTNumericalParameters{W}(;
            dτ_max_threshold = p.numerics.dτ_max_threshold,
            dτ_min_floor = W(floor_mult) * eps(W),
            blas_threads = p.numerics.blas_threads,
            verbose = p.numerics.verbose)
    end
    p.vza = W[vza, vza]
    p.vaz = W[0, 180]
    # `model_from_parameters` builds the quadrature from `params.l_trunc`
    # positionally (`Nquad = (l_trunc + 2) ÷ 2`); `params.nstreams` is only a
    # parse-time alias. Set `l_trunc` to size the OPERATOR, and use
    # `max_m` + `truncation` to set the Legendre/Fourier SUPPORT independently.
    L = lcap === nothing ? 2 * nstr - 1 : lcap
    p.l_trunc = 2 * nstr - 2
    p.nstreams = nstr
    p.stream_l_cap = 2 * nstr - 1
    p.legacy_l_cap_override = nothing
    p.max_m = L + 1
    p.truncation = trunc === :none ? vSmartMOM.Scattering.NoTruncation() :
                                     vSmartMOM.Scattering.δBGE{W}(trunc, zero(W))
    surface = CoreRT.LambertianSurfaceScalar(W(sum(state.coeff[1])))
    p.spec_bands = [W.(ν)]
    p.brdf = [surface]
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

"Solar/SIF source pair at working precision `W`."
function source_prec(ν, solar_T, ::Type{W}) where {W}
    P = W.(planck_spectrum_wn(W(5777), collect(W.(ν))) .* W(2.1629e-5 * π))
    F0 = zeros(W, 3, length(ν))
    F0[1, :] .= W.(solar_T.(ν)) .* P
    return SolarBeam(F₀=F0) + SurfaceSIF(SIF₀=zeros(W, 3, length(ν)))
end

function run_probe(state, ν, solar_T, vza; kwargs...)
    W = get(kwargs, :W, Float32)
    model = probe_model(state, ν, vza; kwargs...)
    src = source_prec(ν, solar_T, W)
    elapsed = @elapsed result = CoreRT.rt_run_toa(model; i_band=1, sources=src)
    r = Float64.(Array(result)[:, 1:3, :])
    return (; I=r[1, 1, 1], Q=r[1, 2, 1], U=r[1, 3, 1], elapsed,
            nstreams=model.quad_points.Nstreams, nquad=model.quad_points.Nquad,
            m_max=model.solver.m_max_bands[1], l_max=model.solver.l_max[1])
end

# ---------------------------------------------------------------------------
# Case matrix
# ---------------------------------------------------------------------------
#
# ref32/ref64   — converged δBGE reference (small operator).
# big32/big64   — LARGE operator, but the SAME small phase support (δBGE 32).
#                 Isolates operator size / precision from Legendre support.
# unt50_*       — NoTruncation, l = 99  (reproduces the 50-stream dataset).
# unt50l32_*    — NoTruncation quadrature size with l capped at 32.
#
const CASES = [
    (name="ref32_n17_l32",     nstr=17,  lcap=32, trunc=32,    W=Float32),
    (name="ref64_n17_l32",     nstr=17,  lcap=32, trunc=32,    W=Float64),
    (name="big32_n50_l32",     nstr=50,  lcap=32, trunc=32,    W=Float32),
    (name="big64_n50_l32",     nstr=50,  lcap=32, trunc=32,    W=Float64),
    (name="big32_n111_l32",    nstr=111, lcap=32, trunc=32,    W=Float32),
    (name="big64_n111_l32",    nstr=111, lcap=32, trunc=32,    W=Float64),
    (name="unt32_n50_l99",     nstr=50,  lcap=nothing, trunc=:none, W=Float32),
    (name="unt64_n50_l99",     nstr=50,  lcap=nothing, trunc=:none, W=Float64),
    # Float32, but with the dτ floor lowered so `dτ_max = 0.001·μ_min` is
    # actually honoured. If these move to the reference, the floor clamp — not
    # generic Float32 roundoff — is the cause.
    (name="f32lo_n50_l32",     nstr=50,  lcap=32, trunc=32, W=Float32, floor_mult=1),
    (name="f32lo_n111_l32",    nstr=111, lcap=32, trunc=32, W=Float32, floor_mult=1),
]

function main_precision()
    state = read_states()[9]
    ν = sparse_grid(1)
    solar_T = solar_interpolator()
    rows = []
    for vza in PREC_VZAS, c in CASES
        @info "precision probe" vza c.name
        r = run_probe(state, ν, solar_T, vza;
                      nstr=c.nstr, lcap=c.lcap, trunc=c.trunc, W=c.W,
                      floor_mult=get(c, :floor_mult, nothing))
        @info "  result" I=r.I Q=r.Q t=r.elapsed nquad=r.nquad m_max=r.m_max l_max=r.l_max
        push!(rows, (; vza, c.name, r...))
        open(PREC_OUT, "w") do io
            println(io, "# vza case I Q U nstreams nquad m_max l_max elapsed_s")
            for x in rows
                @printf(io, "%5.1f %20s %.9e %.9e %.9e %4d %4d %4d %4d %9.2f\n",
                        x.vza, x.name, x.I, x.Q, x.U,
                        x.nstreams, x.nquad, x.m_max, x.l_max, x.elapsed)
            end
        end
        GC.gc()
    end
    println(PREC_OUT)
end

# Guarded so this file can be `include`d as a helper library.
abspath(PROGRAM_FILE) == abspath(@__FILE__) && main_precision()
