#!/usr/bin/env julia

"""Compute isolated and cumulative untruncated Stokes contributions.

Environment controls:
- `M_DIAG_MAX` — last Fourier order (default 99, from the beta-tail study).
- `NSTREAMS_DIAG` — diffuse streams per hemisphere (default 50, the minimum
  satisfying `2*N-1 >= M_DIAG_MAX`).
- `BY_M_OUT` — checkpoint/output path. Use distinct paths when comparing
  quadratures; checkpoints record and validate both resolution settings.
"""

# High-stream untruncated adding-doubling is not numerically reliable in
# Float32 (see truth_map_aerosols/HANDOFF_untruncated_artifact.md). Set this
# before loading generate_truth_map.jl, where the experiment FT is selected.
ENV["RRS_XCO2_FLOAT_TYPE"] = get(ENV, "RRS_XCO2_FLOAT_TYPE", "Float64")
include(joinpath(@__DIR__, "benchmark_aerosol_untruncated.jl"))
using JLD2
using DelimitedFiles

const M_DIAG_MAX = parse(Int, get(ENV, "M_DIAG_MAX", "99"))
const NSTREAMS_DIAG = parse(Int, get(ENV, "NSTREAMS_DIAG",
    string(cld(M_DIAG_MAX + 1, 2))))
const BY_M_OUT = get(ENV, "BY_M_OUT", joinpath(
    ROOT, "truth_map_aerosols",
    "untruncated_stokes_by_m_float64_nstreams$(NSTREAMS_DIAG)_mmax$(M_DIAG_MAX).jld2"))

FT === Float64 || error(
    "Untruncated moment benchmarks require Float64; got $FT. " *
    "Do not override RRS_XCO2_FLOAT_TYPE with Float32.")

function model_for_m_diagnostic(state, ν, vza)
    p = load_parameters()
    set_common!(p, state)
    p.architecture = CPU()
    p.vza = FT[vza, vza]
    p.vaz = copy(TEST_VAZ)
    req = full_mie_requirements(ν)
    p.nstreams = NSTREAMS_DIAG
    # `rt_set_streams` still consumes the normalized Legendre cap stored in
    # `l_trunc`; keep all resolution fields synchronized when mutating a
    # parsed parameter object. NoTruncation preserves the complete Mie Greek
    # series—the cap controls quadrature and the Fourier projection only.
    p.l_trunc = 2 * NSTREAMS_DIAG - 1
    p.stream_l_cap = p.l_trunc
    p.legacy_l_cap_override = nothing
    p.max_m = M_DIAG_MAX + 1     # solver convention: number of m orders
    p.truncation = vSmartMOM.Scattering.NoTruncation()
    surface = CoreRT.LambertianSurfaceScalar(FT(sum(state.coeff[1])))
    select_band!(p, 1, ν, surface)
    # The per-m callback is a diagnostic-only hook and is not exposed by the
    # external-solar TOA driver. Use the algebraically equivalent embedded-μ₀
    # representation here only; production defaults remain external-solar.
    return model_from_parameters(p; external_solar=false), req
end

function checkpoint_by_m(path; m, state, ν, contribution, elapsed_by_vza,
                         completed, requirements, final=false)
    kwargs = (; m, vza=TEST_VZA, vaz=TEST_VAZ, contribution,
              elapsed_by_vza, completed, requirements,
              nstreams=NSTREAMS_DIAG, m_max=M_DIAG_MAX,
              float_type=string(FT),
              external_solar=false,
              wavelength_nm=1e7 / ν[1], wavenumber_cm⁻¹=ν[1],
              state_index=state.index, final)
    tmp = path * ".tmp"
    jldsave(tmp; kwargs...)
    mv(tmp, path; force=true)  # atomic rename on the shared filesystem
end

function main_by_m()
    state = read_states()[9]
    ν = sparse_grid(1)
    solar_T = solar_interpolator()
    nv = length(TEST_VZA)
    requirements = full_mie_requirements(ν)
    m = collect(0:M_DIAG_MAX)
    contribution = zeros(FT, length(m), nv, 2, 3)
    elapsed_by_vza = zeros(Float64, nv)
    completed = falses(nv)

    if isfile(BY_M_OUT)
        saved = load(BY_M_OUT)
        saved["m"] == m || error("checkpoint m grid does not match this run")
        FT.(saved["vza"]) == TEST_VZA || error("checkpoint VZA grid does not match this run")
        saved["nstreams"] == NSTREAMS_DIAG ||
            error("checkpoint stream count does not match this run")
        saved["m_max"] == M_DIAG_MAX ||
            error("checkpoint m_max does not match this run")
        saved["float_type"] == string(FT) ||
            error("checkpoint float type does not match this run")
        size(saved["contribution"]) == size(contribution) ||
            error("checkpoint contribution shape does not match this run")
        contribution .= saved["contribution"]
        elapsed_by_vza .= saved["elapsed_by_vza"]
        completed .= saved["completed"]
        @info "resuming per-m diagnostic checkpoint" completed=count(completed) total=nv
    end

    for (iv, vza) in enumerate(TEST_VZA)
        if completed[iv]
            @info "skipping checkpointed VZA" vza iv nv
            continue
        end
        @info "untruncated per-m Stokes diagnostic" vza iv nv
        model, _ = model_for_m_diagnostic(state, ν, vza)
        model.quad_points.Nstreams == NSTREAMS_DIAG || error(
            "requested $NSTREAMS_DIAG streams, model built $(model.quad_points.Nstreams)")
        model.solver.m_max_bands[1] == M_DIAG_MAX || error(
            "requested m_max=$M_DIAG_MAX, model derived $(model.solver.m_max_bands[1])")
        _, sources = source_set(ν, false, solar_T)
        callback = function (step)
            J = Array(step.composite_layer.J₀⁻)
            for slot in values(step.composite_layer.J₀_by_src)
                J .+= Array(slot.J₀⁻)
            end
            iμ = CoreRT.nearest_point(step.qp_μ, cosd(vza))
            _, first_stokes, last_stokes = CoreRT.get_indices(iμ, step.pol_type)
            raw = @view J[first_stokes:last_stokes, 1, 1]
            for iaz in eachindex(TEST_VAZ)
                φ = TEST_VAZ[iaz]
                c = cosd(step.m * φ)
                s = sind(step.m * φ)
                contribution[step.m + 1, iv, iaz, 1] = step.weight * c * raw[1]
                contribution[step.m + 1, iv, iaz, 2] = step.weight * c * raw[2]
                contribution[step.m + 1, iv, iaz, 3] = step.weight * s * raw[3]
            end
            nothing
        end
        FTtype = CoreRT.float_type(model)
        elapsed_by_vza[iv] = @elapsed CoreRT._rt_run_column(
            vSmartMOM.InelasticScattering.noRS{FTtype}(), model, [1];
            sources, streams_callback=callback)
        completed[iv] = true
        checkpoint_by_m(BY_M_OUT; m, state, ν, contribution, elapsed_by_vza,
                        completed, requirements)
        @info "checkpointed completed VZA" vza elapsed=elapsed_by_vza[iv] completed=count(completed) total=nv
        GC.gc()
    end

    cumulative = cumsum(contribution; dims=1)
    signed_vza = vcat(-reverse(TEST_VZA[2:end]), TEST_VZA)
    # RelAz=180 maps to negative VZA; avoid duplicating nadir.
    isolated_signed = cat(reverse(contribution[:, 2:end, 2, :]; dims=2),
                          contribution[:, :, 1, :]; dims=2)
    cumulative_signed = cumsum(isolated_signed; dims=1)
    checkpoint_by_m(BY_M_OUT; m, state, ν, contribution, elapsed_by_vza,
                    completed, requirements, final=true)
    # Add the derived arrays only to the final file. The base checkpoint is
    # already complete, so interruption here remains safely resumable.
    jldopen(BY_M_OUT, "a+") do file
        file["signed_vza"] = signed_vza
        file["cumulative"] = cumulative
        file["isolated_signed"] = isolated_signed
        file["cumulative_signed"] = cumulative_signed
    end

    table = Matrix{Float64}(undef, length(m) * length(signed_vza), 9)
    row = 1
    for im in eachindex(m), iv in eachindex(signed_vza)
        table[row, :] .= (m[im], signed_vza[iv],
            isolated_signed[im, iv, 1], isolated_signed[im, iv, 2],
            isolated_signed[im, iv, 3], cumulative_signed[im, iv, 1],
            cumulative_signed[im, iv, 2], cumulative_signed[im, iv, 3],
            sum(elapsed_by_vza))
        row += 1
    end
    open(replace(BY_M_OUT, r"\.jld2$" => ".dat"), "w") do io
        println(io, "# m signed_vza dI_m dQ_m dU_m cumulative_I cumulative_Q cumulative_U total_elapsed_s")
        writedlm(io, table)
    end
    println(BY_M_OUT)
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main_by_m()
