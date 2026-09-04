#!/usr/bin/env julia

"""
Validate the elastic truth builder against the linearized retrieval forward model.

The default case is truth state 009: urban surface, all three aerosols present,
SIF off, and uniform CO2=380 ppm. The comparison deliberately uses SIF=0 so it
tests only physics shared by the two models; the retrieval still carries the
two-parameter `(SIF760,mSIF)` source with both values set to exact zero.

Both paths use Float32, nine streams, ABSCO v5.2 O2 plus the archived truth's
rebuilt HITRAN A-band H2O LUT, the same band-specific CO2/H2O ABSCO tables, the same
1000-hPa/16-layer materialized p/T/q profile, the same high-resolution solar
spectrum, and the same full-band aerosol interpolation anchors. The closure
continues through the fixed OCO analyzer, spectral-density conversion,
Gaussian convolution, and synthetic OCO-grid resampling. Every analytic
Jacobian column is independently perturbed through that same instrument
operator.

Environment controls:

- `CLOSURE_STATE` (default `9`; must be aerosol-on and SIF-off)
- `CUDA_DEVICE` (default `1`)
- `CLOSURE_REPORT` (default `inversion/retrieval_setup/ABSCO_CLOSURE.md`)
"""

ENV["RRS_XCO2_FLOAT_TYPE"] = "Float32"
ENV["NLAYERS"] = "16"
ENV["AEROSOL_NSTREAMS"] = "9"
ENV["TRUTH_SZA_DEG"] = "30"
ENV["TRUTH_VZA_DEG"] = "0"
ENV["TRUTH_RELATIVE_AZIMUTH_DEG"] = "0"

using CUDA
using Dates
using LinearAlgebra
using Printf

include(joinpath(@__DIR__, "..", "scripts", "generate_truth_map.jl"))
include(joinpath(@__DIR__, "OptimalEstimation.jl"))
include(joinpath(@__DIR__, "VSmartMOMForward.jl"))
using .VSmartMOMForward

const REPORT = get(ENV, "CLOSURE_REPORT",
    joinpath(@__DIR__, "retrieval_setup", "ABSCO_CLOSURE.md"))

function truth_retrieval_state(state, params)
    x = zeros(Float64, VSmartMOMForward.ACTIVE_STATE_COUNT)
    x[1] = 1000.0
    x[VSmartMOMForward.CO2_RANGE] .= state.xco2_ppm * 1e-6
    scale = VSmartMOMForward.tau_ref_per_aod760()
    aod760 = collect(state.aod550) ./ scale
    x[VSmartMOMForward.LOG_AOD_RANGE] .= log.(aod760)
    heights = [exp(aerosol.profile.μ)
               for aerosol in params.scattering_params.rt_aerosols]
    x[VSmartMOMForward.LOG_HEIGHT_RANGE] .= log.(heights)
    for iband in 1:3
        first_index = first(VSmartMOMForward.SURFACE_RANGE) + 3 * (iband - 1)
        x[first_index:first_index + 2] .= collect(state.coeff[iband])
    end
    # The closure state is intentionally SIF-free. These are physical zeros,
    # not fixed retrieval coordinates; both analytic SIF columns remain active.
    x[VSmartMOMForward.SIF_RANGE] .= 0.0
    return x
end

function difference_metrics(reference, candidate)
    ref = Float64.(reference)
    got = Float64.(candidate)
    size(ref) == size(got) || throw(DimensionMismatch(
        "closure arrays differ in shape: $(size(ref)) versus $(size(got))"))
    delta = got .- ref
    maxabs = maximum(abs, delta)
    rms = norm(delta) / sqrt(length(delta))
    refnorm = norm(ref)
    rel_l2 = iszero(refnorm) ? (iszero(norm(delta)) ? 0.0 : Inf) :
             norm(delta) / refnorm
    rel_max = maxabs / max(maximum(abs, ref), eps(Float64))
    return (; maxabs, rms, rel_l2, rel_max)
end

function exact_profile_checks(truth_model, retrieval_model)
    t, r = truth_model.profile, retrieval_model.profile
    layer_values(value, nlayer) = value isa AbstractVector ? value : fill(value, nlayer)
    nlayer = length(t.T)
    return (
        T = t.T == r.T,
        p_half = t.p_half == r.p_half,
        p_full = t.p_full == r.p_full,
        q = t.q == r.q,
        vmr_h2o = t.vmr_h2o == r.vmr_h2o,
        vcd_dry = t.vcd_dry == r.vcd_dry,
        O2 = layer_values(t.vmr["O2"], nlayer) ==
             layer_values(r.vmr["O2"], nlayer),
        CO2 = layer_values(t.vmr["CO2"], nlayer) ==
              layer_values(r.vmr["CO2"], nlayer),
    )
end

function canonical_global_jacobian(result)
    local_jacobian = Array(result.toa_jacobian)
    global_jacobian = Array(globalize_jacobian(local_jacobian, result.layout))
    size(global_jacobian, 1) == 1 || error(
        "closure test requires exactly one viewing zenith angle")
    size(global_jacobian, 2) >= 3 || error(
        "closure test requires at least I, Q, and U")
    return Array(@view(global_jacobian[1, 1:3, :, :]))
end

"""
Check the fixed instrument operator against central perturbations along every
high-resolution analytic Jacobian column. This isolates Mueller projection,
the per-cm⁻¹ to per-nm density conversion, Gaussian convolution, and spectral
resampling from the RT derivative itself.
"""
function instrument_jacobian_checks(wavelength_nm, stokes, jacobian,
                                    coefficients, spec)
    instrument = VSmartMOMForward.SyntheticOCO2
    processed = instrument.process_stokes_jacobian(
        wavelength_nm, jacobian, coefficients, spec)
    base64 = Float64.(stokes)
    scale = max(maximum(abs, base64), 1.0)
    metrics = NamedTuple[]
    for iparam in axes(jacobian, 3)
        column = Float64.(@view jacobian[:, :, iparam])
        amplitude = maximum(abs, column)
        if iszero(amplitude)
            push!(metrics, difference_metrics(
                zeros(size(processed, 1)), @view(processed[:, iparam])))
            continue
        end
        # Give the instrument-only perturbation a resolvable amplitude. This
        # does not perturb or rerun the RT model, whose analytic derivative is
        # the input under test here.
        step = 0.1 * scale / amplitude
        plus = instrument.process_stokes_spectrum(
            wavelength_nm, base64 .+ step .* column, coefficients, spec)
        minus = instrument.process_stokes_spectrum(
            wavelength_nm, base64 .- step .* column, coefficients, spec)
        finite_difference = (plus .- minus) ./ (2step)
        push!(metrics, difference_metrics(
            @view(processed[:, iparam]), finite_difference))
    end
    return processed, metrics
end

function raman_shoulder_convolution_check(core_wavenumber, core_stokes,
                                          coefficients, spec)
    solve_wavenumber, keep = o2_solve_grid(core_wavenumber)
    length(keep) == length(core_wavenumber) || error(
        "Raman solve grid does not retain exactly the O2 output grid")
    grid_metric = difference_metrics(core_wavenumber, solve_wavenumber[keep])

    # The solve-only shoulders are deliberately filled with a conspicuous
    # value. They are many Gaussian widths outside the measurement support;
    # after the retained core is selected they must be irrelevant to the ILS.
    padded = fill(eltype(core_stokes)(12345),
                  size(core_stokes, 1), length(solve_wavenumber))
    padded[:, keep] .= core_stokes
    instrument = VSmartMOMForward.SyntheticOCO2
    core_processed = instrument.process_stokes_spectrum(
        1.0e7 ./ Float64.(core_wavenumber), core_stokes, coefficients, spec)
    padded_processed = instrument.process_stokes_spectrum(
        1.0e7 ./ Float64.(solve_wavenumber), padded, coefficients, spec)
    return grid_metric, difference_metrics(core_processed, padded_processed)
end

function truth_model_for_band(state, grids, iband)
    iband == 1 && return build_o2(state, grids[1], grids[1])
    return build_co2(state, iband, grids[iband], grids[iband])
end

function main()
    CUDA.functional() || error("CUDA is not functional")
    CUDA.device!(parse(Int, get(ENV, "CUDA_DEVICE", "1")))
    state_index = parse(Int, get(ENV, "CLOSURE_STATE", "9"))
    state = read_states()[state_index]
    any(>(0), state.aod550) || error("closure state must contain aerosols")
    iszero(state.sif_angular_integral760) ||
        error("closure state must have SIF=0")

    retrieval_params = VSmartMOMForward.prepare_retrieval_parameters(;
        architecture=:GPU, float_type=Float32, nstreams=9)
    x = truth_retrieval_state(state, retrieval_params)
    physical = VSmartMOMForward.apply_retrieval_state!(
        retrieval_params, x, VSmartMOMForward.tau_ref_per_aod760();
        fixed_upper_co2_vmr=state.xco2_ppm * 1e-6)
    all(iszero, (physical.SIF760, physical.mSIF)) || error(
        "closure state did not map to exact zero SIF")

    retrieval_model, retrieval_lin = model_from_parameters(
        OCO_RRS_synth(), retrieval_params; external_solar=true)
    grids = RRSXCO2Common.basis_grids(Float32)
    solar_T = RRSXCO2Common.solar_interpolator(Float32)
    instrument = VSmartMOMForward.SyntheticOCO2
    coefficients = instrument.read_representative_coefficients(
        VSmartMOMForward.DEFAULT_COEFFICIENT_PATH)
    rows = NamedTuple[]
    profiles = NamedTuple[]
    shoulder_checks = nothing
    elastic_shoulder_check = nothing
    active_o2_check = nothing

    for iband in 1:3
        truth_model = truth_model_for_band(state, grids, iband)
        push!(profiles, exact_profile_checks(truth_model, retrieval_model))

        truth_sources = source_set(grids[iband], false, solar_T)[2]
        retrieval_sources =
            VSmartMOMForward.RRSXCO2Common.sources_for_band(
                retrieval_params, iband; SIF760=0, mSIF=0,
                solar_T=solar_T)

        truth_toa = toa3(CoreRT.rt_run_toa(
            truth_model; i_band=1, sources=truth_sources))
        if iband == 1
            # Reproduce one actual 256-point aerosol truth chunk. Its noRS
            # radiance must equal the same wavelengths from the canonical
            # full-band solve even though the RT batch also carries ±234 cm⁻¹
            # Raman source shoulders. This is distinct from the convolution
            # check below: it verifies that solve batching does not alter the
            # elastic spectrum retained for the corrected measurement.
            width = min(256, length(grids[iband]))
            first_index = fld(length(grids[iband]) - width, 2) + 1
            chunk_range = first_index:first_index + width - 1
            chunk_wavenumber = grids[iband][chunk_range]
            solve_wavenumber, keep = o2_solve_grid(chunk_wavenumber)
            shoulder_model = build_o2(state, grids[iband], solve_wavenumber)
            shoulder_sources = source_set(solve_wavenumber, false, solar_T)[2]
            shoulder_toa = toa3(CoreRT.rt_run_toa(
                shoulder_model; i_band=1, sources=shoulder_sources), keep)
            elastic_shoulder_check = difference_metrics(
                truth_toa[:, chunk_range], shoulder_toa)
            CUDA.synchronize()
            GC.gc()
            CUDA.reclaim()
        end
        retrieval_toa = Array(CoreRT.rt_run_toa(
            retrieval_model; i_band=iband, sources=retrieval_sources))[1, 1:3, :]
        lin_result = CoreRT.rt_run_lin(
            retrieval_model, retrieval_lin; i_band=iband,
            sources=retrieval_sources)
        lin_toa = Array(lin_result.toa)[1, 1:3, :]
        global_jacobian = canonical_global_jacobian(lin_result)
        spec = instrument.BAND_SPECS[iband]
        wavelength_nm = 1.0e7 ./ Float64.(grids[iband])
        truth_oco = instrument.process_stokes_spectrum(
            wavelength_nm, truth_toa, coefficients[spec.name], spec)
        retrieval_oco = instrument.process_stokes_spectrum(
            wavelength_nm, retrieval_toa, coefficients[spec.name], spec)
        lin_oco = instrument.process_stokes_spectrum(
            wavelength_nm, lin_toa, coefficients[spec.name], spec)
        _, instrument_columns = instrument_jacobian_checks(
            wavelength_nm, lin_toa, global_jacobian,
            coefficients[spec.name], spec)
        if iband == 1
            active_scene = joinpath(
                RRSXCO2Common.ROOT, "truth_map", "aerosol_chunked",
                @sprintf("hiressim_%03d.nc", state.index))
            isfile(active_scene) || error(
                "active A-band truth scene is missing: $active_scene")
            active_o2 = NCDataset(active_scene) do dataset
                get(dataset.attrib, "o2_truth_regenerated", 0) == 1 || error(
                    "active closure scene lacks regenerated-O2 provenance")
                get(dataset.attrib, "o2_absco_regeneration_complete", 0) == 1 ||
                    error("active closure scene has incomplete O2 regeneration")
                get(dataset.attrib, "o2_truth_reused", 1) == 0 || error(
                    "active closure scene is still marked as reused O2 truth")
                Int(dataset.attrib["o2_core_grid_version"]) == 2 || error(
                    "active closure scene does not use exact O2 core nodes")
                Float32.(dataset["radiance_rayleigh_o2a"][:, :])
            end
            active_o2_oco = instrument.process_stokes_spectrum(
                wavelength_nm, active_o2, coefficients[spec.name], spec)
            active_o2_check = (
                high_resolution=difference_metrics(active_o2, retrieval_toa),
                oco=difference_metrics(active_o2_oco, retrieval_oco),
            )
            shoulder_checks = raman_shoulder_convolution_check(
                grids[iband], lin_toa, coefficients[spec.name], spec)
        end

        optical = (
            tau_abs = difference_metrics(
                truth_model.τ_abs[1], retrieval_model.τ_abs[iband]),
            tau_rayleigh = difference_metrics(
                truth_model.τ_rayl[1], retrieval_model.τ_rayl[iband]),
            tau_aerosol = difference_metrics(
                truth_model.τ_aer[1], retrieval_model.τ_aer[iband]),
        )
        push!(rows, (
            band=RRSXCO2Common.BAND_NAMES[iband],
            truth_vs_retrieval=difference_metrics(truth_toa, retrieval_toa),
            retrieval_vs_linearized=difference_metrics(retrieval_toa, lin_toa),
            oco_truth_vs_retrieval=difference_metrics(truth_oco, retrieval_oco),
            oco_retrieval_vs_linearized=difference_metrics(retrieval_oco, lin_oco),
            instrument_columns=instrument_columns,
            optical=optical,
        ))
        CUDA.synchronize()
        GC.gc()
        CUDA.reclaim()
    end

    profile_pass = all(all(values(check)) for check in profiles)
    # Float32 GPU/CPU Mie construction is allowed roundoff, but a physically
    # meaningful mismatch is not. L2 is the primary whole-spectrum criterion;
    # the max metric remains in the report for locating isolated deviations.
    radiance_pass = all(
        row.truth_vs_retrieval.rel_l2 <= 2e-5 &&
        row.retrieval_vs_linearized.rel_l2 <= 1e-5 &&
        row.oco_truth_vs_retrieval.rel_l2 <= 2e-5 &&
        row.oco_retrieval_vs_linearized.rel_l2 <= 1e-5 for row in rows)
    optical_pass = all(
        metric.rel_l2 <= 2e-5
        for row in rows for metric in values(row.optical))
    instrument_pass = all(
        metric.maxabs <= 2e-9 || metric.rel_l2 <= 2e-9
        for row in rows for metric in row.instrument_columns)
    shoulder_grid, shoulder_convolution = shoulder_checks
    shoulder_grid_pass = shoulder_grid.maxabs <= 2eps(Float32) *
                         maximum(abs, grids[1])
    shoulder_convolution_pass = shoulder_convolution.maxabs <= 2e-12
    elastic_shoulder_pass = elastic_shoulder_check.rel_l2 <= 2e-5
    active_o2_pass = active_o2_check.high_resolution.rel_l2 <= 2e-5 &&
                     active_o2_check.oco.rel_l2 <= 2e-5
    shoulder_pass = shoulder_grid_pass && shoulder_convolution_pass &&
                    elastic_shoulder_pass
    passed = profile_pass && optical_pass && radiance_pass &&
             instrument_pass && shoulder_pass && active_o2_pass

    mkpath(dirname(REPORT))
    open(REPORT, "w") do io
        println(io, "# ABSCO truth/retrieval forward-model closure")
        println(io)
        println(io, "- Generated: `$(now(UTC))` UTC")
        println(io, "- Truth state: `$(lpad(state.index, 3, '0'))` " *
                    "($(state.surface), $(state.aerosol_case), SIF off, " *
                    "$(state.xco2_ppm) ppm CO2)")
        println(io, "- Precision/backend: `Float32`, CUDA device " *
                    "`$(get(ENV, "CUDA_DEVICE", "1"))`")
        println(io, "- Spectroscopy: ABSCO `$(RRSXCO2Common.ABSCO_VERSION)`")
        println(io, "- A-band H2O: rebuilt HITRAN LUT retained from the " *
                    "archived ABSCO-O2 truth configuration")
        println(io, "- Retrieval SIF source: `SIF760 + mSIF*(nu-nu760)`; " *
                    "both values are zero for this closure state, while both " *
                    "Jacobian columns remain active")
        println(io, "- Overall result: **$(passed ? "PASS" : "FAIL")**")
        println(io)
        println(io, "Exact profile-array checks:")
        println(io)
        println(io, "| Band | T | p_half | p_full | q | H2O VMR | dry column | O2 VMR | CO2 VMR |")
        println(io, "|---|---:|---:|---:|---:|---:|---:|---:|---:|")
        for (iband, check) in pairs(profiles)
            println(io, "| $(RRSXCO2Common.BAND_NAMES[iband]) | " *
                join(string.(values(check)), " | ") * " |")
        end
        println(io)
        println(io, "All profile arrays matched element-for-element: **$profile_pass**.")
        println(io)
        println(io, "| Band | comparison | max abs | relative L2 | relative max |")
        println(io, "|---|---|---:|---:|---:|")
        for row in rows
            for (label, metric) in (
                    ("truth vs retrieval noRS", row.truth_vs_retrieval),
                    ("retrieval noRS vs linearized forward", row.retrieval_vs_linearized))
                @printf(io, "| %s | %s | %.9e | %.9e | %.9e |\n",
                row.band, label, metric.maxabs,
                        metric.rel_l2, metric.rel_max)
            end
            for (label, metric) in (
                    ("OCO-processed truth vs retrieval noRS", row.oco_truth_vs_retrieval),
                    ("OCO-processed retrieval vs linearized forward",
                     row.oco_retrieval_vs_linearized))
                @printf(io, "| %s | %s | %.9e | %.9e | %.9e |\n",
                        row.band, label, metric.maxabs,
                        metric.rel_l2, metric.rel_max)
            end
        end
        println(io)
        println(io, "| Band | optical field | max abs | relative L2 |")
        println(io, "|---|---|---:|---:|")
        for row in rows, (name, metric) in pairs(row.optical)
            @printf(io, "| %s | %s | %.9e | %.9e |\n",
                    row.band, String(name), metric.maxabs, metric.rel_l2)
        end
        println(io)
        println(io, "## Instrument/Jacobian consistency")
        println(io)
        println(io, "The same fixed operator—`M11*I - M12*Q + M13*U`, " *
                    "per-cm⁻¹ to per-nm conversion, Gaussian ILS, then " *
                    "synthetic OCO sampling—was applied to the forward " *
                    "spectrum and all 30 analytic Jacobian columns.")
        println(io)
        println(io, "| Band | worst column max abs | worst column relative L2 | result |")
        println(io, "|---|---:|---:|---:|")
        for row in rows
            worst_abs = maximum(metric.maxabs for metric in row.instrument_columns)
            worst_rel = maximum(metric.rel_l2 for metric in row.instrument_columns)
            ok = all(metric.maxabs <= 2e-9 || metric.rel_l2 <= 2e-9
                     for metric in row.instrument_columns)
            @printf(io, "| %s | %.9e | %.9e | %s |\n",
                    row.band, worst_abs, worst_rel, ok ? "PASS" : "FAIL")
        end
        println(io)
        @printf(io, "Raman solve-grid retained-node max difference: `%.9e cm⁻¹`.\n\n",
                shoulder_grid.maxabs)
        @printf(io, "OCO convolution difference when solve-only Raman shoulders are present: max abs `%.9e`, relative L2 `%.9e` (**%s**).\n\n",
                shoulder_convolution.maxabs, shoulder_convolution.rel_l2,
                shoulder_convolution_pass ? "PASS" : "FAIL")
        @printf(io, "Elastic noRS difference for a production-size 256-point core solved with ±234 cm⁻¹ Raman shoulders: max abs `%.9e`, relative L2 `%.9e`.\n\n",
                elastic_shoulder_check.maxabs, elastic_shoulder_check.rel_l2)
        println(io, "## Active regenerated A-band product closure")
        println(io)
        println(io, "The active state file's regenerated Rayleigh A-band was " *
                    "compared directly with the retrieval noRS forward model.")
        println(io)
        @printf(io, "High-resolution regenerated-truth/retrieval difference: max abs `%.9e`, relative L2 `%.9e`.\n\n",
                active_o2_check.high_resolution.maxabs,
                active_o2_check.high_resolution.rel_l2)
        @printf(io, "OCO-processed regenerated-truth/retrieval difference: max abs `%.9e`, relative L2 `%.9e` (**%s**).\n\n",
                active_o2_check.oco.maxabs, active_o2_check.oco.rel_l2,
                active_o2_pass ? "PASS" : "FAIL")
        println(io, "Raman shoulders are used by the RRS solve, then discarded " *
                    "with the retained-core index before the Mueller/ILS " *
                    "operator. Their physical Raman redistribution into the " *
                    "core is retained; the shoulder samples themselves never " *
                    "enter convolution.")
        println(io)
        println(io, "The test uses the canonical complete output band as the aerosol " *
                    "interpolation anchor. Raman chunks may extend beyond it, but their " *
                    "aerosol phase/extinction interpolation is anchored to these same " *
                    "band endpoints so chunking cannot change core-band elastic physics.")
    end

    println("closure_report=$REPORT")
    for row in rows
        @printf("%-10s truth/retrieval relL2=%.9e maxabs=%.9e retrieval/lin relL2=%.9e\n",
                row.band, row.truth_vs_retrieval.rel_l2,
                row.truth_vs_retrieval.maxabs,
                row.retrieval_vs_linearized.rel_l2)
    end
    @printf("instrument columns=%s shoulder maxabs=%.9e grid maxabs=%.9e elastic shoulder relL2=%.9e\n",
            instrument_pass ? "PASS" : "FAIL",
            shoulder_convolution.maxabs, shoulder_grid.maxabs,
            elastic_shoulder_check.rel_l2)
    @printf("active O2 product high-resolution relL2=%.9e OCO relL2=%.9e (%s)\n",
            active_o2_check.high_resolution.rel_l2,
            active_o2_check.oco.rel_l2,
            active_o2_pass ? "PASS" : "FAIL")
    passed || error("truth/retrieval closure failed; inspect $REPORT")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
