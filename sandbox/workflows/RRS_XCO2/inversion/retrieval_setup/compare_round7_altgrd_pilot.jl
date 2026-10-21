#!/usr/bin/env julia
# Compare the Round 7_altgrd pilot (FixedAltitudeGrid) with the completed
# Round 7 (FixedPressureGrid) retrievals, scene by scene and member by member.
# Read-only on both output trees. Run with --project=<worktree> and the Round 7
# environment exported (configure_round7_altgrd off). CPU only.
#   julia compare_round7_altgrd_pilot.jl OUT.csv OUT_summary.md [JACOBIAN_PRESSURE.jls JACOBIAN_ALTITUDE.jls]
# The optional .jls files are outputs of vertical_grid_jacobian_states.jl.

using LinearAlgebra
using NCDatasets
using Printf
using Serialization
using Statistics

const W = dirname(Base.active_project())
const INV = joinpath(W, "sandbox", "workflows", "RRS_XCO2", "inversion")
for name in ("OptimalEstimation.jl", "RetrievalCases.jl", "RetrievalState.jl",
             "RetrievalOutput.jl", "Round4SIFTruthConvention.jl",
             "Round5FixedSIF.jl", "Round5RetrievalCampaign.jl")
    include(joinpath(INV, name))
end
using .RetrievalCases, .Round5FixedSIF
import .Round5RetrievalCampaign as R5
include(joinpath(INV, "VSmartMOMForward.jl"))
using .VSmartMOMForward
using vSmartMOM

const BOTTOM = joinpath(ENV["RRS_XCO2_DATA_ROOT"], "bottom_layer_XCO2_retrievals")
const ROUND7 = joinpath(BOTTOM, "round7_fixed_sif_imperfect_correction")
const ALTGRD = joinpath(BOTTOM, "round7_altgrd_fixed_sif_imperfect_correction")

output_path(root, truth, member) = joinpath(root,
    truth.sif_case == :off ? "retrievals_nosif" : "retrievals_sif", "corrected",
    @sprintf("retrieval_state%03d_perturbation%02d.nc", truth.state_index, member))

function read_retrieval(path)
    isfile(path) || return nothing
    NCDataset(path) do ds
        (; xco2=Float64(ds["XCO2"][]), state=Float64.(ds["final_state"][:]),
           chi2=Float64.(ds["final_band_reduced_chi_squared"][:]),
           converged=Int(ds.attrib["converged"]), fit_ok=Int(ds.attrib["fit_quality_ok"]),
           trials=Int(ds.dim["trial"]), iterations=Int(maximum(ds["iteration_index"][:])),
           Shat=Float64.(ds["posterior_covariance"][:, :]),
           outcome=Int(ds.attrib["outcome"]), divergences=Int(ds.attrib["divergence_count"]),
           model=String(ds.attrib["retrieval_state_model"]),
           grid=haskey(ds["XCO2"].attrib, "vertical_grid_convention") ?
                String(ds["XCO2"].attrib["vertical_grid_convention"]) : "FixedPressureGrid")
    end
end

read_names(path) = NCDataset(ds -> split(String(ds.attrib["parameter_names"])), path)

"Gradient of the evaluator's own XCO2 diagnostic (it depends only on psurf and the CO2 coordinates 2:13)."
function xco2_gradient(evaluator, map, state)
    h = zeros(length(state))
    for (j, δ) in ((1, 0.1), ((j, 1e-7) for j in 2:13)...)
        plus, minus = copy(state), copy(state)
        plus[j] += δ; minus[j] -= δ
        h[j] = (column_averaged_co2_ppm(evaluator, expand_round5_state(map, plus)) -
                column_averaged_co2_ppm(evaluator, expand_round5_state(map, minus))) / (2δ)
    end
    return h
end

function main(csv_path, summary_path, jacobian_paths...)
    truths = read_truth_cases(ENV["RETRIEVAL_TRUTH_TABLE"])
    pressure = OCOForwardEvaluator(; architecture=:CPU, float_type=Float32, nstreams=9,
        coefficient_path=ENV["RETRIEVAL_STOKES_COEFFICIENT_PATH"],
        component_path=ENV["RETRIEVAL_SCENE_COMPONENTS_PATH"])
    altitude = OCOForwardEvaluator(pressure.base_parameters, pressure.coefficients,
        pressure.tau_ref_scale, pressure.fixed_upper_co2_vmr, pressure.solar_transmission,
        pressure.synchronize_backend, FixedAltitudeGrid())
    rows = NamedTuple[]
    for truth in truths
        map = R5.resolve_round5_sif_map([truth])
        set_fixed_upper_co2_ppm!(pressure, truth.fixed_upper_co2_ppm)
        set_fixed_upper_co2_ppm!(altitude, truth.fixed_upper_co2_ppm)
        noisy = filter(!isnothing, [read_retrieval(output_path(ROUND7, truth, m)) for m in 1:10])
        scatter_xco2 = length(noisy) > 1 ? std([r.xco2 for r in noisy]) : NaN
        scatter_psurf = length(noisy) > 1 ? std([r.state[1] for r in noisy]) : NaN
        for member in (11, 1)
            r7 = read_retrieval(output_path(ROUND7, truth, member))
            ag = read_retrieval(output_path(ALTGRD, truth, member))
            (r7 === nothing || ag === nothing) && continue
            ag.model == "round7_altgrd_fixed_sif_imperfect_correction" && ag.grid == "FixedAltitudeGrid" &&
                r7.model == "round7_fixed_sif_imperfect_correction" ||
                error("unexpected state model/grid for scene $(truth.state_index) member $member")
            # Each round's XCO2 is recomputed with its own convention as a check.
            x7 = column_averaged_co2_ppm(pressure, expand_round5_state(map, r7.state))
            xa = column_averaged_co2_ppm(altitude, expand_round5_state(map, ag.state))
            h7 = xco2_gradient(pressure, map, r7.state)
            ha = xco2_gradient(altitude, map, ag.state)
            push!(rows, (; scene=truth.state_index, member, surface=truth.surface,
                aerosol=truth.aerosol_case, sif=truth.sif_case, truth_xco2=truth.xco2_ppm,
                xco2_r7=r7.xco2, xco2_ag=ag.xco2, recheck_r7=x7 - r7.xco2, recheck_ag=xa - ag.xco2,
                psurf_r7=r7.state[1], psurf_ag=ag.state[1],
                sigma_xco2_r7=sqrt(dot(h7, r7.Shat * h7)), sigma_xco2_ag=sqrt(dot(ha, ag.Shat * ha)),
                sigma_psurf_r7=sqrt(r7.Shat[1, 1]), sigma_psurf_ag=sqrt(ag.Shat[1, 1]),
                chi2_r7=r7.chi2, chi2_ag=ag.chi2, conv_r7=r7.converged, conv_ag=ag.converged,
                fit_r7=r7.fit_ok, fit_ag=ag.fit_ok, iter_r7=r7.iterations, iter_ag=ag.iterations,
                trials_r7=r7.trials, trials_ag=ag.trials, scatter_xco2, scatter_psurf,
                n_noisy=length(noisy), model_ag=ag.model, grid_ag=ag.grid,
                outcome_ag=ag.outcome, divergences_ag=ag.divergences,
                state_r7=r7.state, state_ag=ag.state, names=read_names(output_path(ROUND7, truth, member))))
        end
    end
    open(csv_path, "w") do io
        println(io, "scene,member,surface,aerosol,sif,truth_xco2,xco2_r7,xco2_altgrd,dxco2,",
            "psurf_r7,psurf_altgrd,dpsurf,sigma_xco2_r7,sigma_xco2_altgrd,sigma_psurf_r7,",
            "sigma_psurf_altgrd,chi2_o2a_r7,chi2_weak_r7,chi2_strong_r7,chi2_o2a_altgrd,",
            "chi2_weak_altgrd,chi2_strong_altgrd,converged_r7,converged_altgrd,fit_ok_r7,",
            "fit_ok_altgrd,iterations_r7,iterations_altgrd,trials_r7,trials_altgrd,",
            "r7_noisy_sd_xco2,r7_noisy_sd_psurf")
        for r in rows
            @printf(io, "%03d,%02d,%s,%s,%s,%.6f,%.6f,%.6f,%.6f,%.4f,%.4f,%.4f,%.5f,%.5f,%.4f,%.4f,%s,%s,%d,%d,%d,%d,%d,%d,%d,%d,%.5f,%.4f\n",
                r.scene, r.member, r.surface, r.aerosol, r.sif, r.truth_xco2, r.xco2_r7, r.xco2_ag,
                r.xco2_ag - r.xco2_r7, r.psurf_r7, r.psurf_ag, r.psurf_ag - r.psurf_r7,
                r.sigma_xco2_r7, r.sigma_xco2_ag, r.sigma_psurf_r7, r.sigma_psurf_ag,
                join((@sprintf("%.4f", c) for c in r.chi2_r7), ","),
                join((@sprintf("%.4f", c) for c in r.chi2_ag), ","),
                r.conv_r7, r.conv_ag, r.fit_r7, r.fit_ag, r.iter_r7, r.iter_ag,
                r.trials_r7, r.trials_ag, r.scatter_xco2, r.scatter_psurf)
        end
    end
    serialize(csv_path * ".jls", rows)
    open(summary_path, "w") do io
        write_summary(io, rows, jacobian_paths)
    end
end

stat(v) = isempty(v) ? (NaN, NaN, NaN) : (mean(v), length(v) > 1 ? std(v) : NaN, maximum(abs, v))

function write_summary(io, rows, jacobian_paths)
    # Coverage is reported, never assumed: missing files would change the
    # populations behind every statistic below.
    println(io, "pairs compared=", length(rows), "/160; solves whose Round 7 scene has all ten noisy members=",
        count(r -> r.n_noisy == 10, rows), "/", length(rows))
    println(io, "rows=", length(rows), " max|XCO2 recheck| r7=",
        maximum(abs, [r.recheck_r7 for r in rows]; init=0.0), " altgrd=",
        maximum(abs, [r.recheck_ag for r in rows]; init=0.0))
    println(io, "\n| subset | n | ΔXCO2 mean | SD | max\\|Δ\\| | Δpsurf mean | SD | max\\|Δ\\| | max \\|ΔXCO2\\|/σ_XCO2 | mean \\|ΔXCO2\\|/σ_XCO2 | max \\|ΔXCO2\\|/SD₁₋₁₀ | mean \\|ΔXCO2\\|/SD₁₋₁₀ |")
    println(io, "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
    subsets = [("member 11 (noiseless)", r -> r.member == 11), ("member 01 (noisy)", r -> r.member == 1),
               ("all", r -> true),
               ("no aerosol", r -> r.aerosol == :none), ("aerosol", r -> r.aerosol != :none),
               ("SIF off", r -> r.sif == :off), ("SIF on", r -> r.sif != :off)]
    for (label, keep) in subsets
        sel = filter(keep, rows)
        isempty(sel) && continue
        dx = [r.xco2_ag - r.xco2_r7 for r in sel]
        dp = [r.psurf_ag - r.psurf_r7 for r in sel]
        zs = abs.(dx) ./ [r.sigma_xco2_r7 for r in sel]
        zn = abs.(dx) ./ [r.scatter_xco2 for r in sel]
        a, b, c = stat(dx); d, e, f = stat(dp)
        @printf(io, "| %s | %d | %+.4f | %.4f | %.4f | %+.3f | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f |\n",
            label, length(sel), a, b, c, d, e, f, maximum(zs), mean(zs), maximum(zn), mean(zn))
    end
    println(io, "\nTruth errors (retrieved − truth XCO2, ppm):\n")
    println(io, "| subset | R7 mean | R7 SD | altgrd mean | altgrd SD |")
    println(io, "|---|---:|---:|---:|---:|")
    for (label, keep) in subsets[1:2]
        sel = filter(keep, rows)
        isempty(sel) && continue
        e7 = [r.xco2_r7 - r.truth_xco2 for r in sel]; ea = [r.xco2_ag - r.truth_xco2 for r in sel]
        @printf(io, "| %s | %+.4f | %.4f | %+.4f | %.4f |\n", label, mean(e7), std(e7), mean(ea), std(ea))
    end
    println(io, "\nTruth errors in psurf (retrieved − 1000 hPa):\n")
    println(io, "| subset | R7 mean | R7 SD | altgrd mean | altgrd SD |")
    println(io, "|---|---:|---:|---:|---:|")
    for (label, keep) in subsets
        sel = filter(keep, rows)
        isempty(sel) && continue
        e7 = [r.psurf_r7 - 1000 for r in sel]; ea = [r.psurf_ag - 1000 for r in sel]
        @printf(io, "| %s | %+.3f | %.3f | %+.3f | %.3f |\n", label, mean(e7),
            length(e7) > 1 ? std(e7) : NaN, mean(ea), length(ea) > 1 ? std(ea) : NaN)
    end
    println(io, "\nOther state differences (altgrd − R7), aerosol scenes, all members:\n")
    println(io, "| parameter | mean Δ | SD Δ | max \\|Δ\\| |")
    println(io, "|---|---:|---:|---:|")
    aer = filter(r -> r.aerosol != :none, rows)
    for j in (13, 14, 15, 16, 17, 18, 19)
        isempty(aer) && break
        d = [r.state_ag[j] - r.state_r7[j] for r in aer]
        @printf(io, "| %s | %+.4g | %.4g | %.4g |\n", aer[1].names[j], mean(d),
            length(d) > 1 ? std(d) : NaN, maximum(abs, d))
    end
    println(io, "\nConvergence/fit/iterations:\n")
    for (label, keep) in subsets[1:2]
        sel = filter(keep, rows)
        isempty(sel) && continue
        @printf(io, "- %s: converged R7=%d altgrd=%d; fit_ok R7=%d altgrd=%d; both converged+fit_ok=%d/%d; iterations equal=%d, altgrd more=%d, fewer=%d; mean iterations R7=%.2f altgrd=%.2f\n",
            label, sum(r.conv_r7 for r in sel), sum(r.conv_ag for r in sel),
            sum(r.fit_r7 for r in sel), sum(r.fit_ag for r in sel),
            count(r -> r.conv_r7 == r.conv_ag == r.fit_r7 == r.fit_ag == 1, sel), length(sel),
            count(r -> r.iter_ag == r.iter_r7, sel), count(r -> r.iter_ag > r.iter_r7, sel),
            count(r -> r.iter_ag < r.iter_r7, sel), mean(r.iter_r7 for r in sel), mean(r.iter_ag for r in sel))
    end
    println(io, "\nChi² per band (altgrd − R7):\n")
    println(io, "| subset | band | mean Δχ² | max \\|Δχ²\\| | mean χ² R7 | mean χ² altgrd |")
    println(io, "|---|---|---:|---:|---:|---:|")
    for (label, keep) in subsets[1:2], (ib, band) in enumerate(("O2A", "weak CO2", "strong CO2"))
        sel = filter(keep, rows)
        isempty(sel) && continue
        d = [r.chi2_ag[ib] - r.chi2_r7[ib] for r in sel]
        @printf(io, "| %s | %s | %+.4f | %.4f | %.4f | %.4f |\n", label, band, mean(d), maximum(abs, d),
            mean(r.chi2_r7[ib] for r in sel), mean(r.chi2_ag[ib] for r in sel))
    end
    println(io, "\nPosterior σ (altgrd/R7 ratio):\n")
    for (label, keep) in subsets[1:2]
        sel = filter(keep, rows)
        isempty(sel) && continue
        rx = [r.sigma_xco2_ag / r.sigma_xco2_r7 for r in sel]
        rp = [r.sigma_psurf_ag / r.sigma_psurf_r7 for r in sel]
        @printf(io, "- %s: σ_XCO2 R7 mean %.4f ppm, altgrd mean %.4f ppm, ratio %.4f–%.4f; σ_psurf R7 mean %.3f hPa, altgrd mean %.3f hPa, ratio %.4f–%.4f\n",
            label, mean(r.sigma_xco2_r7 for r in sel), mean(r.sigma_xco2_ag for r in sel),
            minimum(rx), maximum(rx), mean(r.sigma_psurf_r7 for r in sel),
            mean(r.sigma_psurf_ag for r in sel), minimum(rp), maximum(rp))
    end
    println(io, "\nAgreement by class (fractions of solves with |ΔXCO2| ≤ 0.5σ and ≤ 1σ, σ = Round 7 posterior):\n")
    println(io, "| class | n | mean ΔXCO2 | max \\|ΔXCO2\\| | ≤0.5σ | ≤1σ | max \\|Δ\\|/σ | outcomes agree |")
    println(io, "|---|---:|---:|---:|---:|---:|---:|---:|")
    for aerosol in (:none, :aerosol), class in (:all, :off, :on, :urban, :rural, :desert, :forest)
        sel = filter(rows) do r
            (aerosol == :none) == (r.aerosol == :none) &&
                (class == :all || (class in (:off, :on) ? (class == :off) == (r.sif == :off) : r.surface == class))
        end
        isempty(sel) && continue
        z = abs.([r.xco2_ag - r.xco2_r7 for r in sel]) ./ [r.sigma_xco2_r7 for r in sel]
        agree = count(r -> r.conv_r7 == r.conv_ag && r.fit_r7 == r.fit_ag, sel)
        @printf(io, "| %s %s | %d | %+.3f | %.3f | %.1f %% | %.1f %% | %.2f | %d/%d |\n",
            aerosol == :none ? "no aerosol" : "aerosol", class, length(sel),
            mean(r.xco2_ag - r.xco2_r7 for r in sel), maximum(abs(r.xco2_ag - r.xco2_r7) for r in sel),
            100count(<=(0.5), z) / length(z), 100count(<=(1), z) / length(z), maximum(z), agree, length(sel))
    end
    println(io, "\nTruth XCO2 bias by class (retrieved − truth, ppm, mean ± SD):\n")
    println(io, "| class | member | Round 7 | altgrd |")
    println(io, "|---|---|---:|---:|")
    for aerosol in (:none, :aerosol), member in (11, 1)
        sel = filter(r -> (aerosol == :none) == (r.aerosol == :none) && r.member == member, rows)
        isempty(sel) && continue
        e7 = [r.xco2_r7 - r.truth_xco2 for r in sel]; ea = [r.xco2_ag - r.truth_xco2 for r in sel]
        @printf(io, "| %s | %02d | %+.3f ± %.3f | %+.3f ± %.3f |\n", aerosol == :none ? "no aerosol" : "aerosol",
            member, mean(e7), std(e7), mean(ea), std(ea))
    end
    println(io, "\nSolves that did not converge under altgrd:\n")
    for r in filter(r -> r.conv_ag != 1, rows)
        @printf(io, "- scene %03d member %02d (%s, %s, SIF %s): outcome=%d, divergences=%d, iterations %d (R7 %d), ΔXCO2 %+.3f ppm at stop, fit_ok=%d\n",
            r.scene, r.member, r.surface, r.aerosol, r.sif, r.outcome_ag, r.divergences_ag,
            r.iter_ag, r.iter_r7, r.xco2_ag - r.xco2_r7, r.fit_ag)
    end
    noisy = filter(r -> r.member == 11, rows)
    @printf(io, "\nRound 7 noise-driven scatter (SD over members 01–10, per scene): XCO2 median %.4f ppm (range %.4f–%.4f); psurf median %.3f hPa\n",
        median(r.scatter_xco2 for r in noisy), minimum(r.scatter_xco2 for r in noisy),
        maximum(r.scatter_xco2 for r in noisy), median(r.scatter_psurf for r in noisy))
    length(jacobian_paths) == 2 || return
    jp, ja = deserialize(jacobian_paths[1]), deserialize(jacobian_paths[2])
    println(io, "\nSurface-pressure Jacobian column at truth (instrument space, per hPa):\n")
    println(io, "| state | band | ‖K₁‖ pressure | ‖K₁‖ altitude | ratio | cos angle | ‖ΔK₁‖/‖K₁ pressure‖ |")
    println(io, "|---|---|---:|---:|---:|---:|---:|")
    ranges = (1:934, 935:1741, 1742:2742)
    for key in ("s011_truth", "s036_truth")
        kp, ka = jp[key].jacobian[:, 1], ja[key].jacobian[:, 1]
        for (band, r) in zip(("O2A", "weak CO2", "strong CO2", "all"), (ranges..., 1:2742))
            @printf(io, "| %s | %s | %.4e | %.4e | %.4f | %.6f | %.4f |\n", key, band, norm(kp[r]),
                norm(ka[r]), norm(ka[r]) / norm(kp[r]), dot(kp[r], ka[r]) / (norm(kp[r]) * norm(ka[r])),
                norm(ka[r] - kp[r]) / norm(kp[r]))
        end
    end
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main(ARGS...)
