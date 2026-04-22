# =============================================================================
# compare_rt_EMIT.jl — compare vSmartMOM EMIT .dat outputs against MODTRAN
# (Phase 6 port from sanghavi).
#
# Reads per-scenario .dat spectra produced by emit_modtran_noRS_scenarios.jl
# and compares them against MODTRAN simulations stored in a NetCDF file.
#
# Paths are configurable via environment variables (fall back to the sanghavi
# defaults):
#   EMIT_DAT_DIR       — directory of vSmartMOM .dat outputs
#                        (default: ~/data/EMIT_MODTRANcomp/vSmartMOM_out)
#   EMIT_MODTRAN_NC    — path to MODTRAN NetCDF LUT
#                        (default: ~/data/EMIT_MODTRANcomp/emit_20250404.nc)
#   EMIT_OUT_DIR       — output dir for comparison CSVs + per-scenario .dat
#                        (default: EMIT_DAT_DIR/comparison)
#   EMIT_WRITE_PY      — 1 to write a matplotlib helper script; 0 (default)
#                        skips the external-Python dependency.
#
# Returns `nothing` if EMIT_DAT_DIR or EMIT_MODTRAN_NC does not exist — the
# function is idempotent when the data bundle is absent from disk.
# =============================================================================

using NCDatasets
using Interpolations
using Printf
using Dates
using DelimitedFiles
using Statistics

const EMIT_DAT_DIR    = get(ENV, "EMIT_DAT_DIR",
                             expanduser("~/data/EMIT_MODTRANcomp/vSmartMOM_out"))
const EMIT_MODTRAN_NC = get(ENV, "EMIT_MODTRAN_NC",
                             expanduser("~/data/EMIT_MODTRANcomp/emit_20250404.nc"))
const EMIT_OUT_DIR    = get(ENV, "EMIT_OUT_DIR",
                             joinpath(EMIT_DAT_DIR, "comparison"))

"""
    read_emit_dat(path)

Parse a vSmartMOM EMIT .dat file, returning scenario coordinates and spectra.

Expected filename pattern: `...H2O{val}_AOT{val}_GNDALT{val}_TSZ{val}.dat`.
Column layout (2-line header, then data): wl_nm, wn, R_tot, T_tot, R_ss, T_ss.
"""
function read_emit_dat(path)
    fname = basename(path)
    m = match(r"H2O([\d.]+)_AOT([\d.]+)_GNDALT([\d.]+)_TSZ([\d.]+)\.dat$", fname)
    isnothing(m) && error("Cannot parse filename: $fname")
    h2o    = parse(Float64, m[1])
    aot    = parse(Float64, m[2])
    gndalt = parse(Float64, m[3])
    tsz    = parse(Float64, m[4])

    data = readdlm(path, skipstart = 2)
    return (h2o = h2o, aot = aot, gndalt = gndalt, tsz = tsz,
            wl_nm = data[:, 1], wn = data[:, 2],
            R_tot = data[:, 3], T_tot = data[:, 4],
            R_ss  = data[:, 5], T_ss  = data[:, 6])
end

"""
    load_modtran(nc_path)

Load the MODTRAN LUT NetCDF into a NamedTuple.
"""
function load_modtran(nc_path)
    ds = NCDataset(nc_path)
    result = (
        wl_nm     = Float64.(nomissing(Array(ds["wl"]),           NaN)),
        H2OSTR    = Float64.(nomissing(Array(ds["H2OSTR"]),       NaN)),
        AOT550    = Float64.(nomissing(Array(ds["AOT550"]),       NaN)),
        obs_zen   = Float64.(nomissing(Array(ds["observer_zenith"]), NaN)),
        elev_km   = Float64.(nomissing(Array(ds["surface_elevation_km"]), NaN)),
        solar_irr = Float64.(nomissing(Array(ds["solar_irr"]),    NaN)),
        solzen    = Float64(coalesce(ds["solzen"][], NaN)),
        rhoatm    = Float64.(nomissing(Array(ds["rhoatm"]),       NaN)),
        transm_down_dir = Float64.(nomissing(Array(ds["transm_down_dir"]), NaN)),
        transm_down_dif = Float64.(nomissing(Array(ds["transm_down_dif"]), NaN)),
        transm_up_dir   = Float64.(nomissing(Array(ds["transm_up_dir"]),   NaN)),
        transm_up_dif   = Float64.(nomissing(Array(ds["transm_up_dif"]),   NaN)),
    )
    close(ds)
    return result
end

closest_idx(arr, val) = argmin(abs.(arr .- val))

"""
    match_modtran_indices(mod, h2o, aot, gndalt, tsz)

Return (iH, iA, iZ, iE) into the MODTRAN arrays for the scenario matching
the vSmartMOM coordinates (closest-gridpoint).
"""
function match_modtran_indices(mod, h2o, aot, gndalt, tsz)
    iH = closest_idx(mod.H2OSTR,  h2o)
    iA = closest_idx(mod.AOT550,  aot)
    iZ = closest_idx(mod.obs_zen, tsz)      # TSZ ↔ observer_zenith
    iE = closest_idx(mod.elev_km, gndalt)   # GNDALT ↔ surface_elevation_km
    return iH, iA, iZ, iE
end

"""
    interp_to_grid(src_wl, src_spec, dst_wl)

Linearly interpolate `src_spec` onto `dst_wl`. NaN where out-of-range.
"""
function interp_to_grid(src_wl, src_spec, dst_wl)
    if src_wl[1] > src_wl[end]
        src_wl   = reverse(src_wl)
        src_spec = reverse(src_spec)
    end
    itp = LinearInterpolation(src_wl, src_spec; extrapolation_bc = NaN)
    return itp.(dst_wl)
end

"""
    spectral_stats(vsm, mod_interp)

Return (n, rmse, mae, max_abs, mean_rel, max_rel) over valid (non-NaN) points.
"""
function spectral_stats(vsm, mod_interp)
    mask = .!isnan.(mod_interp) .& .!isnan.(vsm) .& (mod_interp .!= 0)
    n = count(mask)
    n == 0 && return (n = 0, rmse = NaN, mae = NaN, max_abs = NaN,
                      mean_rel = NaN, max_rel = NaN)
    diff    = vsm[mask] .- mod_interp[mask]
    reldiff = diff ./ abs.(mod_interp[mask])
    return (n       = n,
            rmse    = sqrt(mean(diff .^ 2)),
            mae     = mean(abs.(diff)),
            max_abs = maximum(abs.(diff)),
            mean_rel = mean(reldiff),
            max_rel  = maximum(abs.(reldiff)))
end

"""
    write_comparison_dat(path, wl, vsm, mod_interp; header="")

Write a four-column comparison file (wl, vSmartMOM, MODTRAN, rel_diff) for
external plotting.
"""
function write_comparison_dat(path, wl, vsm, mod_interp; header = "")
    open(path, "w") do io
        println(io, "# ", header)
        @printf(io, "# %16s %16s %16s %16s\n",
                "wl_nm", "vSmartMOM", "MODTRAN", "rel_diff")
        for i in eachindex(wl)
            rd = (mod_interp[i] != 0 && !isnan(mod_interp[i])) ?
                 (vsm[i] - mod_interp[i]) / abs(mod_interp[i]) : NaN
            @printf(io, "  %16.6f %16.8e %16.8e %16.8e\n",
                    wl[i], vsm[i], mod_interp[i], rd)
        end
    end
end

"""
    run_comparison(; dat_dir=EMIT_DAT_DIR, modtran_nc=EMIT_MODTRAN_NC, out_dir=EMIT_OUT_DIR)

Main entry point. Returns the stats_rows vector or `nothing` if inputs
are not available.
"""
function run_comparison(; dat_dir     = EMIT_DAT_DIR,
                          modtran_nc  = EMIT_MODTRAN_NC,
                          out_dir     = EMIT_OUT_DIR)
    isdir(dat_dir) || (@info "compare_rt_EMIT: dat_dir not found — skipping" dat_dir; return nothing)
    isfile(modtran_nc) || (@info "compare_rt_EMIT: MODTRAN NC not found — skipping" modtran_nc; return nothing)

    mkpath(out_dir)

    @info "Loading MODTRAN LUT" modtran_nc
    mod = load_modtran(modtran_nc)
    @info "MODTRAN grid" H2OSTR = mod.H2OSTR AOT550 = mod.AOT550 obs_zen = mod.obs_zen elev_km = mod.elev_km wl_range = extrema(mod.wl_nm)

    dat_files = filter(f -> endswith(f, ".dat"), readdir(dat_dir; join = true))
    isempty(dat_files) && (@info "compare_rt_EMIT: no .dat files found — skipping" dat_dir; return nothing)
    @info "Found .dat files" n = length(dat_files)

    stats_rows = NamedTuple[]

    for (i, fpath) in enumerate(sort(dat_files))
        sc = read_emit_dat(fpath)

        iH, iA, iZ, iE = match_modtran_indices(mod, sc.h2o, sc.aot, sc.gndalt, sc.tsz)
        mod_h2o = mod.H2OSTR[iH]; mod_aot  = mod.AOT550[iA]
        mod_zen = mod.obs_zen[iZ]; mod_elev = mod.elev_km[iE]

        mod_rhoatm        = mod.rhoatm[iH, iA, iZ, iE, :]
        mod_rhoatm_interp = interp_to_grid(mod.wl_nm, mod_rhoatm, sc.wl_nm)
        mod_tdown_dir = interp_to_grid(mod.wl_nm,
                            mod.transm_down_dir[iH, iA, iZ, iE, :], sc.wl_nm)
        mod_tdown_dif = interp_to_grid(mod.wl_nm,
                            mod.transm_down_dif[iH, iA, iZ, iE, :], sc.wl_nm)
        mod_tdown = mod_tdown_dir .+ mod_tdown_dif

        st = spectral_stats(sc.R_tot, mod_rhoatm_interp)

        label = @sprintf("H2O=%.2f AOT=%.4f GNDALT=%.1f TSZ=%.0f",
                         sc.h2o, sc.aot, sc.gndalt, sc.tsz)
        @info "[$i/$(length(dat_files))] $label" matched_to = "H2O=$(mod_h2o) AOT=$(mod_aot) zen=$(mod_zen) elev=$(mod_elev)" n_pts = st.n mean_rel = round(st.mean_rel; sigdigits = 3) max_rel = round(st.max_rel; sigdigits = 3)

        push!(stats_rows, (
            h2o = sc.h2o, aot = sc.aot, gndalt = sc.gndalt, tsz = sc.tsz,
            mod_h2o = mod_h2o, mod_aot = mod_aot, mod_zen = mod_zen,
            mod_elev = mod_elev,
            n_pts = st.n, rmse = st.rmse, mae = st.mae,
            max_abs = st.max_abs, mean_rel = st.mean_rel, max_rel = st.max_rel,
        ))

        tag = @sprintf("H2O%.4f_AOT%.4f_GNDALT%.3f_TSZ%.1f",
                       sc.h2o, sc.aot, sc.gndalt, sc.tsz)

        write_comparison_dat(
            joinpath(out_dir, "cmp_rhoatm_$(tag).dat"),
            sc.wl_nm, sc.R_tot, mod_rhoatm_interp;
            header = "rhoatm comparison — $label  |  MODTRAN matched: H2O=$(mod_h2o) AOT=$(mod_aot) zen=$(mod_zen) elev=$(mod_elev)")
        write_comparison_dat(
            joinpath(out_dir, "cmp_tdown_$(tag).dat"),
            sc.wl_nm, sc.T_tot, mod_tdown;
            header = "T_down comparison — $label  |  MODTRAN matched: H2O=$(mod_h2o) AOT=$(mod_aot) zen=$(mod_zen) elev=$(mod_elev)")
    end

    csv_path = joinpath(out_dir, "comparison_summary.csv")
    open(csv_path, "w") do io
        println(io, "h2o,aot,gndalt,tsz,mod_h2o,mod_aot,mod_zen,mod_elev,n_pts,rmse,mae,max_abs,mean_rel,max_rel")
        for r in stats_rows
            @printf(io, "%.4f,%.4f,%.3f,%.1f,%.4f,%.4f,%.1f,%.1f,%d,%.6e,%.6e,%.6e,%.6e,%.6e\n",
                    r.h2o, r.aot, r.gndalt, r.tsz,
                    r.mod_h2o, r.mod_aot, r.mod_zen, r.mod_elev,
                    r.n_pts, r.rmse, r.mae, r.max_abs, r.mean_rel, r.max_rel)
        end
    end
    @info "Summary CSV written" csv_path

    # Optional: write Python/matplotlib helper script (external dep)
    if get(ENV, "EMIT_WRITE_PY", "0") == "1"
        plot_script = joinpath(out_dir, "_plot_comparison.py")
        write_plot_script(plot_script, out_dir)
        @info "Wrote matplotlib helper — run it externally" plot_script
    end

    @info "Comparison complete" out_dir
    return stats_rows
end

function write_plot_script(path, out_dir)
    open(path, "w") do io
        print(io, raw"""
import os, glob, re, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = sys.argv[1] if len(sys.argv) > 1 else """ * "\"$(escape_string(out_dir))\"" * raw"""

rhoatm_files = sorted(glob.glob(os.path.join(OUT, "cmp_rhoatm_*.dat")))
print(f"Plotting {len(rhoatm_files)} rhoatm comparison files")

for fpath in rhoatm_files:
    tag = os.path.basename(fpath).replace("cmp_rhoatm_", "").replace(".dat", "")
    data = np.loadtxt(fpath)
    wl, vsm, mod = data[:, 0], data[:, 1], data[:, 2]
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 7), sharex=True,
                                    gridspec_kw={"height_ratios": [3, 1]})
    ax1.plot(wl, mod, linewidth=0.3, alpha=0.8, label="MODTRAN rhoatm")
    ax1.plot(wl, vsm, linewidth=0.3, alpha=0.8, label="vSmartMOM R_tot")
    ax1.set_ylabel("Reflectance"); ax1.set_title(tag.replace("_", "  "), fontsize=10)
    ax1.legend(fontsize=8); ax1.set_xlim(wl.min(), wl.max())
    with np.errstate(divide="ignore", invalid="ignore"):
        rel = np.where(np.abs(mod) > 0, (vsm - mod) / np.abs(mod), np.nan)
    ax2.plot(wl, rel, linewidth=0.3, color="C2"); ax2.axhline(0, color="k", linewidth=0.5)
    ax2.set_ylabel("Rel. diff"); ax2.set_xlabel("Wavelength (nm)"); ax2.set_ylim(-1, 1)
    plt.tight_layout()
    plt.savefig(os.path.join(OUT, f"cmp_rhoatm_{tag}.png"), dpi=150)
    plt.close(fig)

csv_path = os.path.join(OUT, "comparison_summary.csv")
if os.path.isfile(csv_path):
    csv = np.genfromtxt(csv_path, delimiter=",", skip_header=1)
    mean_rel = csv[:, 12]; max_rel = csv[:, 13]; idx = np.arange(len(mean_rel))
    fig, ax = plt.subplots(figsize=(12, 5))
    ax.bar(idx - 0.15, mean_rel, width=0.3, label="mean rel. diff", alpha=0.8)
    ax.bar(idx + 0.15, max_rel,  width=0.3, label="max |rel. diff|", alpha=0.8)
    ax.axhline(0, color="k", linewidth=0.5)
    ax.set_xlabel("Scenario index"); ax.set_ylabel("Relative difference (rhoatm)")
    ax.set_title("vSmartMOM vs MODTRAN — rhoatm comparison")
    ax.legend(); plt.tight_layout()
    plt.savefig(os.path.join(OUT, "aggregate_rhoatm_comparison.png"), dpi=150)
    plt.close(fig)
print("Done")
""")
    end
end

# Script entry point
if abspath(PROGRAM_FILE) == @__FILE__
    stats = run_comparison()
    if isnothing(stats)
        @info "compare_rt_EMIT: inputs not available; skipped."
    end
end
