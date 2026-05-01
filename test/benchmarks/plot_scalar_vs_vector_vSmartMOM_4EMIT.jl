#=
plot_scalar_vs_vector_vSmartMOM_4EMIT.jl — vSmartMOM intercomparison plotter.

Compares Stokes-I from scalar (Stokes_I) and vector (Stokes_IQU) RT runs of
the EMIT regression for the same atmospheric scenarios. Vector RT is
physically more complete (includes the depolarising effect of multiple
scattering on intensity); the Stokes-I difference between scalar and vector
is the "scalar approximation error."

Reads noRS .dat files from:
  scalar : $HOME/data/EMIT_MODTRANcomp/vSmartMOM_out/scalarI/
  vector : $HOME/data/EMIT_MODTRANcomp/vSmartMOM_out/vectorIQU/

Each .dat has columns: wl_nm wn_cm-1 R_tot T_tot tau_total hem_R hem_T

Output: /home/sanghavi/data/EMIT_MODTRANcomp/plots/scalar_vs_vector/
        — overlay (R, T, τ) and ΔI (vec-sca abs and relative) per scenario,
          plus full-range and per-spectral-slice variants. Filenames carry
          the `_scalar_vs_vector` tag.
        — log-y variants of every linear plot are also written to the
          adjacent `logplots/` subdirectory, with `_log` appended to the
          filename. For diff plots, the log version plots `|Δ|` (signed Δ
          can't be log-scaled) and adjusts the title accordingly.

Sister scripts:
  plot_MODTRAN_vSmartMOM_4EMIT.jl       (scalar  vs MODTRAN)
  plot_MODTRAN_vec_vSmartMOM_4EMIT.jl   (vector  vs MODTRAN, _IQU output tag)
=#
using Plots
using DelimitedFiles
using LaTeXStrings
using Printf
using Statistics

default(alpha = 0.5)

const SCAL_DIR = "/home/sanghavi/data/EMIT_MODTRANcomp/vSmartMOM_out/scalarI"
const VEC_DIR  = "/home/sanghavi/data/EMIT_MODTRANcomp/vSmartMOM_out/vectorIQU"
const OUT_DIR  = "/home/sanghavi/data/EMIT_MODTRANcomp/plots/scalar_vs_vector"
mkpath(OUT_DIR)

# ---------------------------------------------------------------------------
# Scenario set to plot. Adjust here if the regression-run grid changes.
# Values must match files actually present in BOTH subdirs.
const H2O_SET    = (0.05, 2.05)
const AOT_SET    = (0.0001,)
const GNDALT_SET = (0.0, 2.0, 5.0)
const TSZ_SET    = (168.0, 178.0)

# Spectral slices for zoomed comparison plots.
const SLICES = (
    ("uv",    350.0,  400.0),
    ("vis",   500.0,  700.0),
    ("nir",   950.0, 1100.0),
    ("swir1",1500.0, 1700.0),
    ("swir2",2100.0, 2200.0),
)

scen_tag(h2o, aot, gndalt, tsz) =
    @sprintf("H2O%.4f_AOT%.4f_GNDALT%.3f_TSZ%.1f", h2o, aot, gndalt, tsz)

dat_path(dir, h2o, aot, gndalt, tsz) =
    joinpath(dir, "vSmartMOM_$(scen_tag(h2o, aot, gndalt, tsz)).dat")

function load_dat(path)
    isfile(path) || return nothing
    d = readdlm(path, comments=true, comment_char='#')
    (λ = d[:, 1], R = d[:, 3], T = d[:, 4], τ = d[:, 5],
     hem_R = d[:, 6], hem_T = d[:, 7])
end

"""
    log_path(linear_path)

Translate `<dir>/foo.pdf` → `<dir>/logplots/foo_log.pdf`.
"""
function log_path(linear_path::AbstractString)
    dir = dirname(linear_path)
    isempty(dir) && (dir = ".")
    log_dir = joinpath(dir, "logplots")
    mkpath(log_dir)
    stem, ext = splitext(basename(linear_path))
    return joinpath(log_dir, stem * "_log" * ext)
end

function _build_overlay(s, v, tag, slice_label, m, yscale_RT)
    # yscale_RT applies to the R and T panels; τ panel always log-y (already
    # spans many orders of magnitude on a linear plot is unhelpful).
    pR = plot(s.λ[m], max.(s.R[m], yscale_RT === :log10 ? 1e-30 : -Inf),
              label="scalar I", color=:blue, lw=0.8,
              ylabel=L"$R_{TOA}$", legend=:topright, yscale=yscale_RT,
              title="R_TOA  $(tag)  [$slice_label]")
    plot!(pR, v.λ[m], max.(v.R[m], yscale_RT === :log10 ? 1e-30 : -Inf),
          label="vector  I (IQU)", color=:red, lw=0.8)

    pT = plot(s.λ[m], max.(s.T[m], yscale_RT === :log10 ? 1e-30 : -Inf),
              label="scalar I", color=:blue, lw=0.8,
              ylabel=L"$T_{BOA}$", legend=:topright, yscale=yscale_RT,
              title="T_BOA  $(tag)  [$slice_label]")
    plot!(pT, v.λ[m], max.(v.T[m], yscale_RT === :log10 ? 1e-30 : -Inf),
          label="vector  I (IQU)", color=:red, lw=0.8)

    pτ = plot(s.λ[m], max.(s.τ[m], 1e-30), label="scalar τ", color=:blue, lw=0.8,
              ylabel=L"$\tau_\mathrm{total}$", legend=:topright,
              yscale=:log10, ylims=(1e-4, 1e2),
              title="τ_total  $(tag)  [$slice_label]")
    plot!(pτ, v.λ[m], max.(v.τ[m], 1e-30), label="vector τ", color=:red, lw=0.8)

    plot(pR, pT, pτ, layout=(3, 1), size=(1400, 900),
         xlabel="Wavelength [nm]")
end

function plot_overlay(s, v, tag, slice_label, lo, hi)
    m = (s.λ .>= lo) .& (s.λ .<= hi)
    out = joinpath(OUT_DIR,
                   "overlay_$(tag)_$(slice_label)_scalar_vs_vector.pdf")
    pl_lin = _build_overlay(s, v, tag, slice_label, m, :identity)
    savefig(pl_lin, out)
    pl_log = _build_overlay(s, v, tag, slice_label, m, :log10)
    savefig(pl_log, log_path(out))
    return out
end

function _build_diff(s, v, tag, slice_label, m, yscale_diff)
    ΔR = v.R[m] .- s.R[m]
    ΔT = v.T[m] .- s.T[m]
    rel_R = ΔR ./ max.(abs.(s.R[m]), eps())
    rel_T = ΔT ./ max.(abs.(s.T[m]), eps())

    # On log-y, plot |Δ| since Δ is signed.
    yvals_R   = yscale_diff === :log10 ? max.(abs.(ΔR), 1e-30) : ΔR
    yvals_T   = yscale_diff === :log10 ? max.(abs.(ΔT), 1e-30) : ΔT
    yvals_rR  = yscale_diff === :log10 ? max.(abs.(rel_R), 1e-30) : rel_R
    yvals_rT  = yscale_diff === :log10 ? max.(abs.(rel_T), 1e-30) : rel_T
    abs_lbl   = yscale_diff === :log10 ? "|" : ""
    title_a   = yscale_diff === :log10 ?
                "|Stokes-I difference|  $(tag)  [$slice_label]" :
                "Stokes-I difference (vector − scalar)  $(tag)  [$slice_label]"
    title_b   = yscale_diff === :log10 ?
                "|Relative Stokes-I difference|  $(tag)  [$slice_label]" :
                "Relative Stokes-I difference  $(tag)  [$slice_label]"

    pa = plot(s.λ[m], yvals_R, label="$(abs_lbl)ΔR$(abs_lbl)",
              color=:purple, lw=0.7, ylabel="ΔR_TOA", yscale=yscale_diff,
              title=title_a)
    plot!(pa, s.λ[m], yvals_T, label="$(abs_lbl)ΔT$(abs_lbl)",
          color=:darkgreen, lw=0.7)

    pb = plot(s.λ[m], yvals_rR, label="$(abs_lbl)ΔR / R_sca$(abs_lbl)",
              color=:purple, lw=0.7, ylabel="relative ΔI", yscale=yscale_diff,
              title=title_b)
    plot!(pb, s.λ[m], yvals_rT, label="$(abs_lbl)ΔT / T_sca$(abs_lbl)",
          color=:darkgreen, lw=0.7)

    plot(pa, pb, layout=(2, 1), size=(1400, 700),
         xlabel="Wavelength [nm]")
end

function plot_diff(s, v, tag, slice_label, lo, hi)
    m = (s.λ .>= lo) .& (s.λ .<= hi)
    out = joinpath(OUT_DIR,
                   "diff_$(tag)_$(slice_label)_scalar_vs_vector.pdf")
    pl_lin = _build_diff(s, v, tag, slice_label, m, :identity)
    savefig(pl_lin, out)
    pl_log = _build_diff(s, v, tag, slice_label, m, :log10)
    savefig(pl_log, log_path(out))
    return out
end

function summarise(s, v, tag)
    function stats(x, name)
        x = filter(isfinite, x)
        @sprintf("  %-22s med=%9.3e  mean=%9.3e  max=%9.3e",
                 name, median(x), mean(x), maximum(x))
    end
    ΔR = v.R .- s.R
    ΔT = v.T .- s.T
    rel_R = abs.(ΔR) ./ max.(abs.(s.R), eps())
    rel_T = abs.(ΔT) ./ max.(abs.(s.T), eps())
    println("=== ", tag, " ===")
    println(stats(abs.(ΔR), "|ΔR| absolute"))
    println(stats(abs.(ΔT), "|ΔT| absolute"))
    println(stats(rel_R,    "|ΔR|/R relative"))
    println(stats(rel_T,    "|ΔT|/T relative"))
end

# ---------------------------------------------------------------------------
function main()
    n_scenarios = 0
    n_missing   = 0
    for h2o in H2O_SET, aot in AOT_SET, gndalt in GNDALT_SET, tsz in TSZ_SET
        tag = scen_tag(h2o, aot, gndalt, tsz)
        s_path = dat_path(SCAL_DIR, h2o, aot, gndalt, tsz)
        v_path = dat_path(VEC_DIR,  h2o, aot, gndalt, tsz)

        s = load_dat(s_path)
        v = load_dat(v_path)
        if s === nothing || v === nothing
            @warn "missing one side — skipping" tag scalar_present=(s!==nothing) vector_present=(v!==nothing)
            n_missing += 1
            continue
        end
        @assert length(s.λ) == length(v.λ) "λ-grid mismatch in $(tag)"
        n_scenarios += 1

        # Full-range overlay + diff
        plot_overlay(s, v, tag, "full",  minimum(s.λ), maximum(s.λ))
        plot_diff(s, v, tag, "full",  minimum(s.λ), maximum(s.λ))
        # Per-slice
        for (slice, lo, hi) in SLICES
            plot_overlay(s, v, tag, slice, lo, hi)
            plot_diff(s, v, tag, slice, lo, hi)
        end
        summarise(s, v, tag)
    end
    println()
    println("Scenarios plotted: ", n_scenarios)
    println("Scenarios skipped (file missing): ", n_missing)
    println("Output: ", OUT_DIR)
end

main()
