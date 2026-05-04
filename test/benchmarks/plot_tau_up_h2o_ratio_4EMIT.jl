#=
plot_tau_up_h2o_ratio_4EMIT.jl — vertical τ from upwelling-direct transmittance,
MODTRAN vs vSmartMOM at H2OSTR=0.05 and 2.05, plus the τ(2.05)/τ(0.05) ratio.

Vertical τ is recovered from the directly-transmitted upwelling beam via
    τ_up_vert(λ) = -μ_v · ln(tu_dir(λ))     with μ_v = cos(VZA)
matching the convention used by `plot_MODTRAN_vSmartMOM_4EMIT.jl`.

The ratio τ(H2OSTR=2.05) / τ(H2OSTR=0.05) is a clean diagnostic of the H2O
column response: a constant ratio across the full spectrum would imply the
reported τ scales linearly with the H2OSTR multiplier (i.e. file-level scaling
of all gases together — see HANDOFF / project_emit_modtran_swir2_gap.md);
H2O-band morphology with ratios → 1 in the windows is the physically correct
behaviour. The naive 41 = 2.05/0.05 reference line marks "all gases scaled".

Reads:
  MODTRAN  $HOME/data/EMIT_MODTRANcomp/MODTRAN_out/MODTRAN_*.dat
           cols: λ_nm wn_cm-1 rho_atm td_dir td_dif tu_dir tu_dif sphalb
  vSmartMOM $HOME/data/EMIT_MODTRANcomp/vSmartMOM_out/scalar_modtran_equiv/
           same six MODTRAN-equivalent fields.

Output:
  $HOME/data/EMIT_MODTRANcomp/plots/h2o_tau_up_ratio/
     overlay_H2O_<scen>.pdf            both codes, both H2OSTR levels
     tau_up_H2O0.05_<scen>.pdf         MODTRAN vs vSmartMOM at H2OSTR=0.05
     tau_up_H2O2.05_<scen>.pdf         MODTRAN vs vSmartMOM at H2OSTR=2.05
     ratio_H2O_<scen>.pdf              ratio τ(2.05)/τ(0.05) for both codes
  ...plus _<slice>.pdf variants for uv / vis / nir / swir1 / swir2 zooms,
  and adjacent logplots/<…>_log.pdf log-y twins for every linear plot.

Sister scripts:
  plot_MODTRAN_vSmartMOM_4EMIT.jl       (scalar-I, full RT field comparison)
  plot_MODTRAN_vec_vSmartMOM_4EMIT.jl   (vector-IQU, full RT field comparison)
  plot_scalar_vs_vector_vSmartMOM_4EMIT.jl (scalar vs vector, vSmartMOM only)
=#
using Plots
using DelimitedFiles
using LaTeXStrings
using Printf
using Plots.PlotMeasures: mm

default(alpha = 0.5)

const MOD_DIR = "/home/sanghavi/data/EMIT_MODTRANcomp/MODTRAN_out"
const MOM_DIR = "/home/sanghavi/data/EMIT_MODTRANcomp/vSmartMOM_out/scalar_modtran_equiv"
const OUT_DIR = "/home/sanghavi/data/EMIT_MODTRANcomp/plots/h2o_tau_up_ratio"
mkpath(OUT_DIR)

# Geometry — matches the rest of the EMIT regression suite.
const SZA = 39.583965990615404
const VZA = 12.0
const μ_v = cosd(VZA)

# Single scenario set — adjust here if more (AOT, GNDALT, TSZ) combinations
# are needed. The H2O sensitivity is between the two H2OSTR levels below.
const H2O_LO = 0.05
const H2O_HI = 2.05
const AOT    = 0.0001
const GNDALT = 0.0
const TSZ    = 168.0

const SLICES = (
    ("full",  -Inf,    Inf),
    ("uv",     350.0,  400.0),
    ("vis",    500.0,  700.0),
    ("nir",    950.0, 1100.0),
    ("swir1", 1500.0, 1700.0),
    ("swir2", 2100.0, 2200.0),
)

# Wide layout with explicit margins — the LaTeXStrings y-labels and the
# (1400, 500) canvas otherwise clip on the left, which is exactly the
# defect we are correcting from the prior PNG output.
const FIG_SIZE   = (1400, 500)
const FIG_MARGIN = (left=14mm, bottom=10mm, right=8mm, top=6mm)

scen_tag(h2o, aot, gndalt, tsz) =
    @sprintf("H2O%.4f_AOT%.4f_GNDALT%.3f_TSZ%.1f", h2o, aot, gndalt, tsz)

mod_path(h2o)  = joinpath(MOD_DIR, "MODTRAN_$(scen_tag(h2o, AOT, GNDALT, TSZ)).dat")
vsm_path(h2o)  = joinpath(MOM_DIR, "vSmartMOM_modtran_equiv_$(scen_tag(h2o, AOT, GNDALT, TSZ)).dat")

# MODTRAN .dat has two '#'-prefixed header lines for some sets and a single
# `skipstart=2` line for others; using `comments` is robust to both.
load_dat(p) = readdlm(p, comments=true, comment_char='#')

"""
    tau_up_vertical(tu_dir, μ_v)

Vertical optical depth recovered from the directly-transmitted upwelling
beam: τ = -μ_v · ln(tu_dir). Negative or zero `tu_dir` (numerical noise
in saturated bands) is floored at 1e-30 so the log is well-defined.
"""
tau_up_vertical(tu_dir, μ) = -μ .* log.(max.(tu_dir, 1e-30))

# === robust savefig (catches GR's silent 0-byte failure mode) ============
"""
    _savefig_checked(name)

`Plots.savefig` opens the destination via `fopen("wb")` (creates a 0-byte
file) before GR renders and writes the content. If GR fails silently — most
often because the GKS workstation is in a wedged state — the file is left
at 0 bytes with no exception thrown. We detect that, retry once after
`closeall`, and remove the file if it stays empty so callers can see the
failure instead of inheriting a corrupt artifact.
"""
function _savefig_checked(plt, name::AbstractString)
    try
        Plots.savefig(plt, name)
    catch err
        @warn "savefig: threw exception" name exception=err
        try; Plots.closeall(); catch; end
        isfile(name) && filesize(name) == 0 && rm(name; force=true)
        return false
    end
    if isfile(name) && filesize(name) == 0
        @warn "savefig: produced 0-byte file, retrying after closeall" name
        try; Plots.closeall(); catch; end
        try
            Plots.savefig(plt, name)
        catch err
            @warn "savefig: retry threw" name exception=err
        end
    end
    if isfile(name) && filesize(name) == 0
        @warn "savefig: file still 0 bytes after retry — removing" name
        rm(name; force=true)
        return false
    end
    return true
end

"""
    save_both(plt, name)

Save `plt` to `name` (linear y-axis) and to `<dir>/logplots/<stem>_log<ext>`
with `yscale=:log10`. Restores the y-axis to linear afterwards so the
caller's plot is not mutated for any subsequent operations.

Wrapped in try/catch: switching to `:log10` over data containing
non-positive values (signed differences, empty slices) can collapse the
plot layout and trigger `AssertionError: total_plotarea_vertical > 0mm`
on PlotUtils/Plots backends. We skip the affected save rather than
aborting the whole script.
"""
function save_both(plt, name::AbstractString)
    dir = dirname(name)
    isempty(dir) && (dir = ".")
    mkpath(dir)
    saved_linear = _savefig_checked(plt, name)
    log_dir = joinpath(dir, "logplots")
    mkpath(log_dir)
    stem, ext = splitext(basename(name))
    saved_linear || return
    log_name = joinpath(log_dir, stem * "_log" * ext)
    try
        yaxis!(plt, :log10)
        _savefig_checked(plt, log_name)
    catch err
        @warn "save_both: log-axis save failed" name exception=err
        try; Plots.closeall(); catch; end
        isfile(log_name) && filesize(log_name) == 0 && rm(log_name; force=true)
    finally
        try; yaxis!(plt, :identity); catch; end
    end
end

# === plotting helpers =====================================================
const SCEN_BASE = scen_tag(H2O_LO, AOT, GNDALT, TSZ)  # naming anchor; H2O field gets stripped per call

function _slice_mask(λ, lo, hi)
    isfinite(lo) && isfinite(hi) || return trues(length(λ))
    return (λ .>= lo) .& (λ .<= hi)
end

# Suffix used in output filenames to identify the (AOT, GNDALT, TSZ) scenario
# without redundantly carrying the H2O label that varies between panels.
const SCEN_SUFFIX = @sprintf("AOT%.4f_GNDALT%.3f_TSZ%.1f", AOT, GNDALT, TSZ)

scen_tag_short() = SCEN_SUFFIX

"""
    plot_tau_one_h2o(mod_d, vsm_d, h2o_label, slice_label, lo, hi)

Compare MODTRAN vs vSmartMOM vertical τ_up for a single H2OSTR level over
the wavelength window [lo, hi]. Saves linear + log via `save_both`.
"""
function plot_tau_one_h2o(mod_d, vsm_d, h2o::Real, slice_label, lo, hi)
    m_mod = _slice_mask(mod_d[:,1], lo, hi)
    m_vsm = _slice_mask(vsm_d[:,1], lo, hi)
    if !any(m_mod) || !any(m_vsm)
        @info "plot_tau_one_h2o: empty slice" h2o slice_label
        return
    end
    τ_mod = tau_up_vertical(mod_d[m_mod, 6], μ_v)
    τ_vsm = tau_up_vertical(vsm_d[m_vsm, 6], μ_v)

    plt = plot(mod_d[m_mod, 1], τ_mod,
               label=@sprintf("MODTRAN, H₂O=%.2f", h2o),
               color=:blue, lw=0.9,
               xlabel="Wavelength [nm]",
               ylabel=L"$\tau_\mathrm{atm}$ [upward]",
               title=@sprintf("τ_up  H₂O=%.2f  [%s]  (%s)",
                              h2o, slice_label, scen_tag_short()),
               size=FIG_SIZE,
               left_margin=FIG_MARGIN.left, bottom_margin=FIG_MARGIN.bottom,
               right_margin=FIG_MARGIN.right, top_margin=FIG_MARGIN.top,
               legend=:topright)
    plot!(plt, vsm_d[m_vsm, 1], τ_vsm,
          label=@sprintf("vSmartMOM, H₂O=%.2f", h2o),
          color=:red, lw=0.9)

    out = joinpath(OUT_DIR,
                   @sprintf("tau_up_H2O%.2f_%s_%s.pdf",
                            h2o, slice_label, SCEN_SUFFIX))
    save_both(plt, out)
end

"""
    plot_tau_overlay(mod_lo, vsm_lo, mod_hi, vsm_hi, slice_label, lo, hi)

Four-line τ_up overlay: MODTRAN and vSmartMOM at H2OSTR=0.05 and 2.05.
Replicates the existing `tau_up_H2O.pdf` figure but with explicit margins.
"""
function plot_tau_overlay(mod_lo, vsm_lo, mod_hi, vsm_hi, slice_label, lo, hi)
    m_mod = _slice_mask(mod_lo[:,1], lo, hi)
    m_vsm = _slice_mask(vsm_lo[:,1], lo, hi)
    if !any(m_mod) || !any(m_vsm)
        @info "plot_tau_overlay: empty slice" slice_label
        return
    end
    plt = plot(mod_lo[m_mod, 1], tau_up_vertical(mod_lo[m_mod, 6], μ_v),
               label=@sprintf("MODTRAN, H₂O=%.2f", H2O_LO),
               color=:blue, lw=0.9,
               xlabel="Wavelength [nm]",
               ylabel=L"$\tau_\mathrm{atm}$ [upward]",
               title=@sprintf("τ_up overlay (H₂O=%.2f vs %.2f)  [%s]  (%s)",
                              H2O_LO, H2O_HI, slice_label, scen_tag_short()),
               size=FIG_SIZE,
               left_margin=FIG_MARGIN.left, bottom_margin=FIG_MARGIN.bottom,
               right_margin=FIG_MARGIN.right, top_margin=FIG_MARGIN.top,
               legend=:topright)
    plot!(plt, vsm_lo[m_vsm, 1], tau_up_vertical(vsm_lo[m_vsm, 6], μ_v),
          label=@sprintf("vSmartMOM, H₂O=%.2f", H2O_LO),
          color=:cyan, lw=0.9)
    plot!(plt, mod_hi[m_mod, 1], tau_up_vertical(mod_hi[m_mod, 6], μ_v),
          label=@sprintf("MODTRAN, H₂O=%.2f", H2O_HI),
          color=:red, lw=0.9)
    plot!(plt, vsm_hi[m_vsm, 1], tau_up_vertical(vsm_hi[m_vsm, 6], μ_v),
          label=@sprintf("vSmartMOM, H₂O=%.2f", H2O_HI),
          color=:orange, lw=0.9)

    out = joinpath(OUT_DIR,
                   @sprintf("overlay_H2O_%s_%s.pdf", slice_label, SCEN_SUFFIX))
    save_both(plt, out)
end

"""
    plot_tau_ratio(mod_lo, vsm_lo, mod_hi, vsm_hi, slice_label, lo, hi)

Ratio τ(H2OSTR=H2O_HI) / τ(H2OSTR=H2O_LO) for both codes overlaid. The
naive value 41 = 2.05/0.05 is drawn as a horizontal reference — values
near 41 across the full spectrum would imply file-level scaling of all
gases by H2OSTR; physically correct H2O sensitivity should peak in the
H2O bands and approach 1 in atmospheric windows.

Each code is plotted on its own λ grid (no resampling), since both grids
are fine and overlap densely; this avoids interpolation artefacts at
band edges.
"""
function plot_tau_ratio(mod_lo, vsm_lo, mod_hi, vsm_hi, slice_label, lo, hi)
    m_mod = _slice_mask(mod_lo[:,1], lo, hi)
    m_vsm = _slice_mask(vsm_lo[:,1], lo, hi)
    if !any(m_mod) || !any(m_vsm)
        @info "plot_tau_ratio: empty slice" slice_label
        return
    end
    τ_mod_lo = tau_up_vertical(mod_lo[m_mod, 6], μ_v)
    τ_mod_hi = tau_up_vertical(mod_hi[m_mod, 6], μ_v)
    τ_vsm_lo = tau_up_vertical(vsm_lo[m_vsm, 6], μ_v)
    τ_vsm_hi = tau_up_vertical(vsm_hi[m_vsm, 6], μ_v)

    # Floor the denominator to keep the ratio finite at saturated/zero τ.
    eps_τ = 1e-12
    r_mod = τ_mod_hi ./ max.(τ_mod_lo, eps_τ)
    r_vsm = τ_vsm_hi ./ max.(τ_vsm_lo, eps_τ)

    naive = H2O_HI / H2O_LO  # 41.0 for the canonical pair

    plt = plot(mod_lo[m_mod, 1], r_mod,
               label="MODTRAN ratio",
               color=:darkorange, lw=0.9,
               xlabel="Wavelength [nm]",
               ylabel=L"$\tau_\mathrm{up}(H_2O=2.05)\;/\;\tau_\mathrm{up}(H_2O=0.05)$",
               title=@sprintf("τ_up ratio  (H₂O=%.2f / %.2f)  [%s]  (%s)",
                              H2O_HI, H2O_LO, slice_label, scen_tag_short()),
               size=FIG_SIZE,
               left_margin=FIG_MARGIN.left, bottom_margin=FIG_MARGIN.bottom,
               right_margin=FIG_MARGIN.right, top_margin=FIG_MARGIN.top,
               legend=:topright)
    plot!(plt, vsm_lo[m_vsm, 1], r_vsm,
          label="vSmartMOM ratio",
          color=:darkblue, lw=0.9)
    hline!(plt, [naive],
           label=@sprintf("naive H₂O scaling = %g", naive),
           color=:black, ls=:dash, lw=1.0)

    out = joinpath(OUT_DIR,
                   @sprintf("ratio_H2O_%s_%s.pdf", slice_label, SCEN_SUFFIX))
    save_both(plt, out)
end

# === main =================================================================
function main()
    println("Reading MODTRAN/vSmartMOM .dat for H₂O=$(H2O_LO) and $(H2O_HI) …")
    mod_lo = load_dat(mod_path(H2O_LO))
    mod_hi = load_dat(mod_path(H2O_HI))
    vsm_lo = load_dat(vsm_path(H2O_LO))
    vsm_hi = load_dat(vsm_path(H2O_HI))

    @assert size(mod_lo, 2) ≥ 6 "MODTRAN file missing tu_dir column"
    @assert size(vsm_lo, 2) ≥ 6 "vSmartMOM file missing tu_dir column"

    n_plots = 0
    for (slice, lo, hi) in SLICES
        plot_tau_one_h2o(mod_lo, vsm_lo, H2O_LO, slice, lo, hi)
        plot_tau_one_h2o(mod_hi, vsm_hi, H2O_HI, slice, lo, hi)
        plot_tau_overlay(mod_lo, vsm_lo, mod_hi, vsm_hi, slice, lo, hi)
        plot_tau_ratio(mod_lo, vsm_lo, mod_hi, vsm_hi, slice, lo, hi)
        n_plots += 4
    end
    println("Slices plotted: ", length(SLICES))
    println("Linear PDFs written: ", n_plots, " (× 2 with logplots/ twins)")
    println("Output: ", OUT_DIR)
end

main()
