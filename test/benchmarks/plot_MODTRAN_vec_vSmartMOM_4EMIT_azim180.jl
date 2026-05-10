#=
plot_MODTRAN_vec_vSmartMOM_4EMIT_azim180.jl — azim180 convention-probe
variant of the MODTRAN/vSmartMOM EMIT comparison plotter.

This plotter is identical to its parent script except that it
reads vSmartMOM .dat files from
  $HOME/data/EMIT_MODTRANcomp/azim180/vSmartMOM_out/...
and writes plots into
  $HOME/data/EMIT_MODTRANcomp/azim180/plots/...
with `_azim180` injected into every output filename.

Background. vSmartMOM defines vaz = 0° as Sun and satellite on
*opposite* sides of zenith in the principal plane, and vaz = 180°
as Sun directly *behind* the satellite (same side of zenith). Many
other models — MODTRAN among them — define RAA = 0 as
Sun-behind-detector, which corresponds to vSmartMOM vaz = 180°. The azim180 YAML series
(`ParamsEMIT_MODTRANcomp_newLUT*_azim180.yaml`) tests this
convention-mismatch hypothesis by running vSmartMOM with
vaz = (parent vaz) - 180°. If the hypothesis is correct, the
azim180 outputs should match MODTRAN better than the parent
outputs.

MODTRAN data is unchanged between the two comparisons; only the
vSmartMOM model output and the resulting plot directory differ.
=#
using Plots
using DelimitedFiles
using LaTeXStrings

default(alpha = 0.5)

"""
    _savefig_checked(name)

`Plots.savefig` opens the destination via `fopen("wb")` (creates a 0-byte
file) before GR renders and writes the content. If GR fails silently — most
often because the GKS workstation is in a wedged state — the file is left
at 0 bytes with no exception thrown. We detect that, attempt one recovery
pass (closeall + re-save), and remove the file if it stays empty so callers
can see the failure instead of inheriting a corrupt artifact.
"""
function _savefig_checked(name::AbstractString)
    try
        Plots.savefig(name)
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
            Plots.savefig(name)
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
    _LOG_FLOOR

Tiny positive value used to clamp non-positive (and non-finite) y-data
before plotting. Workaround for a Plots/GR PDF rendering bug: when
`yaxis!(:log10)` runs against y-data containing zero or negative values,
GR emits literal `nan`/`inf` tokens in the PDF path operators. The PDF
spec has no such operator, so any reader silently drops the affected
paths — the user sees a plot frame with no curves.

Mutating series y after `plot()` doesn't help because Plots/GR caches
the rendering state independently of `series.plotattributes[:y]`. The
fix has to apply at plot-call time, hence the source-level helpers
below (`_pos_floor!` for matrix columns, `_pos_floor.(...)` for vectors).
1e-30 is visually indistinguishable from 0 on any sensible log10 axis
and keeps every drawing operation well-defined.
"""
const _LOG_FLOOR = 1e-30

"""
    _pos_floor(x) -> y

Element-wise clamp of `x` that replaces non-positive or non-finite
values with `_LOG_FLOOR`. Use this at the *plot call site* to make a
y-data vector log-axis-safe.
"""
_pos_floor(x) = map(v -> isfinite(v) && v > 0 ? float(v) : float(_LOG_FLOOR), x)

"""
    _pos_floor!(d, cols)

In-place clamp of selected columns of a matrix `d`: every non-positive
or non-finite entry becomes `_LOG_FLOOR`. Apply once after `readdlm`
to make any direct-value plot of those columns log-axis-safe; τ-style
plots (`-log(d) · μ`) should still wrap their log argument explicitly.
"""
function _pos_floor!(d::AbstractMatrix, cols)
    for j in cols, i in axes(d, 1)
        v = d[i, j]
        if !(isfinite(v) && v > 0)
            d[i, j] = _LOG_FLOOR
        end
    end
    return d
end

"""
    savefig_both(name)

Save the current plot at `name` (linear y-axis) and also at
`<dir>/logplots/<stem>_log<ext>` with `yscale=:log10`. Restores the
y-axis scale to linear afterwards so subsequent `plot!` calls are
unaffected.

Wraps each save in a try/catch: on PlotUtils/Plots backends, switching
to `:log10` over data containing non-positive values can collapse the
plot layout and trigger `AssertionError: total_plotarea_vertical > 0mm`.
We skip the affected save rather than aborting the whole script. To
avoid the related "0-byte file" failure mode see `_savefig_checked`,
and to avoid the "log10 → `nan` PDF operators" rendering bug see
`_pos_floor` (clamp y-values at the plot call, not here — post-hoc
mutation is not picked up by the GR backend).
"""
function savefig_both(name::AbstractString)
    dir = dirname(name)
    isempty(dir) && (dir = ".")
    mkpath(dir)
    saved_linear = _savefig_checked(name)
    log_dir = joinpath(dir, "logplots")
    mkpath(log_dir)
    stem, ext = splitext(basename(name))
    saved_linear || return
    log_name = joinpath(log_dir, stem * "_log" * ext)
    try
        yaxis!(:log10)
        _savefig_checked(log_name)
    catch err
        @warn "savefig_both: failed to save log-axis plot" name exception=err
        try; Plots.closeall(); catch; end
        isfile(log_name) && filesize(log_name) == 0 && rm(log_name; force=true)
    finally
        try; yaxis!(:identity); catch; end
    end
end

function fmt(x, n)
    s = string(round(x; digits=n))
    i = findfirst('.', s)
    i === nothing ? s * "." * "0"^n : rpad(s, i + n, '0')
end

H2O = collect(0.05:0.5:2.5)
AOT = collect(0.0001:0.1:0.7)
GNDALT = collect(0.0:5.0)
TSZ = [168.0, 178.0]
h2o = H2O[1]
aot = AOT[1]
gndalt = GNDALT[1]
tsz = TSZ[1]


MOD_DIR = "/home/sanghavi/data/EMIT_MODTRANcomp/MODTRAN_out/"
filename1(h2o, aot, gndalt, tsz) =
           "MODTRAN_H2O$(fmt(h2o,4))_AOT$(fmt(aot,4))_GNDALT$(fmt(gndalt,3))_TSZ$(fmt(tsz,1)).dat"
path1 = joinpath(MOD_DIR, filename1(h2o, aot, gndalt, tsz))
data1 = readdlm(path1, skipstart=2)
_pos_floor!(data1, 3:8)

MOM_DIR = "/home/sanghavi/data/EMIT_MODTRANcomp/azim180/vSmartMOM_out/vector_modtran_equiv/"
filename2(h2o, aot, gndalt, tsz) =
           "vSmartMOM_modtran_equiv_H2O$(fmt(h2o,4))_AOT$(fmt(aot,4))_GNDALT$(fmt(gndalt,3))_TSZ$(fmt(tsz,1)).dat"
path2 = joinpath(MOM_DIR, filename2(h2o, aot, gndalt, tsz))
# The vSmartMOM file has two '#'-prefixed header lines; skipping only one would
# leave the tokenised column-name line in data2 and splay sphalb into a phantom
# 9th column of empty strings.
data2 = readdlm(path2, comments=true, comment_char='#')
_pos_floor!(data2, 3:8)

ylabels = [L"$\rho_\mathrm{atm}$", L"$T^\downarrow_\mathrm{dir}$", L"$T^\downarrow_\mathrm{dif}$", L"$T^\uparrow_\mathrm{dir}$", L"$T^\uparrow_\mathrm{dif}$", L"$S_\mathrm{atm}$"]
xlabel = "Wavelength [nm]"
for i=3:8
    plot(data1[:,1], data1[:,i], label="MODTRAN")#, xaxis=:log, yaxis=:log)
    plot!(data2[:,1], data2[:,i], label="vSmartMOM IQU")#, xaxis=:log, yaxis=:log)
    xlabel!(xlabel)
    ylabel!(ylabels[i-2])
    savefig_both("/home/sanghavi/data/EMIT_MODTRANcomp/azim180/plots/MODTRAN_vSmartMOM_4EMIT_azim180_IQU_$(i-2).pdf")
end

vza = 12.0
sza = 39.583965990615404
ylabels = [L"$\tau_\mathrm{atm}$ [downward]", L"$\tau_\mathrm{atm}$ [upward]"]
for i in [4,6]
    if i==4
        plot(data1[:,1], -log.(data1[:,i])*cosd(sza), label="MODTRAN")#, xaxis=:log, yaxis=:log)
        plot!(data2[:,1], -log.(data2[:,i])*cosd(sza), label="vSmartMOM IQU")#, xaxis=:log, yaxis=:log)
        ylabel!(ylabels[1])
        xlabel!(xlabel)
        savefig_both("/home/sanghavi/data/EMIT_MODTRANcomp/azim180/plots/MODTRAN_vSmartMOM_4EMIT_azim180_IQU_tau_down.pdf")
    else
        plot(data1[:,1], -log.(data1[:,i])*cosd(vza), label="MODTRAN")#, xaxis=:log, yaxis=:log)
        plot!(data2[:,1], -log.(data2[:,i])*cosd(vza), label="vSmartMOM IQU")#, xaxis=:log, yaxis=:log)
        ylabel!(ylabels[2])
        xlabel!(xlabel)
        savefig_both("/home/sanghavi/data/EMIT_MODTRANcomp/azim180/plots/MODTRAN_vSmartMOM_4EMIT_azim180_IQU_tau_up.pdf")
    end
end
    
plot(data1[:,1], -log.(data1[:,4])*cosd(sza) + log.(data1[:,6])*cosd(vza), label="MODTRAN")#, xaxis=:log, yaxis=:log)
plot!(data2[:,1], -log.(data2[:,4])*cosd(sza) + log.(data2[:,6])*cosd(vza), label="vSmartMOM IQU")#, xaxis=:log, yaxis=:log)
ylabel!(L"$\Delta\tau_\mathrm{atm}$ [down-up]")
xlabel!(xlabel)
savefig_both("/home/sanghavi/data/EMIT_MODTRANcomp/azim180/plots/MODTRAN_vSmartMOM_4EMIT_azim180_IQU_tau_diff.pdf")

if false  # DISABLED for azim180: not a full-spectrum MODTRAN-equivalent comparison
ind_uv1 = findall(x -> 350<=x<= 400.0, data1[:,1])
ind_uv2 = findall(x -> 350<=x<= 400.0, data2[:,1])

ind_vis1 = findall(x -> 500<=x<=700.0, data1[:,1])
ind_vis2 = findall(x -> 500<=x<=700.0, data2[:,1])

ind_nir1 = findall(x -> 950<=x<=1100.0, data1[:,1])
ind_nir2 = findall(x -> 950<=x<=1100.0, data2[:,1])

ind_swir11 = findall(x -> 1500<=x<=1700.0, data1[:,1])
ind_swir12 = findall(x -> 1500<=x<=1700.0, data2[:,1])

ind_swir21 = findall(x -> 2100<=x<=2200.0, data1[:,1])
ind_swir22 = findall(x -> 2100<=x<=2200.0, data2[:,1])

spec_slice = ["uv", "vis", "nir", "swir1", "swir2"]

for ctr=1:5
    if ctr==1
        ind1 = ind_uv1
        ind2 = ind_uv2
    elseif ctr==2
        ind1 = ind_vis1
        ind2 = ind_vis2
    elseif ctr==3
        ind1 = ind_nir1
        ind2 = ind_nir2
    elseif ctr==4
        ind1 = ind_swir11           
        ind2 = ind_swir12
    elseif ctr==5
        ind1 = ind_swir21
        ind2 = ind_swir22
    end
    if isempty(ind1) || isempty(ind2)
        @info "Skipping spec_slice $(spec_slice[ctr]) — empty index range (data1=$(length(ind1)), data2=$(length(ind2)))"
        continue
    end
    local ylabels = [L"$\rho_\mathrm{atm}$", L"$T^\downarrow_\mathrm{dir}$", L"$T^\downarrow_\mathrm{dif}$", L"$T^\uparrow_\mathrm{dir}$", L"$T^\uparrow_\mathrm{dif}$", L"$S_\mathrm{atm}$"]
    local xlabel = "Wavelength [nm]"
    for i=3:8
        plot(data1[ind1,1], data1[ind1,i], label="MODTRAN")#, xaxis=:log, yaxis=:log)
        plot!(data2[ind2,1], data2[ind2,i], label="vSmartMOM IQU")#, xaxis=:log, yaxis=:log)
        xlabel!(xlabel)
        ylabel!(ylabels[i-2])
        savefig_both("/home/sanghavi/data/EMIT_MODTRANcomp/azim180/plots/MODTRAN_vSmartMOM_4EMIT_azim180_IQU_$(i-2)"*"_"*spec_slice[ctr]*".pdf")
    end

    ylabels = [L"$\tau_\mathrm{atm}$ [downward]", L"$\tau_\mathrm{atm}$ [upward]"]
    for i in [4,6]
        if i==4
            plot(data1[ind1,1], -log.(data1[ind1,i])*cosd(sza), label="MODTRAN")#, xaxis=:log, yaxis=:log)
            plot!(data2[ind2,1], -log.(data2[ind2,i])*cosd(sza), label="vSmartMOM IQU")#, xaxis=:log, yaxis=:log)
            ylabel!(ylabels[1])
            xlabel!(xlabel)
            savefig_both("/home/sanghavi/data/EMIT_MODTRANcomp/azim180/plots/MODTRAN_vSmartMOM_4EMIT_azim180_IQU_tau_down_"*spec_slice[ctr]*".pdf")
        else
            plot(data1[ind1,1], -log.(data1[ind1,i])*cosd(vza), label="MODTRAN")#, xaxis=:log, yaxis=:log)
            plot!(data2[ind2,1], -log.(data2[ind2,i])*cosd(vza), label="vSmartMOM IQU")#, xaxis=:log, yaxis=:log)
            ylabel!(ylabels[2])
            xlabel!(xlabel)
            savefig_both("/home/sanghavi/data/EMIT_MODTRANcomp/azim180/plots/MODTRAN_vSmartMOM_4EMIT_azim180_IQU_tau_up_"*spec_slice[ctr]*".pdf")
        end
    end
    
    plot(data1[ind1,1], -log.(data1[ind1,4])*cosd(sza) + log.(data1[ind1,6])*cosd(vza), label="MODTRAN")#, xaxis=:log, yaxis=:log)
    plot!(data2[ind2,1], -log.(data2[ind2,4])*cosd(sza) + log.(data2[ind2,6])*cosd(vza), label="vSmartMOM IQU")#, xaxis=:log, yaxis=:log)
    ylabel!(L"$\Delta\tau_\mathrm{atm}$ [down-up]")
    xlabel!(xlabel)
    savefig_both("/home/sanghavi/data/EMIT_MODTRANcomp/azim180/plots/MODTRAN_vSmartMOM_4EMIT_azim180_IQU_tau_diff_"*spec_slice[ctr]*".pdf")
end
end  # /azim180-disabled

#=================================================================================================#
# Compare H2O0.05 and H2O2.05 cases to see how the difference between MODTRAN and vSmartMOM changes with H2O column.
#=================================================================================================#

H2O = collect(0.05:0.5:2.5)
AOT = collect(0.0001:0.1:0.7)
GNDALT = collect(0.0:5.0)
TSZ = [168.0, 178.0]
h2o = H2O[end]
aot = AOT[1]
gndalt = GNDALT[1]
tsz = TSZ[1]


MOD_DIR = "/home/sanghavi/data/EMIT_MODTRANcomp/MODTRAN_out/"
filename11(h2o, aot, gndalt, tsz) =
           "MODTRAN_H2O$(fmt(h2o,4))_AOT$(fmt(aot,4))_GNDALT$(fmt(gndalt,3))_TSZ$(fmt(tsz,1)).dat"
path11 = joinpath(MOD_DIR, filename11(h2o, aot, gndalt, tsz))
data11 = readdlm(path11, skipstart=2)
_pos_floor!(data11, 3:8)

MOM_DIR = "/home/sanghavi/data/EMIT_MODTRANcomp/azim180/vSmartMOM_out/vector_modtran_equiv/"
filename21(h2o, aot, gndalt, tsz) =
           "vSmartMOM_modtran_equiv_H2O$(fmt(h2o,4))_AOT$(fmt(aot,4))_GNDALT$(fmt(gndalt,3))_TSZ$(fmt(tsz,1)).dat"
path21 = joinpath(MOM_DIR, filename21(h2o, aot, gndalt, tsz))
# The vSmartMOM file has two '#'-prefixed header lines; skipping only one would
# leave the tokenised column-name line in data2 and splay sphalb into a phantom
# 9th column of empty strings.
data21 = readdlm(path21, comments=true, comment_char='#')
_pos_floor!(data21, 3:8)

ylabels = [L"$\rho_\mathrm{atm}$", L"$T^\downarrow_\mathrm{dir}$", L"$T^\downarrow_\mathrm{dif}$", L"$T^\uparrow_\mathrm{dir}$", L"$T^\uparrow_\mathrm{dif}$", L"$S_\mathrm{atm}$"]
xlabel = "Wavelength [nm]"
for i=3:8
    plot(data1[:,1], data1[:,i], label="MODTRAN, H₂O=0.05")#, xaxis=:log, yaxis=:log)
    plot!(data2[:,1], data2[:,i], label="vSmartMOM IQU, H₂O=0.05")#, xaxis=:log, yaxis=:log)
    plot!(data11[:,1], data11[:,i], label="MODTRAN, H₂O=2.05")#, xaxis=:log, yaxis=:log)
    plot!(data21[:,1], data21[:,i], label="vSmartMOM IQU, H₂O=2.05")
    xlabel!(xlabel)
    ylabel!(ylabels[i-2])
    savefig_both("/home/sanghavi/data/EMIT_MODTRANcomp/azim180/plots/MODTRAN_vSmartMOM_4EMIT_azim180_IQU_$(i-2)_H2O.pdf")
end

vza = 12.0
sza = 39.583965990615404
ylabels = [L"$\tau_\mathrm{atm}$ [downward]", L"$\tau_\mathrm{atm}$ [upward]"]
for i in [4,6]
    if i==4
        plot(data1[:,1], -log.(data1[:,i])*cosd(sza), label="MODTRAN, H₂O=0.05")#, xaxis=:log, yaxis=:log)
        plot!(data2[:,1], -log.(data2[:,i])*cosd(sza), label="vSmartMOM IQU, H₂O=0.05")#, xaxis=:log, yaxis=:log)
        plot!(data11[:,1], -log.(data11[:,i])*cosd(sza), label="MODTRAN, H₂O=2.05")#, xaxis=:log, yaxis=:log
        plot!(data21[:,1], -log.(data21[:,i])*cosd(sza), label="vSmartMOM IQU, H₂O=2.05")#, xaxis=:log, yaxis=:log)
        ylabel!(ylabels[1])
        xlabel!(xlabel)
        savefig_both("/home/sanghavi/data/EMIT_MODTRANcomp/azim180/plots/MODTRAN_vSmartMOM_4EMIT_azim180_IQU_tau_down_H2O.pdf")
    else
        plot(data1[:,1], -log.(data1[:,i])*cosd(vza), label="MODTRAN, H₂O=0.05")#, xaxis=:log, yaxis=:log)
        plot!(data2[:,1], -log.(data2[:,i])*cosd(vza), label="vSmartMOM IQU, H₂O=0.05")#, xaxis=:log, yaxis=:log)
        plot!(data11[:,1], -log.(data11[:,i])*cosd(vza), label="MODTRAN, H₂O=2.05")#, xaxis=:log, yaxis=:log)
        plot!(data21[:,1], -log.(data21[:,i])*cosd(vza), label="vSmartMOM IQU, H₂O=2.05")#, xaxis=:log, yaxis=:log)
        ylabel!(ylabels[2])
        xlabel!(xlabel)
        savefig_both("/home/sanghavi/data/EMIT_MODTRANcomp/azim180/plots/MODTRAN_vSmartMOM_4EMIT_azim180_IQU_tau_up_H2O.pdf")
    end
end
    
plot(data1[:,1], -log.(data1[:,4])*cosd(sza) + log.(data1[:,6])*cosd(vza), label="MODTRAN, H₂O=0.05")#, xaxis=:log, yaxis=:log)
plot!(data2[:,1], -log.(data2[:,4])*cosd(sza) + log.(data2[:,6])*cosd(vza), label="vSmartMOM IQU, H₂O=0.05")#, xaxis=:log, yaxis=:log)
plot!(data11[:,1], -log.(data11[:,4])*cosd(sza) + log.(data11[:,6])*cosd(vza), label="MODTRAN, H₂O=2.05")#, xaxis=:log, yaxis=:log)
plot!(data21[:,1], -log.(data21[:,4])*cosd(sza) + log.(data21[:,6])*cosd(vza), label="vSmartMOM IQU, H₂O=2.05")#, xaxis=:log, yaxis=:log)
ylabel!(L"$\Delta\tau_\mathrm{atm}$ [down-up]")
xlabel!(xlabel)
savefig_both("/home/sanghavi/data/EMIT_MODTRANcomp/azim180/plots/MODTRAN_vSmartMOM_4EMIT_azim180_IQU_tau_diff_H2O.pdf")

if false  # DISABLED for azim180: not a full-spectrum MODTRAN-equivalent comparison
ind_uv1 = findall(x -> 350<=x<= 400.0, data1[:,1])
ind_uv2 = findall(x -> 350<=x<= 400.0, data2[:,1])

ind_vis1 = findall(x -> 500<=x<=700.0, data1[:,1])
ind_vis2 = findall(x -> 500<=x<=700.0, data2[:,1])

ind_nir1 = findall(x -> 950<=x<=1100.0, data1[:,1])
ind_nir2 = findall(x -> 950<=x<=1100.0, data2[:,1])

ind_swir11 = findall(x -> 1500<=x<=1700.0, data1[:,1])
ind_swir12 = findall(x -> 1500<=x<=1700.0, data2[:,1])

ind_swir21 = findall(x -> 2100<=x<=2200.0, data1[:,1])
ind_swir22 = findall(x -> 2100<=x<=2200.0, data2[:,1])

spec_slice = ["uv", "vis", "nir", "swir1", "swir2"]

for ctr=1:5
    if ctr==1
        ind1 = ind_uv1
        ind2 = ind_uv2
    elseif ctr==2
        ind1 = ind_vis1
        ind2 = ind_vis2
    elseif ctr==3
        ind1 = ind_nir1
        ind2 = ind_nir2
    elseif ctr==4
        ind1 = ind_swir11           
        ind2 = ind_swir12
    elseif ctr==5
        ind1 = ind_swir21
        ind2 = ind_swir22
    end
    if isempty(ind1) || isempty(ind2)
        @info "Skipping spec_slice $(spec_slice[ctr]) (H2O comparison) — empty index range"
        continue
    end
    local ylabels = [L"$\rho_\mathrm{atm}$", L"$T^\downarrow_\mathrm{dir}$", L"$T^\downarrow_\mathrm{dif}$", L"$T^\uparrow_\mathrm{dir}$", L"$T^\uparrow_\mathrm{dif}$", L"$S_\mathrm{atm}$"]
    local xlabel = "Wavelength [nm]"
    for i=3:8
        plot(data1[ind1,1], data1[ind1,i], label="MODTRAN, H₂O=0.05")#, xaxis=:log, yaxis=:log)
        plot!(data2[ind2,1], data2[ind2,i], label="vSmartMOM IQU, H₂O=0.05")#, xaxis=:log, yaxis=:log)
        plot!(data11[ind1,1], data11[ind1,i], label="MODTRAN, H₂O=2.05")#, xaxis=:log, yaxis=:log)
        plot!(data21[ind2,1], data21[ind2,i], label="vSmartMOM IQU, H₂O=2.05")#, xaxis=:log, yaxis=:log)
        xlabel!(xlabel)
        ylabel!(ylabels[i-2])
        savefig_both("/home/sanghavi/data/EMIT_MODTRANcomp/azim180/plots/MODTRAN_vSmartMOM_4EMIT_azim180_IQU_$(i-2)"*"_"*spec_slice[ctr]*"_H2O.pdf")
    end

    ylabels = [L"$\tau_\mathrm{atm}$ [downward]", L"$\tau_\mathrm{atm}$ [upward]"]
    for i in [4,6]
        if i==4
            plot(data1[ind1,1], -log.(data1[ind1,i])*cosd(sza), label="MODTRAN, H₂O=0.05")#, xaxis=:log, yaxis=:log)
            plot!(data2[ind2,1], -log.(data2[ind2,i])*cosd(sza), label="vSmartMOM IQU, H₂O=0.05")#, xaxis=:log, yaxis=:log)
            plot!(data11[ind1,1], -log.(data11[ind1,i])*cosd(sza), label="MODTRAN, H₂O=2.05")#, xaxis=:log, yaxis=:log)
            plot!(data21[ind2,1], -log.(data21[ind2,i])*cosd(sza), label="vSmartMOM IQU, H₂O=2.05")#, xaxis=:log, yaxis=:log)
            ylabel!(ylabels[1])
            xlabel!(xlabel)
            savefig_both("/home/sanghavi/data/EMIT_MODTRANcomp/azim180/plots/MODTRAN_vSmartMOM_4EMIT_azim180_IQU_tau_down_"*spec_slice[ctr]*"_H2O.pdf")
        else
            plot(data1[ind1,1], -log.(data1[ind1,i])*cosd(vza), label="MODTRAN, H₂O=0.05")#, xaxis=:log, yaxis=:log)
            plot!(data2[ind2,1], -log.(data2[ind2,i])*cosd(vza), label="vSmartMOM IQU, H₂O=0.05")#, xaxis=:log, yaxis=:log)
            plot!(data11[ind1,1], -log.(data11[ind1,i])*cosd(vza), label="MODTRAN, H₂O=2.05")#, xaxis=:log, yaxis=:log)
            plot!(data21[ind2,1], -log.(data21[ind2,i])*cosd(vza), label="vSmartMOM IQU, H₂O=2.05")#, xaxis=:log, yaxis=:log)
            ylabel!(ylabels[2])
            xlabel!(xlabel)
            savefig_both("/home/sanghavi/data/EMIT_MODTRANcomp/azim180/plots/MODTRAN_vSmartMOM_4EMIT_azim180_IQU_tau_up_"*spec_slice[ctr]*"_H2O.pdf")
        end
    end
    
    plot(data1[ind1,1], -log.(data1[ind1,4])*cosd(sza) + log.(data1[ind1,6])*cosd(vza), label="MODTRAN, H₂O=0.05")#, xaxis=:log, yaxis=:log)
    plot!(data2[ind2,1], -log.(data2[ind2,4])*cosd(sza) + log.(data2[ind2,6])*cosd(vza), label="vSmartMOM IQU, H₂O=0.05")#, xaxis=:log, yaxis=:log)
    plot!(data11[ind1,1], -log.(data11[ind1,4])*cosd(sza) + log.(data11[ind1,6])*cosd(vza), label="MODTRAN, H₂O=2.05")#, xaxis=:log, yaxis=:log)
    plot!(data21[ind2,1], -log.(data21[ind2,4])*cosd(sza) + log.(data21[ind2,6])*cosd(vza), label="vSmartMOM IQU, H₂O=2.05")#, xaxis=:log, yaxis=:log)
    ylabel!(L"$\Delta\tau_\mathrm{atm}$ [down-up]")
    xlabel!(xlabel)
    savefig_both("/home/sanghavi/data/EMIT_MODTRANcomp/azim180/plots/MODTRAN_vSmartMOM_4EMIT_azim180_IQU_tau_diff_"*spec_slice[ctr]*"_H2O.pdf")
end
end  # /azim180-disabled

if false  # DISABLED for azim180: not a full-spectrum MODTRAN-equivalent comparison
#====== only vSmartMOM ======#
ylabels = [L"$\rho_\mathrm{atm}$", L"$T^\downarrow_\mathrm{dir}$", L"$T^\downarrow_\mathrm{dif}$", L"$T^\uparrow_\mathrm{dir}$", L"$T^\uparrow_\mathrm{dif}$", L"$S_\mathrm{atm}$"]
xlabel = "Wavelength [nm]"
for i=3:8
    plot(data21[:,1], data21[:,i], label="vSmartMOM IQU, H₂O=0.05")#, xaxis=:log, yaxis=:log)
    plot!(data21[:,1], data21[:,i], label="vSmartMOM IQU, H₂O=2.05")
    xlabel!(xlabel)
    ylabel!(ylabels[i-2])
    savefig_both("/home/sanghavi/data/EMIT_MODTRANcomp/azim180/plots/h2o_sensitivity/vSmartMOM_4EMIT_azim180_IQU_$(i-2)_H2O.pdf")
end

vza = 12.0
sza = 39.583965990615404
ylabels = [L"$\tau_\mathrm{atm}$ [downward]", L"$\tau_\mathrm{atm}$ [upward]"]
for i in [4,6]
    if i==4
        plot(data2[:,1], -log.(data2[:,i])*cosd(sza), label="vSmartMOM IQU, H₂O=0.05")#, xaxis=:log, yaxis=:log)
        plot!(data21[:,1], -log.(data21[:,i])*cosd(sza), label="vSmartMOM IQU, H₂O=2.05")#, xaxis=:log, yaxis=:log)
        ylabel!(ylabels[1])
        xlabel!(xlabel)
        savefig_both("/home/sanghavi/data/EMIT_MODTRANcomp/azim180/plots/h2o_sensitivity/vSmartMOM_4EMIT_azim180_IQU_tau_down_H2O.pdf")
    else
        plot(data2[:,1], -log.(data2[:,i])*cosd(vza), label="vSmartMOM IQU, H₂O=0.05")#, xaxis=:log, yaxis=:log)
        plot!(data21[:,1], -log.(data21[:,i])*cosd(vza), label="vSmartMOM IQU, H₂O=2.05")#, xaxis=:log, yaxis=:log)
        ylabel!(ylabels[2])
        xlabel!(xlabel)
        savefig_both("/home/sanghavi/data/EMIT_MODTRANcomp/azim180/plots/h2o_sensitivity/vSmartMOM_4EMIT_azim180_IQU_tau_up_H2O.pdf")
    end
end
    
plot(data2[:,1], -log.(data2[:,4])*cosd(sza) + log.(data2[:,6])*cosd(vza), label="vSmartMOM IQU, H₂O=0.05")#, 
plot!(data21[:,1], -log.(data21[:,4])*cosd(sza) + log.(data21[:,6])*cosd(vza), label="vSmartMOM IQU, H₂O=2.05")#, xaxis=:log, yaxis=:log)
ylabel!(L"$\Delta\tau_\mathrm{atm}$ [down-up]")
xlabel!(xlabel)
savefig_both("/home/sanghavi/data/EMIT_MODTRANcomp/azim180/plots/h2o_sensitivity/vSmartMOM_4EMIT_azim180_IQU_tau_diff_H2O.pdf")

ind_uv1 = findall(x -> 350<=x<= 400.0, data1[:,1])
ind_uv2 = findall(x -> 350<=x<= 400.0, data2[:,1])

ind_vis1 = findall(x -> 500<=x<=700.0, data1[:,1])
ind_vis2 = findall(x -> 500<=x<=700.0, data2[:,1])

ind_nir1 = findall(x -> 950<=x<=1100.0, data1[:,1])
ind_nir2 = findall(x -> 950<=x<=1100.0, data2[:,1])

ind_swir11 = findall(x -> 1500<=x<=1700.0, data1[:,1])
ind_swir12 = findall(x -> 1500<=x<=1700.0, data2[:,1])

ind_swir21 = findall(x -> 2100<=x<=2200.0, data1[:,1])
ind_swir22 = findall(x -> 2100<=x<=2200.0, data2[:,1])

spec_slice = ["uv", "vis", "nir", "swir1", "swir2"]

for ctr=1:5
    if ctr==1
        ind1 = ind_uv1
        ind2 = ind_uv2
    elseif ctr==2
        ind1 = ind_vis1
        ind2 = ind_vis2
    elseif ctr==3
        ind1 = ind_nir1
        ind2 = ind_nir2
    elseif ctr==4
        ind1 = ind_swir11           
        ind2 = ind_swir12
    elseif ctr==5
        ind1 = ind_swir21
        ind2 = ind_swir22
    end
    if isempty(ind2)
        @info "Skipping spec_slice $(spec_slice[ctr]) (vSmartMOM-only) — empty data2 index range"
        continue
    end
    local ylabels = [L"$\rho_\mathrm{atm}$", L"$T^\downarrow_\mathrm{dir}$", L"$T^\downarrow_\mathrm{dif}$", L"$T^\uparrow_\mathrm{dir}$", L"$T^\uparrow_\mathrm{dif}$", L"$S_\mathrm{atm}$"]
    local xlabel = "Wavelength [nm]"
    for i=3:8
        plot(data2[ind2,1], data2[ind2,i], label="vSmartMOM IQU, H₂O=0.05")#, xaxis=:log, yaxis=:log)
        plot!(data21[ind2,1], data21[ind2,i], label="vSmartMOM IQU, H₂O=2.05")#, xaxis=:log, yaxis=:log)
        xlabel!(xlabel)
        ylabel!(ylabels[i-2])
        savefig_both("/home/sanghavi/data/EMIT_MODTRANcomp/azim180/plots/h2o_sensitivity/vSmartMOM_4EMIT_azim180_IQU_$(i-2)"*"_"*spec_slice[ctr]*"_H2O.pdf")
    end

    ylabels = [L"$\tau_\mathrm{atm}$ [downward]", L"$\tau_\mathrm{atm}$ [upward]"]
    for i in [4,6]
        if i==4
            plot(data2[ind2,1], -log.(data2[ind2,i])*cosd(sza), label="vSmartMOM IQU, H₂O=0.05")#, xaxis=:log, yaxis=:log)
            plot!(data21[ind2,1], -log.(data21[ind2,i])*cosd(sza), label="vSmartMOM IQU, H₂O=2.05")#, xaxis=:log, yaxis=:log)
            ylabel!(ylabels[1])
            xlabel!(xlabel)
            savefig_both("/home/sanghavi/data/EMIT_MODTRANcomp/azim180/plots/h2o_sensitivity/vSmartMOM_4EMIT_azim180_IQU_tau_down_"*spec_slice[ctr]*"_H2O.pdf")
        else
            plot(data2[ind2,1], -log.(data2[ind2,i])*cosd(vza), label="vSmartMOM IQU, H₂O=0.05")#, xaxis=:log, yaxis=:log)
            plot!(data21[ind2,1], -log.(data21[ind2,i])*cosd(vza), label="vSmartMOM IQU, H₂O=2.05")#, xaxis=:log, yaxis=:log)
            ylabel!(ylabels[2])
            xlabel!(xlabel)
            savefig_both("/home/sanghavi/data/EMIT_MODTRANcomp/azim180/plots/h2o_sensitivity/vSmartMOM_4EMIT_azim180_IQU_tau_up_"*spec_slice[ctr]*"_H2O.pdf")
        end
    end
    
    plot(data2[ind2,1], -log.(data2[ind2,4])*cosd(sza) + log.(data2[ind2,6])*cosd(vza), label="vSmartMOM IQU, H₂O=0.05")#, xaxis=:log, yaxis=:log)
    plot!(data21[ind2,1], -log.(data21[ind2,4])*cosd(sza) + log.(data21[ind2,6])*cosd(vza), label="vSmartMOM IQU, H₂O=2.05")#, xaxis=:log, yaxis=:log)
    ylabel!(L"$\Delta\tau_\mathrm{atm}$ [down-up]")
    xlabel!(xlabel)
    savefig_both("/home/sanghavi/data/EMIT_MODTRANcomp/azim180/plots/h2o_sensitivity/vSmartMOM_4EMIT_azim180_IQU_tau_diff_"*spec_slice[ctr]*"_H2O.pdf")
end
#====== only MODTRAN ======#
ylabels = [L"$\rho_\mathrm{atm}$", L"$T^\downarrow_\mathrm{dir}$", L"$T^\downarrow_\mathrm{dif}$", L"$T^\uparrow_\mathrm{dir}$", L"$T^\uparrow_\mathrm{dif}$", L"$S_\mathrm{atm}$"]
xlabel = "Wavelength [nm]"
for i=3:8
    plot(data11[:,1], data11[:,i], label="MODTRAN, H₂O=0.05")#, xaxis=:log, yaxis=:log)
    plot!(data11[:,1], data11[:,i], label="MODTRAN, H₂O=2.05")#, xaxis=:log, yaxis=:log)

    xlabel!(xlabel)
    ylabel!(ylabels[i-2])
    savefig_both("/home/sanghavi/data/EMIT_MODTRANcomp/azim180/plots/h2o_sensitivity/MODTRAN_4EMIT_azim180_IQU_$(i-2)_H2O.pdf")
end

vza = 12.0
sza = 39.583965990615404
ylabels = [L"$\tau_\mathrm{atm}$ [downward]", L"$\tau_\mathrm{atm}$ [upward]"]
for i in [4,6]
    if i==4
        plot(data1[:,1], -log.(data1[:,i])*cosd(sza), label="MODTRAN, H₂O=0.05")#, xaxis=:log, yaxis=:log)
        plot!(data11[:,1], -log.(data11[:,i])*cosd(sza), label="MODTRAN, H₂O=2.05")#, xaxis=:log, yaxis=:log
        ylabel!(ylabels[1])
        xlabel!(xlabel)
        savefig_both("/home/sanghavi/data/EMIT_MODTRANcomp/azim180/plots/h2o_sensitivity/MODTRAN_4EMIT_azim180_IQU_tau_down_H2O.pdf")
    else
        plot(data1[:,1], -log.(data1[:,i])*cosd(vza), label="MODTRAN, H₂O=0.05")#, xaxis=:log, yaxis=:log)
        plot!(data11[:,1], -log.(data11[:,i])*cosd(vza), label="MODTRAN, H₂O=2.05")#, xaxis=:log, yaxis=:log)
        ylabel!(ylabels[2])
        xlabel!(xlabel)
        savefig_both("/home/sanghavi/data/EMIT_MODTRANcomp/azim180/plots/h2o_sensitivity/MODTRAN_4EMIT_azim180_IQU_tau_up_H2O.pdf")
    end
end
    
plot(data1[:,1], -log.(data1[:,4])*cosd(sza) + log.(data1[:,6])*cosd(vza), label="MODTRAN, H₂O=0.05")#, xaxis=:log, yaxis=:log)
plot!(data11[:,1], -log.(data11[:,4])*cosd(sza) + log.(data11[:,6])*cosd(vza), label="MODTRAN, H₂O=2.05")#, xaxis=:log, yaxis=:log)
ylabel!(L"$\Delta\tau_\mathrm{atm}$ [down-up]")
xlabel!(xlabel)
savefig_both("/home/sanghavi/data/EMIT_MODTRANcomp/azim180/plots/h2o_sensitivity/MODTRAN_4EMIT_azim180_IQU_tau_diff_H2O.pdf")

ind_uv1 = findall(x -> 350<=x<= 400.0, data1[:,1])
ind_vis1 = findall(x -> 500<=x<=700.0, data1[:,1])
ind_nir1 = findall(x -> 950<=x<=1100.0, data1[:,1])
ind_swir11 = findall(x -> 1500<=x<=1700.0, data1[:,1])
ind_swir21 = findall(x -> 2100<=x<=2200.0, data1[:,1])

spec_slice = ["uv", "vis", "nir", "swir1", "swir2"]

for ctr=1:5
    if ctr==1
        ind1 = ind_uv1
    elseif ctr==2
        ind1 = ind_vis1
    elseif ctr==3
        ind1 = ind_nir1
    elseif ctr==4
        ind1 = ind_swir11           
    elseif ctr==5
        ind1 = ind_swir21
    end
    if isempty(ind1)
        @info "Skipping spec_slice $(spec_slice[ctr]) (MODTRAN-only) — empty data1 index range"
        continue
    end
    local ylabels = [L"$\rho_\mathrm{atm}$", L"$T^\downarrow_\mathrm{dir}$", L"$T^\downarrow_\mathrm{dif}$", L"$T^\uparrow_\mathrm{dir}$", L"$T^\uparrow_\mathrm{dif}$", L"$S_\mathrm{atm}$"]
    local xlabel = "Wavelength [nm]"
    for i=3:8
        plot(data1[ind1,1], data1[ind1,i], label="MODTRAN, H₂O=0.05")#, xaxis=:log, yaxis=:log)
        plot!(data11[ind1,1], data11[ind1,i], label="MODTRAN, H₂O=2.05")#, xaxis=:log, yaxis=:log)
        xlabel!(xlabel)
        ylabel!(ylabels[i-2])
        savefig_both("/home/sanghavi/data/EMIT_MODTRANcomp/azim180/plots/h2o_sensitivity/MODTRAN_4EMIT_azim180_IQU_$(i-2)"*"_"*spec_slice[ctr]*"_H2O.pdf")
    end

    ylabels = [L"$\tau_\mathrm{atm}$ [downward]", L"$\tau_\mathrm{atm}$ [upward]"]
    for i in [4,6]
        if i==4
            plot(data1[ind1,1], -log.(data1[ind1,i])*cosd(sza), label="MODTRAN, H₂O=0.05")#, xaxis=:log, yaxis=:log)
            plot!(data11[ind1,1], -log.(data11[ind1,i])*cosd(sza), label="MODTRAN, H₂O=2.05")#, xaxis=:log, yaxis=:log)
            ylabel!(ylabels[1])
            xlabel!(xlabel)
            savefig_both("/home/sanghavi/data/EMIT_MODTRANcomp/azim180/plots/h2o_sensitivity/MODTRAN_4EMIT_azim180_IQU_tau_down_"*spec_slice[ctr]*"_H2O.pdf")
        else
            plot(data1[ind1,1], -log.(data1[ind1,i])*cosd(vza), label="MODTRAN, H₂O=0.05")#, xaxis=:log, yaxis=:log)
            plot!(data11[ind1,1], -log.(data11[ind1,i])*cosd(vza), label="MODTRAN, H₂O=2.05")#, xaxis=:log, yaxis=:log)
            ylabel!(ylabels[2])
            xlabel!(xlabel)
            savefig_both("/home/sanghavi/data/EMIT_MODTRANcomp/azim180/plots/h2o_sensitivity/MODTRAN_4EMIT_azim180_IQU_tau_up_"*spec_slice[ctr]*"_H2O.pdf")
        end
    end
    
    plot(data1[ind1,1], -log.(data1[ind1,4])*cosd(sza) + log.(data1[ind1,6])*cosd(vza), label="MODTRAN, H₂O=0.05")#, xaxis=:log, yaxis=:log)
    plot!(data11[ind1,1], -log.(data11[ind1,4])*cosd(sza) + log.(data11[ind1,6])*cosd(vza), label="MODTRAN, H₂O=2.05")#, xaxis=:log, yaxis=:log)
    ylabel!(L"$\Delta\tau_\mathrm{atm}$ [down-up]")
    xlabel!(xlabel)
    savefig_both("/home/sanghavi/data/EMIT_MODTRANcomp/azim180/plots/h2o_sensitivity/MODTRAN_4EMIT_azim180_IQU_tau_diff_"*spec_slice[ctr]*"_H2O.pdf")
end
end  # /azim180-disabled
