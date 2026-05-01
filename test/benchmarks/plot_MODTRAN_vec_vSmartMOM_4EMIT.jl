#=
plot_MODTRAN_vec_vSmartMOM_4EMIT.jl — vector-RT (Stokes_IQU) MODTRAN/vSmartMOM compare.

Reads:
  MODTRAN .dat from $HOME/data/EMIT_MODTRANcomp/MODTRAN_out/
  vSmartMOM .dat from $HOME/data/EMIT_MODTRANcomp/vSmartMOM_out/vectorIQU/

The vSmartMOM filenames expected here are `vSmartMOM_modtran_equiv_*.dat`
(MODTRAN-equivalent 6S decomposition: rho_atm, td_dir, td_dif, tu_dir, tu_dif,
sphalb) generated from IQU runs of `modtran_equivalent_fields.jl`. The
Stokes-I component of vSmartMOM is compared against MODTRAN here; Q/U are
not consumed (MODTRAN is scalar). Output filenames carry the `_IQU` suffix
to distinguish them from the scalar plotter.

For the scalar (Stokes_I) equivalent see plot_MODTRAN_vSmartMOM_4EMIT.jl.
For a vSmartMOM scalar-vs-vector intercomparison see plot_scalar_vs_vector_vSmartMOM_4EMIT.jl.
=#
using Plots
using DelimitedFiles
using LaTeXStrings

default(alpha = 0.5)

"""
    savefig_both(name)

Save the current plot at `name` (linear y-axis) and also at
`<dir>/logplots/<stem>_log<ext>` with `yscale=:log10`. Restores the
y-axis scale to linear afterwards so subsequent `plot!` calls are unaffected.
"""
function savefig_both(name::AbstractString)
    dir = dirname(name)
    isempty(dir) && (dir = ".")
    mkpath(dir)
    Plots.savefig(name)
    log_dir = joinpath(dir, "logplots")
    mkpath(log_dir)
    stem, ext = splitext(basename(name))
    yaxis!(:log10)
    Plots.savefig(joinpath(log_dir, stem * "_log" * ext))
    yaxis!(:identity)
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

MOM_DIR = "/home/sanghavi/data/EMIT_MODTRANcomp/vSmartMOM_out/vectorIQU/"
filename2(h2o, aot, gndalt, tsz) =
           "vSmartMOM_modtran_equiv_H2O$(fmt(h2o,4))_AOT$(fmt(aot,4))_GNDALT$(fmt(gndalt,3))_TSZ$(fmt(tsz,1)).dat"
path2 = joinpath(MOM_DIR, filename2(h2o, aot, gndalt, tsz))
# The vSmartMOM file has two '#'-prefixed header lines; skipping only one would
# leave the tokenised column-name line in data2 and splay sphalb into a phantom
# 9th column of empty strings.
data2 = readdlm(path2, comments=true, comment_char='#')

ylabels = [L"$\rho_\mathrm{atm}$", L"$T^\downarrow_\mathrm{dir}$", L"$T^\downarrow_\mathrm{dif}$", L"$T^\uparrow_\mathrm{dir}$", L"$T^\uparrow_\mathrm{dif}$", L"$S_\mathrm{atm}$"]
xlabel = "Wavelength [nm]"
for i=3:8
    plot(data1[:,1], data1[:,i], label="MODTRAN")#, xaxis=:log, yaxis=:log)
    plot!(data2[:,1], data2[:,i], label="vSmartMOM IQU")#, xaxis=:log, yaxis=:log)
    xlabel!(xlabel)
    ylabel!(ylabels[i-2])
    savefig_both("MODTRAN_vSmartMOM_4EMIT_IQU_$(i-2).pdf")
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
        savefig_both("MODTRAN_vSmartMOM_4EMIT_IQU_tau_down.pdf")
    else
        plot(data1[:,1], -log.(data1[:,i])*cosd(vza), label="MODTRAN")#, xaxis=:log, yaxis=:log)
        plot!(data2[:,1], -log.(data2[:,i])*cosd(vza), label="vSmartMOM IQU")#, xaxis=:log, yaxis=:log)
        ylabel!(ylabels[2])
        xlabel!(xlabel)
        savefig_both("MODTRAN_vSmartMOM_4EMIT_IQU_tau_up.pdf")
    end
end
    
plot(data1[:,1], -log.(data1[:,4])*cosd(sza) + log.(data1[:,6])*cosd(vza), label="MODTRAN")#, xaxis=:log, yaxis=:log)
plot!(data2[:,1], -log.(data2[:,4])*cosd(sza) + log.(data2[:,6])*cosd(vza), label="vSmartMOM IQU")#, xaxis=:log, yaxis=:log)
ylabel!(L"$\Delta\tau_\mathrm{atm}$ [down-up]")
xlabel!(xlabel)
savefig_both("MODTRAN_vSmartMOM_4EMIT_IQU_tau_diff.pdf")

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
    ylabels = [L"$\rho_\mathrm{atm}$", L"$T^\downarrow_\mathrm{dir}$", L"$T^\downarrow_\mathrm{dif}$", L"$T^\uparrow_\mathrm{dir}$", L"$T^\uparrow_\mathrm{dif}$", L"$S_\mathrm{atm}$"]
    xlabel = "Wavelength [nm]"
    for i=3:8
        plot(data1[ind1,1], data1[ind1,i], label="MODTRAN")#, xaxis=:log, yaxis=:log)
        plot!(data2[ind2,1], data2[ind2,i], label="vSmartMOM IQU")#, xaxis=:log, yaxis=:log)
        xlabel!(xlabel)
        ylabel!(ylabels[i-2])
        savefig_both("/home/sanghavi/data/EMIT_MODTRANcomp/plots/MODTRAN_vSmartMOM_4EMIT_IQU_$(i-2)"*"_"*spec_slice[ctr]*".pdf")
    end

    ylabels = [L"$\tau_\mathrm{atm}$ [downward]", L"$\tau_\mathrm{atm}$ [upward]"]
    for i in [4,6]
        if i==4
            plot(data1[ind1,1], -log.(data1[ind1,i])*cosd(sza), label="MODTRAN")#, xaxis=:log, yaxis=:log)
            plot!(data2[ind2,1], -log.(data2[ind2,i])*cosd(sza), label="vSmartMOM IQU")#, xaxis=:log, yaxis=:log)
            ylabel!(ylabels[1])
            xlabel!(xlabel)
            savefig_both("/home/sanghavi/data/EMIT_MODTRANcomp/plots/MODTRAN_vSmartMOM_4EMIT_IQU_tau_down_"*spec_slice[ctr]*".pdf")
        else
            plot(data1[ind1,1], -log.(data1[ind1,i])*cosd(vza), label="MODTRAN")#, xaxis=:log, yaxis=:log)
            plot!(data2[ind2,1], -log.(data2[ind2,i])*cosd(vza), label="vSmartMOM IQU")#, xaxis=:log, yaxis=:log)
            ylabel!(ylabels[2])
            xlabel!(xlabel)
            savefig_both("/home/sanghavi/data/EMIT_MODTRANcomp/plots/MODTRAN_vSmartMOM_4EMIT_IQU_tau_up_"*spec_slice[ctr]*".pdf")
        end
    end
    
    plot(data1[ind1,1], -log.(data1[ind1,4])*cosd(sza) + log.(data1[ind1,6])*cosd(vza), label="MODTRAN")#, xaxis=:log, yaxis=:log)
    plot!(data2[ind2,1], -log.(data2[ind2,4])*cosd(sza) + log.(data2[ind2,6])*cosd(vza), label="vSmartMOM IQU")#, xaxis=:log, yaxis=:log)
    ylabel!(L"$\Delta\tau_\mathrm{atm}$ [down-up]")
    xlabel!(xlabel)
    savefig_both("/home/sanghavi/data/EMIT_MODTRANcomp/plots/MODTRAN_vSmartMOM_4EMIT_IQU_tau_diff_"*spec_slice[ctr]*".pdf")
end

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

MOM_DIR = "/home/sanghavi/data/EMIT_MODTRANcomp/vSmartMOM_out/vectorIQU/"
filename21(h2o, aot, gndalt, tsz) =
           "vSmartMOM_modtran_equiv_H2O$(fmt(h2o,4))_AOT$(fmt(aot,4))_GNDALT$(fmt(gndalt,3))_TSZ$(fmt(tsz,1)).dat"
path21 = joinpath(MOM_DIR, filename21(h2o, aot, gndalt, tsz))
# The vSmartMOM file has two '#'-prefixed header lines; skipping only one would
# leave the tokenised column-name line in data2 and splay sphalb into a phantom
# 9th column of empty strings.
data21 = readdlm(path21, comments=true, comment_char='#')

ylabels = [L"$\rho_\mathrm{atm}$", L"$T^\downarrow_\mathrm{dir}$", L"$T^\downarrow_\mathrm{dif}$", L"$T^\uparrow_\mathrm{dir}$", L"$T^\uparrow_\mathrm{dif}$", L"$S_\mathrm{atm}$"]
xlabel = "Wavelength [nm]"
for i=3:8
    plot(data11[:,1], data11[:,i], label="MODTRAN, H₂O=0.05")#, xaxis=:log, yaxis=:log)
    plot!(data21[:,1], data21[:,i], label="vSmartMOM IQU, H₂O=0.05")#, xaxis=:log, yaxis=:log)
    plot!(data11[:,1], data11[:,i], label="MODTRAN, H₂O=2.05")#, xaxis=:log, yaxis=:log)
    plot!(data21[:,1], data21[:,i], label="vSmartMOM IQU, H₂O=2.05")
    xlabel!(xlabel)
    ylabel!(ylabels[i-2])
    savefig_both("/home/sanghavi/data/EMIT_MODTRANcomp/plots/MODTRAN_vSmartMOM_4EMIT_IQU_$(i-2)_H2O.pdf")
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
        savefig_both("/home/sanghavi/data/EMIT_MODTRANcomp/plots/MODTRAN_vSmartMOM_4EMIT_IQU_tau_down_H2O.pdf")
    else
        plot(data1[:,1], -log.(data1[:,i])*cosd(vza), label="MODTRAN, H₂O=0.05")#, xaxis=:log, yaxis=:log)
        plot!(data2[:,1], -log.(data2[:,i])*cosd(vza), label="vSmartMOM IQU, H₂O=0.05")#, xaxis=:log, yaxis=:log)
        plot!(data11[:,1], -log.(data11[:,i])*cosd(vza), label="MODTRAN, H₂O=2.05")#, xaxis=:log, yaxis=:log)
        plot!(data21[:,1], -log.(data21[:,i])*cosd(vza), label="vSmartMOM IQU, H₂O=2.05")#, xaxis=:log, yaxis=:log)
        ylabel!(ylabels[2])
        xlabel!(xlabel)
        savefig_both("/home/sanghavi/data/EMIT_MODTRANcomp/plots/MODTRAN_vSmartMOM_4EMIT_IQU_tau_up_H2O.pdf")
    end
end
    
plot(data1[:,1], -log.(data1[:,4])*cosd(sza) + log.(data1[:,6])*cosd(vza), label="MODTRAN, H₂O=0.05")#, xaxis=:log, yaxis=:log)
plot!(data2[:,1], -log.(data2[:,4])*cosd(sza) + log.(data2[:,6])*cosd(vza), label="vSmartMOM IQU, H₂O=0.05")#, xaxis=:log, yaxis=:log)
plot!(data11[:,1], -log.(data11[:,4])*cosd(sza) + log.(data11[:,6])*cosd(vza), label="MODTRAN, H₂O=2.05")#, xaxis=:log, yaxis=:log)
plot!(data21[:,1], -log.(data21[:,4])*cosd(sza) + log.(data21[:,6])*cosd(vza), label="vSmartMOM IQU, H₂O=2.05")#, xaxis=:log, yaxis=:log)
ylabel!(L"$\Delta\tau_\mathrm{atm}$ [down-up]")
xlabel!(xlabel)
savefig_both("/home/sanghavi/data/EMIT_MODTRANcomp/plots/MODTRAN_vSmartMOM_4EMIT_IQU_tau_diff_H2O.pdf")

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
    ylabels = [L"$\rho_\mathrm{atm}$", L"$T^\downarrow_\mathrm{dir}$", L"$T^\downarrow_\mathrm{dif}$", L"$T^\uparrow_\mathrm{dir}$", L"$T^\uparrow_\mathrm{dif}$", L"$S_\mathrm{atm}$"]
    xlabel = "Wavelength [nm]"
    for i=3:8
        plot(data1[ind1,1], data1[ind1,i], label="MODTRAN, H₂O=0.05")#, xaxis=:log, yaxis=:log)
        plot!(data2[ind2,1], data2[ind2,i], label="vSmartMOM IQU, H₂O=0.05")#, xaxis=:log, yaxis=:log)
        plot!(data11[ind1,1], data11[ind1,i], label="MODTRAN, H₂O=2.05")#, xaxis=:log, yaxis=:log)
        plot!(data21[ind2,1], data21[ind2,i], label="vSmartMOM IQU, H₂O=2.05")#, xaxis=:log, yaxis=:log)
        xlabel!(xlabel)
        ylabel!(ylabels[i-2])
        savefig_both("/home/sanghavi/data/EMIT_MODTRANcomp/plots/MODTRAN_vSmartMOM_4EMIT_IQU_$(i-2)"*"_"*spec_slice[ctr]*"_H2O.pdf")
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
            savefig_both("/home/sanghavi/data/EMIT_MODTRANcomp/plots/MODTRAN_vSmartMOM_4EMIT_IQU_tau_down_"*spec_slice[ctr]*"_H2O.pdf")
        else
            plot(data1[ind1,1], -log.(data1[ind1,i])*cosd(vza), label="MODTRAN, H₂O=0.05")#, xaxis=:log, yaxis=:log)
            plot!(data2[ind2,1], -log.(data2[ind2,i])*cosd(vza), label="vSmartMOM IQU, H₂O=0.05")#, xaxis=:log, yaxis=:log)
            plot!(data11[ind1,1], -log.(data11[ind1,i])*cosd(vza), label="MODTRAN, H₂O=2.05")#, xaxis=:log, yaxis=:log)
            plot!(data21[ind2,1], -log.(data21[ind2,i])*cosd(vza), label="vSmartMOM IQU, H₂O=2.05")#, xaxis=:log, yaxis=:log)
            ylabel!(ylabels[2])
            xlabel!(xlabel)
            savefig_both("/home/sanghavi/data/EMIT_MODTRANcomp/plots/MODTRAN_vSmartMOM_4EMIT_IQU_tau_up_"*spec_slice[ctr]*"_H2O.pdf")
        end
    end
    
    plot(data1[ind1,1], -log.(data1[ind1,4])*cosd(sza) + log.(data1[ind1,6])*cosd(vza), label="MODTRAN, H₂O=0.05")#, xaxis=:log, yaxis=:log)
    plot!(data2[ind2,1], -log.(data2[ind2,4])*cosd(sza) + log.(data2[ind2,6])*cosd(vza), label="vSmartMOM IQU, H₂O=0.05")#, xaxis=:log, yaxis=:log)
    plot!(data11[ind1,1], -log.(data11[ind1,4])*cosd(sza) + log.(data11[ind1,6])*cosd(vza), label="MODTRAN, H₂O=2.05")#, xaxis=:log, yaxis=:log)
    plot!(data21[ind2,1], -log.(data21[ind2,4])*cosd(sza) + log.(data21[ind2,6])*cosd(vza), label="vSmartMOM IQU, H₂O=2.05")#, xaxis=:log, yaxis=:log)
    ylabel!(L"$\Delta\tau_\mathrm{atm}$ [down-up]")
    xlabel!(xlabel)
    savefig_both("/home/sanghavi/data/EMIT_MODTRANcomp/plots/MODTRAN_vSmartMOM_4EMIT_IQU_tau_diff_"*spec_slice[ctr]*"_H2O.pdf")
end

#====== only vSmartMOM ======#
ylabels = [L"$\rho_\mathrm{atm}$", L"$T^\downarrow_\mathrm{dir}$", L"$T^\downarrow_\mathrm{dif}$", L"$T^\uparrow_\mathrm{dir}$", L"$T^\uparrow_\mathrm{dif}$", L"$S_\mathrm{atm}$"]
xlabel = "Wavelength [nm]"
for i=3:8
    plot(data21[:,1], data21[:,i], label="vSmartMOM IQU, H₂O=0.05")#, xaxis=:log, yaxis=:log)
    plot!(data21[:,1], data21[:,i], label="vSmartMOM IQU, H₂O=2.05")
    xlabel!(xlabel)
    ylabel!(ylabels[i-2])
    savefig_both("/home/sanghavi/data/EMIT_MODTRANcomp/plots/h2o_sensitivity/vSmartMOM_4EMIT_IQU_$(i-2)_H2O.pdf")
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
        savefig_both("/home/sanghavi/data/EMIT_MODTRANcomp/plots/h2o_sensitivity/vSmartMOM_4EMIT_IQU_tau_down_H2O.pdf")
    else
        plot(data2[:,1], -log.(data2[:,i])*cosd(vza), label="vSmartMOM IQU, H₂O=0.05")#, xaxis=:log, yaxis=:log)
        plot!(data21[:,1], -log.(data21[:,i])*cosd(vza), label="vSmartMOM IQU, H₂O=2.05")#, xaxis=:log, yaxis=:log)
        ylabel!(ylabels[2])
        xlabel!(xlabel)
        savefig_both("/home/sanghavi/data/EMIT_MODTRANcomp/plots/h2o_sensitivity/vSmartMOM_4EMIT_IQU_tau_up_H2O.pdf")
    end
end
    
plot(data2[:,1], -log.(data2[:,4])*cosd(sza) + log.(data2[:,6])*cosd(vza), label="vSmartMOM IQU, H₂O=0.05")#, 
plot!(data21[:,1], -log.(data21[:,4])*cosd(sza) + log.(data21[:,6])*cosd(vza), label="vSmartMOM IQU, H₂O=2.05")#, xaxis=:log, yaxis=:log)
ylabel!(L"$\Delta\tau_\mathrm{atm}$ [down-up]")
xlabel!(xlabel)
savefig_both("/home/sanghavi/data/EMIT_MODTRANcomp/plots/h2o_sensitivity/vSmartMOM_4EMIT_IQU_tau_diff_H2O.pdf")

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
    ylabels = [L"$\rho_\mathrm{atm}$", L"$T^\downarrow_\mathrm{dir}$", L"$T^\downarrow_\mathrm{dif}$", L"$T^\uparrow_\mathrm{dir}$", L"$T^\uparrow_\mathrm{dif}$", L"$S_\mathrm{atm}$"]
    xlabel = "Wavelength [nm]"
    for i=3:8
        plot(data2[ind2,1], data2[ind2,i], label="vSmartMOM IQU, H₂O=0.05")#, xaxis=:log, yaxis=:log)
        plot!(data21[ind2,1], data21[ind2,i], label="vSmartMOM IQU, H₂O=2.05")#, xaxis=:log, yaxis=:log)
        xlabel!(xlabel)
        ylabel!(ylabels[i-2])
        savefig_both("/home/sanghavi/data/EMIT_MODTRANcomp/plots/h2o_sensitivity/vSmartMOM_4EMIT_IQU_$(i-2)"*"_"*spec_slice[ctr]*"_H2O.pdf")
    end

    ylabels = [L"$\tau_\mathrm{atm}$ [downward]", L"$\tau_\mathrm{atm}$ [upward]"]
    for i in [4,6]
        if i==4
            plot(data2[ind2,1], -log.(data2[ind2,i])*cosd(sza), label="vSmartMOM IQU, H₂O=0.05")#, xaxis=:log, yaxis=:log)
            plot!(data21[ind2,1], -log.(data21[ind2,i])*cosd(sza), label="vSmartMOM IQU, H₂O=2.05")#, xaxis=:log, yaxis=:log)
            ylabel!(ylabels[1])
            xlabel!(xlabel)
            savefig_both("/home/sanghavi/data/EMIT_MODTRANcomp/plots/h2o_sensitivity/vSmartMOM_4EMIT_IQU_tau_down_"*spec_slice[ctr]*"_H2O.pdf")
        else
            plot(data2[ind2,1], -log.(data2[ind2,i])*cosd(vza), label="vSmartMOM IQU, H₂O=0.05")#, xaxis=:log, yaxis=:log)
            plot!(data21[ind2,1], -log.(data21[ind2,i])*cosd(vza), label="vSmartMOM IQU, H₂O=2.05")#, xaxis=:log, yaxis=:log)
            ylabel!(ylabels[2])
            xlabel!(xlabel)
            savefig_both("/home/sanghavi/data/EMIT_MODTRANcomp/plots/h2o_sensitivity/vSmartMOM_4EMIT_IQU_tau_up_"*spec_slice[ctr]*"_H2O.pdf")
        end
    end
    
    plot(data2[ind2,1], -log.(data2[ind2,4])*cosd(sza) + log.(data2[ind2,6])*cosd(vza), label="vSmartMOM IQU, H₂O=0.05")#, xaxis=:log, yaxis=:log)
    plot!(data21[ind2,1], -log.(data21[ind2,4])*cosd(sza) + log.(data21[ind2,6])*cosd(vza), label="vSmartMOM IQU, H₂O=2.05")#, xaxis=:log, yaxis=:log)
    ylabel!(L"$\Delta\tau_\mathrm{atm}$ [down-up]")
    xlabel!(xlabel)
    savefig_both("/home/sanghavi/data/EMIT_MODTRANcomp/plots/h2o_sensitivity/vSmartMOM_4EMIT_IQU_tau_diff_"*spec_slice[ctr]*"_H2O.pdf")
end
#====== only MODTRAN ======#
ylabels = [L"$\rho_\mathrm{atm}$", L"$T^\downarrow_\mathrm{dir}$", L"$T^\downarrow_\mathrm{dif}$", L"$T^\uparrow_\mathrm{dir}$", L"$T^\uparrow_\mathrm{dif}$", L"$S_\mathrm{atm}$"]
xlabel = "Wavelength [nm]"
for i=3:8
    plot(data11[:,1], data11[:,i], label="MODTRAN, H₂O=0.05")#, xaxis=:log, yaxis=:log)
    plot!(data11[:,1], data11[:,i], label="MODTRAN, H₂O=2.05")#, xaxis=:log, yaxis=:log)

    xlabel!(xlabel)
    ylabel!(ylabels[i-2])
    savefig_both("/home/sanghavi/data/EMIT_MODTRANcomp/plots/h2o_sensitivity/MODTRAN_4EMIT_IQU_$(i-2)_H2O.pdf")
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
        savefig_both("/home/sanghavi/data/EMIT_MODTRANcomp/plots/h2o_sensitivity/MODTRAN_4EMIT_IQU_tau_down_H2O.pdf")
    else
        plot(data1[:,1], -log.(data1[:,i])*cosd(vza), label="MODTRAN, H₂O=0.05")#, xaxis=:log, yaxis=:log)
        plot!(data11[:,1], -log.(data11[:,i])*cosd(vza), label="MODTRAN, H₂O=2.05")#, xaxis=:log, yaxis=:log)
        ylabel!(ylabels[2])
        xlabel!(xlabel)
        savefig_both("/home/sanghavi/data/EMIT_MODTRANcomp/plots/h2o_sensitivity/MODTRAN_4EMIT_IQU_tau_up_H2O.pdf")
    end
end
    
plot(data1[:,1], -log.(data1[:,4])*cosd(sza) + log.(data1[:,6])*cosd(vza), label="MODTRAN, H₂O=0.05")#, xaxis=:log, yaxis=:log)
plot!(data11[:,1], -log.(data11[:,4])*cosd(sza) + log.(data11[:,6])*cosd(vza), label="MODTRAN, H₂O=2.05")#, xaxis=:log, yaxis=:log)
ylabel!(L"$\Delta\tau_\mathrm{atm}$ [down-up]")
xlabel!(xlabel)
savefig_both("/home/sanghavi/data/EMIT_MODTRANcomp/plots/h2o_sensitivity/MODTRAN_4EMIT_IQU_tau_diff_H2O.pdf")

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
    ylabels = [L"$\rho_\mathrm{atm}$", L"$T^\downarrow_\mathrm{dir}$", L"$T^\downarrow_\mathrm{dif}$", L"$T^\uparrow_\mathrm{dir}$", L"$T^\uparrow_\mathrm{dif}$", L"$S_\mathrm{atm}$"]
    xlabel = "Wavelength [nm]"
    for i=3:8
        plot(data1[ind1,1], data1[ind1,i], label="MODTRAN, H₂O=0.05")#, xaxis=:log, yaxis=:log)
        plot!(data11[ind1,1], data11[ind1,i], label="MODTRAN, H₂O=2.05")#, xaxis=:log, yaxis=:log)
        xlabel!(xlabel)
        ylabel!(ylabels[i-2])
        savefig_both("/home/sanghavi/data/EMIT_MODTRANcomp/plots/h2o_sensitivity/MODTRAN_4EMIT_IQU_$(i-2)"*"_"*spec_slice[ctr]*"_H2O.pdf")
    end

    ylabels = [L"$\tau_\mathrm{atm}$ [downward]", L"$\tau_\mathrm{atm}$ [upward]"]
    for i in [4,6]
        if i==4
            plot(data1[ind1,1], -log.(data1[ind1,i])*cosd(sza), label="MODTRAN, H₂O=0.05")#, xaxis=:log, yaxis=:log)
            plot!(data11[ind1,1], -log.(data11[ind1,i])*cosd(sza), label="MODTRAN, H₂O=2.05")#, xaxis=:log, yaxis=:log)
            ylabel!(ylabels[1])
            xlabel!(xlabel)
            savefig_both("/home/sanghavi/data/EMIT_MODTRANcomp/plots/h2o_sensitivity/MODTRAN_4EMIT_IQU_tau_down_"*spec_slice[ctr]*"_H2O.pdf")
        else
            plot(data1[ind1,1], -log.(data1[ind1,i])*cosd(vza), label="MODTRAN, H₂O=0.05")#, xaxis=:log, yaxis=:log)
            plot!(data11[ind1,1], -log.(data11[ind1,i])*cosd(vza), label="MODTRAN, H₂O=2.05")#, xaxis=:log, yaxis=:log)
            ylabel!(ylabels[2])
            xlabel!(xlabel)
            savefig_both("/home/sanghavi/data/EMIT_MODTRANcomp/plots/h2o_sensitivity/MODTRAN_4EMIT_IQU_tau_up_"*spec_slice[ctr]*"_H2O.pdf")
        end
    end
    
    plot(data1[ind1,1], -log.(data1[ind1,4])*cosd(sza) + log.(data1[ind1,6])*cosd(vza), label="MODTRAN, H₂O=0.05")#, xaxis=:log, yaxis=:log)
    plot!(data11[ind1,1], -log.(data11[ind1,4])*cosd(sza) + log.(data11[ind1,6])*cosd(vza), label="MODTRAN, H₂O=2.05")#, xaxis=:log, yaxis=:log)
    ylabel!(L"$\Delta\tau_\mathrm{atm}$ [down-up]")
    xlabel!(xlabel)
    savefig_both("/home/sanghavi/data/EMIT_MODTRANcomp/plots/h2o_sensitivity/MODTRAN_4EMIT_IQU_tau_diff_"*spec_slice[ctr]*"_H2O.pdf")
end