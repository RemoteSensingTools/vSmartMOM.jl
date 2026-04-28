using Plots
using DelimitedFiles
using LaTeXStrings

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

MOM_DIR = "/home/sanghavi/data/EMIT_MODTRANcomp/vSmartMOM_out/"
filename2(h2o, aot, gndalt, tsz) =
           "vSmartMOM_modtran_equiv_H2O$(fmt(h2o,4))_AOT$(fmt(aot,4))_GNDALT$(fmt(gndalt,3))_TSZ$(fmt(tsz,1)).dat"
path2 = joinpath(MOM_DIR, filename2(h2o, aot, gndalt, tsz))
# The vSmartMOM file has two '#'-prefixed header lines; skipping only one would
# leave the tokenised column-name line in data2 and splay sphalb into a phantom
# 9th column of empty strings.
data2 = readdlm(path2, comments=true, comment_char='#')

ylabels = [L"$\rho_\mathrm{atm}$", L"$T^\downarrow_\mathrm{dir}$", L"$T^\downarrow_\mathrm{dif}$", L"$T^\uparrow_\mathrm{dir}$", L"$T^\uparrow_\mathrm{dif}$", L"$S_\mathrm{atm}$"]
xlabel = L"Wavelength [nm]"
for i=3:8
    plot(data1[:,1], data1[:,i], label="MODTRAN")#, xaxis=:log, yaxis=:log)
    plot!(data2[:,1], data2[:,i], label="vSmartMOM")#, xaxis=:log, yaxis=:log)
    xlabel!(xlabel)
    ylabel!(ylabels[i-2])
    savefig("MODTRAN_vSmartMOM_4EMIT_$(i-2).pdf")
end

vza = 12.0
sza = 39.583965990615404
ylabels = [L"$\tau_\mathrm{atm}$ [downward]", L"$\tau_\mathrm{atm}$ [upward]"]
for i in [4,6]
    if i==4
        plot(-log.(data1[:,1])*cosd(sza), data1[:,i], label="MODTRAN")#, xaxis=:log, yaxis=:log)
        plot!(-log.(data2[:,1])*cosd(sza), data2[:,i], label="vSmartMOM")#, xaxis=:log, yaxis=:log)
        ylabel!(ylabels[1])
        xlabel!(xlabel)
        savefig("MODTRAN_vSmartMOM_4EMIT_tau_down.pdf")
    else
        plot(-log.(data1[:,1])*cosd(vza), data1[:,i], label="MODTRAN")#, xaxis=:log, yaxis=:log)
        plot!(-log.(data2[:,1])*cosd(vza), data2[:,i], label="vSmartMOM")#, xaxis=:log, yaxis=:log)
        ylabel!(ylabels[2])
        xlabel!(xlabel)
        savefig("MODTRAN_vSmartMOM_4EMIT_tau_up.pdf")
    end
end
    
plot(data1[:,1], -log.(data1[:,4])*cosd(sza) + log.(data1[:,6])*cosd(vza), label="MODTRAN")#, xaxis=:log, yaxis=:log)
plot!(data2[:,1], -log.(data2[:,4])*cosd(sza) + log.(data2[:,6])*cosd(vza), label="vSmartMOM")#, xaxis=:log, yaxis=:log)
ylabel!(L"$\Delta\tau_\mathrm{atm}$ [down-up]")
xlabel!(xlabel)
savefig("MODTRAN_vSmartMOM_4EMIT_tau_diff.pdf")
