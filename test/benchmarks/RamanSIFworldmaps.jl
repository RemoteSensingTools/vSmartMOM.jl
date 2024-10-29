using Plots 
using LegendrePolynomials
using DelimitedFiles
using Insolation, Dates
import Insolation.Parameters as IP
import ClimaParams as CP
using NCDatasets
using Interpolations
using Statistics
using InstrumentOperator
using LaTeXStrings
#using PlotlyJS

# Use the following data for wide spectrum reflectances 
# https://speclib.jpl.nasa.gov/library

dpath = "/net/squid/data3/data/FluoData1/group/oco2/"
L1File   = dpath*"L1bSc/oco2_L1bScND_26782a_190715_B10003r_200429202458.h5"
metFile  = dpath*"L2Met/oco2_L2MetND_26782a_190715_B10003r_200429202458.h5"
dictFile = "/home/cfranken/code/gitHub/InstrumentOperator.jl/json/oco2.yaml"
# Load L1 file (could just use filenames here as well)
oco = InstrumentOperator.load_L1(dictFile,L1File, metFile);
oco = InstrumentOperator.load_L1(dictFile,L1File, metFile);
# Pick some bands as tuple (or just one)
bands=(1,); #bands = (1,2,3);
# Indices within that band:
#indices = (92:885,114:845,50:916);
indices = (2:1015,); #indices = (114:845,);
# Geo Index (footprint,sounding):
GeoInd = [5,5000];
# Get data for that sounding:
oco_sounding = InstrumentOperator.getMeasurement(oco, bands, indices, GeoInd);

ocoMM=oco_sounding.mueller
iBand=1
res = 0.001e-3;
off = 0.5e-3
oco_wl = oco_sounding.ils[iBand].ν_out[1]-off:res:oco_sounding.ils[iBand].ν_out[end]+off;
#oco_wl *= 1e3

#Constants for later OCO A-band noise computation
MaxMS = 7.0e20 # photons/m^2/sr/μm/s
c_bkg = 0.0042
c_pht = 0.0089
h_Pl  = 6.62607015e-34 # J.s
c_l   = 299792458 # m/s


fit_window1 = [758.0, 759.2] 
fit_window2 = [770.0, 770.25] #[770.4, 771.6]
#sza_str   = "50"
#psurf_str = "1000"
#alb_str   = "0p2"
FT = Float32
a=[]
for i=0:20
    push!(a, acosd(i/20))  
end
sza=reverse(Int.(ceil.(a[8:21])))
ρ = zeros(FT,15)
ρ_str = []
for iρ = 1:15
    ρ[iρ] = (iρ-1)*0.05
    push!(ρ_str, replace(string(round(ρ[iρ], digits=2)),"."=>"p"))
end
psurf=[1000 750 500]
psurf=reverse(psurf')
#xnSIF_noRS = zeros(length(psurf), length(ρ_str), length(sza))

#coeff = zeros(4,length(psurf), length(ρ_str), length(sza))
#coeff_noRS = zeros(4,length(psurf), length(ρ_str), length(sza))
#coeffa = zeros(5,length(psurf), length(ρ_str), length(sza))
#coeffa_noRS = zeros(5,length(psurf), length(ρ_str), length(sza))

# Find indices for wavelength window:
# Determine ind: start
isurf = 3
iρ = 1
iA = 1
sza_str  = string(sza[iA])
alb_str  = ρ_str[iρ]
psurf_str= string(psurf[isurf])
fname0   = "/home/sanghavi/RamanSIFgrid/raylSIF_sza"*sza_str*"_alb"*alb_str*"_psurf"*psurf_str*"hpa_nors_ABO2.dat"
specNoRS = readdlm(fname0) 
specNoRS = specNoRS[end:-1:1,:]
wl_ = 1e7./specNoRS[:,1]
ind_to_elim = findall(x->x>1e-6, wl_[1:end-2]-wl_[3:end]) #indices of duplicate entries (to be eliminated)
wl = wl_[ filter(x->!(x in ind_to_elim), eachindex(wl_)) ]

ind1 = findall(x->x>fit_window1[1]*1e-3 && x<fit_window1[2]*1e-3, oco_sounding.SpectralGrid);         
ind2 = findall(x->x>fit_window2[1]*1e-3 && x<fit_window2[2]*1e-3, oco_sounding.SpectralGrid);         
ind1_hires = findall(x->x>fit_window1[1] && x<fit_window1[2], wl);         
ind2_hires = findall(x->x>fit_window2[1] && x<fit_window2[2], wl);         

# Determine ind: end

#Tabulate inelastic components
ieI1 = zeros(length(psurf), length(ρ_str), length(sza), length(ind1))
ieI2 = zeros(length(psurf), length(ρ_str), length(sza), length(ind2))
ieI1_hires = zeros(length(psurf), length(ρ_str), length(sza), length(ind1_hires))
ieI2_hires = zeros(length(psurf), length(ρ_str), length(sza), length(ind2_hires))

for isurf = 1:3 #3 # 1:1 #2:2 # 
    for iρ = 1:15 #3 #1:15
        for iA = 1:14
            sza_str  = string(sza[iA])
            alb_str  = ρ_str[iρ]
            psurf_str= string(psurf[isurf])

            #fname0   = "/home/sanghavi/RamanSIFgrid/raylSIF_sza"*sza_str*"_alb"*alb_str*"_psurf"*psurf_str*"hpa_nors_ABO2.dat"
            fname1   = "/home/sanghavi/RamanSIFgrid/raylSIF_sza"*sza_str*"_alb"*alb_str*"_psurf"*psurf_str*"hpa_rrs_ABO2.dat"

            #fname0_0   = "/home/sanghavi/RamanSIFgrid/rayl_sza"*sza_str*"_alb"*alb_str*"_psurf"*psurf_str*"hpa_nors_ABO2.dat"
            #fname1_0   = "/home/sanghavi/RamanSIFgrid/rayl_sza"*sza_str*"_alb"*alb_str*"_psurf"*psurf_str*"hpa_rrs_ABO2.dat"
            specRRS_  = readdlm(fname1)
            specRRS_ = specRRS_[end:-1:1,:]
            specRRS = zeros(length(wl), size(specRRS_,2))
            for i=1:size(specRRS_,2)
                specRRS[:,i] = specRRS_[ filter(x->!(x in ind_to_elim), eachindex(specRRS_[:,i])), i]
            end

            #specNoRS = readdlm(fname0)
            

            #ieI = specRRS[ filter(x->!(x in ind_to_elim), eachindex(specRRS[:,5])), 5]
            #ieQ = specRRS[ filter(x->!(x in ind_to_elim), eachindex(specRRS[:,6])), 6]
            #ieU = specRRS[ filter(x->!(x in ind_to_elim), eachindex(specRRS[:,7])), 7]
            ieI_tmp = (ocoMM[1]*specRRS[:,5] + 
                        ocoMM[2]*specRRS[:,6] +
                        ocoMM[3]*specRRS[:,7]).* 
                        (1e7./wl.^2) #mW/m^2/s/nm

            interp_I = LinearInterpolation((wl)*1e-3, (ieI_tmp));  
            # Re-interpolate I from ν_grid to new grid/resolution
            tmp1 = interp_I(oco_wl);
            # Convolve input spectrum with variable kernel
            tmp2 = InstrumentOperator.conv_spectra(oco_sounding.ils[iBand], oco_wl, tmp1)
            #specNoRS0 = readdlm(fname0_0)
            #specRRS0  = readdlm(fname1_0)

            #wl = 1e7./specRRS[:,1] # converting from wn [cm^{-1}] to wl [nm]
            #ind = findall(x->x>fit_window[1] && x<fit_window[2], wl);
            
            ieI1[isurf, iρ, iA, :] = tmp2[ind1] #ieI_tmp[ind1]
            ieI2[isurf, iρ, iA, :] = tmp2[ind2] #ieI_tmp[ind2]
            ieI1_hires[isurf, iρ, iA, :] = ieI_tmp[ind1_hires]
            ieI2_hires[isurf, iρ, iA, :] = ieI_tmp[ind2_hires]
        end
    end
end

l = @layout [a1 a2 a3; b1 b2 b3]
n=15;
CList = reshape( range(RGB(0.1, 0.1, 0.1), stop=RGB(0.1, 0.1, 0.1),length=n), 1, n );
p1 = plot(oco_sounding.SpectralGrid[ind1]*1e3, ieI1[3,:,7,:]',linecolor=CList, label="")
CList = reshape( range(colorant"blue", stop=colorant"red",length=n), 1, n );
p1 = plot!(wl[ind1_hires], ieI1_hires[3,:,7,:]',linecolor=CList, label="", ylabel=L"$I_\mathrm{ie}$", xtickfontsize=8, xrotation=45)

n=14;
CList = reshape( range(RGB(0.1, 0.1, 0.1), stop=RGB(0.1, 0.1, 0.1),length=n), 1, n );
p2 = plot(oco_sounding.SpectralGrid[ind1]*1e3, ieI1[3,7,:,:]',linecolor=CList, label="")
CList = reshape( range(colorant"red", stop=colorant"blue",length=n), 1, n );
p2 = plot!(wl[ind1_hires], ieI1_hires[3,7,:,:]',linecolor=CList, label="", xtickfontsize=8, xrotation=45)

n=3;
CList = reshape( range(RGB(0.1, 0.1, 0.1), stop=RGB(0.1, 0.1, 0.1),length=n), 1, n );
p3 = plot(oco_sounding.SpectralGrid[ind1]*1e3, ieI1[:,7,7,:]',linecolor=CList, label="")
CList = reshape( range(colorant"blue", stop=colorant"red",length=n), 1, n );
p3 = plot!(wl[ind1_hires], ieI1_hires[:,7,7,:]',linecolor=CList, label="", xtickfontsize=8, xrotation=45)

n=15;
CList = reshape( range(RGB(0.1, 0.1, 0.1), stop=RGB(0.1, 0.1, 0.1),length=n), 1, n );
q1 = plot(oco_sounding.SpectralGrid[ind2]*1e3, ieI2[3,:,7,:]',linecolor=CList, label="")
CList = reshape( range(colorant"blue", stop=colorant"red",length=n), 1, n );
q1 = plot!(wl[ind2_hires], ieI2_hires[3,:,7,:]',linecolor=CList, label="", xlabel=L"$\lambda\,$[nm]", ylabel=L"$I_\mathrm{ie}$", xtickfontsize=8, xrotation=45)

n=14;
CList = reshape( range(RGB(0.1, 0.1, 0.1), stop=RGB(0.1, 0.1, 0.1),length=n), 1, n );
q2 = plot(oco_sounding.SpectralGrid[ind2]*1e3, ieI2[3,7,:,:]',linecolor=CList, label="")
CList = reshape( range(colorant"red", stop=colorant"blue",length=n), 1, n );
q2 = plot!(wl[ind2_hires], ieI2_hires[3,7,:,:]',linecolor=CList, label="", xlabel=L"$\lambda\,$[nm]", xtickfontsize=8, xrotation=45)

n=3;
CList = reshape( range(RGB(0.1, 0.1, 0.1), stop=RGB(0.1, 0.1, 0.1),length=n), 1, n );
q3 = plot(oco_sounding.SpectralGrid[ind2]*1e3, ieI2[:,7,7,:]',linecolor=CList, label="")
CList = reshape( range(colorant"blue", stop=colorant"red",length=n), 1, n );
q3 = plot!(wl[ind2_hires], ieI2_hires[:,7,7,:]',linecolor=CList, label="", xlabel=L"$\lambda\,$[nm]", xtickfontsize=8, xrotation=45)
title1a = L"$\rho=0:0.05:0.7$"
title1b = L"$\mathrm{SZA}=46^\circ, p_\mathrm{surf}=1000\,\mathrm{hPa}$"
title2a = L"$\cos(\mathrm{SZA})\approx 0.35:0.05:1$"
title2b = L"$\rho=0.3, p_\mathrm{surf}=1000\,\mathrm{hPa}$" 
title3a = L"$p_\mathrm{surf}=500:250:1000\,\mathrm{hPa}$"
title3b = L"$\rho=0.3, \mathrm{SZA}=46^\circ$"


plot(p1, p2, p3, q1, q2, q3, layout = l, legend = false, title = [title1a*"\n"*title1b title2a*"\n"*title2b title3a*"\n"*title3b "" "" ""], titlefont = font(8))
savefig("/home/sanghavi/RamanSIFgrid/plots/ieI_wrt_rho_res.png")
#Determine scaled mean inelastic component
I_ie_mean = zeros(length(ind1_hires))
ocoI1_ie_mean = zeros(length(ind1))
ocoI2_ie_mean = zeros(length(ind2))
for isurf = 1:3 #3 # 1:1 #2:2 # 
    for iρ = 1:15 #3 #1:15
        for iA = 1:14
            sza_str  = string(sza[iA])
            alb_str  = ρ_str[iρ]
            psurf_str= string(psurf[isurf])

            #fname0   = "/home/sanghavi/RamanSIFgrid/raylSIF_sza"*sza_str*"_alb"*alb_str*"_psurf"*psurf_str*"hpa_nors_ABO2.dat"
            fname1   = "/home/sanghavi/RamanSIFgrid/raylSIF_sza"*sza_str*"_alb"*alb_str*"_psurf"*psurf_str*"hpa_rrs_ABO2.dat"

            #fname0_0   = "/home/sanghavi/RamanSIFgrid/rayl_sza"*sza_str*"_alb"*alb_str*"_psurf"*psurf_str*"hpa_nors_ABO2.dat"
            #fname1_0   = "/home/sanghavi/RamanSIFgrid/rayl_sza"*sza_str*"_alb"*alb_str*"_psurf"*psurf_str*"hpa_rrs_ABO2.dat"


            #specNoRS = readdlm(fname0)
            specRRS_  = readdlm(fname1)
            specRRS_ = specRRS_[end:-1:1,:]
            specRRS = zeros(length(wl), size(specRRS_,2))
            for i=1:size(specRRS_,2)
                specRRS[:,i] = specRRS_[ filter(x->!(x in ind_to_elim), eachindex(specRRS_[:,i])), i]
            end
            #specRRS = specRRS[end:-1:1,:]
            

            #ieI = specRRS[ filter(x->!(x in ind_to_elim), eachindex(specRRS[:,5])), 5]
            #ieQ = specRRS[ filter(x->!(x in ind_to_elim), eachindex(specRRS[:,6])), 6]
            #ieU = specRRS[ filter(x->!(x in ind_to_elim), eachindex(specRRS[:,7])), 7]
            ieI_tmp = (ocoMM[1]*specRRS[:,5] + 
                        ocoMM[2]*specRRS[:,6] +
                        ocoMM[3]*specRRS[:,7]).* 
                        (1e7./wl.^2) #mW/m^2/s/nm

            interp_I = LinearInterpolation((wl)*1e-3, (ieI_tmp));  
            # Re-interpolate I from ν_grid to new grid/resolution
            tmp1 = interp_I(oco_wl);
            # Convolve input spectrum with variable kernel
            tmp2 = InstrumentOperator.conv_spectra(oco_sounding.ils[iBand], oco_wl, tmp1)
            #specNoRS0 = readdlm(fname0_0)
            #specRRS0  = readdlm(fname1_0)

            #wl = 1e7./specRRS[:,1] # converting from wn [cm^{-1}] to wl [nm]
            #ind = findall(x->x>fit_window[1] && x<fit_window[2], wl);
            
            ocoI1_ie_mean .+= tmp2[ind1]./tmp2[ind1[1]] #ieI_tmp[ind1]
            ocoI2_ie_mean .+= tmp2[ind2]./tmp2[ind2[1]] #ieI_tmp[ind2]
            #ieI1_hires[isurf, iρ, iA, :] = ieI_tmp[ind1_hires]
            #ieI2_hires[isurf, iρ, iA, :] = ieI_tmp[ind2_hires]

            #specNoRS0 = readdlm(fname0_0)
            #specRRS0  = readdlm(fname1_0)

            #wl = 1e7./specRRS[:,1] # converting from wn [cm^{-1}] to wl [nm]
            #ind = findall(x->x>fit_window[1] && x<fit_window[2], wl);
            
            I_ie_mean .+= (specRRS[ind1_hires,5].* (1e7./wl[ind1_hires].^2))/(specRRS[ind1_hires[1],5].* (1e7./wl[ind1_hires[1]].^2)) #converting radiances from mW/m^2/sr/cm^{-1} to mW/m^2/sr/nm
        end
    end
end
I_ie_mean/=(3*15*14)
ocoI1_ie_mean/=(3*15*14)
ocoI2_ie_mean/=(3*15*14)

w1_xSIF0_corr = zeros(length(psurf), length(ρ_str), length(sza))
w1_xnSIF0_corr = zeros(length(psurf), length(ρ_str), length(sza))
w1_xSIF1_corr = zeros(length(psurf), length(ρ_str), length(sza))
w1_xnSIF1_corr = zeros(length(psurf), length(ρ_str), length(sza))

w1_xSIF0 = zeros(length(psurf), length(ρ_str), length(sza))
w1_xnSIF0 = zeros(length(psurf), length(ρ_str), length(sza))
w1_xSIF1 = zeros(length(psurf), length(ρ_str), length(sza))
w1_xnSIF1 = zeros(length(psurf), length(ρ_str), length(sza))

w1_xSIF0_corr_0 = zeros(length(psurf), length(ρ_str), length(sza))
w1_xnSIF0_corr_0 = zeros(length(psurf), length(ρ_str), length(sza))
w1_xSIF1_corr_0 = zeros(length(psurf), length(ρ_str), length(sza))
w1_xnSIF1_corr_0 = zeros(length(psurf), length(ρ_str), length(sza))

w1_xSIF0_0 = zeros(length(psurf), length(ρ_str), length(sza))
w1_xnSIF0_0 = zeros(length(psurf), length(ρ_str), length(sza))
w1_xSIF1_0 = zeros(length(psurf), length(ρ_str), length(sza))
w1_xnSIF1_0 = zeros(length(psurf), length(ρ_str), length(sza))

w2_xSIF0_corr = zeros(length(psurf), length(ρ_str), length(sza))
w2_xnSIF0_corr = zeros(length(psurf), length(ρ_str), length(sza))
w2_xSIF1_corr = zeros(length(psurf), length(ρ_str), length(sza))
w2_xnSIF1_corr = zeros(length(psurf), length(ρ_str), length(sza))

w2_xSIF0 = zeros(length(psurf), length(ρ_str), length(sza))
w2_xnSIF0 = zeros(length(psurf), length(ρ_str), length(sza))
w2_xSIF1 = zeros(length(psurf), length(ρ_str), length(sza))
w2_xnSIF1 = zeros(length(psurf), length(ρ_str), length(sza))

w2_xSIF0_corr_0 = zeros(length(psurf), length(ρ_str), length(sza))
w2_xnSIF0_corr_0 = zeros(length(psurf), length(ρ_str), length(sza))
w2_xSIF1_corr_0 = zeros(length(psurf), length(ρ_str), length(sza))
w2_xnSIF1_corr_0 = zeros(length(psurf), length(ρ_str), length(sza))

w2_xSIF0_0 = zeros(length(psurf), length(ρ_str), length(sza))
w2_xnSIF0_0 = zeros(length(psurf), length(ρ_str), length(sza))
w2_xSIF1_0 = zeros(length(psurf), length(ρ_str), length(sza))
w2_xnSIF1_0 = zeros(length(psurf), length(ρ_str), length(sza))

coAdd_N = 36
for isurf = 1:3 # 1:1 #2:2 # 
    for iρ = 1:15 #3 #1:15
        for iA = 1:14
            sza_str  = string(sza[iA])
            alb_str  = ρ_str[iρ]
            psurf_str= string(psurf[isurf])

            #fname0   = "/home/sanghavi/RamanSIFgrid/raylSIF_sza"*sza_str*"_alb"*alb_str*"_psurf"*psurf_str*"hpa_nors_ABO2.dat"
            fname1   = "/home/sanghavi/RamanSIFgrid/raylSIF_sza"*sza_str*"_alb"*alb_str*"_psurf"*psurf_str*"hpa_rrs_ABO2.dat"

            #fname0_0   = "/home/sanghavi/RamanSIFgrid/rayl_sza"*sza_str*"_alb"*alb_str*"_psurf"*psurf_str*"hpa_nors_ABO2.dat"
            fname1_0   = "/home/sanghavi/RamanSIFgrid/rayl_sza"*sza_str*"_alb"*alb_str*"_psurf"*psurf_str*"hpa_rrs_ABO2.dat"

            #specNoRS = readdlm(fname0)
            #specNoRS = specNoRS[end:-1:1,:]
            specRRS_  = readdlm(fname1)
            specRRS_ = specRRS_[end:-1:1,:]
            specRRS = zeros(length(wl), size(specRRS_,2))
            for i=1:size(specRRS_,2)
                specRRS[:,i] = specRRS_[ filter(x->!(x in ind_to_elim), eachindex(specRRS_[:,i])), i]
            end
            #specRRS = specRRS[end:-1:1,:]

            #specNoRS0 = readdlm(fname0_0)
            #specNoRS0 = specNoRS0[end:-1:1,:]
            specRRS_  = readdlm(fname1_0)
            specRRS_ = specRRS_[end:-1:1,:]
            specRRS0 = zeros(length(wl), size(specRRS_,2))
            for i=1:size(specRRS_,2)
                specRRS0[:,i] = specRRS_[ filter(x->!(x in ind_to_elim), eachindex(specRRS_[:,i])), i]
            end
            #specRRS0 = specRRS0[end:-1:1,:]

            ocoI = specRRS[:,2] + specRRS[:,5]
            ocoQ = specRRS[:,3] + specRRS[:,6]
            ocoU = specRRS[:,4] + specRRS[:,7]
            ocoI_tmp = (ocoMM[1]*ocoI + 
                        ocoMM[2]*ocoQ +
                        ocoMM[3]*ocoU).* 
                        (1e7./wl.^2) #mW/m^2/s/nm

            ocoI_corr = specRRS[:,2]
            ocoQ_corr = specRRS[:,3]
            ocoU_corr = specRRS[:,4]
            ocoI_tmp_corr = (ocoMM[1]*ocoI_corr + 
                            ocoMM[2]*ocoQ_corr +
                            ocoMM[3]*ocoU_corr).* 
                            (1e7./wl.^2) #mW/m^2/s/nm
     

            ocoI0 = specRRS0[:,2] + specRRS0[:,5] 
            ocoQ0 = specRRS0[:,3] + specRRS0[:,6]
            ocoU0 = specRRS0[:,4] + specRRS0[:,7]
            ocoI0_tmp = (ocoMM[1]*ocoI0 + 
                         ocoMM[2]*ocoQ0 +
                         ocoMM[3]*ocoU0).* 
                         (1e7./wl.^2) #mW/m^2/s/nm

            ocoI0_corr = specRRS0[:,2]
            ocoQ0_corr = specRRS0[:,3]
            ocoU0_corr = specRRS0[:,4]
            ocoI0_tmp_corr = (ocoMM[1]*ocoI0_corr + 
                        ocoMM[2]*ocoQ0_corr +
                        ocoMM[3]*ocoU0_corr).* 
                        (1e7./wl.^2) #mW/m^2/s/nm             
     
            ocoF₀ = specRRS[:,8]
            ocoF₀_tmp = #(ocoF₀).*  
                        (ocoMM[1]*ocoF₀).* 
                        (1e7./wl.^2) #mW/m^2/s/nm

            interp_I  = LinearInterpolation((wl)*1e-3, (ocoI_tmp));  
            interp_I0 = LinearInterpolation((wl)*1e-3, (ocoI0_tmp)); 
            interp_I_corr  = LinearInterpolation((wl)*1e-3, (ocoI_tmp_corr));  
            interp_I0_corr = LinearInterpolation((wl)*1e-3, (ocoI0_tmp_corr));  
            interp_F₀ = LinearInterpolation((wl)*1e-3, (ocoF₀_tmp));  
            # Re-interpolate I from ν_grid to new grid/resolution
            tmp1   = interp_I(oco_wl);
            tmp1_0 = interp_I0(oco_wl);
            tmp1_corr   = interp_I_corr(oco_wl);
            tmp1_0_corr = interp_I0_corr(oco_wl);
            tmp1_F₀= interp_F₀(oco_wl);
            # Convolve input spectrum with variable kernel
            tmp2   = InstrumentOperator.conv_spectra(oco_sounding.ils[iBand], oco_wl, tmp1)
            tmp2_0 = InstrumentOperator.conv_spectra(oco_sounding.ils[iBand], oco_wl, tmp1_0)
            tmp2_corr   = InstrumentOperator.conv_spectra(oco_sounding.ils[iBand], oco_wl, tmp1_corr)
            tmp2_0_corr = InstrumentOperator.conv_spectra(oco_sounding.ils[iBand], oco_wl, tmp1_0_corr)
            tmp2_F₀= InstrumentOperator.conv_spectra(oco_sounding.ils[iBand], oco_wl, tmp1_F₀)
            
            #specNoRS0 = readdlm(fname0_0)
            #specRRS0  = readdlm(fname1_0)

            #wl = 1e7./specRRS[:,1] # converting from wn [cm^{-1}] to wl [nm]
            #ind = findall(x->x>fit_window[1] && x<fit_window[2], wl);
            
            oco_y1 = tmp2[ind1] #ieI_tmp[ind1]
            #oco_y2 = tmp2[ind2]
            oco_y1_corr = tmp2_corr[ind1]
            #oco_y2_corr = tmp2_corr[ind2]

            oco_y1_0 = tmp2_0[ind1] #ieI_tmp[ind1]
            #oco_y2_0 = tmp2_0[ind2]
            oco_y1_corr_0 = tmp2_0_corr[ind1]
            #oco_y2_corr_0 = tmp2_0_corr[ind2]

            ocoF₀_1 = tmp2_F₀[ind1]
            #ocoF₀_2 = tmp2_F₀[ind2]

            #radRRS_I_corr = (specRRS[ind1_hires,2]) .* (1e7./wl[ind1_hires].^2) #converting radiances from mW/m^2/sr/cm^{-1} to mW/m^2/sr/nm
            
            #radRRS_I = (specRRS[ind1_hires,2] + specRRS[ind1_hires,5]) .* (1e7./wl[ind1_hires].^2) #converting radiances from mW/m^2/sr/cm^{-1} to mW/m^2/sr/nm
            #radnoRS_I = specNoRS[ind1_hires,2] .* (1e7./wl[ind1_hires].^2) #converting radiances from mW/m^2/sr/cm^{-1} to mW/m^2/sr/nm

            #radRRS_I0_corr = (specRRS0[ind1_hires,2]) .* (1e7./wl[ind1_hires].^2) #converting radiances from mW/m^2/sr/cm^{-1} to mW/m^2/sr/nm
            
            #radRRS_I0 = (specRRS0[ind1_hires,2] + specRRS0[ind1_hires,5]) .* (1e7./wl[ind1_hires].^2) #converting radiances from mW/m^2/sr/cm^{-1} to mW/m^2/sr/nm
            #radnoRS_I0 = specNoRS0[ind1_hires,2] .* (1e7./wl[ind1_hires].^2) #converting radiances from mW/m^2/sr/cm^{-1} to mW/m^2/sr/nm

            #F0 = specNoRS[ind1_hires,5] .* (1e7./wl[ind1_hires].^2) #converting radiances from mW/m^2/sr/cm^{-1} to mW/m^2/sr/nm
            
            # Create synthetic noise
            # Start
            #Sig = (radRRS_I .* wl[ind1_hires]) / (1e9 * h_Pl * c_l)
            #SNR = sqrt.(100*Sig.^2 ./ (MaxMS*(c_bkg^2 * MaxMS/100 .+ c_pht^2*Sig)))
            #noise = zeros(length(radRRS_I))
            #for i=1:length(radRRS_I)
            #    noise[i] = (radRRS_I[i]/SNR[i])*randn()
            #end
            # Window 1
            ocoSig = (oco_y1 .* oco_sounding.SpectralGrid[ind1]*1e3) / (1e9 * h_Pl * c_l) #in Photons/m^2/sr/nm
            ocoSNR = sqrt.(100*ocoSig.^2 ./ (MaxMS*(c_bkg^2 * MaxMS/100 .+ c_pht^2*ocoSig)))
            oco_noise1 = zeros(length(ocoSig))
            for i=1:length(ocoSig)
                oco_noise1[i] = (oco_y1[i]/ocoSNR[i])*randn()
            end
            # Window 2
            #ocoSig = (oco_y2 .* oco_sounding.SpectralGrid[ind2]*1e3) / (1e9 * h_Pl * c_l) #in Photons/m^2/sr/nm
            #ocoSNR = sqrt.(100*ocoSig.^2 ./ (MaxMS*(c_bkg^2 * MaxMS/100 .+ c_pht^2*ocoSig)))
            #oco_noise2 = zeros(length(ocoSig))
            #for i=1:length(ocoSig)
            #    oco_noise2[i] = (oco_y2[i]/ocoSNR[i])*randn()
            #end
            # End
            
            # Fit SIF
            # Grid for Legendre Polynomials:
            #iLeg = range(-1,1,length(ind1_hires))
            # Take noRS run as solar reference (can in principle also use raw solar spectrum, shouldn't matter)
            #Io = specNoRS[ind,5] #radnoRS_I[ind]
            # Define polynomial terms
            #poly = Pl.(iLeg,(0:3)');
            # Multiply polynomial with solar spectrum
            #K_ = F0 .* poly 
            # Add a constant SIF term to the Jacobian (can add Shape later)
            #K0 = [K_ ones(length(ind1_hires))];
            #K = [K_ 1e12*wl[ind].^(-4) ones(length(ind))];
            #K = [K_ I_ie_mean ones(length(ind1_hires))];
            #K1 = [K_ I_ie_mean wl[ind1_hires] ones(length(ind1_hires))];

            iLeg1 = range(-1,1,length(ind1))
            #iLeg2 = range(-1,1,length(ind2))
            # Take noRS run as solar reference (can in principle also use raw solar spectrum, shouldn't matter)
            #Io = specNoRS[ind,5] #radnoRS_I[ind]
            # Define polynomial terms
            oco_poly1 = Pl.(iLeg1,(0:3)');
            #oco_poly2 = Pl.(iLeg2,(0:3)');
            # Multiply polynomial with solar spectrum
            ocoK_ = ocoF₀_1 .* oco_poly1 
            # Add a constant SIF term to the Jacobian (can add Shape later)
            w1_ocoK0 = [ocoK_ ones(length(ind1))];
            w1_ocoK1 = [ocoK_ oco_sounding.SpectralGrid[ind1]*1e3 ones(length(ind1))];
            
            #ocoK_ = ocoF₀_2 .* oco_poly2 
            # Add a constant SIF term to the Jacobian (can add Shape later)
            #w2_ocoK0 = [ocoK_ ones(length(ind2))];
            #w2_ocoK1 = [ocoK_ oco_sounding.SpectralGrid[ind2]*1e3 ones(length(ind2))];
            
            # Fit with Least squares:
            
            #===OCO: window 1===#
            #oco_x = ocoK \ oco_y1;
            #oco_xn = ocoK \ (oco_y1+oco_noise/coAdd_N)
            #x_noRS = K \ radnoRS_I;
            w1_oco_x0 = w1_ocoK0 \ oco_y1;
            w1_oco_xn0 = w1_ocoK0 \ (oco_y1+oco_noise1/coAdd_N)
            w1_oco_x1 = w1_ocoK1 \ oco_y1;
            w1_oco_xn1 = w1_ocoK1 \ (oco_y1+oco_noise1/coAdd_N)

            w1_oco_x0_corr = w1_ocoK0 \ oco_y1_corr;
            w1_oco_xn0_corr = w1_ocoK0 \ (oco_y1_corr+oco_noise1/coAdd_N)
            w1_oco_x1_corr = w1_ocoK1 \ oco_y1_corr; 
            w1_oco_xn1_corr = w1_ocoK1 \ (oco_y1_corr+oco_noise1/coAdd_N)

            
            #x0_noRS = K0 \ radnoRS_I;
            #oco_x_0 = ocoK \ oco_y1_0;
            #oco_xn_0 = ocoK \ (oco_y1_0+oco_noise/coAdd_N)
            #x_noRS = K \ radnoRS_I;
            w1_oco_x0_0 = w1_ocoK0 \ oco_y1_0;
            w1_oco_xn0_0 = w1_ocoK0 \ (oco_y1_0+oco_noise1/coAdd_N)
            w1_oco_x1_0 = w1_ocoK1 \ oco_y1_0;
            w1_oco_xn1_0 = w1_ocoK1 \ (oco_y1_0+oco_noise1/coAdd_N)

            w1_oco_x0_corr0 = w1_ocoK0 \ oco_y1_corr_0; 
            w1_oco_xn0_corr0 = w1_ocoK0 \ (oco_y1_corr_0+oco_noise1/coAdd_N)
            w1_oco_x1_corr0 = w1_ocoK1 \ oco_y1_corr_0; 
            w1_oco_xn1_corr0 = w1_ocoK1 \ (oco_y1_corr_0+oco_noise1/coAdd_N)

            #x0 = K_ \ radRRS_I0;
            #x0_noRS = K_ \ radnoRS_I0;
        

            #===OCO===#
            #xSIF1[isurf, iρ, iA] = mean(oco_x1[end].+oco_x1[end-1]*oco_sounding.SpectralGrid[ind1]*1e3)
            #x1_ie[isurf, iρ, iA] = oco_x1[end-2]
            #xnSIF1[isurf, iρ, iA] = mean(oco_xn1[end].+oco_xn1[end-1]*oco_sounding.SpectralGrid[ind1]*1e3)
            #xn1_ie[isurf, iρ, iA] = oco_xn1[end-2]

            #xSIF1_0[isurf, iρ, iA] = mean(oco_x1_0[end].+oco_x1_0[end-1]*oco_sounding.SpectralGrid[ind1]*1e3)
            #x1_ie_0[isurf, iρ, iA] = oco_x1_0[end-2]
            #xnSIF1_0[isurf, iρ, iA] = mean(oco_xn1_0[end].+oco_xn1_0[end-1]*oco_sounding.SpectralGrid[ind1]*1e3)
            #xn1_ie_0[isurf, iρ, iA] = oco_xn1_0[end-2]


            #xSIF[isurf, iρ, iA] = oco_x[end]
            #x_ie[isurf, iρ, iA] = oco_x[end-1]
            #xnSIF[isurf, iρ, iA] = oco_xn[end]
            #xn_ie[isurf, iρ, iA] = oco_xn[end-1]

            #xSIF_0[isurf, iρ, iA] = oco_x_0[end]
            #x_ie_0[isurf, iρ, iA] = oco_x_0[end-1]
            #xnSIF_0[isurf, iρ, iA] = oco_xn_0[end]
            #xn_ie_0[isurf, iρ, iA] = oco_xn_0[end-1]

            # Method 0
            # Uncorrected, with SIF
            w1_xSIF0[isurf, iρ, iA] = w1_oco_x0[end]/ocoMM[1]
            w1_xnSIF0[isurf, iρ, iA] = w1_oco_xn0[end]/ocoMM[1]
            # Uncorrected, no SIF
            w1_xSIF0_0[isurf, iρ, iA] = w1_oco_x0_0[end]/ocoMM[1]
            w1_xnSIF0_0[isurf, iρ, iA] = w1_oco_xn0_0[end]/ocoMM[1]

            # Corrected, with SIF
            w1_xSIF0_corr[isurf, iρ, iA] = w1_oco_x0_corr[end]/ocoMM[1]
            w1_xnSIF0_corr[isurf, iρ, iA] = w1_oco_xn0_corr[end]/ocoMM[1]
            # Corrected, no SIF
            w1_xSIF0_corr_0[isurf, iρ, iA] = w1_oco_x0_corr0[end]/ocoMM[1]
            w1_xnSIF0_corr_0[isurf, iρ, iA] = w1_oco_xn0_corr0[end]/ocoMM[1]

            #Method 1
            # Uncorrected, with SIF
            w1_xSIF1[isurf, iρ, iA] = mean(w1_oco_x1[end].+w1_oco_x1[end-1]*oco_sounding.SpectralGrid[ind1]*1e3)/ocoMM[1]
            w1_xnSIF1[isurf, iρ, iA] = mean(w1_oco_xn1[end].+w1_oco_xn1[end-1]*oco_sounding.SpectralGrid[ind1]*1e3)/ocoMM[1]
            # Uncorrected, no SIF
            w1_xSIF1_0[isurf, iρ, iA] = mean(w1_oco_x1_0[end].+w1_oco_x1_0[end-1]*oco_sounding.SpectralGrid[ind1]*1e3)/ocoMM[1]
            w1_xnSIF1_0[isurf, iρ, iA] = mean(w1_oco_xn1_0[end].+w1_oco_xn1_0[end-1]*oco_sounding.SpectralGrid[ind1]*1e3)/ocoMM[1]

            # Corrected, with SIF
            w1_xSIF1_corr[isurf, iρ, iA] = mean(w1_oco_x1_corr[end].+w1_oco_x1_corr[end-1]*oco_sounding.SpectralGrid[ind1]*1e3)/ocoMM[1]
            w1_xnSIF1_corr[isurf, iρ, iA] = mean(w1_oco_xn1_corr[end].+w1_oco_xn1_corr[end-1]*oco_sounding.SpectralGrid[ind1]*1e3)/ocoMM[1]
            # Corrected, no SIF
            w1_xSIF1_corr_0[isurf, iρ, iA] = mean(w1_oco_x1_corr0[end].+w1_oco_x1_corr0[end-1]*oco_sounding.SpectralGrid[ind1]*1e3)/ocoMM[1]
            w1_xnSIF1_corr_0[isurf, iρ, iA] = mean(w1_oco_xn1_corr0[end].+w1_oco_xn1_corr0[end-1]*oco_sounding.SpectralGrid[ind1]*1e3)/ocoMM[1]

            #xSIF_noRS[isurf, iρ, iA] = x_noRS[end]
            #coeffa[:, isurf, iρ, iA] = x
            #coeffa_noRS[:, isurf, iρ, iA] = x_noRS
            #coeff[:, isurf, iρ, iA] = x0
            #coeff_noRS[:, isurf, iρ, iA] = x0_noRS
            #xSIF0[isurf, iρ, iA] = x0[end]
            #xSIF0_noRS[isurf, iρ, iA] = x0_noRS[end]

            #===OCO: window 2===#
            #oco_x = ocoK \ oco_y1;
            #oco_xn = ocoK \ (oco_y1+oco_noise/coAdd_N)
            #x_noRS = K \ radnoRS_I;
            #oco_x0 = w2_ocoK0 \ oco_y2;
            #oco_xn0 = w2_ocoK0 \ (oco_y2+oco_noise2/coAdd_N)
            #oco_x1 = w2_ocoK1 \ oco_y2;
            #oco_xn1 = w2_ocoK1 \ (oco_y2+oco_noise2/coAdd_N)

            #oco_x0_corr = w2_ocoK0 \ oco_y2_corr;
            #oco_xn0_corr = w2_ocoK0 \ (oco_y2_corr+oco_noise2/coAdd_N)
            #oco_x1_corr = w2_ocoK1 \ oco_y2_corr; 
            #oco_xn1_corr = w2_ocoK1 \ (oco_y2_corr+oco_noise2/coAdd_N)

            
            #x0_noRS = K0 \ radnoRS_I;
            #oco_x_0 = ocoK \ oco_y1_0;
            #oco_xn_0 = ocoK \ (oco_y1_0+oco_noise/coAdd_N)
            #x_noRS = K \ radnoRS_I;
            #oco_x0_0 = w2_ocoK0 \ oco_y2_0;
            #oco_xn0_0 = w2_ocoK0 \ (oco_y2_0+oco_noise2/coAdd_N)
            #oco_x1_0 = w2_ocoK1 \ oco_y2_0;
            #oco_xn1_0 = w2_ocoK1 \ (oco_y2_0+oco_noise2/coAdd_N)

            #oco_x0_corr0 = w2_ocoK0 \ oco_y2_corr_0; 
            #oco_xn0_corr0 = w2_ocoK0 \ (oco_y2_corr_0+oco_noise2/coAdd_N)
            #oco_x1_corr0 = w2_ocoK1 \ oco_y2_corr_0; 
            #oco_xn1_corr0 = w2_ocoK1 \ (oco_y2_corr_0+oco_noise2/coAdd_N)

            

            #x0 = K_ \ radRRS_I0;
            #x0_noRS = K_ \ radnoRS_I0;
        

            #===OCO===#
            #xSIF1[isurf, iρ, iA] = mean(oco_x1[end].+oco_x1[end-1]*oco_sounding.SpectralGrid[ind1]*1e3)
            #x1_ie[isurf, iρ, iA] = oco_x1[end-2]
            #xnSIF1[isurf, iρ, iA] = mean(oco_xn1[end].+oco_xn1[end-1]*oco_sounding.SpectralGrid[ind1]*1e3)
            #xn1_ie[isurf, iρ, iA] = oco_xn1[end-2]

            #xSIF1_0[isurf, iρ, iA] = mean(oco_x1_0[end].+oco_x1_0[end-1]*oco_sounding.SpectralGrid[ind1]*1e3)
            #x1_ie_0[isurf, iρ, iA] = oco_x1_0[end-2]
            #xnSIF1_0[isurf, iρ, iA] = mean(oco_xn1_0[end].+oco_xn1_0[end-1]*oco_sounding.SpectralGrid[ind1]*1e3)
            #xn1_ie_0[isurf, iρ, iA] = oco_xn1_0[end-2]


            #xSIF[isurf, iρ, iA] = oco_x[end]
            #x_ie[isurf, iρ, iA] = oco_x[end-1]
            #xnSIF[isurf, iρ, iA] = oco_xn[end]
            #xn_ie[isurf, iρ, iA] = oco_xn[end-1]

            #xSIF_0[isurf, iρ, iA] = oco_x_0[end]
            #x_ie_0[isurf, iρ, iA] = oco_x_0[end-1]
            #xnSIF_0[isurf, iρ, iA] = oco_xn_0[end]
            #xn_ie_0[isurf, iρ, iA] = oco_xn_0[end-1]

            # Method 0
            # Uncorrected, with SIF
            #w2_xSIF0[isurf, iρ, iA] = oco_x0[end]/ocoMM[1]
            #w2_xnSIF0[isurf, iρ, iA] = oco_xn0[end]/ocoMM[1]
            # Uncorrected, no SIF
            #w2_xSIF0_0[isurf, iρ, iA] = oco_x0_0[end]/ocoMM[1]
            #w2_xnSIF0_0[isurf, iρ, iA] = oco_xn0_0[end]/ocoMM[1]

            # Corrected, with SIF
            #w2_xSIF0_corr[isurf, iρ, iA] = oco_x0_corr[end]/ocoMM[1]
            #w2_xnSIF0_corr[isurf, iρ, iA] = oco_xn0_corr[end]/ocoMM[1]
            # Corrected, no SIF
            #w2_xSIF0_corr_0[isurf, iρ, iA] = oco_x0_corr0[end]/ocoMM[1]
            #w2_xnSIF0_corr_0[isurf, iρ, iA] = oco_xn0_corr0[end]/ocoMM[1]

            #Method 1
            # Uncorrected, with SIF
            #w2_xSIF1[isurf, iρ, iA] = mean(oco_x1[end].+oco_x1[end-1]*oco_sounding.SpectralGrid[ind2]*1e3)/ocoMM[1]
            #w2_xnSIF1[isurf, iρ, iA] = mean(oco_xn1[end].+oco_xn1[end-1]*oco_sounding.SpectralGrid[ind2]*1e3)/ocoMM[1]
            # Uncorrected, no SIF
            #w2_xSIF1_0[isurf, iρ, iA] = mean(oco_x1_0[end].+oco_x1_0[end-1]*oco_sounding.SpectralGrid[ind2]*1e3)/ocoMM[1]
            #w2_xnSIF1_0[isurf, iρ, iA] = mean(oco_xn1_0[end].+oco_xn1_0[end-1]*oco_sounding.SpectralGrid[ind2]*1e3)/ocoMM[1]

            # Corrected, with SIF
            #w2_xSIF1_corr[isurf, iρ, iA] = mean(oco_x1_corr[end].+oco_x1_corr[end-1]*oco_sounding.SpectralGrid[ind2]*1e3)/ocoMM[1]
            #w2_xnSIF1_corr[isurf, iρ, iA] = mean(oco_xn1_corr[end].+oco_xn1_corr[end-1]*oco_sounding.SpectralGrid[ind2]*1e3)/ocoMM[1]
            # Corrected, no SIF
            #w2_xSIF1_corr_0[isurf, iρ, iA] = mean(oco_x1_corr0[end].+oco_x1_corr0[end-1]*oco_sounding.SpectralGrid[ind2]*1e3)/ocoMM[1]
            #w2_xnSIF1_corr_0[isurf, iρ, iA] = mean(oco_xn1_corr0[end].+oco_xn1_corr0[end-1]*oco_sounding.SpectralGrid[ind2]*1e3)/ocoMM[1]

        end
    end
end


for isurf = 1:3 # 1:1 #2:2 # 
    for iρ = 1:15 #3 #1:15
        for iA = 1:14
            sza_str  = string(sza[iA])
            alb_str  = ρ_str[iρ]
            psurf_str= string(psurf[isurf])

            #fname0   = "/home/sanghavi/RamanSIFgrid/raylSIF_sza"*sza_str*"_alb"*alb_str*"_psurf"*psurf_str*"hpa_nors_ABO2.dat"
            fname1   = "/home/sanghavi/RamanSIFgrid/raylSIF_sza"*sza_str*"_alb"*alb_str*"_psurf"*psurf_str*"hpa_rrs_ABO2.dat"

            #fname0_0   = "/home/sanghavi/RamanSIFgrid/rayl_sza"*sza_str*"_alb"*alb_str*"_psurf"*psurf_str*"hpa_nors_ABO2.dat"
            fname1_0   = "/home/sanghavi/RamanSIFgrid/rayl_sza"*sza_str*"_alb"*alb_str*"_psurf"*psurf_str*"hpa_rrs_ABO2.dat"

            #specNoRS = readdlm(fname0)
            #specNoRS = specNoRS[end:-1:1,:]
            specRRS_  = readdlm(fname1)
            specRRS_ = specRRS_[end:-1:1,:]
            specRRS = zeros(length(wl), size(specRRS_,2))
            for i=1:size(specRRS_,2)
                specRRS[:,i] = specRRS_[ filter(x->!(x in ind_to_elim), eachindex(specRRS_[:,i])), i]
            end
            #specRRS = specRRS[end:-1:1,:]

            #specNoRS0 = readdlm(fname0_0)
            #specNoRS0 = specNoRS0[end:-1:1,:]
            specRRS_  = readdlm(fname1_0)
            specRRS_ = specRRS_[end:-1:1,:]
            specRRS0 = zeros(length(wl), size(specRRS_,2))
            for i=1:size(specRRS_,2)
                specRRS0[:,i] = specRRS_[ filter(x->!(x in ind_to_elim), eachindex(specRRS_[:,i])), i]
            end
            #specRRS0 = specRRS0[end:-1:1,:]

            ocoI = specRRS[:,2] + specRRS[:,5]
            ocoQ = specRRS[:,3] + specRRS[:,6]
            ocoU = specRRS[:,4] + specRRS[:,7]
            ocoI_tmp = (ocoMM[1]*ocoI + 
                        ocoMM[2]*ocoQ +
                        ocoMM[3]*ocoU).* 
                        (1e7./wl.^2) #mW/m^2/s/nm

            ocoI_corr = specRRS[:,2]
            ocoQ_corr = specRRS[:,3]
            ocoU_corr = specRRS[:,4]
            ocoI_tmp_corr = (ocoMM[1]*ocoI_corr + 
                            ocoMM[2]*ocoQ_corr +
                            ocoMM[3]*ocoU_corr).* 
                            (1e7./wl.^2) #mW/m^2/s/nm
     

            ocoI0 = specRRS0[:,2] + specRRS0[:,5] 
            ocoQ0 = specRRS0[:,3] + specRRS0[:,6]
            ocoU0 = specRRS0[:,4] + specRRS0[:,7]
            ocoI0_tmp = (ocoMM[1]*ocoI0 + 
                         ocoMM[2]*ocoQ0 +
                         ocoMM[3]*ocoU0).* 
                         (1e7./wl.^2) #mW/m^2/s/nm

            ocoI0_corr = specRRS0[:,2]
            ocoQ0_corr = specRRS0[:,3]
            ocoU0_corr = specRRS0[:,4]
            ocoI0_tmp_corr = (ocoMM[1]*ocoI0_corr + 
                        ocoMM[2]*ocoQ0_corr +
                        ocoMM[3]*ocoU0_corr).* 
                        (1e7./wl.^2) #mW/m^2/s/nm             
     
            ocoF₀ = specRRS0[:,2]
            ocoF₀_tmp = #(ocoF₀).*  
                        (ocoMM[1]*ocoF₀).* 
                        (1e7./wl.^2) #mW/m^2/s/nm

            interp_I  = LinearInterpolation((wl)*1e-3, (ocoI_tmp));  
            interp_I0 = LinearInterpolation((wl)*1e-3, (ocoI0_tmp)); 
            interp_I_corr  = LinearInterpolation((wl)*1e-3, (ocoI_tmp_corr));  
            interp_I0_corr = LinearInterpolation((wl)*1e-3, (ocoI0_tmp_corr));  
            interp_F₀ = LinearInterpolation((wl)*1e-3, (ocoF₀_tmp));  
            # Re-interpolate I from ν_grid to new grid/resolution
            tmp1   = interp_I(oco_wl);
            tmp1_0 = interp_I0(oco_wl);
            tmp1_corr   = interp_I_corr(oco_wl);
            tmp1_0_corr = interp_I0_corr(oco_wl);
            tmp1_F₀= interp_F₀(oco_wl);
            # Convolve input spectrum with variable kernel
            tmp2   = InstrumentOperator.conv_spectra(oco_sounding.ils[iBand], oco_wl, tmp1)
            tmp2_0 = InstrumentOperator.conv_spectra(oco_sounding.ils[iBand], oco_wl, tmp1_0)
            tmp2_corr   = InstrumentOperator.conv_spectra(oco_sounding.ils[iBand], oco_wl, tmp1_corr)
            tmp2_0_corr = InstrumentOperator.conv_spectra(oco_sounding.ils[iBand], oco_wl, tmp1_0_corr)
            tmp2_F₀= InstrumentOperator.conv_spectra(oco_sounding.ils[iBand], oco_wl, tmp1_F₀)
            
            #specNoRS0 = readdlm(fname0_0)
            #specRRS0  = readdlm(fname1_0)

            #wl = 1e7./specRRS[:,1] # converting from wn [cm^{-1}] to wl [nm]
            #ind = findall(x->x>fit_window[1] && x<fit_window[2], wl);
            
            #oco_y1 = tmp2[ind1] #ieI_tmp[ind1]
            oco_y2 = tmp2[ind2]
            #oco_y1_corr = tmp2_corr[ind1]
            oco_y2_corr = tmp2_corr[ind2]

            #oco_y1_0 = tmp2_0[ind1] #ieI_tmp[ind1]
            oco_y2_0 = tmp2_0[ind2]
            #oco_y1_corr_0 = tmp2_0_corr[ind1]
            oco_y2_corr_0 = tmp2_0_corr[ind2]

            #ocoF₀_1 = tmp2_F₀[ind1]
            ocoF₀_2 = tmp2_F₀[ind2]

            #radRRS_I_corr = (specRRS[ind1_hires,2]) .* (1e7./wl[ind1_hires].^2) #converting radiances from mW/m^2/sr/cm^{-1} to mW/m^2/sr/nm
            
            #radRRS_I = (specRRS[ind1_hires,2] + specRRS[ind1_hires,5]) .* (1e7./wl[ind1_hires].^2) #converting radiances from mW/m^2/sr/cm^{-1} to mW/m^2/sr/nm
            #radnoRS_I = specNoRS[ind1_hires,2] .* (1e7./wl[ind1_hires].^2) #converting radiances from mW/m^2/sr/cm^{-1} to mW/m^2/sr/nm

            #radRRS_I0_corr = (specRRS0[ind1_hires,2]) .* (1e7./wl[ind1_hires].^2) #converting radiances from mW/m^2/sr/cm^{-1} to mW/m^2/sr/nm
            
            #radRRS_I0 = (specRRS0[ind1_hires,2] + specRRS0[ind1_hires,5]) .* (1e7./wl[ind1_hires].^2) #converting radiances from mW/m^2/sr/cm^{-1} to mW/m^2/sr/nm
            #radnoRS_I0 = specNoRS0[ind1_hires,2] .* (1e7./wl[ind1_hires].^2) #converting radiances from mW/m^2/sr/cm^{-1} to mW/m^2/sr/nm

            #F0 = specNoRS[ind1_hires,5] .* (1e7./wl[ind1_hires].^2) #converting radiances from mW/m^2/sr/cm^{-1} to mW/m^2/sr/nm
            
            # Create synthetic noise
            # Start
            #Sig = (radRRS_I .* wl[ind1_hires]) / (1e9 * h_Pl * c_l)
            #SNR = sqrt.(100*Sig.^2 ./ (MaxMS*(c_bkg^2 * MaxMS/100 .+ c_pht^2*Sig)))
            #noise = zeros(length(radRRS_I))
            #for i=1:length(radRRS_I)
            #    noise[i] = (radRRS_I[i]/SNR[i])*randn()
            #end
            # Window 1
            #ocoSig = (oco_y1 .* oco_sounding.SpectralGrid[ind1]*1e3) / (1e9 * h_Pl * c_l) #in Photons/m^2/sr/nm
            #ocoSNR = sqrt.(100*ocoSig.^2 ./ (MaxMS*(c_bkg^2 * MaxMS/100 .+ c_pht^2*ocoSig)))
            #oco_noise1 = zeros(length(ocoSig))
            #for i=1:length(ocoSig)
            #    oco_noise1[i] = (oco_y1[i]/ocoSNR[i])*randn()
            #end
            # Window 2
            ocoSig = (oco_y2 .* oco_sounding.SpectralGrid[ind2]*1e3) / (1e9 * h_Pl * c_l) #in Photons/m^2/sr/nm
            ocoSNR = sqrt.(100*ocoSig.^2 ./ (MaxMS*(c_bkg^2 * MaxMS/100 .+ c_pht^2*ocoSig)))
            oco_noise2 = zeros(length(ocoSig))
            for i=1:length(ocoSig)
                oco_noise2[i] = (oco_y2[i]/ocoSNR[i])*randn()
            end
            # End
            
            # Fit SIF
            # Grid for Legendre Polynomials:
            #iLeg = range(-1,1,length(ind1_hires))
            # Take noRS run as solar reference (can in principle also use raw solar spectrum, shouldn't matter)
            #Io = specNoRS[ind,5] #radnoRS_I[ind]
            # Define polynomial terms
            #poly = Pl.(iLeg,(0:3)');
            # Multiply polynomial with solar spectrum
            #K_ = F0 .* poly 
            # Add a constant SIF term to the Jacobian (can add Shape later)
            #K0 = [K_ ones(length(ind1_hires))];
            #K = [K_ 1e12*wl[ind].^(-4) ones(length(ind))];
            #K = [K_ I_ie_mean ones(length(ind1_hires))];
            #K1 = [K_ I_ie_mean wl[ind1_hires] ones(length(ind1_hires))];

            #iLeg1 = range(-1,1,length(ind1))
            iLeg2 = range(-1,1,length(ind2))
            # Take noRS run as solar reference (can in principle also use raw solar spectrum, shouldn't matter)
            #Io = specNoRS[ind,5] #radnoRS_I[ind]
            # Define polynomial terms
            #oco_poly1 = Pl.(iLeg1,(0:3)');
            oco_poly2 = Pl.(iLeg2,(0:3)');
            # Multiply polynomial with solar spectrum
            #ocoK_ = ocoF₀_1 .* oco_poly1 
            # Add a constant SIF term to the Jacobian (can add Shape later)
            #w1_ocoK0 = [ocoK_ ones(length(ind1))];
            #w1_ocoK1 = [ocoK_ oco_sounding.SpectralGrid[ind1]*1e3 ones(length(ind1))];
            
            ocoK_ = ocoF₀_2 .* oco_poly2 
            # Add a constant SIF term to the Jacobian (can add Shape later)
            w2_ocoK0 = [ocoK_ ones(length(ind2))];
            w2_ocoK1 = [ocoK_ oco_sounding.SpectralGrid[ind2]*1e3 ones(length(ind2))];
            
            # Fit with Least squares:
            
            #===OCO: window 1===#
            #oco_x = ocoK \ oco_y1;
            #oco_xn = ocoK \ (oco_y1+oco_noise/coAdd_N)
            #x_noRS = K \ radnoRS_I;
            #w1_oco_x0 = w1_ocoK0 \ oco_y1;
            #w1_oco_xn0 = w1_ocoK0 \ (oco_y1+oco_noise1/coAdd_N)
            #w1_oco_x1 = w1_ocoK1 \ oco_y1;
            #w1_oco_xn1 = w1_ocoK1 \ (oco_y1+oco_noise1/coAdd_N)

            #w1_oco_x0_corr = w1_ocoK0 \ oco_y1_corr;
            #w1_oco_xn0_corr = w1_ocoK0 \ (oco_y1_corr+oco_noise1/coAdd_N)
            #w1_oco_x1_corr = w1_ocoK1 \ oco_y1_corr; 
            #w1_oco_xn1_corr = w1_ocoK1 \ (oco_y1_corr+oco_noise1/coAdd_N)

            
            #x0_noRS = K0 \ radnoRS_I;
            #oco_x_0 = ocoK \ oco_y1_0;
            #oco_xn_0 = ocoK \ (oco_y1_0+oco_noise/coAdd_N)
            #x_noRS = K \ radnoRS_I;
            #w1_oco_x0_0 = w1_ocoK0 \ oco_y1_0;
            #w1_oco_xn0_0 = w1_ocoK0 \ (oco_y1_0+oco_noise1/coAdd_N)
            #w1_oco_x1_0 = w1_ocoK1 \ oco_y1_0;
            #w1_oco_xn1_0 = w1_ocoK1 \ (oco_y1_0+oco_noise1/coAdd_N)

            #w1_oco_x0_corr0 = w1_ocoK0 \ oco_y1_corr_0; 
            #w1_oco_xn0_corr0 = w1_ocoK0 \ (oco_y1_corr_0+oco_noise1/coAdd_N)
            #w1_oco_x1_corr0 = w1_ocoK1 \ oco_y1_corr_0; 
            #w1_oco_xn1_corr0 = w1_ocoK1 \ (oco_y1_corr_0+oco_noise1/coAdd_N)

            #x0 = K_ \ radRRS_I0;
            #x0_noRS = K_ \ radnoRS_I0;
        

            #===OCO===#
            #xSIF1[isurf, iρ, iA] = mean(oco_x1[end].+oco_x1[end-1]*oco_sounding.SpectralGrid[ind1]*1e3)
            #x1_ie[isurf, iρ, iA] = oco_x1[end-2]
            #xnSIF1[isurf, iρ, iA] = mean(oco_xn1[end].+oco_xn1[end-1]*oco_sounding.SpectralGrid[ind1]*1e3)
            #xn1_ie[isurf, iρ, iA] = oco_xn1[end-2]

            #xSIF1_0[isurf, iρ, iA] = mean(oco_x1_0[end].+oco_x1_0[end-1]*oco_sounding.SpectralGrid[ind1]*1e3)
            #x1_ie_0[isurf, iρ, iA] = oco_x1_0[end-2]
            #xnSIF1_0[isurf, iρ, iA] = mean(oco_xn1_0[end].+oco_xn1_0[end-1]*oco_sounding.SpectralGrid[ind1]*1e3)
            #xn1_ie_0[isurf, iρ, iA] = oco_xn1_0[end-2]


            #xSIF[isurf, iρ, iA] = oco_x[end]
            #x_ie[isurf, iρ, iA] = oco_x[end-1]
            #xnSIF[isurf, iρ, iA] = oco_xn[end]
            #xn_ie[isurf, iρ, iA] = oco_xn[end-1]

            #xSIF_0[isurf, iρ, iA] = oco_x_0[end]
            #x_ie_0[isurf, iρ, iA] = oco_x_0[end-1]
            #xnSIF_0[isurf, iρ, iA] = oco_xn_0[end]
            #xn_ie_0[isurf, iρ, iA] = oco_xn_0[end-1]

            # Method 0
            # Uncorrected, with SIF
            #w1_xSIF0[isurf, iρ, iA] = w1_oco_x0[end]/ocoMM[1]
            #w1_xnSIF0[isurf, iρ, iA] = w1_oco_xn0[end]/ocoMM[1]
            # Uncorrected, no SIF
            #w1_xSIF0_0[isurf, iρ, iA] = w1_oco_x0_0[end]/ocoMM[1]
            #w1_xnSIF0_0[isurf, iρ, iA] = w1_oco_xn0_0[end]/ocoMM[1]

            # Corrected, with SIF
            #w1_xSIF0_corr[isurf, iρ, iA] = w1_oco_x0_corr[end]/ocoMM[1]
            #w1_xnSIF0_corr[isurf, iρ, iA] = w1_oco_xn0_corr[end]/ocoMM[1]
            # Corrected, no SIF
            #w1_xSIF0_corr_0[isurf, iρ, iA] = w1_oco_x0_corr0[end]/ocoMM[1]
            #w1_xnSIF0_corr_0[isurf, iρ, iA] = w1_oco_xn0_corr0[end]/ocoMM[1]

            #Method 1
            # Uncorrected, with SIF
            #w1_xSIF1[isurf, iρ, iA] = mean(w1_oco_x1[end].+w1_oco_x1[end-1]*oco_sounding.SpectralGrid[ind1]*1e3)/ocoMM[1]
            #w1_xnSIF1[isurf, iρ, iA] = mean(w1_oco_xn1[end].+w1_oco_xn1[end-1]*oco_sounding.SpectralGrid[ind1]*1e3)/ocoMM[1]
            # Uncorrected, no SIF
            #w1_xSIF1_0[isurf, iρ, iA] = mean(w1_oco_x1_0[end].+w1_oco_x1_0[end-1]*oco_sounding.SpectralGrid[ind1]*1e3)/ocoMM[1]
            #w1_xnSIF1_0[isurf, iρ, iA] = mean(w1_oco_xn1_0[end].+w1_oco_xn1_0[end-1]*oco_sounding.SpectralGrid[ind1]*1e3)/ocoMM[1]

            # Corrected, with SIF
            #w1_xSIF1_corr[isurf, iρ, iA] = mean(w1_oco_x1_corr[end].+w1_oco_x1_corr[end-1]*oco_sounding.SpectralGrid[ind1]*1e3)/ocoMM[1]
            #w1_xnSIF1_corr[isurf, iρ, iA] = mean(w1_oco_xn1_corr[end].+w1_oco_xn1_corr[end-1]*oco_sounding.SpectralGrid[ind1]*1e3)/ocoMM[1]
            # Corrected, no SIF
            #w1_xSIF1_corr_0[isurf, iρ, iA] = mean(w1_oco_x1_corr0[end].+w1_oco_x1_corr0[end-1]*oco_sounding.SpectralGrid[ind1]*1e3)/ocoMM[1]
            #w1_xnSIF1_corr_0[isurf, iρ, iA] = mean(w1_oco_xn1_corr0[end].+w1_oco_xn1_corr0[end-1]*oco_sounding.SpectralGrid[ind1]*1e3)/ocoMM[1]

            #xSIF_noRS[isurf, iρ, iA] = x_noRS[end]
            #coeffa[:, isurf, iρ, iA] = x
            #coeffa_noRS[:, isurf, iρ, iA] = x_noRS
            #coeff[:, isurf, iρ, iA] = x0
            #coeff_noRS[:, isurf, iρ, iA] = x0_noRS
            #xSIF0[isurf, iρ, iA] = x0[end]
            #xSIF0_noRS[isurf, iρ, iA] = x0_noRS[end]

            #===OCO: window 2===#
            #oco_x = ocoK \ oco_y1;
            #oco_xn = ocoK \ (oco_y1+oco_noise/coAdd_N)
            #x_noRS = K \ radnoRS_I;
            oco_x0 = w2_ocoK0 \ oco_y2;
            oco_xn0 = w2_ocoK0 \ (oco_y2+oco_noise2/coAdd_N)
            oco_x1 = w2_ocoK1 \ oco_y2;
            oco_xn1 = w2_ocoK1 \ (oco_y2+oco_noise2/coAdd_N)

            oco_x0_corr = w2_ocoK0 \ oco_y2_corr;
            oco_xn0_corr = w2_ocoK0 \ (oco_y2_corr+oco_noise2/coAdd_N)
            oco_x1_corr = w2_ocoK1 \ oco_y2_corr; 
            oco_xn1_corr = w2_ocoK1 \ (oco_y2_corr+oco_noise2/coAdd_N)

            
            #x0_noRS = K0 \ radnoRS_I;
            #oco_x_0 = ocoK \ oco_y1_0;
            #oco_xn_0 = ocoK \ (oco_y1_0+oco_noise/coAdd_N)
            #x_noRS = K \ radnoRS_I;
            oco_x0_0 = w2_ocoK0 \ oco_y2_0;
            oco_xn0_0 = w2_ocoK0 \ (oco_y2_0+oco_noise2/coAdd_N)
            oco_x1_0 = w2_ocoK1 \ oco_y2_0;
            oco_xn1_0 = w2_ocoK1 \ (oco_y2_0+oco_noise2/coAdd_N)

            oco_x0_corr0 = w2_ocoK0 \ oco_y2_corr_0; 
            oco_xn0_corr0 = w2_ocoK0 \ (oco_y2_corr_0+oco_noise2/coAdd_N)
            oco_x1_corr0 = w2_ocoK1 \ oco_y2_corr_0; 
            oco_xn1_corr0 = w2_ocoK1 \ (oco_y2_corr_0+oco_noise2/coAdd_N)

            

            #x0 = K_ \ radRRS_I0;
            #x0_noRS = K_ \ radnoRS_I0;
        

            #===OCO===#
            #xSIF1[isurf, iρ, iA] = mean(oco_x1[end].+oco_x1[end-1]*oco_sounding.SpectralGrid[ind1]*1e3)
            #x1_ie[isurf, iρ, iA] = oco_x1[end-2]
            #xnSIF1[isurf, iρ, iA] = mean(oco_xn1[end].+oco_xn1[end-1]*oco_sounding.SpectralGrid[ind1]*1e3)
            #xn1_ie[isurf, iρ, iA] = oco_xn1[end-2]

            #xSIF1_0[isurf, iρ, iA] = mean(oco_x1_0[end].+oco_x1_0[end-1]*oco_sounding.SpectralGrid[ind1]*1e3)
            #x1_ie_0[isurf, iρ, iA] = oco_x1_0[end-2]
            #xnSIF1_0[isurf, iρ, iA] = mean(oco_xn1_0[end].+oco_xn1_0[end-1]*oco_sounding.SpectralGrid[ind1]*1e3)
            #xn1_ie_0[isurf, iρ, iA] = oco_xn1_0[end-2]


            #xSIF[isurf, iρ, iA] = oco_x[end]
            #x_ie[isurf, iρ, iA] = oco_x[end-1]
            #xnSIF[isurf, iρ, iA] = oco_xn[end]
            #xn_ie[isurf, iρ, iA] = oco_xn[end-1]

            #xSIF_0[isurf, iρ, iA] = oco_x_0[end]
            #x_ie_0[isurf, iρ, iA] = oco_x_0[end-1]
            #xnSIF_0[isurf, iρ, iA] = oco_xn_0[end]
            #xn_ie_0[isurf, iρ, iA] = oco_xn_0[end-1]

            # Method 0
            # Uncorrected, with SIF
            w2_xSIF0[isurf, iρ, iA] = oco_x0[end]/ocoMM[1]
            w2_xnSIF0[isurf, iρ, iA] = oco_xn0[end]/ocoMM[1]
            # Uncorrected, no SIF
            w2_xSIF0_0[isurf, iρ, iA] = oco_x0_0[end]/ocoMM[1]
            w2_xnSIF0_0[isurf, iρ, iA] = oco_xn0_0[end]/ocoMM[1]

            # Corrected, with SIF
            w2_xSIF0_corr[isurf, iρ, iA] = oco_x0_corr[end]/ocoMM[1]
            w2_xnSIF0_corr[isurf, iρ, iA] = oco_xn0_corr[end]/ocoMM[1]
            # Corrected, no SIF
            w2_xSIF0_corr_0[isurf, iρ, iA] = oco_x0_corr0[end]/ocoMM[1]
            w2_xnSIF0_corr_0[isurf, iρ, iA] = oco_xn0_corr0[end]/ocoMM[1]

            #Method 1
            # Uncorrected, with SIF
            w2_xSIF1[isurf, iρ, iA] = mean(oco_x1[end].+oco_x1[end-1]*oco_sounding.SpectralGrid[ind2]*1e3)/ocoMM[1]
            w2_xnSIF1[isurf, iρ, iA] = mean(oco_xn1[end].+oco_xn1[end-1]*oco_sounding.SpectralGrid[ind2]*1e3)/ocoMM[1]
            # Uncorrected, no SIF
            w2_xSIF1_0[isurf, iρ, iA] = mean(oco_x1_0[end].+oco_x1_0[end-1]*oco_sounding.SpectralGrid[ind2]*1e3)/ocoMM[1]
            w2_xnSIF1_0[isurf, iρ, iA] = mean(oco_xn1_0[end].+oco_xn1_0[end-1]*oco_sounding.SpectralGrid[ind2]*1e3)/ocoMM[1]

            # Corrected, with SIF
            w2_xSIF1_corr[isurf, iρ, iA] = mean(oco_x1_corr[end].+oco_x1_corr[end-1]*oco_sounding.SpectralGrid[ind2]*1e3)/ocoMM[1]
            w2_xnSIF1_corr[isurf, iρ, iA] = mean(oco_xn1_corr[end].+oco_xn1_corr[end-1]*oco_sounding.SpectralGrid[ind2]*1e3)/ocoMM[1]
            # Corrected, no SIF
            w2_xSIF1_corr_0[isurf, iρ, iA] = mean(oco_x1_corr0[end].+oco_x1_corr0[end-1]*oco_sounding.SpectralGrid[ind2]*1e3)/ocoMM[1]
            w2_xnSIF1_corr_0[isurf, iρ, iA] = mean(oco_xn1_corr0[end].+oco_xn1_corr0[end-1]*oco_sounding.SpectralGrid[ind2]*1e3)/ocoMM[1]

        end
    end
end

l = @layout [a1 a2; b1 b2]
p1 = heatmap(ρ, sza, w1_xSIF0[3,:,:]', xlabel="", ylabel="SZA [ᵒ]")
# title=L"$\mathrm{Uncorr.}\,\mathrm{SIF\,\,retr.} (\hat{\mathrm{SIF}}=0.34\,\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m})$")
p2 = heatmap(ρ, sza, w2_xSIF0[3,:,:]', xlabel="", ylabel="")
# title=L"$\mathrm{Uncorr.}\,\mathrm{SIF\,\,retr.} (\hat{\mathrm{SIF}}=0.20\,\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m})$")
p3 = heatmap(ρ, sza, w1_xSIF0_corr[3,:,:]', xlabel="Albedo", ylabel="SZA [ᵒ]") 

p4 = heatmap(ρ, sza, w2_xSIF0_corr[3,:,:]', xlabel="Albedo", ylabel="") 
        
title1a = L"$\mathrm{Uncorr.}\,\mathrm{SIF\,\,retr.}$" 
title1b = L"$(\hat{\mathrm{SIF}}=0.34\,\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m})$"
title2a = L"$\mathrm{Uncorr.}\,\mathrm{SIF\,\,retr.}$"
title2b = L"$(\hat{\mathrm{SIF}}=0.20\,\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m})$"
title3a = L"$\mathrm{Corr.}\,\mathrm{SIF\,\,retr.}$"
title3b = L"$(\hat{\mathrm{SIF}}=0.34\,\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m})$"
title4a = L"$\mathrm{Corr.}\,\mathrm{SIF\,\,retr.}$"
title4b = L"$(\hat{\mathrm{SIF}}=0.20\,\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m})$"

plot(p1, p2, p3, p4, layout = l, legend = false, title = [title1a*"\n"*title1b title2a*"\n"*title2b title3a*"\n"*title3b title4a*"\n"*title4b], titlefont = font(8))
savefig("/home/sanghavi/RamanSIFgrid/plots/SIF0_retr.png")

l = @layout [a1 a2; b1 b2]
p1 = heatmap(ρ, sza, w1_xnSIF0[3,:,:]', xlabel="", ylabel="SZA [ᵒ]")
# title=L"$\mathrm{Uncorr.}\,\mathrm{SIF\,\,retr.} (\hat{\mathrm{SIF}}=0.34\,\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m})$")
p2 = heatmap(ρ, sza, w2_xnSIF0[3,:,:]', xlabel="", ylabel="")
# title=L"$\mathrm{Uncorr.}\,\mathrm{SIF\,\,retr.} (\hat{\mathrm{SIF}}=0.20\,\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m})$")
p3 = heatmap(ρ, sza, w1_xnSIF0_corr[3,:,:]', xlabel="Albedo", ylabel="SZA [ᵒ]") 

p4 = heatmap(ρ, sza, w2_xnSIF0_corr[3,:,:]', xlabel="Albedo", ylabel="") 
        
title1a = L"$\mathrm{Uncorr.}\,\mathrm{SIF\,\,retr.}$" 
title1b = L"$(\hat{\mathrm{SIF}}=0.34\,\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m})$"
title2a = L"$\mathrm{Uncorr.}\,\mathrm{SIF\,\,retr.}$"
title2b = L"$(\hat{\mathrm{SIF}}=0.20\,\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m})$"
title3a = L"$\mathrm{Corr.}\,\mathrm{SIF\,\,retr.}$"
title3b = L"$(\hat{\mathrm{SIF}}=0.34\,\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m})$"
title4a = L"$\mathrm{Corr.}\,\mathrm{SIF\,\,retr.}$"
title4b = L"$(\hat{\mathrm{SIF}}=0.20\,\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m})$"

plot(p1, p2, p3, p4, layout = l, legend = false, title = [title1a*"\n"*title1b title2a*"\n"*title2b title3a*"\n"*title3b title4a*"\n"*title4b], titlefont = font(8))
savefig("/home/sanghavi/RamanSIFgrid/plots/SIF0_retr_noisy36.png")

l = @layout [a1 a2; b1 b2]
p1 = heatmap(ρ, sza, w1_xSIF1[3,:,:]', xlabel="", ylabel="SZA [ᵒ]")
# title=L"$\mathrm{Uncorr.}\,\mathrm{SIF\,\,retr.} (\hat{\mathrm{SIF}}=0.34\,\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m})$")
p2 = heatmap(ρ, sza, w2_xSIF1[3,:,:]', xlabel="", ylabel="")
# title=L"$\mathrm{Uncorr.}\,\mathrm{SIF\,\,retr.} (\hat{\mathrm{SIF}}=0.20\,\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m})$")
p3 = heatmap(ρ, sza, w1_xSIF1_corr[3,:,:]', xlabel="Albedo", ylabel="SZA [ᵒ]") 

p4 = heatmap(ρ, sza, w2_xSIF1_corr[3,:,:]', xlabel="Albedo", ylabel="") 
        
title1a = L"$\mathrm{Uncorr.}\,\mathrm{SIF\,\,retr.}$" 
title1b = L"$(\hat{\mathrm{SIF}}=0.34\,\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m})$"
title2a = L"$\mathrm{Uncorr.}\,\mathrm{SIF\,\,retr.}$"
title2b = L"$(\hat{\mathrm{SIF}}=0.20\,\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m})$"
title3a = L"$\mathrm{Corr.}\,\mathrm{SIF\,\,retr.}$"
title3b = L"$(\hat{\mathrm{SIF}}=0.34\,\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m})$"
title4a = L"$\mathrm{Corr.}\,\mathrm{SIF\,\,retr.}$"
title4b = L"$(\hat{\mathrm{SIF}}=0.20\,\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m})$"

plot(p1, p2, p3, p4, layout = l, legend = false, title = [title1a*"\n"*title1b title2a*"\n"*title2b title3a*"\n"*title3b title4a*"\n"*title4b], titlefont = font(8))
savefig("/home/sanghavi/RamanSIFgrid/plots/SIF1_retr.png")

l = @layout [a1 a2; b1 b2]
p1 = heatmap(ρ, sza, w1_xnSIF1[3,:,:]', xlabel="", ylabel="SZA [ᵒ]")
# title=L"$\mathrm{Uncorr.}\,\mathrm{SIF\,\,retr.} (\hat{\mathrm{SIF}}=0.34\,\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m})$")
p2 = heatmap(ρ, sza, w2_xnSIF1[3,:,:]', xlabel="", ylabel="")
# title=L"$\mathrm{Uncorr.}\,\mathrm{SIF\,\,retr.} (\hat{\mathrm{SIF}}=0.20\,\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m})$")
p3 = heatmap(ρ, sza, w1_xnSIF1_corr[3,:,:]', xlabel="Albedo", ylabel="SZA [ᵒ]") 

p4 = heatmap(ρ, sza, w2_xnSIF1_corr[3,:,:]', xlabel="Albedo", ylabel="") 
        
title1a = L"$\mathrm{Uncorr.}\,\mathrm{SIF\,\,retr.}$" 
title1b = L"$(\hat{\mathrm{SIF}}=0.34\,\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m})$"
title2a = L"$\mathrm{Uncorr.}\,\mathrm{SIF\,\,retr.}$"
title2b = L"$(\hat{\mathrm{SIF}}=0.20\,\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m})$"
title3a = L"$\mathrm{Corr.}\,\mathrm{SIF\,\,retr.}$"
title3b = L"$(\hat{\mathrm{SIF}}=0.34\,\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m})$"
title4a = L"$\mathrm{Corr.}\,\mathrm{SIF\,\,retr.}$"
title4b = L"$(\hat{\mathrm{SIF}}=0.20\,\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m})$"

plot(p1, p2, p3, p4, layout = l, legend = false, title = [title1a*"\n"*title1b title2a*"\n"*title2b title3a*"\n"*title3b title4a*"\n"*title4b], titlefont = font(8))
savefig("/home/sanghavi/RamanSIFgrid/plots/SIF1_retr_noisy36.png")

l = @layout [a1 a2; b1 b2]
p1 = heatmap(ρ, sza, w1_xSIF0_0[3,:,:]', xlabel="", ylabel="SZA [ᵒ]")
# title=L"$\mathrm{Uncorr.}\,\mathrm{SIF\,\,retr.} (\hat{\mathrm{SIF}}=0.34\,\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m})$")
p2 = heatmap(ρ, sza, w2_xSIF0_0[3,:,:]', xlabel="", ylabel="")
# title=L"$\mathrm{Uncorr.}\,\mathrm{SIF\,\,retr.} (\hat{\mathrm{SIF}}=0.20\,\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m})$")
p3 = heatmap(ρ, sza, w1_xSIF0_corr_0[3,:,:]', xlabel="Albedo", ylabel="SZA [ᵒ]") 

p4 = heatmap(ρ, sza, w2_xSIF0_corr_0[3,:,:]', xlabel="Albedo", ylabel="") 
        
title1a = L"$\mathrm{Uncorr.}\,\mathrm{SIF\,\,retr.}$" 
title1b = L"$(\hat{\mathrm{SIF}}=0.0\,\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m})$"
title2a = L"$\mathrm{Uncorr.}\,\mathrm{SIF\,\,retr.}$"
title2b = L"$(\hat{\mathrm{SIF}}=0.0\,\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m})$"
title3a = L"$\mathrm{Corr.}\,\mathrm{SIF\,\,retr.}$"
title3b = L"$(\hat{\mathrm{SIF}}=0.0\,\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m})$"
title4a = L"$\mathrm{Corr.}\,\mathrm{SIF\,\,retr.}$"
title4b = L"$(\hat{\mathrm{SIF}}=0.0\,\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m})$"

plot(p1, p2, p3, p4, layout = l, legend = false, title = [title1a*"\n"*title1b title2a*"\n"*title2b title3a*"\n"*title3b title4a*"\n"*title4b], titlefont = font(8))
savefig("/home/sanghavi/RamanSIFgrid/plots/SIF0_retr_0.png")

l = @layout [a1 a2; b1 b2]
p1 = heatmap(ρ, sza, w1_xnSIF0_0[3,:,:]', xlabel="", ylabel="SZA [ᵒ]")
# title=L"$\mathrm{Uncorr.}\,\mathrm{SIF\,\,retr.} (\hat{\mathrm{SIF}}=0.34\,\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m})$")
p2 = heatmap(ρ, sza, w2_xnSIF0_0[3,:,:]', xlabel="", ylabel="")
# title=L"$\mathrm{Uncorr.}\,\mathrm{SIF\,\,retr.} (\hat{\mathrm{SIF}}=0.20\,\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m})$")
p3 = heatmap(ρ, sza, w1_xnSIF0_corr_0[3,:,:]', xlabel="Albedo", ylabel="SZA [ᵒ]") 

p4 = heatmap(ρ, sza, w2_xnSIF0_corr_0[3,:,:]', xlabel="Albedo", ylabel="") 
        
title1a = L"$\mathrm{Uncorr.}\,\mathrm{SIF\,\,retr.}$" 
title1b = L"$(\hat{\mathrm{SIF}}=0.0\,\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m})$"
title2a = L"$\mathrm{Uncorr.}\,\mathrm{SIF\,\,retr.}$"
title2b = L"$(\hat{\mathrm{SIF}}=0.0\,\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m})$"
title3a = L"$\mathrm{Corr.}\,\mathrm{SIF\,\,retr.}$"
title3b = L"$(\hat{\mathrm{SIF}}=0.0\,\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m})$"
title4a = L"$\mathrm{Corr.}\,\mathrm{SIF\,\,retr.}$"
title4b = L"$(\hat{\mathrm{SIF}}=0.0\,\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m})$"

plot(p1, p2, p3, p4, layout = l, legend = false, title = [title1a*"\n"*title1b title2a*"\n"*title2b title3a*"\n"*title3b title4a*"\n"*title4b], titlefont = font(8))
savefig("/home/sanghavi/RamanSIFgrid/plots/SIF0_retr_0_noisy.png")

l = @layout [a1 a2; b1 b2]
p1 = heatmap(ρ, sza, w1_xSIF1_0[3,:,:]', xlabel="", ylabel="SZA [ᵒ]")
# title=L"$\mathrm{Uncorr.}\,\mathrm{SIF\,\,retr.} (\hat{\mathrm{SIF}}=0.34\,\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m})$")
p2 = heatmap(ρ, sza, w2_xSIF1_0[3,:,:]', xlabel="", ylabel="")
# title=L"$\mathrm{Uncorr.}\,\mathrm{SIF\,\,retr.} (\hat{\mathrm{SIF}}=0.20\,\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m})$")
p3 = heatmap(ρ, sza, w1_xSIF1_corr_0[3,:,:]', xlabel="Albedo", ylabel="SZA [ᵒ]") 

p4 = heatmap(ρ, sza, w2_xSIF1_corr_0[3,:,:]', xlabel="Albedo", ylabel="") 
        
title1a = L"$\mathrm{Uncorr.}\,\mathrm{SIF\,\,retr.}$" 
title1b = L"$(\hat{\mathrm{SIF}}=0.0\,\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m})$"
title2a = L"$\mathrm{Uncorr.}\,\mathrm{SIF\,\,retr.}$"
title2b = L"$(\hat{\mathrm{SIF}}=0.0\,\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m})$"
title3a = L"$\mathrm{Corr.}\,\mathrm{SIF\,\,retr.}$"
title3b = L"$(\hat{\mathrm{SIF}}=0.0\,\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m})$"
title4a = L"$\mathrm{Corr.}\,\mathrm{SIF\,\,retr.}$"
title4b = L"$(\hat{\mathrm{SIF}}=0.0\,\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m})$"

plot(p1, p2, p3, p4, layout = l, legend = false, title = [title1a*"\n"*title1b title2a*"\n"*title2b title3a*"\n"*title3b title4a*"\n"*title4b], titlefont = font(8))
savefig("/home/sanghavi/RamanSIFgrid/plots/SIF1_retr_0.png")


l = @layout [a1 a2; b1 b2]
p1 = heatmap(ρ, sza, w1_xnSIF1_0[3,:,:]', xlabel="", ylabel="SZA [ᵒ]")
# title=L"$\mathrm{Uncorr.}\,\mathrm{SIF\,\,retr.} (\hat{\mathrm{SIF}}=0.34\,\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m})$")
p2 = heatmap(ρ, sza, w2_xnSIF1_0[3,:,:]', xlabel="", ylabel="")
# title=L"$\mathrm{Uncorr.}\,\mathrm{SIF\,\,retr.} (\hat{\mathrm{SIF}}=0.20\,\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m})$")
p3 = heatmap(ρ, sza, w1_xnSIF1_corr_0[3,:,:]', xlabel="Albedo", ylabel="SZA [ᵒ]") 

p4 = heatmap(ρ, sza, w2_xnSIF1_corr_0[3,:,:]', xlabel="Albedo", ylabel="") 
        
title1a = L"$\mathrm{Uncorr.}\,\mathrm{SIF\,\,retr.}$" 
title1b = L"$(\hat{\mathrm{SIF}}=0.0\,\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m})$"
title2a = L"$\mathrm{Uncorr.}\,\mathrm{SIF\,\,retr.}$"
title2b = L"$(\hat{\mathrm{SIF}}=0.0\,\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m})$"
title3a = L"$\mathrm{Corr.}\,\mathrm{SIF\,\,retr.}$"
title3b = L"$(\hat{\mathrm{SIF}}=0.0\,\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m})$"
title4a = L"$\mathrm{Corr.}\,\mathrm{SIF\,\,retr.}$"
title4b = L"$(\hat{\mathrm{SIF}}=0.0\,\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m})$"

plot(p1, p2, p3, p4, layout = l, legend = false, title = [title1a*"\n"*title1b title2a*"\n"*title2b title3a*"\n"*title3b title4a*"\n"*title4b], titlefont = font(8))
savefig("/home/sanghavi/RamanSIFgrid/plots/SIF1_retr_0_noisy.png")

#=         
# noisy retrieval with tabulated IEF correction
xnSIF0_corr = zeros(length(psurf), length(ρ_str), length(sza))
xnSIF1_corr = zeros(length(psurf), length(ρ_str), length(sza))
for isurf = 1:3 # 1:1 #2:2 # 
    for iρ = 1:15 #3 #1:15
        for iA = 1:14
            sza_str  = string(sza[iA])
            alb_str  = ρ_str[iρ]
            psurf_str= string(psurf[isurf])

            fname0   = "/home/sanghavi/RamanSIFgrid/raylSIF_sza"*sza_str*"_alb"*alb_str*"_psurf"*psurf_str*"hpa_nors_ABO2.dat"
            fname1   = "/home/sanghavi/RamanSIFgrid/raylSIF_sza"*sza_str*"_alb"*alb_str*"_psurf"*psurf_str*"hpa_rrs_ABO2.dat"

            fname0_0   = "/home/sanghavi/RamanSIFgrid/rayl_sza"*sza_str*"_alb"*alb_str*"_psurf"*psurf_str*"hpa_nors_ABO2.dat"
            fname1_0   = "/home/sanghavi/RamanSIFgrid/rayl_sza"*sza_str*"_alb"*alb_str*"_psurf"*psurf_str*"hpa_rrs_ABO2.dat"


            specNoRS = readdlm(fname0)
            specRRS  = readdlm(fname1)

            specNoRS0 = readdlm(fname0_0)
            specRRS0  = readdlm(fname1_0)

            radRRS_I = (specRRS[ind,2] + specRRS[ind,5]) .* (1e7./wl[ind].^2) #converting radiances from mW/m^2/sr/cm^{-1} to mW/m^2/sr/nm
            radnoRS_I = specNoRS[ind,2] .* (1e7./wl[ind].^2) #converting radiances from mW/m^2/sr/cm^{-1} to mW/m^2/sr/nm

            radRRS_I0 = (specRRS0[ind,2] + specRRS0[ind,5]) .* (1e7./wl[ind].^2) #converting radiances from mW/m^2/sr/cm^{-1} to mW/m^2/sr/nm
            radnoRS_I0 = specNoRS0[ind,2] .* (1e7./wl[ind].^2) #converting radiances from mW/m^2/sr/cm^{-1} to mW/m^2/sr/nm

            F0 = specNoRS[ind,5] .* (1e7./wl[ind].^2) #converting radiances from mW/m^2/sr/cm^{-1} to mW/m^2/sr/nm
            
            # Create synthetic noise
            # Start
            Sig = (radRRS_I .* wl[ind]) / (1e9 * h_Pl * c_l)
            SNR = sqrt.(100*Sig.^2 ./ (MaxMS*(c_bkg^2 * MaxMS/100 .+ c_pht^2*Sig)))
            noise = zeros(length(radRRS_I))
            for i=1:length(radRRS_I)
                noise[i] = (radRRS_I[i]/SNR[i])*randn()
            end
            # End
            
            # Fit SIF
            # Grid for Legendre Polynomials:
            iLeg = range(-1,1,length(ind))
            # Take noRS run as solar reference (can in principle also use raw solar spectrum, shouldn't matter)
            #Io = specNoRS[ind,5] #radnoRS_I[ind]
            # Define polynomial terms
            poly = Pl.(iLeg,(0:3)');
            # Multiply polynomial with solar spectrum
            K_ = F0 .* poly 
            # Add a constant SIF term to the Jacobian (can add Shape later)
            K0 = [K_ ones(length(ind))];
            #K = [K_ 1e12*wl[ind].^(-4) ones(length(ind))];
            K = [K_  ones(length(ind))];
            #K = [K_  I_ie_mean ones(length(ind))];
            K1 = [K_ wl[ind] ones(length(ind))];
            
            # Fit with Least squares:
            #x1 = K1 \ radRRS_I;
            #xn1 = K1 \ (radRRS_I+noise)

            #x = K \ radRRS_I;
            xn1_corr = K \ (radRRS_I + noise/4)# - I_ie_mean*x_ie[isurf, iρ, iA])
            #xn1_corr = K1 \ (radRRS_I + noise/4 - I_ie_mean*x_ie[isurf, iρ, iA])
            #x_noRS = K \ radnoRS_I;
            #x0 = K0 \ radRRS_I;
            xn0_corr = K0 \ (radRRS_I + noise/4 - I_ie_mean*x_ie[isurf, iρ, iA])
            #x0_noRS = K0 \ radnoRS_I;
            
            #x0 = K_ \ radRRS_I0;
            #x0_noRS = K_ \ radnoRS_I0;
            #xSIF1[isurf, iρ, iA] = mean(x1[end].+x1[end-1]*wl[ind])
            #x1_ie[isurf, iρ, iA] = x1[end-2]
            xnSIF1_corr[isurf, iρ, iA] = xn1_corr[end]
            #xnSIF1_corr[isurf, iρ, iA] = mean(xn1_corr[end].+xn1_corr[end-1]*wl[ind])
            #xn1_ie[isurf, iρ, iA] = xn1[end-2]

            #xSIF[isurf, iρ, iA] = x[end]
            #x_ie[isurf, iρ, iA] = x[end-1]
            xnSIF0_corr[isurf, iρ, iA] = xn0_corr[end]
            #xn_ie[isurf, iρ, iA] = xn[end-1]

            #xSIF0[isurf, iρ, iA] = x0[end]
            #xnSIF0[isurf, iρ, iA] = xn0[end]
            #xSIF_noRS[isurf, iρ, iA] = x_noRS[end]
            #coeffa[:, isurf, iρ, iA] = x
            #coeffa_noRS[:, isurf, iρ, iA] = x_noRS
            #coeff[:, isurf, iρ, iA] = x0
            #coeff_noRS[:, isurf, iρ, iA] = x0_noRS
            #xSIF0[isurf, iρ, iA] = x0[end]
            #xSIF0_noRS[isurf, iρ, iA] = x0_noRS[end]
        end
    end
end


for i=0:1
    if i==0
        heatmap(ρ, sza, xnSIF0_corr[3,:,:]', xlabel="Albedo", ylabel="SZA", title="Retr SIF (true SIF=0.36), method 0, noisy w/ ie-corr")
        savefig("SIF0_psurf1000_noisy_wIEcorr.png")
    elseif i==1    
        heatmap(ρ, sza, xnSIF1_corr[3,:,:]', xlabel="Albedo", ylabel="SZA", title="Retr SIF (true SIF=0.36), method 0, noisy w/ ie_corr")
        savefig("SIF1_psurf1000_noisy_wIEcorr.png")
    end
end


heatmap(ρ, sza, xSIF[1,:,:]')
heatmap(ρ, sza, xSIF0[1,:,:]')
heatmap(ρ, sza, xnSIF[1,:,:]')
heatmap(ρ, sza, xnSIF0[1,:,:]')
heatmap(ρ, sza, x_ie[1,:,:]')
heatmap(ρ, sza, xn_ie[1,:,:]')

for i=0:2
    if i==0
        heatmap(ρ, sza, xSIF0[3,:,:]', xlabel="Albedo", ylabel="SZA", title="Retrieved SIF (true SIF=0.36), method 0, no noise")
        savefig("SIF0_psurf1000_nonoise.png")
        heatmap(ρ, sza, xnSIF0[3,:,:]', xlabel="Albedo", ylabel="SZA", title="Retrieved SIF (true SIF=0.36), method 0, noisy")
        savefig("SIF0_psurf1000_noisy.png")
    elseif i==1
        # SIF
        heatmap(ρ, sza, xSIF[3,:,:]', xlabel="Albedo", ylabel="SZA", title="Retrieved SIF (true SIF=0.36), method 1, no noise")
        savefig("SIF1_psurf1000_nonoise.png")
        heatmap(ρ, sza, xnSIF[3,:,:]', xlabel="Albedo", ylabel="SZA", title="Retrieved SIF (true SIF=0.36), method 1, noisy")
        savefig("SIF1_psurf1000_noisy.png")
        # IEF
        heatmap(ρ, sza, x_ie[3,:,:]', xlabel="Albedo", ylabel="SZA", title="Inelastic factor, method 1, no noise")
        savefig("IEF1_psurf1000_nonoise.png")
        heatmap(ρ, sza, xn_ie[3,:,:]', xlabel="Albedo", ylabel="SZA", title="Inelastic factor, method 1, noisy")
        savefig("IEF1_psurf1000_noisy.png")
    elseif i==2
        # SIF
        heatmap(ρ, sza, xSIF1[3,:,:]', xlabel="Albedo", ylabel="SZA", title="Retrieved SIF (true SIF=0.36), method 2, no noise")
        savefig("SIF2_psurf1000_nonoise.png")
        heatmap(ρ, sza, xnSIF1[3,:,:]', xlabel="Albedo", ylabel="SZA", title="Retrieved SIF (true SIF=0.36), method 2, noisy")
        savefig("SIF2_psurf1000_noisy.png")
        # IEF
        heatmap(ρ, sza, x1_ie[3,:,:]', xlabel="Albedo", ylabel="SZA", title="Inelastic factor, method 2, no noise")
        savefig("IEF2_psurf1000_nonoise.png")
        heatmap(ρ, sza, xn1_ie[3,:,:]', xlabel="Albedo", ylabel="SZA", title="Inelastic factor, method 2, noisy")
        savefig("IEF2_psurf1000_noisy.png")
    end
end
=#
    #geoSIF_itp = interpolate(xSIF[:,:,end:-1:1], BSpline(Cubic(Line(OnGrid()))))
#geoSIF_sitp = scale(geoSIF_itp, 500.0:250:1000, 0:0.05:0.7, 0.35:0.05:1.0)

#geoSIF_noRS_itp = interpolate(xSIF_noRS[:,:,end:-1:1], BSpline(Cubic(Line(OnGrid()))))
#geoSIF_noRS_sitp = scale(geoSIF_noRS_itp, 500.0:250:1000, 0:0.05:0.7, 0.35:0.05:1.0)


#geo_SIF_itp = []
#geo_SIF_sitp = []
#geo_corrSIF_itp = []
#geo_corrSIF_sitp = []


# Window 1, no noise
w1_geo_SIF0_itp = interpolate(w1_xSIF0[:,:,end:-1:1], BSpline(Cubic(Line(OnGrid()))))    
w1_geo_SIF0_sitp = scale(w1_geo_SIF0_itp, 500.0:250:1000, 0:0.05:0.7, 0.35:0.05:1.0)
  
w1_geo_corrSIF0_itp = interpolate(w1_xSIF0_corr[:,:,end:-1:1], BSpline(Cubic(Line(OnGrid()))))    
w1_geo_corrSIF0_sitp = scale(w1_geo_corrSIF0_itp, 500.0:250:1000, 0:0.05:0.7, 0.35:0.05:1.0)
  
w1_geo_SIF00_itp = interpolate(w1_xSIF0_0[:,:,end:-1:1], BSpline(Cubic(Line(OnGrid()))))    
w1_geo_SIF00_sitp = scale(w1_geo_SIF00_itp, 500.0:250:1000, 0:0.05:0.7, 0.35:0.05:1.0)
  
w1_geo_corrSIF00_itp = interpolate(w1_xnSIF0_corr_0[:,:,end:-1:1], BSpline(Cubic(Line(OnGrid()))))    
w1_geo_corrSIF00_sitp = scale(w1_geo_corrSIF00_itp, 500.0:250:1000, 0:0.05:0.7, 0.35:0.05:1.0)
  
# Window 1, noisy
w1n_geo_SIF0_itp = interpolate(w1_xnSIF0[:,:,end:-1:1], BSpline(Cubic(Line(OnGrid()))))    
w1n_geo_SIF0_sitp = scale(w1n_geo_SIF0_itp, 500.0:250:1000, 0:0.05:0.7, 0.35:0.05:1.0)
  
w1n_geo_corrSIF0_itp = interpolate(w1_xnSIF0_corr[:,:,end:-1:1], BSpline(Cubic(Line(OnGrid()))))    
w1n_geo_corrSIF0_sitp = scale(w1n_geo_corrSIF0_itp, 500.0:250:1000, 0:0.05:0.7, 0.35:0.05:1.0)
  
w1n_geo_SIF00_itp = interpolate(w1_xnSIF0_0[:,:,end:-1:1], BSpline(Cubic(Line(OnGrid()))))    
w1n_geo_SIF00_sitp = scale(w1n_geo_SIF00_itp, 500.0:250:1000, 0:0.05:0.7, 0.35:0.05:1.0)
  
w1n_geo_corrSIF00_itp = interpolate(w1_xnSIF0_corr_0[:,:,end:-1:1], BSpline(Cubic(Line(OnGrid()))))    
w1n_geo_corrSIF00_sitp = scale(w1n_geo_corrSIF00_itp, 500.0:250:1000, 0:0.05:0.7, 0.35:0.05:1.0)
  
# Window 2, no noise
w2_geo_SIF0_itp = interpolate(w2_xSIF0[:,:,end:-1:1], BSpline(Cubic(Line(OnGrid()))))    
w2_geo_SIF0_sitp = scale(w2_geo_SIF0_itp, 500.0:250:1000, 0:0.05:0.7, 0.35:0.05:1.0)
  
w2_geo_corrSIF0_itp = interpolate(w2_xSIF0_corr[:,:,end:-1:1], BSpline(Cubic(Line(OnGrid()))))    
w2_geo_corrSIF0_sitp = scale(w2_geo_corrSIF0_itp, 500.0:250:1000, 0:0.05:0.7, 0.35:0.05:1.0)
  
w2_geo_SIF00_itp = interpolate(w2_xSIF0_0[:,:,end:-1:1], BSpline(Cubic(Line(OnGrid()))))    
w2_geo_SIF00_sitp = scale(w2_geo_SIF00_itp, 500.0:250:1000, 0:0.05:0.7, 0.35:0.05:1.0)
  
w2_geo_corrSIF00_itp = interpolate(w2_xnSIF0_corr_0[:,:,end:-1:1], BSpline(Cubic(Line(OnGrid()))))    
w2_geo_corrSIF00_sitp = scale(w2_geo_corrSIF00_itp, 500.0:250:1000, 0:0.05:0.7, 0.35:0.05:1.0)
  
# Window 1, noisy
w2n_geo_SIF0_itp = interpolate(w2_xnSIF0[:,:,end:-1:1], BSpline(Cubic(Line(OnGrid()))))    
w2n_geo_SIF0_sitp = scale(w2n_geo_SIF0_itp, 500.0:250:1000, 0:0.05:0.7, 0.35:0.05:1.0)
  
w2n_geo_corrSIF0_itp = interpolate(w2_xnSIF0_corr[:,:,end:-1:1], BSpline(Cubic(Line(OnGrid()))))    
w2n_geo_corrSIF0_sitp = scale(w2n_geo_corrSIF0_itp, 500.0:250:1000, 0:0.05:0.7, 0.35:0.05:1.0)
  
w2n_geo_SIF00_itp = interpolate(w2_xnSIF0_0[:,:,end:-1:1], BSpline(Cubic(Line(OnGrid()))))    
w2n_geo_SIF00_sitp = scale(w2n_geo_SIF00_itp, 500.0:250:1000, 0:0.05:0.7, 0.35:0.05:1.0)
  
w2n_geo_corrSIF00_itp = interpolate(w2_xnSIF0_corr_0[:,:,end:-1:1], BSpline(Cubic(Line(OnGrid()))))    
w2n_geo_corrSIF00_sitp = scale(w2n_geo_corrSIF00_itp, 500.0:250:1000, 0:0.05:0.7, 0.35:0.05:1.0)
  

#geoSIF_itp = interpolate((FT(1.0)*psurf', ρ, FT(1.0)*sza), xSIF, BSpline(Linear())); #Gridded(Linear()));
#geoSIF_noRS_itp = interpolate((FT(1.0)*psurf', ρ, FT(1.0)*sza), xSIF_noRS, BSpline(Linear())); #Gridded(Linear()));
#=
geo_coeffa_itp = []
geo_coeffa_sitp = []
geo_coeffa_noRS_itp = []
geo_coeffa_noRS_sitp = []

for ctr=1:5
    atmp = interpolate(coeffa[ctr, :,:,end:-1:1], BSpline(Cubic(Line(OnGrid()))))
    push!(geo_coeffa_itp, atmp)
    btmp = scale(geo_coeffa_itp[ctr], 500.0:250:1000, 0:0.05:0.7, 0.35:0.05:1.0)
    push!(geo_coeffa_sitp, btmp)

    atmp = interpolate(coeffa_noRS[ctr, :,:,end:-1:1], BSpline(Cubic(Line(OnGrid()))))
    push!(geo_coeffa_noRS_itp, atmp)
    btmp = scale(geo_coeffa_noRS_itp[ctr], 500.0:250:1000, 0:0.05:0.7, 0.35:0.05:1.0)
    push!(geo_coeffa_noRS_sitp, btmp)
    #push!(geo_coeff_noRS_itp, interpolate(coeff_noRS[ctr], :,:,end:-1:1], BSpline(Cubic(Line(OnGrid())))))
    #push!(geo_coeff_noRS_sitp, scale(geoSIF_noRS_itp[ctr], 500.0:250:1000, 0:0.05:0.7, 0.35:0.05:1.0))
end
=#
#include("/home/cfranken/code/gitHub/CarbonI/TestSZA.jl")
include("test/benchmarks/TestSZA.jl")


FT = Float64
param_set = IP.InsolationParameters(FT)

# Define local overpass times:
hours = [9, 11, 13]
minutes = [30, 0, 30]

lon = -180:180
lat = -90:90
t_w1_geo_SIF0 = zeros(12,length(lon),length(lat))
t_w2_geo_SIF0 = zeros(12,length(lon),length(lat))
t_w1_geo_corrSIF0 = zeros(12,length(lon),length(lat))
t_w2_geo_corrSIF0 = zeros(12,length(lon),length(lat))

t_w1_geo_SIF00 = zeros(12,length(lon),length(lat))
t_w2_geo_SIF00 = zeros(12,length(lon),length(lat))
t_w1_geo_corrSIF00 = zeros(12,length(lon),length(lat))
t_w2_geo_corrSIF00 = zeros(12,length(lon),length(lat))

for imonth = 1:12
    for ilon = 1:length(lon)
        for ilat = 1:length(lat)
            lat_idx = argmin(abs.(modis_lat .- lat[ilat]));
            lon_idx = argmin(abs.(modis_lon .- lon[ilon]));
            @show ilon, ilat, lon_idx, lat_idx
            plon = lon[ilon]
            plat = lat[ilat]
            palb  = raman_per_month[imonth].alb_850[lon_idx, lat_idx] #alb_850_t[i,j]
            ppsurf = raman_per_month[imonth].modis_psurf[lat_idx, lon_idx]
            pSZA = raman_per_month[imonth].sza_1330[lat_idx] #get_local_time_sza.(plat, month, day, hours[3], minutes[3]);
            if (cosd(pSZA)>=0.35) && (0.7>palb>=0) && (500<=ppsurf<=1000)
                #tSIF[i,j] = geoSIF_sitp(1000., palb, cosd(pSZA))
                #tSIF_noRS[i,j] = geoSIF_noRS_sitp(1000., palb, cosd(pSZA))

                t_w1_geo_SIF0[imonth,ilon,ilat] = w1_geo_SIF0_sitp(ppsurf, palb, cosd(pSZA))
                t_w1_geo_corrSIF0[imonth,ilon,ilat] = w1_geo_corrSIF0_sitp(ppsurf, palb, cosd(pSZA))
                t_w2_geo_SIF0[imonth,ilon,ilat] = w2_geo_SIF0_sitp(ppsurf, palb, cosd(pSZA))
                t_w2_geo_corrSIF0[imonth,ilon,ilat] = w2_geo_corrSIF0_sitp(ppsurf, palb, cosd(pSZA))
                
                t_w1_geo_SIF00[imonth,ilon,ilat] = w1_geo_SIF00_sitp(ppsurf, palb, cosd(pSZA))
                t_w1_geo_corrSIF00[imonth,ilon,ilat] = w1_geo_corrSIF00_sitp(ppsurf, palb, cosd(pSZA))
                t_w2_geo_SIF00[imonth,ilon,ilat] = w2_geo_SIF00_sitp(ppsurf, palb, cosd(pSZA))
                t_w2_geo_corrSIF00[imonth,ilon,ilat] = w2_geo_corrSIF00_sitp(ppsurf, palb, cosd(pSZA))
                #=
                for ictr=1:5
                    tcoeffa[ictr,i,j] = geo_coeffa_sitp[ictr](1000., palb, cosd(pSZA))
                    tcoeffa_noRS[ictr,i,j] = geo_coeffa_noRS_sitp[ictr](1000., palb, cosd(pSZA))
                end
                =#
                #end
            end
        end
    end
end
for t=1:12
    l = @layout [a1 a2; b1 b2]
    x = convert(Vector{Float32}, lon)
    y = convert(Vector{Float32}, lat)

    z = convert(Array{Float32}, 100*(t_w1_geo_SIF0[t,:,:]/0.34 .- 1)')
    p1 = heatmap(x, y, z, size=(800,400), grid=true, xlabel="", ylabel="Lon [ᵒ]")
    #p1 = heatmap(lon, lat, 100*(t_w1_geo_SIF0[t,:,:]/0.34 .- 1)', size=(800,400), grid=true, xlabel="", ylabel="Lon [ᵒ]")
    
    # title=L"$\mathrm{Uncorr.}\,\mathrm{SIF\,\,retr.} (\hat{\mathrm{SIF}}=0.34\,\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m})$")
    z = convert(Array{Float32}, 100*(t_w2_geo_SIF0[t,:,:]/0.2 .- 1)')
    p2 = heatmap(x, y, z, size=(800,400), grid=true, xlabel="", ylabel="")
    #p2 = heatmap(lon, lat, 100*(t_w2_geo_SIF0[t,:,:]/0.20 .- 1)', size=(800,400), grid=true, xlabel="", ylabel="")
    #heatmap(ρ, sza, w2_xnSIF1_0[3,:,:]', xlabel="", ylabel="")
    # title=L"$\mathrm{Uncorr.}\,\mathrm{SIF\,\,retr.} (\hat{\mathrm{SIF}}=0.20\,\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m})$")
    z = convert(Array{Float32}, 100*(t_w1_geo_corrSIF0[t,:,:]/0.34 .- 1)')
    p3 = heatmap(x, y, z, size=(800,400), grid=true, xlabel="Lat [ᵒ]", ylabel="Lon [ᵒ]")
    p3 = heatmap(lon, lat, 100*(t_w1_geo_corrSIF0[t,:,:]/0.34 .- 1)', size=(800,400), grid=true, xlabel="Lat [ᵒ]", ylabel="Lon [ᵒ]")
    #heatmap(ρ, sza, w1_xnSIF1_corr_0[3,:,:]', xlabel="Albedo", ylabel="SZA [ᵒ]") 
    z = convert(Array{Float32}, 100*(t_w2_geo_corrSIF0[t,:,:]/0.2 .- 1)')
    p4 = heatmap(x, y, z, size=(800,400), grid=true, xlabel="Lat [ᵒ]", ylabel="")
    #p4 = heatmap(lon, lat, 100*(t_w2_geo_corrSIF0[t,:,:]/0.20 .- 1)', size=(800,400), grid=true, xlabel="Lat [ᵒ]", ylabel="")
    #heatmap(ρ, sza, w2_xnSIF1_corr_0[3,:,:]', xlabel="Albedo", ylabel="") 
            
    title1a = L"$\mathrm{Uncorr.}\,\Delta\mathrm{SIF\,\,retr.} [\%]$" 
    title1b = L"$(\hat{\mathrm{SIF}}=0.34\,\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m})$"
    title2a = L"$\mathrm{Uncorr.}\,\Delta\mathrm{SIF\,\,retr.} [\%]$"
    title2b = L"$(\hat{\mathrm{SIF}}=0.20\,\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m})$"
    title3a = L"$\mathrm{Corr.}\,\Delta\mathrm{SIF\,\,retr.} [\%]$"
    title3b = L"$(\hat{\mathrm{SIF}}=0.34\,\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m})$"
    title4a = L"$\mathrm{Corr.}\,\Delta\mathrm{SIF\,\,retr.} [\%]$"
    title4b = L"$(\hat{\mathrm{SIF}}=0.20\,\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m})$"

    plot(p1, p2, p3, p4, layout = l, legend = false, title = [title1a*"\n"*title1b title2a*"\n"*title2b title3a*"\n"*title3b title4a*"\n"*title4b], titlefont = font(8))
    savefig("/home/sanghavi/RamanSIFgrid/plots/globmap_SIF0_retr"*string(t)*".png")
end

for t=1:12
    l = @layout [a1 a2; b1 b2]
    p1 = heatmap(lon, lat, (t_w1_geo_SIF00[t,:,:])', size=(800,400), grid=true, xlabel="", ylabel="Lon [ᵒ]")
    # title=L"$\mathrm{Uncorr.}\,\mathrm{SIF\,\,retr.} (\hat{\mathrm{SIF}}=0.34\,\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m})$")
    p2 = heatmap(lon, lat, (t_w2_geo_SIF00[t,:,:])', size=(800,400), grid=true, xlabel="", ylabel="")
    #heatmap(ρ, sza, w2_xnSIF1_0[3,:,:]', xlabel="", ylabel="")
    # title=L"$\mathrm{Uncorr.}\,\mathrm{SIF\,\,retr.} (\hat{\mathrm{SIF}}=0.20\,\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m})$")
    p3 = heatmap(lon, lat, (t_w1_geo_corrSIF00[t,:,:])', size=(800,400), grid=true, xlabel="Lat [ᵒ]", ylabel="Lon [ᵒ]")
    #heatmap(ρ, sza, w1_xnSIF1_corr_0[3,:,:]', xlabel="Albedo", ylabel="SZA [ᵒ]") 
    p4 = heatmap(lon, lat, (t_w2_geo_corrSIF00[t,:,:])', size=(800,400), grid=true, xlabel="Lat [ᵒ]", ylabel="")
    #heatmap(ρ, sza, w2_xnSIF1_corr_0[3,:,:]', xlabel="Albedo", ylabel="") 
            
    title1a = L"$\mathrm{Uncorr.}\,\mathrm{SIF\,\,retr.} $" 
    title1b = L"$(\hat{\mathrm{SIF}}=0.0\,\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m})$"
    title2a = L"$\mathrm{Uncorr.}\,\mathrm{SIF\,\,retr.} $"
    title2b = L"$(\hat{\mathrm{SIF}}=0.0\,\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m})$"
    title3a = L"$\mathrm{Corr.}\,\mathrm{SIF\,\,retr.} $"
    title3b = L"$(\hat{\mathrm{SIF}}=0.0\,\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m})$"
    title4a = L"$\mathrm{Corr.}\,\mathrm{SIF\,\,retr.} $"
    title4b = L"$(\hat{\mathrm{SIF}}=0.0\,\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m})$"

    plot(p1, p2, p3, p4, layout = l, legend = false, title = [title1a*"\n"*title1b title2a*"\n"*title2b title3a*"\n"*title3b title4a*"\n"*title4b], titlefont = font(8))
    savefig("/home/sanghavi/RamanSIFgrid/plots/globmap_SIF00_retr"*string(t)*".png")
end
#heatmap(lon, lat, alb_850_Jan')
heatmap(lon, lat, t_w1_geo_corrSIF0', size=(400,400), grid=true)
            # COmpute fitted spectrum:
            #radRRS_I_fit = K * x;
            # Fitted SIF:
            #SIF_fit = x[end]
            #relative_SIF = SIF_fit / maximum(radRRS_I[ind])

            # Fit SIF without Raman
            # Fit with Least squares:
            #x_noRS = K \ radnoRS_I[ind];

            # COmpute fitted spectrum:
            #radnoRS_I_fit0 = K * x_noRS;
            # Fitted SIF:
            #SIF_fit_noRS = x_noRS[end]
            #relative_SIF_noRS = SIF_fit_noRS / maximum(radnoRS_I[ind])

            #xSIF[iρ, iA] = relative_SIF # SIF_fit
            #xSIF_noRS[iρ, iA] = relative_SIF_noRS #SIF_fit_noRS
#        end
#    end
#end

#heatmap(x, x, f, c = :thermal)
heatmap(ρ, sza, xSIF', c = :thermal)

#N = 100000
plot(scattergl(
    x=ρ[1:13], y=sza, mode="markers",
    marker=attr(color=xSIF, colorscale="Viridis", line_width=1)
))

l = @layout [a ; b]
p1 = plot(wl[ind], radRRS_I[ind], label="Measured")
plot!(wl[ind], radRRS_I_fit, label="Fitted")
p2 = plot(wl[ind], (radRRS_I[ind] .- radRRS_I_fit)./radRRS_I[ind]*1000, label="Residuals", xlabel="Wavelength (nm)", ylabel="Residuals  in permille")
plot(p1, p2,  layout = l)



# Fit without SIF:
#x2 = K_ \ radRRS_I[ind];
#radStupidFit = K_ * x2;
#l = @layout [a ; b]
# Fit with Least squares:



p1 = plot(wl[ind], radRRS_I0[ind], label="Measured")
plot!(wl[ind], radRRS_I_fit0, label="Fit without SIF")
p2 = plot(wl[ind], (radRRS_I0[ind] .- radRRS_I_fit0)./radRRS_I0[ind]*1000, label="Residuals", xlabel="Wavelength (nm)", ylabel="Residuals  in permille")
plot(p1, p2,  layout = l)

#plot!(wl[ind], radStupidFit, label="Fit without SIF")
#p2 = plot(wl[ind], (radRRS_I[ind] .- radStupidFit)./radRRS_I[ind]*1000, label="Residuals", xlabel="Wavelength (nm)", ylabel="Residuals  in permille")
#plot(p1, p2,  layout = l)
 

# Fit without Raman:
#x2 = K_ \ radRRS_I[ind];
#radStupidFit = K_ * x2;
#l = @layout [a ; b]
# Fit with Least squares:
#x_noRS = K \ radnoRS_I[ind];

# COmpute fitted spectrum:
#radnoRS_I_fit = K * x_noRS;


p1 = plot(wl[ind], radnoRS_I[ind], label="Measured")
plot!(wl[ind], radnoRS_I_fit, label="Fit without SIF")
p2 = plot(wl[ind], (radnoRS_I[ind] .- radnoRS_I_fit)./radnoRS_I[ind]*1000, label="Residuals", xlabel="Wavelength (nm)", ylabel="Residuals  in permille")
plot(p1, p2,  layout = l)

#plot!(wl[ind], radStupidFit, label="Fit without SIF")
#p2 = plot(wl[ind], (radRRS_I[ind] .- radStupidFit)./radRRS_I[ind]*1000, label="Residuals", xlabel="Wavelength (nm)", ylabel="Residuals  in permille")
#plot(p1, p2,  layout = l)
 
#x_noRS0 = K \ radnoRS_I0[ind];

# COmpute fitted spectrum:
#radnoRS_I_fit0 = K * x_noRS0;


p1 = plot(wl[ind], radnoRS_I0[ind], label="Measured")
plot!(wl[ind], radnoRS_I_fit0, label="Fit without SIF")
p2 = plot(wl[ind], (radnoRS_I0[ind] .- radnoRS_I_fit0)./radnoRS_I0[ind]*1000, label="Residuals", xlabel="Wavelength (nm)", ylabel="Residuals  in permille")
plot(p1, p2,  layout = l)