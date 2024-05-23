using Plots 
using LegendrePolynomials
using DelimitedFiles
using Insolation, Dates
import Insolation.Parameters as IP
import ClimaParams as CP
using NCDatasets
using Missings
using Interpolations
using Statistics
using InstrumentOperator
using LaTeXStrings
#using PlotlyJS
using HTTP
using JSON

path = "/home/sanghavi/RamanSIFgrid/plots/" #for plot output

function get_elevation(lat, lon)
    # Construct the API URL
    url = "https://api.open-elevation.com/api/v1/lookup?locations=$(lat),$(lon)"
    
    # Make the HTTP GET request
    response = HTTP.get(url)
    
    # Parse the JSON response
    data = JSON.parse(String(response.body))
    
    # Extract the elevation
    elevation = data["results"][1]["elevation"]
    return elevation
end

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
fname0   = "/home/sanghavi/RamanSIFgrid/rayl_sza"*sza_str*"_alb"*alb_str*"_psurf"*psurf_str*"hpa_nors_ABO2.dat"
specNoRS = readdlm(fname0) 
specNoRS = specNoRS[end:-1:1,:]
wl_ = 1e7./specNoRS[:,1]

ind_to_elim = findall(x->x>1e-6, wl_[1:end-2]-wl_[3:end]) #indices of duplicate entries (to be eliminated)
wl = wl_[ filter(x->!(x in ind_to_elim), eachindex(wl_)) ]
refI₀ = specNoRS[ filter(x->!(x in ind_to_elim), eachindex(specNoRS[:,2])), 2]
refI₀ .*= 1e7./wl.^2
refF₀ = specNoRS[ filter(x->!(x in ind_to_elim), eachindex(specNoRS[:,5])), 5]
refF₀ .*= 1e7./wl.^2
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


            #specNoRS = readdlm(fname0)
            specRRS_  = readdlm(fname1)
            specRRS_ = specRRS_[end:-1:1,:]
            specRRS = zeros(length(wl), size(specRRS_,2))
            for i=1:size(specRRS_,2)
                specRRS[:,i] = specRRS_[ filter(x->!(x in ind_to_elim), eachindex(specRRS_[:,i])), i]
            end
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
CList = reshape( range(colorant"blue", stop=colorant"red",length=n), 1, n );
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
CList = reshape( range(colorant"blue", stop=colorant"red",length=n), 1, n );
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
savefig(path*"ieI_wrt_rho_res.png")
#Determine scaled mean inelastic component
I_ie_mean = zeros(length(ind1_hires))
#ocoI1_ie_mean = zeros(length(ind1))
#ocoI2_ie_mean = zeros(length(ind2))
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
            #specNoRS0 = readdlm(fname0_0)
            #specRRS0  = readdlm(fname1_0)

            #wl = 1e7./specRRS[:,1] # converting from wn [cm^{-1}] to wl [nm]
            #ind = findall(x->x>fit_window[1] && x<fit_window[2], wl);
            
            I_ie_mean .+= (specRRS[ind1_hires,5].* (1e7./wl[ind1_hires].^2))/(specRRS[ind1_hires[1],5].* (1e7./wl[ind1_hires[1]].^2)) #converting radiances from mW/m^2/sr/cm^{-1} to mW/m^2/sr/nm
        end
    end
end
I_ie_mean/=(3*15*14)
#ocoI1_ie_mean/=(3*15*14)
#ocoI2_ie_mean/=(3*15*14)

xSIF0_corr = zeros(2,length(psurf), length(ρ_str), length(sza))
xnSIF0_corr = zeros(2,length(psurf), length(ρ_str), length(sza))
xSIF1_corr = zeros(2,length(psurf), length(ρ_str), length(sza))
xnSIF1_corr = zeros(2,length(psurf), length(ρ_str), length(sza))

xSIF0 = zeros(2,length(psurf), length(ρ_str), length(sza))
xnSIF0 = zeros(2,length(psurf), length(ρ_str), length(sza))
xSIF1 = zeros(2,length(psurf), length(ρ_str), length(sza))
xnSIF1 = zeros(2,length(psurf), length(ρ_str), length(sza))

xSIF0_corr_0 = zeros(2,length(psurf), length(ρ_str), length(sza))
xnSIF0_corr_0 = zeros(2,length(psurf), length(ρ_str), length(sza))
xSIF1_corr_0 = zeros(2,length(psurf), length(ρ_str), length(sza))
xnSIF1_corr_0 = zeros(2,length(psurf), length(ρ_str), length(sza))

xSIF0_0 = zeros(2,length(psurf), length(ρ_str), length(sza))
xnSIF0_0 = zeros(2,length(psurf), length(ρ_str), length(sza))
xSIF1_0 = zeros(2,length(psurf), length(ρ_str), length(sza))
xnSIF1_0 = zeros(2,length(psurf), length(ρ_str), length(sza))

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

            specRRS_  = readdlm(fname1_0)
            specRRS_ = specRRS_[end:-1:1,:]
            specRRS0 = zeros(length(wl), size(specRRS_,2))
            for i=1:size(specRRS_,2)
                specRRS0[:,i] = specRRS_[ filter(x->!(x in ind_to_elim), eachindex(specRRS_[:,i])), i]
            end

            # Uncorrected SIF I
            ocoI_tmp = (ocoMM[1]*(specRRS[:,2]+specRRS[:,5]) + 
                        ocoMM[2]*(specRRS[:,3]+specRRS[:,6]) +
                        ocoMM[3]*(specRRS[:,4]+specRRS[:,7])).* 
                        (1e7./wl.^2) #mW/m^2/s/nm
            interp_ocoI = LinearInterpolation((wl)*1e-3, (ocoI_tmp));  
            # Re-interpolate I from ν_grid to new grid/resolution
            ocoI_tmp1 = interp_ocoI(oco_wl);
            # Convolve input spectrum with variable kernel
            ocoI_tmp2 = InstrumentOperator.conv_spectra(oco_sounding.ils[iBand], oco_wl, ocoI_tmp1)

            # Corrected SIF I
            ocoI_corr_tmp = (ocoMM[1]*(specRRS[:,2]) + 
                        ocoMM[2]*(specRRS[:,3]) +
                        ocoMM[3]*(specRRS[:,4])).* 
                        (1e7./wl.^2) #mW/m^2/s/nm
            interp_ocoI_corr = LinearInterpolation((wl)*1e-3, (ocoI_corr_tmp));  
            # Re-interpolate I from ν_grid to new grid/resolution
            ocoI_corr_tmp1 = interp_ocoI_corr(oco_wl);
            # Convolve input spectrum with variable kernel
            ocoI_corr_tmp2 = InstrumentOperator.conv_spectra(oco_sounding.ils[iBand], oco_wl, ocoI_corr_tmp1)

            # Uncorrected noSIF I
            ocoI_0_tmp = (ocoMM[1]*(specRRS0[:,2]+specRRS0[:,5]) + 
                        ocoMM[2]*(specRRS0[:,3]+specRRS0[:,6]) +
                        ocoMM[3]*(specRRS0[:,4]+specRRS0[:,7])).* 
                        (1e7./wl.^2) #mW/m^2/s/nm
            interp_ocoI_0 = LinearInterpolation((wl)*1e-3, (ocoI_0_tmp));  
            # Re-interpolate I from ν_grid to new grid/resolution
            ocoI_0_tmp1 = interp_ocoI_0(oco_wl);
            # Convolve input spectrum with variable kernel
            ocoI_0_tmp2 = InstrumentOperator.conv_spectra(oco_sounding.ils[iBand], oco_wl, ocoI_0_tmp1)
            
            # Corrected noSIF I
            ocoI_0_corr_tmp = (ocoMM[1]*(specRRS0[:,2]) + 
                        ocoMM[2]*(specRRS0[:,3]) +
                        ocoMM[3]*(specRRS0[:,4])).* 
                        (1e7./wl.^2) #mW/m^2/s/nm
            interp_ocoI_0_corr = LinearInterpolation((wl)*1e-3, (ocoI_0_corr_tmp));  
            # Re-interpolate I from ν_grid to new grid/resolution
            ocoI_0_corr_tmp1 = interp_ocoI_0_corr(oco_wl);
            # Convolve input spectrum with variable kernel
            ocoI_0_corr_tmp2 = InstrumentOperator.conv_spectra(oco_sounding.ils[iBand], oco_wl, ocoI_0_corr_tmp1)

            # Reference spectrum: I₀
            oco_refI₀_tmp = (ocoMM[1]*(refI₀)).* 
                        (1e7./wl.^2) #mW/m^2/s/nm
            interp_oco_refI₀ = LinearInterpolation((wl)*1e-3, (oco_refI₀_tmp));  
            # Re-interpolate I from ν_grid to new grid/resolution
            oco_refI₀_tmp1 = interp_oco_refI₀(oco_wl);
            # Convolve input spectrum with variable kernel
            oco_refI₀_tmp2 = InstrumentOperator.conv_spectra(oco_sounding.ils[iBand], oco_wl, oco_refI₀_tmp1)

            # Reference spectrum: F₀
            oco_refF₀_tmp = (ocoMM[1]*(refF₀)).* 
                        (1e7./wl.^2) #mW/m^2/s/nm
            interp_oco_refF₀ = LinearInterpolation((wl)*1e-3, (oco_refF₀_tmp));  
            # Re-interpolate I from ν_grid to new grid/resolution
            oco_refF₀_tmp1 = interp_oco_refF₀(oco_wl);
            # Convolve input spectrum with variable kernel
            oco_refF₀_tmp2 = InstrumentOperator.conv_spectra(oco_sounding.ils[iBand], oco_wl, oco_refF₀_tmp1)
            
            oco_y = ocoI_tmp2
            oco_y_corr = ocoI_corr_tmp2
            oco_y_0 = ocoI_0_tmp2
            oco_y_0_corr = ocoI_0_corr_tmp2

            oco_refI₀ = oco_refI₀_tmp2
            oco_refF₀ = oco_refF₀_tmp2

            #=
            radRRS_I_corr = (specRRS[ind1_hires,2]) .* (1e7./wl[ind1_hires].^2) #converting radiances from mW/m^2/sr/cm^{-1} to mW/m^2/sr/nm
            
            radRRS_I = (specRRS[ind1_hires,2] + specRRS[ind1_hires,5]) .* (1e7./wl[ind1_hires].^2) #converting radiances from mW/m^2/sr/cm^{-1} to mW/m^2/sr/nm
            #radnoRS_I = specNoRS[ind1_hires,2] .* (1e7./wl[ind1_hires].^2) #converting radiances from mW/m^2/sr/cm^{-1} to mW/m^2/sr/nm

            radRRS_I0_corr = (specRRS0[ind1_hires,2]) .* (1e7./wl[ind1_hires].^2) #converting radiances from mW/m^2/sr/cm^{-1} to mW/m^2/sr/nm
            
            radRRS_I0 = (specRRS0[ind1_hires,2] + specRRS0[ind1_hires,5]) .* (1e7./wl[ind1_hires].^2) #converting radiances from mW/m^2/sr/cm^{-1} to mW/m^2/sr/nm
            #radnoRS_I0 = specNoRS0[ind1_hires,2] .* (1e7./wl[ind1_hires].^2) #converting radiances from mW/m^2/sr/cm^{-1} to mW/m^2/sr/nm

            F0 = specRRS[ind1_hires,8] .* (1e7./wl[ind1_hires].^2) #converting radiances from mW/m^2/sr/cm^{-1} to mW/m^2/sr/nm
            
            # Create synthetic noise
            # Start
            Sig = (radRRS_I .* wl[ind1_hires]) / (1e9 * h_Pl * c_l)
            SNR = sqrt.(100*Sig.^2 ./ (MaxMS*(c_bkg^2 * MaxMS/100 .+ c_pht^2*Sig)))
            noise = zeros(length(radRRS_I))
            for i=1:length(radRRS_I)
                noise[i] = (radRRS_I[i]/SNR[i])*randn()
            end

            x = K \ radRRS_I;
            xn = K \ (radRRS_I+noise)
            #x_noRS = K \ radnoRS_I;
            x0 = K0 \ radRRS_I;
            xn0 = K0 \ (radRRS_I+noise)

            x0_corr = K0 \ radRRS_I_corr;
            xn0_corr = K0 \ (radRRS_I_corr+noise)

            x1 = K1 \ radRRS_I;
            xn1 = K1 \ (radRRS_I+noise)
            #x0_noRS = K0 \ radnoRS_I;
            x_0 = K \ radRRS_I0;
            xn_0 = K \ (radRRS_I0+noise)
            #x_noRS = K \ radnoRS_I;
            x0_0 = K0 \ radRRS_I0;
            xn0_0 = K0 \ (radRRS_I0+noise)

            x0_corr0 = K0 \ radRRS_I0_corr;
            xn0_corr0 = K0 \ (radRRS_I0_corr+noise)

            x1_0 = K1 \ radRRS_I0;
            xn1_0 = K1 \ (radRRS_I0+noise)
            =#

            # Create synthetic noise
            ocoSig = (oco_y .* oco_sounding.SpectralGrid*1e3) / (1e9 * h_Pl * c_l)
            ocoSNR = sqrt.(100*ocoSig.^2 ./ (MaxMS*(c_bkg^2 * MaxMS/100 .+ c_pht^2*ocoSig)))
            oco_noise = zeros(length(ocoSig))
            for i=1:length(ocoSig)
                oco_noise[i] = (oco_y[i]/ocoSNR[i])*randn()/coAdd_N
            end

            ocoSig_0 = (oco_y_0 .* oco_sounding.SpectralGrid*1e3) / (1e9 * h_Pl * c_l)
            ocoSNR_0 = sqrt.(100*ocoSig_0.^2 ./ (MaxMS*(c_bkg^2 * MaxMS/100 .+ c_pht^2*ocoSig_0)))
            oco_noise_0 = zeros(length(ocoSig_0))
            for i=1:length(ocoSig_0)
                oco_noise_0[i] = (oco_y_0[i]/ocoSNR_0[i])*randn()/coAdd_N
            end
            # End
            
            for iwin=1:2
                win_ind=[]
                if iwin==1
                    win_ind=ind1
                elseif iwin==2
                    win_ind=ind2
                end
                # Fit SIF
                # Grid for Legendre Polynomials:
                #iLeg = range(-1,1,length(ind1_hires))
                iLeg = range(-1,1,length(win_ind))
                # Take noRS run as solar reference (can in principle also use raw solar spectrum, shouldn't matter)
                #Io = specNoRS[ind,5] #radnoRS_I[ind]
                # Define polynomial terms
                poly = Pl.(iLeg,(0:3)');
                # Multiply polynomial with solar/other reference spectrum
                if iwin==1
                    K_ = oco_refF₀[win_ind] .* poly 
                elseif iwin==2
                    K_ = oco_y_0_corr[win_ind] .* poly
                end 
                # Add a constant SIF term to the Jacobian (can add Shape later)
                K0 = [K_ ones(length(win_ind))];
                #K0 = [K_ ones(length(ind1_hires))];
                #K = [K_ I_ie_mean ones(length(ind1_hires))];
                K1 = [K_ oco_sounding.SpectralGrid[win_ind]*1e3 ones(length(win_ind))];
                #K1 = [K_ I_ie_mean wl[ind1_hires] ones(length(ind1_hires))];

                # Fit with Least squares:
                # SIF
                x0 = K0 \ oco_y[win_ind]
                xn0 = K0 \ (oco_y[win_ind]+oco_noise[win_ind])

                x1 = K1 \ oco_y[win_ind]
                xn1 = K1 \ (oco_y[win_ind]+oco_noise[win_ind])

                x0_corr = K0 \ oco_y_corr[win_ind]
                xn0_corr = K0 \ (oco_y_corr[win_ind]+oco_noise[win_ind])

                x1_corr = K1 \ oco_y_corr[win_ind]
                xn1_corr = K1 \ (oco_y_corr[win_ind]+oco_noise[win_ind])

                # No SIF
                x0_0 = K0 \ oco_y_0[win_ind];
                xn0_0 = K0 \ (oco_y_0[win_ind]+oco_noise_0[win_ind])

                x1_0 = K1 \ oco_y_0[win_ind];
                xn1_0 = K1 \ (oco_y_0[win_ind]+oco_noise_0[win_ind])

                x0_corr_0 = K0 \ oco_y_0_corr[win_ind];
                xn0_corr_0 = K0 \ (oco_y_0_corr[win_ind]+oco_noise_0[win_ind])

                x1_corr_0 = K1 \ oco_y_0_corr[win_ind];
                xn1_corr_0 = K1 \ (oco_y_0_corr[win_ind]+oco_noise_0[win_ind])

                xSIF1[iwin, isurf, iρ, iA] = mean(x1[end].+x1[end-1]*oco_sounding.SpectralGrid[win_ind]*1e3)/ocoMM[1]
                xnSIF1[iwin, isurf, iρ, iA] = mean(xn1[end].+xn1[end-1]*oco_sounding.SpectralGrid[win_ind]*1e3)/ocoMM[1]
                xSIF1_0[iwin, isurf, iρ, iA] = mean(x1_0[end].+x1_0[end-1]*oco_sounding.SpectralGrid[win_ind]*1e3)/ocoMM[1]
                xnSIF1_0[iwin, isurf, iρ, iA] = mean(xn1_0[end].+xn1_0[end-1]*oco_sounding.SpectralGrid[win_ind]*1e3)/ocoMM[1]

                xSIF1_corr[iwin, isurf, iρ, iA] = mean(x1_corr[end].+x1_corr[end-1]*oco_sounding.SpectralGrid[win_ind]*1e3)/ocoMM[1]
                xnSIF1_corr[iwin, isurf, iρ, iA] = mean(xn1_corr[end].+xn1_corr[end-1]*oco_sounding.SpectralGrid[win_ind]*1e3)/ocoMM[1]
                xSIF1_corr_0[iwin, isurf, iρ, iA] = mean(x1_corr_0[end].+x1_corr_0[end-1]*oco_sounding.SpectralGrid[win_ind]*1e3)/ocoMM[1]
                xnSIF1_corr_0[iwin, isurf, iρ, iA] = mean(xn1_corr_0[end].+xn1_corr_0[end-1]*oco_sounding.SpectralGrid[win_ind]*1e3)/ocoMM[1]

                xSIF0[iwin, isurf, iρ, iA] = x0[end]/ocoMM[1]
                xnSIF0[iwin, isurf, iρ, iA] = xn0[end]/ocoMM[1]
                xSIF0_0[iwin, isurf, iρ, iA] = x0_0[end]/ocoMM[1]
                xnSIF0_0[iwin, isurf, iρ, iA] = xn0_0[end]/ocoMM[1]

                xSIF0_corr[iwin, isurf, iρ, iA] = x0_corr[end]/ocoMM[1]
                xnSIF0_corr[iwin, isurf, iρ, iA] = xn0_corr[end]/ocoMM[1]
                xSIF0_corr_0[iwin, isurf, iρ, iA] = x0_corr_0[end]/ocoMM[1]
                xnSIF0_corr_0[iwin, isurf, iρ, iA] = xn0_corr_0[end]/ocoMM[1]
            end
        end
    end
end

# With SIF
l = @layout [a1 a2; b1 b2]
p1 = heatmap(ρ, sza, xSIF0[1,3,:,:]', xlabel="", ylabel="SZA [ᵒ]")
p2 = heatmap(ρ, sza, xSIF0[2,3,:,:]', xlabel="", ylabel="")
p3 = heatmap(ρ, sza, xSIF0_corr[1,3,:,:]', xlabel="Albedo", ylabel="SZA [ᵒ]")
p4 = heatmap(ρ, sza, xSIF0_corr[2,3,:,:]', xlabel="Albedo", ylabel="")

title1a = L"$\mathrm{Uncorr.}\,\mathrm{SIF\,\,retr.}$"
title1b = L"$(\hat{\mathrm{SIF}}=0.33\,\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m})$"
title2a = L"$\mathrm{Uncorr.}\,\mathrm{SIF\,\,retr.}$"
title2b = L"$(\hat{\mathrm{SIF}}=0.20\,\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m})$"
title3a = L"$\mathrm{Corr.}\,\mathrm{SIF\,\,retr.}$"
title3b = L"$(\hat{\mathrm{SIF}}=0.33\,\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m})$"
title4a = L"$\mathrm{Corr.}\,\mathrm{SIF\,\,retr.}$"
title4b = L"$(\hat{\mathrm{SIF}}=0.20\,\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m})$"
plot(p1, p2, p3, p4, layout = l, legend = false, plot_title="Method 0, no noise", title = [title1a*"\n"*title1b title2a*"\n"*title2b title3a*"\n"*title3b title4a*"\n"*title4b], titlefont = font(8))
savefig(path*"SIF0_retr.png")

l = @layout [a1 a2; b1 b2]
p1 = heatmap(ρ, sza, xSIF1[1,3,:,:]', xlabel="", ylabel="SZA [ᵒ]")
p2 = heatmap(ρ, sza, xSIF1[2,3,:,:]', xlabel="", ylabel="")
p3 = heatmap(ρ, sza, xSIF1_corr[1,3,:,:]', xlabel="Albedo", ylabel="SZA [ᵒ]")
p4 = heatmap(ρ, sza, xSIF1_corr[2,3,:,:]', xlabel="Albedo", ylabel="")

title1a = L"$\mathrm{Uncorr.}\,\mathrm{SIF\,\,retr.}$"
title1b = L"$(\hat{\mathrm{SIF}}=0.33\,\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m})$"
title2a = L"$\mathrm{Uncorr.}\,\mathrm{SIF\,\,retr.}$"
title2b = L"$(\hat{\mathrm{SIF}}=0.20\,\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m})$"
title3a = L"$\mathrm{Corr.}\,\mathrm{SIF\,\,retr.}$"
title3b = L"$(\hat{\mathrm{SIF}}=0.33\,\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m})$"
title4a = L"$\mathrm{Corr.}\,\mathrm{SIF\,\,retr.}$"
title4b = L"$(\hat{\mathrm{SIF}}=0.20\,\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m})$"
plot(p1, p2, p3, p4, layout = l, legend = false, plot_title="Method 1, no noise", title = [title1a*"\n"*title1b title2a*"\n"*title2b title3a*"\n"*title3b title4a*"\n"*title4b], titlefont = font(8))
savefig(path*"SIF1_retr.png")

l = @layout [a1 a2; b1 b2]
p1 = heatmap(ρ, sza, xnSIF0[1,3,:,:]', xlabel="", ylabel="SZA [ᵒ]")
p2 = heatmap(ρ, sza, xnSIF0[2,3,:,:]', xlabel="", ylabel="")
p3 = heatmap(ρ, sza, xnSIF0_corr[1,3,:,:]', xlabel="Albedo", ylabel="SZA [ᵒ]")
p4 = heatmap(ρ, sza, xnSIF0_corr[2,3,:,:]', xlabel="Albedo", ylabel="")

title1a = L"$\mathrm{Uncorr.}\,\mathrm{SIF\,\,retr.}$"
title1b = L"$(\hat{\mathrm{SIF}}=0.33\,\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m})$"
title2a = L"$\mathrm{Uncorr.}\,\mathrm{SIF\,\,retr.}$"
title2b = L"$(\hat{\mathrm{SIF}}=0.20\,\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m})$"
title3a = L"$\mathrm{Corr.}\,\mathrm{SIF\,\,retr.}$"
title3b = L"$(\hat{\mathrm{SIF}}=0.33\,\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m})$"
title4a = L"$\mathrm{Corr.}\,\mathrm{SIF\,\,retr.}$"
title4b = L"$(\hat{\mathrm{SIF}}=0.20\,\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m})$"
plot(p1, p2, p3, p4, layout = l, legend = false, plot_title="Method 0, noisy", title = [title1a*"\n"*title1b title2a*"\n"*title2b title3a*"\n"*title3b title4a*"\n"*title4b], titlefont = font(8))
savefig(path*"SIF0_noisy_retr.png")

l = @layout [a1 a2; b1 b2]
p1 = heatmap(ρ, sza, xnSIF1[1,3,:,:]', xlabel="", ylabel="SZA [ᵒ]")
p2 = heatmap(ρ, sza, xnSIF1[2,3,:,:]', xlabel="", ylabel="")
p3 = heatmap(ρ, sza, xnSIF1_corr[1,3,:,:]', xlabel="Albedo", ylabel="SZA [ᵒ]")
p4 = heatmap(ρ, sza, xnSIF1_corr[2,3,:,:]', xlabel="Albedo", ylabel="")

title1a = L"$\mathrm{Uncorr.}\,\mathrm{SIF\,\,retr.}$"
title1b = L"$(\hat{\mathrm{SIF}}=0.33\,\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m})$"
title2a = L"$\mathrm{Uncorr.}\,\mathrm{SIF\,\,retr.}$"
title2b = L"$(\hat{\mathrm{SIF}}=0.20\,\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m})$"
title3a = L"$\mathrm{Corr.}\,\mathrm{SIF\,\,retr.}$"
title3b = L"$(\hat{\mathrm{SIF}}=0.33\,\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m})$"
title4a = L"$\mathrm{Corr.}\,\mathrm{SIF\,\,retr.}$"
title4b = L"$(\hat{\mathrm{SIF}}=0.20\,\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m})$"
plot(p1, p2, p3, p4, layout = l, legend = false, plot_title="Method 1, noisy", title = [title1a*"\n"*title1b title2a*"\n"*title2b title3a*"\n"*title3b title4a*"\n"*title4b], titlefont = font(8))
savefig(path*"SIF1_noisy_retr.png")

# No SIF
l = @layout [a1 a2; b1 b2]
p1 = heatmap(ρ, sza, xSIF0_0[1,3,:,:]', xlabel="", ylabel="SZA [ᵒ]")
p2 = heatmap(ρ, sza, xSIF0_0[2,3,:,:]', xlabel="", ylabel="")
p3 = heatmap(ρ, sza, xSIF0_corr_0[1,3,:,:]', xlabel="Albedo", ylabel="SZA [ᵒ]")
p4 = heatmap(ρ, sza, xSIF0_corr_0[2,3,:,:]', xlabel="Albedo", ylabel="")

title1a = L"$\mathrm{Uncorr.}\,\mathrm{SIF\,\,retr.}$"
title1b = L"$(\hat{\mathrm{SIF}}=0.0\,\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m})$"
title2a = L"$\mathrm{Uncorr.}\,\mathrm{SIF\,\,retr.}$"
title2b = L"$(\hat{\mathrm{SIF}}=0.0\,\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m})$"
title3a = L"$\mathrm{Corr.}\,\mathrm{SIF\,\,retr.}$"
title3b = L"$(\hat{\mathrm{SIF}}=0.0\,\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m})$"
title4a = L"$\mathrm{Corr.}\,\mathrm{SIF\,\,retr.}$"
title4b = L"$(\hat{\mathrm{SIF}}=0.0\,\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m})$"
plot(p1, p2, p3, p4, layout = l, legend = false, plot_title="Method 0, no noise", title = [title1a*"\n"*title1b title2a*"\n"*title2b title3a*"\n"*title3b title4a*"\n"*title4b], titlefont = font(8))
savefig(path*"SIF0_retr_0.png")

l = @layout [a1 a2; b1 b2]
p1 = heatmap(ρ, sza, xSIF1_0[1,3,:,:]', xlabel="", ylabel="SZA [ᵒ]")
p2 = heatmap(ρ, sza, xSIF1_0[2,3,:,:]', xlabel="", ylabel="")
p3 = heatmap(ρ, sza, xSIF1_corr_0[1,3,:,:]', xlabel="Albedo", ylabel="SZA [ᵒ]")
p4 = heatmap(ρ, sza, xSIF1_corr_0[2,3,:,:]', xlabel="Albedo", ylabel="")

title1a = L"$\mathrm{Uncorr.}\,\mathrm{SIF\,\,retr.}$"
title1b = L"$(\hat{\mathrm{SIF}}=0.0\,\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m})$"
title2a = L"$\mathrm{Uncorr.}\,\mathrm{SIF\,\,retr.}$"
title2b = L"$(\hat{\mathrm{SIF}}=0.0\,\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m})$"
title3a = L"$\mathrm{Corr.}\,\mathrm{SIF\,\,retr.}$"
title3b = L"$(\hat{\mathrm{SIF}}=0.0\,\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m})$"
title4a = L"$\mathrm{Corr.}\,\mathrm{SIF\,\,retr.}$"
title4b = L"$(\hat{\mathrm{SIF}}=0.0\,\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m})$"
plot(p1, p2, p3, p4, layout = l, legend = false, plot_title="Method 1, no noise", title = [title1a*"\n"*title1b title2a*"\n"*title2b title3a*"\n"*title3b title4a*"\n"*title4b], titlefont = font(8))
savefig(path*"SIF1_retr_0.png")

l = @layout [a1 a2; b1 b2]
p1 = heatmap(ρ, sza, xnSIF0_0[1,3,:,:]', xlabel="", ylabel="SZA [ᵒ]")
p2 = heatmap(ρ, sza, xnSIF0_0[2,3,:,:]', xlabel="", ylabel="")
p3 = heatmap(ρ, sza, xnSIF0_corr_0[1,3,:,:]', xlabel="Albedo", ylabel="SZA [ᵒ]")
p4 = heatmap(ρ, sza, xnSIF0_corr_0[2,3,:,:]', xlabel="Albedo", ylabel="")

title1a = L"$\mathrm{Uncorr.}\,\mathrm{SIF\,\,retr.}$"
title1b = L"$(\hat{\mathrm{SIF}}=0.0\,\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m})$"
title2a = L"$\mathrm{Uncorr.}\,\mathrm{SIF\,\,retr.}$"
title2b = L"$(\hat{\mathrm{SIF}}=0.0\,\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m})$"
title3a = L"$\mathrm{Corr.}\,\mathrm{SIF\,\,retr.}$"
title3b = L"$(\hat{\mathrm{SIF}}=0.0\,\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m})$"
title4a = L"$\mathrm{Corr.}\,\mathrm{SIF\,\,retr.}$"
title4b = L"$(\hat{\mathrm{SIF}}=0.0\,\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m})$"
plot(p1, p2, p3, p4, layout = l, legend = false, plot_title="Method 0, noisy", title = [title1a*"\n"*title1b title2a*"\n"*title2b title3a*"\n"*title3b title4a*"\n"*title4b], titlefont = font(8))
savefig(path*"SIF0_noisy_retr_0.png")

l = @layout [a1 a2; b1 b2]
p1 = heatmap(ρ, sza, xnSIF1_0[1,3,:,:]', xlabel="", ylabel="SZA [ᵒ]")
p2 = heatmap(ρ, sza, xnSIF1_0[2,3,:,:]', xlabel="", ylabel="")
p3 = heatmap(ρ, sza, xnSIF1_corr_0[1,3,:,:]', xlabel="Albedo", ylabel="SZA [ᵒ]")
p4 = heatmap(ρ, sza, xnSIF1_corr_0[2,3,:,:]', xlabel="Albedo", ylabel="")

title1a = L"$\mathrm{Uncorr.}\,\mathrm{SIF\,\,retr.}$"
title1b = L"$(\hat{\mathrm{SIF}}=0.0\,\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m})$"
title2a = L"$\mathrm{Uncorr.}\,\mathrm{SIF\,\,retr.}$"
title2b = L"$(\hat{\mathrm{SIF}}=0.0\,\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m})$"
title3a = L"$\mathrm{Corr.}\,\mathrm{SIF\,\,retr.}$"
title3b = L"$(\hat{\mathrm{SIF}}=0.0\,\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m})$"
title4a = L"$\mathrm{Corr.}\,\mathrm{SIF\,\,retr.}$"
title4b = L"$(\hat{\mathrm{SIF}}=0.0\,\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m})$"
plot(p1, p2, p3, p4, layout = l, legend = false, plot_title="Method 1, noisy", title = [title1a*"\n"*title1b title2a*"\n"*title2b title3a*"\n"*title3b title4a*"\n"*title4b], titlefont = font(8))
savefig(path*"SIF1_noisy_retr_0.png")


#=
for i=0:1
    if i==0
        heatmap(ρ, sza, xnSIF0_corr[3,:,:]', xlabel="Albedo", ylabel="SZA", title="Retr SIF (true SIF=0.36), method 0, noisy w/ ie-corr")
        savefig("SIF0_psurf1000_noisy_wIEcorr.png")
    elseif i==1    
        heatmap(ρ, sza, xnSIF1_corr[3,:,:]', xlabel="Albedo", ylabel="SZA", title="Retr SIF (true SIF=0.36), method 0, noisy w/ ie_corr")
        savefig("SIF1_psurf1000_noisy_wIEcorr.png")
    end
end



    #geoSIF_itp = interpolate(xSIF[:,:,end:-1:1], BSpline(Cubic(Line(OnGrid()))))
#geoSIF_sitp = scale(geoSIF_itp, 500.0:250:1000, 0:0.05:0.7, 0.35:0.05:1.0)

#geoSIF_noRS_itp = interpolate(xSIF_noRS[:,:,end:-1:1], BSpline(Cubic(Line(OnGrid()))))
#geoSIF_noRS_sitp = scale(geoSIF_noRS_itp, 500.0:250:1000, 0:0.05:0.7, 0.35:0.05:1.0)
=#
geo_SIF0_itp = [] #uncorrected SIF retrievals using method 0
geo_SIF0_sitp = []
geo_corr_SIF0_itp = [] #corrected SIF retrievals using method 0
geo_corr_SIF0_sitp = []
geo_noSIF0_itp = [] #uncorrected noSIF retrievals using method 0
geo_noSIF0_sitp = []
geo_corr_noSIF0_itp = [] #corrected noSIF retrievals using method 0
geo_corr_noSIF0_sitp = []

geo_SIF1_itp = [] #uncorrected SIF retrievals using method 0
geo_SIF1_sitp = []
geo_corr_SIF1_itp = [] #corrected SIF retrievals using method 0
geo_corr_SIF1_sitp = []
geo_noSIF1_itp = [] #uncorrected noSIF retrievals using method 0
geo_noSIF1_sitp = []
geo_corr_noSIF1_itp = [] #corrected noSIF retrievals using method 0
geo_corr_noSIF1_sitp = []

for iwin=1:2
    atmp = interpolate(xSIF0[iwin, :,:,end:-1:1], BSpline(Cubic(Line(OnGrid()))))
    push!(geo_SIF0_itp, atmp)
    btmp = scale(geo_SIF0_itp[iwin], 500.0:250:1000, 0:0.05:0.7, 0.35:0.05:1.0)
    push!(geo_SIF0_sitp, btmp)

    atmp = interpolate(xSIF0_0[iwin, :,:,end:-1:1], BSpline(Cubic(Line(OnGrid()))))
    push!(geo_noSIF0_itp, atmp)
    btmp = scale(geo_noSIF0_itp[iwin], 500.0:250:1000, 0:0.05:0.7, 0.35:0.05:1.0)
    push!(geo_noSIF0_sitp, btmp)

    atmp = interpolate(xSIF1[iwin, :,:,end:-1:1], BSpline(Cubic(Line(OnGrid()))))
    push!(geo_SIF1_itp, atmp)
    btmp = scale(geo_SIF1_itp[iwin], 500.0:250:1000, 0:0.05:0.7, 0.35:0.05:1.0)
    push!(geo_SIF1_sitp, btmp)

    atmp = interpolate(xSIF1_0[iwin, :,:,end:-1:1], BSpline(Cubic(Line(OnGrid()))))
    push!(geo_noSIF1_itp, atmp)
    btmp = scale(geo_noSIF1_itp[iwin], 500.0:250:1000, 0:0.05:0.7, 0.35:0.05:1.0)
    push!(geo_noSIF1_sitp, btmp)

    atmp = interpolate(xSIF0_corr[iwin, :,:,end:-1:1], BSpline(Cubic(Line(OnGrid()))))
    push!(geo_corr_SIF0_itp, atmp)
    btmp = scale(geo_corr_SIF0_itp[iwin], 500.0:250:1000, 0:0.05:0.7, 0.35:0.05:1.0)
    push!(geo_corr_SIF0_sitp, btmp)

    atmp = interpolate(xSIF0_corr_0[iwin, :,:,end:-1:1], BSpline(Cubic(Line(OnGrid()))))
    push!(geo_corr_noSIF0_itp, atmp)
    btmp = scale(geo_corr_noSIF0_itp[iwin], 500.0:250:1000, 0:0.05:0.7, 0.35:0.05:1.0)
    push!(geo_corr_noSIF0_sitp, btmp)

    atmp = interpolate(xSIF1_corr[iwin, :,:,end:-1:1], BSpline(Cubic(Line(OnGrid()))))
    push!(geo_corr_SIF1_itp, atmp)
    btmp = scale(geo_corr_SIF1_itp[iwin], 500.0:250:1000, 0:0.05:0.7, 0.35:0.05:1.0)
    push!(geo_corr_SIF1_sitp, btmp)

    atmp = interpolate(xSIF1_corr_0[iwin, :,:,end:-1:1], BSpline(Cubic(Line(OnGrid()))))
    push!(geo_corr_noSIF1_itp, atmp)
    btmp = scale(geo_corr_noSIF1_itp[iwin], 500.0:250:1000, 0:0.05:0.7, 0.35:0.05:1.0)
    push!(geo_corr_noSIF1_sitp, btmp)
end
#geoSIF_itp = interpolate((FT(1.0)*psurf', ρ, FT(1.0)*sza), xSIF, BSpline(Linear())); #Gridded(Linear()));
#geoSIF_noRS_itp = interpolate((FT(1.0)*psurf', ρ, FT(1.0)*sza), xSIF_noRS, BSpline(Linear())); #Gridded(Linear()));
#Freeing reundant memory
xSIF0_corr = 0
xnSIF0_corr = 0
xSIF1_corr = 0
xnSIF1_corr = 0

xSIF0 = 0
xnSIF0 = 0
xSIF1 = 0
xnSIF1 = 0

xSIF0_corr_0 = 0
xnSIF0_corr_0 = 0
xSIF1_corr_0 = 0
xnSIF1_corr_0 = 0

xSIF0_0 = 0
xnSIF0_0 = 0
xSIF1_0 = 0
xnSIF1_0 = 0

geo_SIF0_itp = 0
geo_corr_SIF0_itp = 0
geo_noSIF0_itp = 0
geo_corr_noSIF0_itp = 0

geo_SIF1_itp = 0
geo_corr_SIF1_itp = 0
geo_noSIF1_itp = 0
geo_corr_noSIF1_itp = 0

FT = Float64
param_set = IP.InsolationParameters(FT)

# Define local overpass times:
hours = [9, 11, 13]
minutes = [30, 0, 30]

#Take January 1st as an example
date0 = DateTime("2000-01-01T11:58:56.816")

lat, lon = [34.15, 0.0]
lat = -90.0:90.0
date = DateTime(2020, 01, 10)
timezone = +0 # Just take local time
od = Insolation.OrbitalData()

#=
datetime = date + Dates.Hour.(hours) + Dates.Minute.(minutes)

S = solar_flux_and_cos_sza.(datetime, date0, [od], lon, lat', [param_set])

# Extract SZAs (can also just take mu):
SZAs = [rad2deg(acos(tup[2])) for tup in S]
SZA_interp = LinearInterpolation(lat, SZAs[2,:])
=#
modis_file = "/net/fluo/data2/DATASERVER/satellite/MODIS/MCD43A4.006/reprocessed/0.0833_lat-lon_8d/global/modis_MCD43A43-006.refl.00833deg_regrid.8d.2020.nc"
m_data = Dataset(modis_file, "r")

# You can check out the time index
time = m_data["time"][:]
lat = m_data["lat"][:]
lon = m_data["lon"][:]  

# Check out bands: https://developers.google.com/earth-engine/datasets/catalog/MODIS_061_MCD43A3
alb_850 = m_data["refl_band2"]
#alb_650 = m_data["refl_band1"]

# Extract data for the first time step (January)
#alb_850_Jan = alb_850[1,:,:]
#alb_650_Jan = alb_650[1,:,:]
for t=1:length(time)
    a=time[t]
    date = DateTime(parse(Int32,string(a)[1:4]),parse(Int32,string(a)[6:7]),parse(Int32,string(a)[9:10]),hours[3],minutes[3])
    tlon = 0.0
    tlat = -90.0:90.0
    datetime = date# + Dates.Hour.(hours[3]) + Dates.Minute.(minutes[3])
    timezone = +0 # Just take local time
    S = solar_flux_and_cos_sza.(datetime, date0, [od], tlon, tlat', [param_set])
    # Extract SZAs (can also just take mu):
    SZAs = [rad2deg(acos(tup[2])) for tup in S]
    SZA_interp = LinearInterpolation(tlat, SZAs[:])

    alb_850_t = alb_850[t,:,:]
    for ilon=1:size(alb_850_t,1)
        for ilat=1:size(alb_850_t,2)
            alb = alb_850_t[ilon,ilat]
            if !ismissing(alb)
                if alb>1.0
                    alb_850_t[ilon,ilat] = 0.0
                end
            end
        end
    end

    tSIF0 = zeros(2,length(lon),length(lat))
    #tSIF0_0 = zeros(2,length(lon),length(lat))
    #tSIF1 = zeros(2,length(lon),length(lat))
    #tSIF1_0 = zeros(2,length(lon),length(lat))
    tcorrSIF0 = zeros(2,length(lon),length(lat))
    #tcorrSIF0_0 = zeros(2,length(lon),length(lat))
    #tcorrSIF1 = zeros(2,length(lon),length(lat))
    #tcorrSIF1_0 = zeros(2,length(lon),length(lat))

    for i = 1:length(lon)
        #ind = findall(x -> string(x)!="missing", alb_850_Jan[i,:])
        ind = findall(x -> (!ismissing(x) && 0.0<=x<0.7), alb_850_t[i,:])

        if length(ind)>0
            for ctr=1:length(ind) # = 1:length(lat)
                j=ind[ctr]
                
                plon = lon[i]
                plat = lat[j]
                #if((i==26 && j==1808)||(i==44 && (j==1816)||(j==1935))||(i==109 && j==1791)||(i==206 && j==1944))
                #    elevation=0.0
                #else 
                try
                    elevation = get_elevation(plat, plon)
                catch
                    elevation=0.0
                end
                ppsurf = 1000.0*exp(-elevation/8000)
                #@show t,i,j, plon,plat, alb_850_t[i,j], ppsurf
                palb  = alb_850_t[i,j]
                pSZA = SZA_interp(plat);
                @show t,i,j, plon,plat, alb_850_t[i,j], ppsurf, pSZA
                #=
                if(cosd(pSZA)>=0.35)
                    #tSIF[i,j] = geoSIF_sitp(1000., palb, cosd(pSZA))
                    #tSIF_noRS[i,j] = geoSIF_noRS_sitp(1000., palb, cosd(pSZA))
                    for iwin=1:2
                        tSIF0[iwin,i,j] = geo_SIF0_sitp[iwin](ppsurf, palb, cosd(pSZA))
                        #tSIF0_0[iwin,i,j] = geo_noSIF0_sitp[iwin](ppsurf, palb, cosd(pSZA))
                        #tSIF1[iwin,i,j] = geo_SIF1_sitp[iwin](ppsurf, palb, cosd(pSZA))
                        #tSIF1_0[iwin,i,j] = geo_noSIF1_sitp[iwin](ppsurf, palb, cosd(pSZA))

                        tcorrSIF0[iwin,i,j] = geo_corr_SIF0_sitp[iwin](ppsurf, palb, cosd(pSZA))
                        #tcorrSIF0_0[iwin,i,j] = geo_corr_noSIF0_sitp[iwin](ppsurf, palb, cosd(pSZA))
                        #tcorrSIF1[iwin,i,j] = geo_corr_SIF1_sitp[iwin](ppsurf, palb, cosd(pSZA))
                        #tcorrSIF1_0[iwin,i,j] = geo_corr_noSIF1_sitp[iwin](ppsurf, palb, cosd(pSZA))
                    end
                end
                =#
            end
        end
    end
    l = @layout [a1 a2; b1 b2]
    p1 = heatmap(lon, lat, tSIF0[1,t,:,:]', size=(800,400), grid=true, xlabel="", ylabel="Lon [ᵒ]")
    # title=L"$\mathrm{Uncorr.}\,\mathrm{SIF\,\,retr.} (\hat{\mathrm{SIF}}=0.34\,\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m})$")
    p2 = heatmap(lon, lat, tSIF0[2,t,:,:]', size=(800,400), grid=true, xlabel="", ylabel="")
    #heatmap(ρ, sza, w2_xnSIF1_0[3,:,:]', xlabel="", ylabel="")
    # title=L"$\mathrm{Uncorr.}\,\mathrm{SIF\,\,retr.} (\hat{\mathrm{SIF}}=0.20\,\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m})$")
    p3 = heatmap(lon, lat, tcorrSIF0[1,t,:,:]', size=(800,400), grid=true, xlabel="Lat [ᵒ]", ylabel="Lon [ᵒ]")
    #heatmap(ρ, sza, w1_xnSIF1_corr_0[3,:,:]', xlabel="Albedo", ylabel="SZA [ᵒ]") 
    p4 = heatmap(lon, lat, tcorrSIF0[2,t,:,:]', size=(800,400), grid=true, xlabel="Lat [ᵒ]", ylabel="")
    #heatmap(ρ, sza, w2_xnSIF1_corr_0[3,:,:]', xlabel="Albedo", ylabel="") 
            
    title1a = L"$\mathrm{Uncorr.}\,\mathrm{SIF\,\,retr.}$" 
    title1b = L"$(\hat{\mathrm{SIF}}=0.33\,\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m})$"
    title2a = L"$\mathrm{Uncorr.}\,\mathrm{SIF\,\,retr.}$"
    title2b = L"$(\hat{\mathrm{SIF}}=0.20\,\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m})$"
    title3a = L"$\mathrm{Corr.}\,\mathrm{SIF\,\,retr.}$"
    title3b = L"$(\hat{\mathrm{SIF}}=0.33\,\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m})$"
    title4a = L"$\mathrm{Corr.}\,\mathrm{SIF\,\,retr.}$"
    title4b = L"$(\hat{\mathrm{SIF}}=0.20\,\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m})$"

    plot(p1, p2, p3, p4, layout = l, legend = false, title = [title1a*"\n"*title1b title2a*"\n"*title2b title3a*"\n"*title3b title4a*"\n"*title4b], titlefont = font(8))
    
    savefig(path*"globmap_SIF0_retr"*string(t)*".png")

    l = @layout [a1 a2; b1 b2]
    p1 = heatmap(lon, lat, tSIF1[1,t,:,:]', size=(800,400), grid=true, xlabel="", ylabel="Lon [ᵒ]")
    # title=L"$\mathrm{Uncorr.}\,\mathrm{SIF\,\,retr.} (\hat{\mathrm{SIF}}=0.34\,\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m})$")
    p2 = heatmap(lon, lat, tSIF1[2,t,:,:]', size=(800,400), grid=true, xlabel="", ylabel="")
    #heatmap(ρ, sza, w2_xnSIF1_0[3,:,:]', xlabel="", ylabel="")
    # title=L"$\mathrm{Uncorr.}\,\mathrm{SIF\,\,retr.} (\hat{\mathrm{SIF}}=0.20\,\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m})$")
    p3 = heatmap(lon, lat, tcorrSIF1[1,t,:,:]', size=(800,400), grid=true, xlabel="Lat [ᵒ]", ylabel="Lon [ᵒ]")
    #heatmap(ρ, sza, w1_xnSIF1_corr_0[3,:,:]', xlabel="Albedo", ylabel="SZA [ᵒ]") 
    p4 = heatmap(lon, lat, tcorrSIF1[2,t,:,:]', size=(800,400), grid=true, xlabel="Lat [ᵒ]", ylabel="")
    #heatmap(ρ, sza, w2_xnSIF1_corr_0[3,:,:]', xlabel="Albedo", ylabel="") 
            
    title1a = L"$\mathrm{Uncorr.}\,\mathrm{SIF\,\,retr.}$" 
    title1b = L"$(\hat{\mathrm{SIF}}=0.33\,\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m})$"
    title2a = L"$\mathrm{Uncorr.}\,\mathrm{SIF\,\,retr.}$"
    title2b = L"$(\hat{\mathrm{SIF}}=0.20\,\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m})$"
    title3a = L"$\mathrm{Corr.}\,\mathrm{SIF\,\,retr.}$"
    title3b = L"$(\hat{\mathrm{SIF}}=0.33\,\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m})$"
    title4a = L"$\mathrm{Corr.}\,\mathrm{SIF\,\,retr.}$"
    title4b = L"$(\hat{\mathrm{SIF}}=0.20\,\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m})$"

    plot(p1, p2, p3, p4, layout = l, legend = false, title = [title1a*"\n"*title1b title2a*"\n"*title2b title3a*"\n"*title3b title4a*"\n"*title4b], titlefont = font(8))
    path = "/home/sanghavi/RamanSIFgrid/plots/"
    savefig(path*"globmap_SIF1_retr"*string(t)*".png")
end


#heatmap(lon, lat, alb_850_Jan')
heatmap(lon, lat, tSIF', size=(400,400), grid=true)
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