using Plots 
using LegendrePolynomials
using DelimitedFiles
using Insolation, Dates
import Insolation.Parameters as IP
import ClimaParams as CP
using NCDatasets
using Interpolations
#using PlotlyJS
using HTTP
using JSON

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


fit_window = [758.0, 759.2]
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
xSIF = zeros(length(psurf), length(ρ_str), length(sza))
xSIF_noRS = zeros(length(psurf), length(ρ_str), length(sza))
isurf = 3
iρ = 1
iA = 1
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

            wl = 1e7./specNoRS[:,1]
            ind = findall(x->x>fit_window[1] && x<fit_window[2], wl);

            radRRS_I = specRRS[ind,2] + specRRS[ind,5] 
            radnoRS_I = specNoRS[ind,2] 

            #radRRS_I0 = specRRS0[ind,2] + specRRS0[ind,5] 
            #radnoRS_I0 = specNoRS0[ind,2] 

            F0 = specNoRS[ind,5] 
            # Find indices for wavelength window:
            
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
            K = [K_ ones(length(ind))];
            
            # Fit with Least squares:
            x = K \ radRRS_I;
            x_noRS = K \ radnoRS_I;
            #x0 = K \ radRRS_I0[ind];
            #x0_noRS = K \ radnoRS_I0[ind];

            xSIF[isurf, iρ, iA] = x[end]
            xSIF_noRS[isurf, iρ, iA] = x_noRS[end]
            #xSIF0[isurf, iρ, iA] = x0[end]
            #xSIF0_noRS[isurf, iρ, iA] = x0_noRS[end]
        end
    end
end

geoSIF_itp = interpolate(xSIF[:,:,end:-1:1], BSpline(Cubic(Line(OnGrid()))))
geoSIF_sitp = scale(geoSIF_itp, 500.0:250:1000, 0:0.05:0.7, 0.35:0.05:1.0)

geoSIF_noRS_itp = interpolate(xSIF_noRS[:,:,end:-1:1], BSpline(Cubic(Line(OnGrid()))))
geoSIF_noRS_sitp = scale(geoSIF_noRS_itp, 500.0:250:1000, 0:0.05:0.7, 0.35:0.05:1.0)

#geoSIF_itp = interpolate((FT(1.0)*psurf', ρ, FT(1.0)*sza), xSIF, BSpline(Linear())); #Gridded(Linear()));
#geoSIF_noRS_itp = interpolate((FT(1.0)*psurf', ρ, FT(1.0)*sza), xSIF_noRS, BSpline(Linear())); #Gridded(Linear()));

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

datetime = date + Dates.Hour.(hours) + Dates.Minute.(minutes)
S = solar_flux_and_cos_sza.(datetime, date0, [od], lon, lat', [param_set])

# Extract SZAs (can also just take mu):
SZAs = [rad2deg(acos(tup[2])) for tup in S]
SZA_interp = LinearInterpolation(lat, SZAs[2,:])

modis_file = "/net/fluo/data2/DATASERVER/satellite/MODIS/MCD43A4.006/reprocessed/0.0833_lat-lon_8d/global/modis_MCD43A43-006.refl.00833deg_regrid.8d.2020.nc"
m_data = Dataset(modis_file, "r")

# You can check out the time index
time = m_data["time"][:]
lat = m_data["lat"][:]
lon = m_data["lon"][:]  

# Check out bands: https://developers.google.com/earth-engine/datasets/catalog/MODIS_061_MCD43A3
alb_850 = m_data["refl_band2"]
alb_650 = m_data["refl_band1"]

# Extract data for the first time step (January)
alb_850_Jan = alb_850[1,:,:]
alb_650_Jan = alb_650[1,:,:]


tSIF = zeros(length(lon),length(lat))
tSIF_noRS = zeros(length(lon),length(lat))
for i = 1:length(lon)
    #ind = findall(x -> string(x)!="missing", alb_850_Jan[i,:])
    ind = findall(x -> !isa(x, Missing) && 0.0<x<0.7, alb_850_Jan[i,:])
    #@show ind
    #@show alb_850_Jan[ind,i]
    #ind = findall(x-> 0.0<x<0.7, alb_850_Jan[ind,i])
    #@show ind
    #@show alb_850_Jan[ind,i]
    if length(ind)>0
        for ctr=1:length(ind) # = 1:length(lat)
            j=ind[ctr]
            @show i,j, alb_850_Jan[i,j]
            plon = lon[i]
            plat = lat[j]
            elevation = get_elevation(plat, plon)
            ppsurf = 1000.0*exp(-elevation/8000)
            palb  = alb_850_Jan[i,j]
            pSZA = SZA_interp(plat);
            if(cosd(pSZA)>=0.35)
                tSIF[i,j] = geoSIF_sitp(ppsurf, palb, cosd(pSZA))
                tSIF_noRS[i,j] = geoSIF_noRS_sitp(ppsurf, palb, cosd(pSZA))
            end
        end
    end
end


#heatmap(lon, lat, alb_850_Jan')
heatmap(lon, lat, tSIF')
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