using Insolation, Dates
import Insolation.Parameters as IP
import ClimaParams as CP
using Plots
using NCDatasets

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
SZA_interp = LinearInterpolation(lat, SZAs)
mus  = [tup[2] for tup in S]
flux = [tup[1] for tup in S]
plot(lat, SZAs[1,:], label="9:30")
plot!(lat, SZAs[2,:], label="11:00")
plot!(lat, SZAs[3,:], label="13:30")
xlabel!("Latitude")
ylabel!("Solar Zenith Angle (deg)")

plot(lat, 1.0 ./mus[1,:], label="9:30")
plot!(lat, 1.0 ./mus[2,:], label="11:00")
plot!(lat, 1.0 ./mus[3,:], label="13:30")
xlabel!("Latitude")
ylabel!("Air mass factor (1/cos(SZA))")
ylims!(0, 10)


plot(lat, flux[1,:], label="9:30")
plot!(lat, flux[2,:], label="11:00")
plot!(lat, flux[3,:], label="13:30")
xlabel!("Latitude")
ylabel!("Solar Flux")
## Load in albedo data
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

heatmap(lon, lat, alb_850_Jan')

for i = 1:length(lon)
    ind = findall(x -> string(x) !="missing", alb_850_Jan[i,:])
    for j in ind # = 1:length(lat)
        plon = lon[i]
        plat = lat[j]
        alb  = alb_850_Jan[j,i]
        SZA = SZA_interp(plat);
    end
end
