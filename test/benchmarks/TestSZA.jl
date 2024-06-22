using Insolation, Dates, Missings, Statistics,NaNMath, Dates, Plots
import Insolation.Parameters as IP
import ClimaParams as CP
FT = Float64
param_set = IP.InsolationParameters(FT)

# Define local overpass times:
hours = [9, 11, 13]
minutes = [30, 0, 30]

#Take January 1st as an example, This is just a baseline year for which the orbital characteristics are computed
date0 = DateTime("2010-01-01T11:58:56.816")

lat, lon = [34.15, 0.0]
lat = -90.0:90.0

timezone = +0 # Just take local time
od = Insolation.OrbitalData()
#datetime = date + Dates.Hour.(hours) + Dates.Minute.(minutes)

#S = solar_flux_and_cos_sza.(datetime, date0, [od], lon, lat', [param_set])

# Extract SZAs (can also just take mu):
#SZAs = [rad2deg(acos(tup[2])) for tup in S]
#SZA_interp = LinearInterpolation(lat, SZAs[2,:])

function get_local_time_sza(lat, month, day, hours, minutes; od=od, date0=date0, param_set=param_set)
    date = DateTime(2020, month, day)
    #@show date
    datetime = date + Dates.Hour(hours) + Dates.Minute(minutes)
    S = solar_flux_and_cos_sza(datetime, date0, od, 0.0, lat, param_set)
    return rad2deg(acos(S[2]))
end

# Add topography as well:
using NCDatasets

nc_path = "/home/cfranken/code/gitHub/CarbonI/data/ETOPO_2022_v1_60s_N90W180_surface.nc"
ds = NCDataset(nc_path);
    
# Extract variables
topo_lats = ds["lat"][:];
topo_lons = ds["lon"][:];
topo_elevation = ds["z"][:];

# Function to get elevation from a netCDF file given lat/lon
function get_elevation(lat, lon; lats=topo_lats, lons=topo_lons, elevation=topo_elevation)
    # Find the closest indices for the given lat/lon
    lat_idx = argmin(abs.(lats .- lat))
    lon_idx = argmin(abs.(lons .- lon))
    
    # Get the elevation value
    elev = elevation[lon_idx, lat_idx]
    return elev
end



# Same with MODIS data:
modis_file = "/net/fluo/data2/DATASERVER/satellite/MODIS/MCD43A4.006/reprocessed/0.0833_lat-lon_8d/global/modis_MCD43A43-006.refl.00833deg_regrid.8d.2020.nc"
m_data = Dataset(modis_file, "r")

# You can check out the time index
time_modis = m_data["time"][:]
modis_lat = m_data["lat"][:]
modis_lon = m_data["lon"][:] 
alb_850 = m_data["refl_band2"][:]
alb_650 = m_data["refl_band1"][:]

#alb_850_ = coalesce.(alb_850, NaN);
#alb_650_ = coalesce.(alb_650, NaN);

i_elev_lat = [argmin(abs.(topo_lats .- iLat)) for iLat in modis_lat];
i_elev_lon = [argmin(abs.(topo_lons .- iLon)) for iLon in modis_lon];

# Get elevation and psurf using 1km elevation database:
modis_elevation = topo_elevation[i_elev_lon, i_elev_lat];
modis_psurf = exp.(-modis_elevation'/8000)*1000;

months_modis = Dates.month.(time_modis);
i1 = findall(months_modis.==1);

function mean_missing(alb, indices) 
    dim1,dim2,dim3 = size(alb)
    alb_mean = zeros(dim2,dim3)
    for i in 1:dim2
        for j in 1:dim3
            alb_mean[i,j] = mean(skipmissing(alb[indices,i,j]))
        end
    end 
    return alb_mean
end

# Define a struct with all important parameters:
struct raman_study_params
    lat
    lon
    index_month 
    alb_850
    alb_650
    modis_psurf
    sza_930
    sza_1130
    sza_1330
end

# Array of all parameters per month
raman_per_month = []
latis = convert.(Float64, modis_lat)
for imonth in 1:12
    
    iModis   = findall(months_modis.==imonth);
    alb_850_ = mean_missing(alb_850, iModis);
    alb_650_ = mean_missing(alb_650, iModis);
    sza_930  = get_local_time_sza.(latis, imonth, 15, 9, 30);
    sza_1130 = get_local_time_sza.(latis, imonth, 15, 11, 30);
    sza_1330 = get_local_time_sza.(latis, imonth, 15, 13, 30);
    raman_params = raman_study_params(modis_lat, modis_lon, imonth, alb_850_, alb_650_, modis_psurf, sza_930, sza_1130, sza_1330);
    append!(raman_per_month, [raman_params])
end


# Test stuff out:
lat = 52.5200; lon = 13.4050 # Berlin
lat = 7.000;   lon = 20.0    # Sahel zone in Africa (somewhere)

lat_idx = argmin(abs.(modis_lat .- lat));
lon_idx = argmin(abs.(modis_lon .- lon));

# get timeseries:
plot(1:12, [raman_per_month[i].sza_1130[lat_idx] for i in 1:12], label="SZA 11:30")
plot!(1:12, [raman_per_month[i].sza_1330[lat_idx] for i in 1:12], label="SZA 13:30")
plot!(1:12, [raman_per_month[i].sza_930[lat_idx] for i in 1:12], label="SZA 9:30")
plot(1:12, [raman_per_month[i].alb_850[lon_idx, lat_idx] for i in 1:12], label="A Band Albedo")
plot!(1:12, [raman_per_month[i].alb_650[lon_idx, lat_idx] for i in 1:12], label="B Band Albedo")