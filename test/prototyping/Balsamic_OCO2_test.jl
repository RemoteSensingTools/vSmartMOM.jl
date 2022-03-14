using Revise
using Plots
using Pkg
# Pkg.activate(".")
using vSmartMOM
using vSmartMOM.Architectures
using vSmartMOM.Absorption
using vSmartMOM.Scattering
using vSmartMOM.CoreRT
using vSmartMOM.SolarModel
using InstrumentOperator
using Interpolations
using Polynomials
using ForwardDiff 
using Distributions
using NCDatasets

#---------------------------------------------<Start Functions>---------------------------------------------#
function getSolar(grid, Tsolar)
    T = 5777.0
    λ_grid = reverse(1e4 ./ grid) #collect(757:0.01:777)
    black_body = planck_spectrum_wl(T, λ_grid) * 2.1629e-05 * π   
    black_body = SolarModel.watts_to_photons(λ_grid, black_body)

    # Get solar transmittance spectrum 
    #solar_transmission = solar_transmission_from_file("/home/rjeyaram/RadiativeTransfer/src/solar_merged_20160127_600_26316_100.out", parameters.spec_bands[1])
    solar_transmission = Tsolar(grid)
    # Get outgoing solar radiation
    sun_out = reverse(solar_transmission) .* black_body;
    return sun_out
end
###################################################


# Runner is used to set AD fields as duals
function runner!(y, x, parameters=parameters, oco_sounding= oco_soundings, Tsolar = Tsolar_interp)
    Nangles = length(oco_sounding);
    # Set parameters fields as the dual numbers
    parameters.brdf =  [CoreRT.LambertianSurfaceLegendre([x[1],x[3],x[4]]),
                        CoreRT.LambertianSurfaceLegendre([x[7],x[8],x[5]]),
                        CoreRT.LambertianSurfaceLegendre([x[9],x[10],x[6]])];

    parameters.scattering_params.rt_aerosols[1].τ_ref = exp(x[2]);
    parameters.scattering_params.rt_aerosols[1].p₀    = x[20]; #800.0; #x[4]
    parameters.scattering_params.rt_aerosols[1].σp    = x[21];
    parameters.scattering_params.rt_aerosols[1].aerosol.size_distribution = LogNormal(x[18], x[19]); #x[4]
    parameters.scattering_params.rt_aerosols[1].aerosol.nᵣ = x[16]
    parameters.scattering_params.rt_aerosols[1].aerosol.nᵣ = x[17]
    #parameters.p   = oco_sounding.p_half
    #parameters.q   = oco_sounding.q 
    #parameters.T   = oco_sounding.T# .+ 1.0 #.+ x[15]
    #parameters.sza = oco_sounding.sza
    #parameters.vza = [oco_sounding.vza]
    for ia=1:Nangles
        if ia==1
            parameters.p   = oco_sounding[ia].p_half
            parameters.q   = oco_sounding[ia].q 
            parameters.T   = oco_sounding[ia].T# .+ 1.0 #.+ x[15]
            parameters.sza = oco_sounding[ia].sza
            parameters.vza = [oco_sounding[ia].vza]   
            parameters.vaz = [oco_sounding[ia].raa]
        else
            parameters.vza = vcat(parameters.vza, oco_sounding[ia].vza)
            parameters.vaz = vcat(parameters.vaz, oco_sounding[ia].raa)
        end
    end
    parameters.absorption_params.vmr["H2O"] = [parameters.q[1:65]*x[11] * 1.8; 
                                               parameters.q[66:end]*x[15] * 1.8];
    a1 = zeros(7) .+ x[12]
    a2 = zeros(7) .+ x[13]
    a3 = zeros(6) .+ x[14]
    parameters.absorption_params.vmr["CO2"] = [a1; a2; a3]
    model = model_from_parameters(parameters);
    for ia = 1:Nangles
        for i = 1:length(oco_sounding[ia].BandID)
            println("Modelling band $(i)")
            # Run the model to obtain reflectance matrix
            #R = rt_run(model, i_band=i)[1];
            R = CoreRT.rt_run_test(CoreRT.noRS(), model, i)[1]
            # Re-interpolate I from ν_grid to new grid/resolution
            λ_grid = reverse(1e4 ./ parameters.spec_bands[i])
            res = 0.001e-3;
            off = 0.5e-3
            # Get sun:
            @time sun_out = getSolar(parameters.spec_bands[i],Tsolar)
        
            RR = oco_sounding[ia].mueller[1]*R[ia,1,:] + 
                    oco_sounding[ia].mueller[2]*R[ia,2,:] + 
                    oco_sounding[ia].mueller[3]*R[ia,3,:];        
        
            # Apply Earth reflectance matrix 
            earth_out = sun_out .* reverse(RR[:])
        
            # Re-interpolate I from ν_grid to new grid/resolution
            @time interp_I = LinearInterpolation(λ_grid, earth_out);        
            #The following holds only for the use of a single ground pixel (because the ils changes with ground pixel index)
            wl = oco_sounding[ia].ils[i].ν_out[1]-off:res:oco_sounding[ia].ils[i].ν_out[end]+off;
            @show wl[1],wl[end], λ_grid[1],λ_grid[end]
            I_wl = interp_I(wl);

            # Convolve input spectrum with variable kernel
            @time I_conv = InstrumentOperator.conv_spectra(oco_sounding[ia].ils[i], wl, I_wl)
            if ia==1
                y[oco_sounding[ia].BandID[i]] = [I_conv]
            else
                y[oco_sounding[ia].BandID[i]] = vcat(y[oco_sounding.BandID[i]], I_conv)
            end
        end
    end
end
#---------------------------------------------<End Functions>----------------------------------------------#

## Atmospheric Radiative Transfer

# Load parameters from file
parameters = parameters_from_yaml("test/test_parameters/3BandParameters.yaml")
#parameters.architecture = CPU()
FT = Float64

# Load OCO Data: 
# File names:

L1File   = "/net/fluo/data1/group/oco3/L1bSc/oco3_L1bScSC_15935a_220226_B10311_220226124444.h5" #"/net/fluo/data1/group/oco2/L1bSc/oco2_L1bScND_26780a_190715_B10003r_200429212407.h5"
metFile  = "/net/fluo/data1/group/oco3/L2Met/oco3_L2MetSC_15935a_220226_B10311_220226124455.h5" #"/net/fluo/data1/group/oco2/L2Met/oco2_L2MetND_26780a_190715_B10003r_200429212406.h5"
dictFile = "/home/cfranken/code/gitHub/InstrumentOperator.jl/json/oco2.yaml"
# Load L1 file (could just use filenames here as well)
oco = InstrumentOperator.load_L1(dictFile,L1File, metFile);
grnd_pxl = 5
aa = Dataset(L1File)
bb = aa.group["SoundingGeometry"];
cc = bb["sounding_pcs_mode"];
pcs = [cc[grnd_pxl,i] for i=1:size(cc,2)];
iSam = findall(pcs .== "AM")

# Pick some bands as tuple (or just one)
bands = (1,2,3);
#bands = (1,3);
# Indices within that band:
indices = (92:885,114:845,50:916);
#indices = (92:885,50:916);

ii = iSam[1:36:end]
#y=[];
#oco_soundings=[];
# Get data for that sounding:
oco_sounding = InstrumentOperator.getMeasurement(oco, bands, indices, GeoInd);
oco_soundings = [oco_sounding];
y = [oco_sounding.SpectralMeasurement];

for i=2:length(ii)
    # Geo Index (footprint,sounding):
    GeoInd = [grnd_pxl,ii[i]];
    # Get data for that sounding:
    oco_sounding = InstrumentOperator.getMeasurement(oco, bands, indices, GeoInd);
    if i==1    
        oco_soundings = [oco_sounding];
        y = [oco_sounding.SpectralMeasurement];
    else
        oco_soundings = vcat(oco_soundings, oco_sounding);
        y = vcat(y, oco_sounding.SpectralMeasurement)
    end
end
###################################################
# Produce black-body in wavenumber range
Tsolar = solar_transmission_from_file("/home/rjeyaram/vSmartMOM/src/SolarModel/solar.out")
Tsolar_interp = LinearInterpolation(Tsolar[:, 1], Tsolar[:, 2])

#Initial guess state vector 
x = FT[0.2377, -3, -0.00244, 0, 0,0, 0.4, 0, 0.4, 0,1,400e-6,400e-6,400e-6, 1.0, 1.5, 0.0, -0.69, 0.3, 690., 0.3]
# State vector box-limits
limnᵣ = [1.3,1.7]
limnᵢ = [0.0,0.2]
limσᵣ = [0.05,1.6]
limlnr₀ = [-2.3,1.6] #r₀ ≈ 0.1, 5.0  
limp₀ = [100.,1100.]
limσₚ = [0.05,1.6]
# Run FW model:
# @time runner(x);
#y = oco_sounding.SpectralMeasurement;
Fx = zeros(length(y));
#ind = 92:885

for i=1:5
    K = ForwardDiff.jacobian(runner!, Fx, x);
    dx = K \ (y-Fx);
    x += dx;

    if x[16]<limnᵣ[1]
        x[16] = limnᵣ[1];
    elseif x[16]>limnᵣ[2]
        x[16] = limnᵣ[2];
    end

    if x[17]<limnᵢ[1]
        x[17] = limnᵢ[1];
    elseif x[17]>limnᵢ[2]
        x[17] = limnᵢ[2];
    end
    
    if x[18]<limlnr₀[1]
        x[18] = limlnr₀[1];
    elseif x[18]>limlnr₀[2]
        x[18] = limlnr₀[2];
    end

    if x[19]<limσᵣ[1]
        x[19] = limσᵣ[1];
    elseif x[19]>limσᵣ[2]
        x[19] = limσᵣ[2];
    end

    if x[20]<limp₀[1]
        x[20] = limp₀[1];
    elseif x[20]>limp₀[2]
        x[20] = limp₀[2];
    end

    if x[21]<limσₚ[1]
        x[21] = limσₚ[1];
    elseif x[21]>limσₚ[2]
        x[21] = limσₚ[2];
    end
    @show x;
end
runner!(Fx,x)

ν = oco_sounding.SpectralGrid*1e3
plot(y/1e20, label="Meas")
plot!(Fx/1e20, label="Mod")
plot!((y-Fx)/1e20, label="Meas-mod")
#=
a = Dataset("/net/fluo/data1/group/oco2/oco2_L2DiaND_26780a_190715_B10004r_200618191413.h5")
g = a.group["RetrievalGeometry"]
s = a.group["SpectralParameters"]
OP_results = a.group["RetrievalResults"]
i = findall(abs.(g["retrieval_latitude"][:].-41.736668).<0.0001)
fp_meas = s["measured_radiance"][:,i]
fp_mod  = s["modeled_radiance"][:,i]
=#
#=
fil = "/net/fluo/data1/group/oco3/L1bSc/oco3_L1bScSC_15935a_220226_B10311_220226124444.h5"
aa = Dataset(fil)
bb = aa.group["SoundingGeometry"];
cc = bb["sounding_pcs_mode"];
pcs = [cc[1,i] for i=1:size(cc,2)];
iSam = findall(pcs .== "AM")
#scatter(bb["sounding_longitude"][1,iSam],bb["sounding_latitude"][1,iSam])
SZA = bb["sounding_solar_zenith"][1,iSam[1:36:end]]
SAz = bb["sounding_solar_azimuth"][1,iSam[1:36:end]]
LZA = bb["sounding_zenith"][1,iSam[1:36:end]]
LAz = bb["sounding_azimuth"][1,iSam[1:36:end]]
plot(LZA)
=#