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
using vSmartMOM.InelasticScattering
using InstrumentOperator
using Interpolations
using Polynomials
using ForwardDiff 
using Distributions
using NCDatasets
using LinearAlgebra
# using Distances
include("/home/sanghavi/code/github/vSmartMOM.jl/test/prototyping/runner.jl")
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

function run_inversion(x,y)
    # State vector box-limits
    limnᵣ = [1.3,1.7]
    limnᵢ = [0.0,0.2]
    limσᵣ = [0.05,1.6]
    limlnr₀ = [-2.3,1.6] #r₀ ≈ 0.1, 5.0  
    limp₀ = [100.,1100.]
    limσₚ = [0.05,1.6]
    for i=1:1#5
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
end

function run_OE_inversion(x, y, xa, Nangles, spec_lengths, convfct)
    Fx = zeros(length(y));
    # State vector box-limits
    limnᵣ = [1.3,1.7]
    limnᵢ = [0.0,0.2]
    limσᵣ = [0.05,1.6]
    limlnr₀ = [-2.3,1.6] #r₀ ≈ 0.1, 5.0  
    limp₀ = [100.,1100.]
    limσₚ = [0.05,1.6]
    Sa = zeros(length(x));
    Sa = [(0.5)^2, #1. A-band LambertianSurfaceLegendre
         (0.5)^2, #AOD
         (0.01)^2, #2. A-band LambertianSurfaceLegendre
         (0.001)^2, #3. A-band LambertianSurfaceLegendre
         (0.001)^2, #3. W-band LambertianSurfaceLegendre
         (0.001)^2, #3. S-band LambertianSurfaceLegendre
         (0.25)^2, #1. W-band LambertianSurfaceLegendre
         (0.01)^2, #2. W-band LambertianSurfaceLegendre
         (0.25)^2, #1. S-band LambertianSurfaceLegendre
         (0.01)^2, #2. S-band LambertianSurfaceLegendre
         (0.75)^2, #H2O 1        
         (600.e-6)^2, #CO2 1
         (475.e-6)^2, #CO2 2
         (250.e-6)^2, #CO2 3
         (0.75)^2, #H2O 2
         (limnᵣ[2]-limnᵣ[1])^2, 
         (limnᵢ[2]-limnᵢ[1])^2, 
         (limlnr₀[2]-limlnr₀[1])^2,
         (limσᵣ[2]-limσᵣ[1])^2,
         (limp₀[2]-limp₀[1])^2,
         (limσₚ[2]-limσₚ[1])^2]
    Sϵ = zeros(length(y));
    #Imax=[8.82e20, 4.33e20, 1.89e20];
    c_ph=[0.0089, 0.0069, 0.0079]; #photon noise
    c_bg=[0.0042, 0.0065, 0.0145]; #background noise

    Nwl = sum(spec_lengths)
    for k=1:Nangles
        s1 = (k-1)*sum(spec_lengths) + 1
        e1 = (k-1)*sum(spec_lengths) + spec_lengths[1]
        I1 = maximum(y[s1:e1])
    
        s2 = s1 + spec_lengths[1]
        e2 = e1 + spec_lengths[2]
        I2 = maximum(y[s2:e2])
    
        s3 = s2 + spec_lengths[2]
        e3 = e2 + spec_lengths[3]
        I3 = maximum(y[s3:e3])
    
        Imax = [I1  I2 I3]

        noise_ph = [];
        noise_bg = [];
        for i=1:length(spec_lengths)
            tmp1 = c_ph[i]^2*ones(spec_lengths[i]);
            tmp2 = 0.01*Imax[i]*ones(spec_lengths[i])*c_bg[i]^2
            push!(noise_ph,tmp1);
            push!(noise_bg,tmp2);
        end
        ϵ1=reduce(vcat,noise_ph);
        ϵ2=reduce(vcat,noise_bg);
    
        Sϵ[(k-1)*Nwl+1 : k*Nwl] .= ϵ1.*y[(k-1)*Nwl+1 : k*Nwl] .+ ϵ2;
    end
    invŜ = zeros(length(x), length(x));
    dx   = zeros(length(x));
    for i=1:1 #5
        K = ForwardDiff.jacobian(runner!, Fx, x);
        @show size(K), size(convfct)
        @show size(Fx)
        K  = K./convfct
        Fx = Fx./convfct
        #Computing invŜ 
        invŜ .= Diagonal( 1.0./Sa ) .+ K'*Diagonal( 1.0./Sϵ )*K

        dx .= invŜ * (K'*Diagonal( 1.0./Sϵ )*(y.-Fx) .- Diagonal( 1.0./Sa )*(x.-xa)); #K\(y-Fx)
        x += dx;
        @show "Unadjusted" x;
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
        @show "Adjusted" x;
    end
end
###################################################



#---------------------------------------------<End Functions>----------------------------------------------#

## Atmospheric Radiative Transfer

# Load parameters from file
parameters = parameters_from_yaml("test/test_parameters/3BandParameters.yaml")
#parameters.architecture = CPU()
FT = Float64

# Load OCO Data: 
# File names:

L1File   = "/net/fluo/data2/groupMembers/sanghavi/oco3_L1bScSC_04379a_200210_B10305r_211019003131.h5" #"/net/fluo/data1/group/oco3/L1bSc/oco3_L1bScSC_15935a_220226_B10311_220226124444.h5" #"/net/fluo/data1/group/oco2/L1bSc/oco2_L1bScND_26780a_190715_B10003r_200429212407.h5"
metFile  = "/net/fluo/data2/groupMembers/sanghavi/oco3_L2MetSC_04379a_200210_B10305r_211018214305.h5" #"/net/fluo/data1/group/oco3/L2Met/oco3_L2MetSC_15935a_220226_B10311_220226124455.h5" #"/net/fluo/data1/group/oco2/L2Met/oco2_L2MetND_26780a_190715_B10003r_200429212406.h5"
dictFile = "/home/cfranken/code/gitHub/InstrumentOperator.jl/json/oco2.yaml"
# Load L1 file (could just use filenames here as well)
oco = InstrumentOperator.load_L1(dictFile,L1File, metFile);
grnd_pxl = 1:1 #8#5
grnd_pxl_tmp=1
aa = Dataset(L1File)
bb = aa.group["SoundingGeometry"];
cc = bb["sounding_pcs_mode"];
pcs = [cc[grnd_pxl_tmp,i] for i=1:size(cc,2)];
iSam = findall(pcs .== "TG")

# Pick some bands as tuple (or just one)
bands = (1,2,3);
#bands = (1,3);
# Indices within that band:
indices = (92:885,114:845,50:916);
#indices = (92:885,50:916);

#ii = iSam[1:36:end]
ii = iSam[1:end]
#y=[];
#oco_soundings=[];
# Get data for that sounding:
#ocal GeoInd = [grnd_pxl,ii[1]];
#oco_sounding = InstrumentOperator.getMeasurement(oco, bands, indices, GeoInd);
oco_soundings0 = [];
tmpy = [];#[];#[oco_soundings[1].SpectralMeasurement];
Nangles = 1 #length(ii)
lat=[];
lon=[];
sza=[];
vza=[];
raz=[];
for i=1:Nangles
    for j=1:length(grnd_pxl)
        # Geo Index (footprint,sounding):
        GeoInd = [grnd_pxl[j],ii[i]];
        # Get data for that sounding:
        oco_sounding = InstrumentOperator.getMeasurement(oco, bands, indices, GeoInd);
        push!(oco_soundings0, oco_sounding)
        push!(lat, oco_sounding.latitude)
        push!(lon, oco_sounding.longitude)
        push!(sza, oco_sounding.sza)
        push!(vza, oco_sounding.vza)
        push!(raz, oco_sounding.vaa-oco_sounding.saa)
        @show  oco_sounding.latitude, oco_sounding.longitude
        @show oco_sounding.vza, oco_sounding.vaa-oco_sounding.saa, oco_sounding.sza
        #if i==1    
        #    #oco_soundings = [oco_sounding];
        #    y = [oco_soundings[i].SpectralMeasurement];
        #else
            #oco_soundings = vcat(oco_soundings, oco_sounding);
        push!(tmpy, oco_sounding.SpectralMeasurement)
    end
end

function dist(p1, p2)
    Rₑ = 6375
    dx = Rₑ*cosd(0.5*(p1[1]+p2[1]))*abs(p1[2]-p2[2])*π/180.
    dy = Rₑ*abs(p1[1]-p2[1])*π/180.
    d  = sqrt(dx^2+dy^2)
end
pts = [lat';lon']
Npts=length(lat)
dd=zeros(Npts, Npts)
for i=1:Npts
    for j=i:Npts
        if(i!=j)
            dd[i,j] = dist(pts[:,i],pts[:,j])
            dd[j,i] = dd[i,j]
        end
    end
end
v_Npts = zeros(Int, Npts)
for i=1:Npts
    nbr = filter(x -> 0<x<1, dd[i,:])
    @show nbr
    j=findall(x-> 0<x<1, dd[i,:])
    @show i, j, length(j)
    v_Npts[i]=length(j) 
end
max_nbrs = maximum(v_Npts)
max_Npts = findall(x -> x==max_nbrs, v_Npts)
select_i=[];
for i=1:length(max_Npts)
    @show i
    push!(select_i, sort([max_Npts[i]; findall(x-> 0<x<1, dd[max_Npts[i],:])]))
end
select_i = unique(select_i)

Nangles_max=100;
i_ctr=1;
#y=reduce(vcat,tmpy[select_i[i_ctr]]);
y=reduce(vcat,tmpy[select_i[i_ctr][1]]);

oco_soundings=[];
v_raz=[];
v_sza=[];
v_vza=[];
if Nangles_max==1
    tmp=[];
    push!(tmp, oco_soundings0[select_i[i_ctr][1]]);
    oco_soundings = tmp;
    v_vza = [vza[select_i[i_ctr][1]]];
    v_raz = [raz[select_i[i_ctr][1]]].+180.;
    v_sza = mean(sza[select_i[i_ctr][1]]);
else
    oco_soundings=oco_soundings0[select_i[i_ctr]];
    v_vza = vza[select_i[i_ctr]];
    v_raz = raz[select_i[i_ctr]].+180.;
    v_sza = mean(sza[select_i[i_ctr]]);
end
Nangles = length(v_vza)
parameters.sza = v_sza
parameters.vza = v_vza
parameters.vaz = v_raz

###################################################
# Produce black-body in wavenumber range
Tsolar = solar_transmission_from_file("/home/rjeyaram/vSmartMOM/src/SolarModel/solar.out")
Tsolar_interp = LinearInterpolation(Tsolar[:, 1], Tsolar[:, 2])

#Initial guess state vector 
x = FT[0.2377,     # 1. Legendre Lambertian Albedo, A-band
         -3.5,       # log(AOD)
         -0.00244, # 2. Legendre Lambertian Albedo, A-band
          0,       # 3. Legendre Lambertian Albedo, A-band
          0,       # 3. Legendre Lambertian Albedo, W-band
          0,       # 3. Legendre Lambertian Albedo, S-band
          0.2,     # 1. Legendre Lambertian Albedo, W-band
          0,       # 2. Legendre Lambertian Albedo, W-band
          0.2,     # 1. Legendre Lambertian Albedo, S-band
          0,       # 2. Legendre Lambertian Albedo, S-band
          0,       # H2O
          400e-6,  # CO2 1
          400e-6,  # CO2 2
          400e-6,  # CO2 3
          0,
          1.33,    # nᵣ
          0.01,    # nᵢ
        -1.69,     # lnr₀
          0.3,     # σ₀
          690.,    #p₀
          50.]     #σₚ


# Run FW model:
# @time runner(x);
#y = oco_sounding.SpectralMeasurement;

λ = oco_soundings[1].SpectralGrid*1e3 #nm
#λ = oco_soundings.SpectralGrid*1e3 #nm
h = 6.62607015e-34 # J.s
c = 2.99792458e8 # m/s
convfct = repeat(h*c*1.e16./λ, Nangles)
#ind = 92:885
#run_inversion(x,y)
#run_OE_inversion(x, y, x, Nangles, [length(indices[1]), length(indices[2]), length(indices[3])], convfct);
Fx = zeros(length(y));
runner!(Fx,x)
#Fx = Fx./convfct
#plot(y, label="Meas") #Ph sec^{-1} m^{-2} sr^{-1} um^{-1}
#plot!(Fx, label="Mod")# mW/m²-str-cm⁻¹ -> Ph sec^{-1} m^{-2} sr^{-1} um^{-1}

#plot!((y-Fx), label="Meas-mod")
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