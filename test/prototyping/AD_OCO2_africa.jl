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
using LinearAlgebra

## Atmospheric Radiative Transfer

# Load parameters from file
parameters = parameters_from_yaml("test/test_parameters/3BandParameters.yaml")
#parameters.architecture = CPU()
FT = Float64

# Load OCO Data: 
# File names:

# Africa
L1File   = "/net/fluo/data1/group/oco2/L1bSc/oco2_L1bScND_31182a_200512_B10206r_210427151019.h5"
metFile  = "/net/fluo/data1/group/oco2/L2Met/oco2_L2MetND_31182a_200512_B10206r_210425233338.h5"
dictFile = "/home/cfranken/code/gitHub/InstrumentOperator.jl/json/oco2.yaml"
# Load L1 file (could just use filenames here as well)
oco = InstrumentOperator.load_L1(dictFile,L1File, metFile);


# Pick some bands as tuple (or just one)
bands = (1,2,3);
#bands = (1,3);
# Indices within that band:
#ind1 = [collect(92:68);  collect(690:870); collect(875:885)]
indices = (92:682,114:845,50:916);

#indices = (92:885,114:845,50:916);
# Geo Index (footprint,sounding):
GeoInd = [3, 3391];

# Get data for that sounding:
oco_sounding = InstrumentOperator.getMeasurement(oco, bands, indices, GeoInd);




###################################################
# Produce black-body in wavenumber range


Tsolar = solar_transmission_from_file("/home/rjeyaram/vSmartMOM/src/SolarModel/solar.out")
Tsolar_interp = LinearInterpolation(Tsolar[:, 1], Tsolar[:, 2])

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
function runner!(y, x, parameters=parameters, oco_sounding= oco_sounding, Tsolar = Tsolar_interp)

    # Set parameters fields as the dual numbers
    parameters.brdf =  [CoreRT.LambertianSurfaceLegendre([x[1],x[2],x[3]]),
                        CoreRT.LambertianSurfaceLegendre([x[4],x[5],x[6]]),
                        CoreRT.LambertianSurfaceLegendre([x[7],x[8],x[9]])];

    parameters.scattering_params.rt_aerosols[1].τ_ref = exp(x[10]);
    parameters.scattering_params.rt_aerosols[1].p₀    = x[11]
    size_distribution = LogNormal(x[12], 0.3)
    parameters.scattering_params.rt_aerosols[1].aerosol.size_distribution =  size_distribution
    parameters.scattering_params.rt_aerosols[1].aerosol.nᵣ = x[13]
    parameters.p   = oco_sounding.p_half
    parameters.q   = oco_sounding.q 
    parameters.T   = oco_sounding.T# .+ 1.0 #.+ x[15]
    parameters.sza = oco_sounding.sza
    parameters.vza = [oco_sounding.vza]
    parameters.absorption_params.vmr["O2"] = x[18]
    parameters.absorption_params.vmr["H2O"] = parameters.q*x[14] * 1.8 #[parameters.q[1:65]*x[11] * 1.8; 
                                               #parameters.q[66:end]*x[15] * 1.8];
    a1 = zeros(25) .+ x[15] * 1e-6
    a2 = zeros(25) .+ x[16] * 1e-6
    a3 = zeros(22) .+ x[17] * 1e-6
    parameters.absorption_params.vmr["CO2"] = [a1; a2; a3]
    model = model_from_parameters(parameters);
    
    for i = 1:length(oco_sounding.BandID)
        println("Modelling band $(i)")
        # Run the model to obtain reflectance matrix
        #R = rt_run(model, i_band=i)[1];
        R = CoreRT.rt_run_test(CoreRT.noRS(), model, i)[1]
        RR = oco_sounding.mueller[1]*R[1,1,:] + oco_sounding.mueller[2]*R[1,2,:] + oco_sounding.mueller[3]*R[1,3,:]
        
        # Get sun:
        @time sun_out = getSolar(parameters.spec_bands[i],Tsolar)
        # Apply Earth reflectance matrix 
        earth_out = sun_out .* reverse(RR[:])
        #ii = findall(earth_out .> 1e50)
        #@show ii
        #@show ii
        #if length(ii==1)
        #    earth_out[ii] =  earth_out[ii+1]
        #end 
        #f = 0.9999908003445712
        #@show findall(sun_out .< 0)
        # Re-interpolate I from ν_grid to new grid/resolution

        # Apply doppler here:
        λ_grid = reverse(1e4 ./ parameters.spec_bands[i]) .* oco_sounding.doppler
        @time interp_I = LinearInterpolation(λ_grid, earth_out);
        res = 0.001e-3;
        off = 0.5e-3
        wl = oco_sounding.ils[i].ν_out[1]-off:res:oco_sounding.ils[i].ν_out[end]+off;
        @show wl[1],wl[end], λ_grid[1],λ_grid[end]
        I_wl = interp_I(wl);

        # Convolve input spectrum with variable kernel
        @time I_conv = InstrumentOperator.conv_spectra(oco_sounding.ils[i], wl, I_wl)
        y[oco_sounding.BandID[i]] = I_conv/π
    end
    @show xco2 = 1e6*(model.profile.vmr["CO2"]' * model.profile.vcd_dry)/sum(model.profile.vcd_dry)

end


# Prior State vector
xₐ = FT[  0.376,   # Albedo band 1, degree 1
        0,         # Albedo band 1, degree 2       
        0,         # Albedo band 1, degree 3       
        0.251,#0.139,     # Albedo band 2, degree 1
        0,         # Albedo band 2, degree 2
        0,         # Albedo band 2, degree 3
        0.1,#0.0512,    # Albedo band 3, degree 1
        0,         # Albedo band 3, degree 2
        0.00,      # Albedo band 3, degree 3
        -1.0,      # log AOD
        800.0,     # Aerosol peak height (hPa)
        -1.0,       # log of size distribution mean (1µm here)
        1.3,       # Real refractive index of aerosol
        1,         # H2O scaling factor for entire profile
        400,    # CO2 Top layers
        400,    # CO2 Middle layers
        400,    # CO2 Bottom layers
        0.2095  # O2 VMR
        ]#,

# Measurement covariance matrix:
Sₑ = Diagonal(oco_sounding.SpectralNoise.^2);
# Prior covariance matrix
n = length(xₐ)
Sₐ = zeros(n,n)
# For albedos
Sₐ[1,1] = 1.2^2
Sₐ[4,4] = 1.2^2
Sₐ[7,7] = 1.2^2
Sₐ[2,2] = 0.2^2
Sₐ[5,5] = 0.2^2
Sₐ[8,8] = 0.2^2
Sₐ[3,3] = 0.2^2
Sₐ[6,6] = 0.2^2
Sₐ[9,9] = 0.2^2
# Aerosols:
Sₐ[10,10] = 0.5^2  # AOD, 2.5 orders of magnitude
Sₐ[11,11] = 10.0^2
Sₐ[12,12] = 0.05^2
Sₐ[13,13] = 0.05^2
Sₐ[14,14] = 0.5^2
Sₐ[15,15] = 100^2
Sₐ[16,16] = 100^2
Sₐ[17,17] = 100^2
Sₐ[18,18] = 0.002^2

# Run FW model:
# @time runner(x);
y = oco_sounding.SpectralMeasurement;
Fx = zeros(length(y));
#ind = 92:885
x = xₐ
xᵢ = []
push!(xᵢ, xₐ)
# LM parameter only for aerosols: 
γ = zeros(n,n);
γ[10,10] = 50;
γ[11,11] = 50;
γ[12,12] = 50;
γ[13,13] = 50;
γ[end,end] = 50;
for i=1:15
    
    K = ForwardDiff.jacobian(runner!, Fx, xᵢ[end]);
    res_1 = (y-Fx)' * inv(Sₑ) * (y-Fx)
    G = inv(K'inv(Sₑ)K + (1.0 .+ γ).*inv(Sₐ))K'inv(Sₑ)
    δy = (y-Fx + K*(xᵢ[end]-xₐ))
    dx = G * δy
    runner!(Fx,xₐ + dx)
    res_2 = (y-Fx)' * inv(Sₑ) * (y-Fx)
    if res_2 < res_1
        # Reduce γ
        @show "Reduce γ", i
        γ *= 0.5
        push!(xᵢ, xₐ + dx);
    else
        #Increase γ
        @show "Increase γ", i
        γ *= 2
    end
    #@show x;
end
Fx2 = similar
runner!(Fx,xᵢ[end])

ν = oco_sounding.SpectralGrid*1e3
plot(y/1e20, label="Meas")
plot!(Fx/1e20, label="Mod")
plot!((y-Fx)/1e20, label="Meas-mod")

a = Dataset("/net/fluo/data1/group/oco2/oco2_L2DiaND_26780a_190715_B10004r_200618191413.h5")
g = a.group["RetrievalGeometry"]
s = a.group["SpectralParameters"]
i = findall(abs.(g["retrieval_latitude"][:].-41.736668).<0.0001)
fp_meas = s["measured_radiance"][:,i]
fp_mod  = s["modeled_radiance"][:,i]
wl = s["wavelength"][:,i][:,1]*1e3
aod_L2 = a.group["AerosolResults"]["aerosol_aod"][i]
RR = a.group["RetrievalResults"]
SV = a.group["RetrievedStateVector"]
sv_L2 = SV["state_vector_result"][:,i]
RR["xco2"][i] * 1e6


l = @layout [a b c; d e f]
ii = oco_sounding.BandID[1]
p1 = plot(ν[ii],y[ii]/1e20, label=nothing)
plot!(ν[ii],Fx[ii]/1e20, label=nothing)
p4 = plot(ν[ii],(Fx[ii] - y[ii])/maximum(Fx[ii]) * 100, label=nothing )
ylims!(-1,1)  
ii = oco_sounding.BandID[2]
p2 = plot(ν[ii],y[ii]/1e20 , label=nothing )
plot!(ν[ii],Fx[ii]/1e20, label=nothing)
p5 = plot(ν[ii],(Fx[ii] - y[ii])/maximum(Fx[ii]) * 100, label=nothing )
ylims!(-1,1) 
ii = oco_sounding.BandID[3]
p3 = plot(ν[ii],y[ii]/1e20, label="Meas/1e20"  )

plot!(ν[ii],Fx[ii]/1e20, label="Mod/1e20")
p6 = plot(ν[ii],(Fx[ii] - y[ii])/maximum(Fx[ii]) * 100, label="(Fx-y)/max(Fx)*100")
ylims!(-1,1)  
plot(p1, p2, p3, p4,p5,p6,layout = l, size = (1000, 700))
savefig("co2_fit_vSmartMOM.pdf")

ii = findall((wl .< 800) .&& (wl.>0))
p1 = plot(wl[ii],fp_meas[ii]/1e20, label=nothing)
plot!(wl[ii],fp_mod[ii]/1e20, label=nothing)
p4 = plot(wl[ii],(fp_mod[ii] - fp_meas[ii])/maximum(fp_mod[ii]) * 100, label=nothing )
ylims!(-1,1)  

ii = findall((wl .> 800) .&& (wl.<1800))
p2 = plot(wl[ii],fp_meas[ii]/1e20 , label=nothing )
plot!(wl[ii],fp_mod[ii]/1e20, label=nothing)
p5 = plot(wl[ii],(fp_mod[ii] - fp_meas[ii])/maximum(fp_mod[ii]) * 100, label=nothing )
ylims!(-1,1) 
ii = findall(((wl.>1900)))
p3 = plot(wl[ii],fp_meas[ii]/1e20, label="Meas/1e20"  )

plot!(wl[ii],fp_mod[ii]/1e20, label="Mod/1e20")
p6 = plot(wl[ii],(fp_mod[ii] - fp_meas[ii])/maximum(fp_mod[ii]) * 100, label="(Fx-y)/max(Fx)*100")
ylims!(-1,1)  
plot(p1, p2, p3, p4,p5,p6,layout = l, size = (1000, 700))
savefig("co2_fit_JPL.pdf")

# EOF test
eof_d = Dataset("/home/cfranken/b10-baseline4_eofsC_oceanG_20141127-20190304_alt5_falt1_sorted_L2.h5")
gg = eof_d.group["Instrument"].group["EmpiricalOrthogonalFunction"].group["Glint"]
plot(gg["EOF_1_waveform_1"])

l = @layout [a b c]
ii = oco_sounding.BandID[1]

p1 = plot(ν[ii],(Fx[ii] - y[ii])/maximum(Fx[ii]) * 100, label=nothing )
plot!(ν[ii],5gg["EOF_1_waveform_1"][indices[1],5], label=nothing)
plot!(ν[ii],5gg["EOF_2_waveform_1"][indices[1],5], label=nothing)

ii = oco_sounding.BandID[2]
p2 = plot(ν[ii],(Fx[ii] - y[ii])/maximum(Fx[ii]) * 100, label=nothing )
plot!(ν[ii],5gg["EOF_1_waveform_2"][indices[2],5], label=nothing)
plot!(ν[ii],5gg["EOF_2_waveform_2"][indices[2],5], label=nothing)

ii = oco_sounding.BandID[3]
p3 = plot(ν[ii],(Fx[ii] - y[ii])/maximum(Fx[ii]) * 100, label="(Fx-y)/max(Fx)*100")
plot!(ν[ii],5gg["EOF_1_waveform_3"][indices[3],5], label="EOF1")
plot!(ν[ii],5gg["EOF_2_waveform_3"][indices[3],5], label="EOF2")
plot(p1, p2, p3,layout = l, size = (1000, 500))
savefig("co2Residual_fit_vSmartMOM.pdf")


# Run O2 band again with Raman!
