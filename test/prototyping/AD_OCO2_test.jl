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

## Atmospheric Radiative Transfer

# Load parameters from file
parameters = parameters_from_yaml("test/test_parameters/3BandParameters.yaml")
#parameters.architecture = CPU()
FT = Float64

# Load OCO Data: 
# File names:
L1File   = "/net/fluo/data1/group/oco2/L1bSc/oco2_L1bScND_26780a_190715_B10003r_200429212407.h5"
metFile  = "/net/fluo/data1/group/oco2/L2Met/oco2_L2MetND_26780a_190715_B10003r_200429212406.h5"
dictFile = "/home/cfranken/code/gitHub/InstrumentOperator.jl/json/oco2.yaml"
# Load L1 file (could just use filenames here as well)
oco = InstrumentOperator.load_L1(dictFile,L1File, metFile);


# Pick some bands as tuple (or just one)
bands = (1,2,3);
#bands = (1,3);
# Indices within that band:
indices = (92:885,114:845,50:916);
#indices = (92:885,50:916);
# Geo Index (footprint,sounding):
GeoInd = [5,5000];

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
    parameters.brdf =  [CoreRT.LambertianSurfaceLegendre([x[1],x[3],x[4],x[16],x[17] ]),
                        CoreRT.LambertianSurfaceLegendre([x[7],x[8],x[5]]),
                        CoreRT.LambertianSurfaceLegendre([x[9],x[10],x[6]])];

    parameters.scattering_params.rt_aerosols[1].τ_ref = 0.1+ 0*(x[2]);
    parameters.scattering_params.rt_aerosols[1].p₀    = 600.0; #x[4]
   
    parameters.p   = oco_sounding.p_half
    parameters.q   = oco_sounding.q 
    parameters.T   = oco_sounding.T# .+ 1.0 #.+ x[15]
    parameters.sza = oco_sounding.sza
    parameters.vza = [oco_sounding.vza]
    parameters.absorption_params.vmr["H2O"] = [parameters.q[1:65]*x[11] * 1.8; 
                                               parameters.q[66:end]*x[15] * 1.8];
    a1 = zeros(25) .+ x[12]
    a2 = zeros(25) .+ x[13]
    a3 = zeros(22) .+ x[14]
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
        
        #@show findall(sun_out .< 0)
        # Re-interpolate I from ν_grid to new grid/resolution
        λ_grid = reverse(1e4 ./ parameters.spec_bands[i])
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


# State vector
x = FT[0.2377, 
        -2, 
        -0.00244, 
        0, 
        0,
        0, 
        0.44, 
        0, 
        0.4, 
        0,
        1,
        410e-6,
        390e-6,
        400e-6, 
        1.0, 0.0, 0.0]#,


# Run FW model:
# @time runner(x);
y = oco_sounding.SpectralMeasurement;
Fx = zeros(length(y));
#ind = 92:885

for i=1:3
    K = ForwardDiff.jacobian(runner!, Fx, x);
    dx = K \ (y-Fx);
    x += dx;
    @show x;
end
runner!(Fx,x)

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
