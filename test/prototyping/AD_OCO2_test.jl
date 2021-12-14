using Revise
using Plots
using Pkg
# Pkg.activate(".")
using RadiativeTransfer
using RadiativeTransfer.Architectures
using RadiativeTransfer.Absorption
using RadiativeTransfer.Scattering
using RadiativeTransfer.vSmartMOM
using RadiativeTransfer.SolarModel
using InstrumentOperator
using Interpolations
using Polynomials
using ForwardDiff 
using Distributions
using NCDatasets

## Atmospheric Radiative Transfer

# Load parameters from file
parameters = vSmartMOM.parameters_from_yaml("test/test_parameters/O2Parameters.yaml")
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
# Indices within that band:
indices = (92:885,50:950,50:950);
# Geo Index (footprint,sounding):
GeoInd = [5,5000];

# Get data for that sounding:
oco_sounding = InstrumentOperator.getMeasurement(oco, bands, indices, GeoInd);

###################################################
# Produce black-body in wavenumber range


Tsolar = solar_transmission_from_file("/home/rjeyaram/RadiativeTransfer/src/SolarModel/solar.out")
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
    parameters.brdf = [vSmartMOM.LambertianSurfaceLegendre([x[1],x[3],x[4],x[5],x[6]]),vSmartMOM.LambertianSurfaceLegendre([x[7],x[8]]),vSmartMOM.LambertianSurfaceLegendre([x[9],x[10]])]

    parameters.scattering_params.rt_aerosols[1].τ_ref = exp(x[2]);
    parameters.scattering_params.rt_aerosols[1].p₀    = 800.0; #x[4]
   
    parameters.p   = oco_sounding.p_half
    parameters.q   = oco_sounding.q 
    parameters.T   = oco_sounding.T .+ 1.0
    parameters.sza = oco_sounding.sza
    parameters.vza = [oco_sounding.vza]
    parameters.absorption_params.vmr["H2O"] = parameters.q*x[11]
    parameters.absorption_params.vmr["CO2"] = x[12]
    model = model_from_parameters(parameters);

    @show sum(model.τ_rayl[1] )
    for i = 1:length(oco_sounding.BandID)
        println("Modelling band $(i)")
        # Run the model to obtain reflectance matrix
        R = vSmartMOM.rt_run(model, i_band=i);
        RR = oco_sounding.mueller[1]*R[1,1,:] + oco_sounding.mueller[2]*R[1,2,:] + oco_sounding.mueller[2]*R[1,3,:]
        
        # Get sun:
        @time sun_out = getSolar(parameters.spec_bands[i],Tsolar)
        # Apply Earth reflectance matrix 
        earth_out = sun_out .* reverse(RR[:])
        
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
        y[oco_sounding.BandID[i]] = I_conv
    end

end


# State vector
x = FT[0.2377, -3, -0.00244, 0, 0,0, 0.4, 0, 0.4, 0,1,400e-6]#,


# Run FW model:
# @time runner(x);
y = oco_sounding.SpectralMeasurement;
Fx = zeros(length(y));
#ind = 92:885

for i=1:5
    K = ForwardDiff.jacobian(runner!, Fx, x);
    dx = K \ (y-Fx);
    x += dx;
    @show x;
end
runner!(Fx,x)

ν = oco_sounding.SpectralGrid*1e3
plot(ν,  y/1e20, label="Meas")
plot!(ν, Fx/1e20, label="Mod")
plot!(ν, (y-Fx)/1e20, label="Meas-mod")

end