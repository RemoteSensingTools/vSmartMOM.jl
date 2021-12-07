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
bands = (1,);
# Indices within that band:
indices = (1:1016,);
# Geo Index (footprint,sounding):
GeoInd = [5,5000];

# Get data for that sounding:
oco_sounding = InstrumentOperator.getMeasurement(oco, bands, indices, GeoInd);

###################################################
# Produce black-body in wavenumber range
T = 5777
λ_grid = reverse(1e4 ./ parameters.spec_bands[1]) #collect(757:0.01:777)
black_body = planck_spectrum_wl(T, λ_grid) * 2.1629e-05 * pi
black_body = SolarModel.watts_to_photons(λ_grid, black_body)

# Get solar transmittance spectrum 
#solar_transmission = solar_transmission_from_file("/home/rjeyaram/RadiativeTransfer/src/solar_merged_20160127_600_26316_100.out", parameters.spec_bands[1])
solar_transmission = solar_transmission_from_file("/home/rjeyaram/RadiativeTransfer/src/SolarModel/solar.out", parameters.spec_bands[1])
# Get outgoing solar radiation
sun_out = reverse(solar_transmission) .* black_body
###################################################


# Runner is used to set AD fields as duals
function runner!(y, x, parameters=parameters, oco_sounding= oco_sounding, sun_out=sun_out)

    # Set parameters fields as the dual numbers
    parameters.brdf = [vSmartMOM.LambertianSurfaceLegendre([x[1],x[3],x[4],x[5],x[6]])]

    parameters.scattering_params.rt_aerosols[1].τ_ref = exp(x[2]);
    parameters.scattering_params.rt_aerosols[1].p₀    = 800.0; #x[4]
   
    parameters.p   = oco_sounding.p_half
    parameters.q   = oco_sounding.q 
    parameters.T   = oco_sounding.T .+ 1.0
    parameters.sza = oco_sounding.sza
    parameters.vza = [oco_sounding.vza]

    model = model_from_parameters(parameters);
    @show sum(model.τ_rayl[1] )
    # Run the model to obtain reflectance matrix
    R = vSmartMOM.rt_run(model, i_band=1);

    RR = oco_sounding.mueller[1]*R[1,1,:] + oco_sounding.mueller[2]*R[1,2,:] + oco_sounding.mueller[2]*R[1,3,:]
    
    # Apply Earth reflectance matrix 
    earth_out = sun_out .* reverse(RR[:])
    # Re-interpolate I from ν_grid to new grid/resolution
    λ_grid = reverse(1e4 ./ parameters.spec_bands[1])
    interp_I = LinearInterpolation(λ_grid, earth_out);
    res = 0.001;
    wl = 757.5:res:771.0;
    I_wl = interp_I(wl/1000);

    # Convolve input spectrum with variable kernel
    I_conv = InstrumentOperator.conv_spectra(oco_sounding.ils[1], wl./1e3, I_wl)
    y[:] = I_conv[:]

end


# State vector
x = FT[0.2377, -3, -0.00244, 0, 0,0]#,


# Run FW model:
# @time runner(x);
I_conv = zeros(1016)
ind = 92:885
y = oco_sounding.SpectralMeasurement[ind]
for i=1:3
    dfdx = ForwardDiff.jacobian(runner!, I_conv, x);
    Fx = I_conv[ind];
    K = dfdx[ind,:];
    dx = K \ (y-Fx);
    x += dx;
    @show x;
end
runner!(I_conv,x)
Fx = I_conv[ind];

ν = oco_sounding.SpectralGrid[ind]*1e3
plot(ν,  y/1e20, label="Meas")
plot!(ν, Fx/1e20, label="Mod")
plot!(ν, (y-Fx)/1e20, label="Meas-mod")

end