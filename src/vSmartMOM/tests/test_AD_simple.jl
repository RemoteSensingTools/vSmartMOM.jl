using Revise
using Plots
using Pkg
# Pkg.activate(".")
using RadiativeTransfer
using RadiativeTransfer.Absorption
using RadiativeTransfer.Scattering
using RadiativeTransfer.vSmartMOM
using RadiativeTransfer.SolarModel
using InstrumentOperator
using Interpolations
using Polynomials
using ForwardDiff 
using Distributions

## Atmospheric Radiative Transfer

# Load parameters from file
parameters = vSmartMOM.parameters_from_yaml("RadiativeTransfer/test/helper/O2Parameters.yaml")
FT = Float64

# Runner is used to set AD fields as duals
function runner!(y, x, parameters=parameters)

    parameters.scattering_params.rt_aerosols[1].τ_ref = x[1];
    parameters.scattering_params.rt_aerosols[1].p₀ = x[2];
    parameters.scattering_params.rt_aerosols[1].aerosol.nᵣ = x[3];
    parameters.scattering_params.rt_aerosols[1].aerosol.nᵢ = x[4];
    parameters.scattering_params.rt_aerosols[1].aerosol.size_distribution = LogNormal(log(x[5]), log(x[6]), check_args=false)

    model = model_from_parameters(parameters);
    R = vSmartMOM.rt_run(model, i_band=1);
    y[:] = R[1,1,:]

end

x = FT[0.1,90001.0,1.3, 1.0e-8, 1.3, 2.0]

# Run FW model:
# @time runner(x);
R = zeros(length(parameters.spec_bands[1]))
@time dfdx = ForwardDiff.jacobian(runner!, R, x);

## Solar Model 

# Produce black-body in wavenumber range
T = 5777
ν_grid = collect((1e7/777):0.015:(1e7/757))
black_body = planck_spectrum_wn(T, ν_grid)

# Get solar transmittance spectrum 
solar_transmission = solar_transmission_from_file("RadiativeTransfer/src/solar_merged_20160127_600_26316_100.out", ν_grid)

# Get outgoing solar radiation
sun_out = solar_transmission .* black_body

## Apply Earth reflectance matrix 

earth_out = sun_out .* R[:]

## Set up Instrument Operator

oco_file = "/home/cfranken/oco2_L1bScND_18688a_180105_B8100r_180206190633.h5"
ils_file = "/home/rjeyaram/RadiativeTransfer/src/vSmartMOM/ils_oco2.json"
ils_Δ, ils_in, dispersion = InstrumentOperator.read_ils_table(oco_file, ils_file);

# Define model grid:
res = 0.001

# Just consider the ILS within ± 0.35nm
grid_x = -0.35e-3:res*1e-3:0.35e-3
extended_dims = [5,1] # Footprint, band

# Re-interpolate I from ν_grid to new grid/resolution
interp_I = CubicSplineInterpolation(range(ν_grid[1], ν_grid[end], length=length(ν_grid)), earth_out);
wl = 757.5:res:771.0
I_wl = interp_I(1e7./wl)

# Pixels to be used
ind_out = collect(1:1016); 

# Eventual grid of OCO-2 for Band 1, FP 5:
dispPoly = Polynomial(view(dispersion, :, extended_dims...))
ν = Float32.(dispPoly.(0:1015))

# Prepare ILS table:
ils_pixel   = InstrumentOperator.prepare_ils_table(grid_x, ils_in, ils_Δ,extended_dims)
oco2_kernel = VariableKernelInstrument(ils_pixel, ν, ind_out)

# Convolve input spectrum with variable kernel
I_conv = conv_spectra(oco2_kernel, wl./1e3, I_wl)