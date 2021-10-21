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
FT = Float64

# oco2_L1bScND_18688a_180105_B8100r_180206190633

oco_file = "/net/fluo/data1/group/oco2/L1bSc/oco2_L1bScGL_15258a_170515_B10003r_200214061601.h5"
oco_met_file = "/net/fluo/data1/group/oco2/L2Met/oco2_L2MetGL_26777a_190715_B10003r_200429213029.h5"
ils_file = "test/prototyping/ils_oco2.json"

fp = 5
iOrbit = 5000
band = 1

sza_ = 32.4436
vza_ = [0.072]


function conv_spectra_local(m::VariableKernelInstrument, ν, spectrum; stride=1)
    # FT = eltype(m.ν_out)
    # Define grid where to perform convolution:
    
    # Padding at both sides required:
    off = ceil(Int, size(m.kernel, 1) / 2)
    ind = off:stride:(length(ν) - off)
    
    # knots where convolution will be applied to
    knots = view(ν, ind)
    te = LinearInterpolation(m.ν_out, Float32.(m.ind_out))
    spec_out = zeros(Real, length(knots));
    for i in eachindex(knots)
        # Simple first, nearest neighbor ILS
        ind_fraction = round(Int, te(knots[i]));
        kernel = view(m.kernel, :, ind_fraction)
        for j in eachindex(kernel)
            spec_out[i] += kernel[j] * spectrum[ind[i] + j] 
        end
    end
    # Change this later to only perform conv around output grid!
    fin = LinearInterpolation(ν[ind], spec_out; extrapolation_bc=Interpolations.Flat())
    return fin(m.ν_out)
end;

# Runner is used to set AD fields as duals
function runner!(y, x, parameters=parameters, oco_file=oco_file, 
                                              oco_met_file=oco_met_file, 
                                              ils_file=ils_file, fp=fp, iOrbit=iOrbit,
                                              sza_=sza_, vza_=vza_)

    # Set parameters fields as the dual numbers
    parameters.brdf = [vSmartMOM.LambertianSurfaceLegendre([x[1],x[3],x[4]])]

    parameters.scattering_params.rt_aerosols[1].τ_ref = exp(x[2]);
    parameters.scattering_params.rt_aerosols[1].p₀    = 70000.0; #x[4]
    # parameters.scattering_params.rt_aerosols[1].aerosol.size_distribution = LogNormal(log(x[3]), log(x[4]), check_args=false)

    # parameters.scattering_params.rt_aerosols[1].aerosol.nᵣ = x[5];
    # parameters.scattering_params.rt_aerosols[1].aerosol.nᵢ = x[6];

    # parameters.scattering_params.rt_aerosols[1].p₀ = x[7];
    # parameters.scattering_params.rt_aerosols[1].σp = x[8];
    ocoData = Dataset(oco_file)
    SoundingGeometry = ocoData.group["SoundingGeometry"]
    ϕ = SoundingGeometry["sounding_polarization_angle"][fp,iOrbit]
    m₁ = 0.5
    m₂ = cosd(2ϕ)/2
    m₃ = sind(2ϕ)/2
    # Stokes:

    # Set profiles properly
    met = Dataset(oco_met_file);
    T_met = met.group["Meteorology"]["temperature_profile_met"][:,fp,iOrbit];
    ak = met.group["MeteorologyDiagnostics"]["ak"][:,fp,iOrbit];
    bk = met.group["MeteorologyDiagnostics"]["bk"][:,fp,iOrbit];
    ak = [0.000000e+00, 4.804826e-02, 6.593752e+00, 1.313480e+01, 1.961311e+01, 2.609201e+01,
3.257081e+01, 3.898201e+01, 4.533901e+01, 5.169611e+01, 5.805321e+01, 6.436264e+01,
7.062198e+01, 7.883422e+01, 8.909992e+01, 9.936521e+01, 1.091817e+02, 1.189586e+02,
1.286959e+02, 1.429100e+02, 1.562600e+02, 1.696090e+02, 1.816190e+02, 1.930970e+02,
2.032590e+02, 2.121500e+02, 2.187760e+02, 2.238980e+02, 2.243630e+02, 2.168650e+02,
2.011920e+02, 1.769300e+02, 1.503930e+02, 1.278370e+02, 1.086630e+02, 9.236572e+01,
7.851231e+01, 6.660341e+01, 5.638791e+01, 4.764391e+01, 4.017541e+01, 3.381001e+01,
2.836781e+01, 2.373041e+01, 1.979160e+01, 1.645710e+01, 1.364340e+01, 1.127690e+01,
9.292942e+00, 7.619842e+00, 6.216801e+00, 5.046801e+00, 4.076571e+00, 3.276431e+00,
2.620211e+00, 2.084970e+00, 1.650790e+00, 1.300510e+00, 1.019440e+00, 7.951341e-01,
6.167791e-01, 4.758061e-01, 3.650411e-01, 2.785261e-01, 2.113490e-01, 1.594950e-01,
1.197030e-01, 8.934502e-02, 6.600001e-02, 4.758501e-02, 3.270000e-02, 2.000000e-02,
1.000000e-02][end:-1:1]

bk = [1.000000e+00, 9.849520e-01, 9.634060e-01, 9.418650e-01, 9.203870e-01, 8.989080e-01,
8.774290e-01, 8.560180e-01, 8.346609e-01, 8.133039e-01, 7.919469e-01, 7.706375e-01,
7.493782e-01, 7.211660e-01, 6.858999e-01, 6.506349e-01, 6.158184e-01, 5.810415e-01,
5.463042e-01, 4.945902e-01, 4.437402e-01, 3.928911e-01, 3.433811e-01, 2.944031e-01,
2.467411e-01, 2.003501e-01, 1.562241e-01, 1.136021e-01, 6.372006e-02, 2.801004e-02,
6.960025e-03, 8.175413e-09, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
0.000000e+00][end:-1:1]

    p_surf = 1.025*met.group["Meteorology"]["surface_pressure_met"][fp,iOrbit];
    p_half = (ak + bk * p_surf);
    #p_half = vcat(p_half, p_surf);
    q = met.group["Meteorology"]["specific_humidity_profile_met"][:,fp,iOrbit];
    parameters.p = p_half / 100
    parameters.q = q # zeros(size(q))
    parameters.T = T_met 
    parameters.T[end-30:end] .-= 10
    parameters.sza = sza_
    parameters.vza = vza_

    @show parameters.sza, parameters.vza

    # @show parameters.sza
    # plot!(parameters.p)
    # return 


    model = model_from_parameters(parameters);
    @show sum(model.τ_rayl[1] )
    # Run the model to obtain reflectance matrix
    R = vSmartMOM.rt_run(model, i_band=1);
    RR = m₁*R[1,1,:] + m₂*R[1,2,:] + m₃*R[1,3,:]
    @show size(R)
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

    # Apply Earth reflectance matrix 
    earth_out = sun_out .* reverse(RR[:])

    # Set up InstrumentOperator
    ils_Δ, ils_in, dispersion = InstrumentOperator.read_ils_table(oco_file, ils_file);

    # Define model grid:
    res = 0.001

    # Just consider the ILS within ± 0.35nm
    grid_x = -0.35e-3:res*1e-3:0.35e-3
    
    extended_dims = [fp,band] # Footprint, band

    # Re-interpolate I from ν_grid to new grid/resolution
    interp_I = LinearInterpolation(λ_grid, earth_out);
    wl = 757.5:res:771.0
    I_wl = interp_I(wl/1000)

    # Pixels to be used
    ind_out = collect(0:1015); 

    # Eventual grid of OCO-2 for Band 1, FP 5:
    dispPoly = Polynomial(view(dispersion, :, extended_dims...))
    ν = Float32.(dispPoly.(1:1016))

    # Prepare ILS table:
    ils_pixel   = InstrumentOperator.prepare_ils_table(grid_x, ils_in, ils_Δ,extended_dims)
    oco2_kernel = VariableKernelInstrument(ils_pixel, ν, ind_out)

    # Convolve input spectrum with variable kernel
    I_conv = conv_spectra_local(oco2_kernel, wl./1e3, I_wl)

    y[:] = I_conv[:]

end

directory = "/net/fluo/data1/group/oco2/L1bSc/"
files_list = filter(x->endswith(x, ".h5") && startswith(x, "oco2_L1bScND_26"), readdir(directory))

met_files_list = filter(x->endswith(x, ".h5"), readdir("/net/fluo/data1/group/oco2/L2Met/"))

I_convs_all = zeros(1016, length(files_list))
oco2_Abands_all = zeros(1016, length(files_list))

albedos = zeros(length(files_list))
cloud_levels = zeros(length(files_list))

optimal_x = zeros(8, length(files_list))

#for i in 1:length(files_list)
i = 1
ind = 92:885
@show i

oco_file = files_list[i]

    # try
ocoData = Dataset(directory * oco_file)
oco_met_file = "/net/fluo/data1/group/oco2/L2Met/" * filter(x->contains(x, oco_file[14:19]), met_files_list)[1]
oco_file = directory * oco_file
met = Dataset(oco_met_file)
cloud_levels[i] = met.group["Meteorology"]["cloud_liquid_water_path_met"][fp,iOrbit]


## Getting an OCO spectrum for fun to fit:
o2AbandSpectra = ocoData.group["SoundingMeasurements"]["radiance_o2"]
SoundingGeometry = ocoData.group["SoundingGeometry"]
ϕ = SoundingGeometry["sounding_polarization_angle"][fp,iOrbit]
m₁ = 0.5
m₂ = cosd(2ϕ)/2
m₃ = sind(2ϕ)/2
vza = SoundingGeometry["sounding_zenith"]
sza = SoundingGeometry["sounding_solar_zenith"]
lat = SoundingGeometry["sounding_latitude"]
lon = SoundingGeometry["sounding_longitude"]
o2AbandSpectra = ocoData.group["SoundingMeasurements"]["radiance_o2"]
oco2_Aband = o2AbandSpectra[:,fp,iOrbit]

println("$(round(Float64(lat[fp,iOrbit]), digits=2)) $(round(Float64(lon[fp,iOrbit]), digits=4))")
# continue
sza_ = sza[fp,iOrbit]
vza_ = [vza[fp,iOrbit]]

# State vector
x = FT[0.2377,
        -3, -0.00244, 0]#,


# Run FW model:
# @time runner(x);
I_conv = zeros(1016)
y = oco2_Aband[ind];
for i=1:3
    dfdx = ForwardDiff.jacobian(runner!, I_conv, x);
    Fx = I_conv[ind];
    K = dfdx[ind,:];
    dx = K \ (y-Fx)
    x += dx
end
runner!(I_conv,x)
Fx = I_conv[ind];

        
# plot fits
ils_Δ, ils_in, dispersion = InstrumentOperator.read_ils_table(oco_file, ils_file);
extended_dims = [fp,band]; # Footprint, band
# Eventual grid of OCO-2 for Band 1, FP 5:
dispPoly = Polynomial(view(dispersion, :, extended_dims...));
ν = Float32.(dispPoly.(1:1016));
plot(ν[ind]*1e3,  y/1e20, label="Meas")
plot!(ν[ind]*1e3, Fx/1e20, label="Mod")
plot!(ν[ind]*1e3, (y-Fx)/1e20, label="Meas-mod")