using Revise
# using Plots
using RadiativeTransfer
using RadiativeTransfer.Absorption
using RadiativeTransfer.Scattering
using RadiativeTransfer.vSmartMOM
# using Distributions
# using BenchmarkTools
# using Test
# using CUDA

# Sets all the "specific" parameters
parameters = vSmartMOM.default_parameters();

# Generates all the derived attributes from above parameters
model = default_model(parameters);

@time R_GPU, T_GPU = vSmartMOM.rt_run(model);

# curr_parameters.

# Quadrature points for RTM
# Nquad, qp_μ, wt_μ = rt_set_streams(vSmartMOM.RadauQuad(), Ltrunc, FT(60.0), FT[0.0, 15.0, 30., 45., 60.])
# Nquad, qp_μ, wt_μ = rt_set_streams(vSmartMOM.GaussQuadFullSphere(), Ltrunc, FT(60.0), FT[0.0, 15.0, 30., 45., 60.])

# Aerosol particle distribution and properties
# μ            = [1.3]    # [0.3,2.0]       # Log mean radius
# σ            = [2.0]    # [2.0,1.8]       # Log stddev of radius
# r_max        = [30.0]   # [30.0,30.0]     # Maximum radius
# nquad_radius = [2500]   # [2500,2500]     # Number of quadrature points for integrating of size dist.
# nᵣ           = [1.3]    # [1.3, 1.66]     # Real part of refractive index
# nᵢ           = [0.00000001]  # [0.001,0.0003]  # Imag part of refractive index

# # Aerosol vertical distribution profiles
# p₀          = FT[90000.]  # [50000., 20000.] # Pressure peak [Pa]
# σp          = FT[5000.]   # [5000., 2000.]   # Pressure peak width [Pa]

# size_distribution = [LogNormal(log(μ[1]), log(σ[1]))] # [LogNormal(log(μ[1]), log(σ[1])), LogNormal(log(μ[2]), log(σ[2]))]

# # Create the aerosols (needs to be generalized through loops):
# aero1 = make_univariate_aerosol(size_distribution[1], r_max[1], nquad_radius[1], nᵣ[1], nᵢ[1])
# aero2 = make_univariate_aerosol(size_distribution[2], r_max[2], nquad_radius[2], nᵣ[2], nᵢ[2])

# Define some details, run aerosol optics
# model_NAI2_aero1 = make_mie_model(NAI2(), aero1, λ, polarization_type, truncation_type)
# aerosol_optics_NAI2_aero1 = compute_aerosol_optical_properties(model_NAI2_aero1, FT);

# # Truncate:
# aerosol_optics_trunc_aero1 = Scattering.truncate_phase(truncation_type, aerosol_optics_NAI2_aero1; reportFit=true)

# Define some details, run aerosol optics
# model_NAI2_aero2 = make_mie_model(NAI2(), aero2, λ, polarization_type, truncation_type)
# aerosol_optics_NAI2_aero2 = compute_aerosol_optical_properties(model_NAI2_aero2);
# Truncate:
# aerosol_optics_trunc_aero2 = Scattering.truncate_phase(truncation_type, aerosol_optics_NAI2_aero2)

# Rayleigh Greek
# GreekRayleigh = Scattering.get_greek_rayleigh(depol)


# vza = FT[60., 45., 30., 15., 0., 15., 30., 45., 60.]
# vaz = FT[180., 180., 180., 180., 0., 0., 0., 0., 0.]
# sza = FT(60.)

# obs_geom = ObsGeometry(FT(1000.0), sza, vza, vaz)

# Nquad, qp_μ, wt_μ = rt_set_streams(vSmartMOM.RadauQuad(), Ltrunc, obs_geom);
# Nquad, qp_μ, wt_μ = rt_set_streams(vSmartMOM.GaussQuadFullSphere(), Ltrunc, sza, vza);
# @show Nquad

# In[ ]:


# " Atmospheric Profiles, basics, needs to be refactore entirely"
# file = "/net/fluo/data1/ftp/XYZT_ESE156/Data/MERRA300.prod.assim.inst6_3d_ana_Nv.20150613.hdf.nc4" 
# # file = "MERRA300.prod.assim.inst6_3d_ana_Nv.20150613.hdf.nc4"  
# timeIndex = 2 # There is 00, 06, 12 and 18 in UTC, i.e. 6 hourly data stacked together

# # What latitude do we want? 
# myLat = 34.1377;
# myLon = -118.1253;

# # Read profile (and generate dry/wet VCDs per layer)
# profile_caltech_hr = vSmartMOM.read_atmos_profile(file, myLat, myLon, timeIndex);
# profile_caltech = vSmartMOM.reduce_profile(20, profile_caltech_hr);

# Compute layer optical thickness for Rayleigh (surface pressure in hPa) 
# τRayl =  vSmartMOM.getRayleighLayerOptProp(profile_caltech.psurf / 100, λ, depol, profile_caltech.vcd_dry);
# ϖRayl = ones(FT, length(τRayl));

# Compute Naer aerosol optical thickness profiles
# τAer_1 = vSmartMOM.getAerosolLayerOptProp(1.0, p₀[1], σp[1], profile_caltech.p_levels)
# # τAer_2 = vSmartMOM.getAerosolLayerOptProp(0.3, p₀[2], σp[2], profile_caltech.p_levels)

# # Can be done with arbitrary length later:
# τAer = FT(0.2) * τAer_1; # [τAer_1 τAer_2]
# @show sum(τAer)# , sum(τAer_2)
# ϖAer = FT[aerosol_optics_NAI2_aero1.ω̃]; # [aerosol_optics_NAI2_aero1.ω̃ aerosol_optics_NAI2_aero2.ω̃];
# # fᵗ   = FT[aerosol_optics_trunc_aero1.fᵗ]; # [aerosol_optics_trunc_aero1.fᵗ aerosol_optics_trunc_aero2.fᵗ];

# aerosol_optics = [aerosol_optics_trunc_aero1] # [aerosol_optics_trunc_aero1 aerosol_optics_trunc_aero2]
# Aer𝐙⁺⁺ = [aero1_Z⁺⁺] # [aero1_Z⁺⁺, aero2_Z⁺⁺];
# Aer𝐙⁻⁺ = [aero1_Z⁻⁺] # [aero1_Z⁻⁺, aero2_Z⁻⁺];

# maxM = 3

# grid = range(1e7 / 774, 1e7 / 757, length=1000);
# τ_abs = zeros(FT, length(grid), length(profile_caltech.p));
# hitran_data = read_hitran(artifact("O2"), iso=1)
# model = make_hitran_model(hitran_data, Voigt(), wing_cutoff=100, CEF=HumlicekWeidemann32SDErrorFunction(), architecture=GPU(), vmr=0.21)

# compute_absorption_profile!(τ_abs, model, grid, profile_caltech);



# @time R_GPU, T_GPU = vSmartMOM.rt_run(polarization_type, obs_geom, τRayl, ϖRayl, τAer, ϖAer, qp_μ, wt_μ, maxM, aerosol_optics, GreekRayleigh, τ_abs, Architectures.GPU());
# @time R_CPU, T_CPU = vSmartMOM.rt_run(polarization_type, obs_geom, τRayl, ϖRayl, τAer, ϖAer, qp_μ, wt_μ, maxM, aerosol_optics, GreekRayleigh, τ_abs, Architectures.CPU());

# @test R_CPU ≈ (R_GPU) 
# @test T_CPU ≈ (T_GPU) 