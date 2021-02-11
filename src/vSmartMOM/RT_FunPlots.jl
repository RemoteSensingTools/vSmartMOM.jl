using Revise
using Plots
using RadiativeTransfer
using RadiativeTransfer.Absorption
using RadiativeTransfer.PhaseFunction
using RadiativeTransfer.RTM
using Distributions
using BenchmarkTools
using Test
using CUDA
using StatsPlots

device!(3)

FT = Float32
"Generate aerosol optical properties"

# Wavelength (just one for now)
λ = FT(0.770)       # Incident wavelength
depol = FT(0.0)
# Truncation 
Ltrunc = 20             # Truncation  
truncation_type   = PhaseFunction.δBGE{Float32}(Ltrunc, 2.0)

# polarization_type
polarization_type = Stokes_IQU{FT}()

# Quadrature points for RTM
# Nquad, qp_μ, wt_μ = rt_set_streams(RTM.RadauQuad(), Ltrunc, FT(60.0), FT[0.0, 15.0, 30., 45., 60.])
# Nquad, qp_μ, wt_μ = rt_set_streams(RTM.GaussQuadFullSphere(), Ltrunc, FT(60.0), FT[0.0, 15.0, 30., 45., 60.])

# Aerosol particle distribution and properties
μ            = [1.3]    # [0.3,2.0]       # Log mean radius
σ            = [2.0]    # [2.0,1.8]       # Log stddev of radius
r_max        = [30.0]   # [30.0,30.0]     # Maximum radius
nquad_radius = [2500]   # [2500,2500]     # Number of quadrature points for integrating of size dist.
nᵣ           = [1.3]    # [1.3, 1.66]     # Real part of refractive index
nᵢ           = [0.00000001]  # [0.001,0.0003]  # Imag part of refractive index

# Aerosol vertical distribution profiles
p₀          = FT[90000.]  # [50000., 20000.] # Pressure peak [Pa]
σp          = FT[5000.]   # [5000., 2000.]   # Pressure peak width [Pa]

size_distribution = [LogNormal(log(μ[1]), log(σ[1]))] # [LogNormal(log(μ[1]), log(σ[1])), LogNormal(log(μ[2]), log(σ[2]))]

# Create the aerosols (needs to be generalized through loops):
aero1 = make_univariate_aerosol(size_distribution[1], r_max[1], nquad_radius[1], nᵣ[1], nᵢ[1])
# aero2 = make_univariate_aerosol(size_distribution[2], r_max[2], nquad_radius[2], nᵣ[2], nᵢ[2])

# Define some details, run aerosol optics
model_NAI2_aero1 = make_mie_model(NAI2(), aero1, λ, polarization_type, truncation_type)
aerosol_optics_NAI2_aero1 = compute_aerosol_optical_properties(model_NAI2_aero1, FT);

# Truncate:
aerosol_optics_trunc_aero1 = PhaseFunction.truncate_phase(truncation_type, aerosol_optics_NAI2_aero1; reportFit=true)

# Define some details, run aerosol optics
# model_NAI2_aero2 = make_mie_model(NAI2(), aero2, λ, polarization_type, truncation_type)
# aerosol_optics_NAI2_aero2 = compute_aerosol_optical_properties(model_NAI2_aero2);
# Truncate:
# aerosol_optics_trunc_aero2 = PhaseFunction.truncate_phase(truncation_type, aerosol_optics_NAI2_aero2)

# Rayleigh Greek
GreekRayleigh = PhaseFunction.get_greek_rayleigh(depol)


# In[ ]:


vza = FT[60., 45., 30., 15., 0., 15., 30., 45., 60.]
vaz = FT[180., 180., 180., 180., 0., 0., 0., 0., 0.]
sza = FT(60.)

obs_geom = ObsGeometry(FT(1000.0), sza, vza, vaz)

Nquad, qp_μ, wt_μ = rt_set_streams(RTM.RadauQuad(), Ltrunc, obs_geom);
# Nquad, qp_μ, wt_μ = rt_set_streams(RTM.GaussQuadFullSphere(), Ltrunc, sza, vza);
@show Nquad

# In[ ]:


" Atmospheric Profiles, basics, needs to be refactore entirely"
file = "/net/fluo/data1/ftp/XYZT_ESE156/Data/MERRA300.prod.assim.inst6_3d_ana_Nv.20150613.hdf.nc4" 
# file = "MERRA300.prod.assim.inst6_3d_ana_Nv.20150613.hdf.nc4"  
timeIndex = 2 # There is 00, 06, 12 and 18 in UTC, i.e. 6 hourly data stacked together

# What latitude do we want? 
myLat = 34.1377;
myLon = -118.1253;

# Read profile (and generate dry/wet VCDs per layer)
profile_caltech_hr = RTM.read_atmos_profile(file, myLat, myLon, timeIndex);
profile_caltech = RTM.reduce_profile(20, profile_caltech_hr);
# Compute layer optical thickness for Rayleigh (surface pressure in hPa) 
τRayl =  RTM.getRayleighLayerOptProp(profile_caltech.psurf / 100, λ, depol, profile_caltech.vcd_dry);
ϖRayl = ones(FT, length(τRayl));

# Compute Naer aerosol optical thickness profiles
τAer_1 = RTM.getAerosolLayerOptProp(1.0, p₀[1], σp[1], profile_caltech.p_levels)
# τAer_2 = RTM.getAerosolLayerOptProp(0.3, p₀[2], σp[2], profile_caltech.p_levels)

# Can be done with arbitrary length later:
τAer = FT(0.2) * τAer_1; # [τAer_1 τAer_2]
@show sum(τAer)# , sum(τAer_2)
ϖAer = FT[aerosol_optics_NAI2_aero1.ω̃]; # [aerosol_optics_NAI2_aero1.ω̃ aerosol_optics_NAI2_aero2.ω̃];
fᵗ   = FT[aerosol_optics_trunc_aero1.fᵗ]; # [aerosol_optics_trunc_aero1.fᵗ aerosol_optics_trunc_aero2.fᵗ];

aerosol_optics = [aerosol_optics_trunc_aero1] # [aerosol_optics_trunc_aero1 aerosol_optics_trunc_aero2]
# Aer𝐙⁺⁺ = [aero1_Z⁺⁺] # [aero1_Z⁺⁺, aero2_Z⁺⁺];
# Aer𝐙⁻⁺ = [aero1_Z⁻⁺] # [aero1_Z⁻⁺, aero2_Z⁻⁺];

maxM = 3

#################################################################
# PLOTTING RUNTIME VS. N_SPEC (CPU vs GPU)
#################################################################

elapsed_absorption_time = []
elapsed_RT_time_cpu = []
elapsed_RT_time_gpu_32 = []

grid_sizes = [100, 500, 1000, 2500, 5000, 7500, 10000]


for grid_size in grid_sizes

    @show grid_size

    grid = range(1e7 / 780, 1e7 / 755, length=grid_size)
    τ_abs = zeros(FT, length(grid), length(profile_caltech.p));

    hitran_data = read_hitran(artifact("O2"), iso=1)
    model = make_hitran_model(hitran_data, Voigt(), wing_cutoff=100, CEF=HumlicekWeidemann32SDErrorFunction(), architecture=Absorption.GPU(), vmr=0.21)

    compute_absorption_profile!(τ_abs, model, grid, profile_caltech);
    # push!(elapsed_absorption_time, @elapsed compute_absorption_profile!(grid, τ_abs, profile_caltech))
    
    push!(elapsed_RT_time_gpu_32, @elapsed R_GPU, T_GPU = RTM.rt_run(polarization_type, obs_geom, τRayl, ϖRayl, τAer, ϖAer, fᵗ, qp_μ, wt_μ, maxM, aerosol_optics, GreekRayleigh, τ_abs, RadiativeTransfer.Architectures.GPU()))

    # push!(elapsed_RT_time_cpu, @elapsed RTM.run_RTM(polarization_type, sza, vza, vaz, τRayl, ϖRayl, τAer, ϖAer, fᵗ, qp_μ, wt_μ, maxM, aerosol_optics, GreekRayleigh, τ_abs, RadiativeTransfer.Architectures.CPU()))

end

plot(grid_sizes, elapsed_absorption_time)
plot(grid_sizes, elapsed_RT_time_gpu)
plot(grid_sizes, elapsed_RT_time_cpu)

groupedbar([elapsed_RT_time_gpu_32 elapsed_RT_time_gpu_64], bar_position = :dodge, bar_width=0.7, xticks=(1:12, grid_sizes), legend=:topleft, label=["32" "64"], ylabel="Time (s)", xlabel="N_Spec", title="Runtime with Float32 vs Float64")

# elapsed_absorption_time = [0.282189253, 0.283802194, 0.307112032, 0.3477742, 0.38201981, 0.372605842, 0.468313961]
# elapsed_RT_time_gpu = [1.368995249, 4.598695304, 30.890917859, 30.341784914, 71.295454828, 107.297501462, 143.935854327]
elapsed_RT_time_gpu_64 = 
elapsed_RT_time_gpu_32 = [1.208653814, 3.607311384 , 29.202259484, 25.388859926, 40.688672364, 74.440302864, 97.265656638]
# elapsed_RT_time_cpu = [14.577631018, 72.8522157, 138.689995059, 349.936594999, 683.309987181, 1013.949079029, 1645.13642589]


#################################################################
# PLOTTING R/T VS. SZA (GPU)
#################################################################

# szas = collect(-90.0:1.0:90.0)
# Rs = zeros(9,3,10000,length(szas))
# Ts = zeros(9,3,10000,length(szas))

# grid = range(1e7 / 780, 1e7 / 755, length=5000)
# τ_abs = zeros(length(grid), length(profile_caltech.p))
# compute_absorption_profile!(grid, τ_abs, profile_caltech)

# for sza_idx in 1:length(szas)

#     @show sza_idx

#     model_NAI2_aero1 = make_mie_model(NAI2(), aero1, λ, polarization_type, truncation_type)
#     aerosol_optics_NAI2_aero1 = compute_aerosol_optical_properties(model_NAI2_aero1);

#     # Truncate:
#     aerosol_optics_trunc_aero1 = PhaseFunction.truncate_phase(truncation_type, aerosol_optics_NAI2_aero1; reportFit=true)

#     fᵗ   = [aerosol_optics_trunc_aero1.fᵗ]; 
#     aerosol_optics = [aerosol_optics_trunc_aero1] 

#     Nquad, qp_μ, wt_μ = rt_set_streams(RTM.RadauQuad(), Ltrunc, szas[sza_idx], vza);
#     R, T = RTM.run_RTM(polarization_type, szas[sza_idx], vza, vaz, τRayl, ϖRayl, τAer, ϖAer, fᵗ, qp_μ, wt_μ, maxM, aerosol_optics, GreekRayleigh, τ_abs, RadiativeTransfer.Architectures.GPU()) 

#     R[:,:,:,sza_idx] = R
#     T[:,:,:,sza_idx] = T
# end


#################################################################
# PLOTTING R/T VS. LTRUNC (GPU)
#################################################################


Ltruncs = collect(2:30)
Rs = zeros((9,3,10000,length(Ltruncs)))
Ts = zeros((9,3,10000,length(Ltruncs)))

grid = range(1e7 / 780, 1e7 / 755, length=10000)
τ_abs = zeros(length(grid), length(profile_caltech.p))
compute_absorption_profile!(grid, τ_abs, profile_caltech)

for ltrunc_idx in 1:length(Ltruncs)

    @show ltrunc_idx
    
    Ltrunc = Ltruncs[ltrunc_idx]  
    truncation_type   = PhaseFunction.δBGE{Float32}(Ltrunc, 2.0)      

    model_NAI2_aero1 = make_mie_model(NAI2(), aero1, λ, polarization_type, truncation_type)
    aerosol_optics_NAI2_aero1 = compute_aerosol_optical_properties(model_NAI2_aero1);

    # Truncate:
    aerosol_optics_trunc_aero1 = PhaseFunction.truncate_phase(truncation_type, aerosol_optics_NAI2_aero1; reportFit=true)

    fᵗ   = [aerosol_optics_trunc_aero1.fᵗ]; 
    aerosol_optics = [aerosol_optics_trunc_aero1] 

    Nquad, qp_μ, wt_μ = rt_set_streams(RTM.RadauQuad(), Ltrunc, sza, vza);
    R, T = RTM.run_RTM(polarization_type, sza, vza, vaz, τRayl, ϖRayl, τAer, ϖAer, fᵗ, qp_μ, wt_μ, maxM, aerosol_optics, GreekRayleigh, τ_abs, RadiativeTransfer.Architectures.GPU()) 

    Rs[:,:,:,ltrunc_idx] = R
    Ts[:,:,:,ltrunc_idx] = T
end

anim = @animate for i = 1:24
    plot(1:nSpec, Rs[1,1,:,i], ylims=(-0.1, 0.2), title="Reflectance, Ltrunc=$(Ltruncs[i])")
end
gif(anim, "reflectance.gif", fps = 5)