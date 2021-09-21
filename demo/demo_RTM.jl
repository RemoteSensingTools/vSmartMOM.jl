using Revise
using RadiativeTransfer
using RadiativeTransfer.Scattering
using RadiativeTransfer.RTM
using Distributions
using BenchmarkTools

"Generate aerosol optical properties"

# Wavelength (just one for now)
λ = 0.770       # Incident wavelength
depol = 0.0
Ltrunc = 72     # Truncation  
truncation_type   = Scattering.δBGE(Ltrunc, 2.0)

# polarization_type
polarization_type = Stokes_IQUV()

# Quadrature points for RTM
Nquad, qp_μ, wt_μ = rt_set_streams(RTM.RadauQuad(), Ltrunc, 60.0, [0.0, 15.0, 30., 45., 60.])

# Aerosol particle distribution and properties
μ            = [0.3]   # Log mean radius
σ            = [2.0]   # Log stddev of radius
r_max        = [30.0]  # Maximum radius
nquad_radius = [2500]  # Number of quadrature points for integrating of size dist.
nᵣ           = [1.3]   # Real part of refractive index
nᵢ           = [0.001] # Imag part of refractive index

#Aerosol vertical distribution profiles
p₀          = [50000.] # Pressure peak [Pa]
σp          = [5000.]  # Pressure peak width [Pa]

size_distribution = [LogNormal(log(μ[1]), log(σ[1]))] 

# Create the aerosols (needs to be generalized through loops):
aero1 = make_univariate_aerosol(size_distribution[1], r_max[1], nquad_radius[1], nᵣ[1], nᵢ[1])

# Define some details, run aerosol optics
model_NAI2_aero1 = make_mie_model(NAI2(), aero1, λ, polarization_type, truncation_type)
aerosol_optics_NAI2_aero1 = compute_aerosol_optical_properties(model_NAI2_aero1);

# Truncate:
aerosol_optics_trunc_aero1 = Scattering.truncate_phase(truncation_type, aerosol_optics_NAI2_aero1)

# Rayleigh Greek
GreekRayleigh = Scattering.get_greek_rayleigh(depol)

###########

# Define viewing and solar angles
vza = [60., 45., 30., 15., 0., 15., 30., 45., 60.]
vaz = [180., 180., 180., 180., 0., 0., 0., 0., 0.]
sza = 60.
Nquad, qp_μ, wt_μ = rt_set_streams(RTM.RadauQuad(), Ltrunc, sza, vza);

###########

" Atmospheric Profiles, basics, needs to be refactore entirely"
file = "/home/rjeyaram/RadiativeTransfer/src/RTM/MERRA300.prod.assim.inst6_3d_ana_Nv.20150613.hdf.nc4"# "/Users/sanghavi/GDrive/code/github/atm_profiles/MERRA300.prod.assim.inst6_3d_ana_Nv.20150613.hdf.nc4"   
timeIndex = 2 # There is 00, 06, 12 and 18 in UTC, i.e. 6 hourly data stacked together

# What latitude do we want? 
myLat = 34.1377
myLon = -118.1253

# Read profile (and generate dry/wet VCDs per layer)
profile_caltech = RTM.read_atmos_profile(file, myLat, myLon, timeIndex);

###########

# Compute layer optical thickness for Rayleigh (surface pressure in hPa) 
τRayl =  RTM.getRayleighLayerOptProp(profile_caltech.psurf/100, λ, depol, profile_caltech.vcd_dry);
ϖRayl = ones(length(τRayl))

#Compute Naer aerosol optical thickness profiles
τAer_1 = RTM.getAerosolLayerOptProp(1.0, p₀[1], σp[1], profile_caltech.p_levels)

# Can be done with arbitrary length later:
τAer = 0*τAer_1 #[τAer_1 τAer_2]
@show sum(τAer_1)#, sum(τAer_2)
ϖAer = [aerosol_optics_NAI2_aero1.ω̃] #[aerosol_optics_NAI2_aero1.ω̃ aerosol_optics_NAI2_aero2.ω̃];
fᵗ   = [aerosol_optics_trunc_aero1.fᵗ] #[aerosol_optics_trunc_aero1.fᵗ aerosol_optics_trunc_aero2.fᵗ];

(τAer[10,:])
profile_caltech.p_levels[73]

m = 0
RaylZ⁺⁺, RaylZ⁻⁺     = Scattering.compute_Z_moments(polarization_type, qp_μ, GreekRayleigh, m);
aero1_Z⁺⁺, aero1_Z⁻⁺ = Scattering.compute_Z_moments(polarization_type, qp_μ, aerosol_optics_trunc_aero1.greek_coefs, m);
#aero2_Z⁺⁺, aero2_Z⁻⁺ = Scattering.compute_Z_moments(polarization_type, qp_μ, aerosol_optics_trunc_aero2.greek_coefs, m);
aerosol_optics = [aerosol_optics_trunc_aero1] #[aerosol_optics_trunc_aero1 aerosol_optics_trunc_aero2]
Aer𝐙⁺⁺ = [aero1_Z⁺⁺] #[aero1_Z⁺⁺, aero2_Z⁺⁺];
Aer𝐙⁻⁺ = [aero1_Z⁻⁺] #[aero1_Z⁻⁺, aero2_Z⁻⁺];
@show size(τAer[1])
iz = 10
τ, ϖ, Z⁺⁺, Z⁻⁺  = RTM.construct_atm_layer(τRayl[iz], τAer[iz,:], ϖRayl[iz], ϖAer, fᵗ, RaylZ⁺⁺, RaylZ⁻⁺, Aer𝐙⁺⁺, Aer𝐙⁻⁺)
@show τ, ϖ
@show τAer[iz], τRayl[iz]


###########

@btime R,T = RTM.run_RTM(polarization_type, sza, vza, vaz, τRayl,ϖRayl, τAer, ϖAer, fᵗ, qp_μ, wt_μ, 3, aerosol_optics, GreekRayleigh)

###########

using Plots
using DelimitedFiles
S1=readdlm("/Users/sanghavi/GDrive/code/github/RadiativeTransfer.jl/src/Notebooks/Rayl_orig.dat", ' ')
vza=[-60. -45. -30. -15. 0. 15. 30. 45. 60.]
I_rayl0=S1[:,1]
Q_rayl0=S1[:,2]
U_rayl0=S1[:,3]
V_rayl0=S1[:,4]
p1=plot(vza', [R[:,1], I_rayl0], title="I, 775 nm", label=["julia" "bolshoi"], lw=2, legend=:topright, size = (600, 400))
p2=plot(vza', [R[:,2], Q_rayl0], title="Q, 775 nm", label=["julia" "bolshoi"], lw=2, legend=:topleft, size = (600, 400))
plot(p1, p2, layout=(2, 1))