using Plots
using OptimalEstimation
using LinearAlgebra
using BenchmarkTools


ny = 1000
nx = 100
Sₐ = (Diagonal(abs.(rand(nx))));
Sₑ = (Diagonal((3rand(ny))));
#Sₐ = (abs.(rand(nx,nx)));
#Sₑ = ((10*rand(ny,ny)));
K = rand(ny,nx);
y = randn(ny);

prob = OptimalEstimation.OptimalEstimationProblem(K, Sₑ, y, Sₐ)


#OptimalEstimation.solve(prob) ≈ OptimalEstimation.solve_stable(prob)
@time OptimalEstimation.solve(prob, errorAnalysis=true) 
@time OptimalEstimation.solve_stable(prob,errorAnalysis=true)



##### 


using Revise
using Plots
using RadiativeTransfer
using RadiativeTransfer.CrossSection
using RadiativeTransfer.PhaseFunction
using RadiativeTransfer.RTM
using Distributions
using BenchmarkTools
using Test
using CUDA
using ForwardDiff

device!(3)


function run_auto(x)

    FT = Float64
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
    μ            = [x[1]]    # [0.3,2.0]       # Log mean radius
    σ            = [x[2]]    # [2.0,1.8]       # Log stddev of radius
    r_max        = [30.0]   # [30.0,30.0]     # Maximum radius
    nquad_radius = [2500]   # [2500,2500]     # Number of quadrature points for integrating of size dist.
    nᵣ           = [x[3]]    # [1.3, 1.66]     # Real part of refractive index
    nᵢ           = [x[4]]  # [0.001,0.0003]  # Imag part of refractive index

    # Aerosol vertical distribution profiles
    p₀          = FT[90000.]  # [50000., 20000.] # Pressure peak [Pa]
    σp          = FT[5000.]   # [5000., 2000.]   # Pressure peak width [Pa]

    size_distribution = [LogNormal(log(μ[1]), log(σ[1]))] # [LogNormal(log(μ[1]), log(σ[1])), LogNormal(log(μ[2]), log(σ[2]))]

    # Create the aerosols (needs to be generalized through loops):
    aero1 = make_univariate_aerosol(size_distribution[1], r_max[1], nquad_radius[1], nᵣ[1], nᵢ[1])
    # aero2 = make_univariate_aerosol(size_distribution[2], r_max[2], nquad_radius[2], nᵣ[2], nᵢ[2])

    # Define some details, run aerosol optics
    model_NAI2_aero1 = make_mie_model(NAI2(), aero1, λ, polarization_type, truncation_type)
    aerosol_optics_NAI2_aero1 = compute_aerosol_optical_properties(model_NAI2_aero1);

    # Truncate:
    aerosol_optics_trunc_aero1 = aerosol_optics_NAI2_aero1;
    # aerosol_optics_trunc_aero1 = PhaseFunction.truncate_phase(truncation_type, aerosol_optics_NAI2_aero1; reportFit=true)

    # Define some details, run aerosol optics
    # model_NAI2_aero2 = make_mie_model(NAI2(), aero2, λ, polarization_type, truncation_type)
    # aerosol_optics_NAI2_aero2 = compute_aerosol_optical_properties(model_NAI2_aero2);
    # Truncate:
    # aerosol_optics_trunc_aero2 = PhaseFunction.truncate_phase(truncation_type, aerosol_optics_NAI2_aero2)

    # Rayleigh Greek
    GreekRayleigh = PhaseFunction.get_greek_rayleigh(depol)

    vza = FT[60., 45., 30., 15., 0., 15., 30., 45., 60.]
    vaz = FT[180., 180., 180., 180., 0., 0., 0., 0., 0.]
    sza = FT(60.)

    obs_geom = ObsGeometry(1000.0, sza, vza, vaz)

    Nquad, qp_μ, wt_μ = rt_set_streams(RTM.RadauQuad(), Ltrunc, obs_geom);
    # Nquad, qp_μ, wt_μ = rt_set_streams(RTM.GaussQuadFullSphere(), Ltrunc, sza, vza);

    # In[ ]:


    " Atmospheric Profiles, basics, needs to be refactore entirely"
    file = "/net/fluo/data1/ftp/XYZT_ESE156/Data/MERRA300.prod.assim.inst6_3d_ana_Nv.20150613.hdf.nc4" 
    # file = "MERRA300.prod.assim.inst6_3d_ana_Nv.20150613.hdf.nc4"  
    timeIndex = 2 # There is 00, 06, 12 and 18 in UTC, i.e. 6 hourly data stacked together

    # What latitude do we want? 
    myLat = 34.1377;
    myLon = -118.1253;

    # Read profile (and generate dry/wet VCDs per layer)
    @time profile_caltech_hr = RTM.read_atmos_profile(file, myLat, myLon, timeIndex);
    profile_caltech = RTM.reduce_profile(20, profile_caltech_hr);
    # Compute layer optical thickness for Rayleigh (surface pressure in hPa) 
    τRayl =  RTM.getRayleighLayerOptProp(profile_caltech.psurf / 100, λ, depol, profile_caltech.vcd_dry);
    ϖRayl = ones(length(τRayl));

    # Compute Naer aerosol optical thickness profiles
    τAer_1 = RTM.getAerosolLayerOptProp(1.0, p₀[1], σp[1], profile_caltech.p_levels)
    # τAer_2 = RTM.getAerosolLayerOptProp(0.3, p₀[2], σp[2], profile_caltech.p_levels)

    # Can be done with arbitrary length later:
    τAer = 0.2 * τAer_1; # [τAer_1 τAer_2]
    @show sum(τAer)# , sum(τAer_2)
    ϖAer = [aerosol_optics_NAI2_aero1.ω̃]; # [aerosol_optics_NAI2_aero1.ω̃ aerosol_optics_NAI2_aero2.ω̃];
    fᵗ   = [aerosol_optics_trunc_aero1.fᵗ]; # [aerosol_optics_trunc_aero1.fᵗ aerosol_optics_trunc_aero2.fᵗ];

    aerosol_optics = [aerosol_optics_trunc_aero1] # [aerosol_optics_trunc_aero1 aerosol_optics_trunc_aero2]
    # Aer𝐙⁺⁺ = [aero1_Z⁺⁺] # [aero1_Z⁺⁺, aero2_Z⁺⁺];
    # Aer𝐙⁻⁺ = [aero1_Z⁻⁺] # [aero1_Z⁻⁺, aero2_Z⁻⁺];

    maxM = Ltrunc

    #grid = range(1e7 / 774, 1e7 / 757, length=100);

    #τ_abs = zeros(length(grid), length(profile_caltech.p));
    #compute_absorption_profile!(grid, τ_abs, profile_caltech);
    #@show aerosol_optics[1].greek_coefs
    @time Aer𝐙⁺⁺_curr, Aer𝐙⁻⁺_curr = PhaseFunction.compute_Z_moments(polarization_type, qp_μ, aerosol_optics[1].greek_coefs, 1)
    @show size(Aer𝐙⁺⁺_curr)
    @show Aer𝐙⁺⁺_curr[1,1]
    Aer𝐙⁺⁺_curr = CuArray(Aer𝐙⁺⁺_curr)
    @show Aer𝐙⁺⁺_curr
    return Aer𝐙⁺⁺_curr
    #@time R_GPU, T_GPU = RTM.run_RTM(polarization_type, obs_geom, τRayl, ϖRayl, τAer, ϖAer, fᵗ, qp_μ, wt_μ, maxM, aerosol_optics, GreekRayleigh, τ_abs, RadiativeTransfer.Architectures.GPU());
end
# @time R_CPU, T_CPU = RTM.run_RTM(polarization_type, obs_geom, τRayl, ϖRayl, τAer, ϖAer, fᵗ, qp_μ, wt_μ, maxM, aerosol_optics, GreekRayleigh, τ_abs, RadiativeTransfer.Architectures.CPU());

x = [1.3, 2.0, 1.3, 0.000001];
run_auto(x)
# g = x -> ForwardDiff.jacobian(run_auto, x);

#@test R_CPU ≈ (R_GPU) 