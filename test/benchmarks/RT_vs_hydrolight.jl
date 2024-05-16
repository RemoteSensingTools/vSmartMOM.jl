##
using CUDA
device!(0)
using Revise
using vSmartMOM, vSmartMOM.CoreRT, vSmartMOM.SolarModel
using vSmartMOM.InelasticScattering
using Statistics
using Interpolations
using InstrumentOperator #for convolution of hires spectrum to instrument grid
using ImageFiltering
using Distributions
using Plots
using DelimitedFiles

# Benchmarks: http://web.gps.caltech.edu/~vijay/Rayleigh_Scattering_Tables/STOKES/
##

# Load YAML files into parameter struct
parameters = 
    parameters_from_yaml("test/test_parameters/BOA_var_test.yaml");
FT = Float64
RS_type = InelasticScattering.noRS(
        fscattRayl  = [FT(1)],
        ϖ_Cabannes  = [FT(1)], 
        bandSpecLim = [],
        iBand       = [1],
        F₀          = zeros(FT,1,1));

# Create model struct (precomputes optical properties) from parameters
model      = model_from_parameters(parameters);
n_bands = length(parameters.spec_bands);
tot_nspec = sum([length(parameters.spec_bands[i]) for i=1:length(parameters.spec_bands)]) - (n_bands-1)
n_cam = length(parameters.vza)
R = zeros(FT, n_cam, parameters.polarization_type.n, tot_nspec)
T = zeros(FT, n_cam, parameters.polarization_type.n, tot_nspec)

iBand = 1
tot_ν=zeros(tot_nspec);
tot_F₀ = zeros(tot_nspec);
spec_end = 1
Fdir = zeros(tot_nspec);
for iBand=1:n_bands
    #### Compute all Raman properties
    ν = model.params.spec_bands[iBand]   
    ν̃ = mean(ν);
    spec_start = spec_end
    spec_end += length(ν)-1
    # Find central reference index for RRS:
    i_ref = argmin(abs.(ν .- ν̃))

    # Effective temperature for Raman calculations
    effT = 300.  #(model.profile.vcd_dry' * model.profile.T) / sum(model.profile.vcd_dry);
    # Define RS type
    # Compute N2 and O2
    n2,o2 = InelasticScattering.getRamanAtmoConstants(ν̃,effT);

    
    T_sun = 5777. # K
    P = planck_spectrum_wn(T_sun, ν)
    Tsolar = solar_transmission_from_file("/home/sanghavi/code/github/vSmartMOM.jl/src/SolarModel/solar.out")
    Tsolar_interp = LinearInterpolation(Tsolar[4:end, 1], Tsolar[4:end, 2])
    F₀ = zeros(length(P));
    RS_type.F₀ = zeros(model.params.polarization_type.n, length(P))
    for i=1:length(P) 
        sol_trans = Tsolar_interp(ν[i]);
        F₀[i] = sol_trans * P[i];
        RS_type.F₀[1,i] = F₀[i];
    end

    @show spec_start, spec_end, length(ν) 
    tot_ν[spec_start:spec_end] = ν 
    tot_F₀[spec_start:spec_end] = F₀ 
    R[:,:,spec_start:spec_end], T[:,:,spec_start:spec_end], _, _ = CoreRT.rt_run_test(RS_type,model,iBand);

    tot_τ = sum(model.τ_abs[iBand], dims=2) + sum(model.τ_rayl[iBand], dims=2) + sum(model.τ_aer[iBand][1,:,:], dims=2) #+ sum(sum(model.τ_aer[iBand], dims=3), dims=1)
    Fdir[spec_start:spec_end] = F₀ .* exp.(-(tot_τ)/cosd(parameters.sza))
end
T_aer_fn0p3_ρ0p5 = T
F_aer_fn0p3_ρ0p5 = Fdir
T_aer_fn1p0_ρ0p5 = T
F_aer_fn1p0_ρ0p5 = Fdir
T_aer_fn10p0_ρ0p5 = T
F_aer_fn10p0_ρ0p5 = Fdir
T_aer_cr0p3_ρ0p5 = T
F_aer_cr0p3_ρ0p5 = Fdir
T_aer_cr1p0_ρ0p5 = T
F_aer_cr1p0_ρ0p5 = Fdir
T_aer_cr10p0_ρ0p5 = T
F_aer_cr10p0_ρ0p5 = Fdir
T_rayl_ρ0p5 = T
F_rayl_ρ0p5 = Fdir

T_aer_fn10p0_ρ0p0 = T
F_aer_fn10p0_ρ0p0 = Fdir
T_aer_fn1p0_ρ0p0 = T
F_aer_fn1p0_ρ0p0 = Fdir
T_aer_fn0p3_ρ0p0 = T
F_aer_fn0p3_ρ0p0 = Fdir
T_aer_cr0p3_ρ0p0 = T
F_aer_cr0p3_ρ0p0 = Fdir
T_aer_cr1p0_ρ0p0 = T
F_aer_cr1p0_ρ0p0 = Fdir
T_aer_cr10p0_ρ0p0 = T
F_aer_cr10p0_ρ0p0 = Fdir
T_rayl_ρ0p0 = T
F_rayl_ρ0p0 = Fdir
# SZA = 30




fname1 = "/home/sanghavi/hydrolight_vs_fullRT/T_SZA30_ρ0p0_fn_0p3.dat"
tmp1 = [ν T_aer_fn0p3_ρ0p0[1,1,:] T_aer_fn0p3_ρ0p0[1,2,:]]# T[i,3,:]]
fname2 = "/home/sanghavi/hydrolight_vs_fullRT/F_SZA30_ρ0p0_fn_0p3.dat"
tmp2 = [ν F_aer_fn0p3_ρ0p0]
writedlm(fname1, tmp1)
writedlm(fname2, tmp2)

fname1 = "/home/sanghavi/hydrolight_vs_fullRT/T_SZA30_ρ0p0_fn_1p0.dat"
tmp1 = [ν T_aer_fn1p0_ρ0p0[1,1,:] T_aer_fn1p0_ρ0p0[1,2,:]]# T[i,3,:]]
fname2 = "/home/sanghavi/hydrolight_vs_fullRT/F_SZA30_ρ0p0_fn_1p0.dat"
tmp2 = [ν F_aer_fn1p0_ρ0p0]
writedlm(fname1, tmp1)
writedlm(fname2, tmp2)

fname1 = "/home/sanghavi/hydrolight_vs_fullRT/T_SZA30_ρ0p0_fn_10p0.dat"
tmp1 = [ν T_aer_fn10p0_ρ0p0[1,1,:] T_aer_fn10p0_ρ0p0[1,2,:]]# T[i,3,:]]
fname2 = "/home/sanghavi/hydrolight_vs_fullRT/F_SZA30_ρ0p0_fn_10p0.dat"
tmp2 = [ν F_aer_fn10p0_ρ0p0]
writedlm(fname1, tmp1)
writedlm(fname2, tmp2)

fname1 = "/home/sanghavi/hydrolight_vs_fullRT/T_SZA30_ρ0p0_cr_0p3.dat"
tmp1 = [ν T_aer_cr0p3_ρ0p0[1,1,:] T_aer_cr0p3_ρ0p0[1,2,:]]# T[i,3,:]]
fname2 = "/home/sanghavi/hydrolight_vs_fullRT/F_SZA30_ρ0p0_cr_0p3.dat"
tmp2 = [ν F_aer_cr0p3_ρ0p0]
writedlm(fname1, tmp1)
writedlm(fname2, tmp2)

fname1 = "/home/sanghavi/hydrolight_vs_fullRT/T_SZA30_ρ0p0_cr_1p0.dat"
tmp1 = [ν T_aer_cr1p0_ρ0p0[1,1,:] T_aer_cr1p0_ρ0p0[1,2,:]]# T[i,3,:]]
fname2 = "/home/sanghavi/hydrolight_vs_fullRT/F_SZA30_ρ0p0_cr_1p0.dat"
tmp2 = [ν F_aer_cr1p0_ρ0p0]
writedlm(fname1, tmp1)
writedlm(fname2, tmp2)

fname1 = "/home/sanghavi/hydrolight_vs_fullRT/T_SZA30_ρ0p0_cr_10p0.dat"
tmp1 = [ν T_aer_cr10p0_ρ0p0[1,1,:] T_aer_cr10p0_ρ0p0[1,2,:]]# T[i,3,:]]
fname2 = "/home/sanghavi/hydrolight_vs_fullRT/F_SZA30_ρ0p0_cr_10p0.dat"
tmp2 = [ν F_aer_cr10p0_ρ0p0]
writedlm(fname1, tmp1)
writedlm(fname2, tmp2)

fname1 = "/home/sanghavi/hydrolight_vs_fullRT/T_SZA30_ρ0p0_rayl.dat"
tmp1 = [ν T_rayl_ρ0p0[1,1,:] T_rayl_ρ0p0[1,2,:]]# T[i,3,:]]
fname2 = "/home/sanghavi/hydrolight_vs_fullRT/F_SZA30_ρ0p0_rayl.dat"
tmp2 = [ν F_rayl_ρ0p0]
writedlm(fname1, tmp1)
writedlm(fname2, tmp2)

#=================#

fname1 = "/home/sanghavi/hydrolight_vs_fullRT/T_SZA30_ρ0p5_fn_0p3.dat"
tmp1 = [ν T_aer_fn0p3_ρ0p5[1,1,:] T_aer_fn0p3_ρ0p5[1,2,:]]# T[i,3,:]]
fname2 = "/home/sanghavi/hydrolight_vs_fullRT/F_SZA30_ρ0p5_fn_0p3.dat"
tmp2 = [ν F_aer_fn0p3_ρ0p5]
writedlm(fname1, tmp1)
writedlm(fname2, tmp2)

fname1 = "/home/sanghavi/hydrolight_vs_fullRT/T_SZA30_ρ0p5_fn_1p0.dat"
tmp1 = [ν T_aer_fn1p0_ρ0p5[1,1,:] T_aer_fn1p0_ρ0p5[1,2,:]]# T[i,3,:]]
fname2 = "/home/sanghavi/hydrolight_vs_fullRT/F_SZA30_ρ0p5_fn_1p0.dat"
tmp2 = [ν F_aer_fn1p0_ρ0p5]
writedlm(fname1, tmp1)
writedlm(fname2, tmp2)

fname1 = "/home/sanghavi/hydrolight_vs_fullRT/T_SZA30_ρ0p5_fn_10p0.dat"
tmp1 = [ν T_aer_fn10p0_ρ0p5[1,1,:] T_aer_fn10p0_ρ0p5[1,2,:]]# T[i,3,:]]
fname2 = "/home/sanghavi/hydrolight_vs_fullRT/F_SZA30_ρ0p5_fn_10p0.dat"
tmp2 = [ν F_aer_fn10p0_ρ0p5]
writedlm(fname1, tmp1)
writedlm(fname2, tmp2)

fname1 = "/home/sanghavi/hydrolight_vs_fullRT/T_SZA30_ρ0p5_cr_0p3.dat"
tmp1 = [ν T_aer_cr0p3_ρ0p5[1,1,:] T_aer_cr0p3_ρ0p5[1,2,:]]# T[i,3,:]]
fname2 = "/home/sanghavi/hydrolight_vs_fullRT/F_SZA30_ρ0p5_cr_0p3.dat"
tmp2 = [ν F_aer_cr0p3_ρ0p5]
writedlm(fname1, tmp1)
writedlm(fname2, tmp2)

fname1 = "/home/sanghavi/hydrolight_vs_fullRT/T_SZA30_ρ0p5_cr_1p0.dat"
tmp1 = [ν T_aer_cr1p0_ρ0p5[1,1,:] T_aer_cr1p0_ρ0p5[1,2,:]]# T[i,3,:]]
fname2 = "/home/sanghavi/hydrolight_vs_fullRT/F_SZA30_ρ0p5_cr_1p0.dat"
tmp2 = [ν F_aer_cr1p0_ρ0p5]
writedlm(fname1, tmp1)
writedlm(fname2, tmp2)

fname1 = "/home/sanghavi/hydrolight_vs_fullRT/T_SZA30_ρ0p5_cr_10p0.dat"
tmp1 = [ν T_aer_cr10p0_ρ0p5[1,1,:] T_aer_cr10p0_ρ0p5[1,2,:]]# T[i,3,:]]
fname2 = "/home/sanghavi/hydrolight_vs_fullRT/F_SZA30_ρ0p5_cr_10p0.dat"
tmp2 = [ν F_aer_cr10p0_ρ0p5]
writedlm(fname1, tmp1)
writedlm(fname2, tmp2)

fname1 = "/home/sanghavi/hydrolight_vs_fullRT/T_SZA30_ρ0p5_rayl.dat"
tmp1 = [ν T_rayl_ρ0p5[1,1,:] T_rayl_ρ0p5[1,2,:]]# T[i,3,:]]
fname2 = "/home/sanghavi/hydrolight_vs_fullRT/F_SZA30_ρ0p5_rayl.dat"
tmp2 = [ν F_rayl_ρ0p5]
writedlm(fname1, tmp1)
writedlm(fname2, tmp2)

aod = ["0p3","1p0","10p0"]
sz = ["fn", "cr"]
surf = ["ρ0p0", "ρ0p5"]

T_rayl = zeros(2,2,length(ν))
T_aer  = zeros(2,2,3,2,length(ν))
F_rayl = zeros(2,length(ν))
F_aer  = zeros(2,2,3,length(ν))
for i=1:2 #surf
    fname1 = "/home/sanghavi/hydrolight_vs_fullRT/T_SZA30_"*surf[i]*"_rayl.dat"
    fname2 = "/home/sanghavi/hydrolight_vs_fullRT/F_SZA30_"*surf[i]*"_rayl.dat"
    a = readdlm(fname1)
    ν = a[:,1]
    T_rayl[i,1,:] = a[:,2]
    T_rayl[i,2,:] = a[:,3]

    b = readdlm(fname2)
    F_rayl[i,:] = b[:,2]

    for j=1:2 #sz
        for k=1:3 #aod
            fname1 = "/home/sanghavi/hydrolight_vs_fullRT/T_SZA30_"*surf[i]*"_"*sz[j]*"_"*aod[k]*".dat"
            fname2 = "/home/sanghavi/hydrolight_vs_fullRT/F_SZA30_"*surf[i]*"_"*sz[j]*"_"*aod[k]*".dat"
            a = readdlm(fname1)
            T_aer[i,j,k,1,:] = a[:,2] 
            T_aer[i,j,k,2,:] = a[:,3]

            b = readdlm(fname2)
            F_aer[i,j,k,:] = b[:,2]
        end
    end
end

#Rayl, bright vs black
l = @layout [a1 a2]
p=plot(1e7./ν, 100*T_rayl[1,1,:]./F_rayl[1,:], label="", ylabel="ΔI/I [%]", xlabel="λ [nm]")
p=plot!(1e7./ν, 100*T_rayl[2,1,:]./F_rayl[1,:], label="")
q=plot(1e7./ν, 100*T_rayl[1,2,:]./F_rayl[1,:], label="ρ=0", ylabel="Q/I [%]", xlabel="λ [nm]")
q=plot!(1e7./ν, 100*T_rayl[2,2,:]./F_rayl[1,:], label="ρ=0.5")
plot(p, q, layout = l, title = "Pure Rayleigh scattering", titlefont = font(10))
savefig("/home/sanghavi/hydrolight_vs_fullRT/Rayl_RT_vs_hydrolight_SZA30_400nm.png")

#Aer, fine, AOD 0.3, bright vs black
l = @layout [a1 a2]
p=plot(1e7./ν, 100*T_aer[1,1,1,1,:]./F_aer[1,1,1,:], label="", ylabel="ΔI/I [%]", xlabel="λ [nm]")
p=plot!(1e7./ν, 100*T_aer[2,1,1,1,:]./F_aer[1,1,1,:], label="")
q=plot(1e7./ν, 100*T_aer[1,1,1,2,:]./F_aer[1,1,1,:], label="ρ=0", ylabel="Q/I [%]", xlabel="λ [nm]")
q=plot!(1e7./ν, 100*T_aer[2,1,1,2,:]./F_aer[1,1,1,:], label="ρ=0.5")
plot(p, q, layout = l, title = "Aerosol AOD=0.3, r₀=0.1 μm", titlefont = font(10))
savefig("/home/sanghavi/hydrolight_vs_fullRT/Aer_fn0p3_RT_vs_hydrolight_SZA30_400nm.png")
#Aer, coarse, AOD 0.3, bright vs black
l = @layout [a1 a2]
p=plot(1e7./ν, 100*T_aer[1,2,1,1,:]./F_aer[1,2,1,:], label="", ylabel="ΔI/I [%]", xlabel="λ [nm]")
p=plot!(1e7./ν, 100*T_aer[2,2,1,1,:]./F_aer[1,2,1,:], label="")
q=plot(1e7./ν, 100*T_aer[1,2,1,2,:]./F_aer[1,2,1,:], label="ρ=0", ylabel="Q/I [%]", xlabel="λ [nm]")
q=plot!(1e7./ν, 100*T_aer[2,2,1,2,:]./F_aer[1,2,1,:], label="ρ=0.5")
plot(p, q, layout = l, title = "Aerosol AOD=0.3, r₀=0.5 μm", titlefont = font(10))
savefig("/home/sanghavi/hydrolight_vs_fullRT/Aer_cr0p3_RT_vs_hydrolight_SZA30_400nm.png")

#Aer, fine, AOD 1.0, bright vs black
l = @layout [a1 a2]
p=plot(1e7./ν, 100*T_aer[1,1,2,1,:]./F_aer[1,1,2,:], label="", ylabel="ΔI/I [%]", xlabel="λ [nm]")
p=plot!(1e7./ν, 100*T_aer[2,1,2,1,:]./F_aer[1,1,2,:], label="")
q=plot(1e7./ν, 100*T_aer[1,1,2,2,:]./F_aer[1,1,2,:], label="ρ=0", ylabel="Q/I [%]", xlabel="λ [nm]")
q=plot!(1e7./ν, 100*T_aer[2,1,2,2,:]./F_aer[1,1,2,:], label="ρ=0.5")
plot(p, q, layout = l, title = "Aerosol AOD=1.0, r₀=0.1 μm", titlefont = font(10))
savefig("/home/sanghavi/hydrolight_vs_fullRT/Aer_fn1p0_RT_vs_hydrolight_SZA30_400nm.png")

#Aer, coarse, AOD 1.0, bright vs black
l = @layout [a1 a2]
p=plot(1e7./ν, 100*T_aer[1,2,2,1,:]./F_aer[1,2,2,:], label="", ylabel="ΔI/I [%]", xlabel="λ [nm]")
p=plot!(1e7./ν, 100*T_aer[2,2,2,1,:]./F_aer[1,2,2,:], label="")
q=plot(1e7./ν, 100*T_aer[1,2,2,2,:]./F_aer[1,2,2,:], label="ρ=0", ylabel="Q/I [%]", xlabel="λ [nm]")
q=plot!(1e7./ν, 100*T_aer[2,2,2,2,:]./F_aer[1,2,2,:], label="ρ=0.5")
plot(p, q, layout = l, title = "Aerosol AOD=1.0, r₀=0.5 μm", titlefont = font(10))
savefig("/home/sanghavi/hydrolight_vs_fullRT/Aer_cr1p0_RT_vs_hydrolight_SZA30_400nm.png")

#Aer, fine, AOD 10.0, bright vs black
l = @layout [a1 a2]
p=plot(1e7./ν, 100*T_aer[1,1,3,1,:]./F_aer[1,1,3,:], label="", ylabel="ΔI/I [%]", xlabel="λ [nm]")
p=plot!(1e7./ν, 100*T_aer[2,1,3,1,:]./F_aer[1,1,3,:], label="")
q=plot(1e7./ν, 100*T_aer[1,1,3,2,:]./F_aer[1,1,3,:], label="ρ=0", ylabel="Q/I [%]", xlabel="λ [nm]")
q=plot!(1e7./ν, 100*T_aer[2,1,3,2,:]./F_aer[1,1,3,:], label="ρ=0.5")
plot(p, q, layout = l, title = "Aerosol AOD=10.0, r₀=0.1 μm", titlefont = font(10))
savefig("/home/sanghavi/hydrolight_vs_fullRT/Aer_fn10p0_RT_vs_hydrolight_SZA30_400nm.png")

#Aer, coarse, AOD 10.0, bright vs black
l = @layout [a1 a2]
p=plot(1e7./ν, 100*T_aer[1,2,3,1,:]./F_aer[1,2,3,:], label="", ylabel="ΔI/I [%]", xlabel="λ [nm]")
p=plot!(1e7./ν, 100*T_aer[2,2,3,1,:]./F_aer[1,2,3,:], label="")
q=plot(1e7./ν, 100*T_aer[1,2,3,2,:]./F_aer[1,2,3,:], label="ρ=0", ylabel="Q/I [%]", xlabel="λ [nm]")
q=plot!(1e7./ν, 100*T_aer[2,2,3,2,:]./F_aer[1,2,3,:], label="ρ=0.5")
plot(p, q, layout = l, title = "Aerosol AOD=10.0, r₀=0.5 μm", titlefont = font(10))
savefig("/home/sanghavi/hydrolight_vs_fullRT/Aer_cr10p0_RT_vs_hydrolight_SZA30_400nm.png")

#==================================================================================#
# writing routine
#==================================================================================#
# black surface, Rayleigh
for i = 1:5
    tmp = [ν T[i,1,:] T[i,2,:]]# T[i,3,:]]
    fname = "/home/sanghavi/hydrolight_vs_fullRT/rayl_SZA30_rho0_"*string(i)*".dat"
    writedlm(fname, tmp)
end

# black surface, fine aerosol (1.48, 0.0, 0.1, log(1.6), log(2), sqrt(log(3/2))
for i = 1:5
    tmp = [ν T[i,1,:] T[i,2,:]]# T[i,3,:]]
    fname = "/home/sanghavi/hydrolight_vs_fullRT/faer_SZA30_rho0_"*string(i)*".dat"
    writedlm(fname, tmp)
end

# black surface, coarse aerosol (1.6, 0.0, 0.5, log(1.3), log(3), sqrt(log(5/3))
for i = 1:5
    tmp = [ν T[i,1,:] T[i,2,:]]# T[i,3,:]]
    fname = "/home/sanghavi/hydrolight_vs_fullRT/faer_SZA30_rho0_"*string(i)*".dat"
    writedlm(fname, tmp)
end

# ρ = 0.2, Rayleigh
for i = 1:5
    tmp = [ν T[i,1,:] T[i,2,:]]# T[i,3,:]]
    fname = "/home/sanghavi/hydrolight_vs_fullRT/rayl_SZA30_rho0p2_"*string(i)*".dat"
    writedlm(fname, tmp)
end

# ρ = 0.2, fine aerosol
for i = 1:5
    tmp = [ν T[i,1,:] T[i,2,:]]# T[i,3,:]]
    fname = "/home/sanghavi/hydrolight_vs_fullRT/faer_SZA30_rho0p2_"*string(i)*".dat"
    writedlm(fname, tmp)
end

# ρ = 0.2, coarse aerosol
for i = 1:5
    tmp = [ν T[i,1,:] T[i,2,:]]# T[i,3,:]]
    fname = "/home/sanghavi/hydrolight_vs_fullRT/faer_SZA30_rho0p2_"*string(i)*".dat"
    writedlm(fname, tmp)
end

# ρ = 0.5, Rayleigh
for i = 1:5
    tmp = [ν T[i,1,:] T[i,2,:]]# T[i,3,:]]
    fname = "/home/sanghavi/hydrolight_vs_fullRT/rayl_SZA30_rho0p5_"*string(i)*".dat"
    writedlm(fname, tmp)
end

# ρ = 0.5, fine aerosol
for i = 1:5
    tmp = [ν T[i,1,:] T[i,2,:]]# T[i,3,:]]
    fname = "/home/sanghavi/hydrolight_vs_fullRT/faer_SZA30_rho0p5_"*string(i)*".dat"
    writedlm(fname, tmp)
end

# ρ = 0.5, coarse aerosol
for i = 1:5
    tmp = [ν T[i,1,:] T[i,2,:]]# T[i,3,:]]
    fname = "/home/sanghavi/hydrolight_vs_fullRT/faer_SZA30_rho0p5_"*string(i)*".dat"
    writedlm(fname, tmp)
end

#==================================================================================#
#reading routine
#==================================================================================#

# black surface, Rayleigh
rayl = [zeros(3,length(ν)) for i=1:5]
for i = 1:5
    fname = fname = "/home/sanghavi/hydrolight_vs_fullRT/rayl_SZA30_rho0_"*string(i)*".dat"
    rayl[i] = readdlm(fname)
end

# black surface, fine aerosol (1.48, 0.0000001, 0.1, log(1.6), log(2), sqrt(log(3/2))
faer = [zeros(3,length(ν)) for i=1:5]
for i = 1:5
    fname = fname = "/home/sanghavi/hydrolight_vs_fullRT/faer_SZA30_rho0_"*string(i)*".dat"
    faer[i] = readdlm(fname)
end

# black surface, coarse aerosol (1.6, 0.001, 0.5, log(1.3), log(3), sqrt(log(5/3))
caer = [zeros(3,length(ν)) for i=1:5]
for i = 1:5
    fname = fname = "/home/sanghavi/hydrolight_vs_fullRT/caer_SZA30_rho0_"*string(i)*".dat"
    caer[i] = readdlm(fname)
end

# ρ = 0.2, Rayleigh
rayl = [zeros(3,length(ν)) for i=1:5]
for i = 1:5
    fname = fname = "/home/sanghavi/hydrolight_vs_fullRT/rayl_SZA30_rho0p2_"*string(i)*".dat"
    rayl[i] = readdlm(fname)
end

# ρ = 0.2, fine aerosol
faer = [zeros(3,length(ν)) for i=1:5]
for i = 1:5
    fname = fname = "/home/sanghavi/hydrolight_vs_fullRT/faer_SZA30_rho0p2_"*string(i)*".dat"
    faer[i] = readdlm(fname)
end

# ρ = 0.2, coarse aerosol
caer = [zeros(3,length(ν)) for i=1:5]
for i = 1:5
    fname = fname = "/home/sanghavi/hydrolight_vs_fullRT/caer_SZA30_rho0p2_"*string(i)*".dat"
    caer[i] = readdlm(fname)
end

# ρ = 0.5, Rayleigh
rayl = [zeros(3,length(ν)) for i=1:5]
for i = 1:5
    fname = fname = "/home/sanghavi/hydrolight_vs_fullRT/rayl_SZA30_rho0p5_"*string(i)*".dat"
    rayl[i] = readdlm(fname)
end

# ρ = 0.5, fine aerosol
faer = [zeros(3,length(ν)) for i=1:5]
for i = 1:5
    fname = fname = "/home/sanghavi/hydrolight_vs_fullRT/faer_SZA30_rho0p5_"*string(i)*".dat"
    faer[i] = readdlm(fname)
end

# ρ = 0.5, coarse aerosol
caer = [zeros(3,length(ν)) for i=1:5]
for i = 1:5
    fname = fname = "/home/sanghavi/hydrolight_vs_fullRT/caer_SZA30_rho0p5_"*string(i)*".dat"
    caer[i] = readdlm(fname)
end


RS_type.F₀ = zeros(model.params.polarization_type.n, length(ν))
RS_type.F₀[1,:].=1.

=#
#RnoRS_ms, TnoRS_ms, _, _ = CoreRT.rt_run_test_ms(noRS([0.0],[1.0], Any[],[1]),[0,3],model,iBand);
#RnoRS_test, TnoRS_test, _, _ = CoreRT.rt_run_test_ms(noRS,[0,3],model,iBand);
#RnoRS, TnoRS, _, _ = CoreRT.rt_run_test(noRS([0.0],[1.0], Any[],[1]),model,iBand);
#sR_ms, T_ms, ieR_ms, ieT_ms = CoreRT.rt_run_test_ms(RS_type,[0,3], model,iBand);
#make F₀ the concatenated solar function across all bands contained in iBand
# R_ms, T_ms, ieR_ms, ieT_ms = 
    # CoreRT.rt_run_test_ms(RS_type,model,iBand);
#@show RS_type.ϖ_Cabannes
#RS_type.ϖ_Cabannes = [1.0]
#@show RS_type.ϖ_Cabannes
#R, T, ieR, ieT = CoreRT.rt_run_test(RS_type,model,iBand);
#R_ss, T_ss, ieR_ss, ieT_ss = CoreRT.rt_run_test_ss(RS_type,model,iBand);

#=
RS_type.F₀ = zeros(model.params.polarization_type.n, length(P))
for i=1:length(P)
    sol_trans = Tsolar_interp(ν[i]);
    F₀[i] = sol_trans * P[i];
    RS_type.F₀[1,i] = F₀[i];
end=# 
#parameters = 
#    parameters_from_yaml("test/test_parameters/FraunhoferMockParameters.yaml");

#parameters.depol = 0.041362343961163395 #0.1032 #(0.028: Cabannes), (0.1032: Rayleigh) 
#model      = model_from_parameters(parameters);

#RnoRS, TnoRS, _, _ = CoreRT.rt_run_test(RS_type,model,iBand);
#RnoRS_ss, TnoRS_ss, _, _ = CoreRT.rt_run_test_ss(RS_type,model,iBand);

#=
R_test, T_test, ieR_test, ieT_test = CoreRT.rt_run_test(RS_type,model,iBand);

# You can now run multiple bands like this (list at the end of band!)
RnoRS_test, TnoRS_test, _, _ = CoreRT.rt_run_test(vSmartMOM.noRS(),model,[1,1]);


R_test, T_test, ieR_test, ieT_test = CoreRT.rt_run_test(RS_type,model,1);
=#

#===Convolution of hires spectral simulations to instrument grid===#
#x = (1e7/415):0.3:(1e7/385)
#ν = (1e7/460):0.3:(1e7/420)
x = -40:0.1:40
#kernel = InstrumentOperator.create_instrument_kernel(Normal(0, 12.5), x) #defining a Gaussian kernel for convulution in wavenumber space
kernel = InstrumentOperator.create_instrument_kernel(Normal(0, 12.5), x)
#I_conv = InstrumentOperator.conv_spectra(kernel, )
I_conv_noRS = imfilter(RnoRS[1,1,:], kernel)
I_conv = imfilter(R[1,1,:], kernel)
ieI_conv = imfilter(ieR[1,1,:], kernel)

Q_conv_noRS = imfilter(RnoRS[1,2,:], kernel)
Q_conv = imfilter(R[1,2,:], kernel)
ieQ_conv = imfilter(ieR[1,2,:], kernel)

convfct = 1e7./ν.^2   # to convert radiance units from mW/m²-str-cm⁻¹ to mW/m²-str-nm
l = @layout [a1 a2 ; b1 b2; c1 c2]
p1 = plot(1e7./ν, RnoRS[1,1,:].*convfct, linecolor=:black, xlims=[385,415])
p1 = plot!(1e7./ν, I_conv_noRS.*convfct, linewidth = 3, linecolor=:red, ylabel="[mW/m²/str/nm]", yguidefontsize=8)
p2 = plot(1e7./ν, (R[1,1,:].-RnoRS[1,1,:].+ieR[1,1,:]).*convfct, linecolor=:black, xlabel = "λ [nm]", xlims=[385,415])
p2 = plot!(1e7./ν, (I_conv.-I_conv_noRS.+ieI_conv).*convfct, linewidth = 3, linecolor=:red, ylabel="[mW/m²/str/nm]", yguidefontsize=8)
p3 = plot(1e7./ν, 100*(R[1,1,:].+ieR[1,1,:].-RnoRS[1,1,:])./RnoRS[1,1,:], linecolor=:black, ylabel="[%]", xlims=[385,415], ylims=[-5,150])
p3 = plot!(1e7./ν, 100*(I_conv.+ieI_conv.-I_conv_noRS)./I_conv_noRS, linewidth = 2,  linecolor=:red)


q1 = plot(1e7./ν, RnoRS[1,2,:].*convfct, linecolor=:black, xlims=[385,415])
q1 = plot!(1e7./ν, Q_conv_noRS.*convfct, linewidth = 3, linecolor=:red)
q2 = plot(1e7./ν, (R[1,2,:].-RnoRS[1,2,:].+ieR[1,2,:]).*convfct, linecolor=:black, xlabel = "λ [nm]", xlims=[385,415])
q2 = plot!(1e7./ν, (Q_conv.-Q_conv_noRS.+ieQ_conv).*convfct, linewidth = 3, linecolor=:red)
q3 = plot(1e7./ν, 100*(R[1,2,:].+ieR[1,2,:].-RnoRS[1,2,:])./RnoRS[1,2,:], linecolor=:black, xlims=[385,415])
q3 = plot!(1e7./ν, 100*(Q_conv.+ieQ_conv.-Q_conv_noRS)./Q_conv_noRS, linewidth = 2, linecolor=:red, ylims=[-0.6,1.5])
plot(p1, q1, p2, q2, p3, q3, layout = l, link=:x, legend = false, title = ["I₀ (no RS)" "Q₀ (no RS)" "I₁-I₀" "Q₁-Q₀" "(I₁-I₀)/I₀ [%]" "(Q₁-Q₀)/Q₀ [%]"], titlefont = font(10))
savefig("RingEffect.png")

l = @layout [a1 a2]
p1 = plot(1e7./ν, 100*ieR[1,1,:]./R[1,1,:], linecolor=:black, xlims=[385,415], ylims=[-5,150])
p1 = plot!(1e7./ν, 100*ieI_conv./I_conv, linewidth = 2, linecolor=:red, xlabel = "λ [nm]")
q1 = plot(1e7./ν, 100*ieR[1,2,:]./R[1,2,:], linecolor=:black, xlims=[385,415], ylims=[-0.1,2])
q1 = plot!(1e7./ν, 100*ieQ_conv./Q_conv, linewidth = 2, linecolor=:red, xlabel = "λ [nm]")
plot(p1, q1, layout = l, legend = false, title = ["Iᵢ/Iₑ [%]" "Qᵢ/Qₑ [%]"], titlefont = font(10))
savefig("RingSpectrum.png")

using DelimitedFiles
raylout_rrs = readdlm("out_ray_rrs_385nm_405nm.dat")
raylout_nors = readdlm("out_ray_nors_385nm_405nm.dat")
raylout_rrs_ss = readdlm("out_ray_rrs_ss_385nm_405nm.dat")
aerout_rrs = readdlm("out_aer_rrs_385nm_405nm.dat")
aerout_nors = readdlm("out_aer_nors_385nm_405nm.dat")
aerout_rrs_ss = readdlm("out_aer_rrs_ss_385nm_405nm.dat")

x = -40:0.6:40
#kernel = InstrumentOperator.create_instrument_kernel(Normal(0, 12.5), x) #defining a Gaussian kernel for convulution in wavenumber space
kernel = InstrumentOperator.create_instrument_kernel(Normal(0, 12.5), x)
#I_conv = InstrumentOperator.conv_spectra(kernel, )

ν = raylout_nors[:,1]
#Rayleigh
rI_noRS = zeros(length(ν))
rQ_noRS = zeros(length(ν))
rI = zeros(length(ν))
rQ = zeros(length(ν))
ierI = zeros(length(ν))
ierQ = zeros(length(ν))
rI_ss = zeros(length(ν))
rQ_ss = zeros(length(ν))
ierI_ss = zeros(length(ν))
ierQ_ss = zeros(length(ν))

rI_conv_noRS = zeros(length(ν))
rQ_conv_noRS = zeros(length(ν))
rI_conv = zeros(length(ν))
rQ_conv = zeros(length(ν))
ierI_conv = zeros(length(ν))
ierQ_conv = zeros(length(ν))
rI_conv_ss = zeros(length(ν))
rQ_conv_ss = zeros(length(ν))
ierI_conv_ss = zeros(length(ν))
ierQ_conv_ss = zeros(length(ν))

rI_noRS = raylout_nors[:,2]
rQ_noRS = raylout_nors[:,3]
rI = raylout_rrs[:,2]
rQ = raylout_rrs[:,3]
ierI = raylout_rrs[:,4]
ierQ = raylout_rrs[:,5]
rI_ss = raylout_rrs_ss[:,2]
rQ_ss = raylout_rrs_ss[:,3]
ierI_ss = raylout_rrs_ss[:,4]
ierQ_ss = raylout_rrs_ss[:,5]

rI_conv_noRS = imfilter(rI_noRS, kernel)
rQ_conv_noRS = imfilter(rQ_noRS, kernel)
rI_conv = imfilter(rI, kernel)
rQ_conv = imfilter(rQ, kernel)
ierI_conv = imfilter(ierI, kernel)
ierQ_conv = imfilter(ierQ, kernel)
rI_conv_ss = imfilter(rI_ss, kernel)
rQ_conv_ss = imfilter(rQ_ss, kernel)
ierI_conv_ss = imfilter(ierI_ss, kernel)
ierQ_conv_ss = imfilter(ierQ_ss, kernel)

#Aerosol
aI_noRS = zeros(length(ν))
aQ_noRS = zeros(length(ν))
aI = zeros(length(ν))
aQ = zeros(length(ν))
ieaI = zeros(length(ν))
ieaQ = zeros(length(ν))
aI_ss = zeros(length(ν))
aQ_ss = zeros(length(ν))
ieaI_ss = zeros(length(ν))
ieaQ_ss = zeros(length(ν))

aI_conv_noRS = zeros(length(ν))
aQ_conv_noRS = zeros(length(ν))
aI_conv = zeros(length(ν))
aQ_conv = zeros(length(ν))
ieaI_conv = zeros(length(ν))
ieaQ_conv = zeros(length(ν))
aI_conv_ss = zeros(length(ν))
aQ_conv_ss = zeros(length(ν))
ieaI_conv_ss = zeros(length(ν))
ieaQ_conv_ss = zeros(length(ν))

aI_noRS = aerout_nors[:,2]
aQ_noRS = aerout_nors[:,3]
aI = aerout_rrs[:,2]
aQ = aerout_rrs[:,3]
ieaI = aerout_rrs[:,4]
ieaQ = aerout_rrs[:,5]
aI_ss = aerout_rrs_ss[:,2]
aQ_ss = aerout_rrs_ss[:,3]
ieaI_ss = aerout_rrs_ss[:,4]
ieaQ_ss = aerout_rrs_ss[:,5]

aI_conv_noRS = imfilter(aI_noRS, kernel)
aQ_conv_noRS = imfilter(aQ_noRS, kernel)
aI_conv = imfilter(aI, kernel)
aQ_conv = imfilter(aQ, kernel)
ieaI_conv = imfilter(ieaI, kernel)
ieaQ_conv = imfilter(ieaQ, kernel)
aI_conv_ss = imfilter(aI_ss, kernel)
aQ_conv_ss = imfilter(aQ_ss, kernel)
ieaI_conv_ss = imfilter(ieaI_ss, kernel)
ieaQ_conv_ss = imfilter(ieaQ_ss, kernel)

convfct = 1e7./ν.^2   # to convert radiance units from mW/m²-str-cm⁻¹ to mW/m²-str-nm

p1 = plot(1e7./ν, rQ_noRS./rI_noRS, linewidth = 3, linecolor=:black, label="Rayleigh", ylims=[0.2,0.5])
p1 = plot!(1e7./ν, (rQ+ierQ)./(rI+ierI), linecolor=:black, label="")
p1 = plot!(1e7./ν, (rQ_conv+ierQ_conv)./(rI_conv+ierI_conv), linewidth = 3, linecolor=:red, label="")
p1 = plot!(1e7./ν, aQ_noRS./aI_noRS, linewidth = 3, linecolor=:grey, label="Aerosol")
p1 = plot!(1e7./ν, (aQ+ieaQ)./(aI+ieaI), linecolor=:grey, label="")
p1 = plot!(1e7./ν, (aQ_conv+ieaQ_conv)./(aI_conv+ieaI_conv), linewidth = 3, linecolor=:red, xlabel = "λ [nm]", ylabel = "DOLP", label="")
plot(p1, title= "Zenith obs, θ₀=63ᵒ, A=0.05")#, legend = false)
savefig("compDOLP.png")


# New
#=
rayl10_rrs_ss = [ν R10_ss[1,1,:] R10_ss[1,2,:] R10_ss[1,3,:] ieR10_ss[1,1,:] ieR10_ss[1,2,:] ieR10_ss[1,3,:]]
rayl70_rrs_ss = [ν R70_ss[1,1,:] R70_ss[1,2,:] R70_ss[1,3,:] ieR70_ss[1,1,:] ieR70_ss[1,2,:] ieR70_ss[1,3,:]]
rayl10_rrs = [ν R10[1,1,:] R10[1,2,:] R10[1,3,:] ieR10[1,1,:] ieR10[1,2,:] ieR10[1,3,:]]
rayl70_rrs = [ν R70[1,1,:] R70[1,2,:] R70[1,3,:] ieR70[1,1,:] ieR70[1,2,:] ieR70[1,3,:]]
rayl10_nors_ss = [ν R10noRS_ss[1,1,:] R10noRS_ss[1,2,:] R10noRS_ss[1,3,:]]
rayl70_nors_ss = [ν R70noRS_ss[1,1,:] R70noRS_ss[1,2,:] R70noRS_ss[1,3,:]]
rayl10_nors = [ν R10noRS[1,1,:] R10noRS[1,2,:] R10noRS[1,3,:]]
rayl70_nors = [ν R70noRS[1,1,:] R70noRS[1,2,:] R70noRS[1,3,:]]

writedlm("out_ray10_rrs_ss_385nm_405nm.dat", rayl10_rrs_ss)
writedlm("out_ray70_rrs_ss_385nm_405nm.dat", rayl70_rrs_ss)
writedlm("out_ray70_rrs_385nm_405nm.dat", rayl70_rrs)
writedlm("out_ray10_rrs_385nm_405nm.dat", rayl10_rrs)
writedlm("out_ray10_nors_385nm_405nm.dat", rayl10_nors)
writedlm("out_ray70_nors_385nm_405nm.dat", rayl70_nors)
writedlm("out_ray70_nors_ss_385nm_405nm.dat", rayl70_nors_ss)
writedlm("out_ray10_nors_ss_385nm_405nm.dat", rayl10_nors_ss)
=#
rayl10_rrs = readdlm("out_ray10_rrs_385nm_405nm.dat")
rayl10_nors = readdlm("out_ray10_nors_385nm_405nm.dat")
rayl10_rrs_ss = readdlm("out_ray10_rrs_ss_385nm_405nm.dat")
rayl10_nors_ss = readdlm("out_ray10_nors_ss_385nm_405nm.dat")
rayl70_rrs = readdlm("out_ray70_rrs_385nm_405nm.dat")
rayl70_nors = readdlm("out_ray70_nors_385nm_405nm.dat")
rayl70_rrs_ss = readdlm("out_ray70_rrs_ss_385nm_405nm.dat")
rayl70_nors_ss = readdlm("out_ray70_nors_ss_385nm_405nm.dat")

ν = rayl10_nors[:,1]

I10 = rayl10_rrs[:,2]
Q10 = rayl10_rrs[:,3]
U10 = rayl10_rrs[:,4]
ieI10 = rayl10_rrs[:,5]
ieQ10 = rayl10_rrs[:,6]
ieU10 = rayl10_rrs[:,7]

I10_ss = rayl10_rrs_ss[:,2]
Q10_ss = rayl10_rrs_ss[:,3]
U10_ss = rayl10_rrs_ss[:,4]
ieI10_ss = rayl10_rrs_ss[:,5]
ieQ10_ss = rayl10_rrs_ss[:,6]
ieU10_ss = rayl10_rrs_ss[:,7]

I10_nors = rayl10_nors[:,2]
Q10_nors = rayl10_nors[:,3]
U10_nors = rayl10_nors[:,4]

I10_nors_ss = rayl10_nors_ss[:,2]
Q10_nors_ss = rayl10_nors_ss[:,3]
U10_nors_ss = rayl10_nors_ss[:,4]

I70 = rayl70_rrs[:,2]
Q70 = rayl70_rrs[:,3]
U70 = rayl70_rrs[:,4]
ieI70 = rayl70_rrs[:,5]
ieQ70 = rayl70_rrs[:,6]
ieU70 = rayl70_rrs[:,7]

I70_ss = rayl70_rrs_ss[:,2]
Q70_ss = rayl70_rrs_ss[:,3]
U70_ss = rayl70_rrs_ss[:,4]
ieI70_ss = rayl70_rrs_ss[:,5]
ieQ70_ss = rayl70_rrs_ss[:,6]
ieU70_ss = rayl70_rrs_ss[:,7]

I70_nors = rayl70_nors[:,2]
Q70_nors = rayl70_nors[:,3]
U70_nors = rayl70_nors[:,4]

I70_nors_ss = rayl70_nors_ss[:,2]
Q70_nors_ss = rayl70_nors_ss[:,3]
U70_nors_ss = rayl70_nors_ss[:,4]

x = -40:0.3:40
#kernel = InstrumentOperator.create_instrument_kernel(Normal(0, 12.5), x) #defining a Gaussian kernel for convulution in wavenumber space
#kernel = InstrumentOperator.create_instrument_kernel(Normal(0, 6.25), x)
kernel = InstrumentOperator.create_instrument_kernel(Normal(0, 7.64), x)

#I_conv = InstrumentOperator.conv_spectra(kernel, )
I10_conv_noRS = imfilter(I10_nors, kernel)
I10_conv_noRS_ss = imfilter(I10_nors_ss, kernel)
I10_conv = imfilter(I10, kernel)
ieI10_conv = imfilter(ieI10, kernel)
I10_ss_conv = imfilter(I10_ss, kernel)
ieI10_ss_conv = imfilter(ieI10_ss, kernel)

Q10_conv_noRS = imfilter(Q10_nors, kernel)
Q10_conv_noRS_ss = imfilter(Q10_nors_ss, kernel)
Q10_conv = imfilter(Q10, kernel)
ieQ10_conv = imfilter(ieQ10, kernel)
Q10_ss_conv = imfilter(Q10_ss, kernel)
ieQ10_ss_conv = imfilter(ieQ10_ss, kernel)

U10_conv_noRS = imfilter(U10_nors, kernel)
U10_conv_noRS_ss = imfilter(U10_nors_ss, kernel)
U10_conv = imfilter(U10, kernel)
ieU10_conv = imfilter(ieU10, kernel)
U10_ss_conv = imfilter(U10_ss, kernel)
ieU10_ss_conv = imfilter(ieU10_ss, kernel)

I70_conv_noRS = imfilter(I70_nors, kernel)
I70_conv_noRS_ss = imfilter(I70_nors_ss, kernel)
I70_conv = imfilter(I70, kernel)
ieI70_conv = imfilter(ieI70, kernel)
I70_ss_conv = imfilter(I70_ss, kernel)
ieI70_ss_conv = imfilter(ieI70_ss, kernel)

Q70_conv_noRS = imfilter(Q70_nors, kernel)
Q70_conv_noRS_ss = imfilter(Q70_nors_ss, kernel)
Q70_conv = imfilter(Q70, kernel)
ieQ70_conv = imfilter(ieQ70, kernel)
Q70_ss_conv = imfilter(Q70_ss, kernel)
ieQ70_ss_conv = imfilter(ieQ70_ss, kernel)

U70_conv_noRS = imfilter(U70_nors, kernel)
U70_conv_noRS_ss = imfilter(U70_nors_ss, kernel)
U70_conv = imfilter(U70, kernel)
ieU70_conv = imfilter(ieU70, kernel)
U70_ss_conv = imfilter(U70_ss, kernel)
ieU70_ss_conv = imfilter(ieU70_ss, kernel)

convfct = 1e7./ν.^2
l = @layout [a1 a2 ; b1 b2; c1 c2]

#p1 = plot(1e7./ν, 100*((I10.+ieI10)./I10_nors.-1), linecolor=:black, ylabel="R₁ [%]", alpha=0.3)
p1 = plot(1e7./ν, 100*((I10_conv.+ieI10_conv)./I10_conv_noRS.-1), linewidth = 3, linecolor=:black)
#p1 = plot!(1e7./ν, 100*((I10_ss.+ieI10_ss)./I10_nors_ss.-1), linecolor=:red, alpha=0.3)
p1 = plot!(1e7./ν, 100*((I10_ss_conv.+ieI10_ss_conv)./I10_conv_noRS_ss.-1), linewidth = 3, linecolor=:red)

#p2 = plot(1e7./ν, 100*((Q10.+ieQ10)./Q10_nors.-1), linecolor=:black, ylabel="R₂ [%]", alpha=0.3)
p2 = plot(1e7./ν, 100*((Q10_conv.+ieQ10_conv)./Q10_conv_noRS.-1), linewidth = 3, linecolor=:black)
#p2 = plot!(1e7./ν, 100*((Q10_ss.+ieQ10_ss)./Q10_nors_ss.-1), linecolor=:red, alpha=0.3)
p2 = plot!(1e7./ν, 100*((Q10_ss_conv.+ieQ10_ss_conv)./Q10_conv_noRS_ss.-1), linewidth = 3, linecolor=:red)

#p3 = plot(1e7./ν, 100*((U10.+ieU10)./U10_nors.-1), linecolor=:black, ylabel="R₃ [%]", alpha=0.3)
p3 = plot(1e7./ν, 100*((U10_conv.+ieU10_conv)./U10_conv_noRS.-1), linewidth = 3, linecolor=:black)
#p3 = plot!(1e7./ν, 100*((U10_ss.+ieU10_ss)./U10_nors_ss.-1), linecolor=:red, alpha=0.3)
p3 = plot!(1e7./ν, 100*((U10_ss_conv.+ieU10_ss_conv)./U10_conv_noRS_ss.-1), linewidth = 3, linecolor=:red)

#q1 = plot(1e7./ν, 100*((I70.+ieI70)./I70_nors.-1), linecolor=:black, ylabel="R₁ [%]", alpha=0.3)
q1 = plot(1e7./ν, 100*((I70_conv.+ieI70_conv)./I70_conv_noRS.-1), linewidth = 3, linecolor=:black)
#q1 = plot!(1e7./ν, 100*((I70_ss.+ieI70_ss)./I70_nors_ss.-1), linecolor=:red, alpha=0.3)
q1 = plot!(1e7./ν, 100*((I70_ss_conv.+ieI70_ss_conv)./I70_conv_noRS_ss.-1), linewidth = 3, linecolor=:red)

#q2 = plot(1e7./ν, 100*((Q70.+ieQ70)./Q70_nors.-1), linecolor=:black, ylabel="R₂ [%]", alpha=0.3)
q2 = plot(1e7./ν, 100*((Q70_conv.+ieQ70_conv)./Q70_conv_noRS.-1), linewidth = 3, linecolor=:black)
#q2 = plot!(1e7./ν, 100*((Q70_ss.+ieQ70_ss)./Q70_nors_ss.-1), linecolor=:red, alpha=0.3)
q2 = plot!(1e7./ν, 100*((Q70_ss_conv.+ieQ70_ss_conv)./Q70_conv_noRS_ss.-1), linewidth = 3, linecolor=:red)

#q3 = plot(1e7./ν, 100*((U70.+ieU70)./U70_nors.-1), linecolor=:black, ylabel="R₃ [%]", alpha=0.3)
q3 = plot(1e7./ν, 100*((U70_conv.+ieU70_conv)./U70_conv_noRS.-1), linewidth = 3, linecolor=:black)
#q3 = plot!(1e7./ν, 100*((U70_ss.+ieU70_ss)./U70_nors_ss.-1), linecolor=:red, alpha=0.3)
q3 = plot!(1e7./ν, 100*((U70_ss_conv.+ieU70_ss_conv)./U70_conv_noRS_ss.-1), linewidth = 3, linecolor=:red)

plot(p1, q1, p2, q2, p3, q3, layout = l, link=:y, legend = false, title = ["SZA=10ᵒ" "SZA=70ᵒ" "" "" "" ""], titlefont = font(10))
savefig("abc.png")

l = @layout [a1 a2 ; b1 b2; c1 c2]

#p1 = plot(1e7./ν, 100*((I10.+ieI10)./I10_nors.-1), linecolor=:black, ylabel="R₁ [%]", alpha=0.3)
p1 = plot(1e7./ν, 100*((I10_conv.+ieI10_conv)./I10_conv_noRS.-1), linewidth = 3, linecolor=:black)
#p1 = plot!(1e7./ν, 100*((I10_ss.+ieI10_ss)./I10_nors_ss.-1), linecolor=:red, alpha=0.3)
p1 = plot!(1e7./ν, 100*((I10_conv.+ieI10_ss_conv)./I10_conv_noRS.-1), linewidth = 3, linecolor=:red)

#p2 = plot(1e7./ν, 100*((Q10.+ieQ10)./Q10_nors.-1), linecolor=:black, ylabel="R₂ [%]", alpha=0.3)
p2 = plot(1e7./ν, 100*((Q10_conv.+ieQ10_conv)./Q10_conv_noRS.-1), linewidth = 3, linecolor=:black)
#p2 = plot!(1e7./ν, 100*((Q10_ss.+ieQ10_ss)./Q10_nors_ss.-1), linecolor=:red, alpha=0.3)
p2 = plot!(1e7./ν, 100*((Q10_conv.+ieQ10_ss_conv)./Q10_conv_noRS.-1), linewidth = 3, linecolor=:red)

#p3 = plot(1e7./ν, 100*((U10.+ieU10)./U10_nors.-1), linecolor=:black, ylabel="R₃ [%]", alpha=0.3)
p3 = plot(1e7./ν, 100*((U10_conv.+ieU10_conv)./U10_conv_noRS.-1), linewidth = 3, linecolor=:black)
#p3 = plot!(1e7./ν, 100*((U10_ss.+ieU10_ss)./U10_nors_ss.-1), linecolor=:red, alpha=0.3)
p3 = plot!(1e7./ν, 100*((U10_conv.+ieU10_ss_conv)./U10_conv_noRS.-1), linewidth = 3, linecolor=:red)

#q1 = plot(1e7./ν, 100*((I70.+ieI70)./I70_nors.-1), linecolor=:black, ylabel="R₁ [%]", alpha=0.3)
q1 = plot(1e7./ν, 100*((I70_conv.+ieI70_conv)./I70_conv_noRS.-1), linewidth = 3, linecolor=:black)
#q1 = plot!(1e7./ν, 100*((I70_ss.+ieI70_ss)./I70_nors_ss.-1), linecolor=:red, alpha=0.3)
q1 = plot!(1e7./ν, 100*((I70_conv.+ieI70_ss_conv)./I70_conv_noRS.-1), linewidth = 3, linecolor=:red)

#q2 = plot(1e7./ν, 100*((Q70.+ieQ70)./Q70_nors.-1), linecolor=:black, ylabel="R₂ [%]", alpha=0.3)
q2 = plot(1e7./ν, 100*((Q70_conv.+ieQ70_conv)./Q70_conv_noRS.-1), linewidth = 3, linecolor=:black)
#q2 = plot!(1e7./ν, 100*((Q70_ss.+ieQ70_ss)./Q70_nors_ss.-1), linecolor=:red, alpha=0.3)
q2 = plot!(1e7./ν, 100*((Q70_conv.+ieQ70_ss_conv)./Q70_conv_noRS.-1), linewidth = 3, linecolor=:red)

#q3 = plot(1e7./ν, 100*((U70.+ieU70)./U70_nors.-1), linecolor=:black, ylabel="R₃ [%]", alpha=0.3)
q3 = plot(1e7./ν, 100*((U70_conv.+ieU70_conv)./U70_conv_noRS.-1), linewidth = 3, linecolor=:black)
#q3 = plot!(1e7./ν, 100*((U70_ss.+ieU70_ss)./U70_nors_ss.-1), linecolor=:red, alpha=0.3)
q3 = plot!(1e7./ν, 100*((U70_conv.+ieU70_ss_conv)./U70_conv_noRS.-1), linewidth = 3, linecolor=:red)

plot(p1, q1, p2, q2, p3, q3, layout = l, link=:y, legend = false, title = ["SZA=10ᵒ" "SZA=70ᵒ" "" "" "" ""], titlefont = font(10))
savefig("abc.png")