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
using Printf 

FT = Float64

# Load YAML files into parameter struct
parameters = 
    parameters_from_yaml("/home/sanghavi/code/github/vSmartMOM.jl/test/test_parameters/noabs_parameters2_SIF_grid.yaml");
#parameters = parameters_from_yaml("test/test_parameters/O2Parameters2.yaml");
# Create model struct (precomputes optical properties) from parameters
model      = model_from_parameters(parameters);
model.τ_aer[1][1,:,:] .= 0;
model.τ_aer[1][1,:,end] .= 100.
model.τ_rayl[1][:,end-1] .+= model.τ_rayl[1][:,end];
model.τ_rayl[1][:,end] .= 0;
a=[]
for i=0:20
    push!(a, acosd(i/20))  
end
sza=reverse(Int.(ceil.(a[8:21])))
ρ = zeros(FT,21) #ρ = zeros(FT,15)
ρ_str = []
for iρ = 1:21 #15
    ρ[iρ] = (iρ-1)*0.05
    for i=1:length(parameters.spec_bands)
        parameters.brdf[i].albedo = ρ[iρ]
    end
    push!(ρ_str, replace(string(round(ρ[iρ], digits=2)),"."=>"p"))
end
psurf=[1000 750 500]
J_SIF = readdlm("/home/sanghavi/code/github/vSmartMOM.jl/src/SIF_emission/sif-spectra.csv", ',') 
#J_SIF = readdlm("/home/sanghavi/code/github/vSmartMOM.jl/src/SIF_emission/PC1_SIFSpectra_allSpecies.csv", ',') 
νSIF = reverse(1e7./J_SIF[2:end,1])
jSIF = reverse(J_SIF[2:end,2]*(0.5*π/maximum(J_SIF[2:end,2]))) # mW/m²/nm
jSIF .*= 1e7./(νSIF.^2) # mW/m²/cm⁻¹
SIF_interp = LinearInterpolation(νSIF, jSIF)
Tsolar = solar_transmission_from_file("/home/sanghavi/code/github/vSmartMOM.jl/src/SolarModel/solar.out")
Tsolar_interp = LinearInterpolation(Tsolar[4:end, 1], Tsolar[4:end, 2])
n_bands = length(parameters.spec_bands);
T_sun = 5777. # K

isurf = 1
iρ    = 1
iA    = 2

parameters.sza = 1.0*sza[iA]
model.obs_geom = CoreRT.ObsGeometry(parameters.sza, parameters.vza, parameters.vaz, parameters.obs_alt)
# Set quadrature points for streams
model.quad_points = CoreRT.rt_set_streams(parameters.quadrature_type, 
                                    parameters.l_trunc, 
                                    model.obs_geom, 
                                    parameters.polarization_type, 
                                    array_type(parameters.architecture))

#model.obs_geom.sza = 0.0;#1.0*sza[iA]
#@show parameters.sza, sza[iA]
#@show parameters.vza


for iBand=1:n_bands
    parameters.brdf[iBand].albedo = ρ[iρ]
    #model.params.brdf[iBand].albedo = ρ[iρ]
    #@show model.params.brdf[iBand].albedo
end
overlap_ν = 250 #230
Δν = parameters.spec_bands[1][2] - parameters.spec_bands[1][1]; #Make sure all bands have the same spectral resolution Δν
n_overlap = Int(floor(overlap_ν/Δν))
tot_nspec = sum([length(parameters.spec_bands[i]) for i=1:length(parameters.spec_bands)]) - 2*(n_bands)*n_overlap
n_cam = length(parameters.vza) #for generating the grid, the SZA and the VZA will be swapped, so that the resulting Rs and ieRs will need to be multiplied by cos(VZA) for the true value 
RnoRS = zeros(FT, n_cam, parameters.polarization_type.n, tot_nspec)
TnoRS = zeros(FT, n_cam, parameters.polarization_type.n, tot_nspec)
R = zeros(FT, n_cam, parameters.polarization_type.n, tot_nspec)
T = zeros(FT, n_cam, parameters.polarization_type.n, tot_nspec)
ieR = zeros(FT, n_cam, parameters.polarization_type.n, tot_nspec)
ieT = zeros(FT, n_cam, parameters.polarization_type.n, tot_nspec)

iBand = 1
tot_ν=zeros(tot_nspec);
tot_F₀ = zeros(tot_nspec);
spec_end = 0
#Fdir = zeros(tot_nspec);

# Common operations for all bands
#### Compute all Raman properties
ν̃ = 0.5*(model.params.spec_bands[1][1]+model.params.spec_bands[end][end])

# Find central reference index for RRS:
#i_ref = argmin(abs.(ν .- ν̃))ν̃ = 0.5*(model.params.spec_bands[1][1]+model.params.spec_bands[end][end])

# TODO_VS: λ_vs_in (get input)
# TODO_VS: ν_vs_in (convert to wavenumbers)
# Effective temperature for Raman calculations
effT = (model.profile.vcd_dry' * model.profile.T) / sum(model.profile.vcd_dry);
n2,o2 = InelasticScattering.getRamanAtmoConstants(ν̃,effT);

RS_type0 = InelasticScattering.noRS(
            fscattRayl  = [FT(1)],
            ϖ_Cabannes  = [FT(1)], 
            bandSpecLim = [],
            iBand       = [1],
            F₀          = zeros(FT,1,1),
            SIF₀        = zeros(FT,1,1));            

RS_type1 = InelasticScattering.RRS(
            n2=n2,
            o2=o2,
            greek_raman = InelasticScattering.GreekCoefs([FT(1)], [FT(1)], [FT(1)], [FT(1)], [FT(1)], [FT(1)]),
            fscattRayl  = [FT(1)],
            ϖ_Cabannes  = [FT(1)], 
            ϖ_λ₁λ₀      = zeros(FT,1),
            i_λ₁λ₀      = zeros(Int,1), 
            Z⁻⁺_λ₁λ₀    = zeros(FT,1,1), 
            Z⁺⁺_λ₁λ₀    = zeros(FT,1,1), 
            i_ref       = 0, # this is a redundant parameter
            n_Raman     = 0,
            F₀          = zeros(FT,1,1),
            SIF₀        = zeros(FT,1,1));
iBAnd =1
ν = model.params.spec_bands[iBand]
parameters.brdf[iBand].albedo = ρ[iρ]
spec_start = spec_end+1
spec_end += length(ν) - 2*n_overlap

# Compute Raman SSA properties:
CoreRT.getRamanSSProp!(RS_type1, 1e7/mean(ν), ν);

P = planck_spectrum_wn(T_sun, ν) * 2.1629e-05 * π  # mW/m²-cm⁻¹

F₀ = zeros(length(P));
SIF₀ = zeros(length(ν))
RS_type0.F₀ = zeros(model.params.polarization_type.n, length(P))
RS_type1.F₀ = zeros(model.params.polarization_type.n, length(P))
RS_type0.SIF₀ = zeros(model.params.polarization_type.n, length(ν))
RS_type1.SIF₀ = zeros(model.params.polarization_type.n, length(ν))
F₀ .= 1.0 # = Tsolar_interp.(ν) .* P;
SIF₀ .= 0.0 # = SIF_interp.(ν); # 
RS_type0.F₀[1,:] = F₀; #1.0 #
RS_type1.F₀[1,:] = F₀;
RS_type0.SIF₀[1,:] = SIF₀; #1.0 #
RS_type1.SIF₀[1,:] = SIF₀;

f_CabN2 = "/home/sanghavi/data/RamanSIFgrid/RamanXS/effCoeff_Cabannes_N2.dat" 
f_CabO2 = "/home/sanghavi/data/RamanSIFgrid/RamanXS/effCoeff_Cabannes_O2.dat" 
f_RRSN2 = "/home/sanghavi/data/RamanSIFgrid/RamanXS/effCoeff_RRS_N2.dat" 
f_RRSO2 = "/home/sanghavi/data/RamanSIFgrid/RamanXS/effCoeff_RRS_O2.dat" 

CabN2 = [0 n2.effCoeff.σ_Rayl_coeff]
CabO2 = [0 o2.effCoeff.σ_Rayl_coeff]
Δν_N2 = []
Δν_O2 = []
kscatt_N2 = []
kscatt_O2 = []

append!(Δν_N2, n2.effCoeff.Δν̃_RoRaman_coeff_JtoJp2)
append!(Δν_N2, n2.effCoeff.Δν̃_RoRaman_coeff_JtoJm2)
append!(kscatt_N2, n2.effCoeff.σ_RoRaman_coeff_JtoJp2)
append!(kscatt_N2, n2.effCoeff.σ_RoRaman_coeff_JtoJm2)

append!(Δν_O2, o2.effCoeff.Δν̃_RoRaman_coeff_JtoJp2)
append!(Δν_O2, o2.effCoeff.Δν̃_RoRaman_coeff_JtoJm2)
append!(kscatt_O2, o2.effCoeff.σ_RoRaman_coeff_JtoJp2)
append!(kscatt_O2, o2.effCoeff.σ_RoRaman_coeff_JtoJm2)

writedlm(f_CabN2, CabN2)
writedlm(f_CabO2, CabO2)
writedlm(f_RRSN2, [Δν_N2 kscatt_N2])
writedlm(f_RRSO2, [Δν_O2 kscatt_O2])


ν₀ = ν̃

plot([ν₀,ν₀], [n2.effCoeff.σ_Rayl_coeff.*ν₀.^4, n2.effCoeff.σ_Rayl_coeff.*ν₀.^4]*1e-3, seriestype=:sticks, lw=2, label="N₂ σₑ(x10⁻³)", color=:red, yformatter=x -> string(@sprintf("%.2e", x)), xlabel="wavenumber [cm⁻¹]", ylabel="Cross sections [cm²]")
plot!(ν₀.+n2.effCoeff.Δν̃_RoRaman_coeff_JtoJp2, n2.effCoeff.σ_RoRaman_coeff_JtoJp2.*(ν₀.+n2.effCoeff.Δν̃_RoRaman_coeff_JtoJp2).^4, seriestype=:sticks, lw=2, color=:blue, label="", yformatter=x -> string(@sprintf("%.2e", x)))
plot!(ν₀.+n2.effCoeff.Δν̃_RoRaman_coeff_JtoJm2, n2.effCoeff.σ_RoRaman_coeff_JtoJm2.*(ν₀.+n2.effCoeff.Δν̃_RoRaman_coeff_JtoJm2).^4, seriestype=:sticks, lw=2, color=:blue, label="N₂ σᵢₑ", yformatter=x -> string(@sprintf("%.2e", x)))
savefig("/home/sanghavi/data/RRS_intercomparison/n2_XS.png")

plot([1e7/ν₀,1e7/ν₀], [n2.effCoeff.σ_Rayl_coeff.*ν₀.^4, n2.effCoeff.σ_Rayl_coeff.*ν₀.^4]*1e-3, seriestype=:sticks, lw=2, label="N₂ σₑ(x10⁻³)", color=:red, yformatter=x -> string(@sprintf("%.2e", x)), xlabel="wavelength [nm]", ylabel="Cross sections [cm²]")
plot!(1e7./(ν₀.+n2.effCoeff.Δν̃_RoRaman_coeff_JtoJp2), n2.effCoeff.σ_RoRaman_coeff_JtoJp2.*(ν₀.+n2.effCoeff.Δν̃_RoRaman_coeff_JtoJp2).^4, seriestype=:sticks, lw=2, color=:blue, label="", yformatter=x -> string(@sprintf("%.2e", x)))
plot!(1e7./(ν₀.+n2.effCoeff.Δν̃_RoRaman_coeff_JtoJm2), n2.effCoeff.σ_RoRaman_coeff_JtoJm2.*(ν₀.+n2.effCoeff.Δν̃_RoRaman_coeff_JtoJm2).^4, seriestype=:sticks, lw=2, color=:blue, label="N₂ σᵢₑ", yformatter=x -> string(@sprintf("%.2e", x)))
savefig("/home/sanghavi/data/RRS_intercomparison/n2_XS_wl.png")

plot([ν₀,ν₀], [o2.effCoeff.σ_Rayl_coeff.*ν₀.^4, o2.effCoeff.σ_Rayl_coeff.*ν₀.^4]*1e-3, seriestype=:sticks, lw=2, label="O₂ σₑ(x10⁻³)", color=:red, yformatter=x -> string(@sprintf("%.2e", x)), xlabel="wavenumber [cm⁻¹]", ylabel="Cross sections [cm²]")
plot!(ν₀.+o2.effCoeff.Δν̃_RoRaman_coeff_JtoJp2, o2.effCoeff.σ_RoRaman_coeff_JtoJp2.*(ν₀.+o2.effCoeff.Δν̃_RoRaman_coeff_JtoJp2).^4, seriestype=:sticks, lw=2, color=:blue, label="", yformatter=x -> string(@sprintf("%.2e", x)))
plot!(ν₀.+o2.effCoeff.Δν̃_RoRaman_coeff_JtoJm2, o2.effCoeff.σ_RoRaman_coeff_JtoJm2.*(ν₀.+o2.effCoeff.Δν̃_RoRaman_coeff_JtoJm2).^4, seriestype=:sticks, lw=2, color=:blue, label="O₂ σᵢₑ", yformatter=x -> string(@sprintf("%.2e", x)))
savefig("/home/sanghavi/data/RRS_intercomparison/o2_XS.png")

plot([1e7/ν₀,1e7/ν₀], [o2.effCoeff.σ_Rayl_coeff.*ν₀.^4, o2.effCoeff.σ_Rayl_coeff.*ν₀.^4]*1e-3, seriestype=:sticks, lw=2, label="O₂ σₑ(x10⁻³)", color=:red, yformatter=x -> string(@sprintf("%.2e", x)), xlabel="wavelength [nm]", ylabel="Cross sections [cm²]")
plot!(1e7./(ν₀.+o2.effCoeff.Δν̃_RoRaman_coeff_JtoJp2), o2.effCoeff.σ_RoRaman_coeff_JtoJp2.*(ν₀.+o2.effCoeff.Δν̃_RoRaman_coeff_JtoJp2).^4, seriestype=:sticks, lw=2, color=:blue, label="", yformatter=x -> string(@sprintf("%.2e", x)))
plot!(1e7./(ν₀.+o2.effCoeff.Δν̃_RoRaman_coeff_JtoJm2), o2.effCoeff.σ_RoRaman_coeff_JtoJm2.*(ν₀.+o2.effCoeff.Δν̃_RoRaman_coeff_JtoJm2).^4, seriestype=:sticks, lw=2, color=:blue, label="O₂ σᵢₑ", yformatter=x -> string(@sprintf("%.2e", x)))
savefig("/home/sanghavi/data/RRS_intercomparison/o2_XS_wl.png")

plot([0], [RS_type1.ϖ_Cabannes]*1e-3, seriestype=:sticks, lw=2, label="Cabannes (x10⁻³)", color=:red, yformatter=x -> string(@sprintf("%.2e", x)), xlabel="Δν [cm⁻¹]", ylabel="Single scattering albedo, ϖ")
plot!(RS_type1.i_λ₁λ₀*0.1, RS_type1.ϖ_λ₁λ₀, seriestype=:sticks, lw=2, label="RRS", color=:blue, yformatter=x -> string(@sprintf("%.2e", x)), xlabel="Δν [cm⁻¹]", ylabel="Single scattering albedo, ϖ")
savefig("/home/sanghavi/data/RRS_intercomparison/RRS_ssa.png")
