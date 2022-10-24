using Revise
using Plots
using Pkg
# Pkg.activate(".")
using vSmartMOM
using vSmartMOM.Architectures
using vSmartMOM.Absorption
using vSmartMOM.Scattering
using vSmartMOM.CoreRT
using vSmartMOM.SolarModel
using InstrumentOperator
using Interpolations
using Polynomials
using ForwardDiff 
using Distributions
using NCDatasets
using Unitful
using CanopyOptics
using TimerOutputs
using Parameters
using LinearAlgebra

include("test/sandbox/canopy_tools.jl")
## Atmospheric Radiative Transfer

# Load parameters from file
parameters = parameters_from_yaml("test/test_parameters/3BandParameters_canopy.yaml")
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
#bands = (1,3);
# Indices within that band:
indices = (92:885,114:845,50:916);
#indices = (92:885,50:916);
# Geo Index (footprint,sounding):
GeoInd = [5,5000];
oco_sounding = InstrumentOperator.getMeasurement(oco, bands, indices, GeoInd);
# Need to force Rayleigh and Aerosols here:
parameters.p   = [oco_sounding.p_half; oco_sounding.p_half[end]+5]
parameters.q   = [oco_sounding.q; oco_sounding.q[end]]
parameters.T   = [oco_sounding.T; oco_sounding.T[end]]# .+ 1.0 #.+ x[15]
parameters.sza = oco_sounding.sza
parameters.vza = [oco_sounding.vza]

model = model_from_parameters(parameters);
parameters.p   = oco_sounding.p_half
parameters.q   = oco_sounding.q
parameters.T   = oco_sounding.T# .+ 1.0 #.+ x[15]
parameters.sza = oco_sounding.sza
parameters.vza = [oco_sounding.vza]
model2 = model_from_parameters(parameters);

## Copied and adapted from rt_run
RS_type = CoreRT.noRS() 
results_canopy = []
results_NoCanopySameP = []
results_NoCanopyLowP = []
for Alb = 0.0:0.1:0.6
    for iBand = 1:3
        #iBand = 3
        lAlb = [0.384,0.211,0.061]
        # Canopy Stuff!
        LD = CanopyOptics.uniform_leaves()
        #LD = CanopyOptics.planophile_leaves2()
        leaf = LeafProspectProProperties{Float64}(Cw=0.006);
        LAI = 5.0
        ν = parameters.spec_bands[iBand]
        range = 1e7/ν[end]:1e7/ν[1];
        opti = createLeafOpticalStruct(range*u"nm");

        T,R = prospect(leaf,opti);
        T_ = mean(T)
        R_ = mean(R)
        ϖ_canopy = T_+R_
        @show T_+R_, iBand
        BiLambMod = CanopyOptics.BiLambertianCanopyScattering(R=R_,T=T_)
        # Run without canopy:
        model3 = deepcopy(model)
        model.params.brdf[iBand] = LambertianSurfaceScalar{Float64}(Alb)
        R_SFI_C, T_SFI_C, _, _       = rt_run_test(RS_type, model, BiLambMod, LAI, LD, ϖ_canopy, iBand)
        @show Alb
        model2.params.brdf[iBand]= LambertianSurfaceScalar{Float64}(lAlb[iBand])
        R_SFI_noC, T_SFI_noC, _, _   = CoreRT.rt_run_test(RS_type, model2, iBand)
        model3.params.brdf[iBand]= LambertianSurfaceScalar{Float64}(lAlb[iBand])
        R_SFI_noC2, T_SFI_noC2, _, _ = CoreRT.rt_run_test(RS_type, model3, iBand)
        
        push!(results_canopy,R_SFI_C[1,1,:])
        push!(results_NoCanopySameP, R_SFI_noC2[1,1,:])
        push!(results_NoCanopyLowP , R_SFI_noC[1,1,:])
    end
end
results_canopy_conv = []
results_NoCanopySameP_conv = []
results_NoCanopyLowP_conv = []
wl = []
for iBand in eachindex(results_canopy)
    @show iBand
    iBand2 = iBand % 3
    iBand2 == 0 ? iBand2 = 3 : nothing
    ils = oco_sounding.ils[iBand2]
    ν = parameters.spec_bands[iBand2]
    push!(results_canopy_conv, convOCO(results_canopy[iBand], ν, ils))
    push!(results_NoCanopySameP_conv, convOCO(results_NoCanopySameP[iBand], ν, ils))
    push!(results_NoCanopyLowP_conv, convOCO(results_NoCanopyLowP[iBand], ν, ils))
    push!(wl, oco_sounding.SpectralGrid[oco_sounding.BandID[iBand2]])
    @show iBand
end
albedos = 0.0:0.1:0.6
for i in eachindex(albedos)
    iAlb = i
    iBand = 1 + (iAlb-1)*3
    #plot(wl[iBand], results_canopy_conv[iBand], label="Canopy")
    #plot!(wl[iBand], results_NoCanopySameP_conv[iBand],label="No Canopy")
    r1 = results_NoCanopySameP_conv[iBand]./results_NoCanopyLowP_conv[iBand]
    r2 = results_canopy_conv[iBand]./results_NoCanopyLowP_conv[iBand]
    plot(wl[iBand]*1e3, r1./maximum(r1), label="NoC_highP/noC_lowP", lw=1.5, alpha=0.5)
    plot!(wl[iBand]*1e3, r2./maximum(r2), label="Integrated Canopy / No Canopy", lw=1.5, alpha=0.5)
    xlabel!("Wavelength (nm)")
    savefig("/home/cfranken/canopy/pdrdf_O2band"*string(albedos[i]*10) * ".pdf")
    savefig("/home/cfranken/canopy/pdrdf_O2band"*string(albedos[i]*10) * ".png")

    iBand = 2+ (iAlb-1)*3
    r1 = results_NoCanopySameP_conv[iBand]./results_NoCanopyLowP_conv[iBand]
    r2 = results_canopy_conv[iBand]./results_NoCanopyLowP_conv[iBand]
    plot(wl[iBand]*1e3, r1./maximum(r1), label="NoC_sameP/noC_lowP", lw=1.5, alpha=0.5)
    plot!(wl[iBand]*1e3, r2./maximum(r2), label="Integrated Canopy / No Canopy", lw=1.5, alpha=0.5)
    xlabel!("Wavelength (nm)")
    savefig("/home/cfranken/canopy/pdrdf_WCO2band"*string(albedos[i]*10) * ".pdf")
    savefig("/home/cfranken/canopy/pdrdf_WCO2band"*string(albedos[i]*10) * ".png")

    iBand = 3+ (iAlb-1)*3
    r1 = results_NoCanopySameP_conv[iBand]./results_NoCanopyLowP_conv[iBand]
    r2 = results_canopy_conv[iBand]./results_NoCanopyLowP_conv[iBand]
    plot(wl[iBand]*1e3, r1./maximum(r1), label="NoC_sameP/noC_lowP", lw=1.5, alpha=0.5)
    plot!(wl[iBand]*1e3, r2./maximum(r2), label="Integrated Canopy / No Canopy", lw=1.5, alpha=0.5)
    xlabel!("Wavelength (nm)")
    savefig("/home/cfranken/canopy/pdrdf_SCO2band."*string(albedos[i]*10) * ".pdf")
    savefig("/home/cfranken/canopy/pdrdf_SCO2band"*string(albedos[i]*10) * ".png")

end

iAlb = 2

using ColorSchemes

CS = ColorSchemes.seaborn_colorblind

iBand = 1 + (iAlb-1)*3
l = @layout [a  b c]
p1 = plot(wl[iBand]*1e3, results_canopy_conv[iBand],color=CS[1], label="Integrated Canopy ", lw=1.5, alpha=0.5)
ylims!(0,0.4)
p2 = plot(wl[iBand+1]*1e3, results_canopy_conv[iBand+1],color=CS[3], label="Integrated Canopy ", lw=1.5, alpha=0.5)
ylims!(0,0.4)
p3 = plot(wl[iBand+2]*1e3, results_canopy_conv[iBand+2],color=CS[4],label="Integrated Canopy ", lw=1.5, alpha=0.5)
ylims!(0,0.4)
plot!(size=(700,400))
plot(p1, p2, p3, layout = l)
savefig("/home/cfranken/canopy/AllBandsCanopyModel.pdf")

l = @layout [a  b c]

iBand = 1 + (iAlb-1)*3
    #plot(wl[iBand], results_canopy_conv[iBand], label="Canopy")
    #plot!(wl[iBand], results_NoCanopySameP_conv[iBand],label="No Canopy")
    r1 = results_NoCanopySameP_conv[iBand]./results_NoCanopyLowP_conv[iBand]
    r2 = results_canopy_conv[iBand]./results_NoCanopyLowP_conv[iBand]
    p1 = plot(wl[iBand]*1e3, r1./maximum(r1), label="NoC_highP/noC_lowP", lw=1.5, alpha=0.5)
    plot!(wl[iBand]*1e3, r2./maximum(r2), label="Integrated Canopy / No Canopy", lw=1.5, alpha=0.5)
    xlabel!("Wavelength (nm)")
    

    iBand = 2+ (iAlb-1)*3
    r1 = results_NoCanopySameP_conv[iBand]./results_NoCanopyLowP_conv[iBand]
    r2 = results_canopy_conv[iBand]./results_NoCanopyLowP_conv[iBand]
    p2 = plot(wl[iBand]*1e3, r1./maximum(r1), label="NoC_sameP/noC_lowP", lw=1.5, alpha=0.5)
    plot!(wl[iBand]*1e3, r2./maximum(r2), label="Integrated Canopy / No Canopy", lw=1.5, alpha=0.5)
    xlabel!("Wavelength (nm)")
    

    iBand = 3+ (iAlb-1)*3
    r1 = results_NoCanopySameP_conv[iBand]./results_NoCanopyLowP_conv[iBand]
    r2 = results_canopy_conv[iBand]./results_NoCanopyLowP_conv[iBand]
    p3 = plot(wl[iBand]*1e3, r1./maximum(r1), label="NoC_sameP/noC_lowP", lw=1.5, alpha=0.5)

    plot!(size=(700,400))
    plot(p1, p2, p3, layout = l)

plot(wl[iBand]*1e3, results_NoCanopySameP_conv[iBand], label="NoC_highP", lw=1.5, alpha=0.5)
plot!(wl[iBand]*1e3, results_canopy_conv[iBand], label="Integrated Canopy ", lw=1.5, alpha=0.5)
xlabel!("Wavelength (nm)")
# savefig("/home/cfranken/canopy/pdrdf_O2band"*string(albedos[i]*10) * ".pdf")
