##
using CUDA
device!(1)
using Revise
using vSmartMOM, vSmartMOM.CoreRT, vSmartMOM.SolarModel
using vSmartMOM.InelasticScattering
using Statistics
using Interpolations
using InstrumentOperator #for convolution of hires spectrum to instrument grid
using ImageFiltering
using Distributions
using Plots

##

# Load YAML files into parameter struct
parameters = 
    parameters_from_yaml("test/test_parameters/O2Parameters2.yaml");
#parameters = parameters_from_yaml("test/test_parameters/O2Parameters2.yaml");
# Create model struct (precomputes optical properties) from parameters
model      = model_from_parameters(parameters);

iBand = 1
FT = Float64

#### Compute all Raman properties
ν = model.params.spec_bands[iBand]
ν̃ = mean(ν);
# Find central reference index for RRS:
i_ref = argmin(abs.(ν .- ν̃))
# TODO_VS: λ_vs_in (get input)
# TODO_VS: ν_vs_in (convert to wavenumbers)
# Effective temperature for Raman calculations
effT = (model.profile.vcd_dry' * model.profile.T) / sum(model.profile.vcd_dry);
# Define RS type
# Compute N2 and O2

n2,o2 = InelasticScattering.getRamanAtmoConstants(ν̃,effT);
#greek_raman = get_greek_raman(RS_type, n2, o2);
RS_type = InelasticScattering.RRS(
            n2=n2,
            o2=o2,
            greek_raman = InelasticScattering.GreekCoefs([FT(1)], [FT(1)], [FT(1)], [FT(1)], [FT(1)], [FT(1)]),
            fscattRayl  = [FT(1)],
            ϖ_Cabannes  = [FT(1)], 
            ϖ_λ₁λ₀      = zeros(FT,1),
            i_λ₁λ₀      = zeros(Int,1), 
            Z⁻⁺_λ₁λ₀    = zeros(FT,1,1), 
            Z⁺⁺_λ₁λ₀    = zeros(FT,1,1), 
            i_ref       = i_ref,
            n_Raman     = 0,
            F₀          = zeros(FT,1,1));


#vSmartMOM.get_greek_raman!(RS_type, n2, o2);
# now compute other optical parameters for Raman:
#fscattRayl = ...
#RS_type = vSmartMOM.RRS(...)
# Compute Raman SSA properties:
CoreRT.getRamanSSProp!(RS_type, 1e7/ν̃, ν);
T_sun = 5777. # K
P = planck_spectrum_wn(T_sun, ν)
Tsolar = solar_transmission_from_file("/home/sanghavi/code/github/vSmartMOM.jl/src/SolarModel/solar.out")
Tsolar_interp = LinearInterpolation(Tsolar[4:end, 1], Tsolar[4:end, 2])
F₀ = zeros(length(P));
RS_type.F₀ = zeros(model.params.polarization_type.n, length(P))
for i=1:length(P)
    sol_trans = Tsolar_interp(ν[i]);
    F₀[i] = sol_trans * P[i];
    RS_type.F₀[1,i] = F₀[i]; #1.0 #
end 


R, T, ieR, ieT = CoreRT.rt_run_test(RS_type,model,iBand);
R_ss, T_ss, ieR_ss, ieT_ss = CoreRT.rt_run_test_ss(RS_type,model,iBand);

RS_type = InelasticScattering.noRS(
    fscattRayl  = [FT(1)],
    ϖ_Cabannes  = [FT(1)], 
    bandSpecLim = [],
    iBand       = [1],
    F₀          = zeros(FT,1,1));
RS_type.F₀ = zeros(model.params.polarization_type.n, length(P))
for i=1:length(P)
    #sol_trans = Tsolar_interp(ν[i]);
    #F₀[i] = sol_trans * P[i];
    RS_type.F₀[1,i] = F₀[i]; #F₀
end 

RnoRS, TnoRS, _, _ = CoreRT.rt_run_test(RS_type,model,iBand);
RnoRS_ss, TnoRS_ss, _, _ = CoreRT.rt_run_test_ss(RS_type,model,iBand);

#=
R_test, T_test, ieR_test, ieT_test = CoreRT.rt_run_test(RS_type,model,iBand);

# You can now run multiple bands like this (list at the end of band!)
RnoRS_test, TnoRS_test, _, _ = CoreRT.rt_run_test(vSmartMOM.noRS(),model,[1,1]);


R_test, T_test, ieR_test, ieT_test = CoreRT.rt_run_test(RS_type,model,1);
=#

#===Convolution of hires spectral simulations to instrument grid===#
#x = (1e7/415):0.3:(1e7/385)
#ν = (1e7/460):0.3:(1e7/420)
x = -40:0.5:40
#kernel = InstrumentOperator.create_instrument_kernel(Normal(0, 12.5), x) #defining a Gaussian kernel for convulution in wavenumber space
#kernel = InstrumentOperator.create_instrument_kernel(Normal(0, 12.5), x)
kernel = InstrumentOperator.create_instrument_kernel(Normal(0, 3.75), x) #SCIAMACHY resolution of FWHM=0.54 nm in the O2 A-band 

I_conv_noRS = zeros(3,length(ν))
I_conv = zeros(3,length(ν))
ieI_conv = zeros(3,length(ν))
Q_conv_noRS = zeros(3,length(ν))
Q_conv = zeros(3,length(ν))
ieQ_conv = zeros(3,length(ν))

I_conv_noRS_ss = zeros(3,length(ν))
I_conv_ss = zeros(3,length(ν))
ieI_conv_ss = zeros(3,length(ν))
Q_conv_noRS_ss = zeros(3,length(ν))
Q_conv_ss = zeros(3,length(ν))
ieQ_conv_ss = zeros(3,length(ν))

for ctr=1:3
    I_conv_noRS[ctr,:] .= imfilter(RnoRS[ctr,1,:], kernel)
    I_conv[ctr,:] .= imfilter(R[ctr,1,:], kernel)
    ieI_conv[ctr,:] .= imfilter(ieR[ctr,1,:], kernel)
    I_conv_noRS_ss[ctr,:] .= imfilter(RnoRS_ss[ctr,1,:], kernel)
    I_conv_ss[ctr,:] .= imfilter(R_ss[ctr,1,:], kernel)
    ieI_conv_ss[ctr,:] .= imfilter(ieR_ss[ctr,1,:], kernel)
    
    Q_conv_noRS[ctr,:] .= imfilter(RnoRS[ctr,2,:], kernel)
    Q_conv[ctr,:] .= imfilter(R[ctr,2,:], kernel)
    ieQ_conv[ctr,:] .= imfilter(ieR[ctr,2,:], kernel)
    Q_conv_noRS_ss[ctr,:] .= imfilter(RnoRS_ss[ctr,2,:], kernel)
    Q_conv_ss[ctr,:] .= imfilter(R_ss[ctr,2,:], kernel)
    ieQ_conv_ss[ctr,:] .= imfilter(ieR_ss[ctr,2,:], kernel)
end
#I_conv = InstrumentOperator.conv_spectra(kernel, )

convfct = 1e7./ν.^2   # to convert radiance units from mW/m²-str-cm⁻¹ to mW/m²-str-nm
#convfct = 1.0#./F₀ # to normalize wrt irradiance 
l = @layout [a1 a2 ; b1 b2; c1 c2]
#p1 = plot(1e7./ν, (R[1,1,:].+ieR[1,1,:]).*convfct, linecolor=:black, xlims=(755,775))
#p1 = plot!(1e7./ν, (R[2,1,:].+ieR[2,1,:]).*convfct, linecolor=:black, xlims=(755,775))
#p1 = plot!(1e7./ν, (R[3,1,:].+ieR[3,1,:]).*convfct, linecolor=:black, xlims=(755,775))
#p1 = plot!(1e7./ν, (I_conv[1,:].+ieI_conv[1,:]).*convfct, linewidth = 2, linecolor=:red, xlims=(755,775))
#p1 = plot!(1e7./ν, (I_conv[2,:].+ieI_conv[2,:]).*convfct, linewidth = 2, linecolor=:blue, xlims=(755,775))
#p1 = plot!(1e7./ν, (I_conv[3,:].+ieI_conv[3,:]).*convfct, linewidth = 2, linecolor=:green, xlims=(755,775))
p1 = plot(1e7./ν, (RnoRS[1,1,:]).*convfct, linecolor=:red, linewidth=0.5, alpha=0.5, xlims=(755,775))
p1 = plot!(1e7./ν, (RnoRS[2,1,:]).*convfct, linecolor=:blue, linewidth=0.5, alpha=0.5, xlims=(755,775))
p1 = plot!(1e7./ν, (RnoRS[3,1,:]).*convfct, linecolor=:green, linewidth=0.5, alpha=0.5, xlims=(755,775))
p1 = plot!(1e7./ν, (I_conv_noRS[1,:]).*convfct, linewidth = 2, linecolor=:red, xlims=(755,775))
p1 = plot!(1e7./ν, (I_conv_noRS[2,:]).*convfct, linewidth = 2, linecolor=:blue, xlims=(755,775))
p1 = plot!(1e7./ν, (I_conv_noRS[3,:]).*convfct, linewidth = 2, linecolor=:green, xlims=(755,775), ylabel="[mW/m²/str/nm]", yguidefontsize=8)


p2 = plot(1e7./ν, (R[1,1,:].-RnoRS[1,1,:].+ieR[1,1,:]).*convfct, linecolor=:red, linewidth=0.5, alpha=0.5, xlims=(755,775))#, ylims=(-2e-6, 1.e-5))#, xlabel = "λ [nm]")
p2 = plot!(1e7./ν, (R[2,1,:].-RnoRS[2,1,:].+ieR[2,1,:]).*convfct, linecolor=:blue, linewidth=0.5, alpha=0.5, xlims=(755,775))#, xlabel = "λ [nm]")
p2 = plot!(1e7./ν, (R[3,1,:].-RnoRS[3,1,:].+ieR[3,1,:]).*convfct, linecolor=:green, linewidth=0.5, alpha=0.5, xlims=(755,775))#, xlabel = "λ [nm]", xlims=(755,775))
p2 = plot!(1e7./ν, (I_conv[1,:].-I_conv_noRS[1,:].+ieI_conv[1,:]).*convfct, linewidth = 2, linecolor=:red, xlims=(755,775))
p2 = plot!(1e7./ν, (I_conv[2,:].-I_conv_noRS[2,:].+ieI_conv[2,:]).*convfct, linewidth = 2, linecolor=:blue, xlims=(755,775))
p2 = plot!(1e7./ν, (I_conv[3,:].-I_conv_noRS[3,:].+ieI_conv[3,:]).*convfct, linewidth = 2, linecolor=:green, xlims=(755,775), ylims=(-1,2.5), ylabel="[mW/m²/str/nm]", yguidefontsize=8)

p3 = plot(1e7./ν, 100*(R[1,1,:].-RnoRS[1,1,:].+ieR[1,1,:])./(RnoRS[1,1,:]), linecolor=:red, linewidth=0.5, alpha=0.5, xlabel = "λ [nm]", xlims=(755,775))
p3 = plot!(1e7./ν, 100*(R[2,1,:].-RnoRS[2,1,:].+ieR[2,1,:])./(RnoRS[2,1,:]), linecolor=:blue, linewidth=0.5, alpha=0.5, xlabel = "λ [nm]", xlims=(755,775))
p3 = plot!(1e7./ν, 100*(R[3,1,:].-RnoRS[3,1,:].+ieR[3,1,:])./(RnoRS[3,1,:]), linecolor=:green, linewidth=0.5, alpha=0.5, xlabel = "λ [nm]", xlims=(755,775))
p3 = plot!(1e7./ν, 100*(I_conv[1,:].-I_conv_noRS[1,:].+ieI_conv[1,:])./(I_conv_noRS[1,:]), linewidth = 2, linecolor=:red, xlims=(755,775))
p3 = plot!(1e7./ν, 100*(I_conv[2,:].-I_conv_noRS[2,:].+ieI_conv[2,:])./(I_conv_noRS[2,:]), linewidth = 2, linecolor=:blue, xlims=(755,775))
p3 = plot!(1e7./ν, 100*(I_conv[3,:].-I_conv_noRS[3,:].+ieI_conv[3,:])./(I_conv_noRS[3,:]), linewidth = 2, linecolor=:green, xlims=(755,775), ylims=(-1,5))

q1 = plot(1e7./ν, (RnoRS[1,2,:]).*convfct, linecolor=:red, linewidth=0.5, alpha=0.5, xlims=(755,775))
q1 = plot!(1e7./ν, (RnoRS[2,2,:]).*convfct, linecolor=:blue, linewidth=0.5, alpha=0.5, xlims=(755,775))
q1 = plot!(1e7./ν, (RnoRS[3,2,:]).*convfct, linecolor=:green, linewidth=0.5, alpha=0.5, xlims=(755,775))
q1 = plot!(1e7./ν, (Q_conv_noRS[1,:]).*convfct, linewidth = 2, linecolor=:red, xlims=(755,775))
q1 = plot!(1e7./ν, (Q_conv_noRS[2,:]).*convfct, linewidth = 2, linecolor=:blue, xlims=(755,775))
q1 = plot!(1e7./ν, (Q_conv_noRS[3,:]).*convfct, linewidth = 2, linecolor=:green, xlims=(755,775))#, ylabel="[mW/m²/str/nm]", yguidefontsize=8)

q2 = plot(1e7./ν, (R[1,2,:].-RnoRS[1,2,:].+ieR[1,2,:]).*convfct, linecolor=:red, linewidth=0.5, alpha=0.5, xlims=(755,775))
q2 = plot!(1e7./ν, (R[2,2,:].-RnoRS[2,2,:].+ieR[2,2,:]).*convfct, linecolor=:blue, linewidth=0.5, alpha=0.5, xlims=(755,775))
q2 = plot!(1e7./ν, (R[3,2,:].-RnoRS[3,2,:].+ieR[3,2,:]).*convfct, linecolor=:green, linewidth=0.5, alpha=0.5, xlims=(755,775))
q2 = plot!(1e7./ν, (Q_conv[1,:].-Q_conv_noRS[1,:].+ieQ_conv[1,:]).*convfct, linewidth = 2, linecolor=:red, xlims=(755,775))
q2 = plot!(1e7./ν, (Q_conv[2,:].-Q_conv_noRS[2,:].+ieQ_conv[2,:]).*convfct, linewidth = 2, linecolor=:blue, xlims=(755,775))
q2 = plot!(1e7./ν, (Q_conv[3,:].-Q_conv_noRS[3,:].+ieQ_conv[3,:]).*convfct, linewidth = 2, linecolor=:green, xlims=(755,775), ylims=(-0.1,0.18))#, ylabel="[mW/m²/str/nm]", yguidefontsize=8)

q3 = plot(1e7./ν, 100*(R[1,2,:].-RnoRS[1,2,:].+ieR[1,2,:])./(R[1,2,:].+ieR[1,2,:]), linecolor=:red, linewidth=0.5, alpha=0.5, xlabel = "λ [nm]", xlims=(755,775))
q3 = plot!(1e7./ν, 100*(R[2,2,:].-RnoRS[2,2,:].+ieR[2,2,:])./(R[2,2,:].+ieR[1,2,:]), linecolor=:blue, linewidth=0.5, alpha=0.5, xlabel = "λ [nm]", xlims=(755,775))
q3 = plot!(1e7./ν, 100*(R[3,2,:].-RnoRS[3,2,:].+ieR[3,2,:])./(R[3,2,:].+ieR[1,2,:]), linecolor=:green, linewidth=0.5, alpha=0.5, xlabel = "λ [nm]", xlims=(755,775))
q3 = plot!(1e7./ν, 100*(Q_conv[1,:].-Q_conv_noRS[1,:].+ieQ_conv[1,:])./(Q_conv_noRS[1,:]), linewidth = 3, linecolor=:red, xlims=(755,775))
q3 = plot!(1e7./ν, 100*(Q_conv[2,:].-Q_conv_noRS[2,:].+ieQ_conv[2,:])./(Q_conv_noRS[2,:]), linewidth = 3, linecolor=:green, xlims=(755,775))
q3 = plot!(1e7./ν, 100*(Q_conv[3,:].-Q_conv_noRS[3,:].+ieQ_conv[3,:])./(Q_conv_noRS[3,:]), linewidth = 3, linecolor=:blue, xlims=(755,775), ylims=(-0.1,0.5))

#Q_conv_noRS = imfilter(RnoRS[1,2,:], kernel)
#Q_conv = imfilter(R[1,2,:], kernel)
#ieQ_conv = imfilter(ieR[1,2,:], kernel)

#Q_conv_noRS = imfilter(RnoRS[2,2,:], kernel)
#Q_conv = imfilter(R[2,2,:], kernel)
#ieQ_conv = imfilter(ieR[2,2,:], kernel)

#Q_conv_noRS = imfilter(RnoRS[3,2,:], kernel)
#Q_conv = imfilter(R[3,2,:], kernel)
#ieQ_conv = imfilter(ieR[3,2,:], kernel)

plot(p1, q1, p2, q2, p3, q3, layout = l, legend = false, title = ["I₀ (no RS)" "Q₀ (no RS)" "I₁-I₀" "Q₁-Q₀" "(I₁-I₀)/I₀ [%]" "(Q₁-Q₀)/Q₀ [%]"], titlefont = font(10))
savefig("RingEffect_O2A_SZA30_wF.png")

l = @layout [a1 a2; b1 b2]
p1 = plot(1e7./ν, 100*(1.0.-(ieR_ss[1,1,:]+R_ss[1,1,:]-RnoRS_ss[1,1,:])./(ieR[1,1,:]+R[1,1,:]-RnoRS[1,1,:])), linecolor=:red, linewidth=0.5, alpha=0.5, xlims=(755,775))
p1 = plot!(1e7./ν, 100*(1.0.-(ieR_ss[2,1,:]+R_ss[2,1,:]-RnoRS_ss[2,1,:])./(ieR[2,1,:]+R[2,1,:]-RnoRS[2,1,:])), linecolor=:blue, linewidth=0.5, alpha=0.5, xlims=(755,775))
p1 = plot!(1e7./ν, 100*(1.0.-(ieR_ss[3,1,:]+R_ss[3,1,:]-RnoRS_ss[3,1,:])./(ieR[3,1,:]+R[3,1,:]-RnoRS[3,1,:])), linecolor=:green, linewidth=0.5, alpha=0.5, xlims=(755,775), ylims=(-100,100))
#p1 = plot!(1e7./ν, 100*(1.0.-(ieI_conv_ss[1,:]+I_conv_ss[1,:]-I_conv_noRS_ss[1,:])./(ieI_conv[1,:]+I_conv[1,:]-I_conv_noRS[1,:])), seriestype=:scatter, ms = 2, mc=:red, xlabel = "λ [nm]", xlims=(755,775))
#p1 = plot!(1e7./ν, 100*(1.0.-(ieI_conv_ss[2,:]+I_conv_ss[2,:]-I_conv_noRS_ss[2,:])./(ieI_conv[2,:]+I_conv[2,:]-I_conv_noRS[2,:])), seriestype=:scatter, ms = 2, mc=:blue, xlabel = "λ [nm]", xlims=(755,775))
#p1 = plot!(1e7./ν, 100*(1.0.-(ieI_conv_ss[3,:]+I_conv_ss[3,:]-I_conv_noRS_ss[3,:])./(ieI_conv[3,:]+I_conv[3,:]-I_conv_noRS[3,:])), seriestype=:scatter, ms = 2, mc=:green, xlabel = "λ [nm]", xlims=(755,775))
p2 = plot(1e7./ν, 100*(1.0.-(RnoRS_ss[1,1,:]./RnoRS[1,1,:])), linecolor=:red, linewidth=0.5, alpha=0.5, xlims=(755,775))
p2 = plot!(1e7./ν, 100*(1.0.-(RnoRS_ss[2,1,:]./RnoRS[2,1,:])), linecolor=:blue, linewidth=0.5, alpha=0.5, xlims=(755,775))
p2 = plot!(1e7./ν, 100*(1.0.-(RnoRS_ss[3,1,:]./RnoRS[3,1,:])), linecolor=:green, linewidth=0.5, alpha=0.5, xlims=(755,775), ylims=(0,5))


q1 = plot(1e7./ν, (ieR_ss[1,2,:]+R_ss[1,2,:]-RnoRS_ss[1,2,:])./(ieR[1,2,:]+R[1,2,:]-RnoRS[1,2,:]), linecolor=:red, linewidth=0.5, alpha=0.5, xlims=(755,775))
q1 = plot!(1e7./ν, (ieR_ss[2,2,:]+R_ss[2,2,:]-RnoRS_ss[2,2,:])./(ieR[2,2,:]+R[2,2,:]-RnoRS[2,2,:]), linecolor=:blue, linewidth=0.5, alpha=0.5, xlims=(755,775))
q1 = plot!(1e7./ν, (ieR_ss[3,2,:]+R_ss[3,2,:]-RnoRS_ss[3,2,:])./(ieR[3,2,:]+R[3,2,:]-RnoRS[3,2,:]), linecolor=:green, linewidth=0.5, alpha=0.5, xlims=(755,775), ylims=(-2,10))
#q1 = plot!(1e7./ν, (ieQ_conv_ss[1,:]+Q_conv_ss[1,:]-Q_conv_noRS_ss[1,:])./(ieQ_conv[1,:]+Q_conv[1,:]-Q_conv_noRS[1,:]), linewidth = 2, linecolor=:red, xlabel = "λ [nm]", xlims=(755,775))
#q1 = plot!(1e7./ν, (ieQ_conv_ss[2,:]+Q_conv_ss[2,:]-Q_conv_noRS_ss[2,:])./(ieQ_conv[2,:]+Q_conv[2,:]-Q_conv_noRS[2,:]), linewidth = 2, linecolor=:blue, xlabel = "λ [nm]", xlims=(755,775))
#q1 = plot!(1e7./ν, (ieQ_conv_ss[3,:]+Q_conv_ss[3,:]-Q_conv_noRS_ss[3,:])./(ieQ_conv[3,:]+Q_conv[3,:]-Q_conv_noRS[3,:]), linewidth = 2, linecolor=:green, xlabel = "λ [nm]", xlims=(755,775))

q2 = plot(1e7./ν, 100*(1.0.-(RnoRS_ss[1,2,:]./RnoRS[1,2,:])), linecolor=:red, linewidth=0.5, alpha=0.5, xlims=(755,775))
q2 = plot!(1e7./ν, 100*(1.0.-(RnoRS_ss[2,2,:]./RnoRS[2,2,:])), linecolor=:blue, linewidth=0.5, alpha=0.5, xlims=(755,775))
q2 = plot!(1e7./ν, 100*(1.0.-(RnoRS_ss[3,2,:]./RnoRS[3,2,:])), linecolor=:green, linewidth=0.5, alpha=0.5, xlims=(755,775), ylims=(0,3))


#q1 = plot(1e7./ν, 100*ieR[1,2,:]./R[1,2,:], linecolor=:black, xlims=(755,775))
#q1 = plot!(1e7./ν, 100*ieR[2,2,:]./R[2,2,:], linecolor=:black, xlims=(755,775))
#q1 = plot!(1e7./ν, 100*ieR[3,2,:]./R[3,2,:], linecolor=:black, xlims=(755,775))
#q1 = plot!(1e7./ν, 100*ieQ_conv[1,:]./Q_conv[1,:], linewidth = 2, linecolor=:red, xlabel = "λ [nm]", xlims=(755,775))
#q1 = plot!(1e7./ν, 100*ieQ_conv[2,:]./Q_conv[2,:], linewidth = 2, linecolor=:blue, xlabel = "λ [nm]", xlims=(755,775))
#q1 = plot!(1e7./ν, 100*ieQ_conv[3,:]./Q_conv[3,:], linewidth = 2, linecolor=:green, xlabel = "λ [nm]", xlims=(755,775))

plot(p2, p1, q2, q1, layout = l, legend = false, title = ["(Iₑ-I'ₑ)/Iₑ [%]" "(Iᵢ-Iₑ)-(I'ᵢ-I'ₑ)/(Iᵢ-Iₑ) [%]"  "(Qₑ-Q'ₑ)/Qₑ [%]" "(Qᵢ-Qₑ)-(Q'ᵢ-Q'ₑ)/(Qᵢ-Qₑ) [%]"], titlefont = font(10))
#savefig("RingSpectrum_O2A_SZA30_wF.png")


l = @layout [a1 a2 ; b1 b2; c1 c2]
#p1 = plot(1e7./ν, (R[1,1,:].+ieR[1,1,:]).*convfct, linecolor=:black, xlims=(755,775))
#p1 = plot!(1e7./ν, (R[2,1,:].+ieR[2,1,:]).*convfct, linecolor=:black, xlims=(755,775))
#p1 = plot!(1e7./ν, (R[3,1,:].+ieR[3,1,:]).*convfct, linecolor=:black, xlims=(755,775))
#p1 = plot!(1e7./ν, (I_conv[1,:].+ieI_conv[1,:]).*convfct, linewidth = 2, linecolor=:red, xlims=(755,775))
#p1 = plot!(1e7./ν, (I_conv[2,:].+ieI_conv[2,:]).*convfct, linewidth = 2, linecolor=:blue, xlims=(755,775))
#p1 = plot!(1e7./ν, (I_conv[3,:].+ieI_conv[3,:]).*convfct, linewidth = 2, linecolor=:green, xlims=(755,775))

p1 = plot(1e7./ν, (R[1,1,:].-RnoRS[1,1,:].+ieR[1,1,:]).*convfct, linecolor=:black, xlims=(755,775))#, ylims=(-2e-6, 1.e-5))#, xlabel = "λ [nm]")
p1 = plot!(1e7./ν, (R_ss[1,1,:].-RnoRS_ss[1,1,:].+ieR_ss[1,1,:]).*convfct, linecolor=:red, xlims=(755,775))#, ylims=(-2e-6, 1.e-5))#, xlabel = "λ [nm]")

p2 = plot(1e7./ν, (R[2,1,:].-RnoRS[2,1,:].+ieR[2,1,:]).*convfct, linecolor=:black, xlims=(755,775))#, xlabel = "λ [nm]")
p2 = plot!(1e7./ν, (R_ss[2,1,:].-RnoRS_ss[2,1,:].+ieR_ss[2,1,:]).*convfct, linecolor=:red, xlims=(755,775))#, xlabel = "λ [nm]")

p3 = plot(1e7./ν, (R[3,1,:].-RnoRS[3,1,:].+ieR[3,1,:]).*convfct, linecolor=:black, xlims=(755,775))#, xlabel = "λ [nm]", xlims=(755,775))
p3 = plot!(1e7./ν, (R_ss[3,1,:].-RnoRS_ss[3,1,:].+ieR_ss[3,1,:]).*convfct, linecolor=:red, xlims=(755,775))#, xlabel = "λ [nm]", xlims=(755,775))

q1 = plot(1e7./ν, (R[1,2,:].-RnoRS[1,2,:].+ieR[1,2,:]).*convfct, linecolor=:black, xlims=(755,775))#, ylims=(-2e-6, 1.e-5))#, xlabel = "λ [nm]")
q1 = plot!(1e7./ν, (R_ss[1,2,:].-RnoRS_ss[1,2,:].+ieR_ss[1,2,:]).*convfct, linecolor=:red, xlims=(755,775))#, ylims=(-2e-6, 1.e-5))#, xlabel = "λ [nm]")

q2 = plot(1e7./ν, (R[2,2,:].-RnoRS[2,2,:].+ieR[2,2,:]).*convfct, linecolor=:black, xlims=(755,775))#, xlabel = "λ [nm]")
q2 = plot!(1e7./ν, (R_ss[2,2,:].-RnoRS_ss[2,2,:].+ieR_ss[2,2,:]).*convfct, linecolor=:red, xlims=(755,775))#, xlabel = "λ [nm]")

q3 = plot(1e7./ν, (R[3,2,:].-RnoRS[3,2,:].+ieR[3,2,:]).*convfct, linecolor=:black, xlims=(755,775))#, xlabel = "λ [nm]", xlims=(755,775))
q3 = plot!(1e7./ν, (R_ss[3,2,:].-RnoRS_ss[3,2,:].+ieR_ss[3,2,:]).*convfct, linecolor=:red, xlims=(755,775))#, xlabel = "λ [nm]", xlims=(755,775))


#p2 = plot!(1e7./ν, (I_conv[1,:].-I_conv_noRS[1,:].+ieI_conv[1,:]).*convfct, linewidth = 2, linecolor=:red, xlims=(755,775))
#p2 = plot!(1e7./ν, (I_conv[2,:].-I_conv_noRS[2,:].+ieI_conv[2,:]).*convfct, linewidth = 2, linecolor=:blue, xlims=(755,775))
#p2 = plot!(1e7./ν, (I_conv[3,:].-I_conv_noRS[3,:].+ieI_conv[3,:]).*convfct, linewidth = 2, linecolor=:green, xlims=(755,775))

#p3 = plot(1e7./ν, 100*(R[1,1,:].-RnoRS[1,1,:].+ieR[1,1,:])./(R[1,1,:].+ieR[1,1,:]), linecolor=:black, xlabel = "λ [nm]", xlims=(755,775))
#p3 = plot!(1e7./ν, 100*(R[2,1,:].-RnoRS[2,1,:].+ieR[2,1,:])./(R[2,1,:].+ieR[2,1,:]), linecolor=:black, xlabel = "λ [nm]", xlims=(755,775))
#p3 = plot!(1e7./ν, 100*(R[3,1,:].-RnoRS[3,1,:].+ieR[3,1,:])./(R[3,1,:].+ieR[3,1,:]), linecolor=:black, xlabel = "λ [nm]", xlims=(755,775))
#p3 = plot!(1e7./ν, 100*(I_conv[1,:].-I_conv_noRS[1,:].+ieI_conv[1,:])./(I_conv[1,:].+ieI_conv[1,:]), linewidth = 2, linecolor=:red, xlims=(755,775))
#p3 = plot!(1e7./ν, 100*(I_conv[2,:].-I_conv_noRS[2,:].+ieI_conv[2,:])./(I_conv[2,:].+ieI_conv[2,:]), linewidth = 2, linecolor=:blue, xlims=(755,775))
#p3 = plot!(1e7./ν, 100*(I_conv[3,:].-I_conv_noRS[3,:].+ieI_conv[3,:])./(I_conv[3,:].+ieI_conv[3,:]), linewidth = 2, linecolor=:green, xlims=(755,775))


#q1 = plot(1e7./ν, (R[1,2,:].+ieR[1,2,:]).*convfct, linecolor=:black, xlims=(755,775))
#q1 = plot!(1e7./ν, (R[2,2,:].+ieR[2,2,:]).*convfct, linecolor=:black, xlims=(755,775))
#q1 = plot!(1e7./ν, (R[3,2,:].+ieR[3,2,:]).*convfct, linecolor=:black, xlims=(755,775))
#q1 = plot!(1e7./ν, (Q_conv[1,:].+ieQ_conv[1,:]).*convfct, linewidth = 2, linecolor=:red, xlims=(755,775))
#q1 = plot!(1e7./ν, (Q_conv[2,:].+ieQ_conv[2,:]).*convfct, linewidth = 2, linecolor=:blue, xlims=(755,775))
#q1 = plot!(1e7./ν, (Q_conv[3,:].+ieQ_conv[3,:]).*convfct, linewidth = 2, linecolor=:green, xlims=(755,775))

#q2 = plot(1e7./ν, (R[1,2,:].-RnoRS[1,2,:].+ieR[1,2,:]).*convfct, linecolor=:black, xlims=(755,775))
#q2 = plot!(1e7./ν, (R[2,2,:].-RnoRS[2,2,:].+ieR[2,2,:]).*convfct, linecolor=:black, xlims=(755,775))
#q2 = plot!(1e7./ν, (R[3,2,:].-RnoRS[3,2,:].+ieR[3,2,:]).*convfct, linecolor=:black, xlims=(755,775))
#q2 = plot!(1e7./ν, (Q_conv[1,:].-Q_conv_noRS[1,:].+ieQ_conv[1,:]).*convfct, linewidth = 2, linecolor=:red, xlims=(755,775))
#q2 = plot!(1e7./ν, (Q_conv[2,:].-Q_conv_noRS[2,:].+ieQ_conv[2,:]).*convfct, linewidth = 2, linecolor=:blue, xlims=(755,775))
#q2 = plot!(1e7./ν, (Q_conv[3,:].-Q_conv_noRS[3,:].+ieQ_conv[3,:]).*convfct, linewidth = 2, linecolor=:green, xlims=(755,775))

#q3 = plot(1e7./ν, 100*(R[1,2,:].-RnoRS[1,2,:].+ieR[1,2,:])./(R[1,2,:].+ieR[1,2,:]), linecolor=:black, xlabel = "λ [nm]", xlims=(755,775))
#q3 = plot!(1e7./ν, 100*(R[2,2,:].-RnoRS[2,2,:].+ieR[2,2,:])./(R[2,2,:].+ieR[1,2,:]), linecolor=:black, xlabel = "λ [nm]", xlims=(755,775))
#q3 = plot!(1e7./ν, 100*(R[3,2,:].-RnoRS[3,2,:].+ieR[3,2,:])./(R[3,2,:].+ieR[1,2,:]), linecolor=:black, xlabel = "λ [nm]", xlims=(755,775))
#q3 = plot!(1e7./ν, 100*(Q_conv[1,:].-Q_conv_noRS[1,:].+ieQ_conv[1,:])./(Q_conv[1,:].+ieQ_conv[1,:]), linewidth = 3, linecolor=:red, xlims=(755,775))
#q3 = plot!(1e7./ν, 100*(Q_conv[2,:].-Q_conv_noRS[2,:].+ieQ_conv[2,:])./(Q_conv[2,:].+ieQ_conv[2,:]), linewidth = 3, linecolor=:green, xlims=(755,775))
#q3 = plot!(1e7./ν, 100*(Q_conv[3,:].-Q_conv_noRS[3,:].+ieQ_conv[3,:])./(Q_conv[3,:].+ieQ_conv[3,:]), linewidth = 3, linecolor=:blue, xlims=(755,775))

#Q_conv_noRS = imfilter(RnoRS[1,2,:], kernel)
#Q_conv = imfilter(R[1,2,:], kernel)
#ieQ_conv = imfilter(ieR[1,2,:], kernel)

#Q_conv_noRS = imfilter(RnoRS[2,2,:], kernel)
#Q_conv = imfilter(R[2,2,:], kernel)
#ieQ_conv = imfilter(ieR[2,2,:], kernel)

#Q_conv_noRS = imfilter(RnoRS[3,2,:], kernel)
#Q_conv = imfilter(R[3,2,:], kernel)
#ieQ_conv = imfilter(ieR[3,2,:], kernel)

plot(p1, q1, p2, q2, p3, q3, layout = l, legend = false, title = ["I₁ (with RS)" "Q₁ (with RS)" "I₁-I₀" "Q₁-Q₀" "(1-I₀/I₁) [%]" "(1-Q₀/Q₁) [%]"], titlefont = font(10))
savefig("RingEffect_O2A_SZA30_wF.png")
