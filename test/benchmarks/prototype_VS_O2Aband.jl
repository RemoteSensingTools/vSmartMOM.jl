##
using Revise
using vSmartMOM, vSmartMOM.CoreRT, vSmartMOM.SolarModel
using vSmartMOM.InelasticScattering
using Statistics
using Interpolations
using InstrumentOperator #for convolution of hires spectrum to instrument grid
using ImageFiltering
using Distributions
using Plots

Fraunhofer=true
Tsolar = solar_transmission_from_file("/home/sanghavi/code/github/vSmartMOM.jl/src/SolarModel/solar.out")
Tsolar_interp = LinearInterpolation(Tsolar[4:end, 1], Tsolar[4:end, 2])
        
ν₀ = (1e7/730.):0.5:(1e7/630.) #Spectral range of sources of VS  
#ν₀ = (1e7/750.):0.5:(1e7/731.5) #Spectral range of sources of VS  
λ₀ = 2e7/(ν₀[1]+ν₀[end]) #400. #nm
ii₀ = findmin(abs.(1e7./ν₀ .- λ₀))[2]
νₐ = (1e7/780.):0.5:(1e7/750.) #Target spectral range (O2 A-band)  
# Computing Fraunhofer Spectrum 
T_sun = 5777. # KRS 
P₀ = planck_spectrum_wn(T_sun, collect(ν₀)) 
Pₐ = planck_spectrum_wn(T_sun, collect(νₐ))    
#F₀ = zeros(length(ν₀));
#ν₀ = 1e7/450.:2.:1e7/375.
R₀ = zeros(length(ν₀)) 
Rₐ = zeros(length(νₐ)) 
#T₀ = zeros(length(ν₀)) 
ieRₐ = zeros(length(νₐ)) 
#ieT₀ = zeros(length(ν₀))
R₀_Q = zeros(length(ν₀)) 
Rₐ_Q = zeros(length(νₐ)) 
#T₀_Q = zeros(length(ν₀)) 
ieRₐ_Q = zeros(length(νₐ)) 
#ieT₀_Q = zeros(length(ν₀))
#R₀_U = zeros(3,length(ν₀)) 
#T₀_U = zeros(length(ν₀)) 
#ieR₀_U = zeros(3,length(ν₀)) 
#ieT₀_U = zeros(length(ν₀))
I₀_conv = zeros(length(ν₀)) 
Iₐ_conv = zeros(length(νₐ)) 
ieIₐ_conv = zeros(length(νₐ)) 
Q₀_conv = zeros(length(ν₀)) 
Qₐ_conv = zeros(length(νₐ)) 
ieQₐ_conv = zeros(length(νₐ)) 

#=
for i=1:length(ν₀)
    #    ##sol_trans = Tsolar_interp(ν₀[i]);
    #    F₀[i] = 1.0##sol_trans * P[i];
    #end 
    λ₀ = 1e7/ν₀[i]        
    FT = Float64
    #λ₀ = 1e7/ν₀[ctr]
    #iBand = 1
    n2,o2 = InelasticScattering.getRamanAtmoConstants(1.7/λ₀, 300.);
    RS_type = InelasticScattering.VS_0to1_plus{FT}(
                n2=n2,
                o2=o2);
    # Load YAML files into parameter struct
    parameters = parameters_from_yaml("test/test_parameters/O2ParametersVS.yaml");
    # Create model struct (precomputes optical properties) from parameters
    model      = model_from_parameters(RS_type, λ₀, parameters, νₐ);
    @show RS_type.bandSpecLim
    Fraunhofer=true    
    RS_type.F₀ = zeros(model.params.polarization_type.n, 
        sum(length(RS_type.grid_in[j]) for j in 1:length(RS_type.iBand)))
    global t_offset = 0
    for iB=1:length(RS_type.iBand)
        ν = collect(RS_type.grid_in[iB])
        #global t_offset
        for iii=1:length(ν)
            sol_trans = Tsolar_interp(ν[iii]);
            #    F₀[i] = 1.0##sol_trans * P[i];
            soltmp = 1.
            if (iB==1)
                soltmp = sol_trans * P₀[i];         
            else
                soltmp = sol_trans * Pₐ[iii]
            end
            #@show iii,  t_offset, soltmp
            RS_type.F₀[1,iii+t_offset] = soltmp; #sol_trans * P[i];#1.0;
        end 
        t_offset += length(ν)
    end

    R, T, ieR, ieT = CoreRT.rt_run_test(RS_type,
            model, RS_type.iBand);

    R₀[i] = R[1,1,1]
    R₀_Q[i] = R[1,2,1]
    Rₐ[:] = R[1,1,2:end]
    Rₐ_Q[:] = R[1,2,2:end]
    ieRₐ[:] += ieR[1,1,2:end]
    ieRₐ_Q[:] += ieR[1,2,2:end]
    @show i, size(ν₀) 
end
=#
using DelimitedFiles
vsout = readdlm("out_VS_O2Aband.dat")
srcout = readdlm("out_src_630nm_730nm.dat")
ν₀  =srcout[:,1]
R₀  =srcout[:,2]
R₀_Q=srcout[:,3]
#O2 A-band
νₐ=vsout[:,1]
Rₐ    = vsout[:,2]
Rₐ_Q  = vsout[:,3]
ieRₐ  = vsout[:,4]
ieRₐ_Q= vsout[:,5]

# =========Computing RRS on O2 A-band  ======== #
Iₐ_RRS_conv   = zeros(length(νₐ)) 
Iₐ_RRS_ss_conv= zeros(length(νₐ)) 
ieIₐ_RRS_conv = zeros(length(νₐ)) 
ieIₐ_RRS_ss_conv = zeros(length(νₐ))
Iₐ_noRS_conv  = zeros(length(νₐ))
Qₐ_RRS_conv   = zeros(length(νₐ)) 
Qₐ_RRS_ss_conv= zeros(length(νₐ)) 
ieQₐ_RRS_conv = zeros(length(νₐ)) 
ieQₐ_RRS_ss_conv = zeros(length(νₐ)) 
Qₐ_noRS_conv  = zeros(length(νₐ))

parameters = 
    parameters_from_yaml("test/test_parameters/O2ParametersVS.yaml");
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
#Tsolar = solar_transmission_from_file("/home/sanghavi/code/github/vSmartMOM.jl/src/SolarModel/solar.out")
#Tsolar_interp = LinearInterpolation(Tsolar[4:end, 1], Tsolar[4:end, 2])
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
# ============================================= #

x = -40.:1.:40.
#kernel = InstrumentOperator.create_instrument_kernel(Normal(0, 12.5), x) #defining a Gaussian kernel for convulution in wavenumber space
kernel = InstrumentOperator.create_instrument_kernel(Normal(0, 12.5), x)
#I_conv = InstrumentOperator.conv_spectra(kernel, )
    
I₀_conv .= imfilter(R₀, kernel)
Iₐ_conv .= imfilter(Rₐ, kernel)
ieIₐ_conv .= imfilter(ieRₐ, kernel)
Iₐ_RRS_conv .= imfilter(R[1,1,:], kernel)
Iₐ_RRS_ss_conv .= imfilter(R_ss[1,1,:], kernel)
ieIₐ_RRS_conv .= imfilter(ieR[1,1,:], kernel)
ieIₐ_RRS_ss_conv .= imfilter(ieR_ss[1,1,:], kernel)
Iₐ_noRS_conv .= imfilter(RnoRS[1,1,:], kernel)
Q₀_conv .= imfilter(R₀_Q, kernel)
Qₐ_conv .= imfilter(Rₐ_Q, kernel)
ieQₐ_conv .= imfilter(ieRₐ_Q, kernel)    
Qₐ_RRS_conv .= imfilter(R[1,2,:], kernel)
Qₐ_RRS_ss_conv .= imfilter(R_ss[1,2,:], kernel)
ieQₐ_RRS_conv .= imfilter(ieR[1,2,:], kernel)
ieQₐ_RRS_ss_conv .= imfilter(ieR_ss[1,2,:], kernel)
Qₐ_noRS_conv .= imfilter(RnoRS[1,2,:], kernel)
#end


convfct0 = 1e7./ν₀.^2   # to convert radiance units from mW/m²-str-cm⁻¹ to mW/m²-str-nm
convfcta = 1e7./νₐ.^2 

l = @layout [a1 a2; b1 b2; c1 c2; d1 d2]
ymax = round(maximum(R₀.*convfct0), sigdigits=2)
p1 = plot(1e7./ν₀, R₀.*convfct0, linecolor=:grey, yticks=0:(ymax/4):ymax)
p1 = plot!(1e7./ν₀, I₀_conv.*convfct0, linecolor=:magenta)#, xticks=385:5:415, ylabel="[mW/m²/str/nm]")
ymax = round(maximum(R₀_Q.*convfct0), sigdigits=1)
q1 = plot(1e7./ν₀, R₀_Q.*convfct0, linecolor=:grey, yticks=0:(ymax/4):ymax)#, xticks=385:5:415) #, xticks=425:3:430)
q1 = plot!(1e7./ν₀, Q₀_conv.*convfct0, linecolor=:magenta)#, xticks=385:5:415) #, xticks=425:3:430)

ymax = round(maximum(Rₐ.*convfcta), sigdigits=2)
p2 = plot(1e7./νₐ, Rₐ.*convfcta, linecolor=:black, yticks=0:(ymax/4):ymax, xlims=(755,775))#)#, xlims=(385,415), xticks=385:5:415, ylabel="[mW/m²/str/nm]")
p2 = plot!(1e7./νₐ, Iₐ_conv.*convfcta, linecolor=:red)#)#, xlims=(385,415), xticks=385:5:415, ylabel="[mW/m²/str/nm]")
ymax = round(maximum(Rₐ_Q.*convfcta), sigdigits=1)
q2 = plot(1e7./νₐ, Rₐ_Q.*convfcta, linecolor=:black, yticks=0:(ymax/4):ymax, xlims=(755,775))#)#, xlims=(385,415), xticks=385:5:415) #, xticks=425:3:430)
q2 = plot!(1e7./νₐ, Qₐ_conv.*convfcta, linecolor=:red)#)#, xlims=(385,415), xticks=385:5:415) #, xticks=425:3:430)

ymax = round(maximum(ieRₐ.*convfcta), sigdigits=2)
p3 = plot(1e7./νₐ, ieRₐ.*convfcta, linecolor=:black, yticks=0:(ymax/4):ymax, xlims=(755,775))#)#, xlims=(385,415), xticks=385:5:415, ylabel="[mW/m²/str/nm]")
p3 = plot!(1e7./νₐ, ieIₐ_conv.*convfcta, linecolor=:red)#)#, xlims=(385,415), xticks=385:5:415, ylabel="[mW/m²/str/nm]")
ymax = round(maximum(ieRₐ_Q.*convfcta), sigdigits=1)
q3 = plot(1e7./νₐ, ieRₐ_Q.*convfcta, linecolor=:black, yticks=0:(ymax/4):ymax, xlims=(755,775))#)#, xlims=(385,415), xticks=385:5:415) #, xticks=425:3:430)
q3 = plot!(1e7./νₐ, ieQₐ_conv.*convfcta, linecolor=:red)#)#, xlims=(385,415), xticks=385:5:415) #, xticks=425:3:430)

ymax = round(maximum((ieR[1,1,:].+R[1,1,:].-RnoRS[1,1,:]).*convfcta), sigdigits=2)
ymin = round(minimum((ieR[1,1,:].+R[1,1,:].-RnoRS[1,1,:]).*convfcta), sigdigits=2)
p4 = plot(1e7./νₐ, (ieR[1,1,:].+R[1,1,:].-RnoRS[1,1,:]).*convfcta, linecolor=:black, yticks=0:((ymax-ymin)/4):ymax, xlims=(755,775), ylims=(0,ymax))#)#, xlims=(385,415), xticks=385:5:415, ylabel="[mW/m²/str/nm]")
p4 = plot!(1e7./νₐ, (ieIₐ_RRS_conv.+Iₐ_RRS_conv.-Iₐ_noRS_conv).*convfcta, linecolor=:red)#)#, xlims=(385,415), xticks=385:5:415, ylabel="[mW/m²/str/nm]")
ymax = round(maximum((ieR[1,2,:].+R[1,2,:].-RnoRS[1,2,:]).*convfcta), sigdigits=2)
ymin = round(minimum((ieR[1,2,:].+R[1,2,:].-RnoRS[1,2,:]).*convfcta), sigdigits=2)
q4 = plot(1e7./νₐ, (ieR[1,2,:].+R[1,2,:].-RnoRS[1,2,:]).*convfcta, linecolor=:black, yticks=ymin:((ymax-ymin)/4):ymax, xlims=(755,775), ylims=(ymin,ymax))#)#, xlims=(385,415), xticks=385:5:415) #, xticks=425:3:430)
q4 = plot!(1e7./νₐ, (ieQₐ_RRS_conv.+Qₐ_RRS_conv.-Qₐ_noRS_conv).*convfcta, linecolor=:red)#)#, xlims=(385,415), xticks=385:5:415) #, xticks=425:3:430)

#q2 = plot(1e7./ν₀, ieR₀₂.*convfct2, linecolor=:black, xlabel = "λ [nm]", xticks=424:2.5:429)
plot(p1, q1, p2, q2, p3, q3, p4, q4, layout = l, legend = false, title = ["Source, I₀ " "Source, Q₀" "No RS, O₂ A-band, Iₐ" "No RS, O₂ A-band, Qₐ" "VRS+RVRS, O₂ A-band, ΔIᵥ" "VRS+RVRS, O₂ A-band, ΔQᵥ" "RRS, O₂ A-band, ΔIᵣ" "RRS, O₂ A-band, ΔQᵣ"], titlefont = font(8))
savefig("VS_O2Aband_new.png")


l = @layout [a1 a2; b1 b2]
#=
ymax = round(maximum(R₀.*convfct0), sigdigits=2)
p1 = plot(1e7./ν₀, R₀.*convfct0, linecolor=:grey, yticks=0:(ymax/4):ymax)
p1 = plot!(1e7./ν₀, I₀_conv.*convfct0, linecolor=:magenta)#, xticks=385:5:415, ylabel="[mW/m²/str/nm]")
ymax = round(maximum(R₀_Q.*convfct0), sigdigits=1)
q1 = plot(1e7./ν₀, R₀_Q.*convfct0, linecolor=:grey, yticks=0:(ymax/4):ymax)#, xticks=385:5:415) #, xticks=425:3:430)
q1 = plot!(1e7./ν₀, Q₀_conv.*convfct0, linecolor=:magenta)#, xticks=385:5:415) #, xticks=425:3:430)

ymax = round(maximum(Rₐ.*convfcta), sigdigits=2)
p2 = plot(1e7./νₐ, Rₐ.*convfcta, linecolor=:black, yticks=0:(ymax/4):ymax, xlims=(755,775))#)#, xlims=(385,415), xticks=385:5:415, ylabel="[mW/m²/str/nm]")
p2 = plot!(1e7./νₐ, Iₐ_conv.*convfcta, linecolor=:red)#)#, xlims=(385,415), xticks=385:5:415, ylabel="[mW/m²/str/nm]")
ymax = round(maximum(Rₐ_Q.*convfcta), sigdigits=1)
q2 = plot(1e7./νₐ, Rₐ_Q.*convfcta, linecolor=:black, yticks=0:(ymax/4):ymax, xlims=(755,775))#)#, xlims=(385,415), xticks=385:5:415) #, xticks=425:3:430)
q2 = plot!(1e7./νₐ, Qₐ_conv.*convfcta, linecolor=:red)#)#, xlims=(385,415), xticks=385:5:415) #, xticks=425:3:430)
=#
ymax = round(maximum(100*ieRₐ./Rₐ), sigdigits=2)
p1 = plot(1e7./νₐ, 100*ieRₐ./Rₐ, linecolor=:black, yticks=0:(ymax/4):ymax, xlims=(755,775))#)#, xlims=(385,415), xticks=385:5:415, ylabel="[mW/m²/str/nm]")
p1 = plot!(1e7./νₐ, 100*ieIₐ_conv./Iₐ_conv, linecolor=:red)#)#, xlims=(385,415), xticks=385:5:415, ylabel="[mW/m²/str/nm]")
ymax = round(maximum(100*ieRₐ_Q./Rₐ_Q), sigdigits=2)
q1 = plot(1e7./νₐ, 100*ieRₐ_Q./Rₐ_Q, linecolor=:black, yticks=0:(ymax/4):ymax, xlims=(755,775))#)#, xlims=(385,415), xticks=385:5:415) #, xticks=425:3:430)
q1 = plot!(1e7./νₐ, 100*ieQₐ_conv./Qₐ_conv, linecolor=:red)#)#, xlims=(385,415), xticks=385:5:415) #, xticks=425:3:430)

ymax = round(maximum(100*(ieR[1,1,:].+R[1,1,:].-RnoRS[1,1,:])./Rₐ), sigdigits=2)
ymin = round(minimum(100*(ieR[1,1,:].+R[1,1,:].-RnoRS[1,1,:])./Rₐ), sigdigits=2)
p2 = plot(1e7./νₐ, 100*(ieR[1,1,:].+R[1,1,:].-RnoRS[1,1,:])./Rₐ, linecolor=:black, yticks=0:((ymax-ymin)/4):ymax, xlims=(755,775), ylims=(0,ymax))#)#, xlims=(385,415), xticks=385:5:415, ylabel="[mW/m²/str/nm]")
p2 = plot!(1e7./νₐ, 100*(ieIₐ_RRS_conv.+Iₐ_RRS_conv.-Iₐ_noRS_conv)./Iₐ_conv, linecolor=:red)#)#, xlims=(385,415), xticks=385:5:415, ylabel="[mW/m²/str/nm]")
ymax = round(maximum(100*(ieR[1,2,:].+R[1,2,:].-RnoRS[1,2,:])./Rₐ_Q), sigdigits=2)
ymin = round(minimum(100*(ieR[1,2,:].+R[1,2,:].-RnoRS[1,2,:])./Rₐ_Q), sigdigits=2)
q2 = plot(1e7./νₐ, 100*(ieR[1,2,:].+R[1,2,:].-RnoRS[1,2,:])./Rₐ_Q, linecolor=:black, yticks=ymin:((ymax-ymin)/4):ymax, xlims=(755,775), ylims=(ymin,ymax))#)#, xlims=(385,415), xticks=385:5:415) #, xticks=425:3:430)
q2 = plot!(1e7./νₐ, 100*(ieQₐ_RRS_conv.+Qₐ_RRS_conv.-Qₐ_noRS_conv)./Qₐ_conv, linecolor=:red)#)#, xlims=(385,415), xticks=385:5:415) #, xticks=425:3:430)

#q2 = plot(1e7./ν₀, ieR₀₂.*convfct2, linecolor=:black, xlabel = "λ [nm]", xticks=424:2.5:429)
plot(p1, q1, p2, q2, layout = l, legend = false, title = ["VRS+RVRS, O₂ A-band, 100*ΔIᵥ/Iₐ [%]" "VRS+RVRS, O₂ A-band, 100*ΔQᵥ/Qₐ [%]" "RRS, O₂ A-band, 100*ΔIᵣ/Iₐ [%]" "RRS, O₂ A-band, 100*ΔQᵣ/Qₐ [%]"], titlefont = font(8))
savefig("VS_O2Aband_new_percentages.png")

l = @layout [a1 a2; b1 b2]
ymax = round(maximum(100*R_ss[1,1,:]./R[1,1,:]), sigdigits=2)
ymin = round(minimum(100*R_ss[1,1,:]./R[1,1,:]), sigdigits=2)
p1 = plot(1e7./νₐ, 100*R_ss[1,1,:]./R[1,1,:], linecolor=:black, yticks=ymin:((ymax-ymin)/4):ymax, xlims=(755,775))#)#, xlims=(385,415), xticks=385:5:415, ylabel="[mW/m²/str/nm]")
p1 = plot!(1e7./νₐ, 100*Iₐ_RRS_ss_conv./Iₐ_RRS_conv, linecolor=:red)#)#, xlims=(385,415), xticks=385:5:415, ylabel="[mW/m²/str/nm]")
ymax = round(maximum(100*R_ss[1,2,:]./R[1,2,:]), sigdigits=2)
ymin = round(minimum(100*R_ss[1,2,:]./R[1,2,:]), sigdigits=2)
q1 = plot(1e7./νₐ, 100*R_ss[1,2,:]./R[1,2,:], linecolor=:black, yticks=ymin:((ymax-ymin)/4):ymax, xlims=(755,775))#)#, xlims=(385,415), xticks=385:5:415) #, xticks=425:3:430)
q1 = plot!(1e7./νₐ, 100*Qₐ_RRS_ss_conv./Qₐ_RRS_conv, linecolor=:red)#)#, xlims=(385,415), xticks=385:5:415) #, xticks=425:3:430)

ymax = 100.#round(maximum(100*ieR_ss[1,1,:]./ieR[1,1,:]), sigdigits=2)
ymin = round(minimum(100*ieR_ss[1,1,:]./ieR[1,1,:]), sigdigits=2)
p2 = plot(1e7./νₐ, 100*(ieR_ss[1,1,:]./ieR[1,1,:]), linecolor=:black, yticks=ymin:((ymax-ymin)/4):ymax, xlims=(755,775), ylims=(ymin,ymax))#)#, xlims=(385,415), xticks=385:5:415, ylabel="[mW/m²/str/nm]")
p2 = plot!(1e7./νₐ, 100*(ieIₐ_RRS_ss_conv./ieIₐ_RRS_conv), linecolor=:red)#)#, xlims=(385,415), xticks=385:5:415, ylabel="[mW/m²/str/nm]")
ymax = 100.#round(maximum(100*ieR_ss[1,2,:]./ieR[1,2,:]), sigdigits=2)
ymin = round(minimum(100*ieR_ss[1,2,:]./ieR[1,2,:]), sigdigits=2)
q2 = plot(1e7./νₐ, 100*(ieR_ss[1,2,:]./ieR[1,2,:]), linecolor=:black, yticks=ymin:((ymax-ymin)/4):ymax, xlims=(755,775), ylims=(ymin,ymax))#)#, xlims=(385,415), xticks=385:5:415) #, xticks=425:3:430)
q2 = plot!(1e7./νₐ, 100*ieQₐ_RRS_ss_conv./ieQₐ_RRS_conv, linecolor=:red)#)#, xlims=(385,415), xticks=385:5:415) #, xticks=425:3:430)

#q2 = plot(1e7./ν₀, ieR₀₂.*convfct2, linecolor=:black, xlabel = "λ [nm]", xticks=424:2.5:429)
plot(p1, q1, p2, q2, layout = l, legend = false, title = ["O₂ A-band, incl. RRS, 100*Iₛₛ/(Iₛₛ+Iₘₛ) [%]" "O₂ A-band, incl. RRS, 100*Qₛₛ/(Qₛₛ+Qₘₛ) [%]" "O₂ A-band, 100*ΔIᵣ₍ₛₛ₎/(ΔIᵣ₍ₛₛ₎+ΔIᵣ₍ₘₛ₎) [%]" "O₂ A-band, 100*ΔQᵣ₍ₛₛ₎/(ΔQᵣ₍ₛₛ₎+ΔQᵣ₍ₘₛ₎) [%]"], titlefont = font(8))
savefig("RRS_O2Aband_ss_percentages.png")
#RnoRS_test, TnoRS_test, _, _ = vSmartMOM.rt_run_test(vSmartMOM.noRS(),model,iBand);

#R_test, T_test, ieR_test, ieT_test = vSmartMOM.rt_run_test(RS_type,model,iBand);

# You can now run multiple bands like this (list at the end of band!)
#RnoRS_test, TnoRS_test, _, _ = vSmartMOM.rt_run_test(vSmartMOM.noRS(),model,[1,1]);


#R_test, T_test, ieR_test, ieT_test = vSmartMOM.rt_run_test(RS_type,model,1);