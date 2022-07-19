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
        

ν₀ = (1e7/445.):0.3:(1e7/335.)
λ₀ = 2e7/(ν₀[1]+ν₀[end]) #400. #nm
ii₀ = findmin(abs.(1e7./ν₀ .- λ₀))[2]
# Computing Fraunhofer Spectrum 
T_sun = 5777. # KRS 
#Tsolar = solar_transmission_from_file("/home/sanghavi/code/github/vSmartMOM.jl/src/SolarModel/solar.out")
#Tsolar_interp = LinearInterpolation(Tsolar[4:end, 1], Tsolar[4:end, 2])
P = planck_spectrum_wn(T_sun, collect(ν₀))    
F₀ = zeros(length(ν₀));
#ν₀ = 1e7/450.:2.:1e7/375.
R₀ = zeros(3,length(ν₀)) 
#T₀ = zeros(length(ν₀)) 
ieR₀ = zeros(3,length(ν₀)) 
#ieT₀ = zeros(length(ν₀))
R₀_Q = zeros(3,length(ν₀)) 
#T₀_Q = zeros(length(ν₀)) 
ieR₀_Q = zeros(3,length(ν₀)) 
#ieT₀_Q = zeros(length(ν₀))
#R₀_U = zeros(3,length(ν₀)) 
#T₀_U = zeros(length(ν₀)) 
#ieR₀_U = zeros(3,length(ν₀)) 
#ieT₀_U = zeros(length(ν₀))
I₀_conv = zeros(3,length(ν₀)) 
ieI₀_conv = zeros(3,length(ν₀)) 
Q₀_conv = zeros(3,length(ν₀)) 
ieQ₀_conv = zeros(3,length(ν₀)) 
for ictr = 1:3
    for i=1:length(ν₀)
        sol_trans = Tsolar_interp(ν₀[i]);
        F₀[i] = sol_trans * P[i];
    end 
        
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
    model      = model_from_parameters(RS_type, λ₀, parameters);
    Fraunhofer=false
    #####TEMPORARY!!#####
    if ictr==2
        RS_type.ϖ_λ₁λ₀_VS_n2[:].=0.0
    elseif ictr==3
        RS_type.ϖ_λ₁λ₀_VS_o2[:].=0.0
    end
    #####################
    RS_type.F₀ = zeros(model.params.polarization_type.n, 
        sum(length(RS_type.grid_in[i]) for i in 1:length(RS_type.iBand)))
    global t_offset = 0
    for iB=1:length(RS_type.iBand)
        ν = collect(RS_type.grid_in[iB])
        #global t_offset
        for i=1:length(ν)
            RS_type.F₀[1,i+t_offset] = 1.0;
        end 
        t_offset += length(ν)
    end

    R, T, ieR, ieT = CoreRT.rt_run_test(RS_type,
            model, RS_type.iBand);

    for ctr=1:length(ν₀)
        
        for iB = 1:length(RS_type.iBand)
            if iB==1
                #@show iB, ctr, R[1,1,1]
                R₀[ictr,ctr] = R[1,1,1]*F₀[ctr]
                #T₀[ctr] = T[1,1,1]*F₀[ctr]
                R₀_Q[ictr,ctr] = R[1,2,1]*F₀[ctr]
                #T₀_Q[ctr] = T[1,2,1]*F₀[ctr]
                #R₀_U[ictr,ctr] = R[1,3,1]*F₀[ctr]
                #T₀_U[ctr] = T[1,3,1]*F₀[ctr]
        #        noRS_R₀[ctr] = RnoRS[1][1]
        #        noRS_T₀[ctr] = TnoRS[1][1]
            else
                jj₀ = findmin(abs.(ν₀.-RS_type.grid_in[iB][1]))[2]
                t_off = ctr-ii₀+jj₀-1
                for i=1:length(RS_type.grid_in[iB])
                    t_idx = t_off+i
                    #if (RS_type.grid_in[iB][1]<=ν₀[i]<=RS_type.grid_in[iB][end])
                    if (1<=t_idx<=length(ν₀))
                        ieR₀[ictr,t_idx] += ieR[1,1,RS_type.bandSpecLim[iB][i]]*F₀[ctr];
                        #ieT₀[t_idx] += ieT[1,1,RS_type.bandSpecLim[iB][i]]*F₀[ctr];
                        ieR₀_Q[ictr,t_idx] += ieR[1,2,RS_type.bandSpecLim[iB][i]]*F₀[ctr];
                        #ieT₀_Q[t_idx] += ieT[1,2,RS_type.bandSpecLim[iB][i]]*F₀[ctr];
                        #ieR₀_U[ictr,t_idx] += ieR[1,3,RS_type.bandSpecLim[iB][i]]*F₀[ctr];
                        #ieT₀_U[t_idx] += ieT[1,3,RS_type.bandSpecLim[iB][i]]*F₀[ctr];
                    end
                end 
            end
        end
        @show ictr, ctr, ν₀[ctr]
    end
    x = -40:0.3:40
    #kernel = InstrumentOperator.create_instrument_kernel(Normal(0, 12.5), x) #defining a Gaussian kernel for convulution in wavenumber space
    kernel = InstrumentOperator.create_instrument_kernel(Normal(0, 12.5), x)
    #I_conv = InstrumentOperator.conv_spectra(kernel, )
    
    I₀_conv[ictr,:] .= imfilter(R₀[ictr,:], kernel)
    ieI₀_conv[ictr,:] .= imfilter(ieR₀[ictr,:], kernel)

    Q₀_conv[ictr,:] = imfilter(R₀_Q[ictr,:], kernel)
    ieQ₀_conv[ictr,:] = imfilter(ieR₀_Q[ictr,:], kernel)    
end


convfct0 = 1e7./ν₀.^2   # to convert radiance units from mW/m²-str-cm⁻¹ to mW/m²-str-nm

l = @layout [a1 a2; b1 b2; c1 c2; d1 d2]
ymax = round(maximum(R₀[1,:].*convfct0), sigdigits=2)
p1 = plot(1e7./ν₀, R₀[1,:].*convfct0, linecolor=:grey, xlims=(380,445), yticks=0:(ymax/4):ymax)
p1 = plot!(1e7./ν₀, I₀_conv[1,:].*convfct0, linecolor=:magenta, xlims=(380,445))#, xticks=385:5:415, ylabel="[mW/m²/str/nm]")
ymax = round(maximum(R₀_Q[1,:].*convfct0), sigdigits=1)
q1 = plot(1e7./ν₀, R₀_Q[1,:].*convfct0, linecolor=:grey, xlims=(380,445), yticks=0:(ymax/4):ymax)#, xticks=385:5:415) #, xticks=425:3:430)
q1 = plot!(1e7./ν₀, Q₀_conv[1,:].*convfct0, linecolor=:magenta, xlims=(380,445))#, xticks=385:5:415) #, xticks=425:3:430)

ymax = round(maximum(ieR₀[1,:].*convfct0), sigdigits=2)
p2 = plot(1e7./ν₀, ieR₀[1,:].*convfct0, linecolor=:black, xlims=(380,445), yticks=0:(ymax/4):ymax)#)#, xlims=(385,415), xticks=385:5:415, ylabel="[mW/m²/str/nm]")
p2 = plot!(1e7./ν₀, ieI₀_conv[1,:].*convfct0, linecolor=:red, xlims=(380,445))#)#, xlims=(385,415), xticks=385:5:415, ylabel="[mW/m²/str/nm]")
ymax = round(maximum(ieR₀_Q[1,:].*convfct0), sigdigits=1)
q2 = plot(1e7./ν₀, ieR₀_Q[1,:].*convfct0, linecolor=:black, xlims=(380,445), yticks=0:(ymax/4):ymax)#)#, xlims=(385,415), xticks=385:5:415) #, xticks=425:3:430)
q2 = plot!(1e7./ν₀, ieQ₀_conv[1,:].*convfct0, linecolor=:red, xlims=(380,445))#)#, xlims=(385,415), xticks=385:5:415) #, xticks=425:3:430)

ymax = round(maximum(ieR₀[2,:].*convfct0), sigdigits=2)
p3 = plot(1e7./ν₀, ieR₀[2,:].*convfct0, linecolor=:black, xlims=(380,445), yticks=0:(ymax/4):ymax, ylabel="[mW/m²/str/nm]", yguidefontvalign = :bottom)#)#, xlims=(385,415), xticks=385:5:415, ylabel="[mW/m²/str/nm]")
p3 = plot!(1e7./ν₀, ieI₀_conv[2,:].*convfct0, linecolor=:red, xlims=(380,445))#)#, xlims=(385,415), xticks=385:5:415, ylabel="[mW/m²/str/nm]")
ymax = round(maximum(ieR₀_Q[2,:].*convfct0), sigdigits=1)
q3 = plot(1e7./ν₀, ieR₀_Q[2,:].*convfct0, linecolor=:black, xlims=(380,445), yticks=0:(ymax/4):ymax)#)#, xlims=(385,415), xticks=385:5:415) #, xticks=425:3:430)
q3 = plot!(1e7./ν₀, ieQ₀_conv[2,:].*convfct0, linecolor=:red, xlims=(380,445))#)#, xlims=(385,415), xticks=385:5:415) #, xticks=425:3:430)

ymax = round(maximum(ieR₀[3,:].*convfct0), sigdigits=2)
p4 = plot(1e7./ν₀, ieR₀[3,:].*convfct0, linecolor=:black, xlabel = "λ [nm]", yticks=0:(ymax/4):ymax, xlims=(380,445))#)#, xlims=(395,405), xticks=395:5:405
p4 = plot!(1e7./ν₀, ieI₀_conv[3,:].*convfct0, linecolor=:red, xlabel = "λ [nm]", xlims=(380,445))#)#, xlims=(395,405), xticks=395:5:405
ymax = round(maximum(ieR₀_Q[3,:].*convfct0), sigdigits=1)
q4 = plot(1e7./ν₀, ieR₀_Q[3,:].*convfct0, linecolor=:black, xlabel = "λ [nm]", xlims=(380,445), yticks=0:(ymax/4):ymax)#)#
q4 = plot!(1e7./ν₀, ieQ₀_conv[3,:].*convfct0, linecolor=:red, xlabel = "λ [nm]", xlims=(380,445))#)#
#q2 = plot(1e7./ν₀, ieR₀₂.*convfct2, linecolor=:black, xlabel = "λ [nm]", xticks=424:2.5:429)
plot(p1, q1, p2, q2, p3, q3, p4, q4, layout = l, legend = false, title = ["Iₑ " "Qₑ" "Total Iᵢ" "Total Qᵢ" "Iᵢ by O₂" "Qᵢ by O₂" "Iᵢ by N₂" "Qᵢ by N₂"], titlefont = font(8))
savefig("VS_ghosting.png")

#RnoRS_test, TnoRS_test, _, _ = vSmartMOM.rt_run_test(vSmartMOM.noRS(),model,iBand);

#R_test, T_test, ieR_test, ieT_test = vSmartMOM.rt_run_test(RS_type,model,iBand);

# You can now run multiple bands like this (list at the end of band!)
#RnoRS_test, TnoRS_test, _, _ = vSmartMOM.rt_run_test(vSmartMOM.noRS(),model,[1,1]);


#R_test, T_test, ieR_test, ieT_test = vSmartMOM.rt_run_test(RS_type,model,1);