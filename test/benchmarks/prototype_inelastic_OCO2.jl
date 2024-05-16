##
using Revise
using Plots
using Statistics
using vSmartMOM
using vSmartMOM.Architectures
using vSmartMOM.Absorption
using vSmartMOM.Scattering
using vSmartMOM.CoreRT
using vSmartMOM.SolarModel
using vSmartMOM.InelasticScattering
using InstrumentOperator
using Interpolations
using DelimitedFiles

# Load OCO Data: 
# File names:
#dpath = "/net/fluo/data1/group/oco2/"
dpath = "/net/squid/data3/data/FluoData1/group/oco2/"
L1File   = dpath*"L1bSc/oco2_L1bScND_26782a_190715_B10003r_200429202458.h5"
metFile  = dpath*"L2Met/oco2_L2MetND_26782a_190715_B10003r_200429202458.h5"
dictFile = "/home/cfranken/code/gitHub/InstrumentOperator.jl/json/oco2.yaml"
# Load L1 file (could just use filenames here as well)
oco = InstrumentOperator.load_L1(dictFile,L1File, metFile);


# Pick some bands as tuple (or just one)
bands=(1,);
#bands = (1,2,3);
#bands = (1,3);
# Indices within that band:
#indices = (92:885,114:845,50:916);
indices = (92:885,);
#indices = (114:845,);
# Geo Index (footprint,sounding):
GeoInd = [5,5000];

# Get data for that sounding:
oco_sounding = InstrumentOperator.getMeasurement(oco, bands, indices, GeoInd);



# Load YAML files into parameter struct
parameters = parameters_from_yaml("test/test_parameters/O2Parameters2.yaml");
#parameters = parameters_from_yaml("test/test_parameters/CO2WParameters2.yaml");
# Create model struct (precomputes optical properties) from parameters
#parameters.depol = 0.028 #0.1032 #0.028 #ρ_Cabannes
#parameters.depol = 0.028 #0.1032 #0.028 #ρ_Cabannes
model      = model_from_parameters(parameters);


iBand = 1
FT = Float64
I_wl=[];
#I_conv=[];
#I_conv_RS=[];
#### Compute all Raman properties
#for iBand in length(bands)
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
        RS_type.F₀[1,i] = F₀[i];
    end 
    # For now, convert these special cases to the right array type:
    #aType = array_type(model.params.architecture);
    #RS_type.Z⁺⁺_λ₁λ₀ = aType(RS_type.Z⁺⁺_λ₁λ₀);
    #RS_type.ϖ_λ₁λ₀   = aType(RS_type.ϖ_λ₁λ₀);
    #RS_type.i_λ₁λ₀   = aType(RS_type.i_λ₁λ₀)
    #RS_type.Z⁻⁺_λ₁λ₀ = aType(RS_type.Z⁻⁺_λ₁λ₀);
    # Add something here, that computes ALL the OP needed for the Raman case.
    #modelRS = ...
    R, T, ieR, ieT =CoreRT.rt_run_test(RS_type, model, iBand)
    #mR, mT, mieR, mieT =CoreRT.rt_run_test(RS_type, model, iBand)
    
    #R_ss, T_ss, ieR_ss, ieT_ss =CoreRT.rt_run_test(RS_type, model, iBand)
    #=R, T, ieR, ieT = rt_run(RS_type,
        model.params.polarization_type,
        model.obs_geom,
        model.τ_rayl[iBand], 
        model.τ_aer[iBand], 
        model.quad_points,
        model.params.max_m,
        model.aerosol_optics[iBand],
        model.greek_rayleigh,
        model.τ_abs[iBand],
        model.params.brdf[iBand],
        model.params.architecture);
        =#
    λ_grid = reverse(1e4 ./ parameters.spec_bands[iBand])
    res = 0.001e-3;
    off = 0.5e-3
    wl = oco_sounding.ils[iBand].ν_out[1]-off:res:oco_sounding.ils[iBand].ν_out[end]+off;
    @show wl[1],wl[end], λ_grid[1],λ_grid[end]
    convfct = 1e7./ν.^2   # to convert radiance units from mW/m²-str-cm⁻¹ to mW/m²-str-nm

    #RR = zeros(length(parameters.vza), length(R[1,1,:]))
    # Full measurement:
    for icam in 1:length(parameters.vza)
        RR = oco_sounding.mueller[1]*(R[icam,1,:]+ieR[icam,1,:]) + 
            oco_sounding.mueller[2]*(R[icam,2,:]+ieR[icam,2,:]) + 
            oco_sounding.mueller[3]*(R[icam,3,:]+ieR[icam,3,:]);
        TT = oco_sounding.mueller[1]*(T[icam,1,:]+ieT[icam,1,:]) + 
            oco_sounding.mueller[2]*(T[icam,2,:]+ieT[icam,2,:]) + 
            oco_sounding.mueller[3]*(T[icam,3,:]+ieT[icam,3,:]);
        earth_out_r = reverse(RR[:].*convfct)
        earth_out_t = reverse(TT[:].*convfct)
        #for RR
        # Re-interpolate I from ν_grid to new grid/resolution
        @time interp_I = LinearInterpolation(λ_grid, earth_out_r);  
        tmp1 = interp_I(wl);
        # Convolve input spectrum with variable kernel
        @time tmp2 = InstrumentOperator.conv_spectra(oco_sounding.ils[iBand], wl, tmp1)
        #R = vSmartMOM.rt_run(model, i_band=1)
        out_R = [λ_grid*1.e3 R[icam,1,:] R[icam,2,:] R[icam,3,:]]
        writedlm("Aband_R_"*string(icam)*".dat", out_R)
        out_R_oco = [oco_sounding.SpectralGrid*1e3 tmp2]
        writedlm("ocoAband_R_"*string(icam)*".dat", out_R_oco)

        #for TT
        # Re-interpolate I from ν_grid to new grid/resolution
        @time interp_I = LinearInterpolation(λ_grid, earth_out_t);  
        tmp1 = interp_I(wl);
        # Convolve input spectrum with variable kernel
        @time tmp2 = InstrumentOperator.conv_spectra(oco_sounding.ils[iBand], wl, tmp1)
        #R = vSmartMOM.rt_run(model, i_band=1)
        out_T = [λ_grid*1.e3 T[icam,1,:] T[icam,2,:] T[icam,3,:]]
        writedlm("Aband_T_"*string(icam)*".dat", out_T)
        out_T_oco = [oco_sounding.SpectralGrid*1e3 tmp2]
        writedlm("ocoAband_T_"*string(icam)*".dat", out_T_oco)
    end
    #push!(I_conv_RS, tmp2);


    #============================================#
    #================no RS=======================#
    #============================================#
    #parameters.depol = 0.1032 #ρ_molec
    #parameters.depol = 0.028 #ρ_Cabannes
    model      = model_from_parameters(parameters);

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
        RS_type.F₀[1,i] = F₀[i];
    end 
    #cRnoRS, cTnoRS, _, _ = CoreRT.rt_run_test(RS_type,model,iBand);
    RnoRS, TnoRS, _, _ = CoreRT.rt_run_test(RS_type,model,iBand);
    
    #λ_grid = reverse(1e4 ./ parameters.spec_bands[iBand])
    #res = 0.001e-3;
    #off = 0.5e-3
    #wl = oco_sounding.ils[iBand].ν_out[1]-off:res:oco_sounding.ils[iBand].ν_out[end]+off;
    #@show wl[1],wl[end], λ_grid[1],λ_grid[end]
    #convfct = 1e7./ν.^2   # to convert radiance units from mW/m²-str-cm⁻¹ to mW/m²-str-nm

    #RR = zeros(length(parameters.vza), length(R[1,1,:]))
    # Full measurement:
    for icam in 1:length(parameters.vza)
        RR_ = oco_sounding.mueller[1]*RnoRS[icam,1,:] + 
            oco_sounding.mueller[2]*RnoRS[icam,2,:] + 
            oco_sounding.mueller[3]*RnoRS[icam,3,:];
        TT_ = oco_sounding.mueller[1]*TnoRS[icam,1,:] + 
            oco_sounding.mueller[2]*TnoRS[icam,2,:] + 
            oco_sounding.mueller[3]*TnoRS[icam,3,:];
        earth_out_r = reverse(RR_[:].*convfct)
        earth_out_t = reverse(TT_[:].*convfct)
        #for RR
        # Re-interpolate I from ν_grid to new grid/resolution
        @time interp_I = LinearInterpolation(λ_grid, earth_out_r);  
        tmp1 = interp_I(wl);
        # Convolve input spectrum with variable kernel
        @time tmp2 = InstrumentOperator.conv_spectra(oco_sounding.ils[iBand], wl, tmp1)
        #R = vSmartMOM.rt_run(model, i_band=1)
        out_R = [λ_grid*1.e3 RnoRS[icam,1,:] RnoRS[icam,2,:] RnoRS[icam,3,:]]
        writedlm("Aband_RnoRS_"*string(icam)*".dat", out_R)
        out_R_oco = [oco_sounding.SpectralGrid*1e3 tmp2]
        writedlm("ocoAband_RnoRS_"*string(icam)*".dat", out_R_oco)

        #for TT
        # Re-interpolate I from ν_grid to new grid/resolution
        @time interp_I = LinearInterpolation(λ_grid, earth_out_t);  
        tmp1 = interp_I(wl);
        # Convolve input spectrum with variable kernel
        @time tmp2 = InstrumentOperator.conv_spectra(oco_sounding.ils[iBand], wl, tmp1)
        #R = vSmartMOM.rt_run(model, i_band=1)
        out_T = [λ_grid*1.e3 TnoRS[icam,1,:] TnoRS[icam,2,:] TnoRS[icam,3,:]]
        writedlm("Aband_TnoRS_"*string(icam)*".dat", out_T)
        out_T_oco = [oco_sounding.SpectralGrid*1e3 tmp2]
        writedlm("ocoAband_TnoRS_"*string(icam)*".dat", out_T_oco)
    end
    #=
    
    RR_ = oco_sounding.mueller[1]*(RnoRS[iBand,1,:]) + 
          oco_sounding.mueller[2]*(RnoRS[iBand,2,:]) + 
          oco_sounding.mueller[3]*(RnoRS[iBand,3,:]);
    earth_out = reverse(RR_[:])
    # Re-interpolate I from ν_grid to new grid/resolution
    λ_grid = reverse(1e4 ./ parameters.spec_bands[iBand])
    @time interp_I = LinearInterpolation(λ_grid, earth_out);

    res = 0.001e-3;
    off = 0.5e-3
    wl = oco_sounding.ils[iBand].ν_out[1]-off:res:oco_sounding.ils[iBand].ν_out[end]+off;
    @show wl[1],wl[end], λ_grid[1],λ_grid[end]
    tmp1 = interp_I(wl);
    push!(I_wl, tmp1);

    # Convolve input spectrum with variable kernel
    @time tmp2 = InstrumentOperator.conv_spectra(oco_sounding.ils[iBand], wl, tmp1)
    =#
    #push!(I_conv, tmp2);
#end
for icam=1:length(parameters.vza)
    in_R = readdlm("Aband_R_"*string(icam)*".dat")
    in_RnoRS = readdlm("Aband_RnoRS_"*string(icam)*".dat")
    in_ocoR = readdlm("ocoAband_R_"*string(icam)*".dat")
    in_ocoRnoRS = readdlm("ocoAband_RnoRS_"*string(icam)*".dat")

    in_T = readdlm("Aband_T_"*string(icam)*".dat")
    in_TnoRS = readdlm("Aband_TnoRS_"*string(icam)*".dat")
    in_ocoT = readdlm("ocoAband_T_"*string(icam)*".dat")
    in_ocoTnoRS = readdlm("ocoAband_TnoRS_"*string(icam)*".dat")
    
    wl = in_R[:,1]
    rI  = in_R[:,2]
    rQ  = in_R[:,3]
    rU  = in_R[:,4]

    #wl = in_R[:,1]
    rInoRS  = in_RnoRS[:,2]
    rQnoRS  = in_RnoRS[:,3]
    rUnoRS  = in_RnoRS[:,4]

    oco_wl = in_ocoR[:,1]
    oco_R = in_ocoR[:,2]
    oco_RnoRS = in_ocoRnoRS[:,2]

    #wl = in_R[:,1]
    tI  = in_T[:,2]
    tQ  = in_T[:,3]
    tU  = in_T[:,4]

    #wl = in_R[:,1]
    tInoRS  = in_TnoRS[:,2]
    tQnoRS  = in_TnoRS[:,3]
    tUnoRS  = in_TnoRS[:,4]

    #oco_wl = in_ocoT[:,1]
    oco_T = in_ocoT[:,2]
    oco_TnoRS = in_ocoTnoRS[:,2]   

#l = @layout [a b]
    l = @layout [a1 a2 ; b]
    #i=1
    p1 = plot(wl, rInoRS, label="Elastic", ylabel = "I [mW/m²/str/nm]")
    p1 = plot!(wl, rI, label="w/ RRS", xlabel = "λ [nm]")

    p2 = plot(wl, rQnoRS, label="Elastic", ylabel = "Q [mW/m²/str/nm]")
    p2 = plot!(wl, rQ, label="w/ RRS", xlabel = "λ [nm]")

    p3 = plot(oco_sounding.SpectralGrid*1e3, oco_RnoRS, label="Elastic", ylabel = "Meas. y  [mW/m²/str/nm]")
    p3 = plot!(oco_sounding.SpectralGrid*1e3, oco_R, label="w/ RRS", xlabel = "λ [nm]")

    plot(p1, p2, p3, layout = l, legend = true, title = ["" "" "θᵥ="*string(parameters.vza[icam])*"ᵒ, Δϕ="*string(parameters.vaz[icam])*"ᵒ"], titlefont = font(10))

    savefig("O2_RRS_SZA30_vgeom"*string(icam)*"_aer0p0_alb0.png")

    l = @layout [a1 a2]
    q1 = plot(oco_sounding.SpectralGrid*1e3, 
        (oco_R - oco_RnoRS)/maximum(oco_RnoRS) * 100, 
        #ylims=(-5,5), 
        label="(RRS-NoRS)/max(noRS)*100")
    q2 = plot(oco_sounding.SpectralGrid*1e3, 
        (oco_R .- oco_RnoRS)./(oco_RnoRS) * 100, 
        #ylims=(-5,5), 
        label="(RRS-NoRS)/noRS*100")
    plot(q1,q2,link=:y,layout=l)
    savefig("O2_RRS_impact_SZA30_vgeom"*string(icam)*"_aer0p0_alb0.png")
end
    #savefig("CO2W_RRS_impact.pdf")
#=
l = @layout [a b c; d e f]
i=1;
p1 = plot(I_wl[i]*1e3, I_conv[i], label="Elastic")
plot!(I_wl[1]*1e3, I_conv_RS[i], label="With RRS")
p4 = plot(I_wl[i]*1e3, (I_conv_RS[i] .- I_conv[i])/maximum(I_conv[i]) * 100, label="(RRS-NoRS)/max(noRS)*100")

i=2;
p2 = plot(I_wl[i]*1e3, I_conv[i], label="Elastic")
plot!(I_wl[1]*1e3, I_conv_RS[i], label="With RRS")
p5 = plot(I_wl[i]*1e3, (I_conv_RS[i] .- I_conv[i])/maximum(I_conv[i]) * 100, label="(RRS-NoRS)/max(noRS)*100")

i=3;
p3 = plot(I_wl[i]*1e3, I_conv[i], label="Elastic")
plot!(I_wl[1]*1e3, I_conv_RS[i], label="With RRS")
p6 = plot(I_wl[i]*1e3, (I_conv_RS[i] .- I_conv[i])/maximum(I_conv[i]) * 100, label="(RRS-NoRS)/max(noRS)*100")

plot(p1, p2, p3, p4, p5, p6, layout=l)
savefig("RRS_impact.pdf")
=#
ip=2
#=plot(mR[1,ip,:], color=:blue, alpha=1)
plot!(R[1,ip,:], color=:red, alpha=1)

plot!(mR[2,ip,:], color=:blue, alpha=0.75)
plot!(R[2,ip,:], color=:red, alpha=0.75)

plot!(mR[3,ip,:], color=:blue, alpha=0.5)
plot!(R[3,ip,:], color=:red, alpha=0.5)

plot!(mR[4,ip,:], color=:blue, alpha=0.25)
plot!(R[4,ip,:], color=:red, alpha=0.25)
=#
plot(mR[1,ip,:]./R[1,ip,:], color=:blue, alpha=1)
plot!(RnoRS[1,ip,:]./cRnoRS[1,ip,:], color=:red, alpha=1)
plot!(mT[1,ip,:]./T[1,ip,:], color=:green, alpha=1)
plot!(TnoRS[1,ip,:]./cTnoRS[1,ip,:], color=:purple, alpha=1)
for icam=2:4
    tra=1-(icam-1)*0.25
    plot!(mR[icam,ip,:]./R[icam,ip,:], color=:blue, alpha=tra)
    plot!(RnoRS[icam,ip,:]./cRnoRS[icam,ip,:], color=:red, alpha=tra)
    plot!(mT[icam,ip,:]./T[icam,ip,:], color=:green, alpha=tra)
    plot!(TnoRS[icam,ip,:]./cTnoRS[icam,ip,:], color=:purple, alpha=tra)
end
plot!(mR[2,ip,:]./R[2,ip,:], color=:blue, alpha=0.75)
plot!(RnoRS[2,ip,:]./cRnoRS[2,ip,:], color=:red, alpha=0.75)

plot!(mR[3,ip,:]./R[3,ip,:], color=:blue, alpha=0.5)
plot!(RnoRS[3,ip,:]./cRnoRS[3,ip,:], color=:red, alpha=0.5)

plot!(mR[4,ip,:]./R[4,ip,:], color=:blue, alpha=0.25)
plot!(RnoRS[4,ip,:]./cRnoRS[4,ip,:], color=:red, alpha=0.25)
