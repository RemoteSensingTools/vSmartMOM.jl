# Runner is used to set AD fields as duals
function runner!(y, x, parameters=parameters, oco_sounding= oco_soundings, Tsolar = Tsolar_interp)
    Nangles = length(oco_sounding);
    #@show Nangles
    # Set parameters fields as the dual numbers
    parameters.brdf =  [CoreRT.LambertianSurfaceLegendre([x[1],x[3],x[4]]),
                        CoreRT.LambertianSurfaceLegendre([x[7],x[8],x[5]]),
                        CoreRT.LambertianSurfaceLegendre([x[9],x[10],x[6]])];

    parameters.scattering_params.rt_aerosols[1].τ_ref = exp(x[2]);
    @show  x[2], exp(x[2])
    parameters.scattering_params.rt_aerosols[1].profile    = Normal(x[20], x[21]); #800.0; #x[4]
    #parameters.scattering_params.rt_aerosols[1].profile.σp    = x[21];
    @show  x[20], x[21]
    parameters.scattering_params.rt_aerosols[1].aerosol.size_distribution = LogNormal(x[18], x[19]); #x[4]
    @show  x[18], x[19]
    parameters.scattering_params.rt_aerosols[1].aerosol.nᵣ = x[16]
    parameters.scattering_params.rt_aerosols[1].aerosol.nᵢ = x[17]
    @show  x[16], x[17]
    #parameters.p   = oco_sounding.p_half
    #parameters.q   = oco_sounding.q 
    #parameters.T   = oco_sounding.T# .+ 1.0 #.+ x[15]
    #parameters.sza = oco_sounding.sza
    #parameters.vza = [oco_sounding.vza]
    for ia=1:Nangles
        #if ia==1
            parameters.p   = oco_sounding[ia].p_half
            parameters.q   = oco_sounding[ia].q 
            parameters.T   = oco_sounding[ia].T# .+ 1.0 #.+ x[15]
            parameters.sza = oco_sounding[ia].sza
            #parameters.vza = [oco_sounding[ia].vza]   
            #parameters.vaz = [oco_sounding[ia].raa]
        #else
        #    parameters.vza = vcat(parameters.vza, oco_sounding[ia].vza)
        #    parameters.vaz = vcat(parameters.vaz, oco_sounding[ia].raa)
        #end
    end
    parameters.absorption_params.vmr["H2O"] = [parameters.q[1:65]*x[11] * 1.8; 
                                               parameters.q[66:end]*x[15] * 1.8];
    a1 = zeros(7) .+ x[12]
    a2 = zeros(7) .+ x[13]
    a3 = zeros(6) .+ x[14]
    parameters.absorption_params.vmr["CO2"] = [a1; a2; a3]
    model = model_from_parameters(parameters);
    RS_type = InelasticScattering.noRS(
        fscattRayl  = [FT(1)],
        ϖ_Cabannes  = [FT(1)], 
        bandSpecLim = [],
        iBand       = [1],
        F₀          = zeros(FT,1,1));
    
    R=[];
    for i = 1:length(oco_sounding[1].BandID)
        println("Modelling band $(i)")
        # Run the model to obtain reflectance matrix
        #R = rt_run(model, i_band=i)[1];
        #<<
        T_sun = 5777. # K
        P = planck_spectrum_wn(T_sun, parameters.spec_bands[i])
        F₀ = zeros(length(P));
        RS_type.F₀ = zeros(model.params.polarization_type.n, length(P))
        for ii=1:length(P)
            sol_trans = Tsolar(parameters.spec_bands[i][ii]);
            F₀[ii] = sol_trans * P[ii];
            RS_type.F₀[i,ii] = F₀[ii];
        end 
        #>>
        #tmpR = CoreRT.rt_run_test(CoreRT.noRS(), model, i)[1]
        #@show  "HAHAHHAHAHHHAHHAH", model
        tmpR = CoreRT.rt_run_test(RS_type, model, i)[1]
        #@show  "TTTTTAHAHHAHAHHHAHHAH"
        push!(R, tmpR);
    end
    for ia = 1:Nangles
        for i = 1:length(oco_sounding[ia].BandID)
            println("Modelling band $(i)")
            # Run the model to obtain reflectance matrix
            #R = rt_run(model, i_band=i)[1];
            #R = CoreRT.rt_run_test(CoreRT.noRS(), model, i)[1]
            # Re-interpolate I from ν_grid to new grid/resolution
            λ_grid = reverse(1e4 ./ parameters.spec_bands[i])
            res = 0.001e-3;
            off = 0.5e-3
            # Get sun:
            #@time sun_out = getSolar(parameters.spec_bands[i],Tsolar)
        
            RR = oco_sounding[ia].mueller[1]*R[i][ia,1,:] + 
                    oco_sounding[ia].mueller[2]*R[i][ia,2,:] + 
                    oco_sounding[ia].mueller[3]*R[i][ia,3,:];        
        
            # Apply Earth reflectance matrix 
            #earth_out = sun_out .* reverse(RR[:])
            earth_out = reverse(RR[:])
            
            # Re-interpolate I from ν_grid to new grid/resolution
            @time interp_I = LinearInterpolation(λ_grid, earth_out);        
            #The following holds only for the use of a single ground pixel (because the ils changes with ground pixel index)
            wl = oco_sounding[ia].ils[i].ν_out[1]-off:res:oco_sounding[ia].ils[i].ν_out[end]+off;
            @show wl[1],wl[end], λ_grid[1],λ_grid[end]
            I_wl = interp_I(wl);

            # Convolve input spectrum with variable kernel
            @time I_conv = InstrumentOperator.conv_spectra(oco_sounding[ia].ils[i], wl, I_wl)
            off = oco_soundings[1].BandID[end][end] * (ia-1)
            y[oco_sounding[ia].BandID[i] .+ off ] = I_conv
            #if ia==1
            #    y[oco_sounding[ia].BandID[i]] = I_conv
            #else
            #    # not sure this will work!
            #    y[oco_sounding[ia].BandID[i]] = vcat(y[oco_sounding.BandID[i]], I_conv)
            #end
        end
    end
end

function lin_runner!(y, x, parameters=parameters, oco_sounding= oco_soundings, Tsolar = Tsolar_interp)
    Nangles = length(oco_sounding);
    #@show Nangles
    # New parameter order:
    # 1: p_surf
    # (2-7) CO2: VMR₀, VMR₁, VMR₂, HCO2, z₀CO2, σCO2
    # (8,9) H2O: a1, a2
    # (10-16) Aerosol1: ln(τ_ref), nᵣ, nᵢ, r₀, σᵣ, z₀, σz
    # (17-23) Aerosol2: ln(τ_ref), nᵣ, nᵢ, r₀, σᵣ, z₀, σz  
    # (24-32) Albedo: ρA₁, ρA₂, ρA₃, ρW₁, ρW₂, ρW₃, ρS₁, ρS₂, ρS₃
    # BRDF: RPV/RossLi/CoxMunk   
    # Set parameters fields as the dual numbers
    parameters.brdf =  [CoreRT.LambertianSurfaceLegendre([x[24],x[25],x[26]]),
                        CoreRT.LambertianSurfaceLegendre([x[27],x[28],x[29]]),
                        CoreRT.LambertianSurfaceLegendre([x[30],x[31],x[32]])];

    parameters.scattering_params.rt_aerosols[1].τ_ref = exp(x[10]);
    @show  x[10], exp(x[10])
    parameters.scattering_params.rt_aerosols[1].aerosol.nᵣ = x[11]
    parameters.scattering_params.rt_aerosols[1].aerosol.nᵢ = x[12]
    @show  x[11], x[12]
    parameters.scattering_params.rt_aerosols[1].aerosol.size_distribution = 
                LogNormal(x[13], x[14]); #x[4]
    @show  x[13], x[14]
    parameters.scattering_params.rt_aerosols[1].profile    = 
                LogNormal(x[15], x[16]); #800.0; #x[4]
    @show  x[15], x[16]
    
    parameters.scattering_params.rt_aerosols[2].τ_ref = exp(x[17]);
    @show  x[17], exp(x[17])
    parameters.scattering_params.rt_aerosols[2].aerosol.nᵣ = x[18]
    parameters.scattering_params.rt_aerosols[2].aerosol.nᵢ = x[19]
    @show  x[18], x[19]
    parameters.scattering_params.rt_aerosols[2].aerosol.size_distribution = 
                LogNormal(x[20], x[21]); #x[4]
    @show  x[20], x[21]
    parameters.scattering_params.rt_aerosols[2].profile    = 
                LogNormal(x[22], x[23]); #800.0; #x[4]
    @show  x[22], x[23]

    nparams=Int(length(x))
    #parameters.p   = oco_sounding.p_half
    #parameters.q   = oco_sounding.q 
    #parameters.T   = oco_sounding.T# .+ 1.0 #.+ x[15]
    #parameters.sza = oco_sounding.sza
    #parameters.vza = [oco_sounding.vza]
    psurf = x[1] 
    for ia=1:Nangles
        #if ia==1
            parameters.p   = oco_sounding[ia].p_half
            #inducing psurf
            parameters.p[end] = psurf
            #@show parameters.q#[end], parameters.p[end-1]
            @assert( parameters.p[end] >= parameters.p[end-1] )
            parameters.q   = oco_sounding[ia].q 
            parameters.T   = oco_sounding[ia].T# .+ 1.0 #.+ x[15]
            parameters.sza = oco_sounding[ia].sza

            #parameters.vza = [oco_sounding[ia].vza]   
            #parameters.vaz = [oco_sounding[ia].raa]
        #else
        #    parameters.vza = vcat(parameters.vza, oco_sounding[ia].vza)
        #    parameters.vaz = vcat(parameters.vaz, oco_sounding[ia].raa)
        #end
    end
    #dVMR_H2O = zeros(FT, 2, length(q)) # x[8], x[9]
    #dVMR_CO2 = zeros(FT, 7, length(q)) # x[1], x[2], x[3], x[4], x[5], x[6], x[7] 

   
    parameters.absorption_params.vmr["H2O"] = [zeros(65).+ x[8]; 
                                               zeros(length(parameters.q)-65) .+ x[9]];
    #a1 = zeros(7) .+ x[12]
    #a2 = zeros(7) .+ x[13]
    #a3 = zeros(6) .+ x[14]
    #H = 8.0  # km (atmospheric scale height)
    #z = H*log.(psurf./parameters.p)
    #prof  = LogNormal(x[6], x[7])
    
    parameters.absorption_params.vmr["CO2"] = (x[2].+zeros(length(parameters.q))) #+ 
                                            #(x[3] * exp.(-z./x[5])) +
                                            #(x[4] * pdf.(prof, z))   #[a1; a2; a3]
    
    #dVMR_H2O[1,:] = [parameters.q[1:65] * 1.8; parameters.q[66:end] * 0.0;] # wrt x[7]
    #dVMR_H2O[2,:] = [parameters.q[1:65] * 0.0; parameters.q[66:end] * 1.8;] # wrt x[8]

    #dVMR_CO2[1,:] = x[3] * exp.(-z./x[5]) * (-H/x[5]) * (1/x[1]); # wrt x[1] 
    #dVMR_CO2[2,:] = 1.0 .+ zeros(length(z)) # wrt x[2]
    #dVMR_CO2[3,:] = exp.(-z./x[5]) # wrt x[3]
    #dVMR_CO2[4,:] = pdf.(prof, z) # wrt x[4]
    #dVMR_CO2[5,:] = x[3] * exp.(-z./x[5]) .* z./(x[5])^2 # wrt x[5]
    #dVMR_CO2[6,:] = x[4] * pdf.(prof, z) .* (log.(z) .- x[6]) / x[7]^2
    #dVMR_CO2[7,:] = (x[4] * pdf.(prof, z) / x[7]) .* 
    #                    (((log.(z) - x[6]) / x[7]).^2 .- 1)

    model, lin = lin_model_from_parameters(parameters, x[1:9])#, dVMR_CO2, dVMR_H2O);
    RS_type = InelasticScattering.noRS(
        fscattRayl  = [FT(1)],
        ϖ_Cabannes  = [FT(1)], 
        bandSpecLim = [],
        iBand       = [1],
        F₀          = zeros(FT,1,1));
    
    R=[];
    for i = 1:length(oco_sounding[1].BandID)
        println("Modelling band $(i)")
        # Run the model to obtain reflectance matrix
        #R = rt_run(model, i_band=i)[1];
        #<<
        T_sun = 5777. # K
        P = planck_spectrum_wn(T_sun, parameters.spec_bands[i])
        F₀ = zeros(length(P));
        RS_type.F₀ = zeros(model.params.polarization_type.n, length(P))
        for ii=1:length(P)
            sol_trans = Tsolar(parameters.spec_bands[i][ii]);
            F₀[ii] = sol_trans * P[ii];
            RS_type.F₀[i,ii] = F₀[ii];
        end 
        #>>
        #tmpR = CoreRT.rt_run_test(CoreRT.noRS(), model, i)[1]
        #@show  "HAHAHHAHAHHHAHHAH", model
        
        tmpR = CoreRT.rt_run_test(RS_type, model, lin, nparams, i)[1]
        
        #@show  "TTTTTAHAHHAHAHHHAHHAH"
        push!(R, tmpR);
    end
    for ia = 1:Nangles
        for i = 1:length(oco_sounding[ia].BandID)
            println("Modelling band $(i)")
            # Run the model to obtain reflectance matrix
            #R = rt_run(model, i_band=i)[1];
            #R = CoreRT.rt_run_test(CoreRT.noRS(), model, i)[1]
            # Re-interpolate I from ν_grid to new grid/resolution
            λ_grid = reverse(1e4 ./ parameters.spec_bands[i])
            res = 0.001e-3;
            off = 0.5e-3
            # Get sun:
            #@time sun_out = getSolar(parameters.spec_bands[i],Tsolar)
        
            RR = oco_sounding[ia].mueller[1]*R[i][ia,1,:] + 
                    oco_sounding[ia].mueller[2]*R[i][ia,2,:] + 
                    oco_sounding[ia].mueller[3]*R[i][ia,3,:];        
        
            # Apply Earth reflectance matrix 
            #earth_out = sun_out .* reverse(RR[:])
            earth_out = reverse(RR[:])
            
            # Re-interpolate I from ν_grid to new grid/resolution
            @time interp_I = LinearInterpolation(λ_grid, earth_out);        
            #The following holds only for the use of a single ground pixel (because the ils changes with ground pixel index)
            wl = oco_sounding[ia].ils[i].ν_out[1]-off:res:oco_sounding[ia].ils[i].ν_out[end]+off;
            @show wl[1],wl[end], λ_grid[1],λ_grid[end]
            I_wl = interp_I(wl);

            # Convolve input spectrum with variable kernel
            @time I_conv = InstrumentOperator.conv_spectra(oco_sounding[ia].ils[i], wl, I_wl)
            off = oco_soundings[1].BandID[end][end] * (ia-1)
            y[oco_sounding[ia].BandID[i] .+ off ] = I_conv
            #if ia==1
            #    y[oco_sounding[ia].BandID[i]] = I_conv
            #else
            #    # not sure this will work!
            #    y[oco_sounding[ia].BandID[i]] = vcat(y[oco_sounding.BandID[i]], I_conv)
            #end
        end
    end
end