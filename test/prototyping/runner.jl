# Runner is used to set AD fields as duals
function runner!(y, x, parameters=parameters, oco_sounding= oco_soundings, Tsolar = Tsolar_interp)
    Nangles = length(oco_sounding);
    # Set parameters fields as the dual numbers
    parameters.brdf =  [CoreRT.LambertianSurfaceLegendre([x[1],x[3],x[4]]),
                        CoreRT.LambertianSurfaceLegendre([x[7],x[8],x[5]]),
                        CoreRT.LambertianSurfaceLegendre([x[9],x[10],x[6]])];

    parameters.scattering_params.rt_aerosols[1].τ_ref = exp(x[2]);
    parameters.scattering_params.rt_aerosols[1].p₀    = x[20]; #800.0; #x[4]
    parameters.scattering_params.rt_aerosols[1].σp    = x[21];
    parameters.scattering_params.rt_aerosols[1].aerosol.size_distribution = LogNormal(x[18], x[19]); #x[4]
    parameters.scattering_params.rt_aerosols[1].aerosol.nᵣ = x[16]
    parameters.scattering_params.rt_aerosols[1].aerosol.nᵣ = x[17]
    #parameters.p   = oco_sounding.p_half
    #parameters.q   = oco_sounding.q 
    #parameters.T   = oco_sounding.T# .+ 1.0 #.+ x[15]
    #parameters.sza = oco_sounding.sza
    #parameters.vza = [oco_sounding.vza]
    for ia=1:Nangles
        if ia==1
            parameters.p   = oco_sounding[ia].p_half
            parameters.q   = oco_sounding[ia].q 
            parameters.T   = oco_sounding[ia].T# .+ 1.0 #.+ x[15]
            parameters.sza = oco_sounding[ia].sza
            parameters.vza = [oco_sounding[ia].vza]   
            parameters.vaz = [oco_sounding[ia].raa]
        else
            parameters.vza = vcat(parameters.vza, oco_sounding[ia].vza)
            parameters.vaz = vcat(parameters.vaz, oco_sounding[ia].raa)
        end
    end
    parameters.absorption_params.vmr["H2O"] = [parameters.q[1:65]*x[11] * 1.8; 
                                               parameters.q[66:end]*x[15] * 1.8];
    a1 = zeros(7) .+ x[12]
    a2 = zeros(7) .+ x[13]
    a3 = zeros(6) .+ x[14]
    parameters.absorption_params.vmr["CO2"] = [a1; a2; a3]
    model = model_from_parameters(parameters);
    for ia = 1:Nangles
        for i = 1:length(oco_sounding[ia].BandID)
            println("Modelling band $(i)")
            # Run the model to obtain reflectance matrix
            #R = rt_run(model, i_band=i)[1];
            R = CoreRT.rt_run_test(CoreRT.noRS(), model, i)[1]
            # Re-interpolate I from ν_grid to new grid/resolution
            λ_grid = reverse(1e4 ./ parameters.spec_bands[i])
            res = 0.001e-3;
            off = 0.5e-3
            # Get sun:
            @time sun_out = getSolar(parameters.spec_bands[i],Tsolar)
        
            RR = oco_sounding[ia].mueller[1]*R[ia,1,:] + 
                    oco_sounding[ia].mueller[2]*R[ia,2,:] + 
                    oco_sounding[ia].mueller[3]*R[ia,3,:];        
        
            # Apply Earth reflectance matrix 
            earth_out = sun_out .* reverse(RR[:])
        
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