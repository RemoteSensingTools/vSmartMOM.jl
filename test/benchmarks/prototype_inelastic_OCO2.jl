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

# Load OCO Data: 
# File names:
L1File   = "/net/fluo/data1/group/oco2/L1bSc/oco2_L1bScND_26780a_190715_B10003r_200429212407.h5"
metFile  = "/net/fluo/data1/group/oco2/L2Met/oco2_L2MetND_26780a_190715_B10003r_200429212406.h5"
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
model      = model_from_parameters(parameters);


iBand = 1
FT = Float64
I_wl=[];
I_conv=[];
I_conv_RS=[];
#### Compute all Raman properties
for iBand in length(bands)
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
                n_Raman=0);

    #vSmartMOM.get_greek_raman!(RS_type, n2, o2);
    # now compute other optical parameters for Raman:
    #fscattRayl = ...
    #RS_type = vSmartMOM.RRS(...)
    # Compute Raman SSA properties:
    CoreRT.getRamanSSProp!(RS_type, 1e7/ν̃, ν);

    # For now, convert these special cases to the right array type:
    #aType = array_type(model.params.architecture);
    #RS_type.Z⁺⁺_λ₁λ₀ = aType(RS_type.Z⁺⁺_λ₁λ₀);
    #RS_type.ϖ_λ₁λ₀   = aType(RS_type.ϖ_λ₁λ₀);
    #RS_type.i_λ₁λ₀   = aType(RS_type.i_λ₁λ₀)
    #RS_type.Z⁻⁺_λ₁λ₀ = aType(RS_type.Z⁻⁺_λ₁λ₀);
    # Add something here, that computes ALL the OP needed for the Raman case.
    #modelRS = ...
    R, T, ieR, ieT =CoreRT.rt_run_test(RS_type, model, iBand)
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

    
    # Full measurement:
    RR = oco_sounding.mueller[1]*(R[iBand,1,:]+ieR[iBand,1,:]) + 
         oco_sounding.mueller[2]*(R[iBand,2,:]+ieR[iBand,2,:]) + 
         oco_sounding.mueller[2]*(R[iBand,3,:]+ieR[iBand,3,:]);

    earth_out = reverse(RR[:])
    # Re-interpolate I from ν_grid to new grid/resolution
    λ_grid = reverse(1e4 ./ parameters.spec_bands[iBand])
    @time interp_I = LinearInterpolation(λ_grid, earth_out);

    res = 0.001e-3;
    off = 0.5e-3
    wl = oco_sounding.ils[iBand].ν_out[1]-off:res:oco_sounding.ils[iBand].ν_out[end]+off;
    @show wl[1],wl[end], λ_grid[1],λ_grid[end]
    tmp1 = interp_I(wl);

    # Convolve input spectrum with variable kernel
    @time tmp2 = InstrumentOperator.conv_spectra(oco_sounding.ils[iBand], wl, tmp1)
    #R = vSmartMOM.rt_run(model, i_band=1)
    push!(I_conv_RS, tmp2);

    RnoRS, TnoRS, _, _ = rt_run(noRS(),
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
    push!(I_conv, tmp2);
end
l = @layout [a b]
i=1
p1 = plot(oco_sounding.SpectralGrid*1e3, I_conv[i], label="Elastic")
plot!(oco_sounding.SpectralGrid*1e3, I_conv_RS[i], label="With RRS")
savefig("O2_RRS_SZA30_VZA45_VAz0_aer0p0_alb0.png")
p2 = plot(oco_sounding.SpectralGrid*1e3, 
    (I_conv_RS[i] - I_conv[i])/maximum(I_conv[i]) * 100, 
    ylims=(-5,5), 
    label="(RRS-NoRS)/max(noRS)*100")
p3 = plot(oco_sounding.SpectralGrid*1e3, 
    (I_conv_RS[i] .- I_conv[i])./(I_conv[i]) * 100, 
    ylims=(-5,5), 
    label="(RRS-NoRS)/noRS*100")
plot(p2,p3,layout=l)
savefig("O2_RRS_impact_SZA30_VZA45_VAz0_aer0p0_alb0.png")
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
