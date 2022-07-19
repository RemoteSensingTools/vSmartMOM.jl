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

##

# Load YAML files into parameter struct
parameters = 
    parameters_from_yaml("test/test_parameters/FraunhoferMockParameters.yaml");
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
    RS_type.F₀[1,i] = F₀[i];
end 

#solar_transmission_from_file(file_name::String) = readdlm(file_name)
# For now, convert these special cases to the right array type:
#aType = array_type(model.params.architecture);
#RS_type.Z⁺⁺_λ₁λ₀ = aType(RS_type.Z⁺⁺_λ₁λ₀);
#RS_type.ϖ_λ₁λ₀   = aType(RS_type.ϖ_λ₁λ₀);
#RS_type.i_λ₁λ₀   = aType(RS_type.i_λ₁λ₀)
#RS_type.Z⁻⁺_λ₁λ₀ = aType(RS_type.Z⁻⁺_λ₁λ₀);
# Add something here, that computes ALL the OP needed for the Raman case.
#modelRS = ...
#=
R, T, ieR, ieT = CoreRT.rt_run_ms(RS_type,
    model.params.polarization_type,
    model.obs_geom,
    [0,3],
    model.τ_rayl[1], 
    model.τ_aer[1], 
    model.quad_points,
    model.params.max_m,
    model.aerosol_optics[1],
    model.greek_rayleigh,
    model.τ_abs[1],
    model.params.brdf[1],
    model.params.architecture);
#R = vSmartMOM.rt_run(model, i_band=1)

RnoRS, TnoRS, _, _ = CoreRT.rt_run_ms(noRS(), 
            model.params.polarization_type,
            model.obs_geom,
            [0,3],
            model.τ_rayl[1], 
            model.τ_aer[1], 
            model.quad_points,
            model.params.max_m,
            model.aerosol_optics[1],
            model.greek_rayleigh,
            model.τ_abs[1],
            model.params.brdf[1],
            model.params.arhitecture);
   =#       
#RnoRS_ms, TnoRS_ms, _, _ = CoreRT.rt_run_test_ms(noRS([0.0],[1.0], Any[],[1]),[0,3],model,iBand);
#RnoRS_test, TnoRS_test, _, _ = CoreRT.rt_run_test_ms(noRS,[0,3],model,iBand);
#RnoRS, TnoRS, _, _ = CoreRT.rt_run_test(noRS([0.0],[1.0], Any[],[1]),model,iBand);
#sR_ms, T_ms, ieR_ms, ieT_ms = CoreRT.rt_run_test_ms(RS_type,[0,3], model,iBand);
#make F₀ the concatenated solar function across all bands contained in iBand
# R_ms, T_ms, ieR_ms, ieT_ms = 
    # CoreRT.rt_run_test_ms(RS_type,model,iBand);

R, T, ieR, ieT = CoreRT.rt_run_test(RS_type,model,iBand);


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

RnoRS, TnoRS, _, _ = CoreRT.rt_run_test(RS_type,model,iBand);
#=
R_test, T_test, ieR_test, ieT_test = CoreRT.rt_run_test(RS_type,model,iBand);

# You can now run multiple bands like this (list at the end of band!)
RnoRS_test, TnoRS_test, _, _ = CoreRT.rt_run_test(vSmartMOM.noRS(),model,[1,1]);


R_test, T_test, ieR_test, ieT_test = CoreRT.rt_run_test(RS_type,model,1);
=#

#===Convolution of hires spectral simulations to instrument grid===#
#x = (1e7/415):0.3:(1e7/385)
#ν = (1e7/460):0.3:(1e7/420)
x = -40:0.3:40
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
p1 = plot(1e7./ν, RnoRS[1,1,:].*convfct, linecolor=:black)
p1 = plot!(1e7./ν, I_conv_noRS.*convfct, linewidth = 3, linecolor=:red)
p2 = plot(1e7./ν, (R[1,1,:].+ieR[1,1,:]).*convfct, linecolor=:black, ylabel="[mW/m²/str/nm]")
p2 = plot!(1e7./ν, (I_conv.+ieI_conv).*convfct, linewidth = 3, linecolor=:red)
p3 = plot(1e7./ν, (R[1,1,:].-RnoRS[1,1,:].+ieR[1,1,:]).*convfct, linecolor=:black, xlabel = "λ [nm]")
p3 = plot!(1e7./ν, (I_conv.-I_conv_noRS.+ieI_conv).*convfct, linewidth = 3, linecolor=:red)

q1 = plot(1e7./ν, RnoRS[1,2,:].*convfct, linecolor=:black)
q1 = plot!(1e7./ν, Q_conv_noRS.*convfct, linewidth = 3, linecolor=:red)
q2 = plot(1e7./ν, (R[1,2,:].+ieR[1,2,:]).*convfct, linecolor=:black)
q2 = plot!(1e7./ν, (Q_conv.+ieQ_conv).*convfct, linewidth = 3, linecolor=:red)
q3 = plot(1e7./ν, (R[1,2,:].-RnoRS[1,2,:].+ieR[1,2,:]).*convfct, linecolor=:black, xlabel = "λ [nm]")
q3 = plot!(1e7./ν, (Q_conv.-Q_conv_noRS.+ieQ_conv).*convfct, linewidth = 3, linecolor=:red)
plot(p1, q1, p2, q2, p3, q3, layout = l, legend = false, title = ["I₀ (no RS)" "Q₀ (no RS)" "I₁ (with RS)" "Q₁ (with RS)" "I₁-I₀" "Q₁-Q₀"], titlefont = font(10))
savefig("RingEffect.png")

l = @layout [a1 a2]
p1 = plot(1e7./ν, ieR[1,1,:]./R[1,1,:].*convfct, linecolor=:black)
p1 = plot!(1e7./ν, ieI_conv./I_conv.*convfct, linewidth = 2, linecolor=:red, xlabel = "λ [nm]")
q1 = plot(1e7./ν, ieR[1,2,:]./R[1,2,:].*convfct, linecolor=:black)
q1 = plot!(1e7./ν, ieQ_conv./Q_conv.*convfct, linewidth = 2, linecolor=:red, xlabel = "λ [nm]")
plot(p1, q1, layout = l, legend = false, title = ["Iᵢ/Iₑ" "Qᵢ/Qₑ"], titlefont = font(10))
savefig("RingSpectrum.png")