using Revise
using vSmartMOM, vSmartMOM.CoreRT, vSmartMOM.SolarModel
using vSmartMOM.InelasticScattering
using Statistics
using Interpolations
using InstrumentOperator #for convolution of hires spectrum to instrument grid
using ImageFiltering
using Distributions
using Plots
using CUDA
device!(1)
using Revise
using DelimitedFiles

Fraunhofer=true
λ₀ = 400. #745. #400. #nm
#iBand = 1
n2,o2 = InelasticScattering.getRamanAtmoConstants(1.7/λ₀, 300.);
FT = Float64
RS_type = InelasticScattering.VS_0to1_plus{FT}(
                   n2=n2,
                   o2=o2);
# Load YAML files into parameter struct
parameters = parameters_from_yaml("test/test_parameters/O2ParametersVS.yaml");
# Create model struct (precomputes optical properties) from parameters
model      = model_from_parameters(RS_type, λ₀, parameters);
τ_abs = model.τ_abs

if Fraunhofer
    T_sun = 5777. # KRS 
    Tsolar = solar_transmission_from_file("/home/sanghavi/code/github/vSmartMOM.jl/src/SolarModel/solar.out")
    Tsolar_interp = LinearInterpolation(Tsolar[4:end, 1], Tsolar[4:end, 2])
    RS_type.F₀ = zeros(model.params.polarization_type.n, 
                        sum(length(RS_type.grid_in[i]) for i in 1:length(RS_type.iBand)))
    local i_offset = 0
    for iB=1:length(RS_type.iBand)
        ν = collect(RS_type.grid_in[iB])
        P = planck_spectrum_wn(T_sun, ν)    
        F₀ = zeros(length(ν));
        
        for i=1:length(ν)
            sol_trans = Tsolar_interp(ν[i]);
            F₀[i] = sol_trans * P[i];
            RS_type.F₀[1,i+i_offset] = F₀[i];
        end 
        i_offset += length(ν)
    end
else
    RS_type.F₀ = zeros(model.params.polarization_type.n, 
        sum(length(RS_type.grid_in[i]) for i in 1:length(RS_type.iBand)))
    local i_offset = 0
    for iB=1:length(RS_type.iBand)
        ν = collect(RS_type.grid_in[iB])

        for i=1:length(ν)
            RS_type.F₀[1,i+i_offset] = 1.0;
        end 
        i_offset += length(ν)
    end
end

R, T, ieR, ieT = CoreRT.rt_run_test(RS_type,
           model, RS_type.iBand);

############# VS with polarized laser light ############
if Fraunhofer
    T_sun = 5777. # KRS 
    Tsolar = solar_transmission_from_file("/home/sanghavi/code/github/vSmartMOM.jl/src/SolarModel/solar.out")
    Tsolar_interp = LinearInterpolation(Tsolar[4:end, 1], Tsolar[4:end, 2])
    RS_type.F₀ = zeros(model.params.polarization_type.n, 
                        sum(length(RS_type.grid_in[i]) for i in 1:length(RS_type.iBand)))
    local i_offset = 0
    for iB=1:length(RS_type.iBand)
        ν = collect(RS_type.grid_in[iB])
        P = planck_spectrum_wn(T_sun, ν)    
        F₀ = zeros(length(ν));
        
        for i=1:length(ν)
            sol_trans = Tsolar_interp(ν[i]);
            F₀[i] = sol_trans * P[i];
            RS_type.F₀[1,i+i_offset] = F₀[i];
            RS_type.F₀[2,i+i_offset] = -F₀[i];
        end 
        i_offset += length(ν)
    end
else
    RS_type.F₀ = zeros(model.params.polarization_type.n, 
        sum(length(RS_type.grid_in[i]) for i in 1:length(RS_type.iBand)))
    local i_offset = 0
    for iB=1:length(RS_type.iBand)
        ν = collect(RS_type.grid_in[iB])

        for i=1:length(ν)
            RS_type.F₀[1,i+i_offset] = 1.0;
            RS_type.F₀[2,i+i_offset] = -1.0;
        end 
        i_offset += length(ν)
    end
end

Rl, Tl, ieRl, ieTl = CoreRT.rt_run_test(RS_type,
           model, RS_type.iBand);            

############# RRS ##############
parameters2 = 
       parameters_from_yaml("test/test_parameters/FraunhoferMockParameters.yaml");
#parameters = parameters_from_yaml("test/test_parameters/O2Parameters2.yaml");
# Create model struct (precomputes optical properties) from parameters
model2      = model_from_parameters(parameters2);

τ_abs2 = model2.τ_abs

iBand2 = 1
#### Compute all Raman properties
ν = model2.params.spec_bands[iBand2]
ν̃ = mean(ν);
# Find central reference index for RRS:
i_ref = argmin(abs.(ν .- ν̃))
# TODO_VS: λ_vs_in (get input)
# TODO_VS: ν_vs_in (convert to wavenumbers)
# Effective temperature for Raman calculations
effT = (model2.profile.vcd_dry' * model2.profile.T) / sum(model2.profile.vcd_dry);
# Define RS type
# Compute N2 and O2
n2,o2 = InelasticScattering.getRamanAtmoConstants(ν̃,effT);
#greek_raman = get_greek_raman(RS_type, n2, o2);
RS_type2 = InelasticScattering.RRS(
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
CoreRT.getRamanSSProp!(RS_type2, 1e7/ν̃, ν);
T_sun = 5777. # K
#P = planck_spectrum_wn(T_sun, ν)
#Tsolar = solar_transmission_from_file("/home/sanghavi/code/github/vSmartMOM.jl/src/SolarModel/solar.out")
#Tsolar_interp = LinearInterpolation(Tsolar[4:end, 1], Tsolar[4:end, 2])
#F₀2 = zeros(length(ν));
RS_type2.F₀ = zeros(model2.params.polarization_type.n, length(ν))
#for i=1:length(ν)
#sol_trans = Tsolar_interp(ν[i]);   
#F₀[i] = sol_trans * P[i];
RRSidx = findmin(abs.(1e7./ν .- 400.))[2]
#RRSidx = findmin(abs.(1e7./ν .- 745.))[2]
#RS_type2.F₀[1,:] .= 1.
RS_type2.F₀[1,RRSidx] = RS_type.F₀[1]
#end 
R_RRS, T_RRS, ieR_RRS, ieT_RRS = CoreRT.rt_run_test(RS_type2,model2,iBand2);
############# RRS with polarized laser light ############
RS_type2.F₀ = zeros(model2.params.polarization_type.n, length(ν))
#for i=1:length(ν)
#sol_trans = Tsolar_interp(ν[i]);   
#F₀[i] = sol_trans * P[i];
RRSidx = findmin(abs.(1e7./ν .- 400.))[2]
#RRSidx = findmin(abs.(1e7./ν .- 745.))[2]
#RS_type2.F₀[1,:] .= 1.
RS_type2.F₀[1,RRSidx] = RS_type.F₀[1]
RS_type2.F₀[2,RRSidx] = -RS_type.F₀[1]
R_RRSl, T_RRSl, ieR_RRSl, ieT_RRSl = CoreRT.rt_run_test(RS_type2,model2,iBand2);

########## RRS w/ Fraunhofer ##########
T_sun = 5777. # KRS 
Tsolar = solar_transmission_from_file("/home/sanghavi/code/github/vSmartMOM.jl/src/SolarModel/solar.out")
Tsolar_interp = LinearInterpolation(Tsolar[4:end, 1], Tsolar[4:end, 2])
RS_type2.F₀ = zeros(model2.params.polarization_type.n, 
                        sum(length(model2.params.spec_bands[i]) for i in 1:length(RS_type2.iBand)))
#i_offset = 0
for iB=1:length(RS_type2.iBand)
    #ν = collect(RS_type2.grid_in[iB])
    ν = model2.params.spec_bands[iB]
    P = planck_spectrum_wn(T_sun, ν)    
    F₀ = zeros(length(ν));
        
    for i=1:length(ν)
        sol_trans = Tsolar_interp(ν[i]);
        F₀[i] = sol_trans * P[i];
        RS_type2.F₀[1,i] = F₀[i];
    end 
    #i_offset += length(ν)
end
R_RRS2, T_RRS2, ieR_RRS2, ieT_RRS2 = CoreRT.rt_run_test(RS_type2,model2,iBand2);


############# RRS with polarized laser light ############
RS_type2.F₀ = zeros(model2.params.polarization_type.n, 
                        sum(length(model2.params.spec_bands[i]) for i in 1:length(RS_type2.iBand)))
#i_offset = 0
for iB=1:length(RS_type2.iBand)
    #ν = collect(RS_type2.grid_in[iB])
    ν = model2.params.spec_bands[iB]
    P = planck_spectrum_wn(T_sun, ν)    
    F₀ = zeros(length(ν));
        
    for i=1:length(ν)
        sol_trans = Tsolar_interp(ν[i]);
        F₀[i] = sol_trans * P[i];
        RS_type2.F₀[1,i] = F₀[i];
        RS_type2.F₀[2,i] = -F₀[i];
    end 
    #i_offset += length(ν)
end
R_RRS2l, T_RRS2l, ieR_RRS2l, ieT_RRS2l = CoreRT.rt_run_test(RS_type2,model2,iBand2);

############## Now starting no-ozone computations ########################
Fraunhofer=true
λ₀ = 400. #745. #400. #nm
#iBand = 1
n2,o2 = InelasticScattering.getRamanAtmoConstants(1.7/λ₀, 300.);
FT = Float64
RS_type = InelasticScattering.VS_0to1_plus{FT}(
                   n2=n2,
                   o2=o2);
# Load YAML files into parameter struct
parameters = parameters_from_yaml("test/test_parameters/O2ParametersVS0.yaml");
# Create model struct (precomputes optical properties) from parameters
model      = model_from_parameters(RS_type, λ₀, parameters);
τ_abs0 = model.τ_abs

if Fraunhofer
    T_sun = 5777. # KRS 
    Tsolar = solar_transmission_from_file("/home/sanghavi/code/github/vSmartMOM.jl/src/SolarModel/solar.out")
    Tsolar_interp = LinearInterpolation(Tsolar[4:end, 1], Tsolar[4:end, 2])
    RS_type.F₀ = zeros(model.params.polarization_type.n, 
                        sum(length(RS_type.grid_in[i]) for i in 1:length(RS_type.iBand)))
    local i_offset = 0
    for iB=1:length(RS_type.iBand)
        ν = collect(RS_type.grid_in[iB])
        P = planck_spectrum_wn(T_sun, ν)    
        F₀ = zeros(length(ν));
        
        for i=1:length(ν)
            sol_trans = Tsolar_interp(ν[i]);
            F₀[i] = sol_trans * P[i];
            RS_type.F₀[1,i+i_offset] = F₀[i];
        end 
        i_offset += length(ν)
    end
else
    RS_type.F₀ = zeros(model.params.polarization_type.n, 
        sum(length(RS_type.grid_in[i]) for i in 1:length(RS_type.iBand)))
    local i_offset = 0
    for iB=1:length(RS_type.iBand)
        ν = collect(RS_type.grid_in[iB])

        for i=1:length(ν)
            RS_type.F₀[1,i+i_offset] = 1.0;
        end 
        i_offset += length(ν)
    end
end

R0, T0, ieR0, ieT0 = CoreRT.rt_run_test(RS_type,
           model, RS_type.iBand);
############# VS with polarized laser light ############
if Fraunhofer
    T_sun = 5777. # KRS 
    Tsolar = solar_transmission_from_file("/home/sanghavi/code/github/vSmartMOM.jl/src/SolarModel/solar.out")
    Tsolar_interp = LinearInterpolation(Tsolar[4:end, 1], Tsolar[4:end, 2])
    RS_type.F₀ = zeros(model.params.polarization_type.n, 
                        sum(length(RS_type.grid_in[i]) for i in 1:length(RS_type.iBand)))
    local i_offset = 0
    for iB=1:length(RS_type.iBand)
        ν = collect(RS_type.grid_in[iB])
        P = planck_spectrum_wn(T_sun, ν)    
        F₀ = zeros(length(ν));
        
        for i=1:length(ν)
            sol_trans = Tsolar_interp(ν[i]);
            F₀[i] = sol_trans * P[i];
            RS_type.F₀[1,i+i_offset] = F₀[i];
            RS_type.F₀[2,i+i_offset] = -F₀[i];
        end 
        i_offset += length(ν)
    end
else
    RS_type.F₀ = zeros(model.params.polarization_type.n, 
        sum(length(RS_type.grid_in[i]) for i in 1:length(RS_type.iBand)))
    local i_offset = 0
    for iB=1:length(RS_type.iBand)
        ν = collect(RS_type.grid_in[iB])

        for i=1:length(ν)
            RS_type.F₀[1,i+i_offset] = 1.0;
            RS_type.F₀[2,i+i_offset] = -1.0;
        end 
        i_offset += length(ν)
    end
end

R0l, T0l, ieR0l, ieT0l = CoreRT.rt_run_test(RS_type,
           model, RS_type.iBand);

############# RRS ##############
parameters2 = 
       parameters_from_yaml("test/test_parameters/FraunhoferMockParameters0.yaml");
#parameters = parameters_from_yaml("test/test_parameters/O2Parameters2.yaml");
# Create model struct (precomputes optical properties) from parameters
model2      = model_from_parameters(parameters2);

τ_abs02 = model2.τ_abs

iBand2 = 1
#### Compute all Raman properties
ν = model2.params.spec_bands[iBand2]
ν̃ = mean(ν);
# Find central reference index for RRS:
i_ref = argmin(abs.(ν .- ν̃))
# TODO_VS: λ_vs_in (get input)
# TODO_VS: ν_vs_in (convert to wavenumbers)
# Effective temperature for Raman calculations
effT = (model2.profile.vcd_dry' * model2.profile.T) / sum(model2.profile.vcd_dry);
# Define RS type
# Compute N2 and O2
n2,o2 = InelasticScattering.getRamanAtmoConstants(ν̃,effT);
#greek_raman = get_greek_raman(RS_type, n2, o2);
RS_type2 = InelasticScattering.RRS(
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
CoreRT.getRamanSSProp!(RS_type2, 1e7/ν̃, ν);
T_sun = 5777. # K
#P = planck_spectrum_wn(T_sun, ν)
#Tsolar = solar_transmission_from_file("/home/sanghavi/code/github/vSmartMOM.jl/src/SolarModel/solar.out")
#Tsolar_interp = LinearInterpolation(Tsolar[4:end, 1], Tsolar[4:end, 2])
#F₀2 = zeros(length(ν));
RS_type2.F₀ = zeros(model2.params.polarization_type.n, length(ν))
#for i=1:length(ν)
#sol_trans = Tsolar_interp(ν[i]);   
#F₀[i] = sol_trans * P[i];
RRSidx = findmin(abs.(1e7./ν .- 400.))[2]
#RRSidx = findmin(abs.(1e7./ν .- 745.))[2]
#RS_type2.F₀[1,:] .= 1.
RS_type2.F₀[1,RRSidx] = RS_type.F₀[1]
#end 
R0_RRS, T0_RRS, ieR0_RRS, ieT0_RRS = CoreRT.rt_run_test(RS_type2,model2,iBand2);
############# RRS with polarized laser light ############
RS_type2.F₀ = zeros(model2.params.polarization_type.n, length(ν))
#for i=1:length(ν)
#sol_trans = Tsolar_interp(ν[i]);   
#F₀[i] = sol_trans * P[i];
RRSidx = findmin(abs.(1e7./ν .- 400.))[2]
#RRSidx = findmin(abs.(1e7./ν .- 745.))[2]
#RS_type2.F₀[1,:] .= 1.
RS_type2.F₀[1,RRSidx] = RS_type.F₀[1]
RS_type2.F₀[2,RRSidx] = -RS_type.F₀[1]
#end 
R0_RRSl, T0_RRSl, ieR0_RRSl, ieT0_RRSl = CoreRT.rt_run_test(RS_type2,model2,iBand2);

########## RRS w/ Fraunhofer ##########
T_sun = 5777. # KRS 
Tsolar = solar_transmission_from_file("/home/sanghavi/code/github/vSmartMOM.jl/src/SolarModel/solar.out")
Tsolar_interp = LinearInterpolation(Tsolar[4:end, 1], Tsolar[4:end, 2])
RS_type2.F₀ = zeros(model2.params.polarization_type.n, 
                        sum(length(model2.params.spec_bands[i]) for i in 1:length(RS_type2.iBand)))
#i_offset = 0
for iB=1:length(RS_type2.iBand)
    #ν = collect(RS_type2.grid_in[iB])
    ν = model2.params.spec_bands[iB]
    P = planck_spectrum_wn(T_sun, ν)    
    F₀ = zeros(length(ν));
        
    for i=1:length(ν)
        sol_trans = Tsolar_interp(ν[i]);
        F₀[i] = sol_trans * P[i];
        RS_type2.F₀[1,i] = F₀[i];
    end 
    #i_offset += length(ν)
end
R0_RRS2, T0_RRS2, ieR0_RRS2, ieT0_RRS2 = CoreRT.rt_run_test(RS_type2,model2,iBand2);

############# RRS with polarized laser light ############
RS_type2.F₀ = zeros(model2.params.polarization_type.n, 
                        sum(length(model2.params.spec_bands[i]) for i in 1:length(RS_type2.iBand)))
#i_offset = 0
for iB=1:length(RS_type2.iBand)
    #ν = collect(RS_type2.grid_in[iB])
    ν = model2.params.spec_bands[iB]
    P = planck_spectrum_wn(T_sun, ν)    
    F₀ = zeros(length(ν));
        
    for i=1:length(ν)
        sol_trans = Tsolar_interp(ν[i]);
        F₀[i] = sol_trans * P[i];
        RS_type2.F₀[1,i] = F₀[i];
        RS_type2.F₀[2,i] = -F₀[i];
    end 
    #i_offset += length(ν)
end
R0_RRS2l, T0_RRS2l, ieR0_RRS2l, ieT0_RRS2l = CoreRT.rt_run_test(RS_type2,model2,iBand2);

################### End of no-ozone computations ###################

########## Plotting ###########
ν₀ = RS_type.grid_in[1]
ν₁ = RS_type.grid_in[2]
ν₂ = RS_type.grid_in[3]

R₀₀ = R[1,1,RS_type.bandSpecLim[1]]
R₀₁ = R[1,1,RS_type.bandSpecLim[2]]
R₀₂ = R[1,1,RS_type.bandSpecLim[3]]

ieR₀₀ = ieR[1,1,RS_type.bandSpecLim[1]]
ieR₀₁ = ieR[1,1,RS_type.bandSpecLim[2]]
ieR₀₂ = ieR[1,1,RS_type.bandSpecLim[3]]

Q₀₀ = R[1,2,RS_type.bandSpecLim[1]]
Q₀₁ = R[1,2,RS_type.bandSpecLim[2]]
Q₀₂ = R[1,2,RS_type.bandSpecLim[3]]

ieQ₀₀ = ieR[1,2,RS_type.bandSpecLim[1]]
ieQ₀₁ = ieR[1,2,RS_type.bandSpecLim[2]]
ieQ₀₂ = ieR[1,2,RS_type.bandSpecLim[3]]

# with polarized laser light
lR₀₀ = Rl[1,1,RS_type.bandSpecLim[1]]
lR₀₁ = Rl[1,1,RS_type.bandSpecLim[2]]
lR₀₂ = Rl[1,1,RS_type.bandSpecLim[3]]

lieR₀₀ = ieRl[1,1,RS_type.bandSpecLim[1]]
lieR₀₁ = ieRl[1,1,RS_type.bandSpecLim[2]]
lieR₀₂ = ieRl[1,1,RS_type.bandSpecLim[3]]

lQ₀₀ = Rl[1,2,RS_type.bandSpecLim[1]]
lQ₀₁ = Rl[1,2,RS_type.bandSpecLim[2]]
lQ₀₂ = Rl[1,2,RS_type.bandSpecLim[3]]

lieQ₀₀ = ieRl[1,2,RS_type.bandSpecLim[1]]
lieQ₀₁ = ieRl[1,2,RS_type.bandSpecLim[2]]
lieQ₀₂ = ieRl[1,2,RS_type.bandSpecLim[3]]
###### No ozone ######

R0₀₀ = R0[1,1,RS_type.bandSpecLim[1]]
R0₀₁ = R0[1,1,RS_type.bandSpecLim[2]]
R0₀₂ = R0[1,1,RS_type.bandSpecLim[3]]

ieR0₀₀ = ieR0[1,1,RS_type.bandSpecLim[1]]
ieR0₀₁ = ieR0[1,1,RS_type.bandSpecLim[2]]
ieR0₀₂ = ieR0[1,1,RS_type.bandSpecLim[3]]

Q0₀₀ = R0[1,2,RS_type.bandSpecLim[1]]
Q0₀₁ = R0[1,2,RS_type.bandSpecLim[2]]
Q0₀₂ = R0[1,2,RS_type.bandSpecLim[3]]

ieQ0₀₀ = ieR0[1,2,RS_type.bandSpecLim[1]]
ieQ0₀₁ = ieR0[1,2,RS_type.bandSpecLim[2]]
ieQ0₀₂ = ieR0[1,2,RS_type.bandSpecLim[3]]

# With polarized laser light
lR0₀₀ = R0l[1,1,RS_type.bandSpecLim[1]]
lR0₀₁ = R0l[1,1,RS_type.bandSpecLim[2]]
lR0₀₂ = R0l[1,1,RS_type.bandSpecLim[3]]

lieR0₀₀ = ieR0l[1,1,RS_type.bandSpecLim[1]]
lieR0₀₁ = ieR0l[1,1,RS_type.bandSpecLim[2]]
lieR0₀₂ = ieR0l[1,1,RS_type.bandSpecLim[3]]

lQ0₀₀ = R0l[1,2,RS_type.bandSpecLim[1]]
lQ0₀₁ = R0l[1,2,RS_type.bandSpecLim[2]]
lQ0₀₂ = R0l[1,2,RS_type.bandSpecLim[3]]

lieQ0₀₀ = ieR0l[1,2,RS_type.bandSpecLim[1]]
lieQ0₀₁ = ieR0l[1,2,RS_type.bandSpecLim[2]]
lieQ0₀₂ = ieR0l[1,2,RS_type.bandSpecLim[3]]

######################

convfct = 1e7./ν.^2
convfct0 = 1e7./ν₀.^2   # to convert radiance units from mW/m²-str-cm⁻¹ to mW/m²-str-nm
convfct1 = 1e7./ν₁.^2
convfct2 = 1e7./ν₂.^2

l = @layout [a1 a2 a3; b1 b2 b3]

p1 = plot(1e7./ν, R_RRS2[1,1,:].*convfct, linecolor=:grey, xlims=(395,405), xticks=395:5:405, ylabel="[mW/m²/str/nm]")
p1 = plot!(1e7./ν, R_RRS[1,1,:].*convfct, linecolor=:black, linewidth=2, xlims=(395,405), xticks=395:5:405, ylabel="[mW/m²/str/nm]")
p1 = plot!(1e7./ν₀, R₀₀.*convfct0, linecolor=:black, xlims=(395,405))
p3 = plot(1e7./ν₁, R₀₁.*convfct1, linecolor=:grey, xlims=(436,446), xticks=436:4:444) #, xticks=425:3:430)
p2 = plot(1e7./ν₂, R₀₂.*convfct2, linecolor=:grey, xticks=424:2.5:429)
q1 = plot(1e7./ν, ieR_RRS[1,1,:].*convfct, linecolor=:black, xlabel = "λ [nm]", xlims=(395,405), xticks=395:5:405, ylabel="[mW/m²/str/nm]")
q3 = plot(1e7./ν₁, ieR₀₁.*convfct1, linecolor=:black, xlabel = "λ [nm]", xlims=(436,446), xticks=436:4:444)
q2 = plot(1e7./ν₂, ieR₀₂.*convfct2, linecolor=:black, xlabel = "λ [nm]", xticks=424:2.5:429)
plot(p1, p2, p3, q1, q2, q3, layout = l, legend = false, title = ["Iₑ (incident)" "Iₑ (O₂)" "Iₑ (N₂)" " Iᵢ (RRS)" "Iᵢ (VRS+RVRS, O₂)" "Iᵢ (VRS+RVRS, N₂)"], titlefont = font(8))
savefig("RVRS_monochromatic.png")

p1 = plot(1e7./ν, R_RRS2l[1,1,:].*convfct, linecolor=:green, xlims=(395,405), xticks=395:5:405, ylabel="[mW/m²/str/nm]")
p1 = plot!(1e7./ν, R_RRSl[1,1,:].*convfct, linecolor=:blue, linewidth=2, xlims=(395,405), xticks=395:5:405, ylabel="[mW/m²/str/nm]")
p1 = plot!(1e7./ν₀, lR₀₀.*convfct0, linecolor=:blue, xlims=(395,405))
p3 = plot(1e7./ν₁, lR₀₁.*convfct1, linecolor=:green, xlims=(436,446), xticks=436:4:444) #, xticks=425:3:430)
p2 = plot(1e7./ν₂, lR₀₂.*convfct2, linecolor=:green, xticks=424:2.5:429)
q1 = plot(1e7./ν, ieR_RRSl[1,1,:].*convfct, linecolor=:blue, xlabel = "λ [nm]", xlims=(395,405), xticks=395:5:405, ylabel="[mW/m²/str/nm]")
q3 = plot(1e7./ν₁, lieR₀₁.*convfct1, linecolor=:blue, xlabel = "λ [nm]", xlims=(436,446), xticks=436:4:444)
q2 = plot(1e7./ν₂, lieR₀₂.*convfct2, linecolor=:blue, xlabel = "λ [nm]", xticks=424:2.5:429)

plot(p1, p2, p3, q1, q2, q3, layout = l, legend = false, title = ["Iₑ (incident)" "Iₑ (O₂)" "Iₑ (N₂)" " Iᵢ (RRS)" "Iᵢ (VRS+RVRS, O₂)" "Iᵢ (VRS+RVRS, N₂)"], titlefont = font(8))
savefig("RVRS_monochromatic_laser.png")

p1 = plot(1e7./ν, R_RRS2[1,1,:].*convfct, linecolor=:grey, xlims=(395,405), xticks=395:5:405, ylabel="[mW/m²/str/nm]")
p1 = plot!(1e7./ν, R_RRSl[1,1,:].*convfct, linecolor=:black, linewidth=2, xlims=(395,405), xticks=395:5:405, ylabel="[mW/m²/str/nm]")
#p1 = plot!(1e7./ν₀, lR₀₀.*convfct0, linecolor=:black, xlims=(395,405))
p3 = plot(1e7./ν₁, R₀₁.*convfct1, linecolor=:grey, xlims=(436,446), xticks=436:4:444) #, xticks=425:3:430)
p2 = plot(1e7./ν₂, R₀₂.*convfct2, linecolor=:grey, xticks=424:2.5:429)
q1 = plot(1e7./ν, ieR_RRSl[1,1,:].*convfct, linecolor=:black, xlabel = "λ [nm]", xlims=(395,405), xticks=395:5:405, ylabel="[mW/m²/str/nm]")
q3 = plot(1e7./ν₁, lieR₀₁.*convfct1, linecolor=:black, xlabel = "λ [nm]", xlims=(436,446), xticks=436:4:444)
q2 = plot(1e7./ν₂, lieR₀₂.*convfct2, linecolor=:black, xlabel = "λ [nm]", xticks=424:2.5:429)

plot(p1, p2, p3, q1, q2, q3, layout = l, legend = false, title = ["Iₑ (incident)" "Iₑ (O₂)" "Iₑ (N₂)" " Iᵢ (RRS)" "Iᵢ (VRS+RVRS, O₂)" "Iᵢ (VRS+RVRS, N₂)"], titlefont = font(8))
savefig("RVRS_monochromatic_sunlasercomp_vv.png")

l = @layout [a1 a2 a3; b1 b2 b3]

p1 = plot(1e7./ν, R_RRS2[1,2,:].*convfct, linecolor=:grey, xlims=(395,405), xticks=395:5:405, ylabel="[mW/m²/str/nm]")
p1 = plot!(1e7./ν, R_RRS[1,2,:].*convfct, linecolor=:black, linewidth=2, xlims=(395,405), xticks=395:5:405, ylabel="[mW/m²/str/nm]")
p1 = plot!(1e7./ν₀, Q₀₀.*convfct0, linecolor=:black, xlims=(395,405))
p3 = plot(1e7./ν₁, Q₀₁.*convfct1, linecolor=:grey, xlims=(436,446), xticks=436:4:444) #, xticks=425:3:430)
p2 = plot(1e7./ν₂, Q₀₂.*convfct2, linecolor=:grey, xticks=424:2.5:429)
q1 = plot(1e7./ν, ieR_RRS[1,2,:].*convfct, linecolor=:black, xlabel = "λ [nm]", xlims=(395,405), xticks=395:5:405, ylabel="[mW/m²/str/nm]")
q3 = plot(1e7./ν₁, ieQ₀₁.*convfct1, linecolor=:black, xlabel = "λ [nm]", xlims=(436,446), xticks=436:4:444)
q2 = plot(1e7./ν₂, ieQ₀₂.*convfct2, linecolor=:black, xlabel = "λ [nm]", xticks=424:2.5:429)
plot(p1, p2, p3, q1, q2, q3, layout = l, legend = false, title = ["Qₑ (incident)" "Qₑ (O₂)" "Qₑ (N₂)" " Qᵢ (RRS)" "Qᵢ (VRS+RVRS, O₂)" "Qᵢ (VRS+RVRS, N₂)"], titlefont = font(8))
savefig("RVRS_Qmonochromatic.png")

p1 = plot(1e7./ν, R_RRS2l[1,2,:].*convfct, linecolor=:green, xlims=(395,405), xticks=395:5:405, ylabel="[mW/m²/str/nm]")
p1 = plot!(1e7./ν, R_RRSl[1,2,:].*convfct, linecolor=:blue, linewidth=2, xlims=(395,405), xticks=395:5:405, ylabel="[mW/m²/str/nm]")
p1 = plot!(1e7./ν₀, lQ₀₀.*convfct0, linecolor=:blue, xlims=(395,405))
p3 = plot(1e7./ν₁, lQ₀₁.*convfct1, linecolor=:green, xlims=(436,446), xticks=436:4:444) #, xticks=425:3:430)
p2 = plot(1e7./ν₂, lQ₀₂.*convfct2, linecolor=:green, xticks=424:2.5:429)
q1 = plot(1e7./ν, ieR_RRSl[1,2,:].*convfct, linecolor=:blue, xlabel = "λ [nm]", xlims=(395,405), xticks=395:5:405, ylabel="[mW/m²/str/nm]")
q3 = plot(1e7./ν₁, lieQ₀₁.*convfct1, linecolor=:blue, xlabel = "λ [nm]", xlims=(436,446), xticks=436:4:444)
q2 = plot(1e7./ν₂, lieQ₀₂.*convfct2, linecolor=:blue, xlabel = "λ [nm]", xticks=424:2.5:429)

plot(p1, p2, p3, q1, q2, q3, layout = l, legend = false, title = ["Qₑ (incident)" "Qₑ (O₂)" "Qₑ (N₂)" " Qᵢ (RRS)" "Qᵢ (VRS+RVRS, O₂)" "Qᵢ (VRS+RVRS, N₂)"], titlefont = font(8))
savefig("RVRS_Qmonochromatic_laser.png")

p1 = plot(1e7./ν, R_RRS2[1,2,:].*convfct, linecolor=:grey, xlims=(395,405), xticks=395:5:405, ylabel="[mW/m²/str/nm]")
p1 = plot!(1e7./ν, R_RRSl[1,2,:].*convfct, linecolor=:black, linewidth=2, xlims=(395,405), xticks=395:5:405, ylabel="[mW/m²/str/nm]")
#p1 = plot!(1e7./ν₀, lQ₀₀.*convfct0, linecolor=:black, xlims=(395,405))
p3 = plot(1e7./ν₁, Q₀₁.*convfct1, linecolor=:grey, xlims=(436,446), xticks=436:4:444) #, xticks=425:3:430)
p2 = plot(1e7./ν₂, Q₀₂.*convfct2, linecolor=:grey, xticks=424:2.5:429)
q1 = plot(1e7./ν, ieR_RRSl[1,2,:].*convfct, linecolor=:black, xlabel = "λ [nm]", xlims=(395,405), xticks=395:5:405, ylabel="[mW/m²/str/nm]")
q3 = plot(1e7./ν₁, lieQ₀₁.*convfct1, linecolor=:black, xlabel = "λ [nm]", xlims=(436,446), xticks=436:4:444)
q2 = plot(1e7./ν₂, lieQ₀₂.*convfct2, linecolor=:black, xlabel = "λ [nm]", xticks=424:2.5:429)

plot(p1, p2, p3, q1, q2, q3, layout = l, legend = false, title = ["Qₑ (incident)" "Qₑ (O₂)" "Qₑ (N₂)" " Qᵢ (RRS)" "Qᵢ (VRS+RVRS, O₂)" "Qᵢ (VRS+RVRS, N₂)"], titlefont = font(8))
savefig("RVRS_Qmonochromatic_sunlasercomp_vv.png")

l = @layout [a1 a2 a3; b1 b2 b3]

p1 = plot(1e7./ν, (R_RRS2[1,1,:]-R0_RRS2[1,1,:]).*convfct, linecolor=:grey, xlims=(395,405), xticks=395:5:405, ylabel="[mW/m²/str/nm]")
p1 = plot!(1e7./ν, (R_RRS[1,1,:]-R0_RRS[1,1,:]).*convfct, linecolor=:black, linewidth=2, xlims=(395,405), xticks=395:5:405, ylabel="[mW/m²/str/nm]")
p1 = plot!(1e7./ν₀, (R₀₀-R0₀₀).*convfct0, linecolor=:black, xlims=(395,405))
p3 = plot(1e7./ν₁, (R₀₁-R0₀₁).*convfct1, linecolor=:grey) #, xticks=425:3:430)
p2 = plot(1e7./ν₂, (R₀₂-R0₀₂).*convfct2, linecolor=:grey, xticks=424:2.5:429)
q1 = plot(1e7./ν, (ieR_RRS[1,1,:]-ieR0_RRS[1,1,:]).*convfct, linecolor=:black, xlabel = "λ [nm]", xlims=(395,405), xticks=395:5:405, ylabel="[mW/m²/str/nm]")
q3 = plot(1e7./ν₁, (ieR₀₁-ieR0₀₁).*convfct1, linecolor=:black, xlabel = "λ [nm]")#
q2 = plot(1e7./ν₂, (ieR₀₂-ieR0₀₂).*convfct2, linecolor=:black, xlabel = "λ [nm]", xticks=424:2.5:429)
plot(p1, p2, p3, q1, q2, q3, layout = l, legend = false, title = ["ΔIₑ (incident)" "ΔIₑ (O₂)" "ΔIₑ (N₂)" " ΔIᵢ (RRS)" "ΔIᵢ (VRS+RVRS, O₂)" "ΔIᵢ (VRS+RVRS, N₂)"], titlefont = font(8))
savefig("DeltaRVRS_monochromatic.png")

p1 = plot(1e7./ν, (R_RRS2[1,1,:]-R0_RRS2[1,1,:]).*convfct, linecolor=:grey, xlims=(395,405), xticks=395:5:405, ylabel="[mW/m²/str/nm]")
p1 = plot!(1e7./ν, (R_RRSl[1,1,:]-R0_RRSl[1,1,:]).*convfct, linecolor=:black, linewidth=2, xlims=(395,405), xticks=395:5:405, ylabel="[mW/m²/str/nm]")
#p1 = plot!(1e7./ν₀, (R₀₀-R0₀₀).*convfct0, linecolor=:black, xlims=(395,405))
p3 = plot(1e7./ν₁, (R₀₁-R0₀₁).*convfct1, linecolor=:grey) #, xticks=425:3:430)
p2 = plot(1e7./ν₂, (R₀₂-R0₀₂).*convfct2, linecolor=:grey, xticks=424:2.5:429)
q1 = plot(1e7./ν, (ieR_RRSl[1,1,:]-ieR0_RRSl[1,1,:]).*convfct, linecolor=:black, xlabel = "λ [nm]", xlims=(395,405), xticks=395:5:405, ylabel="[mW/m²/str/nm]")
q3 = plot(1e7./ν₁, (lieR₀₁-lieR0₀₁).*convfct1, linecolor=:black, xlabel = "λ [nm]")#
q2 = plot(1e7./ν₂, (lieR₀₂-lieR0₀₂).*convfct2, linecolor=:black, xlabel = "λ [nm]", xticks=424:2.5:429)
plot(p1, p2, p3, q1, q2, q3, layout = l, legend = false, title = ["ΔIₑ (incident)" "ΔIₑ (O₂)" "ΔIₑ (N₂)" " ΔIᵢ (RRS)" "ΔIᵢ (VRS+RVRS, O₂)" "ΔIᵢ (VRS+RVRS, N₂)"], titlefont = font(8))
savefig("DeltaRVRS_monochromatic_sunlasercomp_vv.png")

l = @layout [a1 a2 a3; b1 b2 b3]

p1 = plot(1e7./ν, (R_RRS2[1,2,:]-R0_RRS2[1,2,:]).*convfct, linecolor=:grey, xlims=(395,405), xticks=395:5:405, ylabel="[mW/m²/str/nm]")
p1 = plot!(1e7./ν, (R_RRS[1,2,:]-R0_RRS[1,2,:]).*convfct, linecolor=:black, linewidth=2, xlims=(395,405), xticks=395:5:405, ylabel="[mW/m²/str/nm]")
p1 = plot!(1e7./ν₀, (Q₀₀-Q0₀₀).*convfct0, linecolor=:black, xlims=(395,405))
p3 = plot(1e7./ν₁, (Q₀₁-Q0₀₁).*convfct1, linecolor=:grey) #, xticks=425:3:430)
p2 = plot(1e7./ν₂, (Q₀₂-Q0₀₂).*convfct2, linecolor=:grey, xticks=424:2.5:429)
q1 = plot(1e7./ν, (ieR_RRS[1,2,:]-ieR0_RRS[1,2,:]).*convfct, linecolor=:black, xlabel = "λ [nm]", xlims=(395,405), xticks=395:5:405, ylabel="[mW/m²/str/nm]")
q3 = plot(1e7./ν₁, (ieQ₀₁-ieQ0₀₁).*convfct1, linecolor=:black, xlabel = "λ [nm]")#
q2 = plot(1e7./ν₂, (ieQ₀₂-ieQ0₀₂).*convfct2, linecolor=:black, xlabel = "λ [nm]", xticks=424:2.5:429)
plot(p1, p2, p3, q1, q2, q3, layout = l, legend = false, title = ["ΔQₑ (incident)" "ΔQₑ (O₂)" "ΔQₑ (N₂)" " ΔQᵢ (RRS)" "ΔQᵢ (VRS+RVRS, O₂)" "ΔQᵢ (VRS+RVRS, N₂)"], titlefont = font(8))
savefig("DeltaRVRS_Qmonochromatic.png")

l = @layout [a1 a2 a3; b1 b2 b3]

p1 = plot(1e7./ν, (R_RRS2[1,2,:]-R0_RRS2[1,2,:]).*convfct, linecolor=:grey, xlims=(395,405), xticks=395:5:405, ylabel="[mW/m²/str/nm]")
p1 = plot!(1e7./ν, (R_RRSl[1,2,:]-R0_RRSl[1,2,:]).*convfct, linecolor=:black, linewidth=2, xlims=(395,405), xticks=395:5:405, ylabel="[mW/m²/str/nm]")
#p1 = plot!(1e7./ν₀, (Q₀₀-Q0₀₀).*convfct0, linecolor=:black, xlims=(395,405))
p3 = plot(1e7./ν₁, (Q₀₁-Q0₀₁).*convfct1, linecolor=:grey) #, xticks=425:3:430)
p2 = plot(1e7./ν₂, (Q₀₂-Q0₀₂).*convfct2, linecolor=:grey, xticks=424:2.5:429)
q1 = plot(1e7./ν, (ieR_RRSl[1,2,:]-ieR0_RRSl[1,2,:]).*convfct, linecolor=:black, xlabel = "λ [nm]", xlims=(395,405), xticks=395:5:405, ylabel="[mW/m²/str/nm]")
q3 = plot(1e7./ν₁, (lieQ₀₁-lieQ0₀₁).*convfct1, linecolor=:black, xlabel = "λ [nm]")#
q2 = plot(1e7./ν₂, (lieQ₀₂-lieQ0₀₂).*convfct2, linecolor=:black, xlabel = "λ [nm]", xticks=424:2.5:429)
plot(p1, p2, p3, q1, q2, q3, layout = l, legend = false, title = ["ΔQₑ (incident)" "ΔQₑ (O₂)" "ΔQₑ (N₂)" " ΔQᵢ (RRS)" "ΔQᵢ (VRS+RVRS, O₂)" "ΔQᵢ (VRS+RVRS, N₂)"], titlefont = font(8))
savefig("DeltaRVRS_Qmonochromatic_sunlasercomp_vv.png")

norm=RS_type.F₀[1]*convfct0
l = @layout [a1 a2 a3; b1 b2 b3; c1 c2 c3]
o1 = plot(1e7./ν, (sum(τ_abs2[1], dims=2)).*convfct, linecolor=:red, xlims=(395,405), xticks=395:5:405, ylabel="τₐ")
o1 = plot!(1e7./ν, (sum(τ_abs02[1], dims=2)).*convfct, linecolor=:blue, xlims=(395,405), xticks=395:5:405)#, ylabel="\tau_{\mathrm{O}_3}")
o1 = plot!(1e7./ν₀, (sum(τ_abs[1], dims=2)).*convfct0, linecolor=:red, xlims=(395,405), xticks=395:5:405)
o1 = plot!(1e7./ν₀, (sum(τ_abs0[1], dims=2)).*convfct0, linecolor=:blue, xlims=(395,405), xticks=395:5:405)
o3 = plot(1e7./ν₁, (sum(τ_abs[2], dims=2)).*convfct1, linecolor=:red, xlims=(436,446), xticks=436:4:444) #, xticks=425:3:430)
o3 = plot!(1e7./ν₁, (sum(τ_abs0[2], dims=2)).*convfct1, linecolor=:blue) #, xticks=425:3:430)
o2 = plot(1e7./ν₂, (sum(τ_abs[3], dims=2)).*convfct2, linecolor=:red, xticks=424:2.5:429)
o2 = plot!(1e7./ν₂, (sum(τ_abs0[3], dims=2)).*convfct2, linecolor=:blue, xticks=424:2.5:429)

#p1 = plot(1e7./ν, (R_RRS2[1,1,:]-R0_RRS2[1,1,:]).*convfct, linecolor=:grey, xlims=(395,405), xticks=395:5:405, ylabel="[mW/m²/str/nm]")
#p1 = plot!(1e7./ν, (R_RRS[1,1,:]-R0_RRS[1,1,:]).*convfct, linecolor=:black, linewidth=2, xlims=(395,405), xticks=395:5:405, ylabel="[mW/m²/str/nm]")
#p1 = plot!(1e7./ν₀, (R₀₀-R0₀₀).*convfct0, linecolor=:black, xlims=(395,405))
#p3 = plot(1e7./ν₁, (R₀₁-R0₀₁).*convfct1, linecolor=:grey) #, xticks=425:3:430)
#p2 = plot(1e7./ν₂, (R₀₂-R0₀₂).*convfct2, linecolor=:grey, xticks=424:2.5:429)
p1 = plot(1e7./ν, (ieR0_RRS[1,1,:]).*convfct/norm, linecolor=:black, xlabel = "λ [nm]", xlims=(395,405), xticks=395:5:405)#, ylabel="[mW/m²/str/nm]")
p3 = plot(1e7./ν₁, (ieR0₀₁).*convfct1/norm, linecolor=:black, xlabel = "λ [nm]", xlims=(436,446), xticks=436:4:444)
p2 = plot(1e7./ν₂, (ieR0₀₂).*convfct2/norm, linecolor=:black, xlabel = "λ [nm]", xticks=424:2.5:429)


q1 = plot(1e7./ν, (ieR_RRS[1,1,:]-ieR0_RRS[1,1,:]).*convfct/norm, linecolor=:black, xlabel = "λ [nm]", xlims=(395,405), xticks=395:5:405)#, ylabel="[mW/m²/str/nm]")
q3 = plot(1e7./ν₁, (ieR₀₁-ieR0₀₁).*convfct1/norm, linecolor=:black, xlabel = "λ [nm]", xlims=(436,446), xticks=436:4:444)
q2 = plot(1e7./ν₂, (ieR₀₂-ieR0₀₂).*convfct2/norm, linecolor=:black, xlabel = "λ [nm]", xticks=424:2.5:429)
plot(o1, o2, o3, p1, p2, p3, q1, q2, q3, layout = l, legend = false, title = ["" "" "" "Iᵢ/F₀ (RRS)" "Iᵢ/F₀ (VRS+RVRS, O₂)" "Iᵢ/F₀ (VRS+RVRS, N₂)" " ΔIᵢ/F₀ (RRS)" "ΔIᵢ/F₀ (VRS+RVRS, O₂)" "ΔIᵢ/F₀ (VRS+RVRS, N₂)"], titlefont = font(7))
savefig("DeltaRVRS_monochromaticO3.png")

norm=RS_type.F₀[1]*convfct0
l = @layout [a1 a2 a3; b1 b2 b3; c1 c2 c3]
o1 = plot(1e7./ν, (sum(τ_abs2[1], dims=2)).*convfct, linecolor=:red, xlims=(395,405), xticks=395:5:405, ylabel="τₐ")
o1 = plot!(1e7./ν, (sum(τ_abs02[1], dims=2)).*convfct, linecolor=:blue, xlims=(395,405), xticks=395:5:405)#, ylabel="\tau_{\mathrm{O}_3}")
o1 = plot!(1e7./ν₀, (sum(τ_abs[1], dims=2)).*convfct0, linecolor=:red, xlims=(395,405), xticks=395:5:405)
o1 = plot!(1e7./ν₀, (sum(τ_abs0[1], dims=2)).*convfct0, linecolor=:blue, xlims=(395,405), xticks=395:5:405)
o3 = plot(1e7./ν₁, (sum(τ_abs[2], dims=2)).*convfct1, linecolor=:red, xlims=(436,446), xticks=436:4:444) #, xticks=425:3:430)
o3 = plot!(1e7./ν₁, (sum(τ_abs0[2], dims=2)).*convfct1, linecolor=:blue) #, xticks=425:3:430)
o2 = plot(1e7./ν₂, (sum(τ_abs[3], dims=2)).*convfct2, linecolor=:red, xticks=424:2.5:429)
o2 = plot!(1e7./ν₂, (sum(τ_abs0[3], dims=2)).*convfct2, linecolor=:blue, xticks=424:2.5:429)

#p1 = plot(1e7./ν, (R_RRS2[1,1,:]-R0_RRS2[1,1,:]).*convfct, linecolor=:grey, xlims=(395,405), xticks=395:5:405, ylabel="[mW/m²/str/nm]")
#p1 = plot!(1e7./ν, (R_RRS[1,1,:]-R0_RRS[1,1,:]).*convfct, linecolor=:black, linewidth=2, xlims=(395,405), xticks=395:5:405, ylabel="[mW/m²/str/nm]")
#p1 = plot!(1e7./ν₀, (R₀₀-R0₀₀).*convfct0, linecolor=:black, xlims=(395,405))
#p3 = plot(1e7./ν₁, (R₀₁-R0₀₁).*convfct1, linecolor=:grey) #, xticks=425:3:430)
#p2 = plot(1e7./ν₂, (R₀₂-R0₀₂).*convfct2, linecolor=:grey, xticks=424:2.5:429)
p1 = plot(1e7./ν, (ieR0_RRSl[1,1,:]).*convfct/norm, linecolor=:black, xlabel = "λ [nm]", xlims=(395,405), xticks=395:5:405)#, ylabel="[mW/m²/str/nm]")
p3 = plot(1e7./ν₁, (lieR0₀₁).*convfct1/norm, linecolor=:black, xlabel = "λ [nm]", xlims=(436,446), xticks=436:4:444)
p2 = plot(1e7./ν₂, (lieR0₀₂).*convfct2/norm, linecolor=:black, xlabel = "λ [nm]", xticks=424:2.5:429)


q1 = plot(1e7./ν, (ieR_RRSl[1,1,:]-ieR0_RRSl[1,1,:]).*convfct/norm, linecolor=:black, xlabel = "λ [nm]", xlims=(395,405), xticks=395:5:405)#, ylabel="[mW/m²/str/nm]")
q3 = plot(1e7./ν₁, (lieR₀₁-lieR0₀₁).*convfct1/norm, linecolor=:black, xlabel = "λ [nm]", xlims=(436,446), xticks=436:4:444)
q2 = plot(1e7./ν₂, (lieR₀₂-lieR0₀₂).*convfct2/norm, linecolor=:black, xlabel = "λ [nm]", xticks=424:2.5:429)
plot(o1, o2, o3, p1, p2, p3, q1, q2, q3, layout = l, legend = false, title = ["" "" "" "Iᵢ/F₀ (RRS)" "Iᵢ/F₀ (VRS+RVRS, O₂)" "Iᵢ/F₀ (VRS+RVRS, N₂)" " ΔIᵢ/F₀ (RRS)" "ΔIᵢ/F₀ (VRS+RVRS, O₂)" "ΔIᵢ/F₀ (VRS+RVRS, N₂)"], titlefont = font(7))
savefig("DeltaRVRS_monochromaticO3_sunlasercomp_vv.png")

norm=RS_type.F₀[1]*convfct0
l = @layout [a1 a2 a3; b1 b2 b3; c1 c2 c3]
o1 = plot(1e7./ν, (sum(τ_abs2[1], dims=2)).*convfct, linecolor=:red, xlims=(395,405), xticks=395:5:405, ylabel="τₐ")
o1 = plot!(1e7./ν, (sum(τ_abs02[1], dims=2)).*convfct, linecolor=:blue, xlims=(395,405), xticks=395:5:405)#, ylabel="\tau_{\mathrm{O}_3}")
o1 = plot!(1e7./ν₀, (sum(τ_abs[1], dims=2)).*convfct0, linecolor=:red, xlims=(395,405), xticks=395:5:405)
o1 = plot!(1e7./ν₀, (sum(τ_abs0[1], dims=2)).*convfct0, linecolor=:blue, xlims=(395,405), xticks=395:5:405)
o3 = plot(1e7./ν₁, (sum(τ_abs[2], dims=2)).*convfct1, linecolor=:red, xlims=(436,446), xticks=436:4:444)
o3 = plot!(1e7./ν₁, (sum(τ_abs0[2], dims=2)).*convfct1, linecolor=:blue) #, xticks=425:3:430)
o2 = plot(1e7./ν₂, (sum(τ_abs[3], dims=2)).*convfct2, linecolor=:red, xticks=424:2.5:429)
o2 = plot!(1e7./ν₂, (sum(τ_abs0[3], dims=2)).*convfct2, linecolor=:blue, xticks=424:2.5:429)

#p1 = plot(1e7./ν, (R_RRS2[1,1,:]-R0_RRS2[1,1,:]).*convfct, linecolor=:grey, xlims=(395,405), xticks=395:5:405, ylabel="[mW/m²/str/nm]")
#p1 = plot!(1e7./ν, (R_RRS[1,1,:]-R0_RRS[1,1,:]).*convfct, linecolor=:black, linewidth=2, xlims=(395,405), xticks=395:5:405, ylabel="[mW/m²/str/nm]")
#p1 = plot!(1e7./ν₀, (R₀₀-R0₀₀).*convfct0, linecolor=:black, xlims=(395,405))
#p3 = plot(1e7./ν₁, (R₀₁-R0₀₁).*convfct1, linecolor=:grey) #, xticks=425:3:430)
#p2 = plot(1e7./ν₂, (R₀₂-R0₀₂).*convfct2, linecolor=:grey, xticks=424:2.5:429)
p1 = plot(1e7./ν, (ieR0_RRS[1,2,:]).*convfct/norm, linecolor=:black, xlabel = "λ [nm]", xlims=(395,405), xticks=395:5:405) #, ylabel="[mW/m²/str/nm]")
p3 = plot(1e7./ν₁, (ieQ0₀₁).*convfct1/norm, linecolor=:black, xlabel = "λ [nm]", xlims=(436,446), xticks=436:4:444)
p2 = plot(1e7./ν₂, (ieQ0₀₂).*convfct2/norm, linecolor=:black, xlabel = "λ [nm]", xticks=424:2.5:429)


q1 = plot(1e7./ν, (ieR_RRS[1,2,:]-ieR0_RRS[1,2,:]).*convfct/norm, linecolor=:black, xlabel = "λ [nm]", xlims=(395,405), xticks=395:5:405) #, ylabel="[mW/m²/str/nm]")
q3 = plot(1e7./ν₁, (ieQ₀₁-ieQ0₀₁).*convfct1/norm, linecolor=:black, xlabel = "λ [nm]", xlims=(436,446), xticks=436:4:444)
q2 = plot(1e7./ν₂, (ieQ₀₂-ieQ0₀₂).*convfct2/norm, linecolor=:black, xlabel = "λ [nm]", xticks=424:2.5:429)
plot(o1, o2, o3, p1, p2, p3, q1, q2, q3, layout = l, legend = false, title = ["" "" "" "Qᵢ/F₀ (RRS)" "Qᵢ/F₀ (VRS+RVRS, O₂)" "Qᵢ/F₀ (VRS+RVRS, N₂)" " ΔQᵢ/F₀ (RRS)" "ΔQᵢ/F₀ (VRS+RVRS, O₂)" "ΔQᵢ/F₀ (VRS+RVRS, N₂)"], titlefont = font(7))
savefig("DeltaRVRS_QmonochromaticO3.png")

norm=RS_type.F₀[1]*convfct0
l = @layout [a1 a2 a3; b1 b2 b3; c1 c2 c3]
o1 = plot(1e7./ν, (sum(τ_abs2[1], dims=2)).*convfct, linecolor=:red, xlims=(395,405), xticks=395:5:405, ylabel="τₐ")
o1 = plot!(1e7./ν, (sum(τ_abs02[1], dims=2)).*convfct, linecolor=:blue, xlims=(395,405), xticks=395:5:405)#, ylabel="\tau_{\mathrm{O}_3}")
o1 = plot!(1e7./ν₀, (sum(τ_abs[1], dims=2)).*convfct0, linecolor=:red, xlims=(395,405), xticks=395:5:405)
o1 = plot!(1e7./ν₀, (sum(τ_abs0[1], dims=2)).*convfct0, linecolor=:blue, xlims=(395,405), xticks=395:5:405)
o3 = plot(1e7./ν₁, (sum(τ_abs[2], dims=2)).*convfct1, linecolor=:red, xlims=(436,446), xticks=436:4:444)
o3 = plot!(1e7./ν₁, (sum(τ_abs0[2], dims=2)).*convfct1, linecolor=:blue) #, xticks=425:3:430)
o2 = plot(1e7./ν₂, (sum(τ_abs[3], dims=2)).*convfct2, linecolor=:red, xticks=424:2.5:429)
o2 = plot!(1e7./ν₂, (sum(τ_abs0[3], dims=2)).*convfct2, linecolor=:blue, xticks=424:2.5:429)

#p1 = plot(1e7./ν, (R_RRS2[1,1,:]-R0_RRS2[1,1,:]).*convfct, linecolor=:grey, xlims=(395,405), xticks=395:5:405, ylabel="[mW/m²/str/nm]")
#p1 = plot!(1e7./ν, (R_RRS[1,1,:]-R0_RRS[1,1,:]).*convfct, linecolor=:black, linewidth=2, xlims=(395,405), xticks=395:5:405, ylabel="[mW/m²/str/nm]")
#p1 = plot!(1e7./ν₀, (R₀₀-R0₀₀).*convfct0, linecolor=:black, xlims=(395,405))
#p3 = plot(1e7./ν₁, (R₀₁-R0₀₁).*convfct1, linecolor=:grey) #, xticks=425:3:430)
#p2 = plot(1e7./ν₂, (R₀₂-R0₀₂).*convfct2, linecolor=:grey, xticks=424:2.5:429)
p1 = plot(1e7./ν, (ieR0_RRSl[1,2,:]).*convfct/norm, linecolor=:black, xlabel = "λ [nm]", xlims=(395,405), xticks=395:5:405) #, ylabel="[mW/m²/str/nm]")
p3 = plot(1e7./ν₁, (lieQ0₀₁).*convfct1/norm, linecolor=:black, xlabel = "λ [nm]", xlims=(436,446), xticks=436:4:444)
p2 = plot(1e7./ν₂, (lieQ0₀₂).*convfct2/norm, linecolor=:black, xlabel = "λ [nm]", xticks=424:2.5:429)


q1 = plot(1e7./ν, (ieR_RRSl[1,2,:]-ieR0_RRSl[1,2,:]).*convfct/norm, linecolor=:black, xlabel = "λ [nm]", xlims=(395,405), xticks=395:5:405) #, ylabel="[mW/m²/str/nm]")
q3 = plot(1e7./ν₁, (lieQ₀₁-lieQ0₀₁).*convfct1/norm, linecolor=:black, xlabel = "λ [nm]", xlims=(436,446), xticks=436:4:444)
q2 = plot(1e7./ν₂, (lieQ₀₂-lieQ0₀₂).*convfct2/norm, linecolor=:black, xlabel = "λ [nm]", xticks=424:2.5:429)
plot(o1, o2, o3, p1, p2, p3, q1, q2, q3, layout = l, legend = false, title = ["" "" "" "Qᵢ/F₀ (RRS)" "Qᵢ/F₀ (VRS+RVRS, O₂)" "Qᵢ/F₀ (VRS+RVRS, N₂)" " ΔQᵢ/F₀ (RRS)" "ΔQᵢ/F₀ (VRS+RVRS, O₂)" "ΔQᵢ/F₀ (VRS+RVRS, N₂)"], titlefont = font(7))
savefig("DeltaRVRS_QmonochromaticO3_sunlasercomp_vv.png")
#=
l = @layout [a1 a2 a3; b1 b2 b3]

p1 = plot(1e7./ν, R_RRS2[1,1,:].*convfct, linecolor=:grey, ylabel="[mW/m²/str/nm]")
p1 = plot!(1e7./ν, R_RRS[1,1,:].*convfct, linecolor=:black, linewidth=2, ylabel="[mW/m²/str/nm]")
p1 = plot!(1e7./ν₀, R₀₀.*convfct0, linecolor=:black)
p3 = plot(1e7./ν₁, R₀₁.*convfct1, linecolor=:grey) #, xticks=425:3:430)
p2 = plot(1e7./ν₂, R₀₂.*convfct2, linecolor=:grey)
q1 = plot(1e7./ν, ieR_RRS[1,1,:].*convfct, linecolor=:black, xlabel = "λ [nm]", ylabel="[mW/m²/str/nm]")
q3 = plot(1e7./ν₁, ieR₀₁.*convfct1, linecolor=:black, xlabel = "λ [nm]")#
q2 = plot(1e7./ν₂, ieR₀₂.*convfct2, linecolor=:black, xlabel = "λ [nm]")
plot(p1, p2, p3, q1, q2, q3, layout = l, legend = false, title = ["Iₑ (incident)" "Iₑ (O₂)" "Iₑ (N₂)" " Iᵢ (RRS)" "Iᵢ (VRS+RVRS, O₂)" "Iᵢ (VRS+RVRS, N₂)"], titlefont = font(8))
savefig("RVRS_monochromaticRed.png")

l = @layout [a1 a2 a3; b1 b2 b3]

p1 = plot(1e7./ν, R_RRS2[1,2,:].*convfct, linecolor=:grey, ylabel="[mW/m²/str/nm]")
p1 = plot!(1e7./ν, R_RRS[1,2,:].*convfct, linecolor=:black, linewidth=2, ylabel="[mW/m²/str/nm]")
p1 = plot!(1e7./ν₀, Q₀₀.*convfct0, linecolor=:black)
p3 = plot(1e7./ν₁, Q₀₁.*convfct1, linecolor=:grey) #, xticks=425:3:430)
p2 = plot(1e7./ν₂, Q₀₂.*convfct2, linecolor=:grey)
q1 = plot(1e7./ν, ieR_RRS[1,2,:].*convfct, linecolor=:black, xlabel = "λ [nm]", ylabel="[mW/m²/str/nm]")
q3 = plot(1e7./ν₁, ieQ₀₁.*convfct1, linecolor=:black, xlabel = "λ [nm]")#
q2 = plot(1e7./ν₂, ieQ₀₂.*convfct2, linecolor=:black, xlabel = "λ [nm]")
plot(p1, p2, p3, q1, q2, q3, layout = l, legend = false, title = ["Qₑ (incident)" "Qₑ (O₂)" "Qₑ (N₂)" " Qᵢ (RRS)" "Qᵢ (VRS+RVRS, O₂)" "Qᵢ (VRS+RVRS, N₂)"], titlefont = font(8))
savefig("RVRS_QmonochromaticRed.png")
=#