using CUDA
device!(1)
using DelimitedFiles
using Distributions
#using ImageFiltering
using InstrumentOperator
using Interpolations
using JLD2
using LaTeXStrings
using Plots
using Revise
using Statistics
using vSmartMOM
using vSmartMOM.CoreRT
using vSmartMOM.InelasticScattering
using vSmartMOM.SolarModel
# Benchmarks: http://web.gps.caltech.edu/~vijay/Rayleigh_Scattering_Tables/STOKES/
##

# Load YAML files into parameter struct
#parameters = 
#    parameters_from_yaml("test/test_parameters/FraunhoferMockParameters.yaml");
parameters = 
    parameters_from_yaml("test/test_parameters/aerosol_parameters.yaml");
model      = model_from_parameters(parameters);

n_bands = length(parameters.spec_bands)
n_aer = isnothing(parameters.scattering_params) ? 0 : length(parameters.scattering_params.rt_aerosols)
truncation_type = vSmartMOM.Scattering.δBGE{parameters.float_type}(parameters.l_trunc, parameters.Δ_angle)

#vmr = isnothing(parameters.absorption_params) ? Dict() : parameters.absorption_params.vmr
#p_full, p_half, vmr_h2o, vcd_dry, vcd_h2o, new_vmr = 
#    vSmartMOM.CoreRT.compute_atmos_profile_fields(parameters.T, parameters.p, parameters.q, vmr)

#profile = vSmartMOM.CoreRT.AtmosphericProfile(parameters.T, p_full, parameters.q, p_half, vmr_h2o, vcd_dry, vcd_h2o, new_vmr)
    
# Reduce the profile to the number of target layers (if specified)
#if parameters.profile_reduction_n != -1
#    profile = vSmartMOM.CoreRT.reduce_profile(parameters.profile_reduction_n, profile);
#end

# aerosol_optics[iBand][iAer]
aerosol_optics = [Array{vSmartMOM.Scattering.AerosolOptics}(undef, (n_aer)) for i=1:n_bands];
FT2 = isnothing(parameters.scattering_params) ? parameters.float_type : typeof(parameters.scattering_params.rt_aerosols[1].τ_ref)
#FT2 =  parameters.float_type 

# τ_aer[iBand][iAer,iZ]
τ_aer = [zeros(FT2, n_aer, length(model.profile.p_full)) for i=1:n_bands];

# Loop over aerosol type
for i_aer=1:n_aer

    # Get curr_aerosol
    c_aero = parameters.scattering_params.rt_aerosols[i_aer]
    curr_aerosol = c_aero.aerosol
    
    # Create Aerosol size distribution for each aerosol species
    size_distribution = curr_aerosol.size_distribution

    # Create a univariate aerosol distribution
    mie_aerosol = vSmartMOM.Scattering.Aerosol(size_distribution, curr_aerosol.nᵣ, curr_aerosol.nᵢ)
    #@show typeof(curr_aerosol.nᵣ)
    #mie_aerosol = make_mie_aerosol(size_distribution, curr_aerosol.nᵣ, curr_aerosol.nᵢ, parameters.scattering_params.r_max, parameters.scattering_params.nquad_radius) #Suniti: why is the refractive index needed here?

    # Create the aerosol extinction cross-section at the reference wavelength:
    mie_model      = vSmartMOM.Scattering.make_mie_model(parameters.scattering_params.decomp_type, 
                                    mie_aerosol, 
                                    parameters.scattering_params.λ_ref, 
                                    parameters.polarization_type, 
                                    truncation_type, 
                                    parameters.scattering_params.r_max, 
                                    parameters.scattering_params.nquad_radius)       
    k_ref          = vSmartMOM.Scattering.compute_ref_aerosol_extinction(mie_model, parameters.float_type)

    #parameters.scattering_params.rt_aerosols[i_aer].p₀, parameters.scattering_params.rt_aerosols[i_aer].σp
    # Loop over bands
    for i_band=1:n_bands
        
        # i'th spectral band (convert from cm⁻¹ to μm)
        curr_band_λ = 1e4 ./ parameters.spec_bands[i_band]

        # Create the aerosols:
        mie_model      = vSmartMOM.Scattering.make_mie_model(parameters.scattering_params.decomp_type, 
                                        mie_aerosol, 
                                        (maximum(curr_band_λ)+minimum(curr_band_λ))/2, 
                                        parameters.polarization_type, 
                                        truncation_type, 
                                        parameters.scattering_params.r_max, 
                                        parameters.scattering_params.nquad_radius)

        # Compute raw (not truncated) aerosol optical properties (not needed in RT eventually) 
        # @show FT2
        aerosol_optics_raw = vSmartMOM.Scattering.compute_aerosol_optical_properties(mie_model, FT2);

        # Compute truncated aerosol optical properties (phase function and fᵗ), consistent with Ltrunc:
        #@show i_aer, i_band
        aerosol_optics[i_band][i_aer] = vSmartMOM.Scattering.truncate_phase(truncation_type, 
                                                aerosol_optics_raw; reportFit=false)

        # Compute nAer aerosol optical thickness profiles
        τ_aer[i_band][i_aer,:] = 
            parameters.scattering_params.rt_aerosols[i_aer].τ_ref * 
            (aerosol_optics[i_band][i_aer].k/k_ref) * 
            vSmartMOM.CoreRT.getAerosolLayerOptProp(1, c_aero.profile.μ, c_aero.profile.σ, model.profile.p_half)
        @show vSmartMOM.CoreRT.getAerosolLayerOptProp(1, c_aero.profile.μ, c_aero.profile.σ, model.profile.p_half)
        @show aerosol_optics[i_band][i_aer].k, k_ref
        @show τ_aer[i_band][i_aer,:]
        @show parameters.scattering_params.rt_aerosols[i_aer].τ_ref
    end 
end


iBand = 1
FT = Float64

#### Compute all Raman properties
ν = model.params.spec_bands[iBand]
ν̃ = mean(ν);

    #parameters.depol = 0.041362343961163395 #0.028 #(0.028: Cabannes), (0.1032: Rayleigh) 
    #parameters = parameters_from_yaml("test/test_parameters/O2Parameters2.yaml");
# Create model struct (precomputes optical properties) from parameters
sza = [10. 70.]
I_conv = zeros(2,length(ν));
#I_conv_noRS = zeros(2,length(ν));
#ieI_conv = zeros(2,length(ν));

#I_conv_ss = zeros(2,length(ν));
#I_conv_noRS_ss = zeros(2,length(ν));
#ieI_conv_ss = zeros(2,length(ν));

Q_conv = zeros(2,length(ν));
#Q_conv_noRS = zeros(2,length(ν));
#ieQ_conv = zeros(2,length(ν));

#Q_conv_ss = zeros(2,length(ν));
#Q_conv_noRS_ss = zeros(2,length(ν));
#ieQ_conv_ss = zeros(2,length(ν));

U_conv = zeros(2,length(ν));
#U_conv_noRS = zeros(2,length(ν));
#ieU_conv = zeros(2,length(ν));

#U_conv_ss = zeros(2,length(ν));
#U_conv_noRS_ss = zeros(2,length(ν));
#ieU_conv_ss = zeros(2,length(ν));
#for ctr=1:2
    parameters.sza = sza[ctr]
    model      = model_from_parameters(parameters);

    iBand = 1
    FT = Float64

    #### Compute all Raman properties
    ν = model.params.spec_bands[iBand]
    ν̃ = mean(ν);

    # Find central reference index for RRS:
    i_ref = argmin(abs.(ν .- ν̃))
    T_sun = 5777. # K
    P = planck_spectrum_wn(T_sun, ν)
    Tsolar = solar_transmission_from_file("/home/sanghavi/code/github/vSmartMOM.jl/src/SolarModel/solar.out")
    Tsolar_interp = LinearInterpolation(Tsolar[4:end, 1], Tsolar[4:end, 2])
    F₀ = zeros(length(P));

    # TODO_VS: λ_vs_in (get input)
    # TODO_VS: ν_vs_in (convert to wavenumbers)
    # Effective temperature for Raman calculations
    effT = 300.  #(model.profile.vcd_dry' * model.profile.T) / sum(model.profile.vcd_dry);
    # Define RS type
    # Compute N2 and O2

    #n2,o2 = InelasticScattering.getRamanAtmoConstants(ν̃,effT);
    ##greek_raman = get_greek_raman(RS_type, n2, o2);
    #=RS_type = InelasticScattering.RRS(
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
    =#

    # Compute Raman SSA properties:
    #CoreRT.getRamanSSProp!(RS_type, parameters.depol, 1e7/ν̃, ν);

    RS_type.F₀ = zeros(model.params.polarization_type.n, length(P))
    for i=1:length(P)
        sol_trans = Tsolar_interp(ν[i]);
        F₀[i] = sol_trans * P[i];
        RS_type.F₀[1,i] = F₀[i];
    end 
    #RS_type.F₀ = zeros(model.params.polarization_type.n, length(ν))
    #RS_type.F₀[1,:].=1.

    #R, T, ieR, ieT = CoreRT.rt_run_test(RS_type,model,iBand);
    #R_ss, T_ss, ieR_ss, ieT_ss = CoreRT.rt_run_test_ss(RS_type,model,iBand);

    RS_type = InelasticScattering.noRS(
        fscattRayl  = [FT(1)],
        ϖ_Cabannes  = [FT(1)], 
        bandSpecLim = [],
        iBand       = [1],
        F₀          = zeros(FT,1,1));
    RS_type.F₀ = zeros(model.params.polarization_type.n, length(P))
    for i=1:length(P)
        sol_trans = Tsolar_interp(ν[i]);
        F₀[i] = sol_trans * P[i];
        RS_type.F₀[1,i] = F₀[i];
    end 
        #RS_type.F₀ = zeros(model.params.polarization_type.n, length(ν))
    #RS_type.F₀[1,:].=1.

    R, T, _, _ = CoreRT.rt_run_test(RS_type,model,iBand);
    #RnoRS_ss, TnoRS_ss, _, _ = CoreRT.rt_run_test_ss(RS_type,model,iBand);

    #test = jldopen("/home/sanghavi/debugRay2_.jld2")
    #@load "/home/sanghavi/debugRRS2.jld2"
    #RayJ₀ = test["J₀"];
    #test = jldopen("/home/sanghavi/debugRRS2.jld2")
    #ieJ₀ = test["ieJ₀"];
    #test = jldopen("/home/sanghavi/debugCab2.jld2")
    #CabJ₀ = test["J₀"];
    #RayJ₀[:,:,640] == sum(ieJ₀[:,:,640,:], dims=3) + CabJ₀[:,:,640]
    #RayJ₀[:,:,640] ≈ sum(ieJ₀[:,:,640,:], dims=3) + CabJ₀[:,:,640]

    #@load "/home/sanghavi/debugCab4.jld2"
    #@load "/home/sanghavi/debugRRS4.jld2"
    #@load "/home/sanghavi/debugRay3.jld2"

    #=
    a = sum(ieR[:,:,640,:], dims=3) + CabR[:,:,640]
    RayR[:,:,640] .≈ a[:,:,1]
    b = sum(ieT[:,:,640,:], dims=3) + CabT[:,:,640]
    RayR[:,:,640] .≈ a[:,:,1]
    c = sum(ieJ₀m[:,:,640,:], dims=3) + CabJ₀m[:,:,640]
    RayJ₀m[:,:,640] .≈ c[:,:,1]
    d = sum(ieJ₀p[:,:,640,:], dims=3) + CabJ₀p[:,:,640]
    RayJ₀p[:,:,640] .≈ d[:,:,1]
    =#

    x = -40:0.3:40
    #kernel = InstrumentOperator.create_instrument_kernel(Normal(0, 12.5), x) #defining a Gaussian kernel for convulution in wavenumber space
    #kernel = InstrumentOperator.create_instrument_kernel(Normal(0, 7.64), x)
    #kernel = InstrumentOperator.create_instrument_kernel(Normal(0, 3.82), x)
    kernel = InstrumentOperator.create_instrument_kernel(Normal(0, 3.), x)
    
    #I_conv = InstrumentOperator.conv_spectra(kernel, )
    #I_conv_noRS[ctr,:] .= imfilter(RnoRS[1,1,:], kernel)
    I_conv[ctr,:] .= imfilter(R[1,1,:], kernel)
    #ieI_conv[ctr,:] .= imfilter(ieR[1,1,:], kernel)

    #Q_conv_noRS[ctr,:] .= imfilter(RnoRS[1,2,:], kernel)
    Q_conv[ctr,:] .= imfilter(R[1,2,:], kernel)
    #ieQ_conv[ctr,:] .= imfilter(ieR[1,2,:], kernel)

    #U_conv_noRS[ctr,:] .= imfilter(RnoRS[1,3,:], kernel)
    U_conv[ctr,:] .= imfilter(R[1,3,:], kernel)
    #ieU_conv[ctr,:] .= imfilter(ieR[1,3,:], kernel)

    #I_conv_noRS_ss[ctr,:] .= imfilter(RnoRS_ss[1,1,:], kernel)
    #I_conv_ss[ctr,:] .= imfilter(R_ss[1,1,:], kernel)
    #ieI_conv_ss[ctr,:] .= imfilter(ieR_ss[1,1,:], kernel)

    #Q_conv_noRS_ss[ctr,:] .= imfilter(RnoRS_ss[1,2,:], kernel)
    #Q_conv_ss[ctr,:] .= imfilter(R_ss[1,2,:], kernel)
    #ieQ_conv_ss[ctr,:] .= imfilter(ieR_ss[1,2,:], kernel)

    #U_conv_noRS_ss[ctr,:] .= imfilter(RnoRS_ss[1,3,:], kernel)
    #U_conv_ss[ctr,:] .= imfilter(R_ss[1,3,:], kernel)
    #ieU_conv_ss[ctr,:] .= imfilter(ieR_ss[1,3,:], kernel)
#end

convfct = 1e7./ν.^2   # to convert radiance units from mW/m²-str-cm⁻¹ to mW/m²-str-nm
l = @layout [a1 a2 ; b1 b2; c1 c2]

p1 = plot(1e7./ν, 100*(I_conv[1,:]+ieI_conv[1,:]-I_conv_noRS[1,:])./I_conv_noRS[1,:], xlims=(385,405), ylabel=L"R_I [\%]", linecolor=:black, linewidth = 2, alpha=0.7)
p1 = plot!(1e7./ν, 100*(I_conv_ss[1,:]+ieI_conv_ss[1,:]-I_conv_noRS_ss[1,:])./I_conv_noRS_ss[1,:], linecolor=:red, linewidth = 2, alpha=0.7)

q1 = plot(1e7./ν, 100*(I_conv[2,:]+ieI_conv[2,:]-I_conv_noRS[2,:])./I_conv_noRS[2,:], xlims=(385,405), linecolor=:black, linewidth = 2, alpha=0.7)
q1 = plot!(1e7./ν, 100*(I_conv_ss[2,:]+ieI_conv_ss[2,:]-I_conv_noRS_ss[2,:])./I_conv_noRS_ss[2,:], linecolor=:red, linewidth = 2, alpha=0.7)

p2 = plot(1e7./ν, 100*(Q_conv[1,:]+ieQ_conv[1,:]-Q_conv_noRS[1,:])./Q_conv_noRS[1,:], xlims=(385,405), ylabel=L"R_Q [\%]", linecolor=:black, linewidth = 2, alpha=0.7)
p2 = plot!(1e7./ν, 100*(Q_conv_ss[1,:]+ieQ_conv_ss[1,:]-Q_conv_noRS_ss[1,:])./Q_conv_noRS_ss[1,:], linecolor=:red, linewidth = 2, alpha=0.7)

q2 = plot(1e7./ν, 100*(Q_conv[2,:]+ieQ_conv[2,:]-Q_conv_noRS[2,:])./Q_conv_noRS[2,:], xlims=(385,405), linecolor=:black, linewidth = 2, alpha=0.7)
q2 = plot!(1e7./ν, 100*(Q_conv_ss[2,:]+ieQ_conv_ss[2,:]-Q_conv_noRS_ss[2,:])./Q_conv_noRS_ss[2,:], linecolor=:red, linewidth = 2, alpha=0.7)

p3 = plot(1e7./ν, 100*(U_conv[1,:]+ieU_conv[1,:]-U_conv_noRS[1,:])./U_conv_noRS[1,:], xlims=(385,405), ylabel=L"R_U [\%]", linecolor=:black, linewidth = 2, alpha=0.7)
p3 = plot!(1e7./ν, 100*(U_conv_ss[1,:]+ieU_conv_ss[1,:]-U_conv_noRS_ss[1,:])./U_conv_noRS_ss[1,:], xlabel="Wavelength [nm]", linecolor=:red, linewidth = 2, alpha=0.7)

q3 = plot(1e7./ν, 100*(U_conv[2,:]+ieU_conv[2,:]-U_conv_noRS[2,:])./U_conv_noRS[2,:], xlims=(385,405), linecolor=:black, linewidth = 2, alpha=0.7)
q3 = plot!(1e7./ν, 100*(U_conv_ss[2,:]+ieU_conv_ss[2,:]-U_conv_noRS_ss[2,:])./U_conv_noRS_ss[2,:], xlabel="Wavelength [nm]", linecolor=:red, linewidth = 2, alpha=0.7)

plot(p1, q1, p2, q2, p3, q3, layout = l, link=:y, legend = false, title = ["θ₀=10ᵒ" "θ₀=70ᵒ" "" "" "" ""], titlefont = font(10))
savefig("R_comparison_Lelli1.png")


#ϖ_Cabannes, γ_air_Cabannes, γ_air_Rayleigh = vSmartMOM.InelasticScattering.compute_γ_air_Rayleigh!(400.)
ϖ_Cabannes, γ_air_Cabannes, γ_air_Rayleigh = vSmartMOM.InelasticScattering.compute_γ_air_Rayleigh!(395.)
#ϖ_n2_Cabannes, γ_n2_Cabannes, γ_n2_Rayleigh = vSmartMOM.InelasticScattering.compute_γ_mol_Rayleigh!(395.,n2)
γ_air_RRS = 3/4
θ = 0:180

Z11_Cab = (3/(8π*(1+2γ_air_Cabannes)))*(0.5*(1-γ_air_Cabannes)*((cosd.(θ)).^2 .+ 1) .+ 2γ_air_Cabannes)
Z12_Cab = (3/(8π*(1+2γ_air_Cabannes)))*(0.5*(1-γ_air_Cabannes)*((cosd.(θ)).^2 .- 1))

Z11_RRS = (3/(8π*(1+2γ_air_RRS)))*(0.5*(1-γ_air_RRS)*((cosd.(θ)).^2 .+ 1) .+ 2γ_air_RRS)
Z12_RRS = (3/(8π*(1+2γ_air_RRS)))*(0.5*(1-γ_air_RRS)*((cosd.(θ)).^2 .-1 ))

Z11_Ray = (3/(8π*(1+2γ_air_Rayleigh)))*(0.5*(1-γ_air_Rayleigh)*((cosd.(θ)).^2 .+ 1) .+ 2γ_air_Rayleigh)
Z12_Ray = (3/(8π*(1+2γ_air_Rayleigh)))*(0.5*(1-γ_air_Rayleigh)*((cosd.(θ)).^2 .- 1))

plot(θ, Z11_RRS./Z11_Cab, label="j=1", ylabel=L"\mathbf{Z}_{1j, \mathrm{RRS}}/\mathbf{Z}_{1j, \mathrm{e}}")
plot!(θ, Z12_RRS./Z12_Cab, label="j=2", xlabel="Θ[ᵒ]", ylims=(0,1.4))
savefig("Z_comparison_Landgraf.png")

l = @layout [a1; b1]
a=plot(θ, Z11_Ray*4π, label="Rayleigh", color=:blue, ylabel=L"\mathbf{Z}_{11}")#, legend=:bottom)
a=plot!(θ, Z11_Cab*4π, label="Cabannes", color=:green, xlabel="Θ[ᵒ]")#, ylims=(0.7,1.5))
a=plot!(θ, Z11_RRS*4π, label="RRS", color=:red)#, xlabel="Θ[ᵒ]", ylims=(0.7,1.5))

b=plot(θ, -Z12_Ray./Z11_Ray, color=:blue, label="", ylabel=L"-\mathbf{Z}_{12}/\mathbf{Z}_{11}")
b=plot!(θ, -Z12_RRS./Z11_RRS, color=:red, label="")#, xlabel="Θ[ᵒ]", ylims=(0.7,1.5))
b=plot!(θ, -Z12_Cab./Z11_Cab, color=:green, label="", xlabel="Θ[ᵒ]")#, ylims=(0.7,1.5))

plot(a, b, layout = l, link=:x, legend = :top, title = ["normalized (to 4π) phase function" "degree of linear polarization"], titlefont = font(10))

savefig("Z_comparison_Stam.png")