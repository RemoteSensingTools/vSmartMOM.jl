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
using DelimitedFiles
#
## run this code using the following command:
## /net/fluo/data2/software/Julia/julia-1.11.3/bin/julia --project=/home/sanghavi/code/github/vSmartMOM.jl/ /home/sanghavi/code/github/vSmartMOM.jl/test/benchmarks/test_creategrid_noabs_RamanSIF.jl &

FT = Float64

# Load YAML files into parameter struct
# noabs
# parameters = 
#    parameters_from_yaml("/home/sanghavi/code/github/vSmartMOM.jl/test/test_parameters/noabs_parameters2_SIF_grid.yaml");
# w/abs
parameters = 
    parameters_from_yaml("/home/sanghavi/code/github/vSmartMOM.jl/test/test_parameters/abs_parameters2_SIF_grid.yaml");

#parameters = parameters_from_yaml("test/test_parameters/O2Parameters2.yaml");
# Create model struct (precomputes optical properties) from parameters
model      = model_from_parameters(parameters);
#=model.τ_aer[1][1,:,:] .= 0;
model.τ_aer[1][1,:,end] .= 100.
model.τ_rayl[1][:,end-1] .+= model.τ_rayl[1][:,end];
model.τ_rayl[1][:,end] .= 0;=#
a=[]
for i=0:20
    push!(a, acosd(i/20))  
end
sza=reverse(Int.(ceil.(a[8:21])))
ρ = zeros(FT,21) #ρ = zeros(FT,15)
ρ_str = []
for iρ = 1:21 #15
    ρ[iρ] = (iρ-1)*0.05
    for i=1:length(parameters.spec_bands)
        parameters.brdf[i].albedo = ρ[iρ]
    end
    push!(ρ_str, replace(string(round(ρ[iρ], digits=2)),"."=>"p"))
end
psurf=[1000 750 500]
J_SIF = readdlm("/home/sanghavi/code/github/vSmartMOM.jl/src/SIF_emission/sif-spectra.csv", ',') 
#J_SIF = readdlm("/home/sanghavi/code/github/vSmartMOM.jl/src/SIF_emission/PC1_SIFSpectra_allSpecies.csv", ',') 
νSIF = reverse(1e7./J_SIF[2:end,1])
jSIF = reverse(J_SIF[2:end,2]*(0.5*π/maximum(J_SIF[2:end,2]))) # mW/m²/nm
jSIF .*= 1e7./(νSIF.^2) # mW/m²/cm⁻¹
SIF_interp = LinearInterpolation(νSIF, jSIF)
Tsolar = solar_transmission_from_file("/home/sanghavi/code/github/vSmartMOM.jl/src/SolarModel/solar.out")
Tsolar_interp = LinearInterpolation(Tsolar[4:end, 1], Tsolar[4:end, 2])
n_bands = length(parameters.spec_bands);
T_sun = 5777. # K

for isurf = 1:1 # 3:3 # 2:2 # 1:1 # 
    for iρ in [1 7 21 17] # = 21:21 #7:7 #21 #1:15 #3 #1:15
        for iA = 2:2 #14
            parameters.sza = 1.0*sza[iA]
            model.obs_geom = CoreRT.ObsGeometry(parameters.sza, parameters.vza, parameters.vaz, parameters.obs_alt)
            # Set quadrature points for streams
            model.quad_points = CoreRT.rt_set_streams(parameters.quadrature_type, 
                                                parameters.l_trunc, 
                                                model.obs_geom, 
                                                parameters.polarization_type, 
                                                array_type(parameters.architecture))

            #model.obs_geom.sza = 0.0;#1.0*sza[iA]
            #@show parameters.sza, sza[iA]
            #@show parameters.vza

            
            for iBand=1:n_bands
                parameters.brdf[iBand].albedo = ρ[iρ]
                #model.params.brdf[iBand].albedo = ρ[iρ]
                #@show model.params.brdf[iBand].albedo
            end

            #model      = model_from_parameters(parameters);
            overlap_ν = 250 #230
            Δν = parameters.spec_bands[1][2] - parameters.spec_bands[1][1]; #Make sure all bands have the same spectral resolution Δν
            n_overlap = Int(floor(overlap_ν/Δν))
            tot_nspec = sum([length(parameters.spec_bands[i]) for i=1:length(parameters.spec_bands)]) - 2*(n_bands)*n_overlap
            n_cam = length(parameters.vza) #for generating the grid, the SZA and the VZA will be swapped, so that the resulting Rs and ieRs will need to be multiplied by cos(VZA) for the true value 
            RnoRS = zeros(FT, n_cam, parameters.polarization_type.n, tot_nspec)
            TnoRS = zeros(FT, n_cam, parameters.polarization_type.n, tot_nspec)
            R = zeros(FT, n_cam, parameters.polarization_type.n, tot_nspec)
            T = zeros(FT, n_cam, parameters.polarization_type.n, tot_nspec)
            ieR = zeros(FT, n_cam, parameters.polarization_type.n, tot_nspec)
            ieT = zeros(FT, n_cam, parameters.polarization_type.n, tot_nspec)

            iBand = 1
            tot_ν=zeros(tot_nspec);
            tot_F₀ = zeros(tot_nspec);
            spec_end = 0
            #Fdir = zeros(tot_nspec);

            # Common operations for all bands
            #### Compute all Raman properties
            ν̃ = 0.5*(model.params.spec_bands[1][1]+model.params.spec_bands[end][end])
            
            # Find central reference index for RRS:
            #i_ref = argmin(abs.(ν .- ν̃))
            # TODO_VS: λ_vs_in (get input)
            # TODO_VS: ν_vs_in (convert to wavenumbers)
            # Effective temperature for Raman calculations
            effT = (model.profile.vcd_dry' * model.profile.T) / sum(model.profile.vcd_dry);
            # Define RS type
            # Compute N2 and O2

            n2,o2 = InelasticScattering.getRamanAtmoConstants(ν̃,effT);
            #greek_raman = get_greek_raman(RS_type, n2, o2);
            RS_type0 = InelasticScattering.noRS(
                        fscattRayl  = [FT(1)],
                        ϖ_Cabannes  = [FT(1)], 
                        bandSpecLim = [],
                        iBand       = [1],
                        F₀          = zeros(FT,1,1),
                        SIF₀        = zeros(FT,1,1));            

            RS_type1 = InelasticScattering.RRS(
                        n2=n2,
                        o2=o2,
                        greek_raman = InelasticScattering.GreekCoefs([FT(1)], [FT(1)], [FT(1)], [FT(1)], [FT(1)], [FT(1)]),
                        fscattRayl  = [FT(1)],
                        ϖ_Cabannes  = [FT(1)], 
                        ϖ_λ₁λ₀      = zeros(FT,1),
                        i_λ₁λ₀      = zeros(Int,1), 
                        Z⁻⁺_λ₁λ₀    = zeros(FT,1,1), 
                        Z⁺⁺_λ₁λ₀    = zeros(FT,1,1), 
                        i_ref       = 0, # this is a redundant parameter
                        n_Raman     = 0,
                        F₀          = zeros(FT,1,1),
                        SIF₀        = zeros(FT,1,1));

            for iBand=1:n_bands
                #### Compute all Raman properties
                ν = model.params.spec_bands[iBand]
                parameters.brdf[iBand].albedo = ρ[iρ]
                spec_start = spec_end+1
                spec_end += length(ν) - 2*n_overlap
                
                # Compute Raman SSA properties:
                CoreRT.getRamanSSProp!(RS_type1, 1e7/mean(ν), ν);

                P = planck_spectrum_wn(T_sun, ν) * 2.1629e-05 * π  # mW/m²-cm⁻¹
                
                F₀ = zeros(length(P));
                SIF₀ = zeros(length(ν))
                RS_type0.F₀ = zeros(model.params.polarization_type.n, length(P))
                RS_type1.F₀ = zeros(model.params.polarization_type.n, length(P))
                RS_type0.SIF₀ = zeros(model.params.polarization_type.n, length(ν))
                RS_type1.SIF₀ = zeros(model.params.polarization_type.n, length(ν))
                ######################
                #TEMPORARY for w/FR!!#
                ######################
                for i=1:length(P)
                    sol_trans = Tsolar_interp(ν[i]);
                    F₀[i] = sol_trans * P[i];
                    SIF₀[i] = SIF_interp(ν[i]); #0.0
                    RS_type0.F₀[1,i] = F₀[i]; #1.0 #
                    RS_type1.F₀[1,i] = F₀[i];
                    RS_type0.SIF₀[1,i] = SIF₀[i]; #1.0 #
                    RS_type1.SIF₀[1,i] = SIF₀[i];
                end 
                ######################
                #TEMPORARY for w/FR!!#
                ######################
                ######################
                #TEMPORARY for NOFR!!#
                ######################
                #=F₀ .= 1.0 # = Tsolar_interp.(ν) .* P;
                SIF₀ .= 0.0 # = SIF_interp.(ν); # 
                RS_type0.F₀[1,:] = F₀; #1.0 #
                RS_type1.F₀[1,:] = F₀;
                RS_type0.SIF₀[1,:] = SIF₀; #1.0 #
                RS_type1.SIF₀[1,:] = SIF₀;
                =#
                ######################
                #TEMPORARY for NOFR!!#
                ######################

                #@show iBand, spec_start, spec_end, spec_end - spec_start + 1
                #@show iBand, n_overlap + 1, length(ν) - n_overlap,  length(ν) - 2*n_overlap

                R1, T1, ieR1, ieT1 = CoreRT.rt_run_test(RS_type1, model,iBand);
                RnoRS0, TnoRS0, _, _ = CoreRT.rt_run_test(RS_type0, model,iBand);
                tot_ν[spec_start:spec_end] = ν[(n_overlap+1):(end-n_overlap)] 
                tot_F₀[spec_start:spec_end] = F₀[(n_overlap+1):(end-n_overlap)] 
                R[:,:,spec_start:spec_end] = R1[:,:,(n_overlap+1):(end-n_overlap)]
                ieR[:,:,spec_start:spec_end] = ieR1[:,:,(n_overlap+1):(end-n_overlap)]
                RnoRS[:,:,spec_start:spec_end] = RnoRS0[:,:,(n_overlap+1):(end-n_overlap)]
                #@show R[1,:,1], ieR[1,:,1],RnoRS[1,:,1] 
            end

            vza_str = replace(string(round(model.params.vza[1], digits=1)),"."=>"p")
            albedo = ρ_str[iρ]
            sza_str = string(sza[iA])
            for vctr = 1:length(model.params.vaz)
                vaz_str = replace(string(round(model.params.vaz[vctr], digits=1)),"."=>"p")
                
                # With SIF
                #fname0 = "/home/sanghavi/data/RamanSIFgrid/tests/test_raylSIF_sza"*sza_str*"_vza"*vza_str*"_vaz"*vaz_str*"_alb"*albedo*"_psurf"*string(psurf[isurf])*"hpa_nors_ABO2.dat"
                #fname1 = "/home/sanghavi/data/RamanSIFgrid/tests/test_raylSIF_sza"*sza_str*"_vza"*vza_str*"_vaz"*vaz_str*"_alb"*albedo*"_psurf"*string(psurf[isurf])*"hpa_rrs_ABO2.dat"
                # No SIF
                #fname0 = "/home/sanghavi/data/RamanSIFgrid/tests/albtest_nofr_noabs_rayl_sza"*sza_str*"_vza"*vza_str*"_vaz"*vaz_str*"_alb"*albedo*"_psurf"*string(psurf[isurf])*"hpa_nors_ABO2.dat"
                #fname1 = "/home/sanghavi/data/RamanSIFgrid/tests/albtest_nofr_noabs_rayl_sza"*sza_str*"_vza"*vza_str*"_vaz"*vaz_str*"_alb"*albedo*"_psurf"*string(psurf[isurf])*"hpa_rrs_ABO2.dat"
    
                #nofr, noabs
                #fname0 = "/home/sanghavi/data/RamanSIFgrid/tests/test_nofr_noabs_rayl_sza"*sza_str*"_vza"*vza_str*"_vaz"*vaz_str*"_alb"*albedo*"_psurf"*string(psurf[isurf])*"hpa_nors_ABO2.dat"
                #fname1 = "/home/sanghavi/data/RamanSIFgrid/tests/test_nofr_noabs_rayl_sza"*sza_str*"_vza"*vza_str*"_vaz"*vaz_str*"_alb"*albedo*"_psurf"*string(psurf[isurf])*"hpa_rrs_ABO2.dat"
                # w/fr, noabs
                #fname0 = "/home/sanghavi/data/RamanSIFgrid/tests/test_wfr_noabs_rayl_sza"*sza_str*"_vza"*vza_str*"_vaz"*vaz_str*"_alb"*albedo*"_psurf"*string(psurf[isurf])*"hpa_nors_ABO2.dat"
                #fname1 = "/home/sanghavi/data/RamanSIFgrid/tests/test_wfr_noabs_rayl_sza"*sza_str*"_vza"*vza_str*"_vaz"*vaz_str*"_alb"*albedo*"_psurf"*string(psurf[isurf])*"hpa_rrs_ABO2.dat"
                # w/fr, w/abs
                fname0 = "/home/sanghavi/data/RamanSIFgrid/tests/test_wfr_wabs_rayl_sza"*sza_str*"_vza"*vza_str*"_vaz"*vaz_str*"_alb"*albedo*"_psurf"*string(psurf[isurf])*"hpa_nors_ABO2.dat"
                fname1 = "/home/sanghavi/data/RamanSIFgrid/tests/test_wfr_wabs_rayl_sza"*sza_str*"_vza"*vza_str*"_vaz"*vaz_str*"_alb"*albedo*"_psurf"*string(psurf[isurf])*"hpa_rrs_ABO2.dat"
                
                rayl_rrs_ABO2 = []
                rayl_nors_ABO2 = []
                
                rayl_rrs_ABO2 = [tot_ν R[vctr,1,:] R[vctr,2,:] R[vctr,3,:] ieR[vctr,1,:] ieR[vctr,2,:] ieR[vctr,3,:] tot_F₀]
                rayl_nors_ABO2 = [tot_ν RnoRS[vctr,1,:] RnoRS[vctr,2,:] RnoRS[vctr,3,:] tot_F₀]

                writedlm(fname0, rayl_nors_ABO2)
                writedlm(fname1, rayl_rrs_ABO2)
                #@show rayl_nors_ABO2
            end    
            # No SIF
            #fname0 = "/home/sanghavi/data/RamanSIFgrid/rayl_sza"*sza_str*"_alb"*albedo*"_psurf"*string(psurf[isurf])*"hpa_nors_ABO2.dat"
            #fname1 = "/home/sanghavi/data/RamanSIFgrid/rayl_sza"*sza_str*"_alb"*albedo*"_psurf"*string(psurf[isurf])*"hpa_rrs_ABO2.dat"

            # With SIF
            #fname0 = "/home/sanghavi/data/RamanSIFgrid/raylSIF_sza"*sza_str*"_alb"*albedo*"_psurf"*string(psurf[isurf])*"hpa_nors_ABO2.dat"
            #fname1 = "/home/sanghavi/data/RamanSIFgrid/raylSIF_sza"*sza_str*"_alb"*albedo*"_psurf"*string(psurf[isurf])*"hpa_rrs_ABO2.dat"
            
            # Testing
            #fname0 = "/home/sanghavi/data/RamanSIFgrid/test1raylSIF_sza"*sza_str*"_alb"*albedo*"_psurf"*string(psurf[isurf])*"hpa_nors_ABO2.dat"
            #fname1 = "/home/sanghavi/data/RamanSIFgrid/test1raylSIF_sza"*sza_str*"_alb"*albedo*"_psurf"*string(psurf[isurf])*"hpa_rrs_ABO2.dat"
            
            #writedlm(fname0, rayl_nors_ABO2)
            #writedlm(fname1, rayl_rrs_ABO2)
        end
    end
end
#=
rayl_sza70_vza00_nors_ABO2_Δp = [ν RnoRS[1,1,:] RnoRS[1,2,:] RnoRS[1,3,:]]
rayl_sza70_vza00_nors_ABO2_SIF = [ν RnoRS[1,1,:] RnoRS[1,2,:] RnoRS[1,3,:]]

rayl_sza30_vza00_rrs_ABO2 = [ν R[1,1,:] R[1,2,:] R[1,3,:] ieR[1,1,:] ieR[1,2,:] ieR[1,3,:]]
rayl_sza30_vza00_nors_ABO2 = [ν RnoRS[1,1,:] RnoRS[1,2,:] RnoRS[1,3,:]]
rayl_sza30_vza00_nors_ABO2_Δp = [ν RnoRS[1,1,:] RnoRS[1,2,:] RnoRS[1,3,:]]
rayl_sza30_vza00_nors_ABO2_SIF = [ν RnoRS[1,1,:] RnoRS[1,2,:] RnoRS[1,3,:]]


writedlm("/home/sanghavi/CMS_RamanSIF/rayl_sza70_vza00_rrs_ABO2.dat", rayl_sza70_vza00_rrs_ABO2)
writedlm("/home/sanghavi/CMS_RamanSIF/rayl_sza70_vza00_nors_ABO2.dat", rayl_sza70_vza00_nors_ABO2)
writedlm("/home/sanghavi/CMS_RamanSIF/rayl_sza70_vza00_nors_ABO2_dp.dat", rayl_sza70_vza00_nors_ABO2_Δp)
writedlm("/home/sanghavi/CMS_RamanSIF/rayl_sza70_vza00_nors_ABO2_SIF.dat", rayl_sza70_vza00_nors_ABO2_SIF)

writedlm("/home/sanghavi/CMS_RamanSIF/rayl_sza30_vza00_rrs_ABO2.dat", rayl_sza30_vza00_rrs_ABO2)
writedlm("/home/sanghavi/CMS_RamanSIF/rayl_sza30_vza00_nors_ABO2.dat", rayl_sza30_vza00_nors_ABO2)
writedlm("/home/sanghavi/CMS_RamanSIF/rayl_sza30_vza00_nors_ABO2_dp.dat", rayl_sza30_vza00_nors_ABO2_Δp)
writedlm("/home/sanghavi/CMS_RamanSIF/rayl_sza30_vza00_nors_ABO2_SIF.dat", rayl_sza30_vza00_nors_ABO2_SIF)

rayl_sza70_vza00_rrs_ABO2 = readdlm("/home/sanghavi/CMS_RamanSIF/rayl_sza70_vza00_rrs_ABO2.dat")
rayl_sza70_vza00_nors_ABO2 = readdlm("/home/sanghavi/CMS_RamanSIF/rayl_sza70_vza00_nors_ABO2.dat")
rayl_sza70_vza00_nors_ABO2_Δp = readdlm("/home/sanghavi/CMS_RamanSIF/rayl_sza70_vza00_nors_ABO2_dp.dat")
rayl_sza70_vza00_nors_ABO2_SIF = readdlm("/home/sanghavi/CMS_RamanSIF/rayl_sza70_vza00_nors_ABO2_SIF.dat")

rayl_sza30_vza00_rrs_ABO2 = readdlm("/home/sanghavi/CMS_RamanSIF/rayl_sza30_vza00_rrs_ABO2.dat")
rayl_sza30_vza00_nors_ABO2 = readdlm("/home/sanghavi/CMS_RamanSIF/rayl_sza30_vza00_nors_ABO2.dat")
rayl_sza30_vza00_nors_ABO2_Δp = readdlm("/home/sanghavi/CMS_RamanSIF/rayl_sza30_vza00_nors_ABO2_dp.dat")
rayl_sza30_vza00_nors_ABO2_SIF = readdlm("/home/sanghavi/CMS_RamanSIF/rayl_sza30_vza00_nors_ABO2_SIF.dat")

=#


#RnoRS_ss, TnoRS_ss, _, _ = CoreRT.rt_run_test_ss(RS_type,model,iBand);

#=
R_test, T_test, ieR_test, ieT_test = CoreRT.rt_run_test(RS_type,model,iBand);

# You can now run multiple bands like this (list at the end of band!)
RnoRS_test, TnoRS_test, _, _ = CoreRT.rt_run_test(vSmartMOM.noRS(),model,[1,1]);


R_test, T_test, ieR_test, ieT_test = CoreRT.rt_run_test(RS_type,model,1);
=#
#=
rayl_sza70_vza00_rrs_ABO2 = readdlm("/home/sanghavi/CMS_RamanSIF/rayl_sza70_vza00_rrs_ABO2.dat")
rayl_sza70_vza00_nors_ABO2 = readdlm("/home/sanghavi/CMS_RamanSIF/rayl_sza70_vza00_nors_ABO2.dat")
rayl_sza70_vza00_nors_ABO2_Δp = readdlm("/home/sanghavi/CMS_RamanSIF/rayl_sza70_vza00_nors_ABO2_dp.dat")
rayl_sza70_vza00_nors_ABO2_SIF = readdlm("/home/sanghavi/CMS_RamanSIF/rayl_sza70_vza00_nors_ABO2_SIF.dat")

rayl_sza30_vza00_rrs_ABO2 = readdlm("/home/sanghavi/CMS_RamanSIF/rayl_sza30_vza00_rrs_ABO2.dat")
rayl_sza30_vza00_nors_ABO2 = readdlm("/home/sanghavi/CMS_RamanSIF/rayl_sza30_vza00_nors_ABO2.dat")
rayl_sza30_vza00_nors_ABO2_Δp = readdlm("/home/sanghavi/CMS_RamanSIF/rayl_sza30_vza00_nors_ABO2_dp.dat")
rayl_sza30_vza00_nors_ABO2_SIF = readdlm("/home/sanghavi/CMS_RamanSIF/rayl_sza30_vza00_nors_ABO2_SIF.dat")

ν = rayl_sza70_vza00_nors_ABO2[:,1]
I70_rrs = zeros(length(ν),2)
ieI70_rrs = zeros(length(ν),2)
I70_nors = zeros(length(ν),2)
I70_nors_Δp = zeros(length(ν),2)
I70_nors_SIF = zeros(length(ν),2)

I30_rrs = zeros(length(ν),2)
ieI30_rrs = zeros(length(ν),2)
I30_nors = zeros(length(ν),2)
I30_nors_Δp = zeros(length(ν),2)
I30_nors_SIF = zeros(length(ν),2)
=#
#===Convolution of hires spectral simulations to instrument grid===#
#x = (1e7/415):0.3:(1e7/385)
#ν = (1e7/460):0.3:(1e7/420)
#x = -40:0.15:40

#=
x = -10:0.18:10
FWHM_OCO_λ = 0.042 #nm
FWHM_OCO_ν = FWHM_OCO_λ*(mean(ν))^2 * 1e-7 # cm⁻¹
σOCO_λ = FWHM_OCO_λ/2.355
σOCO_ν = FWHM_OCO_ν/2.355
#kernel = InstrumentOperator.create_instrument_kernel(Normal(0, 12.5), x) #defining a Gaussian kernel for convolution in wavenumber space
kernel_ν = InstrumentOperator.create_instrument_kernel(Normal(0, σOCO_ν), x)
#FWHM_OCO = 0.042 #nm
#σOCO = FWHM_OCO/2.355
#kernel = InstrumentOperator.create_instrument_kernel(Normal(0, 12.5), x) #defining a Gaussian kernel for convolution in wavenumber space
#kernel = InstrumentOperator.create_instrument_kernel(Normal(0, σOCO), x)
#kernel = InstrumentOperator.create_instrument_kernel(Normal(0, 3.75), x) #SCIAMACHY resolution of σ=0.54 nm in the O2 A-band 
convfct = ν.^2 / 1e7   # to convert radiance units from mW/m²-str-cm⁻¹ to mW/m²-str-nm

for ctr=1:2
    I70_nors[:,ctr]     .= imfilter(rayl_sza70_vza00_nors_ABO2[:,ctr+1], kernel_ν)
    I70_nors_Δp[:,ctr]  .= imfilter(rayl_sza70_vza00_nors_ABO2_Δp[:,ctr+1], kernel_ν)
    I70_nors_SIF[:,ctr] .= imfilter(rayl_sza70_vza00_nors_ABO2_SIF[:,ctr+1], kernel_ν)
    I70_rrs[:,ctr]      .= imfilter(rayl_sza70_vza00_rrs_ABO2[:,ctr+1], kernel_ν)
    ieI70_rrs[:,ctr]      .= imfilter(rayl_sza70_vza00_rrs_ABO2[:,ctr+4], kernel_ν)

    I30_nors[:,ctr]     .= imfilter(rayl_sza30_vza00_nors_ABO2[:,ctr+1], kernel_ν)
    I30_nors_Δp[:,ctr]  .= imfilter(rayl_sza30_vza00_nors_ABO2_Δp[:,ctr+1], kernel_ν)
    I30_nors_SIF[:,ctr] .= imfilter(rayl_sza30_vza00_nors_ABO2_SIF[:,ctr+1], kernel_ν)
    I30_rrs[:,ctr]      .= imfilter(rayl_sza30_vza00_rrs_ABO2[:,ctr+1], kernel_ν)
    ieI30_rrs[:,ctr]      .= imfilter(rayl_sza30_vza00_rrs_ABO2[:,ctr+4], kernel_ν)
end

#convfct = 1.0#./F₀ # to normalize wrt irradiance 
l = @layout [a1 a2] # ; b1 b2; c1 c2]
p1 = plot(1e7./ν, ((I70_rrs[:,1]+ieI70_rrs[:,1])./I70_nors[:,1] .- 1)*100, linecolor=:red, linewidth=1.5, alpha=0.75, xlims=(755,775), label="RRS")
p1 = plot!(1e7./ν, (I70_nors_Δp[:,1]./I70_nors[:,1] .- 1)*100, linecolor=:green, linewidth=1.5, alpha=0.75, xlims=(755,775), label="Δp")
p1 = plot!(1e7./ν, (I70_nors_SIF[:,1]./I70_nors[:,1] .- 1)*100, linecolor=:blue, linewidth=1.5, alpha=0.75, xlims=(755,775), label = "SIF", ylabel="ΔI [%]", yguidefontsize=8)

q1 = plot(1e7./ν, ((I70_rrs[:,2]+ieI70_rrs[:,2])./I70_nors[:,2] .- 1)*100, linecolor=:red, linewidth=1.5, alpha=0.75, xlims=(755,775), label="RRS")
q1 = plot!(1e7./ν, (I70_nors_Δp[:,2]./I70_nors[:,2] .- 1)*100, linecolor=:green, linewidth=1.5, alpha=0.75, xlims=(755,775), label="Δp")
q1 = plot!(1e7./ν, (I70_nors_SIF[:,2]./I70_nors[:,2] .- 1)*100, linecolor=:blue, linewidth=1.5, alpha=0.75, xlims=(755,775), label="SIF", ylabel="ΔQ [%]", yguidefontsize=8)

plot(p1, q1, layout = l, plot_title = "SZA=70ᵒ, ρ=0.3, SIF=0.5 mW/m²/str/nm", titlefont = font(10))
savefig("/home/sanghavi/CMS_RamanSIF/Relative_RRS_dp_SIF_O2Aband_SZA70.png")

l = @layout [a1 a2] # ; b1 b2; c1 c2]
p1 = plot(1e7./ν, ((I30_rrs[:,1]+ieI30_rrs[:,1])./I30_nors[:,1] .- 1)*100, linecolor=:red, linewidth=1.5, alpha=0.75, xlims=(755,775), label="RRS")
p1 = plot!(1e7./ν, (I30_nors_Δp[:,1]./I30_nors[:,1] .- 1)*100, linecolor=:green, linewidth=1.5, alpha=0.75, xlims=(755,775), label="Δp")
p1 = plot!(1e7./ν, (I30_nors_SIF[:,1]./I30_nors[:,1] .- 1)*100, linecolor=:blue, linewidth=1.5, alpha=0.75, xlims=(755,775), label = "SIF", ylabel="ΔI [%]", yguidefontsize=8)

q1 = plot(1e7./ν, ((I30_rrs[:,2]+ieI30_rrs[:,2])./I30_nors[:,2] .- 1)*100, linecolor=:red, linewidth=1.5, alpha=0.75, xlims=(755,775), label="RRS")
q1 = plot!(1e7./ν, (I30_nors_Δp[:,2]./I30_nors[:,2] .- 1)*100, linecolor=:green, linewidth=1.5, alpha=0.75, xlims=(755,775), label="Δp")
q1 = plot!(1e7./ν, (I30_nors_SIF[:,2]./I30_nors[:,2] .- 1)*100, linecolor=:blue, linewidth=1.5, alpha=0.75, xlims=(755,775), label="SIF", ylabel="ΔQ [%]", yguidefontsize=8)

plot(p1, q1, layout = l, plot_title = "SZA=30ᵒ, ρ=0.3, SIF=0.5 mW/m²/str/nm", titlefont = font(10))
savefig("/home/sanghavi/CMS_RamanSIF/Relative_RRS_dp_SIF_O2Aband_SZA30.png")

# ZOOM
#convfct = 1.0#./F₀ # to normalize wrt irradiance 
l = @layout [a1 a2] # ; b1 b2; c1 c2]
p1 = plot(1e7./ν, ((I70_rrs[:,1]+ieI70_rrs[:,1])./I70_nors[:,1] .- 1)*100, linecolor=:red, linewidth=1.5, alpha=0.75, xlims=(755.5,759.5), label="RRS")
p1 = plot!(1e7./ν, (I70_nors_Δp[:,1]./I70_nors[:,1] .- 1)*100, linecolor=:green, linewidth=1.5, alpha=0.75, xlims=(755.5,759.5), label="Δp")
p1 = plot!(1e7./ν, (I70_nors_SIF[:,1]./I70_nors[:,1] .- 1)*100, linecolor=:blue, linewidth=1.5, alpha=0.75, xlims=(755.5,759.5), label = "SIF", ylabel="ΔI [%]", yguidefontsize=8)

q1 = plot(1e7./ν, ((I70_rrs[:,2]+ieI70_rrs[:,2])./I70_nors[:,2] .- 1)*100, linecolor=:red, linewidth=1.5, alpha=0.75, xlims=(755,759.5), label="RRS")
q1 = plot!(1e7./ν, (I70_nors_Δp[:,2]./I70_nors[:,2] .- 1)*100, linecolor=:green, linewidth=1.5, alpha=0.75, xlims=(755.5,759.5), label="Δp")
q1 = plot!(1e7./ν, (I70_nors_SIF[:,2]./I70_nors[:,2] .- 1)*100, linecolor=:blue, linewidth=1.5, alpha=0.75, xlims=(755.5,759.5), label="SIF", ylabel="ΔQ [%]", yguidefontsize=8)

plot(p1, q1, layout = l, plot_title = "SZA=70ᵒ, ρ=0.3, SIF=0.5 mW/m²/str/nm", titlefont = font(10))
savefig("/home/sanghavi/CMS_RamanSIF/Relative_RRS_dp_SIF_O2Aband_SZA70_zoom.png")

l = @layout [a1 a2] # ; b1 b2; c1 c2]
p1 = plot(1e7./ν, ((I30_rrs[:,1]+ieI30_rrs[:,1])./I30_nors[:,1] .- 1)*100, linecolor=:red, linewidth=1.5, alpha=0.75, xlims=(755.5,759.5), label="RRS")
p1 = plot!(1e7./ν, (I30_nors_Δp[:,1]./I30_nors[:,1] .- 1)*100, linecolor=:green, linewidth=1.5, alpha=0.75, xlims=(755.5,759.5), label="Δp")
p1 = plot!(1e7./ν, (I30_nors_SIF[:,1]./I30_nors[:,1] .- 1)*100, linecolor=:blue, linewidth=1.5, alpha=0.75, xlims=(755.5,759.5), label = "SIF", ylabel="ΔI [%]", yguidefontsize=8)

q1 = plot(1e7./ν, ((I30_rrs[:,2]+ieI30_rrs[:,2])./I30_nors[:,2] .- 1)*100, linecolor=:red, linewidth=1.5, alpha=0.75, xlims=(755.5,759.5), label="RRS")
q1 = plot!(1e7./ν, (I30_nors_Δp[:,2]./I30_nors[:,2] .- 1)*100, linecolor=:green, linewidth=1.5, alpha=0.75, xlims=(755.5,759.5), label="Δp")
q1 = plot!(1e7./ν, (I30_nors_SIF[:,2]./I30_nors[:,2] .- 1)*100, linecolor=:blue, linewidth=1.5, alpha=0.75, xlims=(755.5,759.5), label="SIF", ylabel="ΔQ [%]", yguidefontsize=8)

plot(p1, q1, layout = l, plot_title = "SZA=30ᵒ, ρ=0.3, SIF=0.5 mW/m²/str/nm", titlefont = font(10))
savefig("/home/sanghavi/CMS_RamanSIF/Relative_RRS_dp_SIF_O2Aband_SZA30_zoom.png")

#Only I
#convfct = 1.0#./F₀ # to normalize wrt irradiance 
l = @layout [a1] # ; b1 b2; c1 c2]
p1 = plot(1e7./ν, ((I70_rrs[:,1]+ieI70_rrs[:,1])./I70_nors[:,1] .- 1)*100, linecolor=:red, linewidth=1.5, alpha=0.75, xlims=(755,775), label="RRS")
p1 = plot!(1e7./ν, (I70_nors_Δp[:,1]./I70_nors[:,1] .- 1)*100, linecolor=:green, linewidth=1.5, alpha=0.75, xlims=(755,775), label="Δp")
p1 = plot!(1e7./ν, (I70_nors_SIF[:,1]./I70_nors[:,1] .- 1)*100, linecolor=:blue, linewidth=1.5, alpha=0.75, xlims=(755,775), label = "SIF", ylabel="ΔI [%]", yguidefontsize=8)

#q1 = plot(1e7./ν, ((I70_rrs[:,2]+ieI70_rrs[:,2])./I70_nors[:,2] .- 1)*100, linecolor=:red, linewidth=1.5, alpha=0.75, xlims=(755,775), label="RRS")
#q1 = plot!(1e7./ν, (I70_nors_Δp[:,2]./I70_nors[:,2] .- 1)*100, linecolor=:green, linewidth=1.5, alpha=0.75, xlims=(755,775), label="Δp")
#q1 = plot!(1e7./ν, (I70_nors_SIF[:,2]./I70_nors[:,2] .- 1)*100, linecolor=:blue, linewidth=1.5, alpha=0.75, xlims=(755,775), label="SIF", ylabel="ΔQ [%]", yguidefontsize=8)

plot(p1, layout = l, plot_title = "SZA=70ᵒ, ρ=0.3, SIF=0.5 mW/m²/str/nm", titlefont = font(10))
savefig("/home/sanghavi/CMS_RamanSIF/Relative_RRS_dp_SIF_O2Aband_SZA70_I.png")

l = @layout [a1] # ; b1 b2; c1 c2]
p1 = plot(1e7./ν, ((I30_rrs[:,1]+ieI30_rrs[:,1])./I30_nors[:,1] .- 1)*100, linecolor=:red, linewidth=1.5, alpha=0.75, xlims=(755,775), label="RRS")
p1 = plot!(1e7./ν, (I30_nors_Δp[:,1]./I30_nors[:,1] .- 1)*100, linecolor=:green, linewidth=1.5, alpha=0.75, xlims=(755,775), label="Δp")
p1 = plot!(1e7./ν, (I30_nors_SIF[:,1]./I30_nors[:,1] .- 1)*100, linecolor=:blue, linewidth=1.5, alpha=0.75, xlims=(755,775), label = "SIF", ylabel="ΔI [%]", yguidefontsize=8)

#q1 = plot(1e7./ν, ((I30_rrs[:,2]+ieI30_rrs[:,2])./I30_nors[:,2] .- 1)*100, linecolor=:red, linewidth=1.5, alpha=0.75, xlims=(755,775), label="RRS")
#q1 = plot!(1e7./ν, (I30_nors_Δp[:,2]./I30_nors[:,2] .- 1)*100, linecolor=:green, linewidth=1.5, alpha=0.75, xlims=(755,775), label="Δp")
#q1 = plot!(1e7./ν, (I30_nors_SIF[:,2]./I30_nors[:,2] .- 1)*100, linecolor=:blue, linewidth=1.5, alpha=0.75, xlims=(755,775), label="SIF", ylabel="ΔQ [%]", yguidefontsize=8)

plot(p1, layout = l, plot_title = "SZA=30ᵒ, ρ=0.3, SIF=0.5 mW/m²/str/nm", titlefont = font(10))
savefig("/home/sanghavi/CMS_RamanSIF/Relative_RRS_dp_SIF_O2Aband_SZA30_I.png")

# ZOOM
#convfct = 1.0#./F₀ # to normalize wrt irradiance 
l = @layout [a1] # ; b1 b2; c1 c2]
p1 = plot(1e7./ν, ((I70_rrs[:,1]+ieI70_rrs[:,1])./I70_nors[:,1] .- 1)*100, linecolor=:red, linewidth=1.5, alpha=0.75, xlims=(755.5,759.5), label="RRS")
p1 = plot!(1e7./ν, (I70_nors_Δp[:,1]./I70_nors[:,1] .- 1)*100, linecolor=:green, linewidth=1.5, alpha=0.75, xlims=(755.5,759.5), label="Δp")
p1 = plot!(1e7./ν, (I70_nors_SIF[:,1]./I70_nors[:,1] .- 1)*100, linecolor=:blue, linewidth=1.5, alpha=0.75, xlims=(755.5,759.5), label = "SIF", ylabel="ΔI [%]", yguidefontsize=8)

#q1 = plot(1e7./ν, ((I70_rrs[:,2]+ieI70_rrs[:,2])./I70_nors[:,2] .- 1)*100, linecolor=:red, linewidth=1.5, alpha=0.75, xlims=(755,759.5), label="RRS")
#q1 = plot!(1e7./ν, (I70_nors_Δp[:,2]./I70_nors[:,2] .- 1)*100, linecolor=:green, linewidth=1.5, alpha=0.75, xlims=(755.5,759.5), label="Δp")
#q1 = plot!(1e7./ν, (I70_nors_SIF[:,2]./I70_nors[:,2] .- 1)*100, linecolor=:blue, linewidth=1.5, alpha=0.75, xlims=(755.5,759.5), label="SIF", ylabel="ΔQ [%]", yguidefontsize=8)

plot(p1, layout = l, plot_title = "SZA=70ᵒ, ρ=0.3, SIF=0.5 mW/m²/str/nm", titlefont = font(10))
savefig("/home/sanghavi/CMS_RamanSIF/Relative_RRS_dp_SIF_O2Aband_SZA70_zoom_I.png")

l = @layout [a1] # ; b1 b2; c1 c2]
p1 = plot(1e7./ν, ((I30_rrs[:,1]+ieI30_rrs[:,1])./I30_nors[:,1] .- 1)*100, linecolor=:red, linewidth=1.5, alpha=0.75, xlims=(755.5,759.5), label="RRS")
p1 = plot!(1e7./ν, (I30_nors_Δp[:,1]./I30_nors[:,1] .- 1)*100, linecolor=:green, linewidth=1.5, alpha=0.75, xlims=(755.5,759.5), label="Δp")
p1 = plot!(1e7./ν, (I30_nors_SIF[:,1]./I30_nors[:,1] .- 1)*100, linecolor=:blue, linewidth=1.5, alpha=0.75, xlims=(755.5,759.5), label = "SIF", ylabel="ΔI [%]", yguidefontsize=8)

#q1 = plot(1e7./ν, ((I30_rrs[:,2]+ieI30_rrs[:,2])./I30_nors[:,2] .- 1)*100, linecolor=:red, linewidth=1.5, alpha=0.75, xlims=(755.5,759.5), label="RRS")
#q1 = plot!(1e7./ν, (I30_nors_Δp[:,2]./I30_nors[:,2] .- 1)*100, linecolor=:green, linewidth=1.5, alpha=0.75, xlims=(755.5,759.5), label="Δp")
#q1 = plot!(1e7./ν, (I30_nors_SIF[:,2]./I30_nors[:,2] .- 1)*100, linecolor=:blue, linewidth=1.5, alpha=0.75, xlims=(755.5,759.5), label="SIF", ylabel="ΔQ [%]", yguidefontsize=8)

plot(p1, layout = l, plot_title = "SZA=30ᵒ, ρ=0.3, SIF=0.5 mW/m²/str/nm", titlefont = font(10))
savefig("/home/sanghavi/CMS_RamanSIF/Relative_RRS_dp_SIF_O2Aband_SZA30_zoom_I.png")


# Combined
l = @layout [a1 a2; b1 b2] # ; b1 b2; c1 c2]
p1 = plot(1e7./ν, ((I70_rrs[:,1]+ieI70_rrs[:,1])./I70_nors[:,1] .- 1)*100, linecolor=:red, linewidth=1.5, alpha=0.75, xlims=(755,775), label="RRS")
p1 = plot!(1e7./ν, (I70_nors_Δp[:,1]./I70_nors[:,1] .- 1)*100, linecolor=:green, linewidth=1.5, alpha=0.75, xlims=(755,775), label="Δp")
p1 = plot!(1e7./ν, (I70_nors_SIF[:,1]./I70_nors[:,1] .- 1)*100 .- 0.25, linecolor=:blue, linewidth=1.5, alpha=0.75, xlims=(755,775), label = "SIF", ylabel="ΔI [%]", yguidefontsize=8)

#q1 = plot(1e7./ν, ((I70_rrs[:,2]+ieI70_rrs[:,2])./I70_nors[:,2] .- 1)*100, linecolor=:red, linewidth=1.5, alpha=0.75, xlims=(685,695), label="RRS")
#q1 = plot!(1e7./ν, (I70_nors_Δp[:,2]./I70_nors[:,2] .- 1)*100, linecolor=:green, linewidth=1.5, alpha=0.75, xlims=(685,695), label="Δp")
#q1 = plot!(1e7./ν, (I70_nors_SIF[:,2]./I70_nors[:,2] .- 1)*100, linecolor=:blue, linewidth=1.5, alpha=0.75, xlims=(685,695), label="SIF", ylabel="ΔQ [%]", yguidefontsize=8)

#plot(p1, layout = l, plot_title = "SZA=70ᵒ, ρ=0.05, SIF=0.5 mW/m²/str/nm", titlefont = font(8))
#savefig("/home/sanghavi/CMS_RamanSIF/Relative_RRS_dp_SIF_O2Bband_SZA70_I.png")
#savefig("/home/sanghavi/CMS_RamanSIF/SCIA_Relative_RRS_dp_SIF_O2Bband_SZA70.png")

#l = @layout [a1] # ; b1 b2; c1 c2]
p2 = plot(1e7./ν, ((I30_rrs[:,1]+ieI30_rrs[:,1])./I30_nors[:,1] .- 1)*100, linecolor=:red, linewidth=1.5, alpha=0.75, xlims=(755,775), label="RRS")
p2 = plot!(1e7./ν, (I30_nors_Δp[:,1]./I30_nors[:,1] .- 1)*100, linecolor=:green, linewidth=1.5, alpha=0.75, xlims=(755,775), label="Δp")
p2 = plot!(1e7./ν, (I30_nors_SIF[:,1]./I30_nors[:,1] .- 1)*100 .- 0.1, linecolor=:blue, linewidth=1.5, alpha=0.75, xlims=(755,775), label = "SIF", ylabel="ΔI [%]", yguidefontsize=8)

#q1 = plot(1e7./ν, ((I30_rrs[:,2]+ieI30_rrs[:,2])./I30_nors[:,2] .- 1)*100, linecolor=:red, linewidth=1.5, alpha=0.75, xlims=(685,695), label="RRS")
#q1 = plot!(1e7./ν, (I30_nors_Δp[:,2]./I30_nors[:,2] .- 1)*100, linecolor=:green, linewidth=1.5, alpha=0.75, xlims=(685,695), label="Δp")
#q1 = plot!(1e7./ν, (I30_nors_SIF[:,2]./I30_nors[:,2] .- 1)*100, linecolor=:blue, linewidth=1.5, alpha=0.75, xlims=(685,695), label="SIF", ylabel="ΔQ [%]", yguidefontsize=8)

#plot(p1, layout = l, plot_title = "SZA=30ᵒ, ρ=0.05, SIF=0.5 mW/m²/str/nm", titlefont = font(8))
#savefig("/home/sanghavi/CMS_RamanSIF/Relative_RRS_dp_SIF_O2Bband_SZA30_I.png")
#savefig("/home/sanghavi/CMS_RamanSIF/SCIA_Relative_RRS_dp_SIF_O2Bband_SZA30.png")

# ZOOM
#l = @layout [a1] # ; b1 b2; c1 c2]
p1z = plot(1e7./ν, ((I70_rrs[:,1]+ieI70_rrs[:,1])./I70_nors[:,1] .- 1)*100, linecolor=:red, linewidth=1.5, alpha=0.75, xlims=(755.5,759.5), ylims=(-0.1,0.1), label="RRS")
p1z = plot!(1e7./ν, (I70_nors_Δp[:,1]./I70_nors[:,1] .- 1)*100, linecolor=:green, linewidth=1.5, alpha=0.75, xlims=(755.5,759.5), ylims=(-0.1,0.1), label="Δp")
p1z = plot!(1e7./ν, (I70_nors_SIF[:,1]./I70_nors[:,1] .- 1)*100 .- 0.25, linecolor=:blue, linewidth=1.5, alpha=0.75, xlims=(755.5,759.5), ylims=(-0.1,0.1), label = "SIF", ylabel="ΔI [%]", yguidefontsize=8)

#q1 = plot(1e7./ν, ((I70_rrs[:,2]+ieI70_rrs[:,2])./I70_nors[:,2] .- 1)*100, linecolor=:red, linewidth=1.5, alpha=0.75, xlims=(685.8,686.2), label="RRS")
#q1 = plot!(1e7./ν, (I70_nors_Δp[:,2]./I70_nors[:,2] .- 1)*100, linecolor=:green, linewidth=1.5, alpha=0.75, xlims=(685.8,686.2), label="Δp")
#q1 = plot!(1e7./ν, (I70_nors_SIF[:,2]./I70_nors[:,2] .- 1)*100, linecolor=:blue, linewidth=1.5, alpha=0.75, xlims=(685.8,686.2), label="SIF", ylabel="ΔQ [%]", yguidefontsize=8)

#plot(p1, layout = l, plot_title = "SZA=70ᵒ, ρ=0.05, SIF=0.5 mW/m²/str/nm", titlefont = font(8))
#savefig("/home/sanghavi/CMS_RamanSIF/Relative_RRS_dp_SIF_O2Bband_SZA70_zoom_I.png")
#savefig("/home/sanghavi/CMS_RamanSIF/SCIA_Relative_RRS_dp_SIF_O2Bband_SZA70.png")

#l = @layout [a1] # ; b1 b2; c1 c2]
p2z = plot(1e7./ν, ((I30_rrs[:,1]+ieI30_rrs[:,1])./I30_nors[:,1] .- 1)*100, linecolor=:red, linewidth=1.5, alpha=0.75, xlims=(755.5,759.5), ylims=(-0.1,0.1), label="")
p2z = plot!(1e7./ν, (I30_nors_Δp[:,1]./I30_nors[:,1] .- 1)*100, linecolor=:green, linewidth=1.5, alpha=0.75, xlims=(755.5,759.5), ylims=(-0.1,0.1), label="")
p2z = plot!(1e7./ν, (I30_nors_SIF[:,1]./I30_nors[:,1] .- 1)*100 .- 0.1, linecolor=:blue, linewidth=1.5, alpha=0.75, xlims=(755.5,759.5), ylims=(-0.1,0.1), label = "", ylabel="ΔI [%]", yguidefontsize=8)

#q1 = plot(1e7./ν, ((I30_rrs[:,2]+ieI30_rrs[:,2])./I30_nors[:,2] .- 1)*100, linecolor=:red, linewidth=1.5, alpha=0.75, xlims=(685.8,686.2), label="RRS")
#q1 = plot!(1e7./ν, (I30_nors_Δp[:,2]./I30_nors[:,2] .- 1)*100, linecolor=:green, linewidth=1.5, alpha=0.75, xlims=(685.8,686.2), label="Δp")
#q1 = plot!(1e7./ν, (I30_nors_SIF[:,2]./I30_nors[:,2] .- 1)*100, linecolor=:blue, linewidth=1.5, alpha=0.75, xlims=(685.8,686.2), label="SIF", ylabel="ΔQ [%]", yguidefontsize=8)

plot(p1, p1z, p2, p2z, layout = l, plot_title = "ρ=0.3, Δp=3 hPa, SIF=0.5 mW/m²/str/nm", ploy_titlefont = font(8), title = ["SZA=70ᵒ" "SZA=70ᵒ" "SZA=30ᵒ" "SZA=30ᵒ"], titlefont = font(8))
savefig("/home/sanghavi/CMS_RamanSIF/Relative_RRS_dp_SIF_O2Aband_withzoom_I.pdf")
=#

#=
I_conv_noRS = zeros(3,length(ν))
I_conv = zeros(3,length(ν))
ieI_conv = zeros(3,length(ν))
Q_conv_noRS = zeros(3,length(ν))
Q_conv = zeros(3,length(ν))
ieQ_conv = zeros(3,length(ν))

I_conv_noRS_bin = zeros(3,length(ν))
I_conv_bin = zeros(3,length(ν))
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

convfct = 1e7./ν.^2   # to convert radiance units from mW/m²-str-cm⁻¹ to mW/m²-str-nm
#convfct = 1.0#./F₀ # to normalize wrt irradiance 
l = @layout [a1 a2 ; b1 b2; c1 c2]
p1 = plot(1e7./ν, (RnoRS[1,1,:]).*convfct, linecolor=:red, linewidth=0.5, alpha=0.5, xlims=(755,775))
p1 = plot!(1e7./ν, (RnoRS[2,1,:]).*convfct, linecolor=:blue, linewidth=0.5, alpha=0.5, xlims=(755,775))
p1 = plot!(1e7./ν, (RnoRS[3,1,:]).*convfct, linecolor=:green, linewidth=0.5, alpha=0.5, xlims=(755,775), ylabel="[mW/m²/str/nm]", yguidefontsize=8)
#p1 = plot!(1e7./ν, (I_conv_noRS[1,:]).*convfct, linewidth = 2, linecolor=:red, xlims=(755,775))
#p1 = plot!(1e7./ν, (I_conv_noRS[2,:]).*convfct, linewidth = 2, linecolor=:blue, xlims=(755,775))
#p1 = plot!(1e7./ν, (I_conv_noRS[3,:]).*convfct, linewidth = 2, linecolor=:green, xlims=(755,775), ylabel="[mW/m²/str/nm]", yguidefontsize=8)


p2 = plot(1e7./ν, (R[1,1,:].-RnoRS[1,1,:].+ieR[1,1,:]).*convfct, linecolor=:red, linewidth=0.5, alpha=0.5, xlims=(755,775))#, ylims=(-2e-6, 1.e-5))#, xlabel = "λ [nm]")
p2 = plot!(1e7./ν, (R[2,1,:].-RnoRS[2,1,:].+ieR[2,1,:]).*convfct, linecolor=:blue, linewidth=0.5, alpha=0.5, xlims=(755,775))#, xlabel = "λ [nm]")
p2 = plot!(1e7./ν, (R[3,1,:].-RnoRS[3,1,:].+ieR[3,1,:]).*convfct, linecolor=:green, linewidth=0.5, alpha=0.5, xlims=(755,775), ylims=(-1,2.5), ylabel="[mW/m²/str/nm]", yguidefontsize=8)#, xlabel = "λ [nm]", xlims=(755,775))
#p2 = plot!(1e7./ν, (I_conv[1,:].-I_conv_noRS[1,:].+ieI_conv[1,:]).*convfct, linewidth = 2, linecolor=:red, xlims=(755,775))
#p2 = plot!(1e7./ν, (I_conv[2,:].-I_conv_noRS[2,:].+ieI_conv[2,:]).*convfct, linewidth = 2, linecolor=:blue, xlims=(755,775))
#p2 = plot!(1e7./ν, (I_conv[3,:].-I_conv_noRS[3,:].+ieI_conv[3,:]).*convfct, linewidth = 2, linecolor=:green, xlims=(755,775), ylims=(-1,2.5), ylabel="[mW/m²/str/nm]", yguidefontsize=8)

p3 = plot(1e7./ν, 100*(R[1,1,:].-RnoRS[1,1,:].+ieR[1,1,:])./(RnoRS[1,1,:]), linecolor=:red, linewidth=0.5, alpha=0.5, xlabel = "λ [nm]", xlims=(755,775))
p3 = plot!(1e7./ν, 100*(R[2,1,:].-RnoRS[2,1,:].+ieR[2,1,:])./(RnoRS[2,1,:]), linecolor=:blue, linewidth=0.5, alpha=0.5, xlabel = "λ [nm]", xlims=(755,775))
p3 = plot!(1e7./ν, 100*(R[3,1,:].-RnoRS[3,1,:].+ieR[3,1,:])./(RnoRS[3,1,:]), linecolor=:green, linewidth=0.5, alpha=0.5, xlabel = "λ [nm]", xlims=(755,775), ylims=(-1,5))
#p3 = plot!(1e7./ν, 100*(I_conv[1,:].-I_conv_noRS[1,:].+ieI_conv[1,:])./(I_conv_noRS[1,:]), linewidth = 2, linecolor=:red, xlims=(755,775))
#p3 = plot!(1e7./ν, 100*(I_conv[2,:].-I_conv_noRS[2,:].+ieI_conv[2,:])./(I_conv_noRS[2,:]), linewidth = 2, linecolor=:blue, xlims=(755,775))
#p3 = plot!(1e7./ν, 100*(I_conv[3,:].-I_conv_noRS[3,:].+ieI_conv[3,:])./(I_conv_noRS[3,:]), linewidth = 2, linecolor=:green, xlims=(755,775), ylims=(-1,5))

q1 = plot(1e7./ν, (RnoRS[1,2,:]).*convfct, linecolor=:red, linewidth=0.5, alpha=0.5, xlims=(755,775))
q1 = plot!(1e7./ν, (RnoRS[2,2,:]).*convfct, linecolor=:blue, linewidth=0.5, alpha=0.5, xlims=(755,775))
q1 = plot!(1e7./ν, (RnoRS[3,2,:]).*convfct, linecolor=:green, linewidth=0.5, alpha=0.5, xlims=(755,775))
#q1 = plot!(1e7./ν, (Q_conv_noRS[1,:]).*convfct, linewidth = 2, linecolor=:red, xlims=(755,775))
#q1 = plot!(1e7./ν, (Q_conv_noRS[2,:]).*convfct, linewidth = 2, linecolor=:blue, xlims=(755,775))
#q1 = plot!(1e7./ν, (Q_conv_noRS[3,:]).*convfct, linewidth = 2, linecolor=:green, xlims=(755,775))#, ylabel="[mW/m²/str/nm]", yguidefontsize=8)

q2 = plot(1e7./ν, (R[1,2,:].-RnoRS[1,2,:].+ieR[1,2,:]).*convfct, linecolor=:red, linewidth=0.5, alpha=0.5, xlims=(755,775))
q2 = plot!(1e7./ν, (R[2,2,:].-RnoRS[2,2,:].+ieR[2,2,:]).*convfct, linecolor=:blue, linewidth=0.5, alpha=0.5, xlims=(755,775))
q2 = plot!(1e7./ν, (R[3,2,:].-RnoRS[3,2,:].+ieR[3,2,:]).*convfct, linecolor=:green, linewidth=0.5, alpha=0.5, xlims=(755,775), ylims=(-0.1,0.18))
#q2 = plot!(1e7./ν, (Q_conv[1,:].-Q_conv_noRS[1,:].+ieQ_conv[1,:]).*convfct, linewidth = 2, linecolor=:red, xlims=(755,775))
#q2 = plot!(1e7./ν, (Q_conv[2,:].-Q_conv_noRS[2,:].+ieQ_conv[2,:]).*convfct, linewidth = 2, linecolor=:blue, xlims=(755,775))
#q2 = plot!(1e7./ν, (Q_conv[3,:].-Q_conv_noRS[3,:].+ieQ_conv[3,:]).*convfct, linewidth = 2, linecolor=:green, xlims=(755,775), ylims=(-0.1,0.18))#, ylabel="[mW/m²/str/nm]", yguidefontsize=8)

q3 = plot(1e7./ν, 100*(R[1,2,:].-RnoRS[1,2,:].+ieR[1,2,:])./(RnoRS[1,2,:]), linecolor=:red, linewidth=0.5, alpha=0.5, xlabel = "λ [nm]", xlims=(755,775))
q3 = plot!(1e7./ν, 100*(R[2,2,:].-RnoRS[2,2,:].+ieR[2,2,:])./(RnoRS[2,2,:]), linecolor=:blue, linewidth=0.5, alpha=0.5, xlabel = "λ [nm]", xlims=(755,775))
q3 = plot!(1e7./ν, 100*(R[3,2,:].-RnoRS[3,2,:].+ieR[3,2,:])./(RnoRS[3,2,:]), linecolor=:green, linewidth=0.5, alpha=0.5, xlabel = "λ [nm]", xlims=(755,775), ylims=(-0.1,0.5))
#q3 = plot!(1e7./ν, 100*(Q_conv[1,:].-Q_conv_noRS[1,:].+ieQ_conv[1,:])./(Q_conv_noRS[1,:]), linewidth = 3, linecolor=:red, xlims=(755,775))
#q3 = plot!(1e7./ν, 100*(Q_conv[2,:].-Q_conv_noRS[2,:].+ieQ_conv[2,:])./(Q_conv_noRS[2,:]), linewidth = 3, linecolor=:green, xlims=(755,775))
#q3 = plot!(1e7./ν, 100*(Q_conv[3,:].-Q_conv_noRS[3,:].+ieQ_conv[3,:])./(Q_conv_noRS[3,:]), linewidth = 3, linecolor=:blue, xlims=(755,775), ylims=(-0.1,0.5))


plot(p1, q1, p2, q2, p3, q3, layout = l, legend = false, title = ["I₀ (no RS)" "Q₀ (no RS)" "I₁-I₀" "Q₁-Q₀" "(I₁-I₀)/I₀ [%]" "(Q₁-Q₀)/Q₀ [%]"], titlefont = font(10))
#savefig("RingEffect_O2A_SZA30_wF.png")
savefig("RingEffect_O2A_SZA30_wF_hires.png")


convfct = 1e7./ν.^2   # to convert radiance units from mW/m²-str-cm⁻¹ to mW/m²-str-nm
#convfct = 1.0#./F₀ # to normalize wrt irradiance 
l = @layout [a1 a2; b1 b2]

p4 = plot(1e7./ν, (ieR[1,1,:]).*convfct, linecolor=:red, linewidth=0.5, alpha=0.5, xlims=(755,775))#, ylims=(-2e-6, 1.e-5))#, xlabel = "λ [nm]")
p4 = plot!(1e7./ν, (ieR[2,1,:]).*convfct, linecolor=:blue, linewidth=0.5, alpha=0.5, xlims=(755,775))#, xlabel = "λ [nm]")
p4 = plot!(1e7./ν, (ieR[3,1,:]).*convfct, linecolor=:green, linewidth=0.5, alpha=0.5, xlims=(755,775), ylims=(0,5.5), ylabel="[mW/m²/str/nm]", yguidefontsize=8)#, xlabel = "λ [nm]", xlims=(755,775))
#p4 = plot!(1e7./ν, (ieI_conv[1,:]).*convfct, linewidth = 2, linecolor=:red, xlims=(755,775))
#p4 = plot!(1e7./ν, (ieI_conv[2,:]).*convfct, linewidth = 2, linecolor=:blue, xlims=(755,775))
#p4 = plot!(1e7./ν, (ieI_conv[3,:]).*convfct, linewidth = 2, linecolor=:green, xlims=(755,775), ylims=(0,5.5), ylabel="[mW/m²/str/nm]", yguidefontsize=8)

p5 = plot(1e7./ν, 100*(ieR[1,1,:])./(RnoRS[1,1,:]), linecolor=:red, linewidth=0.5, alpha=0.5, xlabel = "λ [nm]", xlims=(755,775))
p5 = plot!(1e7./ν, 100*(ieR[2,1,:])./(RnoRS[2,1,:]), linecolor=:blue, linewidth=0.5, alpha=0.5, xlabel = "λ [nm]", xlims=(755,775))
p5 = plot!(1e7./ν, 100*(ieR[3,1,:])./(RnoRS[3,1,:]), linecolor=:green, linewidth=0.5, alpha=0.5, xlabel = "λ [nm]", xlims=(755,775), ylims=(0,8))
#p5 = plot!(1e7./ν, 100*(ieI_conv[1,:])./(I_conv_noRS[1,:]), linewidth = 2, linecolor=:red, xlims=(755,775))
#p5 = plot!(1e7./ν, 100*(ieI_conv[2,:])./(I_conv_noRS[2,:]), linewidth = 2, linecolor=:blue, xlims=(755,775))
#p5 = plot!(1e7./ν, 100*(ieI_conv[3,:])./(I_conv_noRS[3,:]), linewidth = 2, linecolor=:green, xlims=(755,775), ylims=(0,8))


q4 = plot(1e7./ν, (ieR[1,2,:]).*convfct, linecolor=:red, linewidth=0.5, alpha=0.5, xlims=(755,775))
q4 = plot!(1e7./ν, (ieR[2,2,:]).*convfct, linecolor=:blue, linewidth=0.5, alpha=0.5, xlims=(755,775))
q4 = plot!(1e7./ν, (ieR[3,2,:]).*convfct, linecolor=:green, linewidth=0.5, alpha=0.5, xlims=(755,775), ylims=(0,0.35))
#q4 = plot!(1e7./ν, (ieQ_conv[1,:]).*convfct, linewidth = 2, linecolor=:red, xlims=(755,775))
#q4 = plot!(1e7./ν, (ieQ_conv[2,:]).*convfct, linewidth = 2, linecolor=:blue, xlims=(755,775))
#q4 = plot!(1e7./ν, (ieQ_conv[3,:]).*convfct, linewidth = 2, linecolor=:green, xlims=(755,775), ylims=(0,0.35))#, ylabel="[mW/m²/str/nm]", yguidefontsize=8)

q5 = plot(1e7./ν, 100*(ieR[1,2,:])./(RnoRS[1,2,:]), linecolor=:red, linewidth=0.5, alpha=0.5, xlabel = "λ [nm]", xlims=(755,775))
q5 = plot!(1e7./ν, 100*(ieR[2,2,:])./(RnoRS[2,2,:]), linecolor=:blue, linewidth=0.5, alpha=0.5, xlabel = "λ [nm]", xlims=(755,775))
q5 = plot!(1e7./ν, 100*(ieR[3,2,:])./(RnoRS[3,2,:]), linecolor=:green, linewidth=0.5, alpha=0.5, xlabel = "λ [nm]", xlims=(755,775), ylims=(0,0.8))
#q5 = plot!(1e7./ν, 100*(ieQ_conv[1,:])./(Q_conv_noRS[1,:]), linewidth = 3, linecolor=:red, xlims=(755,775))
#q5 = plot!(1e7./ν, 100*(ieQ_conv[2,:])./(Q_conv_noRS[2,:]), linewidth = 3, linecolor=:green, xlims=(755,775))
#q5 = plot!(1e7./ν, 100*(ieQ_conv[3,:])./(Q_conv_noRS[3,:]), linewidth = 3, linecolor=:blue, xlims=(755,775), ylims=(0,0.8))


plot(p4, q4, p5, q5, layout = l, legend = false, title = ["ΔI₁" "ΔQ₁" "ΔI₁/I₀ [%]" "ΔQ₁/Q₀ [%]"], titlefont = font(10))
#savefig("RingSpectrum_O2A_SZA30_wF.png")
savefig("RingSpectrum_O2A_SZA30_wF_hires.png")




l = @layout [a1 a2; b1 b2; c1 c2]
p1 = plot(1e7./ν, 100*(1.0.-(RnoRS_ss[1,1,:]./RnoRS[1,1,:])), linecolor=:red, linewidth=0.5, alpha=0.5, xlims=(755,775))
p1 = plot!(1e7./ν, 100*(1.0.-(RnoRS_ss[2,1,:]./RnoRS[2,1,:])), linecolor=:blue, linewidth=0.5, alpha=0.5, xlims=(755,775))
p1 = plot!(1e7./ν, 100*(1.0.-(RnoRS_ss[3,1,:]./RnoRS[3,1,:])), linecolor=:green, linewidth=0.5, alpha=0.5, xlims=(755,775))

p2 = plot(1e7./ν, 100*(1.0.-(ieR_ss[1,1,:]./ieR[1,1,:])), linecolor=:red, linewidth=0.5, alpha=0.5, xlims=(755,775))
p2 = plot!(1e7./ν, 100*(1.0.-(ieR_ss[2,1,:]./ieR[2,1,:])), linecolor=:blue, linewidth=0.5, alpha=0.5, xlims=(755,775))
p2 = plot!(1e7./ν, 100*(1.0.-(ieR_ss[3,1,:]./ieR[3,1,:])), linecolor=:green, linewidth=0.5, alpha=0.5, xlims=(755,775))

#p2 = plot(1e7./ν, ((ieR[1,1,:]+R[1,1,:]-RnoRS[1,1,:]).-(ieR_ss[1,1,:]+R_ss[1,1,:]-RnoRS_ss[1,1,:])), linecolor=:red, linewidth=0.5, alpha=0.5, xlims=(755,775))
#p2 = plot!(1e7./ν, ((ieR[2,1,:]+R[2,1,:]-RnoRS[2,1,:]).-(ieR_ss[2,1,:]+R_ss[2,1,:]-RnoRS_ss[2,1,:])), linecolor=:blue, linewidth=0.5, alpha=0.5, xlims=(755,775))
#p2 = plot!(1e7./ν, ((ieR[3,1,:]+R[3,1,:]-RnoRS[3,1,:]).-(ieR_ss[3,1,:]+R_ss[3,1,:]-RnoRS_ss[3,1,:])), linecolor=:green, linewidth=0.5, alpha=0.5, ylabel="[mW/m²/str/nm]", yguidefontsize=8, xlims=(755,775), ylims=(-2.5,3.5))

p3 = plot(1e7./ν, 100*((ieR[1,1,:]+R[1,1,:]-RnoRS[1,1,:])./(RnoRS[1,1,:])), linecolor=:red, linewidth=0.5, alpha=0.5, xlims=(755,775))
p3 = plot!(1e7./ν, 100*((ieR_ss[1,1,:]+R_ss[1,1,:]-RnoRS_ss[1,1,:])./(RnoRS[1,1,:])), linecolor=:black, linewidth=0.5, alpha=0.5, xlims=(755,775))

p3 = plot!(1e7./ν, 10 .+ 100*((ieR[2,1,:]+R[2,1,:]-RnoRS[2,1,:])./(RnoRS[2,1,:])), linecolor=:blue, linewidth=0.5, alpha=0.5, xlims=(755,775))
p3 = plot!(1e7./ν, 10 .+ 100*((ieR_ss[2,1,:]+R_ss[2,1,:]-RnoRS_ss[2,1,:])./(RnoRS[2,1,:])), linecolor=:black, linewidth=0.5, alpha=0.5, xlims=(755,775))

p3 = plot!(1e7./ν, 20 .+ 100*((ieR[3,1,:]+R[3,1,:]-RnoRS[3,1,:])./(RnoRS[3,1,:])), linecolor=:green, linewidth=0.5, alpha=0.5, xlims=(755,775), xlabel="Wavelength, λ [nm]")
p3 = plot!(1e7./ν, 20 .+ 100*((ieR_ss[3,1,:]+R_ss[3,1,:]-RnoRS_ss[3,1,:])./(RnoRS[3,1,:])), linecolor=:black, linewidth=0.5, alpha=0.5, xlims=(755,775), xlabel="Wavelength, λ [nm]")


#p3 = plot(1e7./ν, 100*((ieR[1,1,:]+R[1,1,:]-RnoRS[1,1,:])./(RnoRS[1,1,:])), linecolor=:red, linewidth=0.5, alpha=0.5, xlims=(759,761), ylims=(-1,30))
#p3 = plot!(1e7./ν, 100*((ieR_ss[1,1,:]+R_ss[1,1,:]-RnoRS_ss[1,1,:])./(RnoRS[1,1,:])), linecolor=:black, linewidth=0.5, alpha=0.5, xlims=(759,761), ylims=(-1,30))

#p3 = plot!(1e7./ν, 10 .+ 100*((ieR[2,1,:]+R[2,1,:]-RnoRS[2,1,:])./(RnoRS[2,1,:])), linecolor=:blue, linewidth=0.5, alpha=0.5, xlims=(759,761), ylims=(-1,30))
#p3 = plot!(1e7./ν, 10 .+ 100*((ieR_ss[2,1,:]+R_ss[2,1,:]-RnoRS_ss[2,1,:])./(RnoRS[2,1,:])), linecolor=:black, linewidth=0.5, alpha=0.5, xlims=(759,761), ylims=(-1,30))

p4 = plot(1e7./ν, 100*((ieR[3,1,:]+R[3,1,:]-RnoRS[3,1,:])./(RnoRS[3,1,:])), linecolor=:green, linewidth=2, alpha=0.5, xlims=(759,761), ylims=(-1,2.8), xlabel="Wavelength, λ [nm]")
p4 = plot!(1e7./ν, 100*((ieR_ss[3,1,:]+R_ss[3,1,:]-RnoRS_ss[3,1,:])./(RnoRS[3,1,:])), linecolor=:black, linewidth=2, alpha=0.5, xlims=(759,761), ylims=(-1,2.8), xlabel="Wavelength, λ [nm]")


#p3 = plot(1e7./ν, 100*((ieR[1,1,:]+R[1,1,:]-RnoRS[1,1,:]).-(ieR_ss[1,1,:]+R_ss[1,1,:]-RnoRS_ss[1,1,:]))./(RnoRS[1,1,:]), linecolor=:red, linewidth=0.5, alpha=0.5, xlims=(755,775))
#p3 = plot!(1e7./ν, 100*((ieR[2,1,:]+R[2,1,:]-RnoRS[2,1,:]).-(ieR_ss[2,1,:]+R_ss[2,1,:]-RnoRS_ss[2,1,:]))./(RnoRS[2,1,:]), linecolor=:blue, linewidth=0.5, alpha=0.5, xlims=(755,775))
#p3 = plot!(1e7./ν, 100*((ieR[3,1,:]+R[3,1,:]-RnoRS[3,1,:]).-(ieR_ss[3,1,:]+R_ss[3,1,:]-RnoRS_ss[3,1,:]))./(RnoRS[3,1,:]), linecolor=:green, linewidth=0.5, alpha=0.5, xlims=(755,775), ylims=(0,20), xlabel="Wavelength, λ [nm]")


q1 = plot(1e7./ν, 100*(1.0.-(RnoRS_ss[1,2,:]./RnoRS[1,2,:])), linecolor=:red, linewidth=0.5, alpha=0.5, xlims=(755,775))
q1 = plot!(1e7./ν, 100*(1.0.-(RnoRS_ss[2,2,:]./RnoRS[2,2,:])), linecolor=:blue, linewidth=0.5, alpha=0.5, xlims=(755,775))
q1 = plot!(1e7./ν, 100*(1.0.-(RnoRS_ss[3,2,:]./RnoRS[3,2,:])), linecolor=:green, linewidth=0.5, alpha=0.5, xlims=(755,775))

q2 = plot(1e7./ν, 100*(1.0.-(ieR_ss[1,2,:]./ieR[1,2,:])), linecolor=:red, linewidth=0.5, alpha=0.5, xlims=(755,775))
q2 = plot!(1e7./ν, 100*(1.0.-(ieR_ss[2,2,:]./ieR[2,2,:])), linecolor=:blue, linewidth=0.5, alpha=0.5, xlims=(755,775))
q2 = plot!(1e7./ν, 100*(1.0.-(ieR_ss[3,2,:]./ieR[3,2,:])), linecolor=:green, linewidth=0.5, alpha=0.5, xlims=(755,775))

q3 = plot(1e7./ν, 100*((ieR[1,2,:]+R[1,2,:]-RnoRS[1,2,:])./(RnoRS[1,2,:])), linecolor=:red, linewidth=0.5, alpha=0.5, xlims=(755,775))
q3 = plot!(1e7./ν, 100*((ieR_ss[1,2,:]+R_ss[1,2,:]-RnoRS_ss[1,2,:])./(RnoRS[1,2,:])), linecolor=:black, linewidth=0.5, alpha=0.5, xlims=(755,775))

q3 = plot!(1e7./ν, 1 .+ 100*((ieR[2,2,:]+R[2,2,:]-RnoRS[2,2,:])./(RnoRS[2,2,:])), linecolor=:blue, linewidth=0.5, alpha=0.5, xlims=(755,775))
q3 = plot!(1e7./ν, 1 .+ 100*((ieR_ss[2,2,:]+R_ss[2,2,:]-RnoRS_ss[2,2,:])./(RnoRS[2,2,:])), linecolor=:black, linewidth=0.5, alpha=0.5, xlims=(755,775))

q3 = plot!(1e7./ν, 2 .+ 100*((ieR[3,2,:]+R[3,2,:]-RnoRS[3,2,:])./(RnoRS[3,2,:])), linecolor=:green, linewidth=0.5, alpha=0.5, xlims=(755,775), xlabel="Wavelength, λ [nm]")
q3 = plot!(1e7./ν, 2 .+ 100*((ieR_ss[3,2,:]+R_ss[3,2,:]-RnoRS_ss[3,2,:])./(RnoRS[3,2,:])), linecolor=:black, linewidth=0.5, alpha=0.5, xlims=(755,775), xlabel="Wavelength, λ [nm]")

#q3 = plot(1e7./ν, 100*((ieR[1,2,:]+R[1,2,:]-RnoRS[1,2,:])./(RnoRS[1,2,:])), linecolor=:red, linewidth=0.5, alpha=0.5, xlims=(759,761), ylims=(-0.1,3))
#q3 = plot!(1e7./ν, 100*((ieR_ss[1,2,:]+R_ss[1,2,:]-RnoRS_ss[1,2,:])./(RnoRS[1,2,:])), linecolor=:black, linewidth=0.5, alpha=0.5, xlims=(759,761), ylims=(-0.1,30))

#q3 = plot!(1e7./ν, 1 .+ 100*((ieR[2,2,:]+R[2,2,:]-RnoRS[2,2,:])./(RnoRS[2,2,:])), linecolor=:blue, linewidth=0.5, alpha=0.5, xlims=(759,761), ylims=(-0.1,3))
#q3 = plot!(1e7./ν, 1 .+ 100*((ieR_ss[2,2,:]+R_ss[2,2,:]-RnoRS_ss[2,2,:])./(RnoRS[2,2,:])), linecolor=:black, linewidth=0.5, alpha=0.5, xlims=(759,761), ylims=(-0.1,3))

q4 = plot(1e7./ν, 100*((ieR[3,2,:]+R[3,2,:]-RnoRS[3,2,:])./(RnoRS[3,2,:])), linecolor=:green, linewidth=2, alpha=0.5, xlims=(759,761), ylims=(-0.1,0.26), xlabel="Wavelength, λ [nm]")
q4 = plot!(1e7./ν, 100*((ieR_ss[3,2,:]+R_ss[3,2,:]-RnoRS_ss[3,2,:])./(RnoRS[3,2,:])), linecolor=:black, linewidth=2, alpha=0.5, xlims=(759,761), ylims=(-0.1,0.26), xlabel="Wavelength, λ [nm]")

#q2 = plot(1e7./ν, ((ieR[1,2,:]+R[1,2,:]-RnoRS[1,2,:]).-(ieR_ss[1,2,:]+R_ss[1,2,:]-RnoRS_ss[1,2,:])), linecolor=:red, linewidth=0.5, alpha=0.5, xlims=(755,775))
#q2 = plot!(1e7./ν, ((ieR[2,2,:]+R[2,2,:]-RnoRS[2,2,:]).-(ieR_ss[2,2,:]+R_ss[2,2,:]-RnoRS_ss[2,2,:])), linecolor=:blue, linewidth=0.5, alpha=0.5, xlims=(755,775))
#q2 = plot!(1e7./ν, ((ieR[3,2,:]+R[3,2,:]-RnoRS[3,2,:]).-(ieR_ss[3,2,:]+R_ss[3,2,:]-RnoRS_ss[3,2,:])), linecolor=:green, linewidth=0.5, alpha=0.5, xlims=(755,775), ylims=(-0.05,0.055))

#q3 = plot(1e7./ν, 100*(1.0.-(ieR_ss[1,2,:]+R_ss[1,2,:]-RnoRS_ss[1,2,:])./(ieR[1,2,:]+R[1,2,:]-RnoRS[1,2,:])), linecolor=:red, linewidth=0.5, alpha=0.5, xlims=(755,775))
#q3 = plot!(1e7./ν, 100*(1.0.-(ieR_ss[2,2,:]+R_ss[2,2,:]-RnoRS_ss[2,2,:])./(ieR[2,2,:]+R[2,2,:]-RnoRS[2,2,:])), linecolor=:blue, linewidth=0.5, alpha=0.5, xlims=(755,775))
#q3 = plot!(1e7./ν, 100*(1.0.-(ieR_ss[3,2,:]+R_ss[3,2,:]-RnoRS_ss[3,2,:])./(ieR[3,2,:]+R[3,2,:]-RnoRS[3,2,:])), linecolor=:green, linewidth=0.5, alpha=0.5, xlims=(755,775), ylims=(-3,10), xlabel="Wavelength, λ [nm]")

plot(p1, q1, p2, q2, p3, q3, layout = l, legend = false, title = ["(I₀-I'₀)/I₀ [%]"  "(Q₀-Q'₀)/Q₀ [%]" "(ΔI₁-ΔI'₁)/ΔI₁ [%]" "(ΔQ₁-ΔQ'₁)/ΔQ₁ [%]" "(I₁-I₀)/I₀ vs. (I'₁-I'₀)/I₀ [%]" "(Q₁-Q₀)/Q₀ vs. (Q'₁-Q'₀)/Q₀ [%]"], titlefont = font(10))
savefig("SingleScatteringApprox_O2A_SZA30_wF.png")

l = @layout [a1 a2]
plot(p4, q4, layout = l, legend = false, title = ["(I₁-I₀)/I₀ vs. (I'₁-I'₀)/I₀ [%]" "(Q₁-Q₀)/Q₀ vs. (Q'₁-Q'₀)/Q₀ [%]"], titlefont = font(10))
savefig("SingleScatteringApprox_O2A_SZA30_wF_zoom.png")

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

=#
#=
#for ctr=1:3
ν = MOM30_nors[:,1] 
νbin = bin30_nors[:,1] 

#No convolution is carried out in the following because the OCO sampling frequency is higher than the wavelength grid
I_noRS = zeros(3,length(ν))
I = zeros(3,length(ν))
ieI = zeros(3,length(ν))
Q_noRS = zeros(3,length(ν))
Q = zeros(3,length(ν))
ieQ = zeros(3,length(ν))

I_noRS_bin = zeros(3,length(νbin))
I_bin = zeros(3,length(νbin))
ieI_bin = zeros(3,length(νbin))
Q_noRS_bin = zeros(3,length(νbin))
Q_bin = zeros(3,length(νbin))
ieQ_bin = zeros(3,length(νbin))

I_noRS[1,:] .= MOM30_nors[:,2]#imfilter(RnoRS[ctr,1,:], kernel)
Q_noRS[1,:] .= MOM30_nors[:,3]
I[1,:] .= MOM30_rrs[:,2]#imfilter(R[ctr,1,:], kernel)
Q[1,:] .= MOM30_rrs[:,3]
ieI[1,:] .= MOM30_rrs[:,5]#imfilter(ieR[ctr,1,:], kernel)
ieQ[1,:] .= MOM30_rrs[:,6]

I_noRS[2,:] .= MOM45_nors[:,2]#imfilter(RnoRS[ctr,1,:], kernel)
Q_noRS[2,:] .= MOM45_nors[:,3]
I[2,:] .= MOM45_rrs[:,2]#imfilter(R[ctr,1,:], kernel)
Q[2,:] .= MOM45_rrs[:,3]
ieI[2,:] .= MOM45_rrs[:,5]#imfilter(ieR[ctr,1,:], kernel)
ieQ[2,:] .= MOM45_rrs[:,6]

I_noRS[3,:] .= MOM60_nors[:,2]#imfilter(RnoRS[ctr,1,:], kernel)
Q_noRS[3,:] .= MOM60_nors[:,3]
I[3,:] .= MOM60_rrs[:,2]#imfilter(R[ctr,1,:], kernel)
Q[3,:] .= MOM60_rrs[:,3]
ieI[3,:] .= MOM60_rrs[:,5]#imfilter(ieR[ctr,1,:], kernel)
ieQ[3,:] .= MOM60_rrs[:,6]


I_noRS_bin[1,:] .= bin30_nors[:,2]#imfilter(RnoRS[ctr,1,:], kernel)
Q_noRS_bin[1,:] .= bin30_nors[:,3]
I_bin[1,:] .= bin30_rrs[:,2]#imfilter(R[ctr,1,:], kernel)
Q_bin[1,:] .= bin30_rrs[:,3]
ieI_bin[1,:] .= bin30_rrs[:,5]#imfilter(ieR[ctr,1,:], kernel)
ieQ_bin[1,:] .= bin30_rrs[:,6]

I_noRS_bin[2,:] .= bin45_nors[:,2]#imfilter(RnoRS[ctr,1,:], kernel)
Q_noRS_bin[2,:] .= bin45_nors[:,3]
I_bin[2,:] .= bin45_rrs[:,2]#imfilter(R[ctr,1,:], kernel)
Q_bin[2,:] .= bin45_rrs[:,3]
ieI_bin[2,:] .= bin45_rrs[:,5]#imfilter(ieR[ctr,1,:], kernel)
ieQ_bin[2,:] .= bin45_rrs[:,6]

I_noRS_bin[3,:] .= bin60_nors[:,2]#imfilter(RnoRS[ctr,1,:], kernel)
Q_noRS_bin[3,:] .= bin60_nors[:,3]
I_bin[3,:] .= bin60_rrs[:,2]#imfilter(R[ctr,1,:], kernel)
Q_bin[3,:] .= bin60_rrs[:,3]
ieI_bin[3,:] .= bin60_rrs[:,5]#imfilter(ieR[ctr,1,:], kernel)
ieQ_bin[3,:] .= bin60_rrs[:,6]
#end
#I_conv = InstrumentOperator.conv_spectra(kernel, )

convfct = 1e7./ν.^2   # to convert radiance units from mW/m²-str-cm⁻¹ to mW/m²-str-nm
#convfct = 1.0#./F₀ # to normalize wrt irradiance 
l = @layout [a1 a2 ; b1 b2; c1 c2]
p1 = plot(1e7./ν, 100*(I[1,:].-I_noRS[1,:].+ieI[1,:])./(I_noRS[1,:]), linecolor=:black, linewidth=0.8, alpha=0.85, xlabel = "λ [nm]", xlims=(755,775), ylims=(-1,5.5))
p1 = plot!(1e7./νbin, 100*(I_bin[1,:].-I_noRS_bin[1,:].+ieI_bin[1,:])./(I_noRS_bin[1,:]), linecolor=:red, linewidth=0.8, alpha=0.85, xlabel = "λ [nm]", xlims=(755,775))

p2 = plot(1e7./ν, 100*(I[2,:].-I_noRS[2,:].+ieI[2,:])./(I_noRS[2,:]), linecolor=:black, linewidth=0.8, alpha=0.85, xlabel = "λ [nm]", xlims=(755,775), ylims=(-1,5.5))
p2 = plot!(1e7./νbin, 100*(I_bin[2,:].-I_noRS_bin[2,:].+ieI_bin[2,:])./(I_noRS_bin[2,:]), linecolor=:red, linewidth=0.8, alpha=0.85, xlabel = "λ [nm]", xlims=(755,775))

p3 = plot(1e7./ν, 100*(I[3,:].-I_noRS[3,:].+ieI[3,:])./(I_noRS[3,:]), linecolor=:black, linewidth=0.8, alpha=0.85, xlabel = "λ [nm]", xlims=(755,775), ylims=(-1,5.5))
p3 = plot!(1e7./νbin, 100*(I_bin[3,:].-I_noRS_bin[3,:].+ieI_bin[3,:])./(I_noRS_bin[3,:]), linecolor=:red, linewidth=0.8, alpha=0.85, xlabel = "λ [nm]", xlims=(755,775))

q1 = plot(1e7./ν, 100*(Q[1,:].-Q_noRS[1,:].+ieQ[1,:])./(Q_noRS[1,:]), linecolor=:black, linewidth=0.8, alpha=0.85, xlabel = "λ [nm]", xlims=(755,775), ylims=(-0.1,1))
q1 = plot!(1e7./νbin, 100*(Q_bin[1,:].-Q_noRS_bin[1,:].+ieQ_bin[1,:])./(Q_noRS_bin[1,:]), linecolor=:red, linewidth=0.8, alpha=0.85, xlabel = "λ [nm]", xlims=(755,775))

q2 = plot(1e7./ν, 100*(Q[2,:].-Q_noRS[2,:].+ieQ[2,:])./(Q_noRS[2,:]), linecolor=:black, linewidth=0.8, alpha=0.85, xlabel = "λ [nm]", xlims=(755,775), ylims=(-0.1,1))
q2 = plot!(1e7./νbin, 100*(Q_bin[2,:].-Q_noRS_bin[2,:].+ieQ_bin[2,:])./(Q_noRS_bin[2,:]), linecolor=:red, linewidth=0.8, alpha=0.85, xlabel = "λ [nm]", xlims=(755,775))

q3 = plot(1e7./ν, 100*(Q[3,:].-Q_noRS[3,:].+ieQ[3,:])./(Q_noRS[3,:]), linecolor=:black, linewidth=0.8, alpha=0.85, xlabel = "λ [nm]", xlims=(755,775), ylims=(-0.1,1))
q3 = plot!(1e7./νbin, 100*(Q_bin[3,:].-Q_noRS_bin[3,:].+ieQ_bin[3,:])./(Q_noRS_bin[3,:]), linecolor=:red, linewidth=0.8, alpha=0.85, xlabel = "λ [nm]", xlims=(755,775))


plot(p1, q1, p2, q2, p3, q3, layout = l, legend = false, title = ["(I₁-I₀)/I₀ [%], SZA=30ᵒ" "(Q₁-Q₀)/Q₀ [%], SZA=30ᵒ" "(I₁-I₀)/I₀ [%], SZA=45ᵒ" "(Q₁-Q₀)/Q₀ [%], SZA=45ᵒ" "(I₁-I₀)/I₀ [%], SZA=60ᵒ" "(Q₁-Q₀)/Q₀ [%], SZA=60ᵒ"], titlefont = font(10))
#savefig("RingEffect_O2A_SZA30_wF.png")
savefig("RingEffect_MOM_vs_bin_2.png")

l = @layout [a1; a2]# ; b1 b2; c1 c2]
p3 = plot(1e7./ν, 100*(I[3,:].-I_noRS[3,:].+ieI[3,:])./(I_noRS[3,:]), linecolor=:black, linewidth=1, alpha=0.85, xlabel = "λ [nm]", xlims=(755,775), ylims=(-1,5.5))
#p3 = plot!(1e7./νbin, 100*(I_bin[3,:].-I_noRS_bin[3,:].+ieI_bin[3,:])./(I_noRS_bin[3,:]), linecolor=:red, linewidth=1, alpha=0.85, xlabel = "λ [nm]", xlims=(755,775))

q3 = plot(1e7./ν, 100*(Q[3,:].-Q_noRS[3,:].+ieQ[3,:])./(Q_noRS[3,:]), linecolor=:black, linewidth=1, alpha=0.85, xlabel = "λ [nm]", xlims=(755,775), ylims=(-0.1,1))
#q3 = plot!(1e7./νbin, 100*(Q_bin[3,:].-Q_noRS_bin[3,:].+ieQ_bin[3,:])./(Q_noRS_bin[3,:]), linecolor=:red, linewidth=1, alpha=0.85, xlabel = "λ [nm]", xlims=(755,775))


plot(p3, q3, layout = l, legend = false, title = ["(I₁-I₀)/I₀ [%], SZA=60ᵒ" "(Q₁-Q₀)/Q₀ [%], SZA=60ᵒ"], titlefont = font(10))
#savefig("RingEffect_O2A_SZA30_wF.png")
savefig("RingEffect_MOM_vs_bin__sza60_1.png")

l = @layout [a1; a2]# ; b1 b2; c1 c2]

p3 = plot(1e7./ν, 100*(I[3,:].-I_noRS[3,:].+ieI[3,:])./(I_noRS[3,:]), linecolor=:black, linewidth=2, alpha=0.85, xlabel = "λ [nm]", xlims=(759,762), ylims=(-1,5.5))
p3 = plot!(1e7./νbin, 100*(I_bin[3,:].-I_noRS_bin[3,:].+ieI_bin[3,:])./(I_noRS_bin[3,:]), linecolor=:red, linewidth=2, alpha=0.85, xlabel = "λ [nm]", xlims=(759,762))

q3 = plot(1e7./ν, 100*(Q[3,:].-Q_noRS[3,:].+ieQ[3,:])./(Q_noRS[3,:]), linecolor=:black, linewidth=2, alpha=0.85, xlabel = "λ [nm]", xlims=(759,762), ylims=(-0.1,1))
q3 = plot!(1e7./νbin, 100*(Q_bin[3,:].-Q_noRS_bin[3,:].+ieQ_bin[3,:])./(Q_noRS_bin[3,:]), linecolor=:red, linewidth=2, alpha=0.85, xlabel = "λ [nm]", xlims=(759,762))

plot(p3, q3, layout = l, legend = false, title = ["(I₁-I₀)/I₀ [%], SZA=60ᵒ" "(Q₁-Q₀)/Q₀ [%], SZA=60ᵒ"], titlefont = font(10))
#savefig("RingEffect_O2A_SZA30_wF.png")
savefig("RingEffect_MOM_vs_bin_sza60_Rbranch_2.png")

l = @layout [a1; a2]# ; b1 b2; c1 c2]

p3 = plot(1e7./ν, 100*(I[3,:].-I_noRS[3,:].+ieI[3,:])./(I_noRS[3,:]), linecolor=:black, linewidth=2, alpha=0.85, xlabel = "λ [nm]", xlims=(769,772), ylims=(-1,5.5))
p3 = plot!(1e7./νbin, 100*(I_bin[3,:].-I_noRS_bin[3,:].+ieI_bin[3,:])./(I_noRS_bin[3,:]), linecolor=:red, linewidth=2, alpha=0.85, xlabel = "λ [nm]", xlims=(769,772))

q3 = plot(1e7./ν, 100*(Q[3,:].-Q_noRS[3,:].+ieQ[3,:])./(Q_noRS[3,:]), linecolor=:black, linewidth=2, alpha=0.85, xlabel = "λ [nm]", xlims=(769,772), ylims=(-0.1,1))
q3 = plot!(1e7./νbin, 100*(Q_bin[3,:].-Q_noRS_bin[3,:].+ieQ_bin[3,:])./(Q_noRS_bin[3,:]), linecolor=:red, linewidth=2, alpha=0.85, xlabel = "λ [nm]", xlims=(769,772))

plot(p3, q3, layout = l, legend = false, title = ["(I₁-I₀)/I₀ [%], SZA=60ᵒ" "(Q₁-Q₀)/Q₀ [%], SZA=60ᵒ"], titlefont = font(10))
#savefig("RingEffect_O2A_SZA30_wF.png")
savefig("RingEffect_MOM_vs_bin_sza60_farPbranch_2.png")

l = @layout [a1; a2]# ; b1 b2; c1 c2]

p3 = plot(1e7./ν, 100*(I[3,:].-I_noRS[3,:].+ieI[3,:])./(I_noRS[3,:]), linecolor=:black, linewidth=2, alpha=0.85, xlabel = "λ [nm]", xlims=(762,767), ylims=(-1,5.5))
p3 = plot!(1e7./νbin, 100*(I_bin[3,:].-I_noRS_bin[3,:].+ieI_bin[3,:])./(I_noRS_bin[3,:]), linecolor=:red, linewidth=2, alpha=0.85, xlabel = "λ [nm]", xlims=(762,767))

q3 = plot(1e7./ν, 100*(Q[3,:].-Q_noRS[3,:].+ieQ[3,:])./(Q_noRS[3,:]), linecolor=:black, linewidth=2, alpha=0.85, xlabel = "λ [nm]", xlims=(762,767), ylims=(-0.1,1))
q3 = plot!(1e7./νbin, 100*(Q_bin[3,:].-Q_noRS_bin[3,:].+ieQ_bin[3,:])./(Q_noRS_bin[3,:]), linecolor=:red, linewidth=2, alpha=0.85, xlabel = "λ [nm]", xlims=(762,767))

plot(p3, q3, layout = l, legend = false, title = ["(I₁-I₀)/I₀ [%], SZA=60ᵒ" "(Q₁-Q₀)/Q₀ [%], SZA=60ᵒ"], titlefont = font(10))
#savefig("RingEffect_O2A_SZA30_wF.png")
savefig("RingEffect_MOM_vs_bin_sza60_Pbranch_2.png")

I_noRS_lp = zeros(length(ν))
I_lp = zeros(length(ν))
ieI_lp = zeros(length(ν))
Q_noRS_lp = zeros(length(ν))
Q_lp = zeros(length(ν))
ieQ_lp = zeros(length(ν))

I_noRS_lp[:] .= MOM60_nors_lp[:,2]#imfilter(RnoRS[ctr,1,:], kernel)
Q_noRS_lp[:] .= MOM60_nors_lp[:,3]
I_lp[:] .= MOM60_rrs_lp[:,2]#imfilter(R[ctr,1,:], kernel)
Q_lp[:] .= MOM60_rrs_lp[:,3]
ieI_lp[:] .= MOM60_rrs_lp[:,5]#imfilter(ieR[ctr,1,:], kernel)
ieQ_lp[:] .= MOM60_rrs_lp[:,6]

l = @layout [a1; a2]# ; b1 b2; c1 c2]
p3 = plot(1e7./ν, 100*(I[3,:].-I_noRS[3,:].+ieI[3,:])./(I_noRS[3,:]), linecolor=:black, linewidth=1, alpha=0.85, xlabel = "λ [nm]", xlims=(755,775), ylims=(-1,5.5))
p3 = plot!(1e7./νbin, 100*(I_bin[3,:].-I_noRS_bin[3,:].+ieI_bin[3,:])./(I_noRS_bin[3,:]), linecolor=:red, linewidth=1, alpha=0.85, xlabel = "λ [nm]", xlims=(755,775))
p3 = plot!(1e7./ν, 100*(-I_lp .+ I[3,:])./I[3,:], linecolor=:blue, linewidth=1, alpha=0.5, xlabel = "λ [nm]", xlims=(755,775))

q3 = plot(1e7./ν, 100*(Q[3,:].-Q_noRS[3,:].+ieQ[3,:])./(Q_noRS[3,:]), linecolor=:black, linewidth=1, alpha=0.85, xlabel = "λ [nm]", xlims=(755,775), ylims=(-0.1,1))
#q3 = plot!(1e7./νbin, 100*(Q_bin[3,:].-Q_noRS_bin[3,:].+ieQ_bin[3,:])./(Q_noRS_bin[3,:]), linecolor=:red, linewidth=1, alpha=0.85, xlabel = "λ [nm]", xlims=(755,775))


plot(p3, q3, layout = l, legend = false, title = ["(I₁-I₀)/I₀ [%], SZA=60ᵒ" "(Q₁-Q₀)/Q₀ [%], SZA=60ᵒ"], titlefont = font(10))

=#




