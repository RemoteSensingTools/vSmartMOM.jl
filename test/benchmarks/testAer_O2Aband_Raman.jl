##
using CUDA
device!(0)
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
## /net/fluo/data2/software/Julia/julia-1.11.3/bin/julia --project=/home/sanghavi/code/github/vSmartMOM.jl/ /home/sanghavi/code/github/vSmartMOM.jl/test/benchmarks/testAer_O2Aband_Raman.jl &

FT = Float64

# Load YAML files into parameter struct
#parameters = 
#    parameters_from_yaml("/home/sanghavi/code/github/vSmartMOM.jl/test/test_parameters/testAer_O2Aband_Raman.yaml");
parameters = 
    parameters_from_yaml("/home/sanghavi/code/github/vSmartMOM.jl/test/test_parameters/testAer_O2Bband_Raman.yaml");
#parameters = parameters_from_yaml("test/test_parameters/O2Parameters2.yaml");
# Create model struct (precomputes optical properties) from parameters
#model      = model_from_parameters(parameters);
#=a=[]
for i=0:20
    push!(a, acosd(i/20))  
end
sza=reverse(Int.(ceil.(a[8:21])))=#

#======================================================================================================================#


aer_r₀ = [0.1;0.5]
τ_ref = [0.1; 0.5]
aer_z₀ = [2;12]
AOD_str= ["0p1";"0p5"]
z_str  = ["2";"12"]
r_str  = ["0p1";"0p5"]
sza=[0; 50; 70]
#ρ = zeros(FT,5) #ρ = zeros(FT,15)
ρ = [0.1]#[0, 0.5, 1.0]
ρ_str = []

for iρ = 1:length(ρ) #15
    #ρ[iρ] = (iρ)*0.1 #ρ[iρ] = (iρ-1)*0.1
    #for i=1:length(parameters.spec_bands)
    #    parameters.brdf[i].albedo = ρ[iρ]
    #end
    push!(ρ_str, replace(string(round(ρ[iρ], digits=2)),"."=>"p"))
end
psurf=[1000]


#psurf=[1000 750 500]
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

for iτ = 1:length(τ_ref)
    parameters.scattering_params.rt_aerosols[1].τ_ref = τ_ref[iτ]
    for iz = 1:length(aer_z₀)
        parameters.scattering_params.rt_aerosols[1].z₀ = aer_z₀[iz]
        for ir = 1:length(aer_r₀)
            parameters.scattering_params.rt_aerosols[1].aerosol.size_distribution = 
                                                            LogNormal(log(aer_r₀[ir]),1.12)
            
            model      = model_from_parameters(parameters);   
            for isurf = 1:1 # 3:3 # 2:2 # 1:1 # 
                for iρ = 1:length(ρ) #21 #1:15 #3 #1:15
                    for iA = 1:length(sza) #14
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
                        overlap_ν = 200 #250 #230
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
                            
                            F₀ = Tsolar_interp.(ν) .* P;
                            SIF₀ .= 0.0 # = SIF_interp.(ν); # 
                            RS_type0.F₀[1,:] = F₀; #1.0 #
                            RS_type1.F₀[1,:] = F₀;
                            RS_type0.SIF₀[1,:] = SIF₀; #1.0 #
                            RS_type1.SIF₀[1,:] = SIF₀;
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
                            #fname0 = "/home/sanghavi/data/RamanSIFgrid/aerosol_test/aer_AOD0p0_sza"*sza_str*"_vza"*vza_str*"_vaz"*vaz_str*"_alb"*albedo*"_psurf"*string(psurf[isurf])*"hpa_nors_ABO2.dat"
                            #fname1 = "/home/sanghavi/data/RamanSIFgrid/aerosol_test/aer_AOD0p0_sza"*sza_str*"_vza"*vza_str*"_vaz"*vaz_str*"_alb"*albedo*"_psurf"*string(psurf[isurf])*"hpa_rrs_ABO2.dat"
                            # No SIF
                            #fname0 = "/home/sanghavi/data/RamanSIFgrid/aerosol_test/aer_AOD"*AOD_str[iτ]*"_z"*z_str[iz]*"_r"*r_str[ir]*"_sza"*sza_str*"_vza"*vza_str*"_vaz"*vaz_str*"_alb"*albedo*"_psurf"*string(psurf[isurf])*"hpa_nors_ABO2.dat"
                            #fname1 = "/home/sanghavi/data/RamanSIFgrid/aerosol_test/aer_AOD"*AOD_str[iτ]*"_z"*z_str[iz]*"_r"*r_str[ir]*"_sza"*sza_str*"_vza"*vza_str*"_vaz"*vaz_str*"_alb"*albedo*"_psurf"*string(psurf[isurf])*"hpa_rrs_ABO2.dat"
                            fname0 = "/home/sanghavi/data/RamanSIFgrid/aerosol_test/aer_AOD"*AOD_str[iτ]*"_z"*z_str[iz]*"_r"*r_str[ir]*"_sza"*sza_str*"_vza"*vza_str*"_vaz"*vaz_str*"_alb"*albedo*"_psurf"*string(psurf[isurf])*"hpa_nors_BBO2.dat"
                            fname1 = "/home/sanghavi/data/RamanSIFgrid/aerosol_test/aer_AOD"*AOD_str[iτ]*"_z"*z_str[iz]*"_r"*r_str[ir]*"_sza"*sza_str*"_vza"*vza_str*"_vaz"*vaz_str*"_alb"*albedo*"_psurf"*string(psurf[isurf])*"hpa_rrs_BBO2.dat"
                            
                            aer_rrs_ABO2 = []
                            aer_nors_ABO2 = []
                            
                            aer_rrs_ABO2 = [tot_ν R[vctr,1,:] R[vctr,2,:] R[vctr,3,:] ieR[vctr,1,:] ieR[vctr,2,:] ieR[vctr,3,:] tot_F₀]
                            aer_nors_ABO2 = [tot_ν RnoRS[vctr,1,:] RnoRS[vctr,2,:] RnoRS[vctr,3,:] tot_F₀]

                            writedlm(fname0, aer_nors_ABO2)
                            writedlm(fname1, aer_rrs_ABO2)
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
        end
    end
end

#=================================================================================================#
#aer_r₀ = [0.1;0.5]
#=
τ_ref = [0.0]
#aer_z₀ = [2;12]
AOD_str= ["0p0"]
#z_str  = ["2";"12"]
#r_str  = ["0p1";"0p5"]
sza=[0; 50; 70]
#ρ = zeros(FT,3) #ρ = zeros(FT,15)
#ρ_str = []
#=for iτ = 1:length(τ_ref)
    parameters.scattering_params.rt_aerosols[1].τ_ref = τ_ref[iτ]
    for iz = 1:length(aer_z₀)
        parameters.scattering_params.rt_aerosols[1].z₀ = aer_z₀[iz]
        for ir = 1:length(aer_r₀) # 3:3 # 2:2 # 1:1 #
            parameters.scattering_params.rt_aerosols[1].aerosol.size_distribution = 
                                                            LogNormal(log(aer_r₀[ir]),1.12)
            model      = model_from_parameters(parameters);                
        end
    end
end=#
#=for iρ = 1:3 #15
    ρ[iρ] = (iρ-1)*0.5
    #for i=1:length(parameters.spec_bands)
    #    parameters.brdf[i].albedo = ρ[iρ]
    #end
    push!(ρ_str, replace(string(round(ρ[iρ], digits=2)),"."=>"p"))
end=#
psurf=[1000]
#psurf=[1000 750 500]
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

for iτ = 1:length(τ_ref)
    parameters.scattering_params.rt_aerosols[1].τ_ref = τ_ref[iτ]
    #for iz = 1:length(aer_z₀)
    #    parameters.scattering_params.rt_aerosols[1].z₀ = aer_z₀[iz]
    #    for ir = 1:length(aer_r₀)
    #        parameters.scattering_params.rt_aerosols[1].aerosol.size_distribution = 
    #                                                        LogNormal(log(aer_r₀[ir]),1.12)
            
    model      = model_from_parameters(parameters);   
    for isurf = 1:1 # 3:3 # 2:2 # 1:1 # 
        for iρ = 1:length(ρ) #21 #1:15 #3 #1:15
            for iA = 1:length(sza) #3 #14
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
                overlap_ν = 200 #250 #230
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
                    #= for i=1:length(P)
                        sol_trans = Tsolar_interp(ν[i]);
                        F₀[i] = sol_trans * P[i];
                        SIF₀[i] = SIF_interp(ν[i]); #0.0
                        RS_type0.F₀[1,i] = F₀[i]; #1.0 #
                        RS_type1.F₀[1,i] = F₀[i];
                        RS_type0.SIF₀[1,i] = SIF₀[i]; #1.0 #
                        RS_type1.SIF₀[1,i] = SIF₀[i];
                    end =#
                    F₀ = Tsolar_interp.(ν) .* P;
                    SIF₀ .= 0.0 # = SIF_interp.(ν); # 
                    RS_type0.F₀[1,:] = F₀; #1.0 #
                    RS_type1.F₀[1,:] = F₀;
                    RS_type0.SIF₀[1,:] = SIF₀; #1.0 #
                    RS_type1.SIF₀[1,:] = SIF₀;
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
                    #fname0 = "/home/sanghavi/data/RamanSIFgrid/aerosol_test/aer_AOD0p0_sza"*sza_str*"_vza"*vza_str*"_vaz"*vaz_str*"_alb"*albedo*"_psurf"*string(psurf[isurf])*"hpa_nors_ABO2.dat"
                    #fname1 = "/home/sanghavi/data/RamanSIFgrid/aerosol_test/aer_AOD0p0_sza"*sza_str*"_vza"*vza_str*"_vaz"*vaz_str*"_alb"*albedo*"_psurf"*string(psurf[isurf])*"hpa_rrs_ABO2.dat"
                    # No SIF
                    #fname0 = "/home/sanghavi/data/RamanSIFgrid/aerosol_test/aer_AOD"*AOD_str[iτ]*"_sza"*sza_str*"_vza"*vza_str*"_vaz"*vaz_str*"_alb"*albedo*"_psurf"*string(psurf[isurf])*"hpa_nors_ABO2.dat"
                    #fname1 = "/home/sanghavi/data/RamanSIFgrid/aerosol_test/aer_AOD"*AOD_str[iτ]*"_sza"*sza_str*"_vza"*vza_str*"_vaz"*vaz_str*"_alb"*albedo*"_psurf"*string(psurf[isurf])*"hpa_rrs_ABO2.dat"
                    fname0 = "/home/sanghavi/data/RamanSIFgrid/aerosol_test/aer_AOD"*AOD_str[iτ]*"_sza"*sza_str*"_vza"*vza_str*"_vaz"*vaz_str*"_alb"*albedo*"_psurf"*string(psurf[isurf])*"hpa_nors_BBO2.dat"
                    fname1 = "/home/sanghavi/data/RamanSIFgrid/aerosol_test/aer_AOD"*AOD_str[iτ]*"_sza"*sza_str*"_vza"*vza_str*"_vaz"*vaz_str*"_alb"*albedo*"_psurf"*string(psurf[isurf])*"hpa_rrs_BBO2.dat"
                    
                    aer_rrs_ABO2 = []
                    aer_nors_ABO2 = []
                    
                    aer_rrs_ABO2 = [tot_ν R[vctr,1,:] R[vctr,2,:] R[vctr,3,:] ieR[vctr,1,:] ieR[vctr,2,:] ieR[vctr,3,:] tot_F₀]
                    aer_nors_ABO2 = [tot_ν RnoRS[vctr,1,:] RnoRS[vctr,2,:] RnoRS[vctr,3,:] tot_F₀]

                    writedlm(fname0, aer_nors_ABO2)
                    writedlm(fname1, aer_rrs_ABO2)
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
end

=#

