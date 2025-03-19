##
#using CUDA
#device!(1)
using Revise
#using vSmartMOM, vSmartMOM.CoreRT, vSmartMOM.SolarModel
#using vSmartMOM.InelasticScattering
using Statistics
using Interpolations
using InstrumentOperator #for convolution of hires spectrum to instrument grid
using ImageFiltering
using Distributions
#using Plots
using DelimitedFiles
#
## run this code using the following command:
## /net/fluo/data2/software/Julia/julia-1.11.1/bin/julia --project=/home/sanghavi/code/github/vSmartMOM.jl/ /home/sanghavi/code/github/vSmartMOM.jl/test/benchmarks/testAer_O2Aband_Raman.jl &

FT = Float64

# Load YAML files into parameter struct
#parameters = 
#    parameters_from_yaml("/home/sanghavi/code/github/vSmartMOM.jl/test/test_parameters/testAer_O2Aband_Raman.yaml");
#parameters = parameters_from_yaml("test/test_parameters/O2Parameters2.yaml");
# Create model struct (precomputes optical properties) from parameters
#model      = model_from_parameters(parameters);
#=a=[]
for i=0:20
    push!(a, acosd(i/20))  
end
sza=reverse(Int.(ceil.(a[8:21])))=#
aer_r₀ = [0.1;0.5]
τ_ref = [0.1; 0.5]
aer_z₀ = [2;12]
AOD_str= ["0p1";"0p5"]
z_str  = ["2";"12"]
r_str  = ["0p1";"0p5"]
sza=[0; 50; 70]
#ρ = zeros(FT,3) #ρ = zeros(FT,15)
#ρ_str = []

#for iρ = 1:3 #15
#    ρ[iρ] = (iρ-1)*0.5
#    push!(ρ_str, replace(string(round(ρ[iρ], digits=2)),"."=>"p"))
#end
#ρ = [0.0, 0.1, 0.5, 1.0]
#ρ_str = ["0p0", "0p1", "0p5", "1.0"]
ρ = [0.0, 0.1, 0.5, 1.0] #[0.0, 0.1, 0.2, 0.3, 0.4, 0.5]
ρ_str = ["0p0", "0p1", "0p5", "1p0"] #["0p0", "0p1", "0p2", "0p3", "0p4", "0p5"]
psurf=[1000]
#psurf=[1000 750 500]
#J_SIF = readdlm("/home/sanghavi/code/github/vSmartMOM.jl/src/SIF_emission/sif-spectra.csv", ',') 
#J_SIF = readdlm("/home/sanghavi/code/github/vSmartMOM.jl/src/SIF_emission/PC1_SIFSpectra_allSpecies.csv", ',') 
#νSIF = reverse(1e7./J_SIF[2:end,1])
#jSIF = reverse(J_SIF[2:end,2]*(0.5*π/maximum(J_SIF[2:end,2]))) # mW/m²/nm
#jSIF .*= 1e7./(νSIF.^2) # mW/m²/cm⁻¹
#SIF_interp = LinearInterpolation(νSIF, jSIF)
#Tsolar = solar_transmission_from_file("/home/sanghavi/code/github/vSmartMOM.jl/src/SolarModel/solar.out")
#Tsolar_interp = LinearInterpolation(Tsolar[4:end, 1], Tsolar[4:end, 2])
#n_bands = length(parameters.spec_bands);
#T_sun = 5777. # K

# O2 A-band
fname0 = "/home/sanghavi/data/RamanSIFgrid/aerosol_test/aer_AOD0p1_z2_r0p1_sza0_vza0p0_vaz0p0_alb0p0_psurf1000hpa_nors_ABO2.dat"
aer_nors_ABO2 = unique(readdlm(fname0), dims=1)
    
ν = aer_nors_ABO2[:,1]
λ = 1e7./reverse(ν)
F₀ = reverse(aer_nors_ABO2[:,5]).* (1e7./λ.^2) #mW/m^2/s/nm

Irrs = zeros(length(τ_ref),length(aer_z₀),length(aer_r₀),length(ρ),length(sza),length(λ))
Qrrs = zeros(length(τ_ref),length(aer_z₀),length(aer_r₀),length(ρ),length(sza),length(λ))
Urrs = zeros(length(τ_ref),length(aer_z₀),length(aer_r₀),length(ρ),length(sza),length(λ))
ieIrrs = zeros(length(τ_ref),length(aer_z₀),length(aer_r₀),length(ρ),length(sza),length(λ))
ieQrrs = zeros(length(τ_ref),length(aer_z₀),length(aer_r₀),length(ρ),length(sza),length(λ))
ieUrrs = zeros(length(τ_ref),length(aer_z₀),length(aer_r₀),length(ρ),length(sza),length(λ))
Inors = zeros(length(τ_ref),length(aer_z₀),length(aer_r₀),length(ρ),length(sza),length(λ))
Qnors = zeros(length(τ_ref),length(aer_z₀),length(aer_r₀),length(ρ),length(sza),length(λ))
Unors = zeros(length(τ_ref),length(aer_z₀),length(aer_r₀),length(ρ),length(sza),length(λ))

for iτ = 1:length(τ_ref)
    for iz = 1:length(aer_z₀)
        for ir = 1:length(aer_r₀)
            for isurf = 1:1 # 3:3 # 2:2 # 1:1 # 
                for iρ = 1:length(ρ) #21 #1:15 #3 #1:15
                    for iA = 1:length(sza) #14

                        vza_str = "0p0"
                        albedo = ρ_str[iρ]
                        sza_str = string(sza[iA])
                        vaz_str = "0p0"
                            
                        # No SIF
                        fname0 = "/home/sanghavi/data/RamanSIFgrid/aerosol_test/aer_AOD"*AOD_str[iτ]*"_z"*z_str[iz]*"_r"*r_str[ir]*"_sza"*sza_str*"_vza"*vza_str*"_vaz"*vaz_str*"_alb"*albedo*"_psurf"*string(psurf[isurf])*"hpa_nors_ABO2.dat"
                        fname1 = "/home/sanghavi/data/RamanSIFgrid/aerosol_test/aer_AOD"*AOD_str[iτ]*"_z"*z_str[iz]*"_r"*r_str[ir]*"_sza"*sza_str*"_vza"*vza_str*"_vaz"*vaz_str*"_alb"*albedo*"_psurf"*string(psurf[isurf])*"hpa_rrs_ABO2.dat"
                            
                        aer_rrs_ABO2  = unique(readdlm(fname1), dims=1)
                        aer_nors_ABO2 = unique(readdlm(fname0), dims=1)
                            
                        #ν = aer_rrs_ABO2[:,1]
                        #λ = 1e7./reverse(ν)
                        Irrs[iτ,iz,ir,iρ,iA,:] = reverse(aer_rrs_ABO2[:,2]).* (1e7./λ.^2) #mW/m^2/s/nm
                        Qrrs[iτ,iz,ir,iρ,iA,:] = reverse(aer_rrs_ABO2[:,3]).* (1e7./λ.^2) #mW/m^2/s/nm
                        Urrs[iτ,iz,ir,iρ,iA,:] = reverse(aer_rrs_ABO2[:,4]).* (1e7./λ.^2) #mW/m^2/s/nm
                        ieIrrs[iτ,iz,ir,iρ,iA,:] = reverse(aer_rrs_ABO2[:,5]).* (1e7./λ.^2) #mW/m^2/s/nm
                        ieQrrs[iτ,iz,ir,iρ,iA,:] = reverse(aer_rrs_ABO2[:,6]).* (1e7./λ.^2) #mW/m^2/s/nm
                        ieUrrs[iτ,iz,ir,iρ,iA,:] = reverse(aer_rrs_ABO2[:,7]).* (1e7./λ.^2) #mW/m^2/s/nm
                        #F₀ = reverse(aer_rrs_ABO2[:,8]).* (1e7./λ.^2) #mW/m^2/s/nm
                        
                        Inors[iτ,iz,ir,iρ,iA,:] = reverse(aer_nors_ABO2[:,2]).* (1e7./λ.^2) #mW/m^2/s/nm
                        Qnors[iτ,iz,ir,iρ,iA,:] = reverse(aer_nors_ABO2[:,3]).* (1e7./λ.^2) #mW/m^2/s/nm
                        Unors[iτ,iz,ir,iρ,iA,:] = reverse(aer_nors_ABO2[:,4]).* (1e7./λ.^2) #mW/m^2/s/nm

                        #@show rayl_nors_ABO2
                            
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
ΔI = Irrs+ieIrrs-Inors
ΔQ = Qrrs+ieQrrs-Qnors
ΔU = Urrs+ieUrrs-Unors

eΔI = Irrs-Inors
eΔQ = Qrrs-Qnors
eΔU = Urrs-Unors

ieI = ieIrrs
ieQ = ieQrrs
ieU = ieUrrs

# O2 B-band
fname0 = "/home/sanghavi/data/RamanSIFgrid/aerosol_test/aer_AOD0p1_z2_r0p1_sza0_vza0p0_vaz0p0_alb0p0_psurf1000hpa_nors_BBO2.dat"
aer_nors_BBO2 = unique(readdlm(fname0), dims=1)
    
ν_BB = aer_nors_BBO2[:,1]
λ_BB = 1e7./reverse(ν_BB)
F₀_BB = reverse(aer_nors_ABO2[:,5]).* (1e7./λ_BB.^2) #mW/m^2/s/nm

Irrs_BB = zeros(length(τ_ref),length(aer_z₀),length(aer_r₀),length(ρ),length(sza),length(λ_BB))
Qrrs_BB = zeros(length(τ_ref),length(aer_z₀),length(aer_r₀),length(ρ),length(sza),length(λ_BB))
Urrs_BB = zeros(length(τ_ref),length(aer_z₀),length(aer_r₀),length(ρ),length(sza),length(λ_BB))
ieIrrs_BB = zeros(length(τ_ref),length(aer_z₀),length(aer_r₀),length(ρ),length(sza),length(λ_BB))
ieQrrs_BB = zeros(length(τ_ref),length(aer_z₀),length(aer_r₀),length(ρ),length(sza),length(λ_BB))
ieUrrs_BB = zeros(length(τ_ref),length(aer_z₀),length(aer_r₀),length(ρ),length(sza),length(λ_BB))
Inors_BB = zeros(length(τ_ref),length(aer_z₀),length(aer_r₀),length(ρ),length(sza),length(λ_BB))
Qnors_BB = zeros(length(τ_ref),length(aer_z₀),length(aer_r₀),length(ρ),length(sza),length(λ_BB))
Unors_BB = zeros(length(τ_ref),length(aer_z₀),length(aer_r₀),length(ρ),length(sza),length(λ_BB))

for iτ = 1:length(τ_ref)
    for iz = 1:length(aer_z₀)
        for ir = 1:length(aer_r₀)
            for isurf = 1:1 # 3:3 # 2:2 # 1:1 # 
                for iρ = 1:length(ρ) #21 #1:15 #3 #1:15
                    for iA = 1:length(sza) #14

                        vza_str = "0p0"
                        albedo = ρ_str[iρ]
                        sza_str = string(sza[iA])
                        vaz_str = "0p0"
                            
                        # No SIF
                        fname0 = "/home/sanghavi/data/RamanSIFgrid/aerosol_test/aer_AOD"*AOD_str[iτ]*"_z"*z_str[iz]*"_r"*r_str[ir]*"_sza"*sza_str*"_vza"*vza_str*"_vaz"*vaz_str*"_alb"*albedo*"_psurf"*string(psurf[isurf])*"hpa_nors_BBO2.dat"
                        fname1 = "/home/sanghavi/data/RamanSIFgrid/aerosol_test/aer_AOD"*AOD_str[iτ]*"_z"*z_str[iz]*"_r"*r_str[ir]*"_sza"*sza_str*"_vza"*vza_str*"_vaz"*vaz_str*"_alb"*albedo*"_psurf"*string(psurf[isurf])*"hpa_rrs_BBO2.dat"
                            
                        aer_rrs_BBO2  = unique(readdlm(fname1), dims=1)
                        aer_nors_BBO2 = unique(readdlm(fname0), dims=1)
                            
                        #ν = aer_rrs_ABO2[:,1]
                        #λ = 1e7./reverse(ν)
                        Irrs_BB[iτ,iz,ir,iρ,iA,:] = reverse(aer_rrs_BBO2[:,2]).* (1e7./λ_BB.^2) #mW/m^2/s/nm
                        Qrrs_BB[iτ,iz,ir,iρ,iA,:] = reverse(aer_rrs_BBO2[:,3]).* (1e7./λ_BB.^2) #mW/m^2/s/nm
                        Urrs_BB[iτ,iz,ir,iρ,iA,:] = reverse(aer_rrs_BBO2[:,4]).* (1e7./λ_BB.^2) #mW/m^2/s/nm
                        ieIrrs_BB[iτ,iz,ir,iρ,iA,:] = reverse(aer_rrs_BBO2[:,5]).* (1e7./λ_BB.^2) #mW/m^2/s/nm
                        ieQrrs_BB[iτ,iz,ir,iρ,iA,:] = reverse(aer_rrs_BBO2[:,6]).* (1e7./λ_BB.^2) #mW/m^2/s/nm
                        ieUrrs_BB[iτ,iz,ir,iρ,iA,:] = reverse(aer_rrs_BBO2[:,7]).* (1e7./λ_BB.^2) #mW/m^2/s/nm
                        #F₀ = reverse(aer_rrs_ABO2[:,8]).* (1e7./λ.^2) #mW/m^2/s/nm
                        
                        Inors_BB[iτ,iz,ir,iρ,iA,:] = reverse(aer_nors_BBO2[:,2]).* (1e7./λ_BB.^2) #mW/m^2/s/nm
                        Qnors_BB[iτ,iz,ir,iρ,iA,:] = reverse(aer_nors_BBO2[:,3]).* (1e7./λ_BB.^2) #mW/m^2/s/nm
                        Unors_BB[iτ,iz,ir,iρ,iA,:] = reverse(aer_nors_BBO2[:,4]).* (1e7./λ_BB.^2) #mW/m^2/s/nm

                        #@show rayl_nors_ABO2
                            
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
ΔI_BB = Irrs_BB+ieIrrs_BB-Inors_BB
ΔQ_BB = Qrrs_BB+ieQrrs_BB-Qnors_BB
ΔU_BB = Urrs_BB+ieUrrs_BB-Unors_BB

eΔI_BB = Irrs_BB-Inors_BB
eΔQ_BB = Qrrs_BB-Qnors_BB
eΔU_BB = Urrs_BB-Unors_BB

ieI_BB = ieIrrs_BB
ieQ_BB = ieQrrs_BB
ieU_BB = ieUrrs_BB

# Rayleigh case
#fname0 = "/home/sanghavi/data/RamanSIFgrid/aerosol_test/aer_AOD0p0_sza0_vza0p0_vaz0p0_alb0p0_psurf1000hpa_nors_ABO2.dat"
#rayl_nors_ABO2 = unique(readdlm(fname0), dims=1)
    
#ν = rayl_nors_ABO2[:,1]
#λ = 1e7./reverse(ν)
I0rrs = zeros(length(ρ),length(sza),length(λ))
Q0rrs = zeros(length(ρ),length(sza),length(λ))
U0rrs = zeros(length(ρ),length(sza),length(λ))
ieI0rrs = zeros(length(ρ),length(sza),length(λ))
ieQ0rrs = zeros(length(ρ),length(sza),length(λ))
ieU0rrs = zeros(length(ρ),length(sza),length(λ))
#F₀ = reverse(aer_nors_ABO2[:,8]).* (1e7./λ.^2) #mW/m^2/s/nm
I0nors = zeros(length(ρ),length(sza),length(λ))
Q0nors = zeros(length(ρ),length(sza),length(λ))
U0nors = zeros(length(ρ),length(sza),length(λ))


for isurf = 1:1 # 3:3 # 2:2 # 1:1 # 
    for iρ = 1:length(ρ) #21 #1:15 #3 #1:15
        for iA = 1:length(sza) #14

            vza_str = "0p0"
            albedo = ρ_str[iρ]
            sza_str = string(sza[iA])
            vaz_str = "0p0"
                
            # No SIF
            fname0 = "/home/sanghavi/data/RamanSIFgrid/aerosol_test/aer_AOD0p0_sza"*sza_str*"_vza"*vza_str*"_vaz"*vaz_str*"_alb"*albedo*"_psurf"*string(psurf[isurf])*"hpa_nors_ABO2.dat"
            fname1 = "/home/sanghavi/data/RamanSIFgrid/aerosol_test/aer_AOD0p0_sza"*sza_str*"_vza"*vza_str*"_vaz"*vaz_str*"_alb"*albedo*"_psurf"*string(psurf[isurf])*"hpa_rrs_ABO2.dat"
                
            aer0_rrs_ABO2  = unique(readdlm(fname1), dims=1)
            aer0_nors_ABO2 = unique(readdlm(fname0), dims=1)
                
            #ν = aer_rrs_ABO2[:,1]
            #λ = 1e7./reverse(ν)
            I0rrs[iρ,iA,:] = reverse(aer0_rrs_ABO2[:,2]).* (1e7./λ.^2) #mW/m^2/s/nm
            Q0rrs[iρ,iA,:] = reverse(aer0_rrs_ABO2[:,3]).* (1e7./λ.^2) #mW/m^2/s/nm
            U0rrs[iρ,iA,:] = reverse(aer0_rrs_ABO2[:,4]).* (1e7./λ.^2) #mW/m^2/s/nm
            ieI0rrs[iρ,iA,:] = reverse(aer0_rrs_ABO2[:,5]).* (1e7./λ.^2) #mW/m^2/s/nm
            ieQ0rrs[iρ,iA,:] = reverse(aer0_rrs_ABO2[:,6]).* (1e7./λ.^2) #mW/m^2/s/nm
            ieU0rrs[iρ,iA,:] = reverse(aer0_rrs_ABO2[:,7]).* (1e7./λ.^2) #mW/m^2/s/nm
            #F₀ = reverse(aer_rrs_ABO2[:,8]).* (1e7./λ.^2) #mW/m^2/s/nm
            
            I0nors[iρ,iA,:] = reverse(aer0_nors_ABO2[:,2]).* (1e7./λ.^2) #mW/m^2/s/nm
            Q0nors[iρ,iA,:] = reverse(aer0_nors_ABO2[:,3]).* (1e7./λ.^2) #mW/m^2/s/nm
            U0nors[iρ,iA,:] = reverse(aer0_nors_ABO2[:,4]).* (1e7./λ.^2) #mW/m^2/s/nm

            #@show rayl_nors_ABO2
                
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

ΔI0 = I0rrs+ieI0rrs-I0nors
ΔQ0 = Q0rrs+ieQ0rrs-Q0nors
ΔU0 = U0rrs+ieU0rrs-U0nors

eΔI0 = I0rrs-I0nors
eΔQ0 = Q0rrs-Q0nors
eΔU0 = U0rrs-U0nors

ieI0 = ieI0rrs
ieQ0 = ieQ0rrs
ieU0 = ieU0rrs

# Convolution to OCO resolution
dpath = "/net/squid/data3/data/FluoData1/group/oco2/"
L1File   = dpath*"L1bSc/oco2_L1bScND_26782a_190715_B10003r_200429202458.h5"
metFile  = dpath*"L2Met/oco2_L2MetND_26782a_190715_B10003r_200429202458.h5"
dictFile = "/home/cfranken/code/gitHub/InstrumentOperator.jl/json/oco2.yaml"
# Load L1 file (could just use filenames here as well)
oco = InstrumentOperator.load_L1(dictFile,L1File, metFile);
#oco = InstrumentOperator.load_L1(dictFile,L1File, metFile);
# Pick some bands as tuple (or just one)
bands=(1,); #bands = (1,2,3);
# Indices within that band:
#indices = (92:885,114:845,50:916);
indices = (2:1015,); #indices = (114:845,);
# Geo Index (footprint,sounding):
GeoInd = [5,5000];
# Get data for that sounding:
oco_sounding = InstrumentOperator.getMeasurement(oco, bands, indices, GeoInd);

ocoMM=oco_sounding.mueller
iBand=1
res = 0.001e-3;
off = 0.5e-3
oco_wl = oco_sounding.ils[iBand].ν_out[1]-off:res:oco_sounding.ils[iBand].ν_out[end]+off;
oco_λ  = oco_sounding.SpectralGrid
oco_ΔI0 = zeros(length(ρ),length(sza),length(oco_λ))
oco_ΔI = zeros(length(τ_ref),length(aer_z₀),length(aer_r₀),length(ρ),length(sza),length(oco_λ))
oco_eΔI0 = zeros(length(ρ),length(sza),length(oco_λ))
oco_eΔI = zeros(length(τ_ref),length(aer_z₀),length(aer_r₀),length(ρ),length(sza),length(oco_λ))
oco_ieI0 = zeros(length(ρ),length(sza),length(oco_λ))
oco_ieI = zeros(length(τ_ref),length(aer_z₀),length(aer_r₀),length(ρ),length(sza),length(oco_λ))

for isurf = 1:1 # 3:3 # 2:2 # 1:1 # 
    for iρ = 1:length(ρ) #21 #1:15 #3 #1:15
        for iA = 1:length(sza) #14
            oco_ΔI0_tmp1 = (ocoMM[1]*ΔI0[iρ, iA, :] + 
                    ocoMM[2]*ΔQ0[iρ, iA, :] +
                    ocoMM[3]*ΔU0[iρ, iA, :]).* 
                    (1e7./λ.^2) #mW/m^2/s/nm           
            interp_ΔI0 = LinearInterpolation(λ*1e-3, oco_ΔI0_tmp1); 
            tmp1_0 = interp_ΔI0(oco_wl);
            oco_ΔI0[iρ, iA, :] = InstrumentOperator.conv_spectra(oco_sounding.ils[iBand], oco_wl, tmp1_0)
            
            oco_eΔI0_tmp1 = (ocoMM[1]*eΔI0[iρ, iA, :] + 
                    ocoMM[2]*eΔQ0[iρ, iA, :] +
                    ocoMM[3]*eΔU0[iρ, iA, :]).* 
                    (1e7./λ.^2) #mW/m^2/s/nm           
            interp_eΔI0 = LinearInterpolation(λ*1e-3, oco_eΔI0_tmp1); 
            tmp1_0 = interp_eΔI0(oco_wl);
            oco_eΔI0[iρ, iA, :] = InstrumentOperator.conv_spectra(oco_sounding.ils[iBand], oco_wl, tmp1_0)
            
            oco_ieI0_tmp1 = (ocoMM[1]*ieI0[iρ, iA, :] + 
                    ocoMM[2]*ieQ0[iρ, iA, :] +
                    ocoMM[3]*ieU0[iρ, iA, :]).* 
                    (1e7./λ.^2) #mW/m^2/s/nm           
            interp_ieI0 = LinearInterpolation(λ*1e-3, oco_ieI0_tmp1); 
            tmp1_0 = interp_ieI0(oco_wl);
            oco_ieI0[iρ, iA, :] = InstrumentOperator.conv_spectra(oco_sounding.ils[iBand], oco_wl, tmp1_0)
            
            for iτ = 1:length(τ_ref)
                for iz = 1:length(aer_z₀)
                    for ir = 1:length(aer_r₀)
                        oco_ΔI_tmp1 = (ocoMM[1]*ΔI[iτ,iz,ir,iρ,iA,:] + 
                            ocoMM[2]*ΔQ[iτ,iz,ir,iρ,iA,:] +
                            ocoMM[3]*ΔU[iτ,iz,ir,iρ,iA,:]).* 
                            (1e7./λ.^2) #mW/m^2/s/nm
                        interp_ΔI  = LinearInterpolation(λ*1e-3, oco_ΔI_tmp1); 
                        tmp1   = interp_ΔI(oco_wl);
                        oco_ΔI[iτ,iz,ir,iρ,iA,:] = InstrumentOperator.conv_spectra(oco_sounding.ils[iBand], oco_wl, tmp1)            
                        
                        oco_eΔI_tmp1 = (ocoMM[1]*eΔI[iτ,iz,ir,iρ,iA,:] + 
                            ocoMM[2]*eΔQ[iτ,iz,ir,iρ,iA,:] +
                            ocoMM[3]*eΔU[iτ,iz,ir,iρ,iA,:]).* 
                            (1e7./λ.^2) #mW/m^2/s/nm
                        interp_eΔI  = LinearInterpolation(λ*1e-3, oco_eΔI_tmp1); 
                        tmp1   = interp_eΔI(oco_wl);
                        oco_eΔI[iτ,iz,ir,iρ,iA,:] = InstrumentOperator.conv_spectra(oco_sounding.ils[iBand], oco_wl, tmp1)            
                        
                        oco_ieI_tmp1 = (ocoMM[1]*ieI[iτ,iz,ir,iρ,iA,:] + 
                            ocoMM[2]*ieQ[iτ,iz,ir,iρ,iA,:] +
                            ocoMM[3]*ieU[iτ,iz,ir,iρ,iA,:]).* 
                            (1e7./λ.^2) #mW/m^2/s/nm
                        interp_ieI  = LinearInterpolation(λ*1e-3, oco_ieI_tmp1); 
                        tmp1   = interp_ieI(oco_wl);
                        oco_ieI[iτ,iz,ir,iρ,iA,:] = InstrumentOperator.conv_spectra(oco_sounding.ils[iBand], oco_wl, tmp1)            
                    end
                end
            end
        end
    end
end

fit_window0 = [684., 686.] #[680.0, 682.0] #[673.0, 675.0] #[680.0, 683.0] 
fit_window1 = [758.0, 759.2] 
fit_window2 = [770.0, 770.25] 

ind1 = findall(x->x>fit_window1[1]*1e-3 && x<fit_window1[2]*1e-3, oco_λ);         
ind2 = findall(x->x>fit_window2[1]*1e-3 && x<fit_window2[2]*1e-3, oco_λ);         
ind1_hires = findall(x->x>fit_window1[1] && x<fit_window1[2], λ);         
ind2_hires = findall(x->x>fit_window2[1] && x<fit_window2[2], λ);         


using CairoMakie#, GeoMakie
#using GLMakie
using Makie.FileIO
CairoMakie.activate!() # hide
using GridLayoutBase
using Colors 

f = Figure(backgroundcolor = RGBf(0.98, 0.98, 0.98),
    size = (1200, 900))#size = (800, 600))#size = (900, 600))
g1 = GridLayout(f[1:8, 1:6], nrows = 4, ncols = 6)
g2 = GridLayout(f[1:8, 7:9], nrows = 1, ncols = 3)
Label(f[1, 1:9, Top()], "Effect of surface brightness",
        fontsize = 15,
        font = :bold,
        padding = (0, 0, 30, 0),
        halign = :center)
  
ax5 = Axis(g2[1,1:2], title = L"$I_\mathrm{i} = I_\mathrm{RRS}$", 
        xlabel = L"$\lambda\,[\mathrm{nm}]$", 
        ylabel = "")  
ax6 = Axis(g2[1,3], title = L"$I_\mathrm{i}\,\,\mathrm{Window 1}$", 
        xlabel = L"$\lambda\,[\mathrm{nm}]$", 
        ylabel = "")               

# plots
ir=1 # r₀=0.1
iz=1 # z₀=2 km 
iA=1 # sza=0
# Effect of aerosol AOD for varying surface brightnesses
iρ=1 # ρ=0
#palette = cgrad(:dense)   # Create a color gradient from "inferno"
#clist0 = palette[LinRange(0, 1, length(ρ)+2)]  # Sample 10 colors
#palette = cgrad(:algae)   # Create a color gradient from "inferno"
#clist1 = palette[LinRange(0, 1, length(ρ)+2)]
#palette = cgrad(:amp)   # Create a color gradient from "inferno"
#clist2 = palette[LinRange(0, 1, length(ρ)+2)]
clist0 = reshape( range(RGBA(0.,0.,1.,0.9), stop=RGBA(0.,0.,1.,0.3),length=length(ρ)), 1, length(ρ) );
clist1 = reshape( range(RGBA(0.,1.,0.,0.9), stop=RGBA(0.,1.,0.,0.3),length=length(ρ)), 1, length(ρ) );
clist2 = reshape( range(RGBA(1.,0.,0.,0.9), stop=RGBA(1.,0.,0.,0.3),length=length(ρ)), 1, length(ρ) );

for iρ=1:length(ρ)
    rev_iρ = (4-iρ)+1
    #i1 = (iρ-1)*2
    #i0 = i1-1
    if(iρ>1)
        ax1 = Axis(g1[rev_iρ,1:2], title = L"$\Delta I = I_\mathrm{Cab}+I_\mathrm{RRS}-I_\mathrm{Rayl}$", 
        xlabel = "", 
        ylabel = L"$I_\mathrm{TOA}\,[\mathrm{mW}/\mathrm{nm}/\mathrm{sr}/\mathrm{m}^2]$")        
    else
        ax1 = Axis(g1[rev_iρ,1:2], title = L"$\Delta I = I_\mathrm{Cab}+I_\mathrm{RRS}-I_\mathrm{Rayl}$", 
        xlabel = L"$\lambda\,[\mathrm{nm}]$", 
        ylabel = L"$I_\mathrm{TOA}\,[\mathrm{mW}/\mathrm{nm}/\mathrm{sr}/\mathrm{m}^2]$")        
    end
    iτ = 2 # AOD=0.5
    lines!(ax1, λ, ΔI[iτ,iz,ir,iρ,iA,:], color=clist2[iρ], linewidth=0.75,  label=L"$\mathrm{AOD}=$"*string(τ_ref[iτ])*L"$\,\,ρ=$"*string(ρ[iρ]))
    iτ = 1
    lines!(ax1, λ, ΔI[iτ,iz,ir,iρ,iA,:], color=clist1[iρ], linewidth=0.75,  label=L"$\mathrm{AOD}=$"*string(τ_ref[iτ])*L"$\,\,ρ=$"*string(ρ[iρ]))
    lines!(ax1, λ, ΔI0[iρ,iA,:], color=clist0[iρ], linewidth=0.75,  label=L"$\mathrm{AOD}=$"*string(τ_ref[iτ])*L"$\,\,ρ=$"*string(ρ[iρ]))
    
    if(iρ>1)
        iτ = 2 # AOD=0.5
        lines!(ax1, λ, ΔI[iτ,iz,ir,1,iA,:], color=clist2[iρ], linewidth=0.75,  label=L"$\mathrm{AOD}=$"*string(τ_ref[iτ])*L"$\,\,ρ=$"*string(ρ[iρ]))
        iτ = 1
        lines!(ax1, λ, ΔI[iτ,iz,ir,1,iA,:], color=clist1[iρ], linewidth=0.75,  label=L"$\mathrm{AOD}=$"*string(τ_ref[iτ])*L"$\,\,ρ=$"*string(ρ[iρ]))
        lines!(ax1, λ, ΔI0[1,iA,:], color=clist0[iρ], linewidth=0.75,  label=L"$\mathrm{AOD}=$"*string(τ_ref[iτ])*L"$\,\,ρ=$"*string(ρ[iρ]))
    end
    ylims!(ax1, -6, 6)

    if(iρ>1)
        ax2 = Axis(g1[rev_iρ,3], title = L"$\Delta I\,\,\mathrm{Window 1}$", 
        xlabel = "", 
        ylabel = "")        
    else
        ax2 = Axis(g1[rev_iρ,3], title = L"$\Delta I\,\,\mathrm{Window 1}$", 
        xlabel = L"$\lambda\,[\mathrm{nm}]$", 
        ylabel = "")        
    end
    
    iτ = 2 # AOD=0.5
    lines!(ax2, λ[ind1_hires], ΔI[iτ,iz,ir,iρ,iA,ind1_hires], color=clist2[iρ], linewidth=0.75,  label=L"$\mathrm{AOD}=$"*string(τ_ref[iτ])*L"$\,\,ρ=$"*string(ρ[iρ]))
    iτ = 1
    lines!(ax2, λ[ind1_hires], ΔI[iτ,iz,ir,iρ,iA,ind1_hires], color=clist1[iρ], linewidth=0.75,  label=L"$\mathrm{AOD}=$"*string(τ_ref[iτ])*L"$\,\,ρ=$"*string(ρ[iρ]))
    lines!(ax2, λ[ind1_hires], ΔI0[iρ,iA,ind1_hires], color=clist0[iρ], linewidth=0.75,  label=L"$\mathrm{AOD}=$"*string(τ_ref[iτ])*L"$\,\,ρ=$"*string(ρ[iρ]))
        
    if(iρ>1)
        iτ = 2 # AOD=0.5
        lines!(ax2, λ[ind1_hires], ΔI[iτ,iz,ir,1,iA,ind1_hires], color=RGBA(0.8, 0.8, 0.8, 0.75), linewidth=0.75,  label=L"$\mathrm{AOD}=$"*string(τ_ref[iτ])*L"$\,\,ρ=$"*string(ρ[1]))
        iτ = 1
        lines!(ax2, λ[ind1_hires], ΔI[iτ,iz,ir,1,iA,ind1_hires], color=RGBA(0.6, 0.6, 0.6, 0.75), linewidth=0.75,  label=L"$\mathrm{AOD}=$"*string(τ_ref[iτ])*L"$\,\,ρ=$"*string(ρ[1]))
        lines!(ax2, λ[ind1_hires], ΔI0[1,iA,ind1_hires], color=RGBA(0.4, 0.4, 0.4, 0.75), linewidth=0.75,  label=L"$\mathrm{AOD}=$"*string(τ_ref[iτ])*L"$\,\,ρ=$"*string(ρ[1]))    
    end
    ylims!(ax2, -0.15, 0.3)
    
    if(iρ>1)
        ax3 = Axis(g1[rev_iρ, 4:5], title = L"$\Delta I_\mathrm{e} = I_\mathrm{Cab}-I_\mathrm{Rayl}$", 
        xlabel = "", 
        ylabel = "")
    else    
        ax3 = Axis(g1[rev_iρ,4:5], title = L"$\Delta I_\mathrm{e} = I_\mathrm{Cab}-I_\mathrm{Rayl}$", 
        xlabel = L"$\lambda\,[\mathrm{nm}]$", 
        ylabel = "")
    end

    iτ = 2 # AOD=0.5
    lines!(ax3, λ, eΔI[iτ,iz,ir,iρ,iA,:], color=clist2[iρ], linewidth=0.75,  label=L"$\mathrm{AOD}=$"*string(τ_ref[iτ])*L"$\,\,ρ=$"*string(ρ[iρ]))
    iτ = 1
    lines!(ax3, λ, eΔI[iτ,iz,ir,iρ,iA,:], color=clist1[iρ], linewidth=0.75,  label=L"$\mathrm{AOD}=$"*string(τ_ref[iτ])*L"$\,\,ρ=$"*string(ρ[iρ]))
    lines!(ax3, λ, eΔI0[iρ,iA,:], color=clist0[iρ], linewidth=0.75,  label=L"$\mathrm{AOD}=$"*string(τ_ref[iτ])*L"$\,\,ρ=$"*string(ρ[iρ]))
        
    if(iρ>1)
        iτ = 2 # AOD=0.5
        lines!(ax3, λ, eΔI[iτ,iz,ir,1,iA,:], color=RGBA(0.8, 0.8, 0.8, 0.75), linewidth=0.75,  label=L"$\mathrm{AOD}=$"*string(τ_ref[iτ])*L"$\,\,ρ=$"*string(ρ[1]))
        iτ = 1
        lines!(ax3, λ, eΔI[iτ,iz,ir,1,iA,:], color=RGBA(0.6, 0.6, 0.6, 0.75), linewidth=0.75,  label=L"$\mathrm{AOD}=$"*string(τ_ref[iτ])*L"$\,\,ρ=$"*string(ρ[1]))
        lines!(ax3, λ, eΔI0[1,iA,:], color=RGBA(0.4, 0.4, 0.4, 0.75), linewidth=0.75,  label=L"$\mathrm{AOD}=$"*string(τ_ref[iτ])*L"$\,\,ρ=$"*string(ρ[1]))    
    end
    ylims!(ax3, -6, 6) 

    if(iρ>1)
        ax4 = Axis(g1[rev_iρ, 6], title = L"$\Delta I_\mathrm{e}\,\,\mathrm{Window 1}$", 
        xlabel = "", 
        ylabel = "")
    else    
        ax4 = Axis(g1[rev_iρ, 6], title = L"$\Delta I_\mathrm{e}\,\,\mathrm{Window 1}$", 
        xlabel = L"$\lambda\,[\mathrm{nm}]$", 
        ylabel = "")
    end

    iτ = 2 # AOD=0.5
    lines!(ax4, λ[ind1_hires], eΔI[iτ,iz,ir,iρ,iA,ind1_hires], color=clist2[iρ], linewidth=0.75,  label=L"$\mathrm{AOD}=$"*string(τ_ref[iτ])*L"$\,\,ρ=$"*string(ρ[iρ]))
    iτ = 1
    lines!(ax4, λ[ind1_hires], eΔI[iτ,iz,ir,iρ,iA,ind1_hires], color=clist1[iρ], linewidth=0.75,  label=L"$\mathrm{AOD}=$"*string(τ_ref[iτ])*L"$\,\,ρ=$"*string(ρ[iρ]))
    lines!(ax4, λ[ind1_hires], eΔI0[iρ,iA,ind1_hires], color=clist0[iρ], linewidth=0.75,  label=L"$\mathrm{AOD}=$"*string(τ_ref[iτ])*L"$\,\,ρ=$"*string(ρ[iρ]))
        
    if(iρ>1)
        iτ = 2 # AOD=0.5
        lines!(ax4, λ[ind1_hires], eΔI[iτ,iz,ir,1,iA,ind1_hires], color=RGBA(0.8, 0.8, 0.8, 0.75), linewidth=0.75,  label=L"$\mathrm{AOD}=$"*string(τ_ref[iτ])*L"$\,\,ρ=$"*string(ρ[1]))
        iτ = 1
        lines!(ax4, λ[ind1_hires], eΔI[iτ,iz,ir,1,iA,ind1_hires], color=RGBA(0.6, 0.6, 0.6, 0.75), linewidth=0.75,  label=L"$\mathrm{AOD}=$"*string(τ_ref[iτ])*L"$\,\,ρ=$"*string(ρ[1]))
        lines!(ax4, λ[ind1_hires], eΔI0[1,iA,ind1_hires], color=RGBA(0.4, 0.4, 0.4, 0.75), linewidth=0.75,  label=L"$\mathrm{AOD}=$"*string(τ_ref[iτ])*L"$\,\,ρ=$"*string(ρ[1]))    
    end
    ylims!(ax4, -0.7,0) 

    iτ = 2 # AOD=0.5
    lines!(ax5, λ, ieI[iτ,iz,ir,iρ,iA,:], color=clist2[iρ], linewidth=0.75,  label=L"$\mathrm{AOD}=$"*string(τ_ref[iτ])*L"$\,\,ρ=$"*string(ρ[iρ]))
    lines!(ax6, λ[ind1_hires], ieI[iτ,iz,ir,iρ,iA,ind1_hires], color=clist2[iρ], linewidth=0.75,  label=L"$\mathrm{AOD}=$"*string(τ_ref[iτ])*L"$\,\,ρ=$"*string(ρ[iρ]))
    
    iτ = 1
    lines!(ax5, λ, ieI[iτ,iz,ir,iρ,iA,:], color=clist1[iρ], linewidth=0.75,  label=L"$\mathrm{AOD}=$"*string(τ_ref[iτ])*L"$\,\,ρ=$"*string(ρ[iρ]))
    lines!(ax6, λ[ind1_hires], ieI[iτ,iz,ir,iρ,iA,ind1_hires], color=clist1[iρ], linewidth=0.75,  label=L"$\mathrm{AOD}=$"*string(τ_ref[iτ])*L"$\,\,ρ=$"*string(ρ[iρ]))
    
    lines!(ax5, λ, ieI0[iρ,iA,:], color=clist0[iρ], linewidth=0.75,  label=L"$\mathrm{AOD}=$"*string(τ_ref[iτ])*L"$\,\,ρ=$"*string(ρ[iρ]))
    lines!(ax6, λ[ind1_hires], ieI0[iρ,iA,ind1_hires], color=clist0[iρ], linewidth=0.75,  label=L"$\mathrm{AOD}=$"*string(τ_ref[iτ])*L"$\,\,ρ=$"*string(ρ[iρ]))    
end
# Manually create a legend entry with a larger marker
legend_line21 = LineElement(color=clist2[1]) #lines!(ax1, λ[1:2], ΔI[1,1,1,1,1,1:2],color=clist2[2],visible=false)
legend_line22 = LineElement(color=clist2[2]) #lines!(ax1, λ[1:2], ΔI[1,1,1,1,1,1:2],color=clist2[3],visible=false)
legend_line23 = LineElement(color=clist2[3]) #lines!(ax1, λ[1:2], ΔI[1,1,1,1,1,1:2],color=clist2[4],visible=false)
legend_line24 = LineElement(color=clist2[4]) #lines!(ax1, λ[1:2], ΔI[1,1,1,1,1,1:2],color=clist2[5],visible=false)

legend_line11 = LineElement(color=clist1[1]) #lines!(ax1, λ[1:2], ΔI[1,1,1,1,1,1:2],color=clist1[2],visible=false)
legend_line12 = LineElement(color=clist1[2]) #lines!(ax1, λ[1:2], ΔI[1,1,1,1,1,1:2],color=clist1[3],visible=false)
legend_line13 = LineElement(color=clist1[3]) #lines!(ax1, λ[1:2], ΔI[1,1,1,1,1,1:2],color=clist1[4],visible=false)
legend_line14 = LineElement(color=clist1[4]) #lines!(ax1, λ[1:2], ΔI[1,1,1,1,1,1:2],color=clist1[5],visible=false)

legend_line01 = LineElement(color=clist0[1]) #lines!(ax1, λ[1:2], ΔI0[1,1,1:2],color=clist0[2],visible=false)
legend_line02 = LineElement(color=clist0[2]) #lines!(ax1, λ[1:2], ΔI0[1,1,1:2],color=clist0[3],visible=false)
legend_line03 = LineElement(color=clist0[3]) #lines!(ax1, λ[1:2], ΔI0[1,1,1:2],color=clist0[4],visible=false)
legend_line04 = LineElement(color=clist0[4]) #lines!(ax1, λ[1:2], ΔI0[1,1,1:2],color=clist0[5],visible=false)

#=legend_labels = [
    L"\mathrm{AOD}=0\,\,\rho=" * LaTeXString(string(ρ[i])) for i in 1:4
] ∪ [
    L"\mathrm{AOD}=" * LaTeXString(string(τ_ref[1])) * L"\,\,\rho=" * LaTeXString(string(ρ[i])) for i in 1:4
] ∪ [
    L"\mathrm{AOD}=" * LaTeXString(string(τ_ref[2])) * L"\,\,\rho=" * LaTeXString(string(ρ[i])) for i in 1:4
]=#
legend_labels = [
    latexstring(L"\mathrm{AOD}=0\,\,\rho=" * latexstring(string(ρ[1]))),
    latexstring(L"\mathrm{AOD}=" * latexstring(string(τ_ref[1])) * L"\,\,\rho=" * latexstring(string(ρ[1]))),
    latexstring(L"\mathrm{AOD}=" * latexstring(string(τ_ref[2])) * L"\,\,\rho=" * latexstring(string(ρ[1]))),
    latexstring(L"\mathrm{AOD}=0\,\,\rho=" * latexstring(string(ρ[2]))),
    latexstring(L"\mathrm{AOD}=" * latexstring(string(τ_ref[1])) * L"\,\,\rho=" * latexstring(string(ρ[2]))),
    latexstring(L"\mathrm{AOD}=" * latexstring(string(τ_ref[2])) * L"\,\,\rho=" * latexstring(string(ρ[2]))),
    latexstring(L"\mathrm{AOD}=0\,\,\rho=" * latexstring(string(ρ[3]))),
    latexstring(L"\mathrm{AOD}=" * latexstring(string(τ_ref[1])) * L"\,\,\rho=" * latexstring(string(ρ[3]))),
    latexstring(L"\mathrm{AOD}=" * latexstring(string(τ_ref[2])) * L"\,\,\rho=" * latexstring(string(ρ[3]))),
    latexstring(L"\mathrm{AOD}=0\,\,\rho=" * latexstring(string(ρ[4]))),
    latexstring(L"\mathrm{AOD}=" * latexstring(string(τ_ref[1])) * L"\,\,\rho=" * latexstring(string(ρ[4]))),
    latexstring(L"\mathrm{AOD}=" * latexstring(string(τ_ref[2])) * L"\,\,\rho=" * latexstring(string(ρ[4])))
]

legend_elements = [
    legend_line01, legend_line11, legend_line21, 
    legend_line02, legend_line12, legend_line22, 
    legend_line03, legend_line13, legend_line23, 
    legend_line04, legend_line14, legend_line24]
legend = Legend(f, legend_elements, legend_labels, 
           orientation = :vertical,  # Ensures elements are stacked properly
           tellwidth = false,  # Avoids stretching legend width
           nbanks = 3,  # Specifies 4 columns (i.e., 3 rows automatically)
           halign = :center,
           valign = :center
       )

f[9, :] = legend
save("newRRScorr_wrt_AOD_z"*z_str[iz]*"_r"*r_str[ir]*"_sza"*string(sza[iA])*".png", f)




f = Figure(backgroundcolor = RGBf(0.98, 0.98, 0.98),
    size = (1200, 900))#size = (800, 600))#size = (900, 600))

Label(f[1, 1:9, Top()], "Effect of surface brightness",
        fontsize = 15,
        font = :bold,
        padding = (0, 0, 30, 0),
        halign = :center)

ax1 = Axis(f[1:4,1:2], title = L"$\Delta I = I_\mathrm{Cab}+I_\mathrm{RRS}-I_\mathrm{Rayl}$", 
        xlabel = L"$\lambda\,[\mathrm{nm}]$", 
        ylabel = L"$I_\mathrm{TOA}\,[\mathrm{mW}/\mathrm{nm}/\mathrm{sr}/\mathrm{m}^2]$")        
ax2 = Axis(f[1:4,3], title = L"$\Delta I\,\,\mathrm{Window 1}$", 
        xlabel = L"$\lambda\,[\mathrm{nm}]$", 
        ylabel = "")   
ax3 = Axis(f[1:4,4:5], title = L"$\Delta I_\mathrm{e} = I_\mathrm{Cab}-I_\mathrm{Rayl}$", 
        xlabel = L"$\lambda\,[\mathrm{nm}]$", 
        ylabel = "")        
ax4 = Axis(f[1:4,6], title = L"$\Delta I_\mathrm{e}\,\,\mathrm{Window 1}$", 
        xlabel = L"$\lambda\,[\mathrm{nm}]$", 
        ylabel = "")    
ax5 = Axis(f[1:4,7:8], title = L"$I_\mathrm{i} = I_\mathrm{RRS}$", 
        xlabel = L"$\lambda\,[\mathrm{nm}]$", 
        ylabel = "")  
ax6 = Axis(f[1:4,9], title = L"$I_\mathrm{i}\,\,\mathrm{Window 1}$", 
        xlabel = L"$\lambda\,[\mathrm{nm}]$", 
        ylabel = "")               

# plots
ir=1 # r₀=0.1
iz=1 # z₀=2 km 
iA=1 # sza=0
# Effect of aerosol AOD for varying surface brightnesses
iρ=1 # ρ=0
#palette = cgrad(:dense)   # Create a color gradient from "inferno"
#clist0 = palette[LinRange(0, 1, length(ρ)+2)]  # Sample 10 colors
#palette = cgrad(:algae)   # Create a color gradient from "inferno"
#clist1 = palette[LinRange(0, 1, length(ρ)+2)]
#palette = cgrad(:amp)   # Create a color gradient from "inferno"
#clist2 = palette[LinRange(0, 1, length(ρ)+2)]
clist0 = reshape( range(RGBA(0.,0.,1.,0.9), stop=RGBA(0.,0.,1.,0.3),length=length(ρ)), 1, length(ρ) );
clist1 = reshape( range(RGBA(0.,1.,0.,0.9), stop=RGBA(0.,1.,0.,0.3),length=length(ρ)), 1, length(ρ) );
clist2 = reshape( range(RGBA(1.,0.,0.,0.9), stop=RGBA(1.,0.,0.,0.3),length=length(ρ)), 1, length(ρ) );

for iρ=1:length(ρ)
    iτ = 2 # AOD=0.5
    lines!(ax1, λ, ΔI[iτ,iz,ir,iρ,iA,:], color=clist2[iρ], linewidth=0.75,  label=L"$\mathrm{AOD}=$"*string(τ_ref[iτ])*L"$\,\,ρ=$"*string(ρ[iρ]))
    lines!(ax2, λ[ind1_hires], ΔI[iτ,iz,ir,iρ,iA,ind1_hires], color=clist2[iρ], linewidth=0.75,  label=L"$\mathrm{AOD}=$"*string(τ_ref[iτ])*L"$\,\,ρ=$"*string(ρ[iρ]))
    lines!(ax3, λ, eΔI[iτ,iz,ir,iρ,iA,:], color=clist2[iρ], linewidth=0.75,  label=L"$\mathrm{AOD}=$"*string(τ_ref[iτ])*L"$\,\,ρ=$"*string(ρ[iρ]))
    lines!(ax4, λ[ind1_hires], eΔI[iτ,iz,ir,iρ,iA,ind1_hires], color=clist2[iρ], linewidth=0.75,  label=L"$\mathrm{AOD}=$"*string(τ_ref[iτ])*L"$\,\,ρ=$"*string(ρ[iρ]))
    lines!(ax5, λ, ieI[iτ,iz,ir,iρ,iA,:], color=clist2[iρ], linewidth=0.75,  label=L"$\mathrm{AOD}=$"*string(τ_ref[iτ])*L"$\,\,ρ=$"*string(ρ[iρ]))
    lines!(ax6, λ[ind1_hires], ieI[iτ,iz,ir,iρ,iA,ind1_hires], color=clist2[iρ], linewidth=0.75,  label=L"$\mathrm{AOD}=$"*string(τ_ref[iτ])*L"$\,\,ρ=$"*string(ρ[iρ]))
    
    iτ = 1
    lines!(ax1, λ, ΔI[iτ,iz,ir,iρ,iA,:], color=clist1[iρ], linewidth=0.75,  label=L"$\mathrm{AOD}=$"*string(τ_ref[iτ])*L"$\,\,ρ=$"*string(ρ[iρ]))
    lines!(ax2, λ[ind1_hires], ΔI[iτ,iz,ir,iρ,iA,ind1_hires], color=clist1[iρ], linewidth=0.75,  label=L"$\mathrm{AOD}=$"*string(τ_ref[iτ])*L"$\,\,ρ=$"*string(ρ[iρ]))
    lines!(ax3, λ, eΔI[iτ,iz,ir,iρ,iA,:], color=clist1[iρ], linewidth=0.75,  label=L"$\mathrm{AOD}=$"*string(τ_ref[iτ])*L"$\,\,ρ=$"*string(ρ[iρ]))
    lines!(ax4, λ[ind1_hires], eΔI[iτ,iz,ir,iρ,iA,ind1_hires], color=clist1[iρ], linewidth=0.75,  label=L"$\mathrm{AOD}=$"*string(τ_ref[iτ])*L"$\,\,ρ=$"*string(ρ[iρ]))
    lines!(ax5, λ, ieI[iτ,iz,ir,iρ,iA,:], color=clist1[iρ], linewidth=0.75,  label=L"$\mathrm{AOD}=$"*string(τ_ref[iτ])*L"$\,\,ρ=$"*string(ρ[iρ]))
    lines!(ax6, λ[ind1_hires], ieI[iτ,iz,ir,iρ,iA,ind1_hires], color=clist1[iρ], linewidth=0.75,  label=L"$\mathrm{AOD}=$"*string(τ_ref[iτ])*L"$\,\,ρ=$"*string(ρ[iρ]))
    
    lines!(ax1, λ, ΔI0[iρ,iA,:], color=clist0[iρ], linewidth=0.75,  label=L"$\mathrm{AOD}=$"*string(τ_ref[iτ])*L"$\,\,ρ=$"*string(ρ[iρ]))
    lines!(ax2, λ[ind1_hires], ΔI0[iρ,iA,ind1_hires], color=clist0[iρ], linewidth=0.75,  label=L"$\mathrm{AOD}=$"*string(τ_ref[iτ])*L"$\,\,ρ=$"*string(ρ[iρ]))
    lines!(ax3, λ, eΔI0[iρ,iA,:], color=clist0[iρ], linewidth=0.75,  label=L"$\mathrm{AOD}=$"*string(τ_ref[iτ])*L"$\,\,ρ=$"*string(ρ[iρ]))
    lines!(ax4, λ[ind1_hires], eΔI0[iρ,iA,ind1_hires], color=clist0[iρ], linewidth=0.75,  label=L"$\mathrm{AOD}=$"*string(τ_ref[iτ])*L"$\,\,ρ=$"*string(ρ[iρ]))
    lines!(ax5, λ, ieI0[iρ,iA,:], color=clist0[iρ], linewidth=0.75,  label=L"$\mathrm{AOD}=$"*string(τ_ref[iτ])*L"$\,\,ρ=$"*string(ρ[iρ]))
    lines!(ax6, λ[ind1_hires], ieI0[iρ,iA,ind1_hires], color=clist0[iρ], linewidth=0.75,  label=L"$\mathrm{AOD}=$"*string(τ_ref[iτ])*L"$\,\,ρ=$"*string(ρ[iρ]))
    
end
# Manually create a legend entry with a larger marker
legend_line21 = LineElement(color=clist2[1]) #lines!(ax1, λ[1:2], ΔI[1,1,1,1,1,1:2],color=clist2[2],visible=false)
legend_line22 = LineElement(color=clist2[2]) #lines!(ax1, λ[1:2], ΔI[1,1,1,1,1,1:2],color=clist2[3],visible=false)
legend_line23 = LineElement(color=clist2[3]) #lines!(ax1, λ[1:2], ΔI[1,1,1,1,1,1:2],color=clist2[4],visible=false)
legend_line24 = LineElement(color=clist2[4]) #lines!(ax1, λ[1:2], ΔI[1,1,1,1,1,1:2],color=clist2[5],visible=false)

legend_line11 = LineElement(color=clist1[1]) #lines!(ax1, λ[1:2], ΔI[1,1,1,1,1,1:2],color=clist1[2],visible=false)
legend_line12 = LineElement(color=clist1[2]) #lines!(ax1, λ[1:2], ΔI[1,1,1,1,1,1:2],color=clist1[3],visible=false)
legend_line13 = LineElement(color=clist1[3]) #lines!(ax1, λ[1:2], ΔI[1,1,1,1,1,1:2],color=clist1[4],visible=false)
legend_line14 = LineElement(color=clist1[4]) #lines!(ax1, λ[1:2], ΔI[1,1,1,1,1,1:2],color=clist1[5],visible=false)

legend_line01 = LineElement(color=clist0[1]) #lines!(ax1, λ[1:2], ΔI0[1,1,1:2],color=clist0[2],visible=false)
legend_line02 = LineElement(color=clist0[2]) #lines!(ax1, λ[1:2], ΔI0[1,1,1:2],color=clist0[3],visible=false)
legend_line03 = LineElement(color=clist0[3]) #lines!(ax1, λ[1:2], ΔI0[1,1,1:2],color=clist0[4],visible=false)
legend_line04 = LineElement(color=clist0[4]) #lines!(ax1, λ[1:2], ΔI0[1,1,1:2],color=clist0[5],visible=false)

#=legend_labels = [
    L"\mathrm{AOD}=0\,\,\rho=" * LaTeXString(string(ρ[i])) for i in 1:4
] ∪ [
    L"\mathrm{AOD}=" * LaTeXString(string(τ_ref[1])) * L"\,\,\rho=" * LaTeXString(string(ρ[i])) for i in 1:4
] ∪ [
    L"\mathrm{AOD}=" * LaTeXString(string(τ_ref[2])) * L"\,\,\rho=" * LaTeXString(string(ρ[i])) for i in 1:4
]=#
legend_labels = [
    latexstring(L"\mathrm{AOD}=0\,\,\rho=" * latexstring(string(ρ[1]))),
    latexstring(L"\mathrm{AOD}=" * latexstring(string(τ_ref[1])) * L"\,\,\rho=" * latexstring(string(ρ[1]))),
    latexstring(L"\mathrm{AOD}=" * latexstring(string(τ_ref[2])) * L"\,\,\rho=" * latexstring(string(ρ[1]))),
    latexstring(L"\mathrm{AOD}=0\,\,\rho=" * latexstring(string(ρ[2]))),
    latexstring(L"\mathrm{AOD}=" * latexstring(string(τ_ref[1])) * L"\,\,\rho=" * latexstring(string(ρ[2]))),
    latexstring(L"\mathrm{AOD}=" * latexstring(string(τ_ref[2])) * L"\,\,\rho=" * latexstring(string(ρ[2]))),
    latexstring(L"\mathrm{AOD}=0\,\,\rho=" * latexstring(string(ρ[3]))),
    latexstring(L"\mathrm{AOD}=" * latexstring(string(τ_ref[1])) * L"\,\,\rho=" * latexstring(string(ρ[3]))),
    latexstring(L"\mathrm{AOD}=" * latexstring(string(τ_ref[2])) * L"\,\,\rho=" * latexstring(string(ρ[3]))),
    latexstring(L"\mathrm{AOD}=0\,\,\rho=" * latexstring(string(ρ[4]))),
    latexstring(L"\mathrm{AOD}=" * latexstring(string(τ_ref[1])) * L"\,\,\rho=" * latexstring(string(ρ[4]))),
    latexstring(L"\mathrm{AOD}=" * latexstring(string(τ_ref[2])) * L"\,\,\rho=" * latexstring(string(ρ[4])))
]

legend_elements = [
    legend_line01, legend_line11, legend_line21, 
    legend_line02, legend_line12, legend_line22, 
    legend_line03, legend_line13, legend_line23, 
    legend_line04, legend_line14, legend_line24]
legend = Legend(f, legend_elements, legend_labels, 
           orientation = :vertical,  # Ensures elements are stacked properly
           tellwidth = false,  # Avoids stretching legend width
           nbanks = 3,  # Specifies 4 columns (i.e., 3 rows automatically)
           halign = :center,
           valign = :center
       )

f[5:6, :] = legend
save("RRScorr_wrt_AOD_z"*z_str[iz]*"_r"*r_str[ir]*"_sza"*string(sza[iA])*".png", f)

f = Figure(backgroundcolor = RGBf(0.98, 0.98, 0.98),
    size = (1200, 900))#size = (800, 600))#size = (900, 600))

Label(f[1, 1:9, Top()], "Effect of surface brightness",
        fontsize = 15,
        font = :bold,
        padding = (0, 0, 30, 0),
        halign = :center)

ax1 = Axis(f[1:4,1:2], title = L"$\Delta I = I_\mathrm{Cab}+I_\mathrm{RRS}-I_\mathrm{Rayl}$", 
        xlabel = L"$\lambda\,[\mathrm{nm}]$", 
        ylabel = L"$I_\mathrm{TOA}\,[\mathrm{mW}/\mathrm{nm}/\mathrm{sr}/\mathrm{m}^2]$")        
ax2 = Axis(f[1:4,3], title = L"$\Delta I\,\,\mathrm{Window 1}$", 
        xlabel = L"$\lambda\,[\mathrm{nm}]$", 
        ylabel = "")   
ax3 = Axis(f[1:4,4:5], title = L"$\Delta I_\mathrm{e} = I_\mathrm{Cab}-I_\mathrm{Rayl}$", 
        xlabel = L"$\lambda\,[\mathrm{nm}]$", 
        ylabel = "")        
ax4 = Axis(f[1:4,6], title = L"$\Delta I_\mathrm{e}\,\,\mathrm{Window 1}$", 
        xlabel = L"$\lambda\,[\mathrm{nm}]$", 
        ylabel = "")    
ax5 = Axis(f[1:4,7:8], title = L"$I_\mathrm{i} = I_\mathrm{RRS}$", 
        xlabel = L"$\lambda\,[\mathrm{nm}]$", 
        ylabel = "")  
ax6 = Axis(f[1:4,9], title = L"$I_\mathrm{i}\,\,\mathrm{Window 1}$", 
        xlabel = L"$\lambda\,[\mathrm{nm}]$", 
        ylabel = "")        

ir=1 # r₀=0.1
iz=1 # z₀=2 km 
iA=1 # sza=0
# Effect of aerosol AOD for varying surface brightnesses
iρ=1 # ρ=0
palette = cgrad(:dense)   # Create a color gradient from "inferno"
clist0 = palette[LinRange(0, 1, length(ρ)+2)]  # Sample 10 colors
palette = cgrad(:algae)   # Create a color gradient from "inferno"
clist1 = palette[LinRange(0, 1, length(ρ)+2)]
palette = cgrad(:amp)   # Create a color gradient from "inferno"
clist2 = palette[LinRange(0, 1, length(ρ)+2)]
for iρ=1:length(ρ)
    iτ = 2 # AOD=0.5
    lines!(ax1, oco_λ, oco_ΔI[iτ,iz,ir,iρ,iA,:], color=clist2[iρ+1], linewidth=0.75,  label=L"$\mathrm{AOD}=$"*string(τ_ref[iτ])*L"$\,\,ρ=$"*string(ρ[iρ]))
    lines!(ax2, oco_λ[ind1_hires], oco_ΔI[iτ,iz,ir,iρ,iA,ind1_hires], color=clist2[iρ+1], linewidth=0.75,  label=L"$\mathrm{AOD}=$"*string(τ_ref[iτ])*L"$\,\,ρ=$"*string(ρ[iρ]))
    lines!(ax3, oco_λ, oco_eΔI[iτ,iz,ir,iρ,iA,:], color=clist2[iρ+1], linewidth=0.75,  label=L"$\mathrm{AOD}=$"*string(τ_ref[iτ])*L"$\,\,ρ=$"*string(ρ[iρ]))
    lines!(ax4, oco_λ[ind1_hires], oco_eΔI[iτ,iz,ir,iρ,iA,ind1_hires], color=clist2[iρ+1], linewidth=0.75,  label=L"$\mathrm{AOD}=$"*string(τ_ref[iτ])*L"$\,\,ρ=$"*string(ρ[iρ]))
    lines!(ax5, oco_λ, oco_ieI[iτ,iz,ir,iρ,iA,:], color=clist2[iρ+1], linewidth=0.75,  label=L"$\mathrm{AOD}=$"*string(τ_ref[iτ])*L"$\,\,ρ=$"*string(ρ[iρ]))
    lines!(ax6, oco_λ[ind1_hires], oco_ieI[iτ,iz,ir,iρ,iA,ind1_hires], color=clist2[iρ+1], linewidth=0.75,  label=L"$\mathrm{AOD}=$"*string(τ_ref[iτ])*L"$\,\,ρ=$"*string(ρ[iρ]))
    
    iτ = 1
    lines!(ax1, oco_λ, oco_ΔI[iτ,iz,ir,iρ,iA,:], color=clist1[iρ+1], linewidth=0.75,  label=L"$\mathrm{AOD}=$"*string(τ_ref[iτ])*L"$\,\,ρ=$"*string(ρ[iρ]))
    lines!(ax2, oco_λ[ind1_hires], oco_ΔI[iτ,iz,ir,iρ,iA,ind1_hires], color=clist1[iρ+1], linewidth=0.75,  label=L"$\mathrm{AOD}=$"*string(τ_ref[iτ])*L"$\,\,ρ=$"*string(ρ[iρ]))
    lines!(ax3, oco_λ, oco_eΔI[iτ,iz,ir,iρ,iA,:], color=clist1[iρ+1], linewidth=0.75,  label=L"$\mathrm{AOD}=$"*string(τ_ref[iτ])*L"$\,\,ρ=$"*string(ρ[iρ]))
    lines!(ax4, oco_λ[ind1_hires], oco_eΔI[iτ,iz,ir,iρ,iA,ind1_hires], color=clist1[iρ+1], linewidth=0.75,  label=L"$\mathrm{AOD}=$"*string(τ_ref[iτ])*L"$\,\,ρ=$"*string(ρ[iρ]))
    lines!(ax5, oco_λ, oco_ieI[iτ,iz,ir,iρ,iA,:], color=clist1[iρ+1], linewidth=0.75,  label=L"$\mathrm{AOD}=$"*string(τ_ref[iτ])*L"$\,\,ρ=$"*string(ρ[iρ]))
    lines!(ax6, oco_λ[ind1_hires], oco_ieI[iτ,iz,ir,iρ,iA,ind1_hires], color=clist1[iρ+1], linewidth=0.75,  label=L"$\mathrm{AOD}=$"*string(τ_ref[iτ])*L"$\,\,ρ=$"*string(ρ[iρ]))
    
    lines!(ax1, oco_λ, oco_ΔI0[iρ,iA,:], color=clist0[iρ+1], linewidth=0.75,  label=L"$\mathrm{AOD}=$"*string(τ_ref[iτ])*L"$\,\,ρ=$"*string(ρ[iρ]))
    lines!(ax2, oco_λ[ind1_hires], oco_ΔI0[iρ,iA,ind1_hires], color=clist0[iρ+1], linewidth=0.75,  label=L"$\mathrm{AOD}=$"*string(τ_ref[iτ])*L"$\,\,ρ=$"*string(ρ[iρ]))
    lines!(ax3, oco_λ, oco_eΔI0[iρ,iA,:], color=clist0[iρ+1], linewidth=0.75,  label=L"$\mathrm{AOD}=$"*string(τ_ref[iτ])*L"$\,\,ρ=$"*string(ρ[iρ]))
    lines!(ax4, oco_λ[ind1_hires], oco_eΔI0[iρ,iA,ind1_hires], color=clist0[iρ+1], linewidth=0.75,  label=L"$\mathrm{AOD}=$"*string(τ_ref[iτ])*L"$\,\,ρ=$"*string(ρ[iρ]))
    lines!(ax5, oco_λ, oco_ieI0[iρ,iA,:], color=clist0[iρ+1], linewidth=0.75,  label=L"$\mathrm{AOD}=$"*string(τ_ref[iτ])*L"$\,\,ρ=$"*string(ρ[iρ]))
    lines!(ax6, oco_λ[ind1_hires], oco_ieI0[iρ,iA,ind1_hires], color=clist0[iρ+1], linewidth=0.75,  label=L"$\mathrm{AOD}=$"*string(τ_ref[iτ])*L"$\,\,ρ=$"*string(ρ[iρ]))
    
    #lines!(ax1, oco_λ, oco_ΔI0[iρ,iA,:], color=clist0[iρ+1], label=L"$\mathrm{AOD}=0\,\,ρ=$"*string(ρ[iρ]))
end
# Manually create a legend entry with a larger marker
legend_line21 = LineElement(color=clist2[2]) #lines!(ax1, λ[1:2], ΔI[1,1,1,1,1,1:2],color=clist2[2],visible=false)
legend_line22 = LineElement(color=clist2[3]) #lines!(ax1, λ[1:2], ΔI[1,1,1,1,1,1:2],color=clist2[3],visible=false)
legend_line23 = LineElement(color=clist2[4]) #lines!(ax1, λ[1:2], ΔI[1,1,1,1,1,1:2],color=clist2[4],visible=false)
legend_line24 = LineElement(color=clist2[5]) #lines!(ax1, λ[1:2], ΔI[1,1,1,1,1,1:2],color=clist2[5],visible=false)

legend_line11 = LineElement(color=clist1[2]) #lines!(ax1, λ[1:2], ΔI[1,1,1,1,1,1:2],color=clist1[2],visible=false)
legend_line12 = LineElement(color=clist1[3]) #lines!(ax1, λ[1:2], ΔI[1,1,1,1,1,1:2],color=clist1[3],visible=false)
legend_line13 = LineElement(color=clist1[4]) #lines!(ax1, λ[1:2], ΔI[1,1,1,1,1,1:2],color=clist1[4],visible=false)
legend_line14 = LineElement(color=clist1[5]) #lines!(ax1, λ[1:2], ΔI[1,1,1,1,1,1:2],color=clist1[5],visible=false)

legend_line01 = LineElement(color=clist0[2]) #lines!(ax1, λ[1:2], ΔI0[1,1,1:2],color=clist0[2],visible=false)
legend_line02 = LineElement(color=clist0[3]) #lines!(ax1, λ[1:2], ΔI0[1,1,1:2],color=clist0[3],visible=false)
legend_line03 = LineElement(color=clist0[4]) #lines!(ax1, λ[1:2], ΔI0[1,1,1:2],color=clist0[4],visible=false)
legend_line04 = LineElement(color=clist0[5]) #lines!(ax1, λ[1:2], ΔI0[1,1,1:2],color=clist0[5],visible=false)

#=legend_labels = [
    L"\mathrm{AOD}=0\,\,\rho=" * LaTeXString(string(ρ[i])) for i in 1:4
] ∪ [
    L"\mathrm{AOD}=" * LaTeXString(string(τ_ref[1])) * L"\,\,\rho=" * LaTeXString(string(ρ[i])) for i in 1:4
] ∪ [
    L"\mathrm{AOD}=" * LaTeXString(string(τ_ref[2])) * L"\,\,\rho=" * LaTeXString(string(ρ[i])) for i in 1:4
]=#
legend_labels = [
    latexstring(L"\mathrm{AOD}=0\,\,\rho=" * latexstring(string(ρ[1]))),
    latexstring(L"\mathrm{AOD}=" * latexstring(string(τ_ref[1])) * L"\,\,\rho=" * latexstring(string(ρ[1]))),
    latexstring(L"\mathrm{AOD}=" * latexstring(string(τ_ref[2])) * L"\,\,\rho=" * latexstring(string(ρ[1]))),
    latexstring(L"\mathrm{AOD}=0\,\,\rho=" * latexstring(string(ρ[2]))),
    latexstring(L"\mathrm{AOD}=" * latexstring(string(τ_ref[1])) * L"\,\,\rho=" * latexstring(string(ρ[2]))),
    latexstring(L"\mathrm{AOD}=" * latexstring(string(τ_ref[2])) * L"\,\,\rho=" * latexstring(string(ρ[2]))),
    latexstring(L"\mathrm{AOD}=0\,\,\rho=" * latexstring(string(ρ[3]))),
    latexstring(L"\mathrm{AOD}=" * latexstring(string(τ_ref[1])) * L"\,\,\rho=" * latexstring(string(ρ[3]))),
    latexstring(L"\mathrm{AOD}=" * latexstring(string(τ_ref[2])) * L"\,\,\rho=" * latexstring(string(ρ[3]))),
    latexstring(L"\mathrm{AOD}=0\,\,\rho=" * latexstring(string(ρ[4]))),
    latexstring(L"\mathrm{AOD}=" * latexstring(string(τ_ref[1])) * L"\,\,\rho=" * latexstring(string(ρ[4]))),
    latexstring(L"\mathrm{AOD}=" * latexstring(string(τ_ref[2])) * L"\,\,\rho=" * latexstring(string(ρ[4])))
]

legend_elements = [
    legend_line01, legend_line11, legend_line21, 
    legend_line02, legend_line12, legend_line22, 
    legend_line03, legend_line13, legend_line23, 
    legend_line04, legend_line14, legend_line24]
legend = Legend(f, legend_elements, legend_labels, 
           orientation = :vertical,  # Ensures elements are stacked properly
           tellwidth = false,  # Avoids stretching legend width
           nbanks = 3,  # Specifies 4 columns (i.e., 3 rows automatically)
           halign = :center,
           valign = :center
       )

f[5:6, :] = legend
save("ocoRRScorr_wrt_AOD_z"*z_str[iz]*"_r"*r_str[ir]*"_sza"*string(sza[iA])*".png", f)

# Effect of z₀
ir=1 # r₀=0.1
iτ=2 # z₀=2 km 
iA=1 # sza=0
# Effect of aerosol AOD for varying surface brightnesses
iρ=1 # ρ=0
for iρ=1:length(ρ)
    iz = 1 # z₀=2
    plot(λ, I[iτ,iz,ir,iρ,iA,:],label=L"$z_0=$"*string(aer_z₀[iz])*L"$\mathrm{km}\,\,ρ=$"*string(ρ[iρ]))
    iz = 2 # z₀=12
    plot!(λ, I[iτ,iz,ir,iρ,iA,:],label=L"$z_0=$"*string(aer_z₀[iz])*L"$\mathrm{km}\,\,ρ=$"*string(ρ[iρ]))
    plot!(λ, I0[iρ,iA,:],label=L"$\mathrm{AOD}=0\,\,ρ=$"*string(ρ[iρ]))

    savefig("RRScorr_wrt_z_AOD"*AOD_str[iτ]*"_r"*r_str[ir]*"_sza"*string(sza[iA])*"_alb"*ρ_str[iρ]*".png")
end

# Effect of z₀
ir=1 # r₀=0.1
iτ=2 # z₀=2 km 
iA=1 # sza=0
# Effect of aerosol AOD for varying surface brightnesses
iρ=1 # ρ=0
for iρ=1:length(ρ)
    iz = 1 # z₀=2
    plot(oco_λ, oco_I[iτ,iz,ir,iρ,iA,:],label=L"$z_0=$"*string(aer_z₀[iz])*L"$\mathrm{km}\,\,ρ=$"*string(ρ[iρ]))
    iz = 2 # z₀=12
    plot!(oco_λ, oco_I[iτ,iz,ir,iρ,iA,:],label=L"$z_0=$"*string(aer_z₀[iz])*L"$\mathrm{km}\,\,ρ=$"*string(ρ[iρ]))
    plot!(oco_λ, oco_I0[iρ,iA,:],label=L"$\mathrm{AOD}=0\,\,ρ=$"*string(ρ[iρ]))

    savefig("ocoRRScorr_wrt_z_AOD"*AOD_str[iτ]*"_r"*r_str[ir]*"_sza"*string(sza[iA])*"_alb"*ρ_str[iρ]*".png")
end

# Effect of r₀
iz=1 # r₀=0.1
iτ=2 # z₀=2 km 
iA=1 # sza=0
# Effect of aerosol AOD for varying surface brightnesses
iρ=1 # ρ=0
for iρ=1:length(ρ)
    ir = 1 # r₀=0.1
    plot(λ, I[iτ,iz,ir,iρ,iA,:],label=L"$r_0=$"*string(aer_r₀[ir])*L"$\mu\mathrm{m}\,\,ρ=$"*string(ρ[iρ]))
    ir = 2 # z₀=12
    plot!(λ, I[iτ,iz,ir,iρ,iA,:],label=L"$r_0=$"*string(aer_r₀[ir])*L"$\mu\mathrm{m}\,\,ρ=$"*string(ρ[iρ]))
    plot!(λ, I0[iρ,iA,:],label=L"$\mathrm{AOD}=0\,\,ρ=$"*string(ρ[iρ]))

    savefig("RRScorr_wrt_r_AOD"*AOD_str[iτ]*"_z"*z_str[iz]*"_sza"*string(sza[iA])*"_alb"*ρ_str[iρ]*".png")
end

# Effect of r₀
iz=1 # r₀=0.1
iτ=2 # z₀=2 km 
iA=1 # sza=0
# Effect of aerosol AOD for varying surface brightnesses
iρ=1 # ρ=0
for iρ=1:length(ρ)
    ir = 1 # r₀=0.1
    plot(oco_λ, oco_I[iτ,iz,ir,iρ,iA,:],label=L"$r_0=$"*string(aer_r₀[ir])*L"$\mu\mathrm{m}\,\,ρ=$"*string(ρ[iρ]))
    ir = 2 # z₀=12
    plot!(oco_λ, oco_I[iτ,iz,ir,iρ,iA,:],label=L"$r_0=$"*string(aer_r₀[ir])*L"$\mu\mathrm{m}\,\,ρ=$"*string(ρ[iρ]))
    plot!(oco_λ, oco_I0[iρ,iA,:],label=L"$\mathrm{AOD}=0\,\,ρ=$"*string(ρ[iρ]))

    savefig("ocoRRScorr_wrt_r_AOD"*AOD_str[iτ]*"_z"*z_str[iz]*"_sza"*string(sza[iA])*"_alb"*ρ_str[iρ]*".png")
end