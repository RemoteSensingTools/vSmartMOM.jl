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

# Benchmarks: http://web.gps.caltech.edu/~vijay/Rayleigh_Scattering_Tables/STOKES/
##

# Load YAML files into parameter struct
parameters = 
    parameters_from_yaml("test/test_parameters/ParamsCarbonI.yaml");
#parameters.depol = 0.041362343961163395 #0.028 #(0.028: Cabannes), (0.1032: Rayleigh) 
    #parameters = parameters_from_yaml("test/test_parameters/O2Parameters2.yaml");
# Create model struct (precomputes optical properties) from parameters
model      = model_from_parameters(parameters);

iBand = 1
FT = Float64
T_sun = 5777. # K 
Tsolar = solar_transmission_from_file("/home/sanghavi/code/github/vSmartMOM.jl/src/SolarModel/solar.out")    
Tsolar_interp = LinearInterpolation(Tsolar[4:end, 1], Tsolar[4:end, 2])
x = -40:0.05:40
#kernel = InstrumentOperator.create_instrument_kernel(Normal(0, 12.5), x) #defining a Gaussian kernel for convulution in wavenumber space
kernel = InstrumentOperator.create_instrument_kernel(Normal(0, 1.43322), x)
I_conv = [];
λCI = [];
I_CI = [];
λ = [];
Δλ = 7.0 #[nm] EMIT sampling 
for iBand=1:length(model.params.spec_bands)
#### Compute all Raman properties
    ν = model.params.spec_bands[iBand]
    ν̃ = mean(ν);
    push!(λCI, 1e7/ν[end]:Δλ:1e7/ν[1])
    
    # Find central reference index for RRS:
    i_ref = argmin(abs.(ν .- ν̃))
    # TODO_VS: λ_vs_in (get input)
    # TODO_VS: ν_vs_in (convert to wavenumbers)
    # Effective temperature for Raman calculations
    effT = 300.  #(model.profile.vcd_dry' * model.profile.T) / sum(model.profile.vcd_dry);
    RS_type = InelasticScattering.noRS(
        fscattRayl  = [FT(1)],
        ϖ_Cabannes  = [FT(1)], 
        bandSpecLim = [],
        iBand       = [1],
        F₀          = zeros(FT,1,1));

    P = planck_spectrum_wn(T_sun, ν)
    F₀ = zeros(length(P));
    RS_type.F₀ = zeros(model.params.polarization_type.n, length(P))
    for i=1:length(P)
        sol_trans = Tsolar_interp(ν[i]);
        F₀[i] = sol_trans * P[i];
        RS_type.F₀[1,i] = F₀[i];
    end 
    R, T, ieR, ieT = CoreRT.rt_run_test(RS_type,model,iBand);
    #R_ss, T_ss, ieR_ss, ieT_ss = CoreRT.rt_run_test_ss(RS_type,model,iBand);

    #===Convolution of hires spectral simulations to instrument grid===#
    #I_conv = InstrumentOperator.conv_spectra(kernel, )
    #I_conv_noRS = imfilter(RnoRS[1,1,:], kernel)
    convfct = 1e7./ν.^2   # to convert radiance units from mW/m²-str-cm⁻¹ to mW/m²-str-nm
    push!(I_conv, imfilter(R[1,1,:].*convfct, kernel))
    #ieI_conv = imfilter(ieR[1,1,:], kernel)
    tmpI_EMIT=[];
    for j=1:length(λEMIT[iBand])
        ii = findall( x->1e7/(λEMIT[iBand][j] + Δλ/2)<x<1e7/(λEMIT[iBand][j] - Δλ/2), ν)
        tmpI=0;
        δλ =0;
        for k=1:length(ii)
            if ii[k]==1
                dλ = (1e7/ν[ii[k]] - 1e7/ν[ii[k]+1])
            elseif ii[k]==length(ν)        
                dλ = (1e7/ν[ii[k]-1] - 1e7/ν[ii[k]])
            else 
                dλ = 0.5*(1e7/ν[ii[k]-1] - 1e7/ν[ii[k]+1])
            end
            tmpI += I_conv[iBand][ii[k]]*dλ
            δλ += dλ
        end
        push!(tmpI_EMIT, tmpI/δλ)
        #push!(tmpI_EMIT, sum(I_conv[iBand][ii]))
    end
    push!(I_EMIT, tmpI_EMIT)
    push!(λ, 1e7./model.params.spec_bands[iBand])
end

#=
for iBand=1:length(model.params.spec_bands)
    writedlm("/home/sanghavi/EMIT/I_EMIT_aer0p0_hiCH4_band"*string(iBand)*".dat", I_EMIT_aer0p0_hiCH4[iBand])
    writedlm("/home/sanghavi/EMIT/I_EMIT_aer0p0_loCH4_band"*string(iBand)*".dat", I_EMIT_aer0p0_loCH4[iBand])
    writedlm("/home/sanghavi/EMIT/I_EMIT_aer0p5_900hPa_hiCH4_band"*string(iBand)*".dat", I_EMIT_aer0p5_900hPa_hiCH4[iBand])
    writedlm("/home/sanghavi/EMIT/I_EMIT_aer5p0_900hPa_hiCH4_band"*string(iBand)*".dat", I_EMIT_aer5p0_900hPa_hiCH4[iBand])
    writedlm("/home/sanghavi/EMIT/I_EMIT_aer0p5_900hPa_loCH4_band"*string(iBand)*".dat", I_EMIT_aer0p5_900hPa_loCH4[iBand])
    writedlm("/home/sanghavi/EMIT/I_EMIT_aer5p0_900hPa_loCH4_band"*string(iBand)*".dat", I_EMIT_aer5p0_900hPa_loCH4[iBand])
    writedlm("/home/sanghavi/EMIT/I_EMIT_aer0p5_500hPa_hiCH4_band"*string(iBand)*".dat", I_EMIT_aer0p5_500hPa_hiCH4[iBand])
    writedlm("/home/sanghavi/EMIT/I_EMIT_aer5p0_500hPa_hiCH4_band"*string(iBand)*".dat", I_EMIT_aer5p0_500hPa_hiCH4[iBand])
    writedlm("/home/sanghavi/EMIT/I_EMIT_aer0p5_500hPa_loCH4_band"*string(iBand)*".dat", I_EMIT_aer0p5_500hPa_loCH4[iBand])
    writedlm("/home/sanghavi/EMIT/I_EMIT_aer5p0_500hPa_loCH4_band"*string(iBand)*".dat", I_EMIT_aer5p0_500hPa_loCH4[iBand])
    writedlm("/home/sanghavi/EMIT/I_conv_aer0p0_hiCH4_band"*string(iBand)*".dat", I_conv_aer0p0_hiCH4[iBand])
    writedlm("/home/sanghavi/EMIT/I_conv_aer0p0_loCH4_band"*string(iBand)*".dat", I_conv_aer0p0_loCH4[iBand])
    writedlm("/home/sanghavi/EMIT/I_conv_aer0p5_900hPa_hiCH4_band"*string(iBand)*".dat", I_conv_aer0p5_900hPa_hiCH4[iBand])
    writedlm("/home/sanghavi/EMIT/I_conv_aer5p0_900hPa_hiCH4_band"*string(iBand)*".dat", I_conv_aer5p0_900hPa_hiCH4[iBand])
    writedlm("/home/sanghavi/EMIT/I_conv_aer0p5_900hPa_loCH4_band"*string(iBand)*".dat", I_conv_aer0p5_900hPa_loCH4[iBand])
    writedlm("/home/sanghavi/EMIT/I_conv_aer5p0_900hPa_loCH4_band"*string(iBand)*".dat", I_conv_aer5p0_900hPa_loCH4[iBand])
    writedlm("/home/sanghavi/EMIT/I_conv_aer0p5_500hPa_hiCH4_band"*string(iBand)*".dat", I_conv_aer0p5_500hPa_hiCH4[iBand])
    writedlm("/home/sanghavi/EMIT/I_conv_aer5p0_500hPa_hiCH4_band"*string(iBand)*".dat", I_conv_aer5p0_500hPa_hiCH4[iBand])
    writedlm("/home/sanghavi/EMIT/I_conv_aer0p5_500hPa_loCH4_band"*string(iBand)*".dat", I_conv_aer0p5_500hPa_loCH4[iBand])
    writedlm("/home/sanghavi/EMIT/I_conv_aer5p0_500hPa_loCH4_band"*string(iBand)*".dat", I_conv_aer5p0_500hPa_loCH4[iBand])
    #writedlm("/home/sanghavi/EMIT/lambda_band"*string(iBand)*".dat", 1e7./model.params.spec_bands[iBand])
    #writedlm("/home/sanghavi/EMIT/lambdaEMIT_band"*string(iBand)*".dat", λEMIT[iBand])
end
=#
##=
aI_EMIT_aer0p5_500hPa_hiCH4 = []; aI_conv_aer0p5_500hPa_hiCH4 = [];
aI_EMIT_aer5p0_500hPa_hiCH4 = []; aI_conv_aer5p0_500hPa_hiCH4 = [];
aI_EMIT_aer0p5_900hPa_hiCH4 = []; aI_conv_aer0p5_900hPa_hiCH4 = [];
aI_EMIT_aer5p0_900hPa_hiCH4 = []; aI_conv_aer5p0_900hPa_hiCH4 = [];
aI_EMIT_aer0p0_hiCH4 = []; aI_conv_aer0p0_hiCH4 = [];
aI_EMIT_aer0p0_loCH4 = []; aI_conv_aer0p0_loCH4 = [];
aI_EMIT_aer0p5_900hPa_loCH4 = []; aI_conv_aer0p5_900hPa_loCH4 = [];
aI_EMIT_aer5p0_900hPa_loCH4 = []; aI_conv_aer5p0_900hPa_loCH4 = [];
aI_EMIT_aer5p0_500hPa_loCH4 = []; aI_conv_aer5p0_500hPa_loCH4 = [];
aI_EMIT_aer0p5_500hPa_loCH4 = []; aI_conv_aer0p5_500hPa_loCH4 = [];

bI_EMIT_aer0p5_500hPa_hiCH4 = []; bI_conv_aer0p5_500hPa_hiCH4 = [];
bI_EMIT_aer5p0_500hPa_hiCH4 = []; bI_conv_aer5p0_500hPa_hiCH4 = [];
bI_EMIT_aer0p5_900hPa_hiCH4 = []; bI_conv_aer0p5_900hPa_hiCH4 = [];
bI_EMIT_aer5p0_900hPa_hiCH4 = []; bI_conv_aer5p0_900hPa_hiCH4 = [];
bI_EMIT_aer0p0_hiCH4 = []; bI_conv_aer0p0_hiCH4 = [];
bI_EMIT_aer0p0_loCH4 = []; bI_conv_aer0p0_loCH4 = [];
bI_EMIT_aer0p5_900hPa_loCH4 = []; bI_conv_aer0p5_900hPa_loCH4 = [];
bI_EMIT_aer5p0_900hPa_loCH4 = []; bI_conv_aer5p0_900hPa_loCH4 = [];
bI_EMIT_aer5p0_500hPa_loCH4 = []; bI_conv_aer5p0_500hPa_loCH4 = [];
bI_EMIT_aer0p5_500hPa_loCH4 = []; bI_conv_aer0p5_500hPa_loCH4 = [];

cI_EMIT_aer0p5_500hPa_hiCH4 = []; cI_conv_aer0p5_500hPa_hiCH4 = [];
cI_EMIT_aer5p0_500hPa_hiCH4 = []; cI_conv_aer5p0_500hPa_hiCH4 = [];
cI_EMIT_aer0p5_900hPa_hiCH4 = []; cI_conv_aer0p5_900hPa_hiCH4 = [];
cI_EMIT_aer5p0_900hPa_hiCH4 = []; cI_conv_aer5p0_900hPa_hiCH4 = [];
cI_EMIT_aer0p0_hiCH4 = []; cI_conv_aer0p0_hiCH4 = [];
cI_EMIT_aer0p0_loCH4 = []; cI_conv_aer0p0_loCH4 = [];
cI_EMIT_aer0p5_900hPa_loCH4 = []; cI_conv_aer0p5_900hPa_loCH4 = [];
cI_EMIT_aer5p0_900hPa_loCH4 = []; cI_conv_aer5p0_900hPa_loCH4 = [];
cI_EMIT_aer5p0_500hPa_loCH4 = []; cI_conv_aer5p0_500hPa_loCH4 = [];
cI_EMIT_aer0p5_500hPa_loCH4 = []; cI_conv_aer0p5_500hPa_loCH4 = [];

λ = [];
λEMIT=[];
for iBand=1:5 #length(model.params.spec_bands)
    tmp = readdlm("/home/sanghavi/EMIT/arid_0p0_0p2_0p4_0p4_0p4/I_EMIT_aer0p0_hiCH4_band"*string(iBand)*".dat")
    push!(aI_EMIT_aer0p0_hiCH4, tmp)
    tmp = readdlm("/home/sanghavi/EMIT/arid_0p0_0p2_0p4_0p4_0p4/I_EMIT_aer0p0_loCH4_band"*string(iBand)*".dat")
    push!(aI_EMIT_aer0p0_loCH4, tmp)
    tmp = readdlm("/home/sanghavi/EMIT/arid_0p0_0p2_0p4_0p4_0p4/I_EMIT_aer0p5_900hPa_hiCH4_band"*string(iBand)*".dat")
    push!(aI_EMIT_aer0p5_900hPa_hiCH4, tmp)
    tmp = readdlm("/home/sanghavi/EMIT/arid_0p0_0p2_0p4_0p4_0p4/I_EMIT_aer5p0_900hPa_hiCH4_band"*string(iBand)*".dat")
    push!(aI_EMIT_aer5p0_900hPa_hiCH4, tmp)
    tmp = readdlm("/home/sanghavi/EMIT/arid_0p0_0p2_0p4_0p4_0p4/I_EMIT_aer0p5_900hPa_loCH4_band"*string(iBand)*".dat")
    push!(aI_EMIT_aer0p5_900hPa_loCH4, tmp)
    tmp = readdlm("/home/sanghavi/EMIT/arid_0p0_0p2_0p4_0p4_0p4/I_EMIT_aer5p0_900hPa_loCH4_band"*string(iBand)*".dat")
    push!(aI_EMIT_aer5p0_900hPa_loCH4, tmp)
    tmp = readdlm("/home/sanghavi/EMIT/arid_0p0_0p2_0p4_0p4_0p4/I_EMIT_aer0p5_500hPa_hiCH4_band"*string(iBand)*".dat")
    push!(aI_EMIT_aer0p5_500hPa_hiCH4, tmp)
    tmp = readdlm("/home/sanghavi/EMIT/arid_0p0_0p2_0p4_0p4_0p4/I_EMIT_aer5p0_500hPa_hiCH4_band"*string(iBand)*".dat")
    push!(aI_EMIT_aer5p0_500hPa_hiCH4, tmp)    
    tmp = readdlm("/home/sanghavi/EMIT/arid_0p0_0p2_0p4_0p4_0p4/I_EMIT_aer0p5_500hPa_loCH4_band"*string(iBand)*".dat")
    push!(aI_EMIT_aer0p5_500hPa_loCH4, tmp)
    tmp = readdlm("/home/sanghavi/EMIT/arid_0p0_0p2_0p4_0p4_0p4/I_EMIT_aer5p0_500hPa_loCH4_band"*string(iBand)*".dat")
    push!(aI_EMIT_aer5p0_500hPa_loCH4, tmp)    
    tmp = readdlm("/home/sanghavi/EMIT/arid_0p0_0p2_0p4_0p4_0p4/I_conv_aer0p0_hiCH4_band"*string(iBand)*".dat")
    push!(aI_conv_aer0p0_hiCH4, tmp)    
    tmp = readdlm("/home/sanghavi/EMIT/arid_0p0_0p2_0p4_0p4_0p4/I_conv_aer0p0_loCH4_band"*string(iBand)*".dat")
    push!(aI_conv_aer0p0_loCH4, tmp)
    tmp = readdlm("/home/sanghavi/EMIT/arid_0p0_0p2_0p4_0p4_0p4/I_conv_aer0p5_900hPa_hiCH4_band"*string(iBand)*".dat")
    push!(aI_conv_aer0p5_900hPa_hiCH4, tmp)
    tmp = readdlm("/home/sanghavi/EMIT/arid_0p0_0p2_0p4_0p4_0p4/I_conv_aer5p0_900hPa_hiCH4_band"*string(iBand)*".dat")
    push!(aI_conv_aer5p0_900hPa_hiCH4, tmp)
    tmp = readdlm("/home/sanghavi/EMIT/arid_0p0_0p2_0p4_0p4_0p4/I_conv_aer0p5_900hPa_loCH4_band"*string(iBand)*".dat")
    push!(aI_conv_aer0p5_900hPa_loCH4, tmp)
    tmp = readdlm("/home/sanghavi/EMIT/arid_0p0_0p2_0p4_0p4_0p4/I_conv_aer5p0_900hPa_loCH4_band"*string(iBand)*".dat")
    push!(aI_conv_aer5p0_900hPa_loCH4, tmp)
    tmp = readdlm("/home/sanghavi/EMIT/arid_0p0_0p2_0p4_0p4_0p4/I_conv_aer0p5_500hPa_hiCH4_band"*string(iBand)*".dat")
    push!(aI_conv_aer0p5_500hPa_hiCH4, tmp)
    tmp = readdlm("/home/sanghavi/EMIT/arid_0p0_0p2_0p4_0p4_0p4/I_conv_aer5p0_500hPa_hiCH4_band"*string(iBand)*".dat")
    push!(aI_conv_aer5p0_500hPa_hiCH4, tmp)
    tmp = readdlm("/home/sanghavi/EMIT/arid_0p0_0p2_0p4_0p4_0p4/I_conv_aer0p5_500hPa_loCH4_band"*string(iBand)*".dat")
    push!(aI_conv_aer0p5_500hPa_loCH4, tmp)
    tmp = readdlm("/home/sanghavi/EMIT/arid_0p0_0p2_0p4_0p4_0p4/I_conv_aer5p0_500hPa_loCH4_band"*string(iBand)*".dat")
    push!(aI_conv_aer5p0_500hPa_loCH4, tmp)

    tmp = readdlm("/home/sanghavi/EMIT/ocean_0p0_0p0_0p0_0p0_0p0/I_EMIT_aer0p0_hiCH4_band"*string(iBand)*".dat")
    push!(bI_EMIT_aer0p0_hiCH4, tmp)
    tmp = readdlm("/home/sanghavi/EMIT/ocean_0p0_0p0_0p0_0p0_0p0/I_EMIT_aer0p0_loCH4_band"*string(iBand)*".dat")
    push!(bI_EMIT_aer0p0_loCH4, tmp)
    tmp = readdlm("/home/sanghavi/EMIT/ocean_0p0_0p0_0p0_0p0_0p0/I_EMIT_aer0p5_900hPa_hiCH4_band"*string(iBand)*".dat")
    push!(bI_EMIT_aer0p5_900hPa_hiCH4, tmp)
    tmp = readdlm("/home/sanghavi/EMIT/ocean_0p0_0p0_0p0_0p0_0p0/I_EMIT_aer5p0_900hPa_hiCH4_band"*string(iBand)*".dat")
    push!(bI_EMIT_aer5p0_900hPa_hiCH4, tmp)
    tmp = readdlm("/home/sanghavi/EMIT/ocean_0p0_0p0_0p0_0p0_0p0/I_EMIT_aer0p5_900hPa_loCH4_band"*string(iBand)*".dat")
    push!(bI_EMIT_aer0p5_900hPa_loCH4, tmp)
    tmp = readdlm("/home/sanghavi/EMIT/ocean_0p0_0p0_0p0_0p0_0p0/I_EMIT_aer5p0_900hPa_loCH4_band"*string(iBand)*".dat")
    push!(bI_EMIT_aer5p0_900hPa_loCH4, tmp)
    tmp = readdlm("/home/sanghavi/EMIT/ocean_0p0_0p0_0p0_0p0_0p0/I_EMIT_aer0p5_500hPa_hiCH4_band"*string(iBand)*".dat")
    push!(bI_EMIT_aer0p5_500hPa_hiCH4, tmp)
    tmp = readdlm("/home/sanghavi/EMIT/ocean_0p0_0p0_0p0_0p0_0p0/I_EMIT_aer5p0_500hPa_hiCH4_band"*string(iBand)*".dat")
    push!(bI_EMIT_aer5p0_500hPa_hiCH4, tmp)    
    tmp = readdlm("/home/sanghavi/EMIT/ocean_0p0_0p0_0p0_0p0_0p0/I_EMIT_aer0p5_500hPa_loCH4_band"*string(iBand)*".dat")
    push!(bI_EMIT_aer0p5_500hPa_loCH4, tmp)
    tmp = readdlm("/home/sanghavi/EMIT/ocean_0p0_0p0_0p0_0p0_0p0/I_EMIT_aer5p0_500hPa_loCH4_band"*string(iBand)*".dat")
    push!(bI_EMIT_aer5p0_500hPa_loCH4, tmp)    
    tmp = readdlm("/home/sanghavi/EMIT/ocean_0p0_0p0_0p0_0p0_0p0/I_conv_aer0p0_hiCH4_band"*string(iBand)*".dat")
    push!(bI_conv_aer0p0_hiCH4, tmp)    
    tmp = readdlm("/home/sanghavi/EMIT/ocean_0p0_0p0_0p0_0p0_0p0/I_conv_aer0p0_loCH4_band"*string(iBand)*".dat")
    push!(bI_conv_aer0p0_loCH4, tmp)
    tmp = readdlm("/home/sanghavi/EMIT/ocean_0p0_0p0_0p0_0p0_0p0/I_conv_aer0p5_900hPa_hiCH4_band"*string(iBand)*".dat")
    push!(bI_conv_aer0p5_900hPa_hiCH4, tmp)
    tmp = readdlm("/home/sanghavi/EMIT/ocean_0p0_0p0_0p0_0p0_0p0/I_conv_aer5p0_900hPa_hiCH4_band"*string(iBand)*".dat")
    push!(bI_conv_aer5p0_900hPa_hiCH4, tmp)
    tmp = readdlm("/home/sanghavi/EMIT/ocean_0p0_0p0_0p0_0p0_0p0/I_conv_aer0p5_900hPa_loCH4_band"*string(iBand)*".dat")
    push!(bI_conv_aer0p5_900hPa_loCH4, tmp)
    tmp = readdlm("/home/sanghavi/EMIT/ocean_0p0_0p0_0p0_0p0_0p0/I_conv_aer5p0_900hPa_loCH4_band"*string(iBand)*".dat")
    push!(bI_conv_aer5p0_900hPa_loCH4, tmp)
    tmp = readdlm("/home/sanghavi/EMIT/ocean_0p0_0p0_0p0_0p0_0p0/I_conv_aer0p5_500hPa_hiCH4_band"*string(iBand)*".dat")
    push!(bI_conv_aer0p5_500hPa_hiCH4, tmp)
    tmp = readdlm("/home/sanghavi/EMIT/ocean_0p0_0p0_0p0_0p0_0p0/I_conv_aer5p0_500hPa_hiCH4_band"*string(iBand)*".dat")
    push!(bI_conv_aer5p0_500hPa_hiCH4, tmp)
    tmp = readdlm("/home/sanghavi/EMIT/ocean_0p0_0p0_0p0_0p0_0p0/I_conv_aer0p5_500hPa_loCH4_band"*string(iBand)*".dat")
    push!(bI_conv_aer0p5_500hPa_loCH4, tmp)
    tmp = readdlm("/home/sanghavi/EMIT/ocean_0p0_0p0_0p0_0p0_0p0/I_conv_aer5p0_500hPa_loCH4_band"*string(iBand)*".dat")
    push!(bI_conv_aer5p0_500hPa_loCH4, tmp)

    tmp = readdlm("/home/sanghavi/EMIT/semiarid_0p1_0p1_0p6_0p2_0p2/I_EMIT_aer0p0_hiCH4_band"*string(iBand)*".dat")
    push!(cI_EMIT_aer0p0_hiCH4, tmp)
    tmp = readdlm("/home/sanghavi/EMIT/semiarid_0p1_0p1_0p6_0p2_0p2/I_EMIT_aer0p0_loCH4_band"*string(iBand)*".dat")
    push!(cI_EMIT_aer0p0_loCH4, tmp)
    tmp = readdlm("/home/sanghavi/EMIT/semiarid_0p1_0p1_0p6_0p2_0p2/I_EMIT_aer0p5_900hPa_hiCH4_band"*string(iBand)*".dat")
    push!(cI_EMIT_aer0p5_900hPa_hiCH4, tmp)
    tmp = readdlm("/home/sanghavi/EMIT/semiarid_0p1_0p1_0p6_0p2_0p2/I_EMIT_aer5p0_900hPa_hiCH4_band"*string(iBand)*".dat")
    push!(cI_EMIT_aer5p0_900hPa_hiCH4, tmp)
    tmp = readdlm("/home/sanghavi/EMIT/semiarid_0p1_0p1_0p6_0p2_0p2/I_EMIT_aer0p5_900hPa_loCH4_band"*string(iBand)*".dat")
    push!(cI_EMIT_aer0p5_900hPa_loCH4, tmp)
    tmp = readdlm("/home/sanghavi/EMIT/semiarid_0p1_0p1_0p6_0p2_0p2/I_EMIT_aer5p0_900hPa_loCH4_band"*string(iBand)*".dat")
    push!(cI_EMIT_aer5p0_900hPa_loCH4, tmp)
    tmp = readdlm("/home/sanghavi/EMIT/semiarid_0p1_0p1_0p6_0p2_0p2/I_EMIT_aer0p5_500hPa_hiCH4_band"*string(iBand)*".dat")
    push!(cI_EMIT_aer0p5_500hPa_hiCH4, tmp)
    tmp = readdlm("/home/sanghavi/EMIT/semiarid_0p1_0p1_0p6_0p2_0p2/I_EMIT_aer5p0_500hPa_hiCH4_band"*string(iBand)*".dat")
    push!(cI_EMIT_aer5p0_500hPa_hiCH4, tmp)    
    tmp = readdlm("/home/sanghavi/EMIT/semiarid_0p1_0p1_0p6_0p2_0p2/I_EMIT_aer0p5_500hPa_loCH4_band"*string(iBand)*".dat")
    push!(cI_EMIT_aer0p5_500hPa_loCH4, tmp)
    tmp = readdlm("/home/sanghavi/EMIT/semiarid_0p1_0p1_0p6_0p2_0p2/I_EMIT_aer5p0_500hPa_loCH4_band"*string(iBand)*".dat")
    push!(cI_EMIT_aer5p0_500hPa_loCH4, tmp)    
    tmp = readdlm("/home/sanghavi/EMIT/semiarid_0p1_0p1_0p6_0p2_0p2/I_conv_aer0p0_hiCH4_band"*string(iBand)*".dat")
    push!(cI_conv_aer0p0_hiCH4, tmp)    
    tmp = readdlm("/home/sanghavi/EMIT/semiarid_0p1_0p1_0p6_0p2_0p2/I_conv_aer0p0_loCH4_band"*string(iBand)*".dat")
    push!(cI_conv_aer0p0_loCH4, tmp)
    tmp = readdlm("/home/sanghavi/EMIT/semiarid_0p1_0p1_0p6_0p2_0p2/I_conv_aer0p5_900hPa_hiCH4_band"*string(iBand)*".dat")
    push!(cI_conv_aer0p5_900hPa_hiCH4, tmp)
    tmp = readdlm("/home/sanghavi/EMIT/semiarid_0p1_0p1_0p6_0p2_0p2/I_conv_aer5p0_900hPa_hiCH4_band"*string(iBand)*".dat")
    push!(cI_conv_aer5p0_900hPa_hiCH4, tmp)
    tmp = readdlm("/home/sanghavi/EMIT/semiarid_0p1_0p1_0p6_0p2_0p2/I_conv_aer0p5_900hPa_loCH4_band"*string(iBand)*".dat")
    push!(cI_conv_aer0p5_900hPa_loCH4, tmp)
    tmp = readdlm("/home/sanghavi/EMIT/semiarid_0p1_0p1_0p6_0p2_0p2/I_conv_aer5p0_900hPa_loCH4_band"*string(iBand)*".dat")
    push!(cI_conv_aer5p0_900hPa_loCH4, tmp)
    tmp = readdlm("/home/sanghavi/EMIT/semiarid_0p1_0p1_0p6_0p2_0p2/I_conv_aer0p5_500hPa_hiCH4_band"*string(iBand)*".dat")
    push!(cI_conv_aer0p5_500hPa_hiCH4, tmp)
    tmp = readdlm("/home/sanghavi/EMIT/semiarid_0p1_0p1_0p6_0p2_0p2/I_conv_aer5p0_500hPa_hiCH4_band"*string(iBand)*".dat")
    push!(cI_conv_aer5p0_500hPa_hiCH4, tmp)
    tmp = readdlm("/home/sanghavi/EMIT/semiarid_0p1_0p1_0p6_0p2_0p2/I_conv_aer0p5_500hPa_loCH4_band"*string(iBand)*".dat")
    push!(cI_conv_aer0p5_500hPa_loCH4, tmp)
    tmp = readdlm("/home/sanghavi/EMIT/semiarid_0p1_0p1_0p6_0p2_0p2/I_conv_aer5p0_500hPa_loCH4_band"*string(iBand)*".dat")
    push!(cI_conv_aer5p0_500hPa_loCH4, tmp)

    tmp = readdlm("/home/sanghavi/EMIT/lambda_band"*string(iBand)*".dat")
    push!(λ, tmp)
    tmp = readdlm("/home/sanghavi/EMIT/lambdaEMIT_band"*string(iBand)*".dat")
    push!(λEMIT, tmp)
end
##=#
l = @layout[a1 a2; a3 a4]
p1 = plot(λ[1], aI_conv_aer0p0_loCH4[1], ylabel="R [mW/m²-str-nm]", linewidth = 1, linecolor=:gray, label=false)
p1 = plot!(λ[1], aI_conv_aer0p5_900hPa_loCH4[1], linewidth = 2, linecolor=:gray, label=false)
p1 = plot!(λ[1], aI_conv_aer5p0_900hPa_loCH4[1], linewidth = 3, linecolor=:gray, label=false)
p1 = scatter!(λEMIT[1], aI_EMIT_aer0p0_loCH4[1], marker = :circle, markersize = 4, markercolor = :red, label=false)
p1 = scatter!(λEMIT[1], aI_EMIT_aer0p5_900hPa_loCH4[1], marker = :diamond, markersize = 4, markercolor = :red, label=false)
p1 = scatter!(λEMIT[1], aI_EMIT_aer5p0_900hPa_loCH4[1], marker = :square, markersize = 4, markercolor = :red, label=false)

p2 = plot(λ[2], aI_conv_aer0p0_loCH4[2], linewidth = 1, linecolor=:gray, label=false)
p2 = plot!(λ[2], aI_conv_aer0p5_900hPa_loCH4[2], linewidth = 2, linecolor=:gray, label=false)
p2 = plot!(λ[2], aI_conv_aer5p0_900hPa_loCH4[2], linewidth = 3, linecolor=:gray, label=false)
p2 = scatter!(λEMIT[2], aI_EMIT_aer0p0_loCH4[2], marker = :circle, markersize = 4, markercolor = :red, label=false)
p2 = scatter!(λEMIT[2], aI_EMIT_aer0p5_900hPa_loCH4[2], marker = :diamond, markersize = 4, markercolor = :red, label=false)
p2 = scatter!(λEMIT[2], aI_EMIT_aer5p0_900hPa_loCH4[2], marker = :square, markersize = 4, markercolor = :red, label=false)

p3 = plot(λ[3], aI_conv_aer0p0_loCH4[3], xlabel = "λ [nm]", ylabel="R [mW/m²-str-nm]", linewidth = 1, linecolor=:gray, label=false)
p3 = plot!(λ[3], aI_conv_aer0p5_900hPa_loCH4[3], linewidth = 2, linecolor=:gray, label=false)
p3 = plot!(λ[3], aI_conv_aer5p0_900hPa_loCH4[3], linewidth = 3, linecolor=:gray, label=false, legend=:bottomright)
p3 = scatter!(λEMIT[3], aI_EMIT_aer0p0_loCH4[3], marker = :circle, markersize = 4, markercolor = :red, label="AOD=0")
p3 = scatter!(λEMIT[3], aI_EMIT_aer0p5_900hPa_loCH4[3], marker = :diamond, markersize = 4, markercolor = :red, label="AOD=0.5")
p3 = scatter!(λEMIT[3], aI_EMIT_aer5p0_900hPa_loCH4[3], marker = :square, markersize = 4, markercolor = :red, label="AOD=5")

p4 = plot(λ[4], cI_conv_aer0p0_loCH4[4], xlabel = "λ [nm]", linewidth = 1, linecolor=:gray, label=false)
p4 = plot!(λ[4], cI_conv_aer0p5_900hPa_loCH4[4], linewidth = 2, linecolor=:gray, label=false)
p4 = plot!(λ[4], cI_conv_aer5p0_900hPa_loCH4[4], linewidth = 3, linecolor=:gray, label=false)
p4 = scatter!(λEMIT[4], cI_EMIT_aer0p0_loCH4[4], marker = :circle, markersize = 4, markercolor = :red, label=false)
p4 = scatter!(λEMIT[4], cI_EMIT_aer0p5_900hPa_loCH4[4], marker = :diamond, markersize = 4, markercolor = :red, label=false)
p4 = scatter!(λEMIT[4], cI_EMIT_aer5p0_900hPa_loCH4[4], marker = :square, markersize = 4, markercolor = :red, label=false)

plot(p1, p2, p3, p4, layout = l, title=["ρ=0" "ρ=0.2" "ρ=0.4" "ρ=0.2"], titlefont = font(10))
#plot(p1, p2, p3, p4, layout = l, titlefont = font(10))
savefig("/home/sanghavi/EMIT/EMIT_aerosol900hPa_loCH4.pdf")
#EMIT_aerosol900hPa_loCH4.png

l = @layout[a1 a2; a3 a4]
p1 = plot(λ[1], aI_conv_aer0p0_loCH4[1], ylabel="R [mW/m²-str-nm]", linewidth = 1, linecolor=:gray, label=false)
p1 = plot!(λ[1], aI_conv_aer0p5_500hPa_loCH4[1], linewidth = 2, linecolor=:gray, label=false)
p1 = plot!(λ[1], aI_conv_aer5p0_500hPa_loCH4[1], linewidth = 3, linecolor=:gray, label=false)
p1 = plot!(λEMIT[1], aI_EMIT_aer0p0_loCH4[1], linewidth = 0, linecolor=:white, marker = :circle, markersize = 4, markercolor = :red, label=false)
p1 = plot!(λEMIT[1], aI_EMIT_aer0p5_500hPa_loCH4[1], linewidth = 0, linecolor=:white, marker = :diamond, markersize = 4, markercolor = :red, label=false)
p1 = plot!(λEMIT[1], aI_EMIT_aer5p0_500hPa_loCH4[1], linewidth = 0, linecolor=:white, marker = :square, markersize = 4, markercolor = :red, label=false)

p2 = plot(λ[2], aI_conv_aer0p0_loCH4[2], linewidth = 1, linecolor=:gray, label=false)
p2 = plot!(λ[2], aI_conv_aer0p5_500hPa_loCH4[2], linewidth = 2, linecolor=:gray, label=false)
p2 = plot!(λ[2], aI_conv_aer5p0_500hPa_loCH4[2], linewidth = 3, linecolor=:gray, label=false)
p2 = plot!(λEMIT[2], aI_EMIT_aer0p0_loCH4[2], linewidth = 0, linecolor=:white, marker = :circle, markersize = 4, markercolor = :red, label=false)
p2 = plot!(λEMIT[2], aI_EMIT_aer0p5_500hPa_loCH4[2], linewidth = 0, linecolor=:white, marker = :diamond, markersize = 4, markercolor = :red, label=false)
p2 = plot!(λEMIT[2], aI_EMIT_aer5p0_500hPa_loCH4[2], linewidth = 0, linecolor=:white, marker = :square, markersize = 4, markercolor = :red, label=false)

p3 = plot(λ[3], aI_conv_aer0p0_loCH4[3], xlabel = "λ [nm]", ylabel="R [mW/m²-str-nm]", linewidth = 1, linecolor=:gray, label=false)
p3 = plot!(λ[3], aI_conv_aer0p5_500hPa_loCH4[3], linewidth = 2, linecolor=:gray, label=false)
p3 = plot!(λ[3], aI_conv_aer5p0_500hPa_loCH4[3], linewidth = 3, linecolor=:gray, label=false, legend=:bottomright)
p3 = plot!(λEMIT[3], aI_EMIT_aer0p0_loCH4[3], linewidth = 0, linecolor=:white, marker = :circle, markersize = 4, markercolor = :red, label="AOD=0")
p3 = plot!(λEMIT[3], aI_EMIT_aer0p5_500hPa_loCH4[3], linewidth = 0, linecolor=:white, marker = :diamond, markersize = 4, markercolor = :red, label="AOD=0.5")
p3 = plot!(λEMIT[3], aI_EMIT_aer5p0_500hPa_loCH4[3], linewidth = 0, linecolor=:white, marker = :square, markersize = 4, markercolor = :red, label="AOD=5")

p4 = plot(λ[4], cI_conv_aer0p0_loCH4[4], xlabel = "λ [nm]", linewidth = 1, linecolor=:gray, label=false)
p4 = plot!(λ[4], cI_conv_aer0p5_500hPa_loCH4[4], linewidth = 2, linecolor=:gray, label=false)
p4 = plot!(λ[4], cI_conv_aer5p0_500hPa_loCH4[4], linewidth = 3, linecolor=:gray, label=false)
p4 = plot!(λEMIT[4], cI_EMIT_aer0p0_loCH4[4], linewidth = 0, linecolor=:white, marker = :circle, markersize = 4, markercolor = :red, label=false)
p4 = plot!(λEMIT[4], cI_EMIT_aer0p5_500hPa_loCH4[4], linewidth = 0, linecolor=:white, marker = :diamond, markersize = 4, markercolor = :red, label=false)
p4 = plot!(λEMIT[4], cI_EMIT_aer5p0_500hPa_loCH4[4], linewidth = 0, linecolor=:white, marker = :square, markersize = 4, markercolor = :red, label=false)

plot(p1, p2, p3, p4, layout = l, title=["ρ=0" "ρ=0.2" "ρ=0.4" "ρ=0.2"], titlefont = font(10))
#plot(p1, p2, p3, p4, layout = l, titlefont = font(10))
savefig("/home/sanghavi/EMIT/EMIT_aerosol500hPa_loCH4.pdf")



#=
l = @layout[a1 a2; a3 a4]
p1 = plot(λ[1], bI_conv_aer0p5_500hPa_loCH4[1], linewidth = 2, linecolor=:gray, label=false, ylabel="R [mW/m²-str-nm]")
p1 = plot!(λ[1], bI_conv_aer5p0_500hPa_loCH4[1], linewidth = 3, linecolor=:gray, label=false)
p1 = plot!(λEMIT[1], bI_EMIT_aer0p5_500hPa_loCH4[1], linewidth = 0, linecolor=:white, marker = :diamond, markersize = 4, markercolor = :red, label=false)
p1 = plot!(λEMIT[1], bI_EMIT_aer5p0_500hPa_loCH4[1], linewidth = 0, linecolor=:white, marker = :square, markersize = 4, markercolor = :red, label=false)

p1 = plot!(λ[1], bI_conv_aer0p5_900hPa_loCH4[1], linewidth = 2, linecolor=:black, label=false)
p1 = plot!(λ[1], bI_conv_aer5p0_900hPa_loCH4[1], linewidth = 3, linecolor=:black, label=false)
p1 = plot!(λEMIT[1], bI_EMIT_aer0p5_900hPa_loCH4[1], linewidth = 0, linecolor=:white, marker = :diamond, markersize = 4, markercolor = :blue, label=false)
p1 = plot!(λEMIT[1], bI_EMIT_aer5p0_900hPa_loCH4[1], linewidth = 0, linecolor=:white, marker = :square, markersize = 4, markercolor = :blue, label=false)


p2 = plot(λ[2], bI_conv_aer0p5_500hPa_loCH4[2], linewidth = 2, linecolor=:gray, label=false, alpha=0.7)
p2 = plot!(λ[2], bI_conv_aer5p0_500hPa_loCH4[2], linewidth = 3, linecolor=:gray, label=false, alpha=0.7)
p2 = plot!(λEMIT[2], bI_EMIT_aer0p5_500hPa_loCH4[2], linewidth = 0, linecolor=:white, marker = :diamond, alpha=0.7, markersize = 4, markercolor = :red, label=false)
p2 = plot!(λEMIT[2], bI_EMIT_aer5p0_500hPa_loCH4[2], linewidth = 0, linecolor=:white, marker = :square, alpha=0.7, markersize = 4, markercolor = :red, label=false)

p2 = plot!(λ[2], bI_conv_aer0p5_900hPa_loCH4[2], linewidth = 2, linecolor=:black, label=false, alpha=0.7)
p2 = plot!(λ[2], bI_conv_aer5p0_900hPa_loCH4[2], linewidth = 3, linecolor=:black, label=false, alpha=0.7)
p2 = plot!(λEMIT[2], bI_EMIT_aer0p5_900hPa_loCH4[2], linewidth = 0, linecolor=:white, marker = :diamond, alpha=0.7, markersize = 4, markercolor = :blue, label=false)
p2 = plot!(λEMIT[2], bI_EMIT_aer5p0_900hPa_loCH4[2], linewidth = 0, linecolor=:white, marker = :square, alpha=0.7, markersize = 4, markercolor = :blue, label=false)
=#
l = @layout[a1 a2; a3 a4]
p3a = plot(λ[3], bI_conv_aer0p5_500hPa_loCH4[3], linewidth = 2, linecolor=:gray, label=false, alpha=0.7, ylabel="R [mW/m²-str-nm]")
p3a = scatter!(λEMIT[3], bI_EMIT_aer0p5_500hPa_loCH4[3], marker = :diamond, alpha=0.7, markersize = 4, markercolor = :red, label=false)
p3a = plot!(λ[3], bI_conv_aer0p5_900hPa_loCH4[3], linewidth = 2, linecolor=:black, label=false, alpha=0.7)
p3a = scatter!(λEMIT[3], bI_EMIT_aer0p5_900hPa_loCH4[3], marker = :diamond, alpha=0.7, markersize = 4, markercolor = :blue, label=false)

p3b = plot(λ[3], bI_conv_aer5p0_500hPa_loCH4[3], linewidth = 3, linecolor=:gray, label=false, alpha=0.7, xlabel = "λ [nm]", ylabel="R [mW/m²-str-nm]", legend=:bottomright)
p3b = scatter!(λEMIT[3], bI_EMIT_aer5p0_500hPa_loCH4[3], marker = :square, alpha=0.7, markersize = 4, markercolor = :red, label="z₀=5000 m (500 hPa)", legendfontsize=6)
p3b = plot!(λ[3], bI_conv_aer5p0_900hPa_loCH4[3], linewidth = 3, linecolor=:black, label=false, alpha=0.7)
p3b = scatter!(λEMIT[3], bI_EMIT_aer5p0_900hPa_loCH4[3], marker = :square, alpha=0.7, markersize = 4, markercolor = :blue, label="z₀=850 m (900 hPa)", legendfontsize=6)

p4a = plot(λ[4], bI_conv_aer0p5_500hPa_loCH4[4], linewidth = 2, linecolor=:gray, label=false, alpha=0.7)
p4a = scatter!(λEMIT[4], bI_EMIT_aer0p5_500hPa_loCH4[4], marker = :diamond, alpha=0.7, markersize = 4, markercolor = :red, label=false)
p4a = plot!(λ[4], bI_conv_aer0p5_900hPa_loCH4[4], linewidth = 2, linecolor=:black, label=false, alpha=0.7)
p4a = scatter!(λEMIT[4], bI_EMIT_aer0p5_900hPa_loCH4[4], marker = :diamond, alpha=0.7, markersize = 4, markercolor = :blue, label=false)

p4b = plot(λ[4], bI_conv_aer5p0_500hPa_loCH4[4], linewidth = 3, linecolor=:gray, label=false, alpha=0.7, xlabel = "λ [nm]")
p4b = scatter!(λEMIT[4], bI_EMIT_aer5p0_500hPa_loCH4[4], marker = :square, alpha=0.7, markersize = 4, markercolor = :red, label=false)
p4b = plot!(λ[4], bI_conv_aer5p0_900hPa_loCH4[4], linewidth = 3, linecolor=:black, label=false, alpha=0.7)
p4b = scatter!(λEMIT[4], bI_EMIT_aer5p0_900hPa_loCH4[4], marker = :square, alpha=0.7, markersize = 4, markercolor = :blue, label=false)

plot(p3a, p4a, p3b, p4b, layout = l, title=["ρ=0, AOD=0.5" "ρ=0, AOD=0.5" "ρ=0, AOD=5" "ρ=0, AOD=5"], titlefont = font(10))
#plot(p1, p2, p3, p4, layout = l, titlefont = font(10))
savefig("/home/sanghavi/EMIT/EMIT_loCH4_ocean_aerht.pdf")
#EMIT_loCH4_ocean_aerht.png

l = @layout[a1 a2 a3; a4 a5 a6]
p1 = plot(λ[5], cI_conv_aer0p0_loCH4[5], linewidth = 2, linecolor=:gray, label=false, alpha=0.7, ylabel="R [mW/m²-str-nm]")
p1 = scatter!(λEMIT[5], cI_EMIT_aer0p0_loCH4[5], marker = :circle, alpha=0.7, markersize = 4, markercolor = :red, label=false)
p1 = plot!(λ[5], cI_conv_aer0p0_hiCH4[5], linewidth = 2, linecolor=:black, label=false, alpha=0.7)
p1 = scatter!(λEMIT[5], cI_EMIT_aer0p0_hiCH4[5], marker = :circle, alpha=0.7, markersize = 4, markercolor = :green, label=false)

p2 = plot(λ[5], cI_conv_aer0p5_900hPa_loCH4[5], linewidth = 2, linecolor=:gray, label=false, alpha=0.7)
p2 = scatter!(λEMIT[5], cI_EMIT_aer0p5_900hPa_loCH4[5], marker = :diamond, alpha=0.7, markersize = 4, markercolor = :red, label=false)
p2 = plot!(λ[5], cI_conv_aer0p5_900hPa_hiCH4[5], linewidth = 2, linecolor=:black, label=false, alpha=0.7)
p2 = scatter!(λEMIT[5], cI_EMIT_aer0p5_900hPa_hiCH4[5], marker = :diamond, alpha=0.7, markersize = 4, markercolor = :green, label=false)

p3 = plot(λ[5], cI_conv_aer5p0_900hPa_loCH4[5], linewidth = 3, linecolor=:gray, label=false, alpha=0.7, legend=:bottomleft)
p3 = scatter!(λEMIT[5], cI_EMIT_aer5p0_900hPa_loCH4[5], marker = :square, alpha=0.7, markersize = 4, markercolor = :red, label="low CH₄", legendfontsize=6)
p3 = plot!(λ[5], cI_conv_aer5p0_900hPa_loCH4[5], linewidth = 3, linecolor=:black, label=false, alpha=0.7)
p3 = scatter!(λEMIT[5], cI_EMIT_aer5p0_900hPa_loCH4[5], marker = :square, alpha=0.7, markersize = 4, markercolor = :green, label="high CH₄", legendfontsize=6)

p4 = plot(λ[5], bI_conv_aer0p0_loCH4[5], linewidth = 2, linecolor=:gray, label=false, alpha=0.7, ylabel="R [mW/m²-str-nm]", xlabel = "λ [nm]")
p4 = scatter!(λEMIT[5], bI_EMIT_aer0p0_loCH4[5], marker = :circle, alpha=0.7, markersize = 4, markercolor = :red, label=false)
p4 = plot!(λ[5], bI_conv_aer0p0_hiCH4[5], linewidth = 2, linecolor=:black, label=false, alpha=0.7)
p4 = scatter!(λEMIT[5], bI_EMIT_aer0p0_hiCH4[5], marker = :circle, alpha=0.7, markersize = 4, markercolor = :green, label=false)

p5 = plot(λ[5], bI_conv_aer0p5_900hPa_loCH4[5], linewidth = 2, linecolor=:gray, label=false, alpha=0.7, xlabel = "λ [nm]")
p5 = scatter!(λEMIT[5], bI_EMIT_aer0p5_900hPa_loCH4[5], marker = :diamond, alpha=0.7, markersize = 4, markercolor = :red, label=false)
p5 = plot!(λ[5], bI_conv_aer0p5_900hPa_hiCH4[5], linewidth = 2, linecolor=:black, label=false, alpha=0.7)
p5 = scatter!(λEMIT[5], bI_EMIT_aer0p5_900hPa_hiCH4[5], marker = :diamond, alpha=0.7, markersize = 4, markercolor = :green, label=false)

p6 = plot(λ[5], bI_conv_aer5p0_900hPa_loCH4[5], linewidth = 3, linecolor=:gray, label=false, alpha=0.7, xlabel = "λ [nm]")
p6 = scatter!(λEMIT[5], bI_EMIT_aer5p0_900hPa_loCH4[5], marker = :square, alpha=0.7, markersize = 4, markercolor = :red, label=false, legendfontsize=6)
p6 = plot!(λ[5], bI_conv_aer5p0_900hPa_loCH4[5], linewidth = 3, linecolor=:black, label=false, alpha=0.7)
p6 = scatter!(λEMIT[5], bI_EMIT_aer5p0_900hPa_loCH4[5], marker = :square, alpha=0.7, markersize = 4, markercolor = :green, label=false, legendfontsize=6)
#p6 = scatter!(λEMIT[5], bI_EMIT_aer5p0_900hPa_loCH4[5], linewidth = 0, linecolor=:white, marker = :square, alpha=0.7, markersize = 4, markercolor = :green, label=false, legendfontsize=6)

plot(p1, p2, p3, p4, p5, p6, layout = l, title=["ρ=0.2, AOD=0" "ρ=0.2, AOD=0.5" "ρ=0.2, AOD=5" "ρ=0, AOD=0" "ρ=0, AOD=0.5" "ρ=0, AOD=5"], titlefont = font(10))
#plot(p1, p2, p3, p4, layout = l, titlefont = font(10))
savefig("/home/sanghavi/EMIT/EMIT_CH4aer.pdf")
#EMIT_CH4aer.png




l = @layout[a1 a2; a3 a4]
p1 = plot(λ[1], bI_conv_aer0p0_loCH4[1], ylabel="R [mW/m²-str-nm]", linewidth = 1, linecolor=:gray, label=false)
p1 = plot!(λ[1], bI_conv_aer0p5_500hPa_loCH4[1], linewidth = 2, linecolor=:gray, label=false)
p1 = plot!(λ[1], bI_conv_aer5p0_500hPa_loCH4[1], linewidth = 3, linecolor=:gray, label=false)
p1 = plot!(λ[1], bI_conv_aer0p5_900hPa_loCH4[1], linewidth = 2, linecolor=:black, label=false)
p1 = plot!(λ[1], bI_conv_aer5p0_900hPa_loCH4[1], linewidth = 3, linecolor=:black, label=false)

p1 = plot!(λEMIT[1], I_EMIT_aer0p0_loCH4[1], linewidth = 0, linecolor=:white, marker = :circle, markersize = 4, markercolor = :red, label=false)
p1 = plot!(λEMIT[1], I_EMIT_aer0p5_500hPa_loCH4[1], linewidth = 0, linecolor=:white, marker = :diamond, markersize = 4, markercolor = :red, label=false)
p1 = plot!(λEMIT[1], I_EMIT_aer5p0_500hPa_loCH4[1], linewidth = 0, linecolor=:white, marker = :square, markersize = 4, markercolor = :red, label=false)
p1 = plot!(λEMIT[1], I_EMIT_aer0p5_900hPa_loCH4[1], linewidth = 0, linecolor=:white, marker = :diamond, markersize = 4, markercolor = :blue, label=false)
p1 = plot!(λEMIT[1], I_EMIT_aer5p0_900hPa_loCH4[1], linewidth = 0, linecolor=:white, marker = :square, markersize = 4, markercolor = :blue, label=false)

p2 = plot(λ[2], I_conv_aer0p0_loCH4[2], linewidth = 1, linecolor=:gray, label=false)
p2 = plot!(λ[2], I_conv_aer0p5_500hPa_loCH4[2], linewidth = 2, linecolor=:gray, label=false)
p2 = plot!(λ[2], I_conv_aer5p0_500hPa_loCH4[2], linewidth = 3, linecolor=:gray, label=false)
p2 = plot!(λEMIT[2], I_EMIT_aer0p0_loCH4[2], linewidth = 0, linecolor=:white, marker = :circle, markersize = 4, markercolor = :red, label=false)
p2 = plot!(λEMIT[2], I_EMIT_aer0p5_500hPa_loCH4[2], linewidth = 0, linecolor=:white, marker = :diamond, markersize = 4, markercolor = :red, label=false)
p2 = plot!(λEMIT[2], I_EMIT_aer5p0_500hPa_loCH4[2], linewidth = 0, linecolor=:white, marker = :square, markersize = 4, markercolor = :red, label=false)

p3 = plot(λ[3], I_conv_aer0p0_loCH4[3], xlabel = "λ [nm]", ylabel="R [mW/m²-str-nm]", linewidth = 1, linecolor=:gray, label=false)
p3 = plot!(λ[3], I_conv_aer0p5_500hPa_loCH4[3], linewidth = 2, linecolor=:gray, label=false)
p3 = plot!(λ[3], I_conv_aer5p0_500hPa_loCH4[3], linewidth = 3, linecolor=:gray, label=false, legend=:bottomright)
p3 = plot!(λEMIT[3], I_EMIT_aer0p0_loCH4[3], linewidth = 0, linecolor=:white, marker = :circle, markersize = 4, markercolor = :red, label="AOD=0")
p3 = plot!(λEMIT[3], I_EMIT_aer0p5_500hPa_loCH4[3], linewidth = 0, linecolor=:white, marker = :diamond, markersize = 4, markercolor = :red, label="AOD=0.5")
p3 = plot!(λEMIT[3], I_EMIT_aer5p0_500hPa_loCH4[3], linewidth = 0, linecolor=:white, marker = :square, markersize = 4, markercolor = :red, label="AOD=5")

p4 = plot(λ[4], I_conv_aer0p0_loCH4[4], xlabel = "λ [nm]", linewidth = 1, linecolor=:gray, label=false)
p4 = plot!(λ[4], I_conv_aer0p5_500hPa_loCH4[4], linewidth = 2, linecolor=:gray, label=false)
p4 = plot!(λ[4], I_conv_aer5p0_500hPa_loCH4[4], linewidth = 3, linecolor=:gray, label=false)
p4 = plot!(λEMIT[4], I_EMIT_aer0p0_loCH4[4], linewidth = 0, linecolor=:white, marker = :circle, markersize = 4, markercolor = :red, label=false)
p4 = plot!(λEMIT[4], I_EMIT_aer0p5_500hPa_loCH4[4], linewidth = 0, linecolor=:white, marker = :diamond, markersize = 4, markercolor = :red, label=false)
p4 = plot!(λEMIT[4], I_EMIT_aer5p0_500hPa_loCH4[4], linewidth = 0, linecolor=:white, marker = :square, markersize = 4, markercolor = :red, label=false)

plot(p1, p2, p3, p4, layout = l, titlefont = font(10))
savefig("/home/sanghavi/EMIT/EMIT_aerosol500hPa_loCH4.pdf")

#Q_conv_noRS = imfilter(RnoRS[1,2,:], kernel)
#Q_conv = imfilter(R[1,2,:], kernel)
#ieQ_conv = imfilter(ieR[1,2,:], kernel)

#convfct = 1e7./ν.^2   # to convert radiance units from mW/m²-str-cm⁻¹ to mW/m²-str-nm


l = @layout [a1 a2 ; b1 b2; c1 c2]
p1 = plot(1e7./ν, RnoRS[1,1,:].*convfct, linecolor=:black, xlims=[385,415])
p1 = plot!(1e7./ν, I_conv_noRS.*convfct, linewidth = 3, linecolor=:red, ylabel="[mW/m²/str/nm]", yguidefontsize=8)
p2 = plot(1e7./ν, (R[1,1,:].-RnoRS[1,1,:].+ieR[1,1,:]).*convfct, linecolor=:black, xlabel = "λ [nm]", xlims=[385,415])
p2 = plot!(1e7./ν, (I_conv.-I_conv_noRS.+ieI_conv).*convfct, linewidth = 3, linecolor=:red, ylabel="[mW/m²/str/nm]", yguidefontsize=8)
p3 = plot(1e7./ν, 100*(R[1,1,:].+ieR[1,1,:].-RnoRS[1,1,:])./RnoRS[1,1,:], linecolor=:black, ylabel="[%]", xlims=[385,415], ylims=[-5,150])
p3 = plot!(1e7./ν, 100*(I_conv.+ieI_conv.-I_conv_noRS)./I_conv_noRS, linewidth = 2,  linecolor=:red)


q1 = plot(1e7./ν, RnoRS[1,2,:].*convfct, linecolor=:black, xlims=[385,415])
q1 = plot!(1e7./ν, Q_conv_noRS.*convfct, linewidth = 3, linecolor=:red)
q2 = plot(1e7./ν, (R[1,2,:].-RnoRS[1,2,:].+ieR[1,2,:]).*convfct, linecolor=:black, xlabel = "λ [nm]", xlims=[385,415])
q2 = plot!(1e7./ν, (Q_conv.-Q_conv_noRS.+ieQ_conv).*convfct, linewidth = 3, linecolor=:red)
q3 = plot(1e7./ν, 100*(R[1,2,:].+ieR[1,2,:].-RnoRS[1,2,:])./RnoRS[1,2,:], linecolor=:black, xlims=[385,415])
q3 = plot!(1e7./ν, 100*(Q_conv.+ieQ_conv.-Q_conv_noRS)./Q_conv_noRS, linewidth = 2, linecolor=:red, ylims=[-0.6,1.5])
plot(p1, q1, p2, q2, p3, q3, layout = l, link=:x, legend = false, title = ["I₀ (no RS)" "Q₀ (no RS)" "I₁-I₀" "Q₁-Q₀" "(I₁-I₀)/I₀ [%]" "(Q₁-Q₀)/Q₀ [%]"], titlefont = font(10))
savefig("RingEffect.png")

l = @layout [a1 a2]
p1 = plot(1e7./ν, 100*ieR[1,1,:]./R[1,1,:], linecolor=:black, xlims=[385,415], ylims=[-5,150])
p1 = plot!(1e7./ν, 100*ieI_conv./I_conv, linewidth = 2, linecolor=:red, xlabel = "λ [nm]")
q1 = plot(1e7./ν, 100*ieR[1,2,:]./R[1,2,:], linecolor=:black, xlims=[385,415], ylims=[-0.1,2])
q1 = plot!(1e7./ν, 100*ieQ_conv./Q_conv, linewidth = 2, linecolor=:red, xlabel = "λ [nm]")
plot(p1, q1, layout = l, legend = false, title = ["Iᵢ/Iₑ [%]" "Qᵢ/Qₑ [%]"], titlefont = font(10))
savefig("RingSpectrum.png")

using DelimitedFiles
raylout_rrs = readdlm("out_ray_rrs_385nm_405nm.dat")
raylout_nors = readdlm("out_ray_nors_385nm_405nm.dat")
raylout_rrs_ss = readdlm("out_ray_rrs_ss_385nm_405nm.dat")
aerout_rrs = readdlm("out_aer_rrs_385nm_405nm.dat")
aerout_nors = readdlm("out_aer_nors_385nm_405nm.dat")
aerout_rrs_ss = readdlm("out_aer_rrs_ss_385nm_405nm.dat")

x = -40:0.6:40
#kernel = InstrumentOperator.create_instrument_kernel(Normal(0, 12.5), x) #defining a Gaussian kernel for convulution in wavenumber space
kernel = InstrumentOperator.create_instrument_kernel(Normal(0, 12.5), x)
#I_conv = InstrumentOperator.conv_spectra(kernel, )

ν = raylout_nors[:,1]
#Rayleigh
rI_noRS = zeros(length(ν))
rQ_noRS = zeros(length(ν))
rI = zeros(length(ν))
rQ = zeros(length(ν))
ierI = zeros(length(ν))
ierQ = zeros(length(ν))
rI_ss = zeros(length(ν))
rQ_ss = zeros(length(ν))
ierI_ss = zeros(length(ν))
ierQ_ss = zeros(length(ν))

rI_conv_noRS = zeros(length(ν))
rQ_conv_noRS = zeros(length(ν))
rI_conv = zeros(length(ν))
rQ_conv = zeros(length(ν))
ierI_conv = zeros(length(ν))
ierQ_conv = zeros(length(ν))
rI_conv_ss = zeros(length(ν))
rQ_conv_ss = zeros(length(ν))
ierI_conv_ss = zeros(length(ν))
ierQ_conv_ss = zeros(length(ν))

rI_noRS = raylout_nors[:,2]
rQ_noRS = raylout_nors[:,3]
rI = raylout_rrs[:,2]
rQ = raylout_rrs[:,3]
ierI = raylout_rrs[:,4]
ierQ = raylout_rrs[:,5]
rI_ss = raylout_rrs_ss[:,2]
rQ_ss = raylout_rrs_ss[:,3]
ierI_ss = raylout_rrs_ss[:,4]
ierQ_ss = raylout_rrs_ss[:,5]

rI_conv_noRS = imfilter(rI_noRS, kernel)
rQ_conv_noRS = imfilter(rQ_noRS, kernel)
rI_conv = imfilter(rI, kernel)
rQ_conv = imfilter(rQ, kernel)
ierI_conv = imfilter(ierI, kernel)
ierQ_conv = imfilter(ierQ, kernel)
rI_conv_ss = imfilter(rI_ss, kernel)
rQ_conv_ss = imfilter(rQ_ss, kernel)
ierI_conv_ss = imfilter(ierI_ss, kernel)
ierQ_conv_ss = imfilter(ierQ_ss, kernel)

#Aerosol
aI_noRS = zeros(length(ν))
aQ_noRS = zeros(length(ν))
aI = zeros(length(ν))
aQ = zeros(length(ν))
ieaI = zeros(length(ν))
ieaQ = zeros(length(ν))
aI_ss = zeros(length(ν))
aQ_ss = zeros(length(ν))
ieaI_ss = zeros(length(ν))
ieaQ_ss = zeros(length(ν))

aI_conv_noRS = zeros(length(ν))
aQ_conv_noRS = zeros(length(ν))
aI_conv = zeros(length(ν))
aQ_conv = zeros(length(ν))
ieaI_conv = zeros(length(ν))
ieaQ_conv = zeros(length(ν))
aI_conv_ss = zeros(length(ν))
aQ_conv_ss = zeros(length(ν))
ieaI_conv_ss = zeros(length(ν))
ieaQ_conv_ss = zeros(length(ν))

aI_noRS = aerout_nors[:,2]
aQ_noRS = aerout_nors[:,3]
aI = aerout_rrs[:,2]
aQ = aerout_rrs[:,3]
ieaI = aerout_rrs[:,4]
ieaQ = aerout_rrs[:,5]
aI_ss = aerout_rrs_ss[:,2]
aQ_ss = aerout_rrs_ss[:,3]
ieaI_ss = aerout_rrs_ss[:,4]
ieaQ_ss = aerout_rrs_ss[:,5]

aI_conv_noRS = imfilter(aI_noRS, kernel)
aQ_conv_noRS = imfilter(aQ_noRS, kernel)
aI_conv = imfilter(aI, kernel)
aQ_conv = imfilter(aQ, kernel)
ieaI_conv = imfilter(ieaI, kernel)
ieaQ_conv = imfilter(ieaQ, kernel)
aI_conv_ss = imfilter(aI_ss, kernel)
aQ_conv_ss = imfilter(aQ_ss, kernel)
ieaI_conv_ss = imfilter(ieaI_ss, kernel)
ieaQ_conv_ss = imfilter(ieaQ_ss, kernel)

convfct = 1e7./ν.^2   # to convert radiance units from mW/m²-str-cm⁻¹ to mW/m²-str-nm

p1 = plot(1e7./ν, rQ_noRS./rI_noRS, linewidth = 3, linecolor=:black, label="Rayleigh", ylims=[0.2,0.5])
p1 = plot!(1e7./ν, (rQ+ierQ)./(rI+ierI), linecolor=:black, label="")
p1 = plot!(1e7./ν, (rQ_conv+ierQ_conv)./(rI_conv+ierI_conv), linewidth = 3, linecolor=:red, label="")
p1 = plot!(1e7./ν, aQ_noRS./aI_noRS, linewidth = 3, linecolor=:grey, label="Aerosol")
p1 = plot!(1e7./ν, (aQ+ieaQ)./(aI+ieaI), linecolor=:grey, label="")
p1 = plot!(1e7./ν, (aQ_conv+ieaQ_conv)./(aI_conv+ieaI_conv), linewidth = 3, linecolor=:red, xlabel = "λ [nm]", ylabel = "DOLP", label="")
plot(p1, title= "Zenith obs, θ₀=63ᵒ, A=0.05")#, legend = false)
savefig("compDOLP.png")


# New
#=
rayl10_rrs_ss = [ν R10_ss[1,1,:] R10_ss[1,2,:] R10_ss[1,3,:] ieR10_ss[1,1,:] ieR10_ss[1,2,:] ieR10_ss[1,3,:]]
rayl70_rrs_ss = [ν R70_ss[1,1,:] R70_ss[1,2,:] R70_ss[1,3,:] ieR70_ss[1,1,:] ieR70_ss[1,2,:] ieR70_ss[1,3,:]]
rayl10_rrs = [ν R10[1,1,:] R10[1,2,:] R10[1,3,:] ieR10[1,1,:] ieR10[1,2,:] ieR10[1,3,:]]
rayl70_rrs = [ν R70[1,1,:] R70[1,2,:] R70[1,3,:] ieR70[1,1,:] ieR70[1,2,:] ieR70[1,3,:]]
rayl10_nors_ss = [ν R10noRS_ss[1,1,:] R10noRS_ss[1,2,:] R10noRS_ss[1,3,:]]
rayl70_nors_ss = [ν R70noRS_ss[1,1,:] R70noRS_ss[1,2,:] R70noRS_ss[1,3,:]]
rayl10_nors = [ν R10noRS[1,1,:] R10noRS[1,2,:] R10noRS[1,3,:]]
rayl70_nors = [ν R70noRS[1,1,:] R70noRS[1,2,:] R70noRS[1,3,:]]

writedlm("out_ray10_rrs_ss_385nm_405nm.dat", rayl10_rrs_ss)
writedlm("out_ray70_rrs_ss_385nm_405nm.dat", rayl70_rrs_ss)
writedlm("out_ray70_rrs_385nm_405nm.dat", rayl70_rrs)
writedlm("out_ray10_rrs_385nm_405nm.dat", rayl10_rrs)
writedlm("out_ray10_nors_385nm_405nm.dat", rayl10_nors)
writedlm("out_ray70_nors_385nm_405nm.dat", rayl70_nors)
writedlm("out_ray70_nors_ss_385nm_405nm.dat", rayl70_nors_ss)
writedlm("out_ray10_nors_ss_385nm_405nm.dat", rayl10_nors_ss)
=#
rayl10_rrs = readdlm("out_ray10_rrs_385nm_405nm.dat")
rayl10_nors = readdlm("out_ray10_nors_385nm_405nm.dat")
rayl10_rrs_ss = readdlm("out_ray10_rrs_ss_385nm_405nm.dat")
rayl10_nors_ss = readdlm("out_ray10_nors_ss_385nm_405nm.dat")
rayl70_rrs = readdlm("out_ray70_rrs_385nm_405nm.dat")
rayl70_nors = readdlm("out_ray70_nors_385nm_405nm.dat")
rayl70_rrs_ss = readdlm("out_ray70_rrs_ss_385nm_405nm.dat")
rayl70_nors_ss = readdlm("out_ray70_nors_ss_385nm_405nm.dat")

ν = rayl10_nors[:,1]

I10 = rayl10_rrs[:,2]
Q10 = rayl10_rrs[:,3]
U10 = rayl10_rrs[:,4]
ieI10 = rayl10_rrs[:,5]
ieQ10 = rayl10_rrs[:,6]
ieU10 = rayl10_rrs[:,7]

I10_ss = rayl10_rrs_ss[:,2]
Q10_ss = rayl10_rrs_ss[:,3]
U10_ss = rayl10_rrs_ss[:,4]
ieI10_ss = rayl10_rrs_ss[:,5]
ieQ10_ss = rayl10_rrs_ss[:,6]
ieU10_ss = rayl10_rrs_ss[:,7]

I10_nors = rayl10_nors[:,2]
Q10_nors = rayl10_nors[:,3]
U10_nors = rayl10_nors[:,4]

I10_nors_ss = rayl10_nors_ss[:,2]
Q10_nors_ss = rayl10_nors_ss[:,3]
U10_nors_ss = rayl10_nors_ss[:,4]

I70 = rayl70_rrs[:,2]
Q70 = rayl70_rrs[:,3]
U70 = rayl70_rrs[:,4]
ieI70 = rayl70_rrs[:,5]
ieQ70 = rayl70_rrs[:,6]
ieU70 = rayl70_rrs[:,7]

I70_ss = rayl70_rrs_ss[:,2]
Q70_ss = rayl70_rrs_ss[:,3]
U70_ss = rayl70_rrs_ss[:,4]
ieI70_ss = rayl70_rrs_ss[:,5]
ieQ70_ss = rayl70_rrs_ss[:,6]
ieU70_ss = rayl70_rrs_ss[:,7]

I70_nors = rayl70_nors[:,2]
Q70_nors = rayl70_nors[:,3]
U70_nors = rayl70_nors[:,4]

I70_nors_ss = rayl70_nors_ss[:,2]
Q70_nors_ss = rayl70_nors_ss[:,3]
U70_nors_ss = rayl70_nors_ss[:,4]

x = -40:0.3:40
#kernel = InstrumentOperator.create_instrument_kernel(Normal(0, 12.5), x) #defining a Gaussian kernel for convulution in wavenumber space
#kernel = InstrumentOperator.create_instrument_kernel(Normal(0, 6.25), x)
kernel = InstrumentOperator.create_instrument_kernel(Normal(0, 7.64), x)

#I_conv = InstrumentOperator.conv_spectra(kernel, )
I10_conv_noRS = imfilter(I10_nors, kernel)
I10_conv_noRS_ss = imfilter(I10_nors_ss, kernel)
I10_conv = imfilter(I10, kernel)
ieI10_conv = imfilter(ieI10, kernel)
I10_ss_conv = imfilter(I10_ss, kernel)
ieI10_ss_conv = imfilter(ieI10_ss, kernel)

Q10_conv_noRS = imfilter(Q10_nors, kernel)
Q10_conv_noRS_ss = imfilter(Q10_nors_ss, kernel)
Q10_conv = imfilter(Q10, kernel)
ieQ10_conv = imfilter(ieQ10, kernel)
Q10_ss_conv = imfilter(Q10_ss, kernel)
ieQ10_ss_conv = imfilter(ieQ10_ss, kernel)

U10_conv_noRS = imfilter(U10_nors, kernel)
U10_conv_noRS_ss = imfilter(U10_nors_ss, kernel)
U10_conv = imfilter(U10, kernel)
ieU10_conv = imfilter(ieU10, kernel)
U10_ss_conv = imfilter(U10_ss, kernel)
ieU10_ss_conv = imfilter(ieU10_ss, kernel)

I70_conv_noRS = imfilter(I70_nors, kernel)
I70_conv_noRS_ss = imfilter(I70_nors_ss, kernel)
I70_conv = imfilter(I70, kernel)
ieI70_conv = imfilter(ieI70, kernel)
I70_ss_conv = imfilter(I70_ss, kernel)
ieI70_ss_conv = imfilter(ieI70_ss, kernel)

Q70_conv_noRS = imfilter(Q70_nors, kernel)
Q70_conv_noRS_ss = imfilter(Q70_nors_ss, kernel)
Q70_conv = imfilter(Q70, kernel)
ieQ70_conv = imfilter(ieQ70, kernel)
Q70_ss_conv = imfilter(Q70_ss, kernel)
ieQ70_ss_conv = imfilter(ieQ70_ss, kernel)

U70_conv_noRS = imfilter(U70_nors, kernel)
U70_conv_noRS_ss = imfilter(U70_nors_ss, kernel)
U70_conv = imfilter(U70, kernel)
ieU70_conv = imfilter(ieU70, kernel)
U70_ss_conv = imfilter(U70_ss, kernel)
ieU70_ss_conv = imfilter(ieU70_ss, kernel)

convfct = 1e7./ν.^2
l = @layout [a1 a2 ; b1 b2; c1 c2]

#p1 = plot(1e7./ν, 100*((I10.+ieI10)./I10_nors.-1), linecolor=:black, ylabel="R₁ [%]", alpha=0.3)
p1 = plot(1e7./ν, 100*((I10_conv.+ieI10_conv)./I10_conv_noRS.-1), linewidth = 3, linecolor=:black)
#p1 = plot!(1e7./ν, 100*((I10_ss.+ieI10_ss)./I10_nors_ss.-1), linecolor=:red, alpha=0.3)
p1 = plot!(1e7./ν, 100*((I10_ss_conv.+ieI10_ss_conv)./I10_conv_noRS_ss.-1), linewidth = 3, linecolor=:red)

#p2 = plot(1e7./ν, 100*((Q10.+ieQ10)./Q10_nors.-1), linecolor=:black, ylabel="R₂ [%]", alpha=0.3)
p2 = plot(1e7./ν, 100*((Q10_conv.+ieQ10_conv)./Q10_conv_noRS.-1), linewidth = 3, linecolor=:black)
#p2 = plot!(1e7./ν, 100*((Q10_ss.+ieQ10_ss)./Q10_nors_ss.-1), linecolor=:red, alpha=0.3)
p2 = plot!(1e7./ν, 100*((Q10_ss_conv.+ieQ10_ss_conv)./Q10_conv_noRS_ss.-1), linewidth = 3, linecolor=:red)

#p3 = plot(1e7./ν, 100*((U10.+ieU10)./U10_nors.-1), linecolor=:black, ylabel="R₃ [%]", alpha=0.3)
p3 = plot(1e7./ν, 100*((U10_conv.+ieU10_conv)./U10_conv_noRS.-1), linewidth = 3, linecolor=:black)
#p3 = plot!(1e7./ν, 100*((U10_ss.+ieU10_ss)./U10_nors_ss.-1), linecolor=:red, alpha=0.3)
p3 = plot!(1e7./ν, 100*((U10_ss_conv.+ieU10_ss_conv)./U10_conv_noRS_ss.-1), linewidth = 3, linecolor=:red)

#q1 = plot(1e7./ν, 100*((I70.+ieI70)./I70_nors.-1), linecolor=:black, ylabel="R₁ [%]", alpha=0.3)
q1 = plot(1e7./ν, 100*((I70_conv.+ieI70_conv)./I70_conv_noRS.-1), linewidth = 3, linecolor=:black)
#q1 = plot!(1e7./ν, 100*((I70_ss.+ieI70_ss)./I70_nors_ss.-1), linecolor=:red, alpha=0.3)
q1 = plot!(1e7./ν, 100*((I70_ss_conv.+ieI70_ss_conv)./I70_conv_noRS_ss.-1), linewidth = 3, linecolor=:red)

#q2 = plot(1e7./ν, 100*((Q70.+ieQ70)./Q70_nors.-1), linecolor=:black, ylabel="R₂ [%]", alpha=0.3)
q2 = plot(1e7./ν, 100*((Q70_conv.+ieQ70_conv)./Q70_conv_noRS.-1), linewidth = 3, linecolor=:black)
#q2 = plot!(1e7./ν, 100*((Q70_ss.+ieQ70_ss)./Q70_nors_ss.-1), linecolor=:red, alpha=0.3)
q2 = plot!(1e7./ν, 100*((Q70_ss_conv.+ieQ70_ss_conv)./Q70_conv_noRS_ss.-1), linewidth = 3, linecolor=:red)

#q3 = plot(1e7./ν, 100*((U70.+ieU70)./U70_nors.-1), linecolor=:black, ylabel="R₃ [%]", alpha=0.3)
q3 = plot(1e7./ν, 100*((U70_conv.+ieU70_conv)./U70_conv_noRS.-1), linewidth = 3, linecolor=:black)
#q3 = plot!(1e7./ν, 100*((U70_ss.+ieU70_ss)./U70_nors_ss.-1), linecolor=:red, alpha=0.3)
q3 = plot!(1e7./ν, 100*((U70_ss_conv.+ieU70_ss_conv)./U70_conv_noRS_ss.-1), linewidth = 3, linecolor=:red)

plot(p1, q1, p2, q2, p3, q3, layout = l, link=:y, legend = false, title = ["SZA=10ᵒ" "SZA=70ᵒ" "" "" "" ""], titlefont = font(10))
savefig("abc.png")

l = @layout [a1 a2 ; b1 b2; c1 c2]

#p1 = plot(1e7./ν, 100*((I10.+ieI10)./I10_nors.-1), linecolor=:black, ylabel="R₁ [%]", alpha=0.3)
p1 = plot(1e7./ν, 100*((I10_conv.+ieI10_conv)./I10_conv_noRS.-1), linewidth = 3, linecolor=:black)
#p1 = plot!(1e7./ν, 100*((I10_ss.+ieI10_ss)./I10_nors_ss.-1), linecolor=:red, alpha=0.3)
p1 = plot!(1e7./ν, 100*((I10_conv.+ieI10_ss_conv)./I10_conv_noRS.-1), linewidth = 3, linecolor=:red)

#p2 = plot(1e7./ν, 100*((Q10.+ieQ10)./Q10_nors.-1), linecolor=:black, ylabel="R₂ [%]", alpha=0.3)
p2 = plot(1e7./ν, 100*((Q10_conv.+ieQ10_conv)./Q10_conv_noRS.-1), linewidth = 3, linecolor=:black)
#p2 = plot!(1e7./ν, 100*((Q10_ss.+ieQ10_ss)./Q10_nors_ss.-1), linecolor=:red, alpha=0.3)
p2 = plot!(1e7./ν, 100*((Q10_conv.+ieQ10_ss_conv)./Q10_conv_noRS.-1), linewidth = 3, linecolor=:red)

#p3 = plot(1e7./ν, 100*((U10.+ieU10)./U10_nors.-1), linecolor=:black, ylabel="R₃ [%]", alpha=0.3)
p3 = plot(1e7./ν, 100*((U10_conv.+ieU10_conv)./U10_conv_noRS.-1), linewidth = 3, linecolor=:black)
#p3 = plot!(1e7./ν, 100*((U10_ss.+ieU10_ss)./U10_nors_ss.-1), linecolor=:red, alpha=0.3)
p3 = plot!(1e7./ν, 100*((U10_conv.+ieU10_ss_conv)./U10_conv_noRS.-1), linewidth = 3, linecolor=:red)

#q1 = plot(1e7./ν, 100*((I70.+ieI70)./I70_nors.-1), linecolor=:black, ylabel="R₁ [%]", alpha=0.3)
q1 = plot(1e7./ν, 100*((I70_conv.+ieI70_conv)./I70_conv_noRS.-1), linewidth = 3, linecolor=:black)
#q1 = plot!(1e7./ν, 100*((I70_ss.+ieI70_ss)./I70_nors_ss.-1), linecolor=:red, alpha=0.3)
q1 = plot!(1e7./ν, 100*((I70_conv.+ieI70_ss_conv)./I70_conv_noRS.-1), linewidth = 3, linecolor=:red)

#q2 = plot(1e7./ν, 100*((Q70.+ieQ70)./Q70_nors.-1), linecolor=:black, ylabel="R₂ [%]", alpha=0.3)
q2 = plot(1e7./ν, 100*((Q70_conv.+ieQ70_conv)./Q70_conv_noRS.-1), linewidth = 3, linecolor=:black)
#q2 = plot!(1e7./ν, 100*((Q70_ss.+ieQ70_ss)./Q70_nors_ss.-1), linecolor=:red, alpha=0.3)
q2 = plot!(1e7./ν, 100*((Q70_conv.+ieQ70_ss_conv)./Q70_conv_noRS.-1), linewidth = 3, linecolor=:red)

#q3 = plot(1e7./ν, 100*((U70.+ieU70)./U70_nors.-1), linecolor=:black, ylabel="R₃ [%]", alpha=0.3)
q3 = plot(1e7./ν, 100*((U70_conv.+ieU70_conv)./U70_conv_noRS.-1), linewidth = 3, linecolor=:black)
#q3 = plot!(1e7./ν, 100*((U70_ss.+ieU70_ss)./U70_nors_ss.-1), linecolor=:red, alpha=0.3)
q3 = plot!(1e7./ν, 100*((U70_conv.+ieU70_ss_conv)./U70_conv_noRS.-1), linewidth = 3, linecolor=:red)

plot(p1, q1, p2, q2, p3, q3, layout = l, link=:y, legend = false, title = ["SZA=10ᵒ" "SZA=70ᵒ" "" "" "" ""], titlefont = font(10))
savefig("abc.png")