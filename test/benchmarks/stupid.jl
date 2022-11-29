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
using JLD2
using DelimitedFiles

# Benchmarks: http://web.gps.caltech.edu/~vijay/Rayleigh_Scattering_Tables/STOKES/
##

# Load YAML files into parameter struct
parameters = 
    parameters_from_yaml("test/test_parameters/FraunhoferMockParameters.yaml");
#parameters.depol = 0.041362343961163395 #0.028 #(0.028: Cabannes), (0.1032: Rayleigh) 
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
effT = 300.  #(model.profile.vcd_dry' * model.profile.T) / sum(model.profile.vcd_dry);
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
CoreRT.getRamanSSProp!(RS_type, parameters.depol, 1e7/ν̃, ν);

RS_type.F₀ = zeros(model.params.polarization_type.n, length(ν))
RS_type.F₀[1,:].=1.

oR, oT, oieR, oieT = CoreRT.rt_run_test(RS_type,model,iBand);

RS_type = InelasticScattering.noRS(
    fscattRayl  = [FT(1)],
    ϖ_Cabannes  = [FT(1)], 
    bandSpecLim = [],
    iBand       = [1],
    F₀          = zeros(FT,1,1));
RS_type.F₀ = zeros(model.params.polarization_type.n, length(ν))
RS_type.F₀[1,:].=1.

oRnoRS, oTnoRS, _, _ = CoreRT.rt_run_test(RS_type,model,iBand);

#test = jldopen("/home/sanghavi/debugRay2_.jld2")
#@load "/home/sanghavi/debugRRS2.jld2"
#RayJ₀ = test["J₀"];
#test = jldopen("/home/sanghavi/debugRRS2.jld2")
#ieJ₀ = test["ieJ₀"];
#test = jldopen("/home/sanghavi/debugCab2.jld2")
#CabJ₀ = test["J₀"];
#RayJ₀[:,:,640] == sum(ieJ₀[:,:,640,:], dims=3) + CabJ₀[:,:,640]
#RayJ₀[:,:,640] ≈ sum(ieJ₀[:,:,640,:], dims=3) + CabJ₀[:,:,640]
@load "/home/sanghavi/debugCab4.jld2"
@load "/home/sanghavi/debugRRS4.jld2"
@load "/home/sanghavi/debugRay3.jld2"

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