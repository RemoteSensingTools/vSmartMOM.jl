##
using Revise
using RadiativeTransfer, RadiativeTransfer.vSmartMOM
using Statistics

# Load YAML files into parameter struct
parameters = parameters_from_yaml("test/test_parameters/O2Parameters.yaml");
# Create model struct (precomputes optical properties) from parameters
model      = model_from_parameters(parameters);

iBand = 1
FT = Float64

#### Compute all Raman properties
ν = model.params.spec_bands[iBand]
ν̃ = mean(ν);
# Find central reference index for RRS:
i_ref = argmin(abs.(ν .- ν̃))
# TODO_VS: incident wavelength for VS
# TODO_VS: λ_vs_in (get input)
# TODO_VS: ν_vs_in (convert to wavenumbers)
# Effective temperature for Raman calculations
effT = (model.profile.vcd_dry' * model.profile.T) / sum(model.profile.vcd_dry);
# Define RS type
# Compute N2 and O2
n2,o2 = vSmartMOM.getRamanAtmoConstants(ν̃,effT);
#greek_raman = get_greek_raman(RS_type, n2, o2);
RS_type = vSmartMOM.RRS(
            n2=n2,
            o2=o2,
            greek_raman = vSmartMOM.Scattering.GreekCoefs([FT(1)], [FT(1)], [FT(1)], [FT(1)], [FT(1)], [FT(1)]),
            fscattRayl=FT(1),
            ϖ_λ₁λ₀=zeros(FT,1),
            i_λ₁λ₀=zeros(Int,1), 
            Z⁻⁺_λ₁λ₀=zeros(FT,1,1), 
            Z⁺⁺_λ₁λ₀=zeros(FT,1,1), 
            i_ref=i_ref,
            n_Raman=0);

vSmartMOM.get_greek_raman!(RS_type, n2, o2);
# now compute other optical parameters for Raman:
#fscattRayl = ...
#RS_type = vSmartMOM.RRS(...)
# Compute Raman SSA properties:
vSmartMOM.getRamanSSProp!(RS_type, 1e7/ν̃, ν);
# Add something here, that computes ALL the OP needed for the Raman case.
#modelRS = ...

R,ieR = rt_run(RS_type,
    model.params.polarization_type,
    model.obs_geom,
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

RnoRS = rt_run(vSmartMOM.noRS(),
            model.params.polarization_type,
            model.obs_geom,
            model.τ_rayl[1], 
            model.τ_aer[1], 
            model.quad_points,
            model.params.max_m,
            model.aerosol_optics[1],
            model.greek_rayleigh,
            model.τ_abs[1],
            model.params.brdf[1],
            model.params.architecture);



