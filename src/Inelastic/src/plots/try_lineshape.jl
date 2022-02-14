using InelasticScattering
using Plots
λ = 440*1.e-7 #cm
ν̃ = 1/λ
T = 250 #k
RRS = true
VS_0to1 = false
VS_1to0 = false

o2 = InelasticScattering.getMolecularConstants(InelasticScattering.O₂(), (0.2));
compute_effective_coefficents!(ν̃, T, o2);
compute_energy_levels!(o2);
compute_σ_Rayl_coeff!(o2);
compute_σ_Rayl_VibRaman_coeff_hires!(T, o2);
compute_σ_VibRaman_coeff!(T, o2);
compute_σ_RoVibRaman_coeff!(T, o2);
n2 = InelasticScattering.getMolecularConstants(InelasticScattering.N₂(), (0.8));
compute_effective_coefficents!(ν̃, T, n2);
compute_energy_levels!(n2);
compute_σ_Rayl_coeff!(n2);
compute_σ_Rayl_VibRaman_coeff_hires!(T, n2);
compute_σ_VibRaman_coeff!(T, n2);
compute_σ_RoVibRaman_coeff!(T, n2);

# Create grid:
grid_out = -250:0.002:250
σ_out = similar(collect(grid_out));
σ_out .= 0;

apply_lineshape!(0, n2.effCoeff.σ_Rayl_coeff, λ, collect(grid_out), σ_out, 1, 300.0, 28);
atmo_σ_Rayl = n2.vmr * σ_out #cross section in cm^2

if (RRS)
    # N2
    apply_lineshape!(n2.effCoeff.Δν̃_RoRaman_coeff_JtoJp2, n2.effCoeff.σ_RoRaman_coeff_JtoJp2,  λ, collect(grid_out), σ_out, 1, 300.0, 28);
    atmo_σ_RRS_JtoJp2 = n2.vmr * σ_out #cross section in cm^2

    apply_lineshape!(n2.effCoeff.Δν̃_RoRaman_coeff_JtoJm2, n2.effCoeff.σ_RoRaman_coeff_JtoJm2, λ, collect(grid_out), σ_out, 1, 300.0, 40);
    atmo_σ_RRS_JtoJm2 = n2.vmr * σ_out #cross section in cm^2

    # O2
    apply_lineshape!(o2.effCoeff.Δν̃_RoRaman_coeff_JtoJp2, o2.effCoeff.σ_RoRaman_coeff_JtoJp2, λ, collect(grid_out), σ_out, 1, 300.0, 28);
    atmo_σ_RRS_JtoJp2 += o2.vmr * σ_out #cross section in cm^2

    apply_lineshape!(o2.effCoeff.Δν̃_RoRaman_coeff_JtoJm2, o2.effCoeff.σ_RoRaman_coeff_JtoJm2, λ, collect(grid_out), σ_out, 1, 300.0, 40);
    atmo_σ_RRS_JtoJm2 += o2.vmr * σ_out #cross section in cm^2


elseif (VS_0to1)
    # N2
    xin = [n2.effCoeff.Δν̃_RoVibRaman_coeff_0to1_JtoJp2 n2.effCoeff.Δν̃_RoVibRaman_coeff_0to1_JtoJm2]
    yin = [n2.effCoeff.σ_RoVibRaman_coeff_0to1_JtoJp2 n2.effCoeff.σ_RoVibRaman_coeff_0to1_JtoJm2];
    apply_lineshape!(xin, yin, λ, collect(grid_out), σ_out, 1, 300.0, 40);
    atmo_σ_RVRS_0to1 = n2.vmr * σ_out #cross section in cm^2

    xin = n2.effCoeff.Δν̃_VibRaman_coeff_0to1_hires
    yin = n2.effCoeff.σ_VibRaman_coeff_0to1_hires
    apply_lineshape!(xin, yin, λ, collect(grid_out), σ_out, 1, 300.0, 40);
    atmo_σ_VRS_0to1 = n2.vmr * σ_out #cross section in cm^2

    # O2
    xin = [o2.effCoeff.Δν̃_RoVibRaman_coeff_0to1_JtoJp2 o2.effCoeff.Δν̃_RoVibRaman_coeff_0to1_JtoJm2]
    yin = [o2.effCoeff.σ_RoVibRaman_coeff_0to1_JtoJp2 o2.effCoeff.σ_RoVibRaman_coeff_0to1_JtoJm2];
    apply_lineshape!(xin, yin, λ, collect(grid_out), σ_out, 1, 300.0, 40);
    atmo_σ_RVRS_0to1 += o2.vmr * σ_out #cross section in cm^2

    xin = o2.effCoeff.Δν̃_VibRaman_coeff_0to1_hires
    yin = o2.effCoeff.σ_VibRaman_coeff_0to1_hires
    apply_lineshape!(xin, yin, λ, collect(grid_out), σ_out, 1, 300.0, 40);
    atmo_σ_VRS_0to1 += o2.vmr * σ_out #cross section in cm^2
elseif (VS_1to0)
    # N2
    xin = [n2.effCoeff.Δν̃_RoVibRaman_coeff_1to0_JtoJp2 n2.effCoeff.Δν̃_RoVibRaman_coeff_1to0_JtoJm2]
    yin = [n2.effCoeff.σ_RoVibRaman_coeff_1to0_JtoJp2 n2.effCoeff.σ_RoVibRaman_coeff_1to0_JtoJm2];
    apply_lineshape!(xin, yin, λ, collect(grid_out), σ_out, 1, 300.0, 40);
    atmo_σ_RVRS_1to0 = n2.vmr * σ_out #cross section in cm^2

    xin = n2.effCoeff.Δν̃_VibRaman_coeff_1to0_hires
    yin = n2.effCoeff.σ_VibRaman_coeff_1to0_hires
    apply_lineshape!(xin, yin, λ, collect(grid_out), σ_out, 1, 300.0, 40);
    atmo_σ_VRS_1to0 = n2.vmr * σ_out #cross section in cm^2

    # O2
    xin = [o2.effCoeff.Δν̃_RoVibRaman_coeff_1to0_JtoJp2 o2.effCoeff.Δν̃_RoVibRaman_coeff_1to0_JtoJm2]
    yin = [o2.effCoeff.σ_RoVibRaman_coeff_1to0_JtoJp2 o2.effCoeff.σ_RoVibRaman_coeff_1to0_JtoJm2];
    apply_lineshape!(xin, yin, λ, collect(grid_out), σ_out, 1, 300.0, 40);
    atmo_σ_RVRS_1to0 += o2.vmr * σ_out #cross section in cm^2
    
    xin = o2.effCoeff.Δν̃_VibRaman_coeff_1to0_hires
    yin = o2.effCoeff.σ_VibRaman_coeff_1to0_hires
    apply_lineshape!(xin, yin, λ, collect(grid_out), σ_out, 1, 300.0, 40);
    atmo_σ_VRS_1to0 += o2.vmr * σ_out #cross section in cm^2
end
plotly()
plot(grid_out, σ_out*1e50, yscale=:log10)