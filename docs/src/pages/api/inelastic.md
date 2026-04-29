# Inelastic Scattering API

`InelasticScattering` defines Raman/Cabannes mode types and helper routines used
by the inelastic source paths in `CoreRT`.

## Mode Types

```@docs
vSmartMOM.InelasticScattering.AbstractRamanType
vSmartMOM.InelasticScattering.noRS
vSmartMOM.InelasticScattering.noRS_plus
vSmartMOM.InelasticScattering.RRS
vSmartMOM.InelasticScattering.RRS_plus
vSmartMOM.InelasticScattering.VS_0to1
vSmartMOM.InelasticScattering.VS_1to0
vSmartMOM.InelasticScattering.VS_0to1_plus
vSmartMOM.InelasticScattering.VS_1to0_plus
```

## Molecular and Source Helpers

```@docs
vSmartMOM.InelasticScattering.compute_γ_air_Cabannes!
vSmartMOM.InelasticScattering.compute_γ_air_Rayleigh!
vSmartMOM.InelasticScattering.compute_σ_Rayl_coeff!
vSmartMOM.InelasticScattering.compute_σ_Rayl_VibRaman_coeff_hires!
vSmartMOM.InelasticScattering.compute_σ_VibRaman_coeff!
vSmartMOM.InelasticScattering.apply_lineshape!
vSmartMOM.InelasticScattering.compute_energy_levels!
vSmartMOM.InelasticScattering.compute_stellar_Rayl
vSmartMOM.InelasticScattering.get_greek_raman
```
