# Developer API Coverage

This page makes exported developer-facing helpers canonical for Documenter.
Most users should start from the [Library](api_reference.md), task pages, or
module overview pages instead.

## Package Modules

```@docs
vSmartMOM
vSmartMOM.CoreRT
vSmartMOM.IO
vSmartMOM.InelasticScattering
vSmartMOM.Aerosols
```

## Absorption Internals

```@docs
vSmartMOM.Absorption.AbstractCrossSectionModel
vSmartMOM.Absorption.AbstractComplexErrorFunction
vSmartMOM.Absorption.save_interpolation_model
vSmartMOM.Absorption.load_interpolation_model
```

## Scattering Internals

```@docs
vSmartMOM.Scattering.compute_aerosol_optical_properties_gpu
vSmartMOM.Scattering.MiePrecisionPolicy
vSmartMOM.Scattering.NativeFloat64
vSmartMOM.Scattering.DSEmulated
vSmartMOM.Scattering.DoubleSingle
vSmartMOM.Scattering.ComplexDS
vSmartMOM.Scattering.NeumaierAccum
vSmartMOM.Scattering.neumaier_add
vSmartMOM.Scattering.neumaier_sum
vSmartMOM.Scattering.compute_Sl
vSmartMOM.Scattering.compute_mie_π_τ
```

## Core RT Internals

```@docs
vSmartMOM.CoreRT.RTModelLin
vSmartMOM.CoreRT.OpticalPropertyJacobian
vSmartMOM.CoreRT.RawAerosolJacobian
vSmartMOM.CoreRT.rt_run_ss
vSmartMOM.CoreRT.lin_added_layer_all_params!
vSmartMOM.CoreRT.delta_m_forward
vSmartMOM.CoreRT.delta_m_truncation_lin
```

## Inelastic Internals

```@docs
vSmartMOM.InelasticScattering.get_greek_raman_VS
vSmartMOM.InelasticScattering.has_inelastic
vSmartMOM.InelasticScattering.uses_cabannes_phase
vSmartMOM.InelasticScattering.needs_interaction_workspace
vSmartMOM.InelasticScattering.needs_rayleigh_expansion
vSmartMOM.InelasticScattering.normalize_raman_weights!
vSmartMOM.InelasticScattering.compute_effective_coefficents!
vSmartMOM.InelasticScattering.compute_σ_RoVibRaman_coeff!
vSmartMOM.InelasticScattering.get_n₀_n₁
vSmartMOM.InelasticScattering.getRamanSSProp!
vSmartMOM.InelasticScattering.sol_RRS
vSmartMOM.InelasticScattering.sol_VS_0to1
vSmartMOM.InelasticScattering.sol_VS_1to0
vSmartMOM.InelasticScattering.sol_VS_0to1_plus
vSmartMOM.InelasticScattering.sol_VS_1to0_plus
```

## Aerosol Internals

```@docs
vSmartMOM.Aerosols.compute_optical_properties
```

## Architecture Internals

```@docs
vSmartMOM.Architectures.AbstractArchitecture
vSmartMOM.Architectures.@hascuda
```

## IO Internals

```@docs
vSmartMOM.IO.Formats.load_config
```
