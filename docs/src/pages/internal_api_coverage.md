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
vSmartMOM.StandaloneSS
```

## Absorption Internals

```@docs
vSmartMOM.Absorption.AbstractCrossSectionModel
vSmartMOM.Absorption.AbstractComplexErrorFunction
vSmartMOM.Absorption.save_interpolation_model
vSmartMOM.Absorption.load_interpolation_model
vSmartMOM.Absorption.CIATable
vSmartMOM.Absorption.MTCKDTable
vSmartMOM.Absorption.MTCKDBand
vSmartMOM.Absorption.parse_cia_file
vSmartMOM.Absorption.build_cia_table
vSmartMOM.Absorption.load_cia_table
vSmartMOM.Absorption.compute_τ_cia!
vSmartMOM.Absorption.cia_σ_at_T!
vSmartMOM.Absorption.load_mtckd
vSmartMOM.Absorption.build_mtckd_band
vSmartMOM.Absorption.compute_τ_h2o_continuum!
```

## Scattering Internals

```@docs
vSmartMOM.Scattering.greek_coefficients_from_scattering_matrix
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

## StandaloneSS Internals

```@docs
vSmartMOM.StandaloneSS.ExactSSConfig
vSmartMOM.StandaloneSS.SSGeometry
vSmartMOM.StandaloneSS.LambertianSSSurface
vSmartMOM.StandaloneSS.CoxMunkSSSurface
vSmartMOM.StandaloneSS.GreekCoefsSSContributor
vSmartMOM.StandaloneSS.TruncatedAndExactScatteringOpticalProperties
vSmartMOM.StandaloneSS.SSMeasurementSelector
vSmartMOM.StandaloneSS.run_exact_ss
vSmartMOM.StandaloneSS.run_exact_ss_with_jacobians
vSmartMOM.StandaloneSS.exact_ss_config_from_model
vSmartMOM.StandaloneSS.selected_measurements
vSmartMOM.StandaloneSS.selected_measurement_jacobian
vSmartMOM.StandaloneSS.surface_brdf_wind_jacobian
vSmartMOM.StandaloneSS.determine_required_l_aerosol
vSmartMOM.StandaloneSS.determine_required_l_from_moments
vSmartMOM.StandaloneSS.determine_required_nbrdf
vSmartMOM.StandaloneSS.determine_required_nbrdf_coxmunk
vSmartMOM.StandaloneSS.determine_required_nquad
vSmartMOM.StandaloneSS.determine_required_nquad_inner
vSmartMOM.StandaloneSS.determine_required_nstreams
vSmartMOM.StandaloneSS.chain_rule_combine_dτ
vSmartMOM.StandaloneSS.chain_rule_combine_dϖ
vSmartMOM.StandaloneSS.chain_rule_combine_dP
vSmartMOM.StandaloneSS.chain_rule_combine_surface_brdf
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
