# Documenter Warning Inventory

Created: 2026-04-28

Source command:

```bash
julia --project=docs --startup-file=no docs/make.jl 2>&1 | tee /tmp/docs_build.log
```

## Docs Block

0 entries.

## Missing Docs

48 entries from Documenter's `CheckDocument` pass:

1. [open] `vSmartMOM.Absorption.AbstractCrossSectionModel`
2. [open] `vSmartMOM.Scattering.compute_aerosol_optical_properties_gpu :: Union{Tuple{FDT}, Tuple{MieModel{FDT}, Any}} where FDT<:NAI2`
3. [open] `vSmartMOM.Scattering.MiePrecisionPolicy`
4. [open] `vSmartMOM.CoreRT.OpticalPropertyJacobian`
5. [open] `vSmartMOM.CoreRT.RawAerosolJacobian`
6. [open] `vSmartMOM.InelasticScattering.get_greek_raman_VS :: Tuple{Union{vSmartMOM.InelasticScattering.VS_0to1, vSmartMOM.InelasticScattering.VS_0to1_plus, vSmartMOM.InelasticScattering.VS_1to0, vSmartMOM.InelasticScattering.VS_1to0_plus}, Any}`
7. [open] `vSmartMOM.InelasticScattering.needs_rayleigh_expansion :: Tuple{vSmartMOM.InelasticScattering.AbstractRamanType}`
8. [open] `vSmartMOM.Absorption.AbstractComplexErrorFunction`
9. [open] `vSmartMOM.InelasticScattering.has_inelastic :: Tuple{vSmartMOM.InelasticScattering.AbstractRamanType}`
10. [open] `vSmartMOM.Aerosols.compute_optical_properties :: Union{Tuple{FT}, Tuple{vSmartMOM.Aerosols.AerosolData{vSmartMOM.Aerosols.TOMAS15Scheme{FT}}, AbstractVector, vSmartMOM.Aerosols.RefractiveIndexDatabase{FT}}} where FT`
11. [open] `vSmartMOM.Aerosols.compute_optical_properties :: Union{Tuple{FT}, Tuple{vSmartMOM.Aerosols.AerosolData{vSmartMOM.Aerosols.TwoMomentScheme{FT}}, AbstractVector, vSmartMOM.Aerosols.RefractiveIndexDatabase{FT}}} where FT`
12. [open] `vSmartMOM.Absorption.load_interpolation_model :: Tuple{String}`
13. [open] `vSmartMOM.Architectures.AbstractArchitecture`
14. [open] `vSmartMOM.Scattering.DSEmulated`
15. [open] `vSmartMOM.InelasticScattering.compute_σ_RoVibRaman_coeff! :: Union{Tuple{FT}, Tuple{Any, vSmartMOM.InelasticScattering.MolecularConstants{FT}}} where FT`
16. [open] `vSmartMOM.InelasticScattering.get_n₀_n₁ :: Tuple{Any, Any}`
17. [open] `vSmartMOM.Scattering.compute_Sl :: Tuple{Integer, Integer, Integer, Bool, Vararg{Any, 7}}`
18. [open] `vSmartMOM.InelasticScattering`
19. [open] `vSmartMOM.InelasticScattering.sol_VS_1to0`
20. [open] `vSmartMOM.InelasticScattering.needs_interaction_workspace :: Tuple{vSmartMOM.InelasticScattering.AbstractRamanType}`
21. [open] `vSmartMOM.Scattering.NeumaierAccum`
22. [open] `vSmartMOM.CoreRT.RTModelLin`
23. [open] `vSmartMOM.Scattering.DoubleSingle`
24. [open] `vSmartMOM.CoreRT.rt_run_ss :: Tuple{Any}`
25. [open] `vSmartMOM.CoreRT.rt_run_ss :: Tuple{vSmartMOM.InelasticScattering.AbstractRamanType, Any, Any}`
26. [open] `vSmartMOM.Aerosols`
27. [open] `vSmartMOM.IO.Formats.load_config :: Tuple{vSmartMOM.IO.Formats.FileSource}`
28. [open] `vSmartMOM.IO`
29. [open] `vSmartMOM.CoreRT`
30. [open] `vSmartMOM.InelasticScattering.sol_RRS`
31. [open] `vSmartMOM.Scattering.NativeFloat64`
32. [open] `vSmartMOM.InelasticScattering.compute_effective_coefficents! :: Union{Tuple{FT}, Tuple{Any, Any, vSmartMOM.InelasticScattering.MolecularConstants{FT}}} where FT`
33. [open] `vSmartMOM.Absorption.save_interpolation_model :: Tuple{InterpolationModel, String}`
34. [open] `vSmartMOM.InelasticScattering.sol_VS_0to1`
35. [open] `vSmartMOM.vSmartMOM`
36. [open] `vSmartMOM.CoreRT.delta_m_truncation_lin :: Tuple{Any, Any, Any, Any, Any, Any, Any, Any, Any, Any, Int64, Any, Any, Any}`
37. [open] `vSmartMOM.InelasticScattering.getRamanSSProp! :: Tuple{Union{vSmartMOM.InelasticScattering.VS_0to1, vSmartMOM.InelasticScattering.VS_1to0}, Any, Any, Any}`
38. [open] `vSmartMOM.Scattering.ComplexDS`
39. [open] `vSmartMOM.InelasticScattering.sol_VS_0to1_plus`
40. [open] `vSmartMOM.InelasticScattering.uses_cabannes_phase :: Tuple{vSmartMOM.InelasticScattering.AbstractRamanType}`
41. [open] `vSmartMOM.Scattering.compute_mie_π_τ :: Tuple{Any, Any}`
42. [open] `vSmartMOM.InelasticScattering.sol_VS_1to0_plus`
43. [open] `vSmartMOM.InelasticScattering.normalize_raman_weights! :: Tuple{vSmartMOM.InelasticScattering.AbstractRamanType, Any, Any}`
44. [open] `vSmartMOM.Architectures.@hascuda :: Tuple{Any}`
45. [open] `vSmartMOM.Scattering.neumaier_add :: Union{Tuple{T}, Tuple{NeumaierAccum{T}, T}} where T`
46. [open] `vSmartMOM.CoreRT.lin_added_layer_all_params! :: Union{Tuple{FT}, Tuple{vSmartMOM.InelasticScattering.noRS{FT}, Any, Any, Any, Any, vSmartMOM.CoreRT.AddedLayerLin{FT}, Any, Int64}} where FT<:Real`
47. [open] `vSmartMOM.CoreRT.delta_m_forward :: NTuple{5, Any}`
48. [open] `vSmartMOM.Scattering.neumaier_sum :: Union{Tuple{NeumaierAccum{T}}, Tuple{T}} where T`

## Cross References

0 entries.
