# Documenter Warning Inventory

Created: 2026-04-28

Source command:

```bash
julia --project=docs --startup-file=no docs/make.jl 2>&1 | tee /tmp/docs_build.log
```

## Docs Block

0 entries. Resolved for Commit B by dropping `:docs_block` from `warnonly`.

## Missing Docs

48 entries from Documenter's `CheckDocument` pass. Resolved for Commit C by adding `docs/src/pages/internal_api_coverage.md` and fixing broken internal docstring links. No symbols were de-exported:

1. [resolved] `vSmartMOM.Absorption.AbstractCrossSectionModel`
2. [resolved] `vSmartMOM.Scattering.compute_aerosol_optical_properties_gpu :: Union{Tuple{FDT}, Tuple{MieModel{FDT}, Any}} where FDT<:NAI2`
3. [resolved] `vSmartMOM.Scattering.MiePrecisionPolicy`
4. [resolved] `vSmartMOM.CoreRT.OpticalPropertyJacobian`
5. [resolved] `vSmartMOM.CoreRT.RawAerosolJacobian`
6. [resolved] `vSmartMOM.InelasticScattering.get_greek_raman_VS :: Tuple{Union{vSmartMOM.InelasticScattering.VS_0to1, vSmartMOM.InelasticScattering.VS_0to1_plus, vSmartMOM.InelasticScattering.VS_1to0, vSmartMOM.InelasticScattering.VS_1to0_plus}, Any}`
7. [resolved] `vSmartMOM.InelasticScattering.needs_rayleigh_expansion :: Tuple{vSmartMOM.InelasticScattering.AbstractRamanType}`
8. [resolved] `vSmartMOM.Absorption.AbstractComplexErrorFunction`
9. [resolved] `vSmartMOM.InelasticScattering.has_inelastic :: Tuple{vSmartMOM.InelasticScattering.AbstractRamanType}`
10. [resolved] `vSmartMOM.Aerosols.compute_optical_properties :: Union{Tuple{FT}, Tuple{vSmartMOM.Aerosols.AerosolData{vSmartMOM.Aerosols.TOMAS15Scheme{FT}}, AbstractVector, vSmartMOM.Aerosols.RefractiveIndexDatabase{FT}}} where FT`
11. [resolved] `vSmartMOM.Aerosols.compute_optical_properties :: Union{Tuple{FT}, Tuple{vSmartMOM.Aerosols.AerosolData{vSmartMOM.Aerosols.TwoMomentScheme{FT}}, AbstractVector, vSmartMOM.Aerosols.RefractiveIndexDatabase{FT}}} where FT`
12. [resolved] `vSmartMOM.Absorption.load_interpolation_model :: Tuple{String}`
13. [resolved] `vSmartMOM.Architectures.AbstractArchitecture`
14. [resolved] `vSmartMOM.Scattering.DSEmulated`
15. [resolved] `vSmartMOM.InelasticScattering.compute_σ_RoVibRaman_coeff! :: Union{Tuple{FT}, Tuple{Any, vSmartMOM.InelasticScattering.MolecularConstants{FT}}} where FT`
16. [resolved] `vSmartMOM.InelasticScattering.get_n₀_n₁ :: Tuple{Any, Any}`
17. [resolved] `vSmartMOM.Scattering.compute_Sl :: Tuple{Integer, Integer, Integer, Bool, Vararg{Any, 7}}`
18. [resolved] `vSmartMOM.InelasticScattering`
19. [resolved] `vSmartMOM.InelasticScattering.sol_VS_1to0`
20. [resolved] `vSmartMOM.InelasticScattering.needs_interaction_workspace :: Tuple{vSmartMOM.InelasticScattering.AbstractRamanType}`
21. [resolved] `vSmartMOM.Scattering.NeumaierAccum`
22. [resolved] `vSmartMOM.CoreRT.RTModelLin`
23. [resolved] `vSmartMOM.Scattering.DoubleSingle`
24. [resolved] `vSmartMOM.CoreRT.rt_run_ss :: Tuple{Any}`
25. [resolved] `vSmartMOM.CoreRT.rt_run_ss :: Tuple{vSmartMOM.InelasticScattering.AbstractRamanType, Any, Any}`
26. [resolved] `vSmartMOM.Aerosols`
27. [resolved] `vSmartMOM.IO.Formats.load_config :: Tuple{vSmartMOM.IO.Formats.FileSource}`
28. [resolved] `vSmartMOM.IO`
29. [resolved] `vSmartMOM.CoreRT`
30. [resolved] `vSmartMOM.InelasticScattering.sol_RRS`
31. [resolved] `vSmartMOM.Scattering.NativeFloat64`
32. [resolved] `vSmartMOM.InelasticScattering.compute_effective_coefficents! :: Union{Tuple{FT}, Tuple{Any, Any, vSmartMOM.InelasticScattering.MolecularConstants{FT}}} where FT`
33. [resolved] `vSmartMOM.Absorption.save_interpolation_model :: Tuple{InterpolationModel, String}`
34. [resolved] `vSmartMOM.InelasticScattering.sol_VS_0to1`
35. [resolved] `vSmartMOM.vSmartMOM`
36. [resolved] `vSmartMOM.CoreRT.delta_m_truncation_lin :: Tuple{Any, Any, Any, Any, Any, Any, Any, Any, Any, Any, Int64, Any, Any, Any}`
37. [resolved] `vSmartMOM.InelasticScattering.getRamanSSProp! :: Tuple{Union{vSmartMOM.InelasticScattering.VS_0to1, vSmartMOM.InelasticScattering.VS_1to0}, Any, Any, Any}`
38. [resolved] `vSmartMOM.Scattering.ComplexDS`
39. [resolved] `vSmartMOM.InelasticScattering.sol_VS_0to1_plus`
40. [resolved] `vSmartMOM.InelasticScattering.uses_cabannes_phase :: Tuple{vSmartMOM.InelasticScattering.AbstractRamanType}`
41. [resolved] `vSmartMOM.Scattering.compute_mie_π_τ :: Tuple{Any, Any}`
42. [resolved] `vSmartMOM.InelasticScattering.sol_VS_1to0_plus`
43. [resolved] `vSmartMOM.InelasticScattering.normalize_raman_weights! :: Tuple{vSmartMOM.InelasticScattering.AbstractRamanType, Any, Any}`
44. [resolved] `vSmartMOM.Architectures.@hascuda :: Tuple{Any}`
45. [resolved] `vSmartMOM.Scattering.neumaier_add :: Union{Tuple{T}, Tuple{NeumaierAccum{T}, T}} where T`
46. [resolved] `vSmartMOM.CoreRT.lin_added_layer_all_params! :: Union{Tuple{FT}, Tuple{vSmartMOM.InelasticScattering.noRS{FT}, Any, Any, Any, Any, vSmartMOM.CoreRT.AddedLayerLin{FT}, Any, Int64}} where FT<:Real`
47. [resolved] `vSmartMOM.CoreRT.delta_m_forward :: NTuple{5, Any}`
48. [resolved] `vSmartMOM.Scattering.neumaier_sum :: Union{Tuple{NeumaierAccum{T}}, Tuple{T}} where T`

## Cross References

0 entries. Resolved for Commit B by dropping `:cross_references` from
`warnonly`.
