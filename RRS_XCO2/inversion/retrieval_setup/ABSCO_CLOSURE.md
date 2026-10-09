# ABSCO truth/retrieval forward-model closure

- Generated: `2026-09-01T16:46:17.628` UTC
- Truth state: `009` (urban, aod760_0p28, SIF off, 380 ppm CO2)
- Precision/backend: `Float32`, CUDA device `0`
- Spectroscopy: ABSCO `5.2`
- A-band H2O: rebuilt HITRAN LUT retained from the archived ABSCO-O2 truth configuration
- Retrieval SIF source: `SIF760 + mSIF*(nu-nu760)`; both values are zero for this closure state, while both Jacobian columns remain active
- Overall result: **PASS**

Exact profile-array checks:

| Band | T | p_half | p_full | q | H2O VMR | dry column | O2 VMR | CO2 VMR |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| o2a | true | true | true | true | true | true | true | true |
| weak_co2 | true | true | true | true | true | true | true | true |
| strong_co2 | true | true | true | true | true | true | true | true |

All profile arrays matched element-for-element: **true**.

| Band | comparison | max abs | relative L2 | relative max |
|---|---|---:|---:|---:|
| o2a | truth vs retrieval noRS | 2.384185791e-04 | 1.258199829e-06 | 4.249144690e-05 |
| o2a | retrieval noRS vs linearized forward | 1.420974731e-04 | 7.753829030e-06 | 2.532490450e-05 |
| o2a | OCO-processed truth vs retrieval noRS | 2.738865180e-04 | 4.396222255e-07 | 5.707761140e-06 |
| o2a | OCO-processed retrieval vs linearized forward | 8.078926991e-04 | 7.310832174e-06 | 1.683638027e-05 |
| weak_co2 | truth vs retrieval noRS | 1.287460327e-05 | 1.866892896e-07 | 3.504811197e-06 |
| weak_co2 | retrieval noRS vs linearized forward | 5.960464478e-06 | 5.441202768e-07 | 1.622597776e-06 |
| weak_co2 | OCO-processed truth vs retrieval noRS | 7.019420621e-06 | 9.633828594e-08 | 9.738601187e-07 |
| weak_co2 | OCO-processed retrieval vs linearized forward | 1.076725066e-05 | 5.283116067e-07 | 1.493826511e-06 |
| strong_co2 | truth vs retrieval noRS | 7.867813110e-06 | 2.315382582e-07 | 3.404218797e-06 |
| strong_co2 | retrieval noRS vs linearized forward | 1.072883606e-06 | 1.525556464e-07 | 4.642117020e-07 |
| strong_co2 | OCO-processed truth vs retrieval noRS | 3.505557562e-06 | 1.313848181e-07 | 1.290553979e-06 |
| strong_co2 | OCO-processed retrieval vs linearized forward | 1.147334677e-06 | 1.401294878e-07 | 4.223856818e-07 |

| Band | optical field | max abs | relative L2 |
|---|---|---:|---:|
| o2a | tau_abs | 3.814697266e-06 | 3.036772923e-08 |
| o2a | tau_rayleigh | 0.000000000e+00 | 0.000000000e+00 |
| o2a | tau_aerosol | 5.587935448e-09 | 2.927094433e-08 |
| weak_co2 | tau_abs | 1.490116119e-08 | 3.646703743e-08 |
| weak_co2 | tau_rayleigh | 0.000000000e+00 | 0.000000000e+00 |
| weak_co2 | tau_aerosol | 4.656612873e-09 | 1.451900080e-07 |
| strong_co2 | tau_abs | 4.768371582e-07 | 4.206417040e-08 |
| strong_co2 | tau_rayleigh | 0.000000000e+00 | 0.000000000e+00 |
| strong_co2 | tau_aerosol | 2.793967724e-09 | 1.334981707e-07 |

## Instrument/Jacobian consistency

The same fixed operator—`M11*I - M12*Q + M13*U`, per-cm⁻¹ to per-nm conversion, Gaussian ILS, then synthetic OCO sampling—was applied to the forward spectrum and all 30 analytic Jacobian columns.

| Band | worst column max abs | worst column relative L2 | result |
|---|---:|---:|---:|
| o2a | 4.064304449e-12 | 2.819983735e-15 | PASS |
| weak_co2 | 2.193800697e-12 | 1.711023956e-14 | PASS |
| strong_co2 | 9.166001291e-13 | 3.243968052e-15 | PASS |

Raman solve-grid retained-node max difference: `0.000000000e+00 cm⁻¹`.

OCO convolution difference when solve-only Raman shoulders are present: max abs `0.000000000e+00`, relative L2 `0.000000000e+00` (**PASS**).

Elastic noRS difference for a production-size 256-point core solved with ±234 cm⁻¹ Raman shoulders: max abs `1.764297485e-05`, relative L2 `1.404181889e-06`.

## Active regenerated A-band product closure

The active state file's regenerated Rayleigh A-band was compared directly with the retrieval noRS forward model.

High-resolution regenerated-truth/retrieval difference: max abs `2.503395081e-04`, relative L2 `2.091371141e-06`.

OCO-processed regenerated-truth/retrieval difference: max abs `3.116540350e-04`, relative L2 `1.011285894e-06` (**PASS**).

Raman shoulders are used by the RRS solve, then discarded with the retained-core index before the Mueller/ILS operator. Their physical Raman redistribution into the core is retained; the shoulder samples themselves never enter convolution.

The test uses the canonical complete output band as the aerosol interpolation anchor. Raman chunks may extend beyond it, but their aerosol phase/extinction interpolation is anchored to these same band endpoints so chunking cannot change core-band elastic physics.
