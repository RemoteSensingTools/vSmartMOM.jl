# vSmartMOM nSpec=1 multi-layer Stokes_I numerical issue

**Status**: Surfaced 2026-05-05 by `test/vlidort_baseline/cases/case_B_solar_tester.jl`. Branch `SS-exact`. Open.

## Symptom

When all of the following hold simultaneously:

- `polarization_type = Stokes_I()`
- 23 atmospheric layers (multi-layer)
- Per-layer optical properties injected post-`model_from_parameters` via
  `model.τ_rayl[1][:, n] .= …`, `model.τ_abs[1][:, n] .= …`,
  `model.τ_aer[1][1, n] = …`, plus `model.aerosol_optics[1][1] = …`
- `LambertianSurfaceScalar(0.05)`
- `radiative_transfer.spec_bands` resolves to `nSpec = 1`
  (e.g. `"[18867.92]"`)

`rt_run(model, i_band=1)[1]` returns intensity ~50× smaller than the VLIDORT
2.8.3 reference for the same atmosphere + geometry.

Switching only to `nSpec = 2` (e.g. `"[18867.92 18867.93]"`) — same physics,
just a duplicated spectral point — restores agreement to ~5e-4.

Concrete numbers at `sza=35°`, `raz=0°`, `vza=10°` (VLIDORT
`results_solar_tester.all` Task 1, geom 1, level 1, dir 1):

| nSpec | modeled `R[1,1,1]` | truth | rel err |
| ----- | ------------------ | ----- | ------- |
| 1     | 1.309e-3           | 6.460e-2 | 0.98 |
| 2     | 6.457e-2           | 6.460e-2 | 4.6e-4 |

Same per-layer τ_rayl / τ_abs / τ_aer values are confirmed in both runs (the
injection prints them — see `inject_solar_tester_optics!` in the case file).

## Why it doesn't trip Case A (Siewert)

Case A is single-layer + `Stokes_IQUV()` and passes at ~1e-6 with `nSpec=1`.
The bug is specific to the conjunction above; either the multi-layer adding
loop, the Stokes_I scalar code path, or the per-layer injection produces a
different result for `nSpec=1` than for `nSpec=2`.

## Workaround in current Case B

`test/vlidort_baseline/configs/solar_tester.yaml` carries
`spec_bands: ["[18867.92 18867.93]"]` (two-point spectral band). The case
file compares `R[:, 1, 1]` (first spec point); both points produce identical
intensity by construction. This is documented in the case file header and
this dev note.

## Triage starting points

- Compare per-layer added-layer matrices `(r⁻⁺, t⁺⁺, j₀⁺)` between nSpec=1
  and nSpec=2 builds — they should differ only by trivial broadcasting.
- Inspect whether the `Stokes_I` scalar code path collapses a `(1,1,nSpec)`
  array to a scalar incorrectly when `nSpec=1`.
- Check whether `compEffectiveLayerProperties.jl` or the doubling kernel
  uses `nSpec` to size temporaries that misbehave at `nSpec=1`.
- Try the same Case B with `Stokes_IQU` to confirm/reject the
  Stokes_I-specific hypothesis.
