# StandaloneSS Session Handoff

Last updated: 2026-05-06

## Repo State

- Branch: `SS-exact`
- Latest commit: `0e65251 Document Stokes IQ polarization option`
- Tracked working tree state at handoff: clean
- This handoff file is still untracked, along with the existing local
  `dev_notes/`, `docs/dev_notes/`, `.codex`, `sandbox/`, and
  `test/vlidort_baseline/` work areas.
- `test/vlidort_baseline/` appears to be Claude's parallel work area. Avoid
  editing it unless explicitly asked.

## Latest Commit Stack To Remember

- `0e65251 Document Stokes IQ polarization option`
- `8f10f32 Add config surface wind Jacobian helper`
- `19af355 Expose Cox-Munk SS wind Jacobian`
- `dfac47c Cover Stokes IQ linearized Z moments`
- `94c7cf3 Support vector phase chain rule`
- `787c167 Support vector surface BRDF chain rule`
- `f1db12e Support vector StandaloneSS seam Jacobians`
- `436b700 Cover Stokes IQ RTModel SS adapter`
- `4fa6036 Support Stokes IQ Cox-Munk paths`
- `015f136 Add Stokes IQ polarization option`
- `ad63082 Restore docs coverage for new exports`
- `5703120 Cover IQUV Greek single-scatter path`
- `da93141 Cover analytic phase-function parser variants`
- `03cbf43 Document analytic phase-function aerosols`

## Latest Validation

Targeted tests that passed during this work:

```julia
include("StandaloneSS/runtests.jl")
```

Result after vector phase chain-rule work: `243/243` passing.

```julia
include("test_Scattering.jl")
include("test_coxmunk.jl")
```

Both passed after the Stokes-IQ and Cox-Munk changes.

Focused checks also passed for:

- `surface_brdf_wind_jacobian(surface, geometry, n_spec; polarization_type=...)`
- `surface_brdf_wind_jacobian(config)`
- Cox-Munk wind chain rule vs `ForwardDiff.jacobian`
- trimmed docs export coverage with Documenter `checkdocs=:exports`

Full `include("runtests.jl")` was attempted. First attempt used the root
project and failed immediately on missing `WignerSymbols`; after activating
`test/Project.toml`, the full harness ran past the tool's 120s transport
limit without usable output, so the stale Julia worker was stopped. Treat this
as "not completed", not as a regression from the SS changes.

## Implemented On `SS-exact`

- Standalone exact single-scatter paths 1-4 for scalar Lambertian cases.
- Vector Stokes support for the important FO target paths:
  - path 1: sun -> atmosphere -> sensor, using Rayleigh or
    `GreekCoefsSSContributor`
  - path 2: direct-beam sun -> surface -> sensor for Lambertian and Cox-Munk
- `Stokes_IQ` is now a public polarization option in Scattering, parser, docs,
  Cox-Munk, and StandaloneSS adapter coverage.
- `SSMeasurementSelector` and retrieval-vector flattening are implemented.
  Default Stokes selection keeps all Stokes components; use
  `stokes_indices=1` for I-only retrieval vectors.
- `run_exact_ss_with_jacobians` now supports vector Stokes f2 seam Jacobians
  for paths 1 and 2.
- Chain-rule combiners now cover:
  - `chain_rule_combine_dτ`
  - `chain_rule_combine_dϖ`
  - `chain_rule_combine_dP` with scalar 4D and vector-Stokes 5D `dP_dp`
  - `chain_rule_combine_surface_brdf` with scalar 3D and vector-Stokes 4D
    `dρ_dp`
- Public Cox-Munk SS wind derivative helper:
  - `surface_brdf_wind_jacobian(config)`
  - `surface_brdf_wind_jacobian(surface, geometry, n_spec; polarization_type)`
  - Reuses CoreRT `coxmunk_brdf_mueller_and_deriv`.
- Truncated first-order diagnostic helpers:
  - `truncated_ss_path1(config, l_trunc)` moment-limits atmospheric phase
    functions while preserving optical depths and single-scattering albedos.
  - `truncated_ss_path2(config, max_m)` reconstructs the direct-beam surface
    path from BRDF Fourier moments; Lambertian is exact at `max_m=1`,
    Cox-Munk mirrors the CoreRT surface-correction Fourier convention.
  - `apply_back_correction!(R_SFI, config; l_trunc, max_m)` and the RTModel
    wrapper add `exact(path1+path2) - truncated(path1+path2)` in place.

## Important Design State

The main target is still FO-style exact single scatter:

- path 1: one exact atmospheric scattering angle
- path 2: direct solar beam reflected by the surface into the sensor

Paths 3/4 are intentionally scalar/Lambertian scaffold for now. They remain
blocked for vector Stokes by validation in `src/StandaloneSS/solver.jl` and by
the Stokes-I-only path 3/4 kernels in `src/StandaloneSS/kernels.jl`.

The retrieval Jacobian seam is:

- `τ_layer[iz, ispec]`
- `ϖ_eff[iz, ispec]`
- scalar phase: `P_eff[igeom, iz, ispec]`
- vector phase: `P_eff[igeom, istokes, iz, ispec]`
- scalar surface BRDF: `surface_brdf[igeom, ispec]`
- vector surface BRDF: `surface_brdf[igeom, istokes, ispec]`

Upstream f1 should provide derivatives of those seam variables with respect to
retrieval/state parameters. The final retrieval Jacobian is the contraction
against the f2 Jacobians and then selection through `SSMeasurementSelector`.

## Recommended Next Steps

1. Let Claude finish or stabilize `test/vlidort_baseline/`; do not overlap
   there unless explicitly coordinating. Status after continuation: scalar
   solar_tester Case B exists locally and matches VLIDORT Task 1 at about
   `4e-4`-`7e-4` relative error in the external smoke.
2. Run a broader test pass in a clean shell with the test project active:

   ```bash
   cd test
   julia --project=. -e 'include("runtests.jl")'
   ```

3. Add a small user-facing doc/example for the vector SS retrieval seam:
   `run_exact_ss_with_jacobians` + `SSMeasurementSelector` +
   `surface_brdf_wind_jacobian(config)`.
4. When ready for the next code slice, prioritize VLIDORT/Siewert comparison
   around path 1/path 2 polarized FO, not vector paths 3/4.
5. Later, if needed, add a formal f1 helper for Rayleigh depolarization or
   aerosol Greek coefficient derivatives so callers do not have to assemble
   vector `dP_dp` manually.

## Continuation Notes

Additional work completed after this handoff:

- Added `examples/standalone_ss_vector_jacobian.jl`, a runnable Stokes-IQ
  Cox-Munk path-2 wind-Jacobian example with a `ForwardDiff.jacobian`
  cross-check.
- Documented the StandaloneSS retrieval seam in `docs/src/pages/jacobians.md`.
- Added the example to `docs/test_examples.jl` and added `ForwardDiff` to the
  docs project.
- Trimmed VLIDORT baseline Siewert and solar_tester configs to one spectral
  point because their case files inject spectrally constant optical properties
  and only compare the first point.
- Verified:
  - `docs/test_examples.jl` passes, including the new example.
  - `test_Scattering.jl` passes when run with the same imports as
    `test/runtests.jl`.
  - `test_coxmunk.jl` passes.
  - The Cox-Munk StandaloneSS wind chain rule matches `ForwardDiff` to
    `~1e-18` absolute difference in the example.
  - VLIDORT Case A and Case B returned control without failure from the test
    project; Case B's external smoke printed TOA/BOA relative errors below
    `7e-4`.

## Continuation Notes — 2026-05-06

- A clean broader test pass was started with the test project explicitly
  activated:

  ```julia
  import Pkg
  Pkg.activate("/home/cfranken/code/gitHub/vSmartMOM.jl/test")
  cd("/home/cfranken/code/gitHub/vSmartMOM.jl/test") do
      @time include("runtests.jl")
  end
  ```

  The first attempt used the root project and failed before tests on
  `WignerSymbols`; explicit activation fixed the environment issue. The
  uninstrumented full run was stopped after ~25 minutes with no section-level
  output. A focused SS-adjacent broad pass was then run with section markers:

  ```julia
  Scattering                  20.05s   pass
  Cox-Munk Surface             6.22s   pass
  IO Exports                   7.13s   pass
  Parameter Parser            23.28s   pass
  Phase 1c SS driver         394.65s   pass (28/28)
  StandaloneSS               104.88s   pass (247/247)
  ```

  Total SS-adjacent broad-pass wall time: `556.29s`. The Julia worker stayed
  CPU-active after returning control and was restarted to avoid leaving a
  runaway session.

- First VLIDORT/Siewert FO probes:
  - Siewert Problem IIA full tables are not a clean oracle for raw
    `run_exact_ss(...; paths=:path1)`. They include multiple scattering, so
    the raw exact-SS path is only a sign/scale smoke, not a hard validation
    gate.
  - Solar tester vector Task 3 minus Task 1 is also not comparable to raw
    `run_exact_ss(...; paths=:paths_1_2)`. It is an FO correction delta, while
    `run_exact_ss` returns absolute first-order radiance. At sza=35 deg,
    raz=0 deg, TOA-up, the exact path-1+2 I values were about
    `[2.69e-2, 2.56e-2, 2.27e-2]`, while the Task3-Task1 I deltas were about
    `[7.77e-4, -4.73e-4, 2.77e-4]`.

  Conclusion at that point: the shipped VLIDORT fixtures were not a hard gate
  for raw `exact_ss`; they needed a truncated-path reconstruction and
  FO-equivalent back-correction layer first.

- Added the truncated-path reconstruction and post-hoc back-correction layer:
  - `src/StandaloneSS/truncated.jl`
  - exported from both `StandaloneSS` and top-level `vSmartMOM`
  - covered in `test/StandaloneSS/runtests.jl`

  Validation after this slice:

  ```julia
  include("StandaloneSS/runtests.jl")
  ```

  Result: `262/262` passing in `49.0s`.

  Remaining deferred item: the raw VLIDORT/Siewert comparison is still skipped
  by choice. The next validation slice should compare the new back-corrected
  FO-equivalent quantity against the appropriate VLIDORT FO fixture, not raw
  `run_exact_ss`.
