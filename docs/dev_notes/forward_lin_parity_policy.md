# Forward / Linearized Parity Policy

**Status:** binding for all new solver features (adopted 2026-08-23, after
the multisensor-merge review found three silent forward/lin divergences).

The forward radiance returned by `rt_run(model, lin_model, ...)` must be
**identical** to `rt_run(model)` for the same configured model, and the
Jacobians must differentiate **exactly the model the forward solver
evaluates** — not a near-copy. Historical violations this policy exists to
prevent (all found by review, all shipped silently):

- the δBGE fit ran on different angular domains in the forward and
  linearized truncation (Jacobians of a different phase function);
- `update_model!` recomputed τ_aer with a different floating-point
  association and a missing spectral branch (batch updates drifted from
  fresh builds);
- the linearized driver ignored the configured Fourier-convergence
  strategy (different moment counts for the same model).

## Rules

1. **Single definition.** Any physics expression consumed by more than one
   path (forward build / linearized build / batch update / cache replay /
   correction term) lives in exactly one function that all paths call.
   Never re-type a formula at a second call site — extract it.
   Standing examples: `_aerosol_τ_slice` (all five τ_aer producers),
   `_layer_absorption_cross_section` (forward + lin absorption),
   `_exact_ss_accumulate!` (TMS correction + `rt_run_ss_exact`),
   `phase_matrix_first_column` (owns scattering-angle & rotation
   conventions).

2. **Match or reject — never silently differ.** A solver feature
   (convergence strategy, single-scattering correction, source type,
   observer mode, ...) is either implemented in *both* the forward and
   linearized drivers with the same semantics, or the driver that lacks it
   **throws an `ArgumentError` naming the missing support**. Silent
   fallback to different behavior is a correctness bug by definition.
   Standing examples: `rt_run_lin` applies the identical
   `IntensityConvergence` test on its own forward accumulators; it rejects
   `TMSCorrection` and external-solar interior heights outright.

3. **Same loop bounds, same criteria.** Fourier loop ceilings come from
   one derivation (`component_m_max` traits); runtime early-exit uses one
   criterion (`_fourier_full_passes`) on the azimuth-weighted accumulators
   in both drivers.

4. **Strategy types, not booleans.** Switchable behavior is expressed as a
   small abstract-type hierarchy (`AbstractFourierConvergence`,
   `AbstractSingleScatteringCorrection`) selected via
   `RTNumericalParameters` / YAML, with the default subtype bit-identical
   to the historical solver. New variants are new subtypes — never a flag
   argument threaded through call chains.

5. **Parity is tested.** Every feature that touches both drivers ships a
   test asserting either bit-level forward parity
   (`test_lin_forward_loop_parity.jl` pattern) or the explicit rejection
   of rule 2. A feature without its parity test is not done.

## Review checklist for PRs touching solver paths

- [ ] No formula duplicated across forward/lin/update/replay (rule 1)
- [ ] Every unsupported combination throws with a message naming the gap (rule 2)
- [ ] Loop bounds/criteria shared, not re-derived (rule 3)
- [ ] Behavior switches are strategy subtypes with bit-identical defaults (rule 4)
- [ ] Parity or rejection test included (rule 5)
