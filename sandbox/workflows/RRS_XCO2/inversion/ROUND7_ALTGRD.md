# Round 7_altgrd: fixed-altitude vertical grid for surface pressure

Status: implementation and pilot set-up on branch `experiment/round7-altgrd`
(local only, never pushed), based on `7acab57000dae259207a6760faae156cfa1734f6`,
the exact code Round 7 ran on. **Scientific choices below await Suniti's
approval**, in particular the aerosol rule (section 3).

## 1. What changed and why

Until now vSmartMOM linearized surface pressure `p_surf` on a grid whose
pressure interfaces above the surface are fixed: only `p_half[end] = p_surf`
moves, so only the bottom layer's column, thickness and mid-layer pressure
respond (`dp_full[end]/dp_surf = 1/2`). That convention is kept, unchanged and
as the default, as `FixedPressureGrid()`.

The new `FixedAltitudeGrid()` holds every interface fixed in geometric
altitude above the surface, with layer `T`, `q` and VMR fixed on those
altitude layers. This is the natural convention for a fixed-height observer
(flux tower, aircraft at fixed altitude): the layering must not move around
the observer when `p_surf` changes.

Both are subtypes of `AbstractVerticalGridConvention` and are selected by
multiple dispatch at the two places that must agree:

* `surface_pressure_interfaces(convention, p_half, p_surf)` rebuilds the
  forward profile at a trial `p_surf` (used by `apply_retrieval_state!` and the
  XCO2 diagnostic in `VSmartMOMForward.jl`);
* `model_from_parameters(LinMode(), params; vertical_grid=convention)` (and the
  planned `OCO_RRS_synth()` constructor, which forwards keywords) selects the
  matching analytic surface-pressure tangents.

`OCOForwardEvaluator(; vertical_grid=...)` carries the convention through every
evaluation. The default everywhere is `FixedPressureGrid()`.

## 2. Formulation (exact, not first order)

vSmartMOM's layer hydrostatics (`compute_atmos_profile_fields`, constant
`g₀ = 9.8032465 m s⁻²`, moist-air molar mass `Mᵢ` fixed by `qᵢ`) give

    Δzᵢ = (R Tᵢ / (g₀ Mᵢ)) ln(p_half[i+1]/p_half[i]),
    N_dry,i = x_dry,i Nₐ (p_half[i+1] - p_half[i]) / (Mᵢ g₀)   (and N_H₂O likewise).

Fixing every `Δzᵢ`, `Tᵢ`, `Mᵢ` fixes every pressure ratio, so the whole grid
scales: `p_half[i](p_surf) = p_half_ref[i] · p_surf/p_surf_ref`. This is exact
for the model's isothermal-layer, constant-gravity hydrostatics (no
linearization). Hence, for every layer,

    dp_full,i/dp_surf = p_full,i/p_surf,   dN_i/dp_surf = N_i/p_surf,   dΔzᵢ/dp_surf = 0.

Optical-depth tangents per absorber/scatterer:

| Term | FixedPressureGrid (unchanged) | FixedAltitudeGrid (new) |
|---|---|---|
| Line/LUT gas, cross section | bottom layer: `(dσ/dp)·N·vmr/2` | every layer: `(dσ/dp)(p_full,i)·N_i·vmr_i·p_full,i/p_surf` |
| Line/LUT gas, column | bottom layer: `τ·(dN/dp_surf)/N` | every layer: `τ_i/p_surf` |
| CIA, MT_CKD (`∝ p_full·N`) | bottom layer: `τ(dN/N + 1/(2p_full))` | every layer: `2τ_i/p_surf` |
| Rayleigh `τ_tot(p_surf)·N_i/ΣN` | quotient rule, bottom column only | `τ_i/p_surf` (normalized weights invariant) |
| Aerosol, altitude `LogNormal` | interfaces above shift with bottom `Δz` | **zero** (no interface moves in altitude) |
| Aerosol, pressure `Normal(p₀,σp)` | bottom mid-pressure moves | every mid-pressure moves, `Δz` fixed |
| XCO2 weights `N_i/ΣN` | bottom weight changes | independent of `p_surf` |

At the reference `p_surf` both conventions return the bit-identical profile
(`surface_pressure_interfaces` returns an exact copy), so the Round 7
observations, generated at 1000 hPa, remain valid without regeneration.

Known limitations (both conventions, documented, not changed here): the
Rayleigh depolarization/Cabannes albedo use the column-weighted temperature,
which shifts very slightly with the bottom weight under `FixedPressureGrid`;
its tangent is not carried (relative effect ≈2e-6 on τ̇_Rayleigh in the unit
test). Under `FixedAltitudeGrid` that weight is invariant, so this is exact.
Pressure-node kinks of ABSCO/regular LUTs use the existing one-sided
convention. `update_model!`/`BatchContext` callers supply their own `p_half`
and must build it with `surface_pressure_interfaces` to follow a convention.
Profile reduction (`profile_reduction > 0`) and interior observer interfaces
are not covered by the test suite for `FixedAltitudeGrid` (only one-off
checks by the Codex reviewer, section 7); RRS_XCO2 uses neither.

## 3. Aerosol rule — for approval

The RRS_XCO2 state retrieves `ln z₀` of `LogNormal(ln z₀, σ₀)` profiles
defined in **geometric altitude above the surface** (km). Under
`FixedAltitudeGrid` I keep the aerosol fixed in altitude: the layer
allocation is unchanged when `p_surf` changes, so its `p_surf` tangent is
exactly zero. Pressure-form `Normal(p₀, σp)` aerosols (not used by RRS_XCO2)
keep `p₀, σp` fixed in pressure and are redistributed as layer pressures
scale. The alternative — aerosol attached to the air mass, i.e. fixed in
pressure — would make an altitude-parameterized aerosol move under a
fixed-altitude grid; I did not choose it.

## 4. Validation

Tests (Float64, CPU unless noted) — see section 6 for commands:

* `test/test_vertical_grid_psurf.jl` (new, 174 passed): profile tangents vs
  central FD for both conventions; per-component optics vs full rebuild for
  line-by-line, native ABSCO and regular LUT absorbers, with fixed + variable
  gases and q-driven H₂O, including **every layer's** `dσ/dp` under
  `FixedAltitudeGrid` (relative error < 2e-6 per layer); Rayleigh; both
  aerosol forms; CIA through the constructor (synthetic O₂–O₂ file); MT_CKD
  through the constructor helper on the real continuum code; TOA radiance
  Jacobian vs FD (embedded and external solar, rtol 5e-5); planned-constructor
  forwarding; reference-state identity of forward optics, non-pressure
  tangents and radiance; a Float32 `FixedAltitudeGrid` case; and a **golden
  file** (`test/reference/vertical_grid_fixed_pressure_7acab570.txt`, written
  by the unmodified 7acab570 code via `reference/generate_vertical_grid_golden.jl`)
  that the default convention reproduces bit for bit on curry with Julia
  1.12.5 (`VSMARTMOM_EXACT_GOLDEN=1`; elsewhere rtol 1e-12), and the radiance
  and its Jacobian to 1e-12.
* Existing tests, all passing (Broken = CUDA-only cases skipped on CPU):
  `test_absorption_pressure` (205 + 1 broken), `test_absorption_column_integration`
  (35), `test_selective_jacobians` (51), `test_jacobian_upstream_contract` (31),
  `test_jacobians_unit` (36), `test_local_jacobian` (5164), `test_update_model`
  (62 + 1 broken), `test_update_model_transaction` (69), `test_copy_parameters`
  (53), `test_aqua` (10, includes method ambiguities), `test_io_exports` (28).
  Sandbox: `test_forward_state_mapping` (37), new
  `test_vertical_grid_forward_mapping` (13), new
  `test_round7_altgrd_retrieval_campaign` (139 passed; its full 80-scene
  preflight needs a clean pinned checkout and was skipped on the development
  tree; it is rerun on the committed checkpoint before the pilot).
* **Round 7 Jacobian bit-identity** (GPU1, Float32, 9 streams, production
  `OCO_RRS_synth` evaluator): at five states (scene 011 and 036 truth at
  1000 hPa and at 1003.7/993.8 hPa, scene 001 prior) the default and explicit
  `FixedPressureGrid` measurement, all 30 Jacobian columns and XCO2 are
  bit-identical to the unmodified `7acab570` checkout. With
  `FixedAltitudeGrid` at 1000 hPa the measurement and columns 2–30 are
  bit-identical; only the `p_surf` column differs (13–14 % relative).
* RRS-level FD of the `p_surf` column on GPU: see section 5.

## 5. RRS-level finite differences

Production RRS_XCO2 evaluator (`OCOForwardEvaluator`, `OCO_RRS_synth`,
external solar, 9 streams, three bands through the synthetic OCO instrument),
GPU1, central difference of the full instrument-space measurement in
`p_surf` against the analytic column, at the truth states of scenes 011
(urban, aerosol, SIF off) and 036 (rural, aerosol, SIF on). Relative L2
error ‖K₁ − FD‖/‖FD‖:

| Precision | Step | FixedAltitudeGrid 011 / 036 | FixedPressureGrid 011 / 036 |
|---|---|---|---|
| Float64 | 0.5 hPa | 1.4e-7 / 1.4e-7 | 5.2e-7 / 1.9e-7 |
| Float64 | 0.1 hPa | 6.0e-9 / 6.1e-9 | 2.0e-7 / 5.6e-8 |
| Float64 | 0.02 hPa | 1.4e-8 / 1.4e-8 | 2.0e-7 / 5.3e-8 |
| Float32 | 2 hPa | 1.2e-2 / 1.2e-2 | 2.3e-2 / 2.3e-2 |
| Float32 | 0.5 hPa | 3.1e-2 / 3.1e-2 | 6.9e-2 / 6.9e-2 |

Float64 confirms the FixedAltitudeGrid column to ≈1e-8. The
FixedPressureGrid floor of ≈2e-7 is plausibly the untracked
depolarization-temperature effect noted above (not separately verified).

In Float32 the error grows as the step shrinks, which is consistent with
Float32 cancellation in the FD (it does not by itself exclude a small Float32
tangent error); noise-weighted errors at 2 hPa are
0.5–0.6 % (altitude) and 0.8–1.0 % (pressure).

## 5b. Pilot (run at checkpoint `4075fa16a0adfa5bf524b32d44b075854ab4e2d5`)

Sealed release `RRS_XCO2/bottom_layer_XCO2_retrievals/round7_altgrd_fixed_sif_imperfect_correction/deployment/round7_altgrd_release.json`;
clean-checkout campaign test 183/183 including the full 80-scene preflight
(44); `run_round7_altgrd_pilot.sh --check-only`, then `--smoke-only` (states
001 and 036, member 11: both converged, fit_ok), then the full pilot: 160/160
solves (members 01 and 11 of all 80 scenes) on curry physical GPU1, 0 failed
solves. The SIF-on block needed two relaunches after a GPU-compiler
`InvalidIRError` at its first kernel compile. The suspected, circumstantial
cause is another session rewriting the shared `~/.julia/compiled` caches;
the root cause was not isolated, and the private first depot does not shield
the runs from that depot. See the report's execution log. Results and the comparison with Round 7 are in
`ROUND7_ALTGRD_PILOT_REPORT.md` (copy of the report in the output tree), made
with `retrieval_setup/compare_round7_altgrd_pilot.jl`.

## 6. Reproduction

From `test/`, with Julia 1.12.5:

    env JULIA_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1 CUDA_VISIBLE_DEVICES='' \
        julia --startup-file=no --project=. test_vertical_grid_psurf.jl

Sandbox tests use `--project=..` with the Round 7 environment exported
(`round7_altgrd_environment.sh` sets it).

RRS-level drivers (GPU1, `--project=<worktree>`, environment as above):

    OUT=base.jls julia retrieval_setup/vertical_grid_jacobian_states.jl     # on 7acab570
    OUT=new.jls [GRID=pressure|altitude] julia retrieval_setup/vertical_grid_jacobian_states.jl
    julia retrieval_setup/compare_vertical_grid_jacobians.jl base.jls new.jls
    OUT=fd.jls GRID=altitude FTYPE=Float64 STEPS=0.5,0.1,0.02 julia retrieval_setup/vertical_grid_psurf_fd.jl

The bit-identity baseline in section 4 was produced by the first command run
in this worktree before any change (HEAD = 7acab570, clean); the drivers were
then promoted into the branch with path parameters taken from the environment.

## 7. Reviews

Commit 1 (implementation, tests, campaign set-up):

* Independent Claude reviewer (fresh context): no blockers; default path
  numerically unchanged; derivation exact; campaign isolated. Findings
  addressed: RRS-level FD missing (section 5 added); exports undocumented
  (added to `docs/src/pages/api/core_rt.md`, concept page example updated);
  `profile_reduction > 0` caveat for the existing convention (documented);
  runner GPU guard (added); bash messages and `${1:?}` (fixed); vacuous default
  test (replaced by the 7acab570 golden file); Float32 case (added). Not
  changed: generic non-Normal/LogNormal aerosol tangent (unreachable, LinMode
  rejects other distributions before the tangent is formed).
* Codex (`gpt-6.1-sol`, xhigh, read-only). The requested
  `codex exec review --uncommitted ... "<prompt>"` form is rejected by the
  installed CLI ("--uncommitted cannot be used with [PROMPT]", and `-s` must
  precede `review`), so two runs were made: the built-in `codex exec -s
  read-only review --uncommitted` and `codex exec -s read-only` with the full
  checklist prompt. No numerical or physics defect found (Codex independently
  re-checked the golden tangents and fixed-altitude FD, and in one-off checks
  not added to the suite, profile reduction and an interior observer). Three P2 campaign findings, all fixed:
  launch-specific member selection was in resume-validated provenance (moved
  to write-time provenance); spectroscopy table contents were outside the
  identity (now hashed into the input set and sealed in the release record);
  the runner accepted any GPU UUID (now must equal `nvidia-smi -i 1`).
  (The installed CLI needs `-s read-only` before the subcommand. One
  checklist run had its working directory in another checkout and reviewed
  this worktree through absolute paths; nothing was modified.)
* Re-review of the fixes. Claude reviewer: all fixes verified; remaining
  items fixed — bitwise golden comparison made opt-in
  (`VSMARTMOM_EXACT_GOLDEN=1`, generator platform recorded), rtol 1e-12
  otherwise; note wording on one-off checks; validation drivers committed
  under `retrieval_setup/`. Not changed: repeated hashing of the 3.7 GB of
  spectroscopy tables per launch (minutes, kept for safety); a global bash
  array name. Codex (`exec -C <worktree> -s read-only`, xhigh): no P1/P2;
  two P3 note-accuracy items (skipped full preflight not disclosed; Float32
  FD attribution too strong), both fixed.

Commit 2 (pilot report, comparison script, this section 5b):

* Codex (`exec -C <worktree> -s read-only`, xhigh): numbers verified against
  all 160 CSV/NetCDF pairs and the 800 Round 7 noisy-member files; three P2
  findings (interpretation wrongly described the change as tangent-only after
  the first step; InvalidIRError attribution stated too strongly; script
  silently tolerated missing files) and one P3 (timing attribution). Coverage
  is now reported (160/160 pairs, all scenes with ten members).
* Independent Claude reviewer: numbers and script confirmed; findings fixed:
  non-converged solves described by their actual trial histories (two
  divergences, γ → 10⁴, still drifting); truth bias reported by class (the
  pooled improvement was a cancellation; aerosol bias reverses sign);
  χ²-based "not a fit failure" claim softened and the convergence-path
  caveat added; maxima attributed to individual solves; cache-churn timeline
  corrected and the private-depot limitation stated; 006/11 provenance
  (diagnostic solve) and Round 7 host split stated; Float32 FD figures added;
  Rounds 1–2 timing argument; script now asserts state model/grid and emits
  the class tables and non-converged list.
* Re-review (both reviewers): fixes verified against the data; remaining
  wording items fixed — desert no-aerosol O₂A χ² change (up to −0.13) stated
  instead of "χ² unchanged"; the attribution hedge restored (it had been lost
  when the report was regenerated) and repeated here; Round 7 host split
  corrected (curry 440 no-aerosol + 165 SIF-on aerosol, wurst 275 aerosol);
  Float32 FD numbers given with their step; the timing statement limited to
  measured totals (a final Codex check caught a remaining attribution of the
  0.11 s model-build increase and a misquoted χ² range; both corrected); Rounds 1–2 described as unpinned; the alternative
  (optimizer-damping) reading of the non-converged solves added; smoke
  outputs identified as official pilot outputs; hand-computed numbers in §1/§5
  labelled as such.
