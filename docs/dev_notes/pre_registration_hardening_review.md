# Pre-Registration Hardening Plan — Review & Additions

Created: 2026-04-24
Reviewer: Claude (Opus 4.7, 1M context)
Companion to: `pre_registration_hardening_plan.md`
Audit scope: `sanghavi-unified` branch at commit `858a408`

This document reviews the existing hardening plan against a full audit of the
tree. The goals are (a) confirm which claims are still true, (b) flag claims
that need revision, (c) add items the plan misses, and (d) propose a concrete
"In Progress" block for the README.

---

## 1. Overall Assessment

The existing plan is the right shape. Categories (test failures, compat, dispatch
cleanup, docs truth) are correctly prioritized and the AI working protocol is
genuinely useful. Keep it.

Two systemic issues to fix before executing the plan as written:

- **Scope creep risk in Week 3.** Sections 5, 6, 9 each propose dispatch
  refactors that touch hot-loop kernels. `elemental_inelastic.jl` (735 LOC),
  `elemental_inelastic_plus.jl` (697 LOC), and `interaction_inelastic.jl`
  (705 LOC) together are >2k lines. Consolidating them in Week 3 is unrealistic
  with a 2–3 week runway. Recommend punting #9 (Inelastic consolidation) to a
  post-1.0 milestone and keeping #5/#6 conservative (source-mode + Raman
  dispatch are fine, but stop at the public boundary; do not rewrite kernels).
- **"CoreRT reference failures" as the first task is risky as framed.** The plan
  cites specific numerical deltas (e.g. `I` = 0.0213 vs tol 0.002). Those numbers
  may already be stale. Step 0 should be: re-run tests, record current deltas
  verbatim in this dev note, and only then decide "reference vs solver."

---

## 2. Audit-Verified Claims (plan is correct)

| Plan claim | Verified | Evidence |
|---|---|---|
| `SFI::Bool` threading through RT paths | ✅ | 12 occurrences; marker types `SFI`/`DNI` **already defined** at `src/CoreRT/types.jl:883-884` under `AbstractSourceType`. The refactor is mostly a migration, not a new abstraction. |
| `RS_type isa noRS` scatter | ✅ | 12 occurrences; critical path at `src/CoreRT/LayerOpticalProperties/compEffectiveLayerProperties.jl:23` and `src/CoreRT/rt_run.jl:93,182,210,384,447`. |
| Export duplication/missing | ✅ | `src/Inelastic/InelasticScattering.jl:48` — `sol_VS_0to1_plus` listed twice; `sol_VS_1to0_plus` missing. |
| SolarModel path bug | ✅ | `src/SolarModel/SolarModel.jl:140` — `pathof(RadiativeTransfer)` (stale module name). |
| RAMI cwd-sensitive path | ✅ | `test/test_rami_smoke.jl:12` — `"../sandbox/rami/RamiNoGas.yaml"`. |
| Manifests untracked | ✅ | `git ls-files '*Manifest.toml'` returns empty; `.gitignore` covers them. |
| `@assert` overuse in IO | ✅ — **worse than stated** | `src/IO/Parameters.jl` has **43 `@assert`** occurrences (plan implies "many"). |
| Dead file: `rpv_rahman.jl` | ✅ | `src/CoreRT/Surfaces/rpv_rahman.jl` is 0 bytes. |
| `types_inelastic.jl`, `rt_run_bck.jl` deleted | ✅ | Confirmed absent. |

## 3. Claims That Need Revision

### 3.1 VS-plus workspace bug — unconfirmed as described

The plan says VS-plus code reads `tmpieJ...` / `tmpieR...` fields that do not
exist on `InteractionWorkspace`. Audit of `src/CoreRT/CoreKernel/interaction_inelastic.jl`
shows `InteractionWorkspace` has fields `tmp_inv, tmpR⁻⁺, tmpR⁺⁻, tmpT⁻⁻,
tmpT⁺⁺, tmpJ₀⁻, tmpJ₀⁺, gpu_ie_mat_{A,B}, gpu_ie_src, cpu_ie{R,T,J}{⁻⁺,⁺⁻,⁻⁻,⁺⁺,₀⁻,₀⁺},
staged`. No `tmpieJ`/`tmpieR` access was found in `elemental_inelastic_plus.jl`.

**Action:** Before implementing the `interaction_workspace(::VS_0to1_plus, ...) = nothing`
gating, first reproduce a concrete VS-plus failure (run a VS-plus YAML end-to-end).
If VS-plus actually runs today, the bug may have been fixed already, and gating
it out would be a regression. If it fails, capture the error verbatim here.

### 3.2 CI compat mismatch — add to plan

Plan says "Current package compat allows Julia 1.9, 1.10, 1.11, and 1.12" and
moves everything to 1.12. Audit shows `.github/workflows/AutomatedTests.yml`
**only tests Julia 1.11.2**. So compat is wider than what CI actually verifies.

Two viable options:
- **(A)** Require Julia 1.12, test only 1.12, advertise only 1.12. Simple.
- **(B)** Require Julia 1.10+ (LTS-ish), test 1.10 + 1.12 in a matrix. Friendlier
  to users on HPC systems that lag.

My recommendation: **(B)**. Caltech/JPL users commonly have pinned Julia
versions. A single supported version in a new registered package will push
users off. If 1.10 holds up, keep it; drop it only if something forces 1.12.

---

## 4. Items the Plan Misses

### 4.1 Registration-blocking compat hygiene

- `test/Project.toml` has **no `[compat]` section at all**. For registration
  this is fine for the main `Project.toml`, but `test/Project.toml` should pin
  at least Aqua/JET when they are added.
- `docs/Project.toml` has **no `[compat]` at all** (10 deps uncapped). Fix
  before first tag — Documenter breakages silently kill docs CI.
- Main `Project.toml` has some loose ranges worth tightening:
  - `DataInterpolations = "4, 5, 6, 7, 8"` — 4 major versions is extreme
  - `Polynomials = "2, 3"` — wide
  - `Downloads`, `TOML`, `UnPack` — no compat declared
  - `ForwardDiff = "0.10"` — implicitly includes 0.10.x only; OK but make explicit

### 4.2 Repository size / data policy

This is a **registration blocker** that the plan only touches in passing.

- `sandbox/` is **73 MB, 103 tracked files**. When a user runs
  `pkg> add vSmartMOM`, Pkg does **not** fetch `sandbox/` (it follows `Project.toml`),
  but `git clone` does. For a registered package, 73 MB of scratch is noise that
  balloons the repo forever.
- `test/benchmarks/baseline_output/` tracks ~36 `.jld2` baseline files (six
  commit hashes × three scenarios). This is a ratchet — it grows monotonically.
- `test/` is 30 MB total.

**Recommended policy:**
- Move `sandbox/` out of the main repo. Options: separate `vSmartMOM-examples`
  repo, or a `sandbox` branch that never merges.
- Replace `test/benchmarks/baseline_output/*.jld2` with an **Artifact**
  (`Artifacts.toml`). Keep one canonical baseline at HEAD; historical baselines
  are git-tagged, not committed.
- Any file >1 MB tracked in git should have a documented justification or be an
  Artifact.

### 4.3 Docs truthfulness — deeper than plan suggests

Plan flags `warnonly=true` and `"path/to/..."` placeholders. Additional issues:

- **Literate tutorials are generated but never executed.** `docs/make.jl` does
  not pass `execute=true` to `Literate.markdown`. So tutorials are just text —
  a broken tutorial ships without CI noticing. Minimum viable fix: add
  `execute=true` for tutorials that don't need GPU; GPU tutorials can stay
  schematic with an explicit banner.
- **API naming inconsistency across docs.** README uses `parameters_from_yaml`,
  `docs/src/pages/tutorials/Tutorial_IO.jl:12` uses `read_parameters`, and
  `docs/src/IO/Examples.md:7` also uses `read_parameters`. Both exist, both
  work, but this will confuse first-time users. Pick one, deprecate the other
  (or alias with a comment).
- **`docs/make.jl` should set `checkdocs = :exports`** after docstring pass, so
  every exported symbol has a docstring or is explicitly opted out.
- `examples/geoschem_integration.jl` and `examples/aerosol_integration_example.jl`
  reference external NetCDF files not in the repo; they use `isfile()` guards
  so they don't error, but they **silently do nothing**. Add a one-line banner
  telling users where to get the data (or an `@warn` when the guard is hit).
- `examples/` contains `.f90` Fortran files. These aren't Julia examples; move
  them to `sandbox/` or a `legacy/` folder before registration — they'll
  confuse readers of the shipping repo.

### 4.4 Dead code beyond `rpv_rahman.jl`

- `src/CoreRT/CoreKernel/raman_kernel_test.jl` — a file named `*_test.jl`
  living inside `src/`. Runtime code path should not have scratch/test files.
  Move or rename.
- `src/Absorption/local_only/` is `.gitignore`d — fine, but ensure nothing in
  it is actually required at runtime (would break on fresh clones).

### 4.5 Public surface area & accessor magic

`src/CoreRT/types.jl:886-899` defines `Base.getproperty(::RTModel, ...)` with 11
forwarded symbols (`:τ_abs`, `:obs_geom`, `:profile`, `:max_m`, `:l_max`, etc.).
This is convenient during refactor but costs clarity for first-time users and
breaks IDE autocompletion. Before registration:

- Decide: are these forwards part of the **public API** forever, or a
  compatibility shim?
- If public: document each one in the `RTModel` docstring with a table.
- If shim: add `@deprecate` hooks pointing users at `model.atmosphere.profile`
  etc., and plan removal in 2.x.

### 4.6 Branch naming / release hygiene

- Current state: branch is `sanghavi-unified`, memory says dev branch is
  `unified-vsmartmom`. CI in `.github/workflows/Documentation.yml` deploys docs
  on both `main`, `develop`, `unified-vsmartmom`. Before registration, decide
  canonical branch names and clean up. The fewer "live" branches, the less
  confusion for new contributors.
- Add a `CHANGELOG.md`. For a 2.0.0 release (Project.toml says 2.0.0 already),
  users and tooling expect a changelog.
- Pick a versioning discipline. If the first registered version is 2.0.0, every
  breaking change after registration triggers 3.0.0 under SemVer — this is
  expensive. Consider releasing as `0.9.0` first if the API is still moving.

### 4.7 Public API consolidation (plan §12 expanded)

Concrete candidates to review:
- `parameters_from_yaml` vs `read_parameters` — pick one.
- `model_from_parameters(params)` vs `model_from_parameters(LinMode(), params)`
  — keep both, but document LinMode as the jacobian entry point, not a variant.
- Every `export` in each module's entry file should be justified. Spot-check
  `src/CoreRT/CoreRT.jl` and `src/Scattering/Scattering.jl` for internals that
  slipped into `export`.

---

## 5. Priority Re-ordering (Concrete 3-Week Schedule)

### Week 1 — Make the suite green, lock compat
- **Day 1:** Run full `Pkg.test()`, capture actual failures verbatim in this
  file. Do not trust the audit-cited numbers.
- **Day 1–2:** Fix SolarModel `pathof(RadiativeTransfer)` bug. Trivial; unblocks
  `default_solar_transmission`.
- **Day 1–2:** Fix export duplication in `InelasticScattering.jl:48`.
- **Day 2:** Fix `test_rami_smoke.jl:12` to use `pkgdir(vSmartMOM)`.
- **Day 2–4:** Diagnose 6SV1 + Natraj reference deltas. Decision gate: reference
  regen vs solver fix. Commit the decision.
- **Day 4–5:** VS-plus reality check. Run a VS-plus YAML. If it works, skip
  the workspace gating; if it fails, capture the exact error before gating.
- **Day 5:** Pick Julia compat strategy (recommend matrix 1.10 + 1.12). Update
  `Project.toml`, `.github/workflows/AutomatedTests.yml`, docs.

### Week 2 — Quality tooling + dispatch cleanup (low-risk subset)
- **Day 1–2:** Wire Aqua + JET. Start with `Aqua.test_all(; ambiguities=false)`
  and one smoke JET check on `rt_run`. Document any suppressions.
- **Day 2–3:** Replace `@assert` in `src/IO/Parameters.jl` with
  `ArgumentError` / `DomainError`. 43 sites; batch with a helper like
  `require_key(cfg, key, ctx)`.
- **Day 3–4:** Source-mode dispatch (plan §5). Scope: only the public boundary
  (`rt_run` entry, `elemental!` signature). Do not touch inner kernels.
- **Day 4–5:** Raman-mode dispatch (plan §6). Replace `RS_type isa noRS`
  checks with `has_inelastic(RS_type)` trait. Keep the existing struct
  hierarchy; just add methods.

### Week 3 — Docs, data, and release prep
- **Day 1–2:** Add `execute=true` to non-GPU Literate tutorials. Fix whatever
  breaks. Remove `warnonly=true` after docstring pass.
- **Day 2–3:** Docstring sweep for every `export`-ed symbol. Add
  `checkdocs=:exports` to `docs/make.jl`.
- **Day 3:** Data hygiene: move `sandbox/` out, convert baseline `.jld2` files
  to Artifacts, remove `.f90` from `examples/`.
- **Day 4:** Compat tightening in `Project.toml`, `test/Project.toml`,
  `docs/Project.toml`. Add `CHANGELOG.md` seeded from `git log`.
- **Day 4:** Add "In Progress" block to `README.md` (see §6 below).
- **Day 5:** Dress rehearsal: clean clone, `Pkg.test()`, docs build with no
  warnings, tag candidate.

Post-registration (explicitly descoped from Week 3):
- Plan §9 (Inelastic elemental consolidation) — 2k LOC touch; defer.
- Plan §13 (performance/allocation pass) — needs benchmark infra first.
- Surface parsing cleanup (plan §8) — valuable but not blocking.

---

## 6. Proposed README "In Progress" Block

Replace the plan's draft bullets with this — it's calibrated to what's actually
in motion and omits items that aren't actually staffed.

```markdown
## 🚧 In Progress (pre-2.0 registration)

Active work on the `sanghavi-unified` branch, targeting a registered-package
release:

- **Julia 1.12 readiness.** Compat bounds, CI matrix, and docs are being
  updated. Current CI pins Julia 1.11; targeting a 1.10 + 1.12 matrix.
- **Quality gates in CI.** Adding [Aqua.jl] and targeted [JET.jl] checks.
- **Executable documentation.** All Literate tutorials will run as part of the
  docs build. Broken examples will fail CI, not silently ship.
- **Unified Raman / source-mode dispatch.** Replacing scattered `isa`/`Bool`
  checks with a small set of dispatch traits (`has_inelastic`, source-mode
  markers) so adding new modes is local.
- **RAMI + Raman regression validation.** Ongoing validation of RRS, vibrational
  Raman (VS), and RAMI-style scenarios against external references.
- **GPU memory/workspace cleanup** for inelastic RT paths, so larger spectral
  windows fit on a single GPU.
- **Public API tightening.** Exports and docstrings are being audited; a few
  internal helpers currently re-exported will be moved to `internal` status.
- **Artifact-backed reference data.** Large test baselines and sample scenes
  are moving to Julia Artifacts so the repo stays slim.

Not in scope for the first registered release (tracked post-2.0):
- Full unification of forward / linearized / multisensor / single-scatter kernels.
- End-to-end `JET.report_package` cleanliness.
- Consolidation of `elemental_inelastic*.jl` driver code.

[Aqua.jl]: https://github.com/JuliaTesting/Aqua.jl
[JET.jl]: https://github.com/aviatesk/JET.jl
```

The "not in scope" sub-block is important. Users and contributors read roadmaps
as promises; being explicit about what won't land in 2.0 prevents disappointment
and prevents the maintainer from feeling obligated to chase stretch goals.

---

## 7. Definition of Ready — Additions to Plan's Checklist

The plan's DoR checklist is good. Add:

- `CHANGELOG.md` exists and lists breaking changes since last tagged release.
- `docs/Project.toml` has `[compat]` for at least `Documenter`, `Literate`,
  and plotting deps.
- `test/Project.toml` has `[compat]` for Aqua and JET.
- No tracked file >1 MB without a justification in `docs/dev_notes/data_inventory.md`
  (new file, short inventory).
- `examples/` contains only Julia files that either run unconditionally or print
  a clear "needs external data from X" message when guarded.
- `checkdocs=:exports` is enabled in `docs/make.jl` and the build passes.
- Branch strategy is declared: one `main`, one dev branch, clean up the rest.

---

## 8. What I'd Skip From the Plan (Opinion, not Blocker)

- **Plan §5 migration plan step "Keep compatibility wrappers for old Boolean
  call sites."** Before registration, just change the signatures. No external
  user depends on them yet. Compatibility wrappers become debt.
- **Plan §9 (Inelastic consolidation).** Defer. It's the right idea but the
  wrong timing. Doing it before strong regression coverage risks quietly
  changing Raman numerics.
- **Plan §12 "Document the public workflow rather than every internal helper"**
  — I'd go further and **hide internal helpers** via `internal` submodules or
  unexport them. Documenting internals invites users to depend on them.

---

## 9. Quick Wins (<2 hours each)

Pick these up first — visible cleanup for low effort:

1. Delete `src/CoreRT/Surfaces/rpv_rahman.jl` (0 bytes).
2. Move `src/CoreRT/CoreKernel/raman_kernel_test.jl` to `sandbox/` or delete.
3. Fix `src/Inelastic/InelasticScattering.jl:48` exports.
4. Fix `src/SolarModel/SolarModel.jl:140` (`pathof(RadiativeTransfer)` →
   `pkgdir(vSmartMOM)`).
5. Fix `test/test_rami_smoke.jl:12` to use `pkgdir`.
6. Remove `.f90` files from `examples/` (move to `legacy/` or `sandbox/`).
7. Add a one-line `@info` banner to `examples/geoschem_integration.jl` and
   `examples/aerosol_integration_example.jl` when their `isfile()` guards fail.

Stack these in a single PR titled "Pre-registration housekeeping" so review is
trivial.

---

## 10. Continuation Snapshot — 2026-04-27

Current integration branch: `sanghavi-unified`.

Current low-risk follow-up branch: `cfranken/refactor-low-risk-cleanup`
at `3b2679e` (`Make tests cwd-independent`), pushed to origin.

Recent completed hardening work:

- numerical Float32 hardening: `expm1` / stable exponential differences,
  tolerance helpers, and CPU-side precision promotion where appropriate;
- Raman dispatch cleanup: trait/method-based inelastic handling in place of
  scattered type checks;
- parser hardening: required config/source checks now throw validation errors;
- SolarModel relocation: default solar transmission no longer writes into the
  package tree;
- test/tooling cleanup: RAMI smoke removed from the unit-test dependency path,
  Aqua/schema quality gates added, top-level IO exports smoke-tested;
- compat pruning/bumping, including CUDA 6 support;
- aerosol config readers hardened;
- current sub-branch made the main test runner and aerosol tests less dependent
  on caller working directory.

Known current blockers / follow-ups before `main` + registration:

- SIF loader tests currently fail because `src/SIF_emission/sif-spectra.csv`
  and `src/SIF_emission/ficus_refl_600to800nm.dat` are referenced but not
  tracked on this branch. Decide whether to commit the small canonical data
  files, move them to artifacts, or gate the loader smoke tests when fixtures
  are absent.
- Sanghavi has uncommitted CIA/MT_CKD and parameter reorganization WIP in her
  home checkout. Do not overwrite or rebase that work until it is preserved as
  a patch or temporary branch.
- Full test pass should be rerun after the SIF fixture decision and after
  merging any Sanghavi WIP back into `sanghavi-unified`.
- Compat changes still deserve one upper-bound verification pass for the major
  jumps that were allowed.
- `parameters_from_yaml` now rejects non-YAML/TOML misuse earlier; document the
  intended migration to `parameters_from_file` before tagging.

Safe next commits while waiting on Sanghavi WIP:

1. Resolve SIF data/test policy without touching absorption or parameter WIP.
2. Finish test cwd cleanup for any remaining standalone test files.
3. Normalize test top-level execution so smoke tests run inside named
   `@testset`s.
4. Add a short `CHANGELOG.md` / registration notes section covering breaking
   parser/API behavior and compat baseline.
5. Review exported names and docstrings for the top-level IO and Raman helper
   surface, without moving implementation files that overlap Sanghavi WIP.
