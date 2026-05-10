# Audit of `docs_change_plan_for_claude.md`

Reviewer: Claude (Opus 4.7)
Date: 2026-04-28
Branch: `cfranken/refactor-low-risk-cleanup`
Source plan: `docs/dev_notes/docs_change_plan_for_claude.md`

This audit is structured as red-team / blue-team passes against the current
codebase, plus completeness gaps and an "examples actually run?" matrix.

## Summary

The plan covers the right surface area (landing page, IO naming, API
reference, CoreRT theory, tutorials, SIF, Documenter strictness, release
notes), but it has **factual drift against the current codebase**, **internally
inconsistent recommendations**, and **omits ~half of the existing docs
surface**.

---

## Red Team — Pass 1 (factual / scope errors)

1. **The plan ignores a whole class of existing pages.** `docs/make.jl` ships
   pages the plan never names: `docs/src/design.md`,
   `pages/vSmartMOM/Overview.md`, `Example.md`, `InputParametersGuide.md`,
   `Types.md`, `Principles.md`, `References.md`,
   `pages/Absorption/HITRAN_Data.md`, `pages/geoschem_integration.md`, and
   the full Absorption/Scattering folders. `design.md:40` and
   `pages/Absorption/HITRAN_Data.md:204-210` both still call
   `parameters_from_yaml`. The "naming policy" pass is silently incomplete.

2. **Plan's IO policy contradicts what the code already provides.**
   `src/IO/IO.jl:106` defines
   `parameters_from_yaml(src::Formats.IOSource) = parameters_from_source(src)`,
   and `pages/geoschem_integration.md:44,252` relies on
   `parameters_from_yaml(src)` for NetCDF sources. The plan's "YAML-only"
   framing would either deprecate that overload or force a doc lie. The
   codebase is already doing dispatch-based IO; the plan is treating it as
   string-extension routing.

3. **Plan's "Open Question 1" (`read_parameters` vs `parameters_from_file`) is
   already answered inconsistently in the docs.**
   `pages/IO/Overview.md:18-26` teaches `read_parameters` as the headline API,
   while the new `docs/src/index.md:14` teaches `parameters_from_file`. Asking
   Claude to pick later means landing the rename twice.

4. **SIF "expected data files not tracked" is understated.** Only
   `sif_loader.jl` exists in `src/SIF_emission/`. The loader's docstring claims
   data is "bundled" but `sif-spectra.csv` and `ficus_refl_600to800nm.dat` are
   absent on this branch. `test/test_sif.jl:11` hits
   `parameters_from_yaml("test_parameters/Phase1b_RRS_761-764nm.yaml")` and
   assumes the bundled fixtures via the loader defaults — so the loader itself
   is broken at the default call site, not just the test.

5. **Acceptance check `rg "1\.9"` will produce false positives.** That regex
   matches `0.1.9`, `julia-1.9`, the changelog entry "dropped 1.9 support,"
   etc. Use `rg -F "Julia 1.9"` or `rg "Julia\s+1\.9"`.

6. **Plan's quick-start acceptance ignores the architecture override.**
   Existing tutorials all do `params.architecture = vSmartMOM.Architectures.CPU()`
   between load and build to override the YAML's GPU default. The landing
   page glosses over this; if a user copies the example as-is from a YAML
   that defaults to GPU, it will try to load CUDA. Add the architecture
   override to the acceptance shape.

## Blue Team — Pass 1 (fixes)

A. **Enumerate the actual page set explicitly.** Add a "Pages touched / pages
   not touched" table at the top of the plan, derived from `docs/make.jl`.
   Anything not on the touched list is an explicit deferral, not an oversight.

B. **Pick the IO API now, not later.** Recommend: `read_parameters` is the
   user-facing API (one verb, dispatch on type — `String`, `Dict`,
   `IOSource`). `parameters_from_file` / `_dict` / `_source` are explicit
   aliases for code that wants its intent in the call site.
   `parameters_from_yaml` becomes a deprecated alias with a `@deprecate`
   pointing at `read_parameters`. This matches what the IO Overview page
   already teaches and lets the dispatch overload at `IO.jl:106` keep working
   coherently.

C. **Treat SIF as a registration blocker, not a documentation question.**
   Either commit the small text fixtures (Option A) before tagging 2.0.0, or
   guard `load_sif_spectrum` / `load_ficus_reflectance` against missing files
   with a clear error pointing at the data policy. The current state —
   exported public functions whose default arg points at non-existent files —
   should not ship.

D. **Tighten acceptance regexes.** Use word boundaries; pin to
   `docs/src/index.md` and the tutorial folder explicitly; run greps with
   `--type md --type jl`.

E. **Drop the `params.architecture = CPU()` landline in the quick-start** by
   shipping a `config/quickstart.yaml` next to the existing `config/`
   examples whose default is CPU. Then the landing-page snippet works
   literally.

---

## Red Team — Pass 2 (architectural / risk)

7. **"Code Path" snippet in §4 is partially fictional.**
   `src/CoreRT/CoreKernel/rt_kernel.jl` does not call
   `constructCoreOpticalProperties` itself — that wiring lives in
   `rt_run.jl`. Documenting a call sequence the reader cannot grep for is
   the exact "aspirational, not truthful" failure the plan's own Goal 1
   forbids.

8. **The plan over-promises that elemental/doubling/interaction "tells a
   reader exactly where to inspect" the implementation, but
   `src/CoreRT/CoreKernel/` has 21 files** including `_lin`, `_inelastic`,
   `_canopy`, `_ss`, `_multisensor`, `_hdrf` variants. Pointing only at the
   `noRS` elastic kernel hides where the actual algorithmic complexity now
   lives. A reader following the doc map will be surprised when they
   discover six dispatch tables.

9. **"Documenter strictness" plan mixes two unrelated concerns:** removing
   `warnonly` (broken refs) and turning on `checkdocs = :exports` (missing
   docstrings). These have very different blast radii. With ~30 newly added
   Inelastic/Aerosols/SIF exports and no doc pages, flipping
   `checkdocs = :exports` will surface dozens of warnings. Stage them
   separately so a regression in one does not block the other.

10. **Plan says "Execute Literate tutorials" but does not budget for it.**
    Tutorial_QuickStart through Tutorial_Jacobians load HITRAN artifacts and
    run RT — that's minutes per tutorial in CI, plus artifact downloads on
    cold cache. Tutorial_GPU and Tutorial_HybridAD need CUDA. This needs a
    separate "smoke" subset, not a blanket `execute=true`.

11. **No mention of `default_parameters()`**, which is exported from
    `src/vSmartMOM.jl:99` and is the *zero-config* entry point. The simplest
    correct quick start is `model_from_parameters(default_parameters())`,
    not loading a path.

12. **`vSmartMOM_Lin` / `LinMode` story is missing from the doc structure.**
    `src/vSmartMOM.jl:99-106` exports `model_from_parameters_lin`,
    `rt_run_lin`, `OpticsLin`, but the plan's API reference §3 doesn't list
    them.

13. **Plan calls `parameters_from_yaml` "compatibility/legacy" but it has no
    `@deprecate`.** Either deprecate it (with a Julia logger warning) or stop
    calling it legacy in the docs. Right now both code and docs send mixed
    signals.

## Blue Team — Pass 2 (fixes)

F. **Replace §4's invented call-path block with a verbatim transcription**
   copied from `rt_run.jl` and `rt_kernel.jl`, each line annotated with the
   exact `file.jl:LINE`. Add an "Other dispatch arms" subsection naming the
   inelastic/canopy/multisensor/single-scatter variants and pointing at where
   to go for each.

G. **Split §7 into two acceptance gates:**
   - Gate A: `warnonly = false` for cross-references (no broken `[link]`
     syntax). Land first.
   - Gate B: `checkdocs = :exports`. Land last, after every new doc page
     exists.

H. **Define a "tutorial smoke set" of 2-3 fast tutorials** (Quick Start with a
   small spectral grid, Absorption demo, IO demo) that run in CI; mark the
   rest `execute = false` with a banner "Run locally, not in CI."

I. **Add `default_parameters()` as Step 0 of the Quick Start** — it removes
   the path-handling entirely.

J. **Add a §3.5 "Linearized RT" subsection** to the API reference that
   documents `LinMode`, `model_from_parameters_lin`, `rt_run_lin`, and the
   `(R, T, dR, dT)` return shape.

K. **Either deprecate `parameters_from_yaml` for real or stop calling it
   legacy.** Recommend
   `@deprecate parameters_from_yaml(path::AbstractString) read_parameters(path)`
   (keep the `IOSource` overload as a normal method).

---

## Completeness — pages & sections to add

The current `docs/make.jl` nav lists Absorption, Scattering, vSmartMOM. The
package exports six more user-facing surface areas with **no dedicated doc
page**:

| Module | Pages exist? | Exports docstring'd? | Recommended |
|---|---|---|---|
| `InelasticScattering` | no | partial | Add `pages/Inelastic/Overview.md` (RRS vs VS theory, when to use what) + `Types.md` (noRS/RRS/VS_0to1/VS_1to0 + traits) |
| `Aerosols` | no | unknown (WIP per `src/Aerosols/Aerosols.jl:1`) | Add `pages/Aerosols/Overview.md` clearly marked "API may evolve" until WIP cleanup lands |
| `SolarModel` | no | mixed | Add `pages/SolarModel/Overview.md` + Types — covers `planck_spectrum_wn/wl`, `solar_transmission_from_file` |
| `SIF_emission` | no | yes (loaders) | Either delete loaders from public exports OR add `pages/SIF/Overview.md` once data policy is decided |
| Surfaces | only buried in API reference | yes | Add `pages/Surfaces/Overview.md` — the `BRDF_MAP` types are user-facing config strings; users will look here first |
| Architectures | only mentioned | yes | Short `pages/Architectures.md` covering CPU/GPU selection, CUDA extension behavior, threading note |

Plus content-shaped gaps:

- **Migration guide v1.x → 2.0.0** (RTModel sub-struct refactor, removed
  `vSmartMOM_Model`, new `Optics`/`Atmosphere`/`SolverConfig` accessors).
  Without this, every existing user's code breaks silently.
- **HITRAN data management** has `pages/Absorption/HITRAN_Data.md` but doesn't
  cover the new artifact-vs-scratch dispatch in
  `src/Artifacts/artifact_helper.jl` or the
  `set_hitran_edition!` / `available_hitran_editions` workflow.
- **GPU section in design or its own page** — currently spread across
  `Tutorial_GPU.jl` and CLAUDE.md.
- **CHANGELOG entry must call out:** removed `vSmartMOM_Model` (breaking),
  new `RTModel` shape, accessor functions (`get_surface`, `get_spec_bands`,
  `polarization_type`, etc.), `LinMode()` construction path, `read_parameters`
  policy, and CUDA 6 baseline.
- **Glossary / convention page** — sign convention (+ down / − up), uppercase
  R, T, J vs lowercase r, t, j, units (cm⁻¹ vs μm), profile direction
  (TOA→BOA). This is in CLAUDE.md but not user-visible.
- **Performance/benchmark page** — there is implicit data in
  `Tutorial_GPU.jl:92-93` but no canonical benchmark page.
- **Single canonical sample config** in `config/` referenced from the landing
  page, instead of pointing users at `test/test_parameters/`.

---

## Do the examples actually run?

| Example / file | Verdict | Reason |
|---|---|---|
| `docs/src/index.md:12-17` Quick example | **Won't run literally** — `"path/to/your/params.yaml"` is a placeholder. Fine as illustrative, but no shipped real path. |
| `pages/IO/Overview.md:18-26` | **Should run** — `DefaultParameters.yaml` exists at the joined path. Validate the inline Dict example actually parses. |
| `pages/IO/Overview.md:30-50` inline Dict | **Unverified** — schema may have drifted; needs an exec test |
| `Tutorial_QuickStart.jl` | **Likely runs** — `test/test_parameters/ParamsEMIT_fast.yaml` exists; needs HITRAN artifact |
| `Tutorial_Jacobians.jl` | **Likely runs** — same config, needs LinMode pathway intact |
| `Tutorial_HybridAD.jl` | **Partial** — `JacobianTestFast.yaml` exists; second block uses GPU and will fail without CUDA |
| `Tutorial_GPU.jl` | **GPU-only** — `PureRayleighParameters.yaml` exists, but YAML defaults to GPU |
| `Tutorial_CoreRT.jl` | Needs verification — uses `parameters_from_yaml(joinpath(...))` with full path |
| `Tutorial_Surfaces.jl` | Needs verification |
| `Tutorial_Canopy.jl` | **Likely runs** — `CanopyTest.yaml` exists |
| `Tutorial_IO.jl` | **Should run** — uses `read_parameters` / `DefaultParameters.yaml` |
| `test/test_sif.jl` | **Will fail** — depends on missing SIF/ficus fixtures |
| `pages/Absorption/HITRAN_Data.md` examples | Needs HITRAN artifact + uses `parameters_from_yaml` (fine, but inconsistent naming) |

**Concrete suggestion:** add a `docs/test_examples.jl` script that
`Pkg.activate("docs")`, then `include`s each tutorial `.jl` (without
rendering) under a `Test.@testset`. Run it as a separate CI job from docs
build. That gives you a real "does it run" gate without coupling Documenter
speed to RT speed.

---

## Items the plan should explicitly add to its scope

1. Update `docs/src/design.md` along with everything else (it's in nav and
   uses the old loader name).
2. Decide what to do with the unmentioned `pages/vSmartMOM/*.md` set — they
   likely contain stale type/method names from the pre-RTModel era.
3. Add migration guide section to release notes covering RTModel refactor
   (it broke `vSmartMOM_Model`).
4. Add explicit deprecation strategy for `parameters_from_yaml`.
5. Add a `docs/test_examples.jl` smoke-runner.
6. Resolve SIF before tagging — the loader's default args point at missing
   files.
7. Add a Surfaces overview page (`BRDF_MAP` types are user-facing).
8. Strip the WIP markdown (`AEROSOL_*.md`, `RESOLUTION_SUMMARY.md`, etc.)
   inside `src/Aerosols/` from any docs `include` paths.

---

## Recommended commit re-sequencing

The plan's six commits are reasonable but the order has an inversion:
tightening Documenter (#5) before the API page (#2) is finished can fail CI
on its own warnings. Suggested order:

1. `docs: align landing page, IO naming, deprecate parameters_from_yaml`
   (covers index + design.md + IO pages + tutorial loader sweep + `@deprecate`)
2. `docs: add module overview pages (Inelastic, Aerosols, SolarModel, Surfaces)`
3. `docs: expand API reference + linearized RT section`
4. `docs: map CoreRT theory to solver code` (with verbatim file:LINE refs)
5. `docs: stabilize tutorial paths + add docs/test_examples.jl`
6. `docs: tighten Documenter — warnonly off for refs only`
7. `docs: tighten Documenter — checkdocs = :exports`
8. `docs: add release notes + migration guide`

Two-step Documenter tightening protects against the "one warning blocks the
merge" pattern.
