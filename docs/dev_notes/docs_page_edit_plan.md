# vSmartMOM Docs Overhaul — Page-by-Page Edit Plan

Created: 2026-04-28
Branch: `cfranken/refactor-low-risk-cleanup`

Source inputs:
- `docs_change_plan_for_claude.md`
- `docs_change_plan_audit.md`
- `docs_story_scaffold.md`
- `theory_references.md`

This plan converts the story scaffold into concrete file edits. The goal is to
make implementation commits small enough that Claude or a human can review each
one independently.

## Guiding Decisions

- Teach `read_parameters` as the broad user-facing loader because it dispatches
  on `String`, `Dict`, and `IOSource`.
- Present `parameters_from_file`, `parameters_from_dict`, and
  `parameters_from_source` as explicit lower-level aliases.
- Keep `parameters_from_yaml` supported for this pass. Do not add `@deprecate`.
- Treat SIF as a product/data-policy blocker. Do not add a SIF docs page until
  the missing fixtures are committed, artifact-backed, or explicitly optional.
- Mark Aerosols as WIP until that API is stabilized.
- Do not claim thermal emission or fully offline source-function integration as
  implemented.
- Prefer task-based entry pages first, module reference pages second.

## Commit 1 — Planning Artifacts

Files:
- `docs/dev_notes/docs_change_plan_for_claude.md`
- `docs/dev_notes/docs_change_plan_audit.md`
- `docs/dev_notes/docs_story_scaffold.md`
- `docs/dev_notes/theory_references.md`
- `docs/dev_notes/docs_page_edit_plan.md`

Intent:
- Preserve the Claude/Codex planning trail before implementation.

Validation:
- `git diff --check`

## Commit 2 — Landing Page and Navigation Shell

Files:
- `docs/src/index.md`
- `docs/make.jl`
- `docs/src/pages/quickstart.md` (new shell)
- `docs/src/pages/jacobians.md` (new shell)
- `docs/src/pages/gpu.md` (new shell)
- `docs/src/pages/extending/surfaces.md` (new shell)
- `docs/src/pages/extending/raman.md` (new shell)

Edits:
- Make `index.md` a 60-second landing page:
  - what vSmartMOM solves;
  - one runnable zero-config example using `default_parameters()`;
  - three next-step links: Quick Start, Configure a Scene, Core RT Theory.
- Rewire `docs/make.jl` to expose the task-based top layer first:
  - Getting Started
  - Quick Start
  - Configure a Scene
  - Compute Jacobians
  - Run on GPU
  - Extend vSmartMOM
- Add placeholder task pages with "For" and "Next" lines. Keep prose minimal
  until the later commits fill them.
- Do not remove existing tutorials from the docs. Move them lower in nav as
  long-form examples.

Validation:
- `julia --startup-file=no --project=docs docs/make.jl`

## Commit 3 — IO Naming and Scene Configuration

Files:
- `docs/src/pages/IO/Overview.md`
- `docs/src/pages/IO/Examples.md`
- `docs/src/pages/IO/Schema.md`
- `docs/src/design.md`
- `docs/src/pages/geoschem_integration.md`
- `docs/src/pages/Absorption/HITRAN_Data.md`
- all tutorial `.jl` and generated `.md` files that teach parameter loading.

Edits:
- Teach `read_parameters` first.
- Explain:
  - `read_parameters(path)`
  - `read_parameters(dict)`
  - `read_parameters(src::IOSource)`
  - explicit aliases `parameters_from_file`, `parameters_from_dict`,
    `parameters_from_source`
  - supported `parameters_from_yaml` wrapper.
- Make all examples use shipped files or `default_parameters()`.
- Replace user-facing `dirname(pathof(vSmartMOM))` path hacks with
  `pkgdir(vSmartMOM)` where paths are still needed.
- Verify that Schema keys match `src/IO/Parameters.jl`.

Validation:
- `rg "dirname\\(pathof\\(vSmartMOM\\)\\)" docs/src`
- `julia --startup-file=no --project=docs docs/make.jl`

## Commit 4 — Quick Start and Example Smoke Runner

Files:
- `docs/src/pages/quickstart.md`
- `docs/test_examples.jl` (new)
- optional `config/quickstart.yaml` if we decide a real file is better than
  `default_parameters()`.

Edits:
- Write a 5-minute CPU-only quick start:
  - load default or quickstart config;
  - force CPU if using a file that could default to GPU;
  - construct model;
  - run `rt_run`;
  - inspect `R[:, 1, :]` and output shapes.
- Add a small docs smoke runner with a `@testset` for the quick-start code.
  Keep this separate from Documenter execution.

Validation:
- `julia --startup-file=no --project=docs docs/test_examples.jl`
- `julia --startup-file=no --project=docs docs/make.jl`

## Commit 5 — Jacobians Page

Files:
- `docs/src/pages/jacobians.md`
- `docs/src/pages/tutorials/Tutorial_Jacobians.jl`
- `docs/src/pages/tutorials/Tutorial_HybridAD.jl`
- `docs/src/pages/api_reference.md`

Edits:
- Add a task page for linearized RT:
  - `LinMode()`;
  - `model_from_parameters_lin`;
  - `rt_run_lin`;
  - return shapes;
  - `ParameterLayout` slicing.
- Keep the finite-difference cross-check small and CPU-only.
- Link to long-form Jacobians and HybridAD tutorials.
- Add API reference subsection for linearized RT.

Validation:
- Add one docs smoke test for the small Jacobian example if runtime is
  acceptable.

## Commit 6 — Core RT Theory Code Map

Files:
- `docs/src/pages/vSmartMOM/CoreRTTheory.md`
- `docs/src/pages/vSmartMOM/References.md`

Edits:
- Rewrite theory around the scaffold and `theory_references.md`:
  1. Vector RTE.
  2. Elemental layer.
  3. Doubling and adding.
  4. Solar source-function integration as currently implemented.
  5. Truncation.
  6. Quadrature.
  7. Linearization.
  8. Inelastic extension.
  9. Other dispatch arms.
- Every equation gets:
  - paper citation;
  - source file path;
  - line reference if stable enough, or function name if line churn is likely.
- Use current code conventions:
  - lowercase `r,t,j` for `AddedLayer`;
  - uppercase `R,T,J` for `CompositeLayer`.
- Avoid claiming fully offline unified SFI or thermal emission is implemented.

Validation:
- `julia --startup-file=no --project=docs docs/make.jl`
- Manual spot-check that all referenced symbols/files exist.

## Commit 7 — Module Overview Pages

Files:
- `docs/src/pages/Inelastic/Overview.md` (new)
- `docs/src/pages/Aerosols/Overview.md` (new, WIP callout)
- `docs/src/pages/SolarModel/Overview.md` (new)
- `docs/src/pages/Surfaces/Overview.md` (new)
- `docs/make.jl`
- `docs/src/pages/api_reference.md`

Edits:
- Add concise module pages:
  - Inelastic: `noRS`, `RRS`, `VS_0to1`, `VS_1to0`, `_plus` variants.
  - Aerosols: TOMAS-15, two-moment, refractive-index DB, WIP caveat.
  - SolarModel: spectra, transmission, artifact/scratch behavior.
  - Surfaces: one row per `BRDF_MAP` entry, parameters, use case.
- Keep SIF out until product/data policy is resolved.

Validation:
- `julia --startup-file=no --project=docs docs/make.jl`

## Commit 8 — Extending Guides

Files:
- `docs/src/pages/extending/surfaces.md`
- `docs/src/pages/extending/raman.md`

Edits:
- Surface guide:
  - subtype/constructor;
  - `create_surface_layer!`;
  - BRDF parser registration in `BRDF_MAP`;
  - YAML example;
  - focused test.
- Raman guide:
  - subtype `AbstractRamanType`;
  - optical-property/cross-section hooks;
  - RT dispatch arms;
  - warning that inelastic linearization is not currently complete.

Validation:
- Docs build only; do not add unimplemented code examples as runnable blocks.

## Commit 9 — Tutorial Path and Execution Cleanup

Files:
- `docs/src/pages/tutorials/*.jl`
- generated tutorial `.md` files
- `docs/test_examples.jl`

Edits:
- Replace brittle paths.
- Add a smoke subset:
  - Quick Start;
  - IO;
  - one small Jacobian example.
- Mark heavier tutorials as run-locally with explicit text.

Validation:
- `julia --startup-file=no --project=docs docs/test_examples.jl`
- `julia --startup-file=no --project=docs docs/make.jl`

## Commit 10 — Documenter Strictness, Step A

Files:
- `docs/make.jl`

Edits:
- Remove broad `warnonly = true` or narrow it to an explicit list if needed.
- Do not enable `checkdocs = :exports` yet.

Validation:
- `julia --startup-file=no --project=docs docs/make.jl`

## Commit 11 — Public API Reference and `checkdocs = :exports`

Files:
- `docs/src/pages/api_reference.md`
- docstrings only where needed, avoiding Sanghavi WIP files.
- `docs/make.jl`

Edits:
- Expand API reference to cover all intentionally public exports.
- Decide which inelastic helper exports are intentionally public; document
  only those, and create a follow-up issue for export cleanup if needed.
- Enable `checkdocs = :exports`.

Validation:
- `julia --startup-file=no --project=docs docs/make.jl`
- export/docstring audit command saved in final notes.

## Commit 12 — Release Notes and Migration

Files:
- `docs/src/pages/release_notes.md`
- `CHANGELOG.md` if preferred for registry-facing release notes.
- `docs/make.jl`

Edits:
- Document:
  - Julia 1.10+ support;
  - CUDA 6;
  - RTModel hierarchy;
  - loader policy;
  - `LinMode`/linearized path;
  - RAMI unit-test decision;
  - unresolved optional data policy if SIF is still unresolved.

Validation:
- Docs build.

## Legacy vSmartMOM Page Triage — Completed 2026-04-28

The legacy `pages/vSmartMOM/*` cleanup was completed as a follow-up to the
strict Documenter work:

- Trimmed `Overview.md` to a reference-layer landing page.
- Folded `Principles.md` into `CoreRTTheory.md`, `gpu.md`, and `IO/Schema.md`,
  then removed it from navigation.
- Folded the missing Cox-Munk and canopy details from `InputParametersGuide.md`
  into `IO/Schema.md`, then deleted the duplicate parameter guide.
- Deleted `Example.md`; the forward and linearized examples now live in the
  task pages and tutorials.
- Refreshed `References.md` so the JOSS software paper is the primary package
  citation and the method papers are grouped separately.

## Deferred Until Product Decision

- SIF docs page.
- SIF default-loader behavior.
- Any `@deprecate parameters_from_yaml` decision.
- Thermal emission / unified offline source-function docs.
- Inelastic linearized RT docs beyond a "not yet supported" statement.
