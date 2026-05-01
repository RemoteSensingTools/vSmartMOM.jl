# Documentation Change Plan for Claude Review

Created: 2026-04-28
Branch: `cfranken/refactor-low-risk-cleanup`
Purpose: plan the documentation cleanup before merging `sanghavi-unified` to
`main` and registering a new official vSmartMOM version.

This is a planning document, not the implementation. Use it to iterate with
Claude before editing the documentation in bulk.

## Current State

- `docs/src/` contains the source documentation pages. `docs/build/` is not
  tracked, so audit and edit `docs/src/` and `docs/make.jl`.
- `docs/make.jl` generates Literate tutorial markdown but does not execute the
  tutorials. It also runs Documenter with `warnonly = true`.
- The exported public symbols appear to have docstrings mechanically, but the
  public API reference page does not list every public area.
- The docs mix three parameter-loading names:
  `parameters_from_yaml`, `parameters_from_file`, and `read_parameters`.
- Current package compat is Julia `1.10, 1.11, 1.12`, while the landing page
  still advertised Julia 1.9 before the current uncommitted candidate edit.
- `src/SIF_emission/sif_loader.jl` says data samples are bundled under
  `src/SIF_emission/`, but the expected data files are not tracked on this
  branch. This also blocks the full test suite at `test/test_sif.jl`.
- There is one current uncommitted candidate edit in `docs/src/index.md`:
  it updates Julia requirements, uses `parameters_from_file` in the quick
  start, and expands the module overview. Keep, revise, or revert after review.

## Goals

1. Make docs truthful for the current code, not aspirational.
2. Make the public entry points clear for new users.
3. Explain where the solver heart lives in code:
   elemental, doubling, and adding/interaction.
4. Separate user-facing API from internal implementation details.
5. Make docs build failures meaningful before registration.

## Non-Goals for This Pass

- Do not refactor solver code.
- Do not touch Sanghavi WIP overlap files unless the WIP has been preserved and
  merged. In particular, avoid source edits in absorption, parameter parsing,
  model construction, and core type files unless necessary.
- Do not promise the offline unified source-function architecture as implemented
  unless the code actually has it.
- Do not document every internal helper as public API.

## Proposed Documentation Structure

### 1. Landing Page

File: `docs/src/index.md`

Change:
- State Julia `1.10+`.
- Prefer `parameters_from_file("...yaml")` in first example.
- Mention `parameters_from_yaml` only as YAML-only compatibility.
- Name the real submodules: `CoreRT`, `Absorption`, `Scattering`,
  `InelasticScattering`, `Aerosols`, `SolarModel`, `IO`.
- Add one sentence pointing advanced readers to "Core RT Theory".

Acceptance checks:
- `rg "1\\.9|1\\.9\\+" docs/src/index.md` returns no hits.
- First code block constructs `params`, `model`, and `R, T`.

### 2. Parameter Loading Policy

Files:
- `docs/src/pages/IO/Overview.md`
- `docs/src/pages/IO/Examples.md`
- `docs/src/pages/IO/Schema.md`
- tutorials that currently call `parameters_from_yaml`

Proposed policy:
- Canonical file loader: `parameters_from_file(path)`.
- Canonical in-memory loader: `parameters_from_dict(dict)`.
- Convenience/legacy aliases:
  - `parameters_from_yaml(path)` for `.yaml` and `.yml` only.
  - `read_parameters` as a broad convenience wrapper if we decide to keep it
    user-facing.

Open question for Claude:
- Should the docs teach `parameters_from_file` everywhere, or teach
  `read_parameters` as the one friendly API and present `parameters_from_file`
  as the explicit extension-point API?

Acceptance checks:
- `parameters_from_yaml` remains in compatibility notes and YAML-specific
  examples only.
- TOML support is documented where file loading is introduced.
- The docs mention that custom formats extend the `IO.Formats` registry.

### 3. Public API Reference

File: `docs/src/pages/api_reference.md`

Add sections for exported areas currently missing or underrepresented:
- IO sources:
  `GeosChemSource`, `NetCDFGridSource`, `NetCDFSource`,
  `geoschem_to_dict`, `read_geoschem_profile`.
- SIF helpers:
  `sif_data_path`, `load_sif_spectrum`, `load_ficus_reflectance`,
  `build_sif_source`.
- SolarModel public helpers:
  `planck_spectrum_wn`, `planck_spectrum_wl`,
  `solar_transmission_from_file`, and any other exported SolarModel symbols
  after checking `names(vSmartMOM.SolarModel)`.
- Inelastic/Raman mode types:
  `noRS`, `RRS`, `RRS_plus`, `VS_0to1`, `VS_1to0`, plus traits such as
  `has_inelastic`, `needs_interaction_workspace`, and
  `normalize_raman_weights!` if they are intentionally public.
- Aerosols:
  schemes, data containers, and the main loader/optical-property functions.
- CoreRT Jacobian helpers:
  `ParameterLayout`, `n_total`, `aerosol_range`, `gas_range`,
  `surface_range`, `surface_index`, `canopy_range`, `n_layer_params`,
  `lin_added_layer_all_params!`, and `rt_run_ss` if it remains public.

Open question for Claude:
- Which exported symbols should be demoted from public docs before
  registration instead of documented? The Inelastic helper exports are the
  main candidates.

Acceptance checks:
- Temporarily set `checkdocs = :exports` locally and record remaining missing
  or intentionally skipped public docs.
- `julia --project=docs docs/make.jl` should complete without missing-doc
  surprises once `warnonly` is removed or narrowed.

### 4. Core RT Theory and Code Map

File: `docs/src/pages/vSmartMOM/CoreRTTheory.md`

Add a short "Code Path" section before or after the equations:

```text
rt_run(model)
  -> constructCoreOpticalProperties(...)
  -> extractEffectiveProps(...)
  -> rt_kernel!(...)
      -> elemental!(...)
      -> doubling!(...)
      -> interaction!(...)
  -> create_surface_layer!(...)
  -> postprocessing_vza!(...)
```

Add a "Layer State in Code" section:

```text
AddedLayer:
  r-+, r+-, t++, t--, j0+, j0-

CompositeLayer:
  R-+, R+-, T++, T--, J0+, J0-
```

Explain the convention in plain terms:
- lowercase fields belong to the current homogeneous layer;
- uppercase fields belong to the accumulated composite atmosphere;
- `j0+`/`J0+` and `j0-`/`J0-` are solar source-function terms.

Add exact code snippets, but keep them short:
- `src/CoreRT/CoreKernel/rt_kernel.jl`: the call sequence
  `elemental!`, `doubling!`, `interaction!`.
- `src/CoreRT/CoreKernel/elemental.jl`: where `get_elem_rt!` fills `r-+`
  and `t++`, and where `get_elem_rt_SFI!` fills `J0+` and `J0-`.
- `src/CoreRT/CoreKernel/rt_helpers.jl`: three helper functions:
  `compute_geometric_progression!`, `doubling_source_update!`,
  `doubling_rt_update!`.
- `src/CoreRT/CoreKernel/interaction.jl`:
  `interaction_helper!(::ScatteringInterface_11, ...)`.

Add equations that match those snippets:
- elemental reflection and transmission for finite optical thickness;
- SFI source terms, including stable `-expm1(-x)` and
  `expdiff_neg(a, b)` motivation in the thin limit;
- doubling:
  `G = (I - r r)^-1`, `T G`, source update, and `expk = expk^2`;
- adding:
  the full-scattering update with two inverses:
  `(I - r R)^-1` and `(I - R r)^-1`.

Important wording:
- Say that current SFI is solar-source integration inside the solver path.
- Do not imply the proposed fully offline unified source function exists yet.
- If mentioning thermal emission, label it as planned/future design.

Acceptance checks:
- Every code snippet names the source file path.
- The page tells a reader exactly where to inspect elemental, doubling, and
  adding implementation.
- Equations use the same uppercase/lowercase convention as the code map.

### 5. Tutorial Truth Pass

Files:
- `docs/src/pages/tutorials/*.jl`
- generated `*.md` files after running Documenter/Literate

Changes:
- Replace repo-relative path hacks with robust paths:
  `pkgdir(vSmartMOM)` or `joinpath(pkgdir(vSmartMOM), "test", ...)`.
- Prefer `parameters_from_file` unless a tutorial is specifically about the
  YAML compatibility wrapper.
- Identify tutorials that can be executed in docs CI.
- Leave GPU-only sections guarded and explicit: print why they are skipped.
- Do not depend on missing SIF fixtures until the data policy is resolved.

Acceptance checks:
- Non-GPU tutorials run with `Literate.markdown(...; execute=true)` or via a
  dedicated docs smoke test.
- `rg "dirname\\(pathof\\(vSmartMOM\\)\\)" docs/src/pages/tutorials` returns no
  hits unless justified.

### 6. SIF Documentation/Data Decision

Files:
- `src/SIF_emission/sif_loader.jl`
- `docs/src/pages/api_reference.md`
- `test/test_sif.jl`

Decision required before final docs:
- Option A: commit the small canonical SIF and ficus data files if provenance
  is acceptable.
- Option B: move them to an artifact and document artifact fetch behavior.
- Option C: keep loaders public but document that the default data bundle is
  optional, and gate smoke tests when fixtures are absent.

Recommendation:
- Do not claim "bundled data" until the files are tracked or artifact-backed.

### 7. Documenter Strictness

File: `docs/make.jl`

Staged approach:
1. Keep `warnonly = true` during initial cleanup to collect warnings.
2. Add `checkdocs = :exports` once the API page and docstrings are updated.
3. Remove or narrow `warnonly` before release so broken references fail docs CI.
4. Execute a subset of Literate tutorials once they are stable.

Acceptance checks:
- `julia --project=docs docs/make.jl` passes locally.
- Docs CI uses the same source tree and does not rely on stale generated files.

### 8. Release-Facing Docs

New file candidate:
- `CHANGELOG.md` or `docs/src/pages/release_notes.md`

Must cover:
- vSmartMOM 2.0.0 status.
- Julia compat baseline.
- `parameters_from_file` and TOML support.
- `parameters_from_yaml` behavior for non-YAML extensions.
- CUDA 6 support.
- RAMI removed from unit-test critical path.
- Known optional-data items, especially SIF if unresolved.

## Suggested Commit Sequence

1. `docs: align landing page and IO naming`
   - index page, IO overview/examples/schema, quick-start tutorial.

2. `docs: expand public API reference`
   - API page only, with any docstring-only additions needed outside WIP files.

3. `docs: map CoreRT theory to solver code`
   - CoreRTTheory page, focused on elemental/doubling/adding.

4. `docs: make tutorials path-stable`
   - tutorial path cleanup, no solver changes.

5. `docs: tighten Documenter checks`
   - `checkdocs`, selected Literate execution, docs build validation.

6. `docs: add release notes`
   - changelog/release notes for registration.

## Review Questions for Claude

1. Should `read_parameters` or `parameters_from_file` be the one taught first?
2. Which Inelastic/Raman exports are intentionally public?
3. Should `rt_run_ss` remain public and documented?
4. Should SIF default loaders be public before the fixture policy is resolved?
5. Is the CoreRT theory page too implementation-heavy for user docs, or should
   the code map live there because it is the solver heart?
6. Should `warnonly = true` be removed before registration, or kept with a
   documented list of allowed warnings?

## Claude Prompt

Use this prompt for a second-pass review:

```text
Review docs/dev_notes/docs_change_plan_for_claude.md against the current
vSmartMOM.jl sanghavi-unified codebase. Focus on whether the plan is truthful,
whether it over-documents internals, whether the proposed public API policy is
coherent, and whether any registration-blocking documentation risks are
missing. Do not edit code. Return a prioritized list of changes to the plan.
```

