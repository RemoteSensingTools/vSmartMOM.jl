# vSmartMOM.jl — Documentation Story Scaffold

Narrative IA for the docs overhaul. Not prose — this is the through-line and
page-by-page job-to-be-done, organized so Codex can turn it into concrete
file edits in `docs/src/`. Companion to:

- [docs_change_plan_for_claude.md](docs_change_plan_for_claude.md) — the plan
- [docs_change_plan_audit.md](docs_change_plan_audit.md) — Claude's audit + Codex agreement
- [theory_references.md](theory_references.md) — paper↔code crib sheet

## Premise

The current docs are organized by **source-tree module** (vSmartMOM,
Absorption, Scattering). That is fine as a *reference* layer. It fails as a
*journey* layer because:

- A new user lands on a module page that doesn't tell them what to do first.
- Tutorials use YAML before YAML is introduced, and use `pathof(vSmartMOM)`
  joinpath hacks that look like a hostile UI.
- `model_from_parameters_lin` / `LinMode` / `rt_run_lin` — the Jacobian path,
  which is half the reason the package exists per S2014 — is invisible above
  the API reference.
- Inelastic, Aerosols, SolarModel, SIF have no doc pages at all even though
  they have public exports.
- Theory is referenced in nav (`CoreRTTheory.md`) but disconnected from the
  solver source code, so a reader cannot follow an equation into a file.

Fix: **add a task-based journey layer on top of the existing module layer**,
and make every page state who it is for and where to go next.

## Personas

| # | Persona | Wants | Does not want |
|---|---|---|---|
| P1 | **First-time user** (grad student, atmospheric scientist) | Get a plotted reflectance in 5 minutes from `add vSmartMOM` | YAML schema, type hierarchies, theory |
| P2 | **Scene configurator** | Define their own atmosphere/surface/spectral grid in YAML or a Dict | Editing Julia source |
| P3 | **Retrieval / inversion developer** | Compute Jacobians, understand parameter ordering, plug into Optim/inversion code | Forward-only documentation |
| P4 | **Method developer** | Add a new surface BRDF / Raman mode / aerosol scheme; understand kernel dispatch | A black box |
| P5 | **Theorist / benchmarker** | Map equations in the paper to code; reproduce S2014 benchmarks; cite correctly | Hand-waving prose |
| P6 | **GPU/HPC user** | Run on CUDA, understand what is GPU-safe, see speedups | CPU-only assumptions |

P1 and P3 are the volume users. P5 is the credibility user (without them the
package looks like undocumented research code). P4 is the future-contributor
pipeline.

## Journey arc

The docs should move a reader along this arc, *in this order*, with each
page knowing which step it serves:

```
1. Why this exists           (60 seconds)         — landing page
2. Smallest correct example  (5 minutes)          — Quick Start tutorial
3. Configure your scene      (20 minutes)         — IO Overview + Schema
4. Understand the physics    (1 hour, optional)   — CoreRT Theory page
5. Get sensitivities         (30 minutes)         — Linearized RT tutorial
6. Extend the model          (open-ended)         — Surface / Raman / GPU guides
7. Validate and cite         (as needed)          — Benchmarks + References + CITATION.bib
```

Each step has a corresponding doc artifact. Steps 1–3 are mandatory reading
for new users; 4–7 are reference / specialist material.

## Information architecture (proposed)

Replace the current `docs/make.jl` `pages` block with a two-layer
organization:

### Top layer — "How do I…" (task-based)

```
Getting Started           index.md
Quick Start (5 min)       pages/quickstart.md         (was Tutorial_QuickStart)
Configure a Scene         pages/IO/Overview.md        (rewritten as journey, not reference)
Compute Jacobians         pages/jacobians.md          (new — was buried in Tutorial_Jacobians)
Run on GPU                pages/gpu.md                (new — extracted from Tutorial_GPU)
Add a Surface BRDF        pages/extending/surfaces.md (new)
Add a Raman Mode          pages/extending/raman.md    (new — depends on inelastic page)
```

### Reference layer — "What is X" (module-based, current organization)

```
Architecture & Design     design.md                   (rewrite for RTModel hierarchy)
Core RT Theory            pages/vSmartMOM/CoreRTTheory.md  (S2014 + SF2023-II equations + file:LINE refs)
API Reference             pages/api_reference.md
Modules
  CoreRT                  pages/vSmartMOM/Overview.md
  Absorption              pages/Absorption/Overview.md
  Scattering              pages/Scattering/Overview.md
  Inelastic               pages/Inelastic/Overview.md       (NEW)
  Aerosols                pages/Aerosols/Overview.md        (NEW, marked WIP)
  SolarModel              pages/SolarModel/Overview.md      (NEW)
  Surfaces                pages/Surfaces/Overview.md        (NEW — extracted from API ref)
  IO                      pages/IO/Schema.md
HITRAN Data               pages/Absorption/HITRAN_Data.md
GEOSChem Integration      pages/geoschem_integration.md
Tutorials (full)          pages/tutorials/*
References                pages/vSmartMOM/References.md     (rewritten with theory_references.md)
Release Notes             pages/release_notes.md            (NEW)
```

### Tutorials become "long-form" — not the front door

The current tutorials are the only place a user can land that runs code, so
they have implicitly become the quick-start. That is wrong. Promote a
genuine `quickstart.md` (curated, fast, tested) to the front door, and
demote the existing tutorials to "long-form, run-locally" examples.

## Page-by-page job-to-be-done

For each page: **who** it serves, **what** it must accomplish, **where the
reader arrived from**, **where they should go next**, and **what it must
NOT do**.

### `index.md` — Landing

- **Who:** P1, P5
- **Job:** In 60 seconds, convince a reader that this package solves
  polarized atmospheric RT with Jacobians and Raman support, and show them
  one runnable line.
- **Arrives from:** GitHub, JuliaHub, JOSS paper
- **Goes next:** Quick Start (P1), Citation (P5), Modules (P3/P4)
- **Must NOT:** introduce YAML schema, type hierarchies, GPU detail,
  inelastic theory. One worked example, three "next-step" links, done.
- **Acceptance:** the example must literally run as written. Use either
  `default_parameters()` or a shipped `config/quickstart.yaml` whose
  default architecture is CPU.

### `quickstart.md` — 5-minute tutorial

- **Who:** P1
- **Job:** First successful R, T plot. Single spectral band, Lambertian
  surface, no aerosol, no Raman, CPU.
- **Arrives from:** index.md
- **Goes next:** "Configure a Scene" if they want to change anything;
  "Compute Jacobians" if they're a retrieval person.
- **Must NOT:** explain the matrix operator method, parameter layout,
  polarization Stokes basis. Defer to theory page.
- **Acceptance:** runnable in CI in <60s on CPU.

### `pages/IO/Overview.md` — Configure a Scene

- **Who:** P2
- **Job:** Teach the user how to describe an atmosphere + surface + geometry
  + spectral grid as a YAML file or in-memory Dict. Resolve the
  `read_parameters` / `parameters_from_file` / `parameters_from_yaml`
  inconsistency by **picking one for teaching** (recommend
  `read_parameters` because it dispatches on String/Dict/IOSource and
  matches what the page already teaches), and naming the others as
  explicit aliases.
- **Arrives from:** quickstart.md (when the user says "but I want to change
  the surface" or similar)
- **Goes next:** Schema page for the full grammar; Surfaces page for BRDF
  options; Theory page for what `max_m` and `l_trunc` do.
- **Must NOT:** deprecate `parameters_from_yaml`. Per Codex, that is an API
  decision; in docs we present `read_parameters` as preferred and
  `parameters_from_yaml` as supported YAML-specific.
- **Acceptance:** every code block runs with a config that ships in `config/`.

### `pages/IO/Schema.md` — Reference grammar

- **Who:** P2 (lookup mode)
- **Job:** Be the authoritative grammar of the YAML/TOML/Dict schema.
- **Acceptance:** every key documented matches what
  `parameters_from_dict()` actually accepts in `src/IO/Parameters.jl`.

### `pages/jacobians.md` — Compute Jacobians (NEW)

- **Who:** P3
- **Job:** Show the linearized path: `LinMode()` →
  `model_from_parameters_lin` → `rt_run_lin` returning `(R, T, dR, dT)`.
  Explain `ParameterLayout` and how to slice `dR` for a specific aerosol /
  gas / surface parameter.
- **Arrives from:** quickstart.md, index.md, Tutorial_Jacobians (which
  becomes the long-form companion)
- **Goes next:** Theory page §Linearization for derivation; ParameterLayout
  reference for index arithmetic; HybridAD tutorial for ForwardDiff overlap.
- **Must NOT:** re-derive S2014 Appendix C. Point at it.
- **Acceptance:** worked example produces dR for known parameters and a FD
  cross-check.

### `pages/vSmartMOM/CoreRTTheory.md` — Core RT Theory

- **Who:** P5, advanced P3, P4
- **Job:** Map every equation in the elastic + inelastic solver to a
  source file and line range, citing S2014 / SF2023-II / S2015. Use
  [theory_references.md](theory_references.md) as the source of truth for
  the mapping.
- **Arrives from:** any page that uses a phrase like "matrix operator
  method," "elemental layer," "doubling," "Cabannes vs RRS."
- **Goes next:** API reference for the function signatures; Extending pages
  for how to add a kernel dispatch arm.
- **Sections (in order):**
  1. Vector RTE (S2014 Eq. 4)
  2. Elemental layer (S2014 Eqs. 19–22) → `elemental.jl`
  3. Doubling and adding (S2014 Eqs. 23–32) → `doubling.jl`, `interaction.jl`
  4. Why no solar SFI in J (S2014 §2.2) — design-decision callout
  5. Truncation (S2014 App. A, S2015) → `LayerOpticalProperties/`
  6. Quadrature (S2014 App. B) — block-Radau direct raytracing
  7. Linearization (S2014 App. C) → `_lin.jl` family
  8. Inelastic extension (SF2023-II §3) → `*_inelastic.jl`
  9. Constant-N_doubl elemental trick (SF2023-II §3.2) — design callout
  10. Single-scattering inelastic correction (SF2023-II Eq. 32) — explains `rt_run_ss`
  11. Other dispatch arms — `_canopy`, `_multisensor`, `_hdrf`, `_ss` map
- **Must NOT:** invent equations or call paths. Every block cites a paper
  equation AND a `file.jl:LINE` range.
- **Acceptance:** every equation block has a paper citation and a file link;
  every file link points at code that compiles.

### `pages/Inelastic/Overview.md` — Inelastic Scattering (NEW)

- **Who:** P3 (UV/Vis retrieval), P4 (extending)
- **Job:** Explain noRS / RRS / VS_0to1 / VS_1to0 (and `_plus` variants)
  as user choices. When to use which. Reference S2022-I tables.
- **Arrives from:** Configure a Scene (when user picks an inelastic mode);
  Theory page; Ring-effect / O2 A-band examples.
- **Goes next:** Theory page §inelastic; SF2023-II for full RT derivation.
- **Must NOT:** document every internal trait (`has_inelastic`,
  `needs_interaction_workspace`, etc.) as public API. Decide per export.

### `pages/Aerosols/Overview.md` — Aerosols (NEW, WIP-flagged)

- **Who:** P2, P3
- **Job:** Document `TOMAS15Scheme`, `TwoMomentScheme`, the refractive
  index database, and the `read_aerosol_data` loader. Mark as WIP per
  the source header.
- **Must:** state explicitly that the API may change.

### `pages/SolarModel/Overview.md` — Solar irradiance (NEW)

- **Who:** P2, P3
- **Job:** `planck_spectrum_wn/wl`, solar transmission file loader,
  artifact behavior.

### `pages/Surfaces/Overview.md` — Surface BRDF (NEW, extracted)

- **Who:** P2, P4
- **Job:** Currently the BRDF types are buried in API reference but they
  are the user-facing config strings (BRDF_MAP keys). Promote to a real
  page with a row per type: parameters, when to use, citation.
- **Acceptance:** every type in `src/IO/Parameters.jl` `BRDF_MAP` is
  documented here.

### `pages/extending/surfaces.md` — Add a Surface BRDF (NEW)

- **Who:** P4
- **Job:** Cookbook: subtype `AbstractSurfaceType`, implement
  `create_surface_layer!()`, register in `BRDF_MAP`, export, write a YAML
  example, write a test.
- **Source:** translates the `## Adding a New Surface Model` workflow in
  CLAUDE.md to user-facing prose.

### `pages/extending/raman.md` — Add a Raman Mode (NEW)

- **Who:** P4 (rare)
- **Job:** Adding to InelasticScattering: subtype `AbstractRamanType`,
  define cross-sections per S2022-I, implement the `noRS`/RRS-style
  dispatch hooks. Pointer to SF2023-II for the RT side.

### `pages/gpu.md` — Run on GPU (NEW, extracted)

- **Who:** P6
- **Job:** CUDA setup, weak-dependency mechanism (`vSmartMOMCUDAExt`),
  `array_type(model)`, what is and is not GPU-safe (Mie is CPU per
  Tutorial_GPU comments), measured speedups.
- **Arrives from:** quickstart, jacobians (when user wants speed)

### `pages/release_notes.md` — Release notes / migration (NEW)

- **Who:** anyone upgrading from v1
- **Job:** Per Claude's audit, document what broke (vSmartMOM_Model
  removal, RTModel sub-struct refactor, accessor functions, CUDA 6
  baseline, `read_parameters` policy, optional SIF data status).

## Cross-cutting story rules

1. **Every page header has a "For" line** stating the persona(s) and a
   "Next" line listing 1–3 onward links. Without these, the docs read as a
   bag of pages instead of a path.

2. **Every code block runs.** `docs/test_examples.jl` includes them under
   `@testset`. CI fails if a block stops running.

3. **Every theory equation has both a paper citation and a `file.jl:LINE`
   link.** No equation without provenance.

4. **No `pathof(vSmartMOM)` joinpath hacks in user-facing code.** Configs
   live in `config/` and are referenced by package-relative paths via a
   helper, or the user's own path.

5. **No deprecation in this pass.** Per Codex: present `read_parameters` as
   preferred, keep `parameters_from_yaml` as supported. API change is a
   separate decision.

6. **WIP areas wear WIP labels** — Aerosols module, anything still being
   sanghavi-merged, gets an "API may evolve" callout. No silent shipping of
   provisional surface area as if it were stable.

7. **"Not yet supported" is a feature.** Explicitly call out: thermal
   emission, fully-offline SFI, vector linearization with Raman (per
   `RAMAN_CODE_HANDOFF.md` note that linearized path is elastic-only).
   Refusing to claim things we don't have is half the truthfulness fight.

## Through-line summary

A new user reads **index.md → quickstart.md → IO Overview → Configure a
Scene** in 30 minutes and has a working scene. A retrieval developer
branches at quickstart into **Jacobians → Theory App. C → ParameterLayout
reference**. A theorist branches into **Theory page → References →
CITATION.bib**. A method developer branches into **Theory page → Extending
guides**. Every branch is signposted at every junction.

The current docs do not signpost. That is the single biggest fix.

## What this scaffold leaves to Codex

- Concrete diff per file (this scaffold names pages and roles, not bytes).
- Wiring `docs/make.jl` `pages` block to the new IA.
- Building `docs/test_examples.jl` as the runtime gate.
- Resolving SIF default-loader issue at the code level (Codex's call: ship
  fixtures, switch to artifact, or add a clear "data not available" error).
- Two-step Documenter strictness rollout per the audit.

## What this scaffold leaves to product / human

- SIF data policy: ship fixtures vs artifact vs gate the loaders. Until
  this is decided, Inelastic/Solar pages can ship but SIF page cannot
  honestly ship.
- Whether `parameters_from_yaml` ever gets `@deprecate`d. Codex says no for
  this pass; agreed.
- Whether `vSmartMOM_Lin` gets renamed for naming consistency with
  `RTModel`. Independent of docs.
- Whether to drop `pages/vSmartMOM/Principles.md`, `Types.md`, etc. or
  rewrite them — they need a read-through to decide if they have material
  the new IA hasn't absorbed.
