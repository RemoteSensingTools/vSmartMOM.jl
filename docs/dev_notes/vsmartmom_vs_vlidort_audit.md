# vSmartMOM vs vLIDORT Audit

**Date:** 2026-05-05  
**Scope:** local `vSmartMOM.jl` checkout at `8e267b3` vs sibling `../PyVLIDORT` checkout at `20d97e9`, which vendors `VLIDORT2p8p3` and `LIDORT3p7`.

This is a technical audit, not a diplomacy exercise. vSmartMOM is the better
modern software project. vLIDORT is the more mature radiative-transfer product.
Those are not the same thing.

## Executive Verdict

vLIDORT wins on community credibility, feature coverage, and edge-case physics.
It has the baggage of a long-lived Fortran codebase, but it also has the
validation history and feature knobs that atmospheric-RT users expect: thermal
emission, surface leaving, water-leaving/fluorescence, BRDF/BPDF catalogs,
linearized profile/column/surface weighting functions, pseudo-spherical
geometry, doublet/observation geometry, and a pile of shipped reference outputs.

vSmartMOM wins on architecture, packaging, composability, GPU design, integrated
line-by-line absorption/scattering, documentation quality, and sane user-facing
examples. It is a real Julia package with CI, Documenter/Vitepress docs, weak
CUDA/Metal dependencies, typed polarization, and a coherent optical-property
algebra. The solver design is strong. The missing piece is not elegance; it is
coverage and validation. Right now vSmartMOM looks like a high-potential modern
code that still needs a hardening campaign before the broader community will
treat it as a drop-in alternative to vLIDORT.

The blunt version:

| Area | Winner | Reason |
|---|---|---|
| Core RT maturity | vLIDORT | Two decades of feature accretion, validation cases, and operational use. |
| Modern software design | vSmartMOM | One Julia package, typed APIs, package extensions, docs, CI, GPU path. |
| GPU/hyperspectral design | vSmartMOM | Spectral axis is the batch axis; CUDA/Metal are designed in. |
| Feature completeness | vLIDORT | Thermal, sphericity, surface-leaving, BRDF/BPDF zoo, weighting functions. |
| Raman / Ring-effect modeling | vSmartMOM | Package-native RRS/VRS/RVRS machinery; local vLIDORT tree has fluorescence/SLEAVE, not equivalent atmospheric Raman redistribution. |
| Testing discipline | Split | vSmartMOM has modern automated tests; vLIDORT has better physics reference fixtures. |
| User friendliness | vSmartMOM package, vLIDORT domain workflows | vSmartMOM is easier to install/run; PyVLIDORT is more mission-pipeline-specific. |
| Community trust today | vLIDORT | Users trust old ugly validated code more than clean code without enough benchmark receipts. |

## Evidence Base

vSmartMOM evidence:

- Public package metadata and weak GPU extensions: `Project.toml`.
- CI: `.github/workflows/AutomatedTests.yml`, `.github/workflows/Documentation.yml`.
- Test harness: `test/runtests.jl`.
- Concepts/docs: `docs/src/pages/concepts/*.md`, especially `01_overview.md`, `04_mom_solver.md`, `06_linearization.md`, `07_architecture.md`.
- Release/maturity caveats: `docs/src/pages/release_notes.md`.
- Architecture implementation: `src/Architectures.jl`, `ext/vSmartMOMCUDAExt.jl`, `ext/vSmartMOMMetalExt.jl`.
- Roadmap/validation seed already present: `dev_notes/vlidort_baseline_suite.md`.

vLIDORT / PyVLIDORT evidence:

- Local package root: `../PyVLIDORT`.
- Vendored solver: `../PyVLIDORT/src/Components/rtms/RTSI/VLIDORT2p8p3`.
- vLIDORT version and feature summary: `vlidort_def/vlidort_pars.f90_save`.
- Inputs and flags: `vlidort_def/vlidort_inputs_def.f90`, `vlidort_def/vlidort_lin_inputs_def.f90`.
- Test scripts: `vlidort_run.bash`, `vlidort_check2.bash`, `vlidort_run_SpherCorr.bash`.
- Shipped references: `vlidort_{s,v,sc}_test/saved_results/gfortran/**`.
- Python/f2py wrapper: `src/Components/rtms/py_vlidort/VLIDORT_POLAR_py.F90`, `py_vlidort/vlidort.py`.
- Mission wrapper: `src/Components/missions/SBG/sbg_vlidort_pyexample.py`.

## Approach and Numerical Model

### vLIDORT

vLIDORT is a vectorized linearized discrete-ordinate solver. It solves the
polarized RTE through the DISORT-family eigen/BVP machinery, with a mature
set of correction paths and switches around it. Its strength is not that the
design is beautiful. Its strength is that almost every ugly operational
condition has a flag, a test driver, or a historical patch.

The local vendored solver is `VLIDORT2p8p3`, released 2021-03-31. Its own
header lists recent capabilities: total-column Jacobians, BPDF kernels,
thermal emission, consolidated BRDF handling, surface-leaving/BRDF scaling,
Taylor-series fixes, OpenMP thread safety, water-leaving, planetary/media
properties, doublet geometry, dynamic memory, Green's-function solutions,
sphericity corrections using multiple-scattering source terms, snow BRDF,
and extended surface-leaving.

That is the feature envelope vSmartMOM is being judged against.

### vSmartMOM

vSmartMOM uses the Matrix Operator Method: elemental single-scattering layer,
doubling to full layer optical depth, and adding to compose layers and
surface. Its design is much cleaner for modern retrieval workloads:

- layer optics reduce to `(τ, ϖ, Z⁺⁺, Z⁻⁺)`;
- spectral dimension is the third array axis, so batched matrix operations run
  over wavelengths;
- CPU, CUDA, and Metal are selected by architecture/array type dispatch;
- polarization is a type, not a runtime flag;
- linearization is intended as operator-level analytic chain rule, with AD only
  upstream at optical-property boundaries.

Two vSmartMOM design choices are genuinely strong:

- The elemental layer uses exact finite-δ single-scatter formulas, not only the
  infinitesimal linear limit.
- Doubling count is sized from scattering optical depth rather than total gas
  optical depth, which keeps line-by-line spectral batches uniform.

These are not cosmetic. They are the kind of choices that make hyperspectral
GPU radiative transfer plausible.

## Technical Maturity

| Criterion | vSmartMOM | vLIDORT / PyVLIDORT |
|---|---|---|
| Solver history | Young compared with vLIDORT; JOSS paper and active Julia refactor line. | Old, widely recognized, explicitly versioned solver lineage from 2003-2021 in local source. |
| Feature envelope | Strong elastic solar polarized RT, absorption, Mie/Rayleigh, surfaces, newer RRS/VRS/RVRS Raman forward path, analytic elastic Jacobian path. Missing or incomplete thermal, sphericity, some surface-leaving/ocean coupling, full inelastic linearization. | Huge envelope: solar, thermal, surface emission, surface leaving, water leaving, fluorescence, pseudo-spherical, observation/doublet geometry, BRDF/BPDF catalog, profile/column/surface linearization. |
| Numerical validation | Has unit/regression tests and some references, but not yet enough public gold-standard comparisons against vLIDORT/Siewert/MODTRAN/6SV. | Ships many saved results and dedicated tester programs. Some are old-school diff files, but they are real reference artifacts. |
| Performance story | Best-in-class direction: batched spectral axis, CUDA/Metal, GPU Mie/absorption paths, weak GPU deps. | CPU/OpenMP-oriented. Py wrapper loops channels/observations and uses multiprocessing. Solid but not modern GPU-native. |
| API stability | v2.0 registration-oriented cleanup in progress. Some public/internal boundary is still loose. | Fortran API is stable by inertia. Usability suffers, but operational callers know what they are getting. |
| Operational readiness | Not yet a drop-in replacement for vLIDORT in conservative production pipelines. | Solver is operationally credible; local PyVLIDORT repo packaging is rough. |

Scorecard, deliberately harsh:

| Area | vSmartMOM | vLIDORT solver | Local PyVLIDORT wrapper |
|---|---:|---:|---:|
| Solver physics coverage | 6/10 | 9/10 | N/A |
| Solver validation depth | 5/10 | 9/10 | 4/10 |
| Modern package quality | 8/10 | 3/10 | 4/10 |
| GPU/hyperspectral design | 9/10 | 2/10 | 2/10 |
| Raman / inelastic atmospheric scattering | 8/10 forward, 4/10 linearized | 3/10 | 2/10 |
| Retrieval Jacobian ergonomics | 7/10 | 8/10 | 4/10 |
| New-user experience | 7/10 | 3/10 | 4/10 |
| Community confidence today | 5/10 | 9/10 | 5/10 |

## Testing Framework

### vSmartMOM Testing

vSmartMOM has a recognizably modern Julia test suite:

- `test/runtests.jl` centralizes the test sets.
- Tests cover absorption, scattering, CoreRT, batched kernels, forward no-Raman
  runs, Jacobian units, type stability, Float32 consistency, quality gates,
  IO, parser validation, canopy, truncation, Cox-Munk, aerosol wiring, Raman
  phase/regression, and script-port parsing.
- CI runs package tests on Ubuntu, Windows, and macOS for Julia 1.11.2.
- Documentation builds in CI through Documenter/Vitepress.

This is good software discipline. The weak spots are also clear:

- GPU tests are conditional and not a hard CI gate.
- Some important end-to-end workflows are disabled because fixtures are missing
  or machine-local.
- The gold-standard comparison layer is still mostly a plan, not a merged
  enforcement mechanism.
- Performance tests exist, but they are not yet an operational dashboard with
  regression thresholds.

vSmartMOM's test framework is clean but underfed. It needs more authoritative
physics fixtures.

### vLIDORT Testing

vLIDORT's testing is old-school but domain-serious:

- `vlidort_run.bash` compiles and runs scalar/vector tester programs.
- `vlidort_check2.bash` compares generated outputs against
  `saved_results/$compiler`.
- Separate drivers exist for solar, thermal, BRDF, surface-leaving,
  Siewert-2000 validation, planetary, longwave coupling, OpenMP, and
  spherical-correction cases.
- The distribution ships reference output files for scalar, vector, Stokes-3,
  observation geometry, doublet geometry, and sphericity tests.

The problems:

- The local `../PyVLIDORT` repo has GitHub workflows for labels/changelog/YAML,
  not build/test execution.
- There is no visible Python test suite for the wrapper.
- The root README is basically a template stub.
- The test scripts mutate active makefiles and `vlidort_pars.f90`, which is
  exactly the kind of stateful workflow modern CI hates.

So the split is blunt: vLIDORT has the better physics fixtures; vSmartMOM has
the better testing infrastructure. vSmartMOM should steal the fixtures, not the
workflow.

## Software Design

### vSmartMOM Strengths

vSmartMOM is designed like a modern scientific package:

- Structured Julia modules with public docs and tutorials.
- Type-driven polarization and architecture dispatch.
- Optional CUDA/Metal package extensions instead of hard GPU dependencies.
- Layer optical properties compose with algebraic operations.
- CPU/GPU/Metal share kernel source where practical.
- Absorption, scattering, surfaces, and RT live in one language.
- Documentation explains the theory and points back to code anchors.

The architecture is credible. The main design debt is API boundary cleanup:
some internals are still exported, some linearization handoff structures are
still being stabilized, and several experimental modules are visible enough to
confuse users.

### vLIDORT Strengths

vLIDORT is a feature machine:

- Explicit input structures for dozens of RT modes.
- Long-standing linearized structures for column/profile/surface Jacobians.
- BRDF, SLEAVE, and single-scatter correction supplements.
- Many special-purpose geometries and correction modes.
- Static dimensioning makes old Fortran builds predictable.

But the design is not friendly:

- Important limits live in parameter files like `MAXLAYERS`, `MAXSTREAMS`, and
  `MAX_USER_VZANGLES`.
- Feature activation often means juggling flags whose valid combinations are
  non-obvious.
- Build and tests rely on makefile/script conventions and source-file swaps.
- The Python wrapper exposes mission-specific shapes and outputs only `I,Q,U`
  in the common paths.
- The local wrapper is not packaged as a normal Python project.

vLIDORT's design optimizes continuity with legacy atmospheric RT practice.
vSmartMOM optimizes maintainability and hardware portability. For new
development, vSmartMOM is the better base.

## User Friendliness

vSmartMOM is easier for a new technical user:

- `Pkg.add`, `using vSmartMOM`, `read_parameters`, `model_from_parameters`,
  `rt_run`.
- YAML scenes with defaults and comments.
- Tutorials for quick start, scattering, absorption, surfaces, GPU, IO, and
  Jacobians.
- Package docs are built and published.

But vSmartMOM still has friction:

- Many examples assume Julia comfort.
- Some scenes require knowing internal conventions like wavenumber vs wavelength
  and relative azimuth convention.
- Error messages and schema validation need to be ruthless, because RT users
  will misconfigure scenes constantly.
- Python users have no first-class low-level package entry point.

vLIDORT/PyVLIDORT is easier for a very specific class of user: someone already
inside the GEOS/SBG/pyobs ecosystem. For everyone else, the root README and
packaging are poor. The Fortran user guide may be extensive, but a `.docx`
manual and shell scripts are not a modern onboarding path.

## Licensing and Redistribution

vSmartMOM is Apache-2.0 at the package level. That is straightforward for
academic, government, and commercial downstream use.

The local `../PyVLIDORT` repository has an Apache-2.0 root `LICENSE`, but the
vendored `VLIDORT2p8p3` solver itself declares GPL-3.0 in its bundled license
statement and source headers. That is not a small footnote. Anyone distributing
PyVLIDORT-linked binaries or derivative wrappers needs to treat the solver as a
GPL component unless separate licensing exists elsewhere. For community adoption,
vSmartMOM's license is cleaner.

## Raman and Inelastic Scattering

This is the major place where the first draft under-credited vSmartMOM.
vLIDORT has mature surface-leaving and fluorescence support, but fluorescence
is an additive surface/source term. It is not the same thing as atmospheric
rotational/vibrational Raman redistribution through a multiple-scattering
solver. In the local `../PyVLIDORT` tree, I found fluorescence hooks in
`vsleave`, but no equivalent package-native RRS/VRS atmospheric Raman path.

vSmartMOM's inelastic work is newer and more interesting:

- `RRS`, `VS_0to1`, `VS_1to0`, and `_plus` modes are dispatch types under
  `src/Inelastic/types.jl`.
- Raman/Cabannes optical properties are precomputed in `src/Inelastic/`.
- Inelastic source terms propagate through parallel elemental, doubling, and
  interaction kernels: `elemental_inelastic*.jl`, `doubling_inelastic.jl`,
  `interaction_inelastic.jl`.
- The implementation follows the linear-in-inelastic-scattering approximation:
  one Raman event per photon path, with multiple elastic scattering before and
  after. That is the right practical approximation for Ring-effect/O2-band work.
- `rt_run_ss` is not just a debug helper; it is the single-scatter inelastic
  correction path described in the docs for fast O2 A-band Raman corrections.
- Tests include Raman smoke/physics checks, GPU-conditional Raman forward tests,
  and a Phase 1b RRS regression against a frozen Sanghavi reference
  (`test/test_forward_raman_phase1b.jl`).

The caveats are real:

- Raman memory scales badly because redistribution couples wavelengths; the
  tests explicitly call out `O(nλ^2)` memory pressure.
- GPU Raman exists, but it is conditional and not a hard CI gate.
- Inelastic linearization is not complete; the mature Jacobian path is still
  elastic/noRS.
- SIF helpers and data policy are not settled, so Raman/SIF workflows should
  not be marketed as fully productized yet.

Bottom line: for Ring-effect and atmospheric Raman physics, vSmartMOM is not
behind vLIDORT; it is ahead of the local vLIDORT target. The catch-up work is
not "add Raman." It is "harden Raman": reduce memory pressure, validate against
published/Sanghavi cases, make GPU coverage a real gate, and finish or clearly
scope inelastic Jacobians.

## Feature Gap: What vSmartMOM Is Missing

These are the gaps that matter for community adoption, not academic elegance.

| Gap | Why it matters | vLIDORT evidence | vSmartMOM status |
|---|---|---|---|
| Gold-standard validation suite | Nobody trusts a new RT solver without benchmark receipts. | Siewert, solar/thermal/BRDF/sphericity saved results. | Design note exists; not yet a hard CI gate. |
| Pseudo-spherical and sphericity corrections | High-SZA and limb-ish geometries are common enough to matter. | `DO_FOCORR_*`, refractive geometry, Chapman, MSST sphericity. | Mostly plane-parallel. |
| First-order/single-scatter correction modes | Community expects FO correction comparisons against VLIDORT. | FO code and correction flags. | Exact SS and correction design work exists, not complete as product. |
| Thermal emission and surface emission | Needed for IR, thermal, and longwave coupling use cases. | Thermal tester, thermal flags, LBBF Jacobians. | Release notes say thermal is a design topic, not user-facing. |
| Surface leaving, water-leaving, fluorescence | Ocean color and SIF workflows need this. | SLEAVE, water-leaving, fluorescence flags and tests. | Raman/SIF pieces exist, data/workflow policy incomplete; ocean coupling future. |
| Raman hardening, not Raman existence | Ring-effect users need validated and affordable RRS/VRS, not only a working research path. | Local vLIDORT tree does not appear to provide equivalent atmospheric Raman redistribution. | vSmartMOM is ahead here, but needs memory/performance work, public examples, and inelastic Jacobian scope. |
| BRDF/BPDF catalog completeness | Users compare by named kernels. | 19 BRDF indices: Hapke, Roujean, BPDF, Cox-Munk variants, hotspot, modified Fresnel, snow. | Lambertian, RPV, Ross-Li, Cox-Munk, canopy. Missing several expected named kernels. |
| Linearized surface/profile completeness | Retrieval users need weighting functions, not only radiances. | Column/profile/surface/SLEAVE weighting functions. | Analytic framework exists; some parameters and inelastic paths incomplete. |
| Observation/doublet geometry modes | Operational products often use compact observation-geometry inputs. | Dedicated flags and saved results. | Geometry handling is simpler and less feature-complete. |
| Python-facing API | Much of the atmospheric community works in Python. | PyVLIDORT exists, mission-specific. | No first-class Python package/wrapper. |
| Operational examples | Users want "run this scene and compare to known output". | Test drivers and saved outputs. | Tutorials are better, but fewer gold references. |

## Roadmap: Fast Catch-Up Plan

The fastest path is not to reimplement every vLIDORT feature in random order.
The fastest path is to close trust gaps first, then close feature gaps in the
order users will actually ask for them.

### Phase 0: Stop Guessing, Build the vLIDORT Baseline Gate

**Timeline:** 1 week.

Implement the already-planned `dev_notes/vlidort_baseline_suite.md` as a real
test asset:

- Parse shipped VLIDORT saved results.
- Add Siewert 2000 Problem IIA as a committed vSmartMOM comparison.
- Add scalar/vector `solar_tester` Task 1 and Task 2 comparisons.
- Put fixtures under `test/vlidort_baseline/`.
- Add a normal Julia test set, skipped only when explicitly marked expensive.
- Generate a Markdown/JSON report with max error, RMS error, tolerances, and
  convention notes.

This is the highest-leverage move. Without it, every feature claim sounds like
"trust us". With it, the discussion changes to numbers.

### Phase 1: Add a Raw Optical-Properties Scene API

**Timeline:** 1-2 weeks.

The baseline work will be painful unless vSmartMOM can ingest layer optical
properties directly. Add a documented low-level API for:

- layer `τ`, `ϖ`, Greek/F-matrix moments or `Z⁺⁺/Z⁻⁺`;
- geometry;
- surface;
- truncation and Fourier controls;
- optional Jacobian perturbations.

This should bypass HITRAN/Mie when the user already has optics. vLIDORT users
think in optical property tables; vSmartMOM needs to meet them there.

Deliverable: `LayerOpticsScene` or equivalent, with one tutorial:
"Run vSmartMOM from VLIDORT-style layer optics."

### Phase 2: First-Order Correction and Sphericity Track

**Timeline:** 3-5 weeks.

Implement and productize the exact-SS / FO correction path described in the
local dev notes:

- Standalone exact single-scatter solver with clear path decomposition.
- Back-correction adapter for full-MOM results.
- Plane-parallel FO comparisons against custom PyVLIDORT runs.
- Regular pseudo-spherical solar-path correction as the first sphericity mode.
- Document the scope: what matches VLIDORT, what deliberately does not.

Do not try to land enhanced outgoing sphericity first. Get the regular case
validated, then expand.

### Phase 3: Surface Feature Parity That Users Recognize

**Timeline:** 4-6 weeks.

Close the named-kernel gap:

- Hapke.
- Roujean.
- BPDF soil/vegetation/NDVI.
- Ross-thick hotspot.
- Modified Fresnel.
- Snow BRDF.
- Surface linearization for every shipped surface except canopy if necessary.

The implementation burden is moderate compared with solver work, and the user
payoff is high because feature checklists mention these names.

### Phase 4: Thermal and Longwave Minimum Product

**Timeline:** 6-8 weeks.

Add a minimal but validated thermal path:

- Planck/source term layer handling.
- Surface emission.
- Thermal-only/transmittance-only mode if it fits the architecture.
- Thermal Jacobian targets: temperature/profile and surface emissivity after
  forward mode is stable.
- Compare against vLIDORT thermal tester fixtures.

This is a serious feature, not a weekend patch. But without thermal, vSmartMOM
will always look solar-only to part of the community.

### Phase 5: Harden Raman, Then Settle Surface-Leaving/SIF Policy

**Timeline:** 4-8 weeks, can overlap with Phase 4 if ownership is separate.

Do not frame Raman as a missing vLIDORT-parity item. Frame it as a vSmartMOM
lead that needs production hardening:

- Add a public Raman tutorial for O2 A-band/Ring-effect with a committed
  reference output.
- Add a small always-on Raman regression case that runs in normal CI, plus a
  larger scheduled/optional benchmark.
- Reduce or tile the `O(nλ^2)` wavelength-coupling memory pressure.
- Make GPU Raman a scheduled gate where hardware exists.
- Decide whether inelastic Jacobians are in scope for the near-term release;
  if not, document the forward-only boundary clearly.

Then turn the adjacent experimental source-term pieces into product:

- Decide fixture/artifact policy for SIF data.
- Add water-leaving and fluorescence examples with committed reference output.
- Add SLEAVE-style weighting function coverage where feasible.
- Validate against vLIDORT SLEAVE/water/fluor saved results or custom runs.

This matters for ocean color and vegetation/SIF users. Half-finished SIF helpers
are worse than no SIF story because they create false confidence.

### Phase 6: Python Entry Point

**Timeline:** 2-4 weeks for a minimal useful bridge.

Do not build another mission-specific wrapper first. Build a boring low-level
Python interface:

- Accept numpy arrays for layer optics and geometry.
- Return `R`, `T`, and optional Jacobians.
- Expose a simple scene-file runner.
- Ship wheels only after the Julia package and artifact story are stable; until
  then document `juliacall` usage.

This is adoption work, not solver work. It is still necessary.

### Phase 7: CI and Release Hardening

**Timeline:** continuous; first useful version in 1 week.

Minimum gates:

- Keep normal Julia test matrix.
- Add vLIDORT baseline tests on Linux.
- Add a scheduled CUDA job if hardware is available; otherwise keep a manually
  triggered GPU workflow.
- Track benchmark numbers for representative CPU/GPU line-by-line scenes.
- Fail docs on missing exported docs, as already intended.
- Clean public/internal API exports after v2.0 registration.

## What Not To Copy From vLIDORT

Do not copy the Fortran workflow:

- Do not use source-file swaps to change dimensions or modes.
- Do not expose dozens of boolean flags without typed configurations.
- Do not make users learn compiler-specific test scripts.
- Do not bury the main onboarding path in a document-office manual.
- Do not make Python support mission-specific before the low-level API exists.

Copy the validation assets and the feature semantics. Leave the ergonomics
behind.

## Practical Priority List

If there are only three months, do this:

1. vLIDORT baseline suite in CI.
2. Raw optical-properties API.
3. FO/back-correction and regular pseudo-spherical validation.
4. Raman hardening: public O2/Ring tutorial, small CI reference, memory plan.
5. Missing BRDF/BPDF kernels that are cheap to add.
6. Minimal thermal forward path.
7. Python low-level wrapper.

If there is only one month, do this:

1. vLIDORT baseline suite.
2. Raw optical-properties API.
3. FO/back-correction MVP.
4. One "vSmartMOM vs vLIDORT" tutorial with numbers.
5. One "vSmartMOM Raman/Ring effect" tutorial with a committed reference.

The community will forgive missing features if the roadmap is honest and the
benchmarks are reproducible. It will not forgive vague claims of correctness.

## Bottom Line

vSmartMOM does not need to become vLIDORT in Julia. It needs to become the
modern GPU-native RT package that can prove, case by case, where it agrees with
vLIDORT and where it intentionally differs.

Today, vSmartMOM is architecturally ahead, Raman-ahead, and validation/feature
parity behind everywhere else that conservative RT users care about. The
catch-up path is straightforward: turn vLIDORT's shipped references into
vSmartMOM's regression gate, add the optical-property ingestion layer, harden
the Raman lead with public benchmarks, then prioritize sphericity, FO
correction, thermal, surface-leaving, and named BRDF parity. That is the
shortest path from "promising" to "credible".
