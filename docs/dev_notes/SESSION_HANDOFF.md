# Resume v2.2 integration and Jacobian work

Saved 2026-09-08 at the user's request before logging off. This is the current
session handoff; older dated audits retain their historical findings.

## Start here

- **Working checkout:** `/home/cfranken/code/gitHub/vSmartMOM-release`
- **Branch:** `integration/release-candidate`
- **Last implementation/documentation commits before this handoff:**
  - `3c1add99` — execute documentation examples and record deep audit (cfranken).
  - `7ae92577` — align scientific guides with supported v2.2 contracts (Suniti).
  - `326c1e91` — preceding integrated-release validation checkpoint.
- Package version is **2.2.0, an unpublished local release candidate**.
- No push, tag, release publication, registration, or external message was sent.
- The user asked to preserve context and resume further changes after returning.
  No additional scientific implementation is running in the background.

Read `AGENTS.md`, `CLAUDE.md`, then:

1. [Deep documentation audit](documentation_audit_2026-09-08.md): latest findings,
   what was fixed, what remains, coverage limits and validation.
2. [Release validation](release_validation_2026-09-07.md): full integration's
   numerical checks, exact CPU caveat and outstanding publication gates.
3. [Release consolidation](release_consolidation_2026-09-07.md): branch inclusion,
   workflow ports, provenance, release plan.
4. [Selective Jacobian contract](selective_jacobian_plans.md) and
   [Jacobian investigation](jacobian_batched/overall_progress.md): architecture,
   coordinate maps, performance/precision evidence.

The original checkout `/home/cfranken/code/gitHub/vSmartMOM.jl` on
`feat/surface-split` has an unrelated unfinished `.gitignore` conflict plus
staged `.githooks` files. **Leave it intact; work in the release checkout.**
The separate performance worktree `/home/cfranken/code/gitHub/vSmartMOM-jacobian-perf`
was also preserved. Do not reset branches, rewrite shared/published history,
or interfere with Suniti's running workflows or private data/products.

## User intent and attribution

The user wants consolidated, future-extensible Jacobians and a properly
validated v2.2 release: numerical accuracy after convolution, retrieval effects,
documentation/docstrings, clean API definitions, schema/Taplo and CI.

The user explicitly clarified: **it is Suniti Sanghavi's code; Christian is
helping**. Scientific formulations, linearization, Raman and retrieval
methodology belong to Suniti. Christian's supporting work covers engineering,
performance, integration and release plumbing. Preserve this distinction.

- Scientific Git identity verified from history:
  `sunitisanghavi <suniti.sanghavi@gmail.com>`.
- Supporting local identity: `cfranken <fronge@gmail.com>`.
- **No co-author trailers.**
- Fourteen scientific author records on the unpublished candidate history were
  corrected in the previous session, preserving commit trees/messages and
  committer information. Other branches/remotes were not rewritten.
- Eight subsequent Suniti OCO workflow commits were ported with her authorship.
  JSON maps are in this directory; do not repeat those operations.
- Scientific docs corrections are in `7ae92577`, maintenance checks in `3c1add99`.

## Highest-priority next implementation

These are recommendations for the next work session, **not completed fixes**.

### A. Complete pressure-dependent line-absorption derivatives

`src/CoreRT/tools/lin_model_from_parameters.jl` around the comment
"The pressure tangent holds cross sections" constructs `τ̇_abs_psurf` from the
bottom molecular-column change while holding ordinary line cross sections fixed.
A full forward rebuild changes pressure-dependent line shapes/LUT interpolation.
The missing upstream response matters for a full surface-pressure retrieval;
it does not invalidate the analytic RT operator chain rule itself.

Start with a centered full-forward pressure finite-difference reproduction for
small direct-line and ABSCO cases. Establish the pressure/profile coordinate
contract, then add the missing cross-section/LUT pressure response at the
upstream optical-property boundary. Keep optical and detector/convolution tests
separate. CIA/MT_CKD have separate pressure-scaling terms and incomplete
abundance tangents. No new ppm bias for this omission has been measured yet.

### B. Reject source requests that currently lose physics silently

The tracked reproduction is
[documentation_audit_2026-09-08/source_contract_probe.jl](documentation_audit_2026-09-08/source_contract_probe.jl).
It demonstrated on CPU:

- `SolarBeam(sza=75)` ignores that source angle and uses model geometry.
- `SolarBeam()+SolarBeam()` yields exactly one beam: **0.5 times** the radiance
  of a single beam with correctly summed irradiance.
- Nonzero prescribed `SurfaceSIF` over an RPV surface is silently ignored.

Current source docs explain these limits. Runtime guards are **not** added yet.
Implement shared validation across forward, linearized, split and TOA paths;
test unsupported combinations explicitly. Multiple source geometries are a
larger feature; same-geometry irradiance can be summed into one beam today.

### C. Other valuable follow-ups

- Formalize instrument/state metadata: physical parameter identity, units,
  native/retrieval coordinates, band mapping, convolution/resampling/shift/width
  derivatives. Preserve `ParameterKey`/`JacobianPlan` instead of hard-coded offsets.
- Add durable retrieval-level precision acceptance after convolution and fitted
  grid shifts, not only raw radiance/per-sample-noise thresholds.
- Metal: `CoreRT._to_aa_arch` has CPU/CUDA mappings only; production absorption
  paths through it lack Metal support. Real Apple hardware is needed to validate
  end-to-end support; CPU macOS CI is not Metal validation.
- Replace three clearly labeled synthetic documentation figures (temperature
  absorption, surface angular comparison, canopy red edge) with actual tutorial
  output if useful. They are now labeled; no new physical calculation was added.
- Refresh stale numeric code anchors using symbol-aware mapping. A valid file
  link does not establish that an old line still points to the intended equation.
- The GCHP sectional-AOD branch has valuable work but an incomplete sectional-to-RT
  bridge and obsolete RT changes. Port a bounded adapter; do not blindly merge it.
- Thermal/Raman Jacobians remain unsupported; AD traits alone do not implement them.

## What the latest documentation pass completed

Manual-wide text/path scan: 63 Markdown pages, 12,017 lines, 202 ordinary Julia
fences. Targeted semantic review and execution found/fixed retired Mie AD calls,
wrong constructor/import/field examples, source/BRDF normalization descriptions,
radiance labels, source/backend support claims and stale legacy SIF behavior.

The Jacobian guide now explains instrument convolution and its shift/width chain
rule, noise covariance whitening, matched/native grids and the local retrieval
state-displacement estimate. The pressure tangent limitation is prominent.

**Coverage is not universal scientific certification:** not every equation,
historical note or all 999 docstrings have been independently reviewed. Ordinary
fences/comments can be illustrative or depend on external datasets. Documentation
CI now executes all ten Literate scripts and the small source/API contracts.

Latest checks (Julia 1.12.6, CPU, four Julia threads, one BLAS thread):

- All ten tutorials in one final run: **26/26**, exit 0, 3m50.8s.
- Documentation contracts: **20/20**, exit 0.
- Strict Documenter/VitePress build plus generated assets: exit 0; no deployment.
- API inventory: **475 exports, 999 docstrings, zero structural failures**.
- Repository blob-link scan: zero missing paths; heading/line semantics not certified.
- Workflow YAML and `git diff --check`: pass. Hosted CI has not run these commits.
- CUDA tutorial sections and optional GIF generation were excluded.

Logs: `/home/cfranken/vsmartmom-release-evidence/docs-audit-2026-09-08/`.
Failed attempts have separate filenames and were retained honestly.
Tracked summary: [inventory.json](documentation_audit_2026-09-08/inventory.json).

## Earlier integration evidence — preserve exact scope

The preceding full release work included package-owned batched multiplication
(no NNlib piracy), Cox–Munk absolute/source/tangent fixes, GPU inverse dispatch,
explicit rejection of placeholder experimental aerosol optics, source-only OCO
workflow tests, parser/schema/Taplo, docs dependencies and release CI updates.

- CUDA A100 full suite: **20,301/20,301** assertions, exit 0.
- OCO workflow harness with separately supplied read-only ancillary data:
  **1,139/1,139**; portable mode 120 passes/14 skips. No new production inversion.
- Full CPU run: **12,429 passes, one failure, one error, 14 skips**. Both issues
  were test-harness defects; fixes passed all **96 affected-group assertions**.
  **There was no complete green CPU rerun after those two fixes.** Do not claim one.
- Fresh Julia 1.10/1.12 resolution and smoke checks passed.
- JET reported **285 advisory findings** versus historical 232; not all reviewed.
- npm audit was zero findings after pinned dependency updates; Taplo/schema checks
  passed. See the release report rather than rerunning unrelated checks blindly.
- Registration/version-sequence handling, hosted CI, deployment credentials and
  publication remain outstanding. No credentials/publication should be assumed.

## Earlier precision and speed findings

Keep the baseline definitions and acceptance caveats attached to all numbers:

- Archived Suniti inversions versus isolated optimized replay: about **12× clear**,
  **19–20× aerosol**. This is historical evidence, not a fresh matched-worker A/B.
- Controlled comparison within the already optimized implementation: **3.25–3.33×**.
  Do not multiply unrelated stage speedups.
- Full Float32/Float64 retrieval displacement: **−0.0184 ppm on native grids**,
  **+0.00249 ppm with matched grid locations**; posterior sigma about **0.314 ppm**.
  Iteration decisions matched for that fixture.
- Original strict 0.01-noise-sigma diagnostic failed 3/6 cases, worst 0.02076 sigma.
  The user accepted the small physical scale; the original failed diagnostic was
  not relabeled as passing. Results and convolution/grid-shift scripts live under
  `docs/dev_notes/jacobian_batched/`.
- These are fixture-specific results, not guarantees for all retrievals.

## Environment and commands

The default `julia` launcher is broken on this machine. Use:

```text
/home/cfranken/.julia/juliaup/julia-1.12.6+0.x64.linux.gnu/bin/julia
/home/cfranken/.julia/juliaup/julia-1.10.12+0.x64.linux.gnu/bin/julia
```

Tests must run from `test/` because fixtures use relative paths. Typical commands:

```bash
cd /home/cfranken/code/gitHub/vSmartMOM-release/test
JULIA_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1 CUDA_VISIBLE_DEVICES='' \
  /home/cfranken/.julia/juliaup/julia-1.12.6+0.x64.linux.gnu/bin/julia \
  --project=../docs ../tools/check_tutorials.jl

/home/cfranken/.julia/juliaup/julia-1.12.6+0.x64.linux.gnu/bin/julia \
  --project=../docs ../docs/test_examples.jl
```

For docs, run from `docs/` with `CI=false` to disable local deployment.
Use `CUDA_VISIBLE_DEVICES=0` for the available A100 when GPU checks are needed.
`/tmp` was nearly full; prefer
`TMPDIR=/home/cfranken/vsmartmom-release-evidence/tmp`.
Use `python3.12` rather than the system Python 3.6. YAML/jsonschema environment:
`/home/cfranken/vsmartmom-release-evidence/checks-venv/bin/python`.
Correct GitHub CLI: `/home/cfranken/.local/bin/gh-cli` (the `gh` command is a
different Python program). Current permissions do not permit sandbox escalation
arguments; execute normal commands. Do not spawn agents unless explicitly asked.
