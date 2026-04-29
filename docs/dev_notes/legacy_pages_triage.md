# Legacy `pages/vSmartMOM/*` Triage

Created: 2026-04-28
Status: **next Codex task** — start after [codex_next_task.md](codex_next_task.md) (warnonly retirement) is signed off.

The current nav still includes four legacy pages under the `vSmartMOM`
section that pre-date the task-based / module-reference split. Each file
was inspected; per-file decisions and rationale below.

The audit and the story scaffold both flagged this as deferred work. This
doc converts that deferral into a concrete commit plan.

## Per-file decisions

| File | Action | Rationale |
|---|---|---|
| [pages/vSmartMOM/Overview.md](../src/pages/vSmartMOM/Overview.md) | **KEEP, trim** | Already updated for `read_parameters` and `RTModel` hierarchy. Acts as the module-landing for the reference layer. Trim to ~15 lines. |
| [pages/vSmartMOM/Example.md](../src/pages/vSmartMOM/Example.md) | **DELETE** | Already duplicated by the new [`pages/quickstart.md`](../src/pages/quickstart.md) (forward) and [`pages/jacobians.md`](../src/pages/jacobians.md) (linearized) task pages. Two examples in two places diverge over time. |
| [pages/vSmartMOM/InputParametersGuide.md](../src/pages/vSmartMOM/InputParametersGuide.md) | **FOLD into [IO/Schema.md](../src/pages/IO/Schema.md), then DELETE** | Same schema documented in both files. `IO/Schema.md` is more authoritative (covers safe parsing, `eval` scope, all four loaders). |
| [pages/vSmartMOM/Principles.md](../src/pages/vSmartMOM/Principles.md) | **FOLD into [CoreRTTheory.md](../src/pages/vSmartMOM/CoreRTTheory.md), then DELETE** | Has substantive truncation equations from Sanghavi & Stephens (2015) that are *missing* from the current CoreRTTheory page. CoreRTTheory should be the single canonical theory page. |
| [pages/vSmartMOM/Types.md](../src/pages/vSmartMOM/Types.md) | **KEEP** (no change) | Already cleaned up in 8c8273a. Stays as the module's type/method reference. |
| [pages/vSmartMOM/References.md](../src/pages/vSmartMOM/References.md) | **KEEP, refresh** | Should be the bibliographic ground truth for the package. Verify against [theory_references.md](theory_references.md) — should cite S2014, S2015, S2022-I, SF2023-II, JSF2022 at minimum. |

## Commit plan

Three small commits, each independently reviewable.

### Commit T1 — Trim Overview.md

File: `docs/src/pages/vSmartMOM/Overview.md`

Reduce to module landing. Target structure:

```markdown
# vSmartMOM Module Overview

The vSmartMOM module wires parameters, model construction, and the solver into
end-to-end forward and linearized RT runs.

For getting-started workflows see [Quick Start](../quickstart.md).
For Jacobians see [Compute Jacobians](../jacobians.md).
For equation-level solver mapping see [Core RT Theory](CoreRTTheory.md).

## Architecture

![Architecture diagram](vSmartMOMDiagram-vSmartMOM.drawio.png)
```

Drop the existing 3-step "perform an RT simulation" block (covered by
Quick Start) and the bullet list of features (covered by index hero
once Vitepress lands; for now keep one short paragraph).

### Commit T2 — Fold Principles.md into CoreRTTheory.md, delete Principles.md

Files:
- `docs/src/pages/vSmartMOM/CoreRTTheory.md` (edited)
- `docs/src/pages/vSmartMOM/Principles.md` (deleted)
- `docs/make.jl` (remove Principles entry from nav)

Fold targets within CoreRTTheory.md:

| Source section in Principles.md | Target section in CoreRTTheory.md |
|---|---|
| "Truncation strategies and accuracy" + the three equation blocks (truncation factor, modified Greek coefficients, modified optical properties) | The existing Truncation section. Currently very thin (2 mentions of "truncation"); these equations make it concrete. Cite Sanghavi & Stephens 2015 Eqs. 26, 27a-27d. |
| "Polarization handling" | New short subsection under the Vector RTE intro |
| "Linearization (aerosol Jacobians)" + the two Sanghavi quotes | The existing Linearization section. Use the quotes for color, then point at S2014 App. C for the actual derivation |
| "Layer optics and composition" | Already covered by Elemental + Doubling sections; only fold what isn't already there |
| "Surface interaction and HDRF" | New subsection at the end of CoreRTTheory before "Other dispatch arms" |
| "Practical setup tips" | Move to [Run on GPU](../gpu.md) for the architecture/perf bits and [IO Schema](../IO/Schema.md) for the accuracy-knobs bits |

Anything that does not have a target in CoreRTTheory.md or another existing
page does not get folded; it gets dropped. Don't preserve content for
preservation's sake.

After this commit, CoreRTTheory.md should mention `truncation` more than
twice (currently 2; aim for ~10 with one full equation block).

### Commit T3 — Delete Example.md, fold InputParametersGuide.md into Schema, delete

Files:
- `docs/src/pages/vSmartMOM/Example.md` (deleted)
- `docs/src/pages/vSmartMOM/InputParametersGuide.md` (deleted)
- `docs/src/pages/IO/Schema.md` (edited if any field is documented in InputParametersGuide.md but not in Schema.md — verify before deleting)
- `docs/make.jl` (remove both entries from nav)
- `docs/src/pages/vSmartMOM/References.md` (refresh against [theory_references.md](theory_references.md))

Fold check for InputParametersGuide.md → Schema.md:

| Field documented in InputParametersGuide.md | In Schema.md? |
|---|---|
| spec_bands | yes ✓ |
| surface (all 6 BRDF types) | yes for 5; verify CoxMunk note about whitecaps/shadowing kwargs |
| quadrature_type / polarization_type / max_m / Δ_angle / l_trunc / depol / float_type / architecture | yes ✓ |
| sza / vza / vaz / obs_alt | yes ✓ |
| T / p / q / profile_reduction | yes ✓ |
| absorption: molecules / vmr / broadening_function / CEF / wing_cutoff | spot-check |
| scattering: aerosols (τ_ref, μ, σ, nᵣ, nᵢ, p₀, σp), r_max, nquad_radius, λ_ref, decomp_type | spot-check |

For any field documented in InputParametersGuide.md but missing from
Schema.md, port the description over before deleting.

References.md refresh: verify bibliographic entries match
[theory_references.md](theory_references.md). Add JSF2022 (the JOSS
software paper) prominently — it is the primary citation for the package.

## Validation per commit

After each of T1, T2, T3:
- `julia --project=docs --startup-file=no docs/make.jl` exits 0.
- The docs build emits no broken-link warnings (`:cross_references` is
  retired by then).
- `docs/test_examples.jl` still passes if it touched any covered code.

## What stays in the `vSmartMOM` nav section after triage

```
vSmartMOM
  Overview            (trimmed, kept)
  Methods & Types     (kept as-is)
  References          (kept, refreshed)
```

`CoreRTTheory.md` stays at top-level nav (already promoted out of this
subsection). Everything else is folded or deleted.

## Out of scope for this triage

- The `Absorption/{Overview,Example,Types,References,HITRAN_Data}.md` and
  `Scattering/{Overview,Example,Types,References}.md` subsections — those
  follow the same pattern but were not flagged as stale by the audit.
  Defer until someone flags concrete drift.
- Module overview pages added in commit d8ba1a0 (`Inelastic/Overview.md`,
  `Aerosols/Overview.md`, `SolarModel/Overview.md`, `Surfaces/Overview.md`)
  — those are new and well-scoped.

## Hand-off

When done:
- Push T1, T2, T3.
- Update [docs_page_edit_plan.md](docs_page_edit_plan.md) to mark this
  task complete.
- Reply with the three commit SHAs and a summary of any content dropped
  rather than folded.

Then the human reviews and unblocks the [Vitepress
cutover](docs_vitepress_migration.md).
