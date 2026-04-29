# `docs/dev_notes/` — Handoff Runbook

Created: 2026-04-28

One-page index of every planning artifact under `docs/dev_notes/`. Read top
to bottom if you are picking this up cold; jump to the **Active** section
if you already know the context.

## Status legend

| Symbol | Meaning |
|---|---|
| ✅ | **Consumed** — the work this artifact planned has been done. Kept for the historical record and audit trail. Do not edit. |
| ▶ | **Active** — currently being executed or in immediate use. |
| ⏸ | **Deferred** — ready to consume but waiting on an upstream gate (Codex completion, product decision, registration cutover). |
| 📖 | **Reference** — long-lived lookup material that doesn't get "consumed." |

## Active (Codex is working from these)

| File | Purpose |
|---|---|
| (none — see Deferred for next pickup) | |

## Deferred (next pickups, in order)
| ⏸ [design_brief.md](design_brief.md) | **Next pickup.** Specifies the visual system for the hero block, icons, and palette. The follow-up commit moves the four staged SVGs into `docs/src/assets/icons/`, adds the hero frontmatter to `docs/src/index.md`, and verifies the build. |
| ⏸ [design/index_hero_draft.md](design/index_hero_draft.md) | Exact frontmatter to paste into `docs/src/index.md` on Vitepress hero day. |
| ⏸ [design/icons/](design/icons/) | Four production-ready SVGs (`logo`, `scattering`, `absorption`, `radiative_transfer`). Move to `docs/src/assets/icons/` on Vitepress hero day. |
| ⏸ [design/preview.html](design/preview.html) | Self-contained browser preview of all four icons at multiple sizes in both color schemes. Open anytime for visual sanity check. |
| ⏸ [design/CITATION.bib](design/CITATION.bib) | Six canonical citations. Two uses: (1) drop at repo root for GitHub's "Cite this repository" widget, (2) input for `pages/vSmartMOM/References.md` refresh during legacy_pages_triage T3. |
| ⏸ [internal_api_cleanup.md](internal_api_cleanup.md) | Post-v2.0.0 follow-up. Triages the 48 internal exports into four action buckets with a six-commit execution plan. Not blocking. |

## Reference (read anytime; no consumption)

| File | What it is |
|---|---|
| 📖 [theory_references.md](theory_references.md) | Paper ↔ code ↔ doc crib sheet. Maps Sanghavi 2014 / 2015 / 2022-I / Sanghavi-Frankenberg 2023-II equations to source files in `src/CoreRT/CoreKernel/`. **The single source of truth for citations and equation provenance.** |
| 📖 [warning_inventory.md](warning_inventory.md) | Records the 48 missing-doc warnings that the warnonly retirement resolved. Useful as the input set for `internal_api_cleanup.md`. |

## Consumed (historical / audit trail)

| File | What it planned | When consumed |
|---|---|---|
| ✅ [docs_change_plan_for_claude.md](docs_change_plan_for_claude.md) | Original docs-overhaul plan from the user. | Audited and superseded by the implementation plan below. |
| ✅ [docs_change_plan_audit.md](docs_change_plan_audit.md) | Claude's red-team audit of the original plan. | Decisions baked into all subsequent planning docs. |
| ✅ [docs_story_scaffold.md](docs_story_scaffold.md) | Narrative IA — personas, journey arc, per-page job-to-be-done. | Current docs `pages/` structure built from this. |
| ✅ [docs_page_edit_plan.md](docs_page_edit_plan.md) | Codex's 12-commit execution plan for the docs overhaul. | Commits 6494ab8 → 2575a9d (plus polish commits a554675, 8c8273a). |
| ✅ [codex_next_task.md](codex_next_task.md) | The warnonly retirement task: three commits A/B/C. | Commits 491e30a, eaa8fd9, 42d9366. |
| ✅ [legacy_pages_triage.md](legacy_pages_triage.md) | Per-file decisions for the four legacy `pages/vSmartMOM/*.md` files. | Commits f455de6, 38d88e0, 776eb97. |
| ✅ [docs_vitepress_migration.md](docs_vitepress_migration.md) | DocumenterVitepress renderer cutover. | Commit 164c98d. `docs/src/.vitepress/config.mts` not committed — DocumenterVitepress 0.3 generates at build time; revisit if sidebar grouping needs customization. |

## Recommended reading order for cold pickup

1. `handoff_runbook.md` (this file) — orient.
2. `theory_references.md` — understand what cites what.
3. `docs_story_scaffold.md` — understand why the docs are organized this way.
4. `legacy_pages_triage.md` (if active) or whichever `⏸` is next.
5. `internal_api_cleanup.md` — for context on the 48 internal-export situation.

Skim the consumed files only if you need to understand a past decision.

## Maintenance rule

When a deferred artifact is consumed:
1. Move its row from the **Deferred** section to **Consumed**.
2. Add the consuming commit SHA(s) in the "When consumed" column.
3. Do not delete the file — the audit trail matters.

When a new planning artifact is added:
1. Add it to the appropriate section (usually **Deferred** or **Reference**).
2. Add a one-line purpose statement.
3. If it has an upstream gate, note that gate in the "When to consume" column.
