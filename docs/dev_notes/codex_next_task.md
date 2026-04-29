# Codex Next Task — Retire `warnonly` Categories

Created: 2026-04-28
Owner: Codex (GPT-5.5)
Companion docs (DO NOT START THESE YET):
- [docs_vitepress_migration.md](docs_vitepress_migration.md) — deferred until this task is done
- [design_brief.md](design_brief.md) — deferred until Vitepress cutover

## Why

The current `docs/make.jl` carries:

```julia
warnonly = [:docs_block, :missing_docs, :cross_references]
checkdocs = :exports
```

That suppression is hiding three different classes of issue:

- `:docs_block` — duplicate `@docs` entries (already cleaned by 8c8273a, may be empty now)
- `:missing_docs` — exported symbols without docstrings (since `checkdocs = :exports` is on)
- `:cross_references` — broken `[link]` references (most dangerous: docs that lie about navigation)

Until these are retired, "the docs build cleanly" is not a meaningful claim. This task makes it meaningful.

## Out of scope (do not touch)

- DocumenterVitepress migration — separate plan, deferred.
- Logo / icon / hero block — separate plan, deferred.
- Legacy `pages/vSmartMOM/{Overview,Example,InputParametersGuide,Principles}.md` rewrite/delete — separate triage task, deferred.
- SIF data policy — product decision, not docs.
- `@deprecate parameters_from_yaml` — API decision, not docs.
- Tutorial heavy-execution wiring — already partially handled.

If during the work you find a broken cross-reference that points at one of the legacy `pages/vSmartMOM/*.md` pages, fix the link side, not the legacy page.

## Task

Land **three small commits** in this order. Each commit must be reviewable in isolation and the docs build must pass after each.

### Commit A — Build inventory

Files: none (informational only).

Action: run

```bash
julia --project=docs --startup-file=no docs/make.jl 2>&1 | tee /tmp/docs_build.log
```

then extract:

```bash
grep -E "Warning|missing docstring|cross-reference|@docs" /tmp/docs_build.log > docs/dev_notes/warning_inventory.md
```

Add `docs/dev_notes/warning_inventory.md` to the repo so the cleanup is bisectable. Group warnings under three headings (`## Docs Block`, `## Missing Docs`, `## Cross References`) and number each entry. This file is the work list for commits B and C.

Validation:
- `docs/dev_notes/warning_inventory.md` exists and has all three sections.
- The file shows zero or N entries per section, with N being concrete.

### Commit B — Retire `:docs_block` and `:cross_references`

Files: `docs/make.jl` plus whatever `.md` files surface in commit A's inventory under those two headings.

Actions:

1. **`:docs_block`**: 8c8273a already removed the duplicate blocks Claude flagged. If the inventory shows zero remaining, just drop `:docs_block` from the `warnonly` list. If any survive, dedupe them per the 8c8273a pattern (canonical block in `api_reference.md`, prose link from elsewhere) before dropping.
2. **`:cross_references`**: For each broken link, choose one of:
   - Fix the link target if the destination exists under a different path.
   - Remove the link and replace with prose if the destination was retired.
   - Add the missing anchor on the destination page if the link is correct but the anchor is gone.
   Do **not** redirect via shim files; fix the link side.
3. Drop `:cross_references` from `warnonly`.

After this commit, `warnonly = [:missing_docs]` is the only remaining suppression.

Validation:
- `julia --project=docs --startup-file=no docs/make.jl` exits 0 with no warnings in the `:docs_block` or `:cross_references` classes.
- Update `docs/dev_notes/warning_inventory.md` to mark resolved entries.

### Commit C — Retire `:missing_docs`

Files: `docs/make.jl`, plus selected `.md` files in `docs/src/pages/api_reference.md`, plus possibly `src/*/*.jl` for de-export decisions.

For each missing-docstring entry from commit A, choose:

- **(a) Add a docstring** if the symbol is genuinely user-facing. Place the docstring on the source definition. Link it into the API reference if not already linked.
- **(b) De-export** the symbol if it is internal-only. Remove from the `export` line in the relevant module. The audit specifically called out the Inelastic helper exports (`needs_interaction_workspace`, `normalize_raman_weights!`, traits) as de-export candidates — confirm with the source author before removing if unclear, otherwise default to (a) with a one-line "internal helper" docstring.

Drop `:missing_docs` from `warnonly`.

After this commit, `make.jl` should have **no** `warnonly` argument (or `warnonly = false` if Documenter requires explicit).

Validation:
- `julia --project=docs --startup-file=no docs/make.jl` exits 0 with no warnings.
- `docs/dev_notes/warning_inventory.md` shows all entries resolved (or explicitly justified if any deferred).

## Gate before next phase

Once commit C is in and the build is fully clean, the next planned work is the **Vitepress cutover** (see [docs_vitepress_migration.md](docs_vitepress_migration.md)). That work touches `docs/make.jl` heavily and would conflict if interleaved with this cleanup.

Do not start the Vitepress migration until this task signs off.

## Hand-off

When done:
- Push commits A, B, C.
- Update `docs/dev_notes/warning_inventory.md` with the final state.
- Reply with the three commit SHAs and a one-line summary of what was de-exported (if anything in commit C).

Then the human reviews the de-export choices and unblocks the Vitepress phase.
