# Codex Next Task — Hero Block + Icon Set

Created: 2026-04-29
Branch: `cfranken/refactor-low-risk-cleanup`
Latest commit at handoff: `d4fd752`
Most recent CI run: #674 (triggered by commit `6d32d35` adding NCDatasets +
Interpolations to test/Project.toml — verify it's green before merging
anything new).

## Why

The Vitepress renderer cutover landed in commit `164c98d` but ships with the
default landing page. A proper hero block + three feature cards + master mark
are the next visual upgrade. All inputs are already on disk; this is a
mechanical execution.

## Inputs (all in `docs/dev_notes/`)

| File | What it gives you |
|---|---|
| [design_brief.md](design_brief.md) | Palette, stroke weights, geometry rules, sizing targets. Read first. |
| [design/icons/logo.svg](design/icons/logo.svg) | Master mark — 3-layer atmospheric stack with photon dot. |
| [design/icons/scattering.svg](design/icons/scattering.svg) | Phase-function lobe with forward-asymmetric ray fan. |
| [design/icons/absorption.svg](design/icons/absorption.svg) | Spectrum continuum with two ink notches. |
| [design/icons/radiative_transfer.svg](design/icons/radiative_transfer.svg) | Layered atmosphere with bidirectional flux arrows. |
| [design/index_hero_draft.md](design/index_hero_draft.md) | Exact Vitepress frontmatter to paste into `docs/src/index.md`. |
| [design/preview.html](design/preview.html) | Open in a browser to sanity-check icon proportions and theming. |

Each SVG uses embedded `prefers-color-scheme` so it self-themes as `<img>` —
no need for separate light/dark variants.

## Task

Single commit. Mechanical execution.

### Steps

1. **Move SVGs into the served assets directory.**
   ```
   mkdir -p docs/src/assets/icons
   git mv docs/dev_notes/design/icons/logo.svg docs/src/assets/icons/logo.svg
   git mv docs/dev_notes/design/icons/scattering.svg docs/src/assets/icons/scattering.svg
   git mv docs/dev_notes/design/icons/absorption.svg docs/src/assets/icons/absorption.svg
   git mv docs/dev_notes/design/icons/radiative_transfer.svg docs/src/assets/icons/radiative_transfer.svg
   ```
   Use `git mv` so history follows.

2. **Update `docs/make.jl` `assets` field** so Documenter copies the new icons:
   The current `MarkdownVitepress(...)` call has
   `assets = ["assets/favicon.ico", "assets/logo.png"]`. Either:
   - Drop `assets/logo.png` if no longer wanted (replaced by SVG logo), or
   - Keep both (PNG for social cards, SVG for the hero/favicon)
   - Add `"assets/icons/scattering.svg"`, `"assets/icons/absorption.svg"`,
     `"assets/icons/radiative_transfer.svg"`, `"assets/icons/logo.svg"` if
     Documenter doesn't already auto-include the assets folder.

3. **Replace `docs/src/index.md`** with the content in
   [`docs/dev_notes/design/index_hero_draft.md`](design/index_hero_draft.md).
   The frontmatter block is the hero; the body underneath is the public
   modules summary + citation pointer. Strip the wrapper triple-fence
   (\`\`\`markdown ... \`\`\`) — that's only there to render in the dev_notes
   preview.

4. **Verify the link paths** in `actions:` and `features:` in the hero
   frontmatter resolve to real pages. The draft uses `/pages/quickstart`,
   `/pages/IO/Overview`, `/pages/Scattering/Overview`, `/pages/Absorption/Overview`,
   `/pages/vSmartMOM/CoreRTTheory` — confirm against the `pages` array in
   `docs/make.jl`. DocumenterVitepress sometimes wants `.md` extensions and
   sometimes wants leading slashes — check the Lux.jl or Oceananigans landing
   pages for the exact convention if unsure.

5. **Build locally and eyeball it.**
   ```
   julia --project=docs --startup-file=no docs/make.jl
   cd docs && npm run docs:preview
   ```
   Open the preview URL. Verify the hero renders with the master mark, the
   three feature cards have icons, and clicking each action button goes
   somewhere real. Check dark mode toggle.

6. **Commit.** Suggested message:
   ```
   Add hero block and icon set to landing page

   Promotes the SVG icons from dev_notes/design/icons/ into
   docs/src/assets/icons/ and replaces docs/src/index.md with the
   Vitepress hero frontmatter from design/index_hero_draft.md. Per
   design_brief.md.

   - logo.svg: master mark, 3-layer stack with photon dot
   - scattering.svg: forward-asymmetric phase-function lobe
   - absorption.svg: spectrum continuum with two ink notches
   - radiative_transfer.svg: layered atmosphere with up/down arrows

   All four use embedded prefers-color-scheme so they self-theme as
   <img> tags. Hero block has three actions (Quick Start, Core RT
   Theory, GitHub) and three feature cards (Scattering, Absorption,
   Radiative Transfer).
   ```

## Acceptance criteria

- `julia --project=docs --startup-file=no docs/make.jl` exits 0 with no new warnings.
- The four SVGs render correctly in the built docs at the feature-card sizes.
- The master mark renders at the hero size and as the favicon.
- All three action buttons link to real pages.
- All three feature card links go to real module overview pages.
- Dark mode (Vitepress's default toggle) shows the SVGs in their dark variant
  via `prefers-color-scheme: dark` (the system mode the browser reports —
  Vitepress's manual toggle doesn't trigger media-query switches; if that's an
  issue, follow-up commit can add `class="dark"` overrides).
- CI run #675 (or whichever number this commit triggers) is green on Linux.

## Out of scope for this commit

- DocumenterVitepress sidebar customization (`docs/src/.vitepress/config.mts`).
  The default sidebar from `pages = Any[...]` is fine; iterate later.
- Algolia DocSearch upgrade.
- Custom 404 page.
- Wordmark variant of the logo (mentioned as optional in design_brief.md).
- Social card PNG export (1200×630 for GitHub previews).
- Moving `docs/dev_notes/design/CITATION.bib` to repo root — separate
  release-prep step.

## What's still red / open after this commit

- **SIF data policy** — disabled in `test/runtests.jl` line 83. No CI cost
  but `test_sif.jl` stays disabled until fixtures land or loaders gain a
  fallback path. Product decision, not engineering.
- **Internal API cleanup** — see [internal_api_cleanup.md](internal_api_cleanup.md).
  48 currently-public internals that should be triaged before v2.0.1 / v2.1.0.
  Not blocking registration.
- **DocumenterVitepress 0.3 sidebar tweaks** — the auto-generated nav from
  `pages = Any[...]` may want grouping. Iterate when convenient.

## Pipeline state after this commit

Refer to [handoff_runbook.md](handoff_runbook.md) and update it: move
[design_brief.md](design_brief.md) and [design/index_hero_draft.md](design/index_hero_draft.md)
from Deferred to Consumed with the hero-commit SHA. The icon SVGs and
`preview.html` stay in the dev_notes audit trail (they're production-source,
not consumable specs).

## Hand-off back

When done, push and reply with:
- The hero commit SHA
- A screenshot of the rendered hero in light + dark mode (optional, but nice)
- Whether CI on the new commit is green
- Any deviations from this brief
