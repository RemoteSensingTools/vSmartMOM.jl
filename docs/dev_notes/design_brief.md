# vSmartMOM.jl — Visual Design Brief

Created: 2026-04-28
Status: **deferred** — to be implemented after Vitepress cutover ([docs_vitepress_migration.md](docs_vitepress_migration.md)).

Specifies the visual system for the docs landing page, icon set, and master
mark. Written so an illustrator (or Claude/Codex) can produce the SVG files
without further design iteration.

## Design system

### Palette

| Role | Light | Dark |
|---|---|---|
| Ink (lines, geometry) | `#1a1a1a` | `#e6e6e6` |
| Accent (the photon) | `#4a90c2` Rayleigh blue | `#7ab8e0` lighter blue |
| Background | `#ffffff` | `#0f1115` |
| Muted (secondary text) | `#666` | `#999` |

The accent color is **always** the photon — the light being acted on. Ink is
always the medium (particle, molecule, atmosphere, surface). This convention
is what makes the three icons read as one family.

### Geometry

- Stroke weight: **1.5px** on a 24×24 grid; scales linearly (2px at 32px,
  4px at 64px). Use `stroke-width="1.5"` and a `viewBox="0 0 24 24"`.
- Stroke linecap: `round`. Stroke linejoin: `round`.
- Safe margin: 2px on all sides.
- No fill except for the scattering particle dot and the master-mark photon.
- Geometric primitives only: circles, straight lines, arcs. No bezier
  curves except for the absorption notch.

### Sizing

Each glyph is designed at 24×24 and tested at:

- 16px (favicon)
- 24px (inline)
- 64px (Vitepress feature card)
- 256px (social card / README header)

If it doesn't read at 16px, the design is wrong.

## Three glyphs

### 1. Scattering — phase-function lobe

```
Concept: central particle with a forward-asymmetric ray fan
SVG file: docs/src/assets/icons/scattering.svg
```

Construction:
- Filled ink circle at `(12, 12)`, radius `2` — the particle.
- 8 accent strokes radiating from the particle, lengths varying to suggest a
  Henyey-Greenstein phase function with `g ≈ 0.5`:
  - Forward (3 o'clock direction): length `9` (longest)
  - 1:30 / 4:30: length `7`
  - 12 / 6: length `5`
  - 10:30 / 7:30: length `4`
  - Backscatter (9 o'clock): length `3` (shortest)
- All strokes start at the particle's edge, not its center.

Reads as: "burst of light off a particle, with a preferred direction."

### 2. Absorption — line in a continuum

```
Concept: flat horizontal spectrum with two sharp downward notches
SVG file: docs/src/assets/icons/absorption.svg
```

Construction:
- Horizontal accent stroke from `(2, 12)` to `(22, 12)` — the continuum.
- Two ink "notches" (V-shapes) dipping below the line:
  - First notch: peak at `(8, 12)`, valley at `(9, 18)`, peak at `(10, 12)`.
  - Second notch: peak at `(15, 12)`, valley at `(16, 18)`, peak at `(17, 12)`.
- Notches are ink color, drawn after (on top of) the accent line so they
  visually break it.

Reads as: "missing light at specific frequencies."

### 3. Radiative Transfer — layered solver

```
Concept: layered atmosphere with bidirectional flux arrows
SVG file: docs/src/assets/icons/radiative_transfer.svg
```

Construction:
- Three thin ink horizontal lines:
  - `(4, 6)` to `(20, 6)`
  - `(4, 11)` to `(20, 11)`
  - `(4, 16)` to `(20, 16)`
- One thicker (2.5px) ink line at the bottom — the surface:
  - `(3, 21)` to `(21, 21)`
- Two accent arrows, side-by-side near the right side, going through all
  layers:
  - Down-arrow: `(15, 4)` → `(15, 19)` with arrowhead at the bottom
  - Up-arrow: `(18, 19)` → `(18, 4)` with arrowhead at the top

Reads as: "light propagates up and down through atmospheric layers."

## Master mark

Used for: favicon, README header, Vitepress hero, social card.

```
Concept: 3-layer RT stack with a single accent dot in the middle layer
SVG file: docs/src/assets/logo.svg
```

Construction:
- Same 3 ink horizontal lines + thicker bottom surface line as the RT glyph.
- One filled accent circle at `(12, 11)`, radius `1.5` — the photon being
  transported. Positioned exactly on the middle layer line.
- No arrows. The single dot is what distinguishes it from the RT glyph and
  what makes it instantly recognizable at favicon size.

Optional alternative (more "branded" but heavier): a wordmark version with
"vSmartMOM" set in a monospaced sans-serif (JetBrains Mono or IBM Plex Mono)
to the right of the glyph, ink color, sentence case. Use the wordmark in the
hero block only; the icon-only version goes everywhere else.

## Hero block (Vitepress frontmatter)

Lives at the top of `docs/src/index.md` once the Vitepress cutover is done.

```yaml
---
layout: home

hero:
  name: "vSmartMOM.jl"
  text: "Polarized atmospheric radiative transfer"
  tagline: "Vector matrix-operator method with analytic Jacobians, Raman scattering, and CUDA support."
  image:
    src: /assets/logo.svg
    alt: vSmartMOM logo
  actions:
    - theme: brand
      text: Quick Start
      link: /pages/quickstart
    - theme: alt
      text: Core RT Theory
      link: /pages/vSmartMOM/CoreRTTheory
    - theme: alt
      text: View on GitHub
      link: https://github.com/RemoteSensingTools/vSmartMOM.jl

features:
  - icon:
      src: /assets/icons/scattering.svg
    title: Scattering
    details: Mie theory, Greek-coefficient phase matrices, NAI2 / PCW Fourier decomposition, vector δ-m truncation.
    link: /pages/Scattering/Overview
  - icon:
      src: /assets/icons/absorption.svg
    title: Absorption
    details: HITRAN line-by-line cross sections with Voigt, Doppler, and Lorentz line shapes; lookup-table interpolation for hot loops.
    link: /pages/Absorption/Overview
  - icon:
      src: /assets/icons/radiative_transfer.svg
    title: Radiative Transfer
    details: Adding-doubling matrix operator method, polarized solver, analytic Jacobians, Raman/Cabannes inelastic path, GPU-ready.
    link: /pages/vSmartMOM/CoreRTTheory
---
```

The current `index.md` body stays under the frontmatter — it becomes the
"below the fold" content for users who scroll past the hero.

## File deliverables

Produced in this order, all under `docs/src/assets/`:

1. `logo.svg` — master mark, ~50 lines
2. `icons/scattering.svg` — ~30 lines
3. `icons/absorption.svg` — ~25 lines
4. `icons/radiative_transfer.svg` — ~35 lines
5. `favicon.svg` — same as logo.svg, browser auto-handles sizing
6. `social-card.png` — 1200×630 PNG generated from the master mark + wordmark, for GitHub/Twitter previews

Total: ~150 lines of hand-written SVG plus one PNG export.

## Pre-cutover staging

While the Vitepress cutover is pending, draft SVGs live under
`docs/dev_notes/design/icons/` so they are not picked up by the current
Documenter build. On cutover day, move them to `docs/src/assets/`.

## Acceptance

- All four SVGs render correctly at 16px, 24px, 64px, 256px.
- Light and dark mode both look correct (test in Vitepress preview).
- The three feature icons read as a family — same stroke weight, same
  accent rule, same compositional footprint.
- The master mark is recognizable at favicon size in a browser tab.
- No hand-written SVG file exceeds 60 lines (forces minimalism).

## Inspiration / not-inspiration

- **Look like:** Oceananigans landing, Lux.jl landing, Documenter itself.
  Restrained, geometric, scientific.
- **Don't look like:** generic SaaS landing page, Material Design icon
  pack, weather app glyphs (no smiling suns), particle-physics stock art
  (no glowing atom-orbital nonsense).
