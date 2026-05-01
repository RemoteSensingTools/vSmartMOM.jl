# Draft: `index.md` after Vitepress cutover

This is the exact frontmatter block + body that replaces the current
`docs/src/index.md` on Vitepress cutover day. Saved here so the cutover
commit is mechanical.

The hero block uses Vitepress's native `layout: home` schema. See
[design_brief.md](../design_brief.md) for icon / palette specs.

Below the `---` fence, the body content is reduced — the hero already
covers tagline, install, and "where to go." Keep the body for the things
the hero cannot: public modules summary, citation pointer, contributing
note. Anything longer belongs on its own page.

---

```markdown
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
    details: Adding-doubling matrix operator method, polarized solver, analytic Jacobians, Raman / Cabannes inelastic path, GPU-ready.
    link: /pages/vSmartMOM/CoreRTTheory
---

## Install

vSmartMOM supports Julia 1.10 or later.

```julia
pkg> add vSmartMOM
```

## Public Modules

- **CoreRT** — adding-doubling solver, model types, optical-property assembly, surface coupling, Jacobian kernels.
- **IO** — YAML, TOML, Dict, NetCDF, and GEOS-Chem inputs.
- **Absorption** — HITRAN line-by-line and lookup-table gas absorption.
- **Scattering** — Mie calculations, phase functions, Greek coefficients, truncation inputs.
- **InelasticScattering** — Raman / Cabannes mode types and optical-property helpers.
- **Aerosols** — TOMAS-15 and two-moment aerosol input support. *API still being stabilized.*
- **SolarModel** — solar / stellar spectra and transmission helpers.

## Cite

If you use vSmartMOM.jl in published work, please cite the JOSS software paper and the underlying methodology papers — see [References](pages/vSmartMOM/References.md).
```

---

## Notes for the cutover commit

- Drop the existing `# vSmartMOM.jl` H1 — Vitepress hero replaces it.
- Drop the existing "For:" / "Next:" header — those belong on internal pages, not the home layout.
- Keep the existing index body's "Smallest Entry Point" code block? **No** — the hero's "Quick Start" action covers it. The home layout should be marketing-grade, not tutorial-grade. Tutorial belongs on `pages/quickstart.md`.
- The current `index.md` "Where To Go" bullet list is replaced by the three feature cards plus the three hero action buttons.
- Net change: index.md goes from ~45 lines of regular markdown to ~50 lines of frontmatter + ~20 lines of body.

## Open question for cutover day

- The link paths in `actions:` and `features:` use `/pages/...` style. Vitepress sometimes wants `/pages/quickstart` (no `.md`, no leading dot, no trailing slash). DocumenterVitepress documents the exact convention; verify against a Lux.jl or Oceananigans landing page before pasting.
