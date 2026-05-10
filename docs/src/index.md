```@raw html
---
layout: home

hero:
  name: "vSmartMOM.jl"
  text: "Polarized atmospheric radiative transfer"
  tagline: "<strong style=\"color:#dc2626\">v</strong>ector <strong style=\"color:#7c3aed\">S</strong>imulated <strong style=\"color:#dc2626\">m</strong>easurements of the <strong style=\"color:#dc2626\">a</strong>tmosphere using <strong style=\"color:#dc2626\">r</strong>adiative <strong style=\"color:#dc2626\">t</strong>ransfer based on the <strong style=\"color:#7c3aed\">M</strong>atrix <strong style=\"color:#7c3aed\">O</strong>perator <strong style=\"color:#7c3aed\">M</strong>ethod — analytic Jacobians, Raman scattering, CUDA-ready."
  image:
    src: /assets/icons/logo.png
    alt: vSmartMOM logo
  actions:
    - theme: brand
      text: Quick Start
      link: /pages/quickstart
    - theme: alt
      text: RT basics (the story)
      link: /pages/concepts/01_overview
    - theme: alt
      text: View on GitHub
      link: https://github.com/RemoteSensingTools/vSmartMOM.jl

features:
  - icon:
      src: /assets/icons/scattering.png
    title: Scattering
    details: Mie theory, Greek-coefficient phase matrices, NAI2 / PCW Fourier decomposition, vector δ-m truncation.
    link: /pages/concepts/03b_scattering
  - icon:
      src: /assets/icons/absorption.png
    title: Absorption
    details: HITRAN line-by-line cross sections with Voigt, Doppler, and Lorentz line shapes; lookup-table interpolation for hot loops.
    link: /pages/concepts/03a_absorption
  - icon:
      src: /assets/icons/radiative_transfer.png
    title: Radiative Transfer
    details: Adding-doubling matrix operator method, polarized solver, analytic Jacobians, Raman / Cabannes inelastic path, GPU-ready.
    link: /pages/concepts/04_mom_solver
---
```

## Install

vSmartMOM supports Julia 1.10 or later.

```julia
pkg> add vSmartMOM
```

## Public Modules

- **[CoreRT](pages/api/core_rt.md)** — adding-doubling solver, model types, optical-property assembly, surface coupling, Jacobian kernels.
- **[IO](pages/api/io.md)** — YAML, TOML, Dict, NetCDF, and GEOS-Chem inputs.
- **[Absorption](pages/api/absorption.md)** — HITRAN line-by-line and lookup-table gas absorption.
- **[Scattering](pages/api/scattering.md)** — Mie calculations, phase functions, Greek coefficients, truncation inputs.
- **[InelasticScattering](pages/api/inelastic.md)** — Raman / Cabannes mode types and optical-property helpers.
- **[Aerosols](pages/api/aerosols.md)** — TOMAS-15 and two-moment aerosol input support. *API still being stabilized.*
- **[SolarModel](pages/api/solar_model.md)** — solar / stellar spectra and transmission helpers.

## Cite

If you use vSmartMOM.jl in published work, please cite the JOSS software paper and the underlying methodology papers — see [References](pages/vSmartMOM/References.md).
