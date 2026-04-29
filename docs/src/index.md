```@raw html
---
layout: home

hero:
  name: "vSmartMOM.jl"
  text: "Polarized atmospheric radiative transfer"
  tagline: "Vector matrix-operator method with analytic Jacobians, Raman scattering, and CUDA support."
  image:
    src: /assets/icons/logo.svg
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
```

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
