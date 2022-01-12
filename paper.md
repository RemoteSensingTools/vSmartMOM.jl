---
title: 'RadiativeTransfer.jl: an Open-Source Julia Package for Atmospheric Remote Sensing Tools'

tags:
  - Julia
  - radiative transfer
  - atmospheric radiation
authors:
  - name: Rupesh Jeyaram
    orcid: 0000-0003-0142-7367
    affiliation: 1
  - name: Suniti Sanghavi 
    orcid: 0000-0003-0754-9154
    affiliation: "1, 2"
  - name: Christian Frankenberg
    orcid: 0000-0002-0546-5857
    affiliation: 1
affiliations:
 - name: California Institute of Technology 
   index: 1
 - name: Jet Propulsion Laboratory 
   index: 2
date: 10 January 2022
bibliography: paper.bib
---

# Summary

Remote sensing researchers use satellite data and radiative transfer modeling to study Earth's atmospheric and surface properties. The field plays a key role in how scientists understand many aspects of our rapidly changing planet – from climate change and pollution to the carbon and water cycles.

**RadiativeTransfer.jl** is a [Julia](https://julialang.org) package that enables the fast computation of full-polarized radiative transfer simulations and atmospheric optical properties, based on the Matrix Operator Method (`@Sanghavi:2013`). Users can fully customize simulation parameters and atmospheric properties, including aerosol distributions, surface reflectance, and quadrature schemes. Independent submodules can also be imported individually; for example, **Absorption.jl** can be used for computing absorption cross-sections and **Scattering.jl** for computing scattering phase-functions. 

The Julia language provides many exciting opportunities to modernize radiative transfer software. Using the ForwardDiff.jl package (`@Revels:2016`), direct Jacobians can be calculated alongside computations using automatic differentiation, allowing for an elegant and straightforward parameter-fitting interface. Julia's multiple dispatch paradigm enables the software architecture to be clean, flexible and reusable. Additionally, optimized techniques have been implemented to speed up the package’s performance on both CPU and GPU by orders-of-magnitude compared to existing RT codes. 

**RadiativeTransfer.jl** has already been used in research projects, ranging from methane-plume simulation to atmospheric profile fitting. It has also been used in graduate-level remote sensing coursework. Ultimately, **RadiativeTransfer.jl** aims to accelerate the pace of atmospheric research through efficient software while lowering the barrier-of-entry for researchers and students in remote sensing. 

# Statement of need

For historical reasons, much of the scientific work in remote sensing is based on legacy code, written in Fortran or C/C++, mixed with “glue languages” such as Python. Researchers who developed these codes also placed greater emphasis on science results than software engineering best practices. As a result, many parts of key codebases are aging, convoluted, and hard to improve by both incoming graduate students and experienced researchers. 

Rather than simply *porting* these codes to a new language, **RadiativeTransfer.jl** entirely redesigns the radiative transfer code from the ground up to include new functionalities like GPU acceleration and automatic differentiation – features that have become computationally feasible and widespread only in the last decade. 

# Overview of functionality

The package has a modular architecture, allowing users to import just the specific module(s) that they need.

![Sample atmospheric reflectance under default atmospheric parameters, calculated using RadiativeTransfer.jl](joss_1.png)

**RadiativeTransfer.jl** is the top-level module that uses absorption and scattering submodules to compute RT simulations. Specifically, it: 
- Enables 1D vectorized plane-parallel RT modeling based on the Matrix Operator Method (`@Sanghavi:2013`)
- Incorporates fast, high fidelity simulations of scattering atmospheres containing haze and clouds, including pressure- and temperature-resolved absorption profiles of gaseous species in the atmosphere
- Enables GPU-accelerated computations of the resulting hyperspectral reflectances/transmittances
- Enables auto-differentiation of the output spectrum with respect to various input parameters

![Sample absorption spectrum of CO2 with 0.01 step size resolution, calculated using Absorption.jl](joss_2.png)

**Absorption.jl** enables absorption cross-section calculations of atmospheric gases at different pressures, temperatures, wavelengths, and broadeners (Doppler, Lorentzian, Voigt). It uses the HITRAN (`@Gordon:2017`) energy transition database for calculations. While it enables lineshape calculations from scratch, the module also allows users to create and save an interpolator object at specified wavelength, pressure, and temperature grids. The module also supports auto-differentiation (AD) of the profile, with respect to pressure and temperature. Calculations can be computed either on CPU or GPU (CUDA).

![Sample scattering phase functions of aerosols, calculated using Scattering.jl (μ = 0.3 , σ = 2.0, nᵣ = 1.3, nᵢ = 0.0, λ = 0.40 μm)](joss_3.png)

**Scattering.jl** is used for calculating Mie scattering phase-functions for aerosols with specified size distributions and refractive indices. This module enables scattering phase-function calculation of atmospheric aerosols with different size distributions, incident wavelengths, and refractive indices. It can perform the calculation using either the Siewert NAI-2 or Domke PCW methods (`@Sanghavi:2017`). The module also supports auto-differentiation (AD) of the phase function, with respect to the aerosol's size distribution parameters and its refractive index. 

# Benchmarks

Standard reference tables from the literature are used to validate **RadiativeTransfer.jl** simulation output. 

| μ     | ϕ = 0°                 | ϕ = 30°                | ϕ = 60°                | ϕ = 90°                | ϕ = 120°               | ϕ = 150°               | ϕ = 180°               |
| :---  | :---                   | :---                   | :---                   | :---                   | :---                   | :---                   | :---                   |
| 0.02  | **0.4410** <br> 0.4413 | **0.3942** <br> 0.3944 | **0.3007** <br> 0.3009 | **0.2545** <br> 0.2547 | **0.3023** <br> 0.3025 | **0.3970** <br> 0.3973 | **0.4442** <br> 0.4445 |
| 0.06  | **0.3921** <br> 0.3925 | **0.3519** <br> 0.3522 | **0.2719** <br> 0.2721 | **0.2332** <br> 0.2334 | **0.2762** <br> 0.2764 | **0.3594** <br> 0.3597 | **0.4007** <br> 0.4011 |
| 0.10  | **0.3505** <br> 0.3508 | **0.3154** <br> 0.3158 | **0.2460** <br> 0.2462 | **0.2133** <br> 0.2135 | **0.2524** <br> 0.2526 | **0.3264** <br> 0.3268 | **0.3632** <br> 0.3635 |
| 0.16  | **0.2983** <br> 0.2986 | **0.2694** <br> 0.2697 | **0.2126** <br> 0.2128 | **0.1869** <br> 0.1870 | **0.2212** <br> 0.2214 | **0.2844** <br> 0.2847 | **0.3155** <br> 0.3159 |
| 0.20  | **0.2690** <br> 0.2694 | **0.2435** <br> 0.2438 | **0.1935** <br> 0.1937 | **0.1716** <br> 0.1717 | **0.2032** <br> 0.2034 | **0.2604** <br> 0.2607 | **0.2885** <br> 0.2889 |
| 0.28  | **0.2212** <br> 0.2215 | **0.2011** <br> 0.2014 | **0.1620** <br> 0.1622 | **0.1461** <br> 0.1462 | **0.1733** <br> 0.1734 | **0.2206** <br> 0.2209 | **0.2437** <br> 0.2441 |
| 0.32  | **0.2016** <br> 0.2019 | **0.1837** <br> 0.1839 | **0.1491** <br> 0.1492 | **0.1355** <br> 0.1356 | **0.1608** <br> 0.1610 | **0.2041** <br> 0.2043 | **0.2251** <br> 0.2254 |
| 0.40  | **0.1687** <br> 0.1689 | **0.1545** <br> 0.1547 | **0.1274** <br> 0.1275 | **0.1178** <br> 0.1179 | **0.1398** <br> 0.1399 | **0.1759** <br> 0.1761 | **0.1934** <br> 0.1936 |
| 0.52  | **0.1308** <br> 0.1310 | **0.1210** <br> 0.1212 | **0.1027** <br> 0.1028 | **0.0976** <br> 0.0976 | **0.1152** <br> 0.1153 | **0.1427** <br> 0.1429 | **0.1558** <br> 0.1560 |
| 0.64  | **0.1022** <br> 0.1023 | **0.0959** <br> 0.0960 | **0.0844** <br> 0.0844 | **0.0824** <br> 0.0824 | **0.0962** <br> 0.0963 | **0.1164** <br> 0.1165 | **0.1259** <br> 0.1261 |
| 0.72  | **0.0868** <br> 0.0868 | **0.0824** <br> 0.0824 | **0.0746** <br> 0.0746 | **0.0742** <br> 0.0742 | **0.0856** <br> 0.0856 | **0.1014** <br> 0.1015 | **0.1088** <br> 0.1089 |
| 0.84  | **0.0678** <br> 0.0678 | **0.0658** <br> 0.0658 | **0.0627** <br> 0.0627 | **0.0640** <br> 0.0639 | **0.0716** <br> 0.0716 | **0.0812** <br> 0.0812 | **0.0856** <br> 0.0856 |
| 0.92  | **0.0578** <br> 0.0578 | **0.0572** <br> 0.0571 | **0.0565** <br> 0.0564 | **0.0582** <br> 0.0581 | **0.0630** <br> 0.0630 | **0.0684** <br> 0.0684 | **0.0709** <br> 0.0709 |
| 0.96  | **0.0539** <br> 0.0539 | **0.0538** <br> 0.0537 | **0.0540** <br> 0.0539 | **0.0555** <br> 0.0555 | **0.0586** <br> 0.0586 | **0.0619** <br> 0.0619 | **0.0633** <br> 0.0632 |
| 0.98  | **0.0524** <br> 0.0524 | **0.0525** <br> 0.0525 | **0.0530** <br> 0.0529 | **0.0543** <br> 0.0542 | **0.0563** <br> 0.0563 | **0.0583** <br> 0.0583 | **0.0591** <br> 0.0591 |
| 1.00  | **0.0530** <br> 0.0530 | **0.0530** <br> 0.0530 | **0.0530** <br> 0.0530 | **0.0530** <br> 0.0530 | **0.0530** <br> 0.0530 | **0.0530** <br> 0.0530 | **0.0530** <br> 0.0530 |
: I (Upwelling at TOA) for τ = 0.5, μ0 = 0.2, and A = 0.0

This validation demonstrates that simulation output using RadiativeTransfer.jl closely matches the published standard values. 

Runtime duration for a given simulation is also compared between using CPU and GPU architectures. (CPU architecture is single-threaded, AMD EPYC 7H12 64-Core Processor; GPU is parallel on an NVIDIA A100 Tensor Core (40Gb))

![](agu_4.png)

A nearly 100x speedup is observed when using the GPU architecture on the same simulation, suggesting that RadiativeTransfer.jl has the potential to greatly accelerate the pace of scientific discovery in the remote sensing field. 

# Acknowledgements

We thank Frankenberg lab members for their enthusiastic support and guidance throughout this project. We also acknowledge support from Caltech’s Schmidt Academy for Software Engineering.

# References

