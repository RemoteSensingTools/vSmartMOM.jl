## Introduction

vSmartMOM.jl aims to revamp and modernize key atmospheric remote sensing tools. Specifically, it will enable the fast computation of atmospheric optical properties, full-polarized radiative transfer simulations, and commonly-used inversion routines. 

By taking advantage of modern software tools, such as GPU acceleration and HPC computing, the software suite significantly accelerates computationally intensive calculations and models, while keeping the interface easy to use for researchers and students.

## Installation

1. Install [Julia](https://julialang.org/downloads/) (1.6+)
2. Start the Julia REPL, enter the Pkg REPL by pressing `]`, and run:  
```julia
add vSmartMOM
```

## Submodules

vSmartMOM.jl has a modular architecture, allowing users to import just the specific module(s) that they need: 

**vSmartMOM** is the top-level module that uses absorption and scattering submodules to compute RT simulations. 
- **Absorption** is used for computing absorption cross-sections for gases over specified pressure, temperature, and wavelength grids. 
- **Scattering** is used for calculating scattering phase-functions for aerosols with specified size distributions and refractive indices. 

## How to use the documentation

Each module has a set of documentation pages, as seen on the navigation bar. Each module's documentation contains a high-level **overview** of the module, an **example** of the module in use, a list of key **methods and types**, and a list of academic **references**. 

The vSmartMOM module's documentation additionally has a guide to all the user-defined RT parameters that can be modified. 

Finally, there are tutorials that combine text, equations, code, and plots to demonstrate uses of the submodules. 