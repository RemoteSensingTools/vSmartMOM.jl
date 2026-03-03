# Introduction

vSmartMOM.jl is a polarized radiative transfer solver that uses the adding-doubling method to compute reflectance and transmittance for atmospheric remote sensing. It supports gas absorption, aerosol scattering, and analytic Jacobians for inversion. The interface is designed for researchers and students, with optional GPU acceleration via CUDA.jl.

## System requirements

- **Julia** 1.9 or later
- **Optional:** NVIDIA GPU with CUDA for `vSmartMOM.Architectures.GPU()` support

## Quick example

```julia
using vSmartMOM
params = parameters_from_yaml("path/to/your/params.yaml")
model = model_from_parameters(params)
R, T = rt_run(model)   # R = TOA reflectance, T = BOA transmittance
```

See [Tutorial: Quick Start](pages/tutorials/Tutorial_QuickStart.md) for a complete walkthrough.

## Installation

1. Install [Julia](https://julialang.org/downloads/) (1.9+)
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