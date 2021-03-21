# RadiativeTransfer.jl


| **Documentation**    | [![dev][docs-latest-img]][docs-latest-url]       |
|----------------------|--------------------------------------------------|
| **Unit Tests**       | [![unit tests][unit-tests-img]][unit-tests-url]  |
| **Code Coverage**    | [![codecov][codecov-img]][codecov-url]           |

[docs-latest-img]: https://img.shields.io/badge/docs-latest-blue.svg
[docs-latest-url]: https://rupeshjey.github.io/RadiativeTransfer.jl/dev/

[unit-tests-img]: https://travis-ci.org/RupeshJey/RadiativeTransfer.jl.svg?branch=master
[unit-tests-url]: https://travis-ci.org/github/RupeshJey/RadiativeTransfer.jl

[codecov-img]: https://codecov.io/gh/RupeshJey/RadiativeTransfer.jl/branch/master/graph/badge.svg
[codecov-url]: https://codecov.io/gh/RupeshJey/RadiativeTransfer.jl

**An end-to-end modular software suite for radiative transfer calculations, written in <a href="https://julialang.org">Julia</a>**

<img src='docs/images/CrossSectionGIF.gif' class='center'></img>

This project aims to revamp and modernize key atmospheric remote sensing tools. Specifically, it will enable the fast computation of atmospheric optical properties, full-polarized radiative transfer simulations, and commonly-used inversion routines.

By taking advantage of modern software tools, such as GPU acceleration and HPC computing, the software suite significantly accelerates computationally-intensive calculations and models, while keeping the interface easy-to-use for researchers and students.

## Table of Contents

- [Installation](#installation)
- [Modules](#modules)
- [Support](#support)
- [License](#license)

## Installation

1. <a href=https://julialang.org/downloads/>Install Julia</a> (1.3+)
2. Start Julia and run:  
```
] add https://github.com/RupeshJey/RadiativeTransfer.jl
```

## Modules

**Note: This section provides only a quick overview of the available modules in RadiativeTransfer.jl.**

For in-depth examples, tutorials, and implementation details, please see the complete <a href="https://rupeshjey.github.io/RadiativeTransfer.jl/dev/">Documentation</a>.

### Ready to use:

- **Absorption**: This module enables absorption cross-section calculations of atmospheric gases at different pressures, temperatures, and broadeners (Doppler, Lorentzian, Voigt). It uses the <a href=https://hitran.org>HITRAN</a> energy transition database for calculations. While it enables lineshape calculations from scratch, it also allows users to create and save an interpolator object at specified wavelength, pressure, and temperature grids. It can perform these computations either on CPU or GPU. <br><br> Key functions:

  - `read_hitran(filepath::String)`: Creates a HitranTable struct from the fixed-width HITRAN file with transitions.
  - `make_hitran_model(hitran::HitranTable, broadening::AbstractBroadeningFunction, ...)`: Create a HitranModel struct that holds all of the model parameters needed to perform a absorption cross-section (transitions, broadening type, wing_cutoff, etc.)
  - `make_interpolation_model(hitran::HitranTable, broadening::AbstractBroadeningFunction, )`: Similar to creating a HitranModel, but this will perform the interpolation at the given wavelength, pressure, and temperature grids and store the interpolator in InterpolationModel.
  - `absorption_cross_section(model::AbstractCrossSectionModel, grid::AbstractRange{<:Real}, pressure::Real, temperature::Real, ...)`: Performs an absorption cross-section calculation with the given model (HitranModel or InterpolationModel), at a given wavelength grid, pressure and temperature

- **PhaseFunction**: This module enables scattering phase-function calculation of atmospheric aerosols with different size distributions, incident wavelengths, and refractive indices. It can perform the calculation using either the Siewert NAI-2 or Domke PCW methods ([Suniti Sanghavi 2014](https://www.sciencedirect.com/science/article/pii/S0022407313004962)). <br><br> Key functions:

  - `make_univariate_aerosol(size_distribution::ContinuousUnivariateDistribution, r_max, nquad_radius::Int, nᵣ, nᵢ`: Create an aerosol object with size distribution and complex refractive index. 
  - `make_mie_model(computation_type::AbstractFourierDecompositionType, aerosol::AbstractAerosolType, λ::Real, polarization::AbstractPolarizationType, truncation_type::AbstractTruncationType, ...)`: Create a MieModel struct that holds all of the model parameters needed to perform a phase function calculation (computation type, aerosol, incident wavelength, etc. )
  - `compute_aerosol_optical_properties(model::MieModel)`: Compute the aerosol optical properties using the specified model parameters

### In development:

- **RTSimulation**

## Support

This project is being developed in the Christian Frankenberg and Paul Wennberg labs at Caltech and Suniti Sanghavi, with support from the Schmidt Academy for Software Engineering (SASE).

Please <a href="mailto:cfranken@caltech.edu,wennberg@gps.caltech.edu?cc=rjeyaram@caltech.edu"> email us</a> if you have any questions, suggestions, or contributions!

## License

MIT License

Copyright (c) 2020 Rupesh Jeyaram

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
