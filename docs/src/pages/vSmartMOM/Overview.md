# vSmartMOM Module Overview

The vSmartMOM module allows end-to-end simulation of radiative transfer (RT) throughout Earth's atmosphere as well as atmospheres of stars and substellar objects like brown dwarfs and exoplanets. Specifically, the module:

- Enables 1D vectorized plane-parallel RT modeling based on the Matrix Operator Method, also known as Discrete Space Theory (please see the [references](https://remotesensingtools.github.io/vSmartMOM.jl/dev/pages/vSmartMOM/References/) for source papers). 
- Incorporates fast, high-fidelity simulations of scattering atmospheres containing haze and clouds, including pressure- and temperature-resolved absorption profiles of gaseous species in the atmosphere
- Enables GPU-accelerated computations of the resulting hyperspectral reflectances/transmittances
- Enables auto-differentiation of the output spectrum with respect to various input parameters

You can perform an RT simulation in a few short steps: 

1. Use [`read_parameters`](@ref) or [`default_parameters`](@ref) to load a set of RT input parameters. You can modify any parameter in the returning [`CoreRT.vSmartMOM_Parameters`](@ref) struct (it is mutable). Please see the parameters [guide](https://remotesensingtools.github.io/vSmartMOM.jl/dev/pages/vSmartMOM/InputParametersGuide/) for more information on what each field specifies.
2. Use [`model_from_parameters`](@ref) to calculate derived fields (ex. aerosol optical parameters, quadrature points, etc.). The output [`CoreRT.RTModel`](@ref) contains hierarchical sub-structs for solver configuration, atmosphere, optics, and surfaces.
3. Use [`rt_run`](@ref) to perform the radiative transfer calculation using the defined model settings. 

For a full demo of how to use this module, please see the [example](https://remotesensingtools.github.io/vSmartMOM.jl/dev/pages/vSmartMOM/Example/) page. 

For equation-level details of the core solver workflow (`elemental`, `doubling`, and `interaction`), see the **Core RT Theory (Doubling/Adding)** page in this section.

## Architecture

![ArchitectureDiagram](vSmartMOMDiagram-vSmartMOM.drawio.png)
