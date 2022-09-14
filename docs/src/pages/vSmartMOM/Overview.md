# vSmartMOM Module Overview

The vSmartMOM module allows end-to-end simulation of radiative transfer (RT) throughout Earth's atmosphere as well as atmospheres of stars and substellar objects like brown dwarfs and exoplanets. Specifically, the module:

- Enables 1D vectorized plane-parallel RT modeling based on the Matrix Operator Method
- Incorporates fast, high-fidelity simulations of scattering atmospheres containing haze and clouds, including pressure- and temperature-resolved absorption profiles of gaseous species in the atmosphere
- Enables GPU-accelerated computations of the resulting hyperspectral reflectances/transmittances
- Enables auto-differentiation of the output spectrum with respect to various input parameters

You can perform an RT simulation in a few short steps: 

1. Use [`parameters_from_yaml`](@ref) or [`default_parameters`](@ref) to load a set of RT input parameters. You can modify any parameter in the returning [`vSmartMOM_Parameters`](@ref) struct (it is mutable). Please see the parameters [guide](https://remotesensingtools.github.io/vSmartMOM.jl/dev/pages/vSmartMOM/InputParametersGuide/) for more information on what each field specifies. 
2. Use [`model_from_parameters`](@ref) to calculate derived fields (ex. aerosol optical parameters, quadrature points, etc.). Again, any parameter in the output [`vSmartMOM_Model`](@ref) can be modified after the struct is created. 
3. Use [`rt_run`](@ref) to perform the radiative transfer calculation using the defined model settings. 

For a full demo of how to use this module, please see the [example](https://remotesensingtools.github.io/vSmartMOM.jl/dev/pages/vSmartMOM/Example/) page. 

## Architecture

![ArchitectureDiagram](vSmartMOMDiagram-vSmartMOM.drawio.png)