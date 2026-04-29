# Compute Jacobians

**For:** retrieval and inversion developers who need analytic sensitivities.

**Next:** [Core RT Theory](vSmartMOM/CoreRTTheory.md), [Library](api_reference.md), [Tutorial: Jacobians](tutorials/Tutorial_Jacobians.md).

vSmartMOM has a linearized radiative-transfer path for Jacobians. It propagates derivative optical properties through the same adding-doubling solver used for the forward model.

## Build A Linearized Model

Use `model_from_parameters_lin(params)` when the call site should make the Jacobian path explicit. It is a convenience alias for `model_from_parameters(LinMode(), params)`.

```julia
using vSmartMOM
using vSmartMOM.CoreRT

params = read_parameters(joinpath(pkgdir(vSmartMOM),
                                  "test", "test_parameters",
                                  "JacobianTestFast.yaml"))
params.architecture = vSmartMOM.Architectures.CPU()

model, lin_model = model_from_parameters_lin(params)
```

`model` is the forward model used by `rt_run(model)`. `lin_model` carries derivative optical properties, including gas absorption derivatives, aerosol optical-property derivatives, and layer optical-depth derivatives.

## Run Linearized RT

The linearized solver needs the number of aerosol, gas, and surface retrieval parameters so it can size the final Jacobian dimension.

```julia
NAer = isnothing(params.scattering_params) ? 0 : length(params.scattering_params.rt_aerosols)
NGas = size(lin_model.τ̇_abs[1], 1)
NSurf = 1

R, T, dR, dT = rt_run_lin(model, lin_model, NAer, NGas, NSurf)
```

The outputs are:

- `R`: reflected Stokes field, shaped `nVZA × nStokes × nSpec`
- `T`: transmitted Stokes field, shaped like `R`
- `dR`: derivative of `R`, shaped `nVZA × nStokes × nSpec × nParams`
- `dT`: derivative of `T`, shaped like `dR`

## Slice The Parameter Dimension

Use [`ParameterLayout`](@ref) rather than hard-coded offsets.

```julia
layout = CoreRT.ParameterLayout(aerosol_params = 7,
                                n_aerosols = NAer,
                                n_gases = NGas,
                                n_surface = NSurf)

surface_idx = CoreRT.surface_index(layout)
gas_block = CoreRT.gas_range(layout)
aerosol1_block = NAer > 0 ? CoreRT.aerosol_range(layout, 1) : 1:0

dR_surface = dR[:, 1, :, surface_idx]
```

The current ordering is aerosol parameters first, then gas VMR parameters, then surface parameters. Each aerosol contributes seven derivative slots: `τ_ref`, `nᵣ`, `nᵢ`, mean radius, geometric width, vertical-location parameter, and vertical-width parameter.

For a longer walkthrough with plots and finite-difference checks, use [Tutorial: Jacobians](tutorials/Tutorial_Jacobians.md). The hybrid AD notes are in [Tutorial: HybridAD](tutorials/Tutorial_HybridAD.md).
