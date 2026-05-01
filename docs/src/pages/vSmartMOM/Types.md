# vSmartMOM Module Methods & Types

## Defining Parameters for RT Simulation

Use [`read_parameters`](@ref) as the general loader for YAML/TOML files, in-memory
dictionaries, and typed IO sources. Use [`default_parameters`](@ref) for the built-in
default scene.

The explicit helpers [`parameters_from_file`](@ref),
[`parameters_from_dict`](@ref), [`parameters_from_source`](@ref), and
[`parameters_from_yaml`](@ref) are also available when a call site should state the
input kind directly.

## Using the Parameters Object to Create a Model 

Use `model_from_parameters(params)` to construct an `RTModel` from a parsed
scene. See the [Library](../api_reference.md) for the canonical docstring
and return-shape details.

## Performing the RT Simulation

Use `rt_run(model)` for the forward solver. Linearized runs use
`model_from_parameters_lin` and `rt_run_lin`; see [Compute Jacobians](../jacobians.md).

## Types

### Parameters Type

`CoreRT.vSmartMOM_Parameters` is the parsed scene container returned by
`read_parameters`, `parameters_from_file`, and related loader helpers.

### Model Type

`RTModel` is the hierarchical model state used by the solver. New code should
prefer the `solver`, `geometry`, `quad_points`, `atmosphere`, `optics`, and
`surfaces` fields over legacy flat-field access.

### Surface Types

Surface models are subtypes of `CoreRT.AbstractSurfaceType`. The public surface
set includes Lambertian scalar/spectrum/Legendre/spline forms, RPV, Ross-Li,
Cox-Munk ocean, and canopy surfaces. See [Surfaces](../Surfaces/Overview.md)
and [Add a Surface BRDF](../extending/surfaces.md) for configuration and
extension guidance.

```@docs
CoreRT.invalidate_canopy_cache!
```

### Quadrature Types
```@docs
CoreRT.AbstractQuadratureType
CoreRT.RadauQuad
CoreRT.GaussQuadHemisphere
CoreRT.GaussQuadFullSphere

```

### Atmospheric Profile Type

`CoreRT.AtmosphericProfile` stores the layer temperature, pressure, humidity,
column-density, and gas-mixing-ratio fields consumed by `model_from_parameters`.
