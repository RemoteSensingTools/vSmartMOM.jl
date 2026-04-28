# Add a Surface BRDF

**For:** method developers adding a new lower-boundary reflectance model.

**Next:** [Surfaces](../Surfaces/Overview.md), [Surfaces tutorial](../tutorials/Tutorial_Surfaces.md), [Core RT Theory](../vSmartMOM/CoreRTTheory.md), [API Reference](../api_reference.md).

Surface models are part of `CoreRT` because the lower boundary is represented as an `AddedLayer` and then interacts with the atmospheric composite layer. A new BRDF usually needs three pieces:

1. a concrete `AbstractSurfaceType`;
2. a `reflectance` or `create_surface_layer!` method that builds the surface layer;
3. a parser registration if the surface should be available from YAML/TOML scene files.

## Choose The Implementation Shape

Use the simplest hook that matches the model:

| Surface behavior | Implement | Examples |
| --- | --- | --- |
| Scalar BRDF that can be integrated over azimuth | `reflectance(surface, n, mu_i, mu_r, dphi)` plus the generic `reflectance(surface, pol_type, mu, m)` path | `rpvSurfaceScalar`, `RossLiSurfaceScalar` |
| Lambertian or spectral albedo with a simple closed form | specialized `create_surface_layer!` | `LambertianSurfaceScalar`, `LambertianSurfaceLegendre`, `LambertianSurfaceSpline` |
| Polarized ocean or other Mueller BRDF | specialized `reflectance(surface, pol_type, mu, m)` and helper kernels | `CoxMunkSurface` |
| Composite lower boundary | specialized `create_surface_layer!` that internally runs sub-layers and interacts with a soil layer | `CanopySurface` |

The generic surface layer stores upward surface reflection in `added_layer.r⁻⁺`, identity transmission in `t⁺⁺`/`t⁻⁻`, and optional direct solar source terms in `j₀⁺`/`j₀⁻` when SFI is active.

## Define The Type

Add the type near the existing surface definitions in `src/CoreRT/types.jl`.

```julia
struct MySurface{FT} <: AbstractSurfaceType
    strength::FT
    shape::FT
end
```

Keep fields concrete when possible. Surface parameters are often part of Jacobian layouts, parser conversion, and GPU array setup, so preserving the configured floating-point type matters.

## Implement Reflectance

For scalar BRDFs, follow the `rpvSurfaceScalar` pattern:

```julia
function reflectance(surface::MySurface{FT}, n, mu_i::FT, mu_r::FT, dphi::FT) where FT
    n == 1 || return zero(FT)
    return surface.strength * my_angular_shape(mu_i, mu_r, dphi, surface.shape)
end
```

The generic Fourier-moment method in `src/CoreRT/Surfaces/rpv_surface.jl` can then integrate over azimuth and return the quadrature-space reflectance matrix.

For a model with a closed-form quadrature matrix, implement `create_surface_layer!` directly. Use existing methods in `src/CoreRT/Surfaces/lambertian_surface.jl` as the reference for:

- moving arrays through `array_type(architecture)`;
- filling `r⁻⁺`, `r⁺⁻`, `t⁺⁺`, and `t⁻⁻`;
- applying the `m == 0` rule for isotropic contributions;
- adding direct-beam source terms only when `SFI` is true.

## Register Scene Parsing

Scene files are parsed without `eval`. To expose `MySurface(...)` in YAML/TOML configs:

1. add an entry to `BRDF_MAP` in `src/IO/Parameters.jl`;
2. add a `_construct_surface(::Val{:MySurface}, FT, args)` method;
3. validate arity and argument types with `_require_nargs` and `_require_config`.

Skeleton:

```julia
const BRDF_MAP = Dict{String, Function}(
    # existing entries...
    "MySurface" => (FT, args) -> _construct_surface(Val(:MySurface), FT, args),
)

function _construct_surface(::Val{:MySurface}, FT, args)
    _require_nargs("MySurface", args, 2, "strength, shape")
    return CoreRT.MySurface(FT(args[1]), FT(args[2]))
end
```

Then a scene can use:

```yaml
radiative_transfer:
  surface:
    - MySurface(0.2, 1.5)
```

## Add Focused Tests

Keep new tests small and local before adding an end-to-end RT case:

- parser: `parse_surface_str("MySurface(...)", Float32)` returns `MySurface{Float32}`;
- parser failure: bad arity throws `ArgumentError`;
- reflectance: finite, non-negative values over representative angles;
- Fourier moment: `m = 0` and at least one higher moment behave as expected;
- RT smoke: one tiny CPU scene runs with finite `R` and `T`.

Use `test/test_parameters_parser.jl`, `test/test_coxmunk.jl`, and `test/test_canopy.jl` as the closest current patterns.

## Linearized RT

If the new surface should contribute analytic Jacobians, add or extend the matching method under `src/CoreRT/Surfaces/*_lin.jl` and update `ParameterLayout` assumptions as needed. Otherwise document that the surface is forward-only and make sure linearized tests either skip it or compare through an explicit finite-difference path.
