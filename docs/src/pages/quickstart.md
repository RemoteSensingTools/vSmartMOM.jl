# Quick Start

**For:** first-time users who want one successful forward radiative-transfer run.

**Next:** [Configure a Scene](IO/Overview.md), [Compute Jacobians](jacobians.md), [RT basics](concepts/01_overview.md).

This page is the 5-minute CPU path. It uses a tiny shipped scene so the first run does not require external line-data downloads or aerosol optics setup. The long-form tutorial remains available at [Tutorial: Quick Start](tutorials/Tutorial_QuickStart.md).

## Run The Shipped Quickstart Scene

The quickstart scene is deliberately small: one spectral point, one Stokes component, two atmospheric layers, a Lambertian surface, and no absorption or aerosols. That keeps the first run focused on the core Rayleigh adding-doubling path.

```julia
using vSmartMOM

scene = joinpath(pkgdir(vSmartMOM), "config", "quickstart.yaml")
params = read_parameters(scene)

model = model_from_parameters(params)
R, T = rt_run(model)
```

`params` is a mutable scene description. For example, this file already sets `architecture: CPU()`, but an interactive session can still override it before model construction:

```julia
params.architecture = vSmartMOM.Architectures.CPU()
```

## Inspect The Output

`R` and `T` are three-dimensional arrays:

```julia
size(R)
size(T)
R[:, 1, :]
```

The dimensions are:

- view zenith angle
- Stokes component
- spectral point

For `config/quickstart.yaml`, both arrays have shape `(1, 1, 1)`. Larger scenes increase these dimensions by adding more viewing angles, polarization components, or spectral samples.

## Next Steps

- Change surface, geometry, atmosphere, or spectral settings in [Configure a Scene](IO/Overview.md).
- Use the same `model_from_parameters` / `rt_run` workflow on a richer example in [Tutorial: Quick Start](tutorials/Tutorial_QuickStart.md).
- Read [The MOM Solver](concepts/04_mom_solver.md) for the elemental, doubling, and adding operations behind the run.
