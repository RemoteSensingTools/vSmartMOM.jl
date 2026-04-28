# vSmartMOM.jl

**For:** first-time users, retrieval developers, and method reviewers.

**Next:** [Quick Start](pages/quickstart.md), [Configure a Scene](pages/IO/Overview.md), [Core RT Theory](pages/vSmartMOM/CoreRTTheory.md).

vSmartMOM.jl is a polarized atmospheric radiative-transfer solver based on the vector matrix-operator adding-doubling method. It supports gas absorption, aerosol scattering, Raman/inelastic source terms, flexible surface models, optional CUDA acceleration, and analytic Jacobians for retrieval work.

## Install

vSmartMOM supports Julia 1.10 or later.

```julia
pkg> add vSmartMOM
```

## Smallest Complete Run

```julia
using vSmartMOM

params = read_parameters(joinpath(pkgdir(vSmartMOM), "config", "quickstart.yaml"))
model = model_from_parameters(params)
R, T = rt_run(model)
```

`R` is top-of-atmosphere reflectance and `T` is bottom-of-atmosphere transmittance, both shaped as view zenith angle × Stokes component × spectral point. The [Quick Start](pages/quickstart.md) walks through the same CPU scene step by step.

## Where To Go

- New to the package: start with [Quick Start](pages/quickstart.md).
- Changing a scene: use [Configure a Scene](pages/IO/Overview.md), then the [schema reference](pages/IO/Schema.md).
- Computing retrieval sensitivities: go to [Compute Jacobians](pages/jacobians.md).
- Running on CUDA: go to [Run on GPU](pages/gpu.md).
- Checking the method against the papers: go to [Core RT Theory](pages/vSmartMOM/CoreRTTheory.md) and [References](pages/vSmartMOM/References.md).
- Extending internals: start with [Add a Surface BRDF](pages/extending/surfaces.md) or [Add a Raman Mode](pages/extending/raman.md).

## Public Modules

- **CoreRT**: adding-doubling solver, model types, optical-property assembly, surface coupling, and Jacobian kernels.
- **IO**: YAML, TOML, Dict, NetCDF, and GEOS-Chem inputs.
- **[Absorption](pages/Absorption/Overview.md)**: HITRAN and lookup-table gas absorption.
- **[Scattering](pages/Scattering/Overview.md)**: Mie calculations, phase functions, Greek coefficients, and truncation inputs.
- **[Surfaces](pages/Surfaces/Overview.md)**: lower-boundary BRDF models and canopy coupling.
- **[InelasticScattering](pages/Inelastic/Overview.md)**: Raman/Cabannes mode types and optical-property helpers.
- **[Aerosols](pages/Aerosols/Overview.md)**: TOMAS-15 and two-moment aerosol input support. This API is still being stabilized.
- **[SolarModel](pages/SolarModel/Overview.md)**: solar/stellar spectra and transmission helpers.
