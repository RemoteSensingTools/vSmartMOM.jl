# Library

**For:** users who already know which vSmartMOM object or function they need.

**Next:** start with [Top-Level API](api/top_level.md) for ordinary runs, or jump
directly to the module that owns the feature you are using.

The task pages explain workflows. The library pages below group the public API
by module and keep the long docstring reference out of the main reading path.

## Core Interfaces

- [Top-Level API](api/top_level.md): scene loading, model construction, solver entry points, architectures.
- [CoreRT](api/core_rt.md): solver modes, model containers, parsed parameter containers, and Jacobian layout helpers.
- [IO](api/io.md): file, dictionary, GEOS-Chem, and NetCDF input sources.

## Physics Components

- [Absorption](api/absorption.md): HITRAN data, line-shape models, and absorption cross sections.
- [Scattering](api/scattering.md): Mie models, Fourier decomposition, phase reconstruction, and aerosol optics.
- [Surfaces](api/surfaces.md): Lambertian, land BRDF, ocean glint, and canopy lower boundaries.
- [Inelastic Scattering](api/inelastic.md): Raman and Cabannes mode types and source helpers.
- [Aerosols](api/aerosols.md): TOMAS-15 / two-moment aerosol input helpers. This API is still stabilizing.
- [SolarModel](api/solar_model.md): solar spectra and transmission helpers.

## Appendix

- [Experimental Helpers](api/experimental.md): exported helpers that are available but not yet a documented end-to-end product workflow.
- [Developer Coverage](internal_api_coverage.md): exported developer-facing helpers kept documented for `checkdocs = :exports`.
- [Function Index](api/function_index.md): alphabetical index generated from the documented library pages.
