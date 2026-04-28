# Add a Surface BRDF

**For:** method developers adding a new surface model.

**Next:** [Surfaces tutorial](../tutorials/Tutorial_Surfaces.md), [Core RT Theory](../vSmartMOM/CoreRTTheory.md), [API Reference](../api_reference.md).

This page will become the surface-extension cookbook. It will cover:

- which subtype or constructor pattern to follow;
- how `create_surface_layer!` participates in the solver;
- how surface strings are parsed from scene configuration;
- where to register a new BRDF for YAML/TOML input;
- what focused tests are expected.

The implementation pass must keep this guide synchronized with `BRDF_MAP` in `src/IO/Parameters.jl`.
