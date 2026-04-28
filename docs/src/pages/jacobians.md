# Compute Jacobians

**For:** retrieval and inversion developers who need analytic sensitivities.

**Next:** [Core RT Theory](vSmartMOM/CoreRTTheory.md), [API Reference](api_reference.md), [Tutorial: Jacobians](tutorials/Tutorial_Jacobians.md).

vSmartMOM has a linearized radiative-transfer path for Jacobians. The detailed task page will document:

- `LinMode()`;
- `model_from_parameters_lin`;
- `rt_run_lin`;
- `ParameterLayout`;
- how to slice `dR` and `dT` by aerosol, gas, surface, and canopy parameter blocks.

Until this page is filled in, use the long-form [Jacobian tutorial](tutorials/Tutorial_Jacobians.md) and the [`ParameterLayout`](@ref) API docs.
