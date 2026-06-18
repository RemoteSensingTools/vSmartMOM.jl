# Setup-dependent tests & configs

These are **not** part of the CI unit suite (`test/runtests.jl`). They need
local data, external LUTs, or features not yet released — i.e. they can only be
run by someone who has the prerequisites on their own machine.

## Run

```bash
julia --project=test test/setup/runtests.jl
```

Each config is parse + forward + linearized-built; whatever builds in your
environment is checked for forward/lin Fourier-loop parity, and the rest report
what they need and are skipped.

## `test_parameters/` here

| Group | Configs | Needs |
|-------|---------|-------|
| H2O override | `2BandParameters`, `3BandParameters`, `3BandParameters_canopy`, `S2BandParameters`, `ParamsEMIT_MODTRANcomp`, `ParamsEMIT_MODTRANcomp_alb02` | a **config-level H2O override** (H2O in the molecule list + explicit VMR). Not yet implemented — H2O is currently driven by `atmospheric_profile.q`, and listing H2O is rejected. Parked here until the override lands. |
| External LUTs | `ParamsEMIT_MODTRANcomp_newLUT`, `…_newLUT_alb02`, `…_newLUT_alb02_azim180`, `…_newLUT_azim180` | the ABSCO LUTs via the `VSMARTMOM_HITRAN_LUT_DIR` environment variable (data not in the repo). |
| Legacy / reference | `ThreeBandsParameters`, `O2ParametersVS`, `O2_parameters2_SIF_grid`, `WCO2_parameters_SIF_grid` | kept as reference scenes; some predate the current schema. |

Why they moved here: the main `test_lin_forward_loop_parity.jl` used to **silently
skip** any config that failed to parse/build, which let these rot unnoticed. The
main parity test now **fails loudly** on any unbuildable config in
`test/test_parameters/`, so anything that can't build in CI lives here instead.

> Note: AtmosphericAbsorption.jl carries its own LBL/absorption unit tests, so
> the LUT-heavy absorption scenes are intentionally not duplicated in the main
> vSmartMOM CI suite.
