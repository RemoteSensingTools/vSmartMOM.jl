```@meta
EditURL = "Tutorial_Canopy.jl"
```

# Canopy-Coupled Surface Tutorial

vSmartMOM supports vegetation canopy as a lower boundary condition via `CanopySurface`.
The canopy replaces the usual surface BRDF and internally runs canopy sub-layers
through the adding-doubling method before combining with the soil reflectance.

## 1) How it works

`CanopySurface` is an `AbstractSurfaceType` that wraps:
- A **soil BRDF** (e.g., `LambertianSurfaceScalar`)
- **Leaf area index (LAI)** controlling canopy density
- **Leaf angle distribution (LAD)** from `CanopyOptics.jl`
- **Leaf reflectance/transmittance** (scalar or spectral vector)
- Optional **spectral grid** for wavelength-dependent leaf optics
- Optional **within-canopy atmosphere** via `canopy_dp` pressure thickness
- Optional **multi-layer** decomposition (1 = big-leaf, >1 = sub-layers)

When `rt_run()` reaches the surface, dispatch on `CanopySurface` internally:
1. Precomputes azimuthal Z-matrices from the canopy scattering model
2. Runs canopy sub-layers through elemental → doubling
3. Applies interaction with the soil BRDF
4. Returns the effective canopy+soil reflectance to the atmospheric RT

```julia
using vSmartMOM
using vSmartMOM.CoreRT
using CanopyOptics
using CairoMakie
```

## 2) Scalar leaf optics (simplest setup)

A big-leaf canopy with constant leaf reflectance and transmittance over a
Lambertian soil:

```julia
soil = LambertianSurfaceScalar(0.1)
canopy_scalar = CanopySurface(
    soil = soil,
    LAI = 3.0,
    n_layers = 1,
    leaf_reflectance = 0.45,
    leaf_transmittance = 0.05,
)
```

## 3) Spectral leaf optics from PROSPECT

The PROSPECT leaf model computes reflectance and transmittance as a function
of leaf biochemistry.  `CanopySurface_from_prospect` wraps this into a
spectrally-resolved canopy surface.

```julia
leaf = CanopyOptics.LeafProspectProProperties(
    N   = 1.4,      # leaf structure parameter
    Ccab = 40.0,    # chlorophyll a+b [μg/cm²]
    Ccar = 8.0,     # carotenoids [μg/cm²]
    Canth = 0.0,    # anthocyanins
    Cbrown = 0.0,   # brown pigments
    Cw = 0.01,      # equivalent water thickness [cm]
    Cm = 0.009,     # dry matter content [g/cm²]
    Cprot = 0.0,    # protein content [g/cm²]
    Ccbc = 0.0,     # carbon-based constituents [g/cm²]
)

canopy_prospect = CanopySurface_from_prospect(
    leaf, 400.0:1.0:2500.0;
    soil = LambertianSurfaceScalar(0.1),
    LAI = 3.0,
    n_layers = 2,
)

println("Spectral grid points: ", length(canopy_prospect.leaf_optics_grid))
println("Grid unit: ", canopy_prospect.grid_unit)
println("Sample R (at grid point 200): ", canopy_prospect.leaf_reflectance[200])
println("Sample T (at grid point 200): ", canopy_prospect.leaf_transmittance[200])
```

## 4) Custom spectral vectors

You can provide your own spectral leaf R/T on any wavelength or wavenumber grid:

```julia
canopy_custom = CanopySurface(
    soil = LambertianSurfaceScalar(0.1),
    LAI = 3.0,
    n_layers = 1,
    leaf_reflectance   = [0.05, 0.08, 0.10, 0.45, 0.48],
    leaf_transmittance = [0.02, 0.03, 0.05, 0.40, 0.42],
    leaf_optics_grid   = [550.0, 660.0, 680.0, 750.0, 780.0],  # nm
    grid_unit = :nm,
)
```

## 5) Within-canopy atmosphere

For tall canopies (e.g. forests), there is an atmospheric sub-column within the
canopy where gas absorption matters.  Set `include_atm=true` and `canopy_dp` to
the pressure thickness (in hPa) of the canopy air column:

```julia
canopy_atm = CanopySurface(
    soil = LambertianSurfaceScalar(0.1),
    LAI = 5.0,
    n_layers = 4,
    leaf_reflectance = 0.45,
    leaf_transmittance = 0.05,
    include_atm = true,
    canopy_dp = 3.0,   # ~30m forest canopy ≈ 3 hPa
)
```

The within-canopy optical depth `_within_canopy_τ` is automatically computed
from the bottom-of-atmosphere conditions by `rt_run()`, using the same
absorption models as the atmospheric RT.

## 6) Running with the forward model

```julia
yaml_path = joinpath(dirname(dirname(pathof(vSmartMOM))),
                     "test", "test_parameters", "PureRayleighParameters.yaml")
params = parameters_from_yaml(yaml_path)
params.architecture = vSmartMOM.Architectures.CPU()
params.max_m = 2
params.l_trunc = 20

params.brdf[1] = canopy_scalar

model = model_from_parameters(params)
R_canopy, T_canopy = rt_run(model)

println("R shape: ", size(R_canopy))
println("R(nadir, I) with canopy: ", R_canopy[1, 1, 1])
```

Compare with bare Lambertian surface:

```julia
params2 = parameters_from_yaml(yaml_path)
params2.architecture = vSmartMOM.Architectures.CPU()
params2.max_m = 2
params2.l_trunc = 20
model2 = model_from_parameters(params2)
R_bare, T_bare = rt_run(model2)
println("R(nadir, I) bare soil:   ", R_bare[1, 1, 1])
println("Canopy effect on TOA R:  ",
    round((R_canopy[1,1,1] - R_bare[1,1,1]) / R_bare[1,1,1] * 100, digits=1), "%")
```

Compare the reflectance spectra:

```julia
fig = Figure(size=(700, 450))
ax = Axis(fig[1,1],
    xlabel = "Spectral index",
    ylabel = "TOA Reflectance (Stokes I)")
lines!(ax, R_canopy[1, 1, :], label="Canopy (LAI=3)")
lines!(ax, R_bare[1, 1, :],   label="Bare soil (α=0.1)")
axislegend(ax, position=:rt)
fig
```

## 7) Effect of within-canopy atmosphere

The within-canopy atmosphere adds gas absorption between canopy sub-layers.
This matters for tall canopies and strong absorption bands. Here we compare
the TOA reflectance with and without the canopy atmosphere, using a multi-layer
canopy so the interleaved atmospheric layers have an effect.

```julia
canopy_no_atm = CanopySurface(
    soil = LambertianSurfaceScalar(0.1),
    LAI = 5.0,
    n_layers = 4,
    leaf_reflectance = 0.45,
    leaf_transmittance = 0.05,
    include_atm = false,
)

canopy_with_atm = CanopySurface(
    soil = LambertianSurfaceScalar(0.1),
    LAI = 5.0,
    n_layers = 4,
    leaf_reflectance = 0.45,
    leaf_transmittance = 0.05,
    include_atm = true,
    canopy_dp = 3.0,   # ~30m forest ≈ 3 hPa
)
```

Run without canopy atmosphere:

```julia
params_na = parameters_from_yaml(yaml_path)
params_na.architecture = vSmartMOM.Architectures.CPU()
params_na.max_m = 2
params_na.l_trunc = 20
params_na.brdf[1] = canopy_no_atm
model_na = model_from_parameters(params_na)
R_no_atm, _ = rt_run(model_na)
invalidate_canopy_cache!(model_na.params.brdf[1])
```

Run with canopy atmosphere:

```julia
params_wa = parameters_from_yaml(yaml_path)
params_wa.architecture = vSmartMOM.Architectures.CPU()
params_wa.max_m = 2
params_wa.l_trunc = 20
params_wa.brdf[1] = canopy_with_atm
model_wa = model_from_parameters(params_wa)
R_with_atm, _ = rt_run(model_wa)
invalidate_canopy_cache!(model_wa.params.brdf[1])

println("R(nadir, I) without canopy atm: ", R_no_atm[1, 1, 1])
println("R(nadir, I) with canopy atm:    ", R_with_atm[1, 1, 1])
diff_pct = (R_with_atm[1,1,1] - R_no_atm[1,1,1]) / R_no_atm[1,1,1] * 100
println("Within-canopy atmosphere effect: ", round(diff_pct, digits=3), "%")
```

For a pure Rayleigh atmosphere (no gas absorption), the effect is negligible.
In bands with strong molecular absorption (O₂ A-band, CO₂), the canopy
atmosphere can change the effective surface reflectance by several percent,
which matters for trace-gas retrievals.

## 8) YAML configuration

The canopy can also be configured via YAML, including spectral properties:

```yaml
radiative_transfer:
  surface:
    - LambertianSurfaceScalar(0.1)

canopy:
  LAI: 3.0
  n_layers: 4
  leaf_reflectance: [0.05, 0.10, 0.45, 0.48]
  leaf_transmittance: [0.02, 0.05, 0.40, 0.42]
  leaf_optics_grid: [660.0, 680.0, 750.0, 780.0]
  grid_unit: nm
  include_atm: true
  canopy_dp: 3.0
```

The parser wraps the surface entry as the soil BRDF inside a `CanopySurface`.

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

