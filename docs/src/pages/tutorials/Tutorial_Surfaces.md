```@meta
EditURL = "Tutorial_Surfaces.jl"
```

# Surface BRDF Models Tutorial

vSmartMOM supports several surface Bidirectional Reflectance Distribution Function (BRDF) models.
This tutorial covers the available types and when to use each.

## 1) Overview

Surface models define how the lower boundary reflects radiation. Available types:

| Model | Parameters | Use case |
|-------|------------|----------|
| Lambertian | scalar or spectral albedo | Simple, isotropic surfaces |
| RPV (Rahman-Pinty-Verstraete) | ρ₀, k, Θ, ρ_c | Vegetation, angular effects |
| Ross-Li | fiso, fvol, fgeo | Vegetation, kernel-based |

## 2) Lambertian Surface

The simplest model: reflectance is independent of viewing angle. Use `LambertianSurfaceScalar`
for a single albedo per band, or `LambertianSurfaceSpectrum` for wavelength-dependent albedo.

```julia
using vSmartMOM
using vSmartMOM.CoreRT
```

Scalar albedo (same for all wavelengths in the band)

```julia
lam_scalar = LambertianSurfaceScalar(0.15)
```

Spectral albedo (vector, one value per spectral point in the band)

```julia
lam_spectrum = LambertianSurfaceSpectrum([0.12, 0.13, 0.14])
```

## 3) RPV (Rahman-Pinty-Verstraete) Model

Four-parameter BRDF: ρ = ρ₀ · M(μᵢ, μᵣ, k) · F(Θ, cos g) · H(ρ_c, G)

- **ρ₀**: Overall amplitude (isotropic scaling)
- **k**: Minnaert limb-darkening exponent; k < 1 → bowl shape, k > 1 → bell shape
- **Θ**: Hot-spot parameter; Θ < 0 → backward scattering, Θ > 0 → forward scattering
- **ρ_c**: Geometric term amplitude; controls bowl shape in phase angle

```julia
rpv = rpvSurfaceScalar(0.15, 0.1, 0.7, -0.3)  # ρ₀, ρ_c, k, Θ
```

## 4) Ross-Li Model

Three-kernel linear combination: ρ = fiso·K_iso + fvol·K_vol + fgeo·K_geo

- **fiso**: Isotropic kernel (K_iso = 1) — baseline reflectance
- **fvol**: Volumetric (Ross Thick) — dense canopy scattering
- **fgeo**: Geometric (Li Sparse) — shadowing by sparse objects

```julia
rossli = RossLiSurfaceScalar(0.05, 0.03, 0.1)  # fvol, fgeo, fiso
```

## 5) Switching Surface Types in YAML

In your config, set `radiative_transfer.surface` to an array of constructor strings.
One entry per spectral band. Examples:

```yaml
radiative_transfer:
  surface:
    - LambertianSurfaceScalar(0.15)
    - LambertianSurfaceScalar(0.2)
    - rpvSurfaceScalar(0.12, 0.08, 0.75, -0.25)
    - RossLiSurfaceScalar(0.04, 0.02, 0.08)
```

Load a YAML with Lambertian and run:

```julia
yaml_path = joinpath(dirname(dirname(pathof(vSmartMOM))),
                     "test", "test_parameters", "PureRayleighParameters.yaml")
params = parameters_from_yaml(yaml_path)
params.brdf[1] = LambertianSurfaceScalar(0.2)  # Override first band
model = model_from_parameters(params)
R, T = rt_run(model)
```

## 6) When to Use Which Model

- **Lambertian**: Oceans, snow, simple land; fast; no angular dependence
- **RPV**: Vegetation, soils; captures hotspot and limb darkening; 4 params
- **Ross-Li**: Vegetation, MODIS-style kernels; widely used in remote sensing; 3 params

RPV and Ross-Li are scalar-only (no polarization). For polarized RT, use Lambertian.

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

