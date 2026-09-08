# Source configuration

Sources are configured programmatically through `sources=` on
`model_from_parameters` or `rt_run`. A YAML/TOML `sources:` list is not
implemented. The default is `SolarBeam()` with unit Stokes-I irradiance.

## Runnable example

```julia
using vSmartMOM
using vSmartMOM.CoreRT

params = read_parameters(joinpath(pkgdir(vSmartMOM), "config", "quickstart.yaml"))
model = model_from_parameters(params)
ν = collect(model.atmosphere.spec_bands[1])
pol_n = params.polarization_type.n

# A collimated external source with a Planck-shaped spectrum.
beam = BlackbodySource(1500.0, ν; pol_n=pol_n)
R, T = rt_run(model; sources=beam)

# A prescribed hemispheric surface-emission spectrum over a Lambertian surface.
sif = SurfaceSIF(SIF₀=fill(0.01, pol_n, length(ν)))
R_sif, T_sif = rt_run(model; sources=beam + sif)
```

`BlackbodySource` requires a temperature and a vector of wavenumbers as
positional arguments. Set `pol_n` to match the scene (its default is 3).
It returns a `SolarBeam`; it does not emit from the atmosphere or surface.
`F₀ = factor * scale * B(ν,T)`, with `factor=π` and `scale=1` by default.

## Supported vocabulary

| Type | Physical meaning | Current scope |
|---|---|---|
| `SolarBeam(; F₀=nothing, sza=nothing)` | Collimated incident Stokes irradiance, matrix `(nStokes, nSpec)` | Uses model geometry; `sza` is reserved metadata and does not change it |
| `BlackbodySource(T, ν; pol_n=3, factor=π, scale=1)` | Planck-shaped collimated `SolarBeam` | Same restrictions as `SolarBeam` |
| `SurfaceSIF(; SIF₀=nothing)` | Prescribed hemispheric surface emission | Lambertian scalar/spectrum/Legendre/spline surfaces; no emission is injected on other surfaces |
| `SurfaceSIF(; SIF760, mSIF=0, wavenumber_cm1)` | Retrievable radiance `SIF760 + mSIF*(ν-1e7/760)` | Two source columns; see [Jacobians](../../jacobians.md) |
| `SurfaceSIF(; SIF755, slope=0, wavelength_nm)` | Retrievable radiance `SIF755 + slope*(λ_nm-755)` | Alternate coordinate convention; do not interchange the two slopes |
| `ThermalEmission(T_layers, ν)` | Atmospheric volume Planck emission | Forward endpoint RT; thermal Jacobians and interior observers are unsupported |
| `NoSource()` | No incident/emitted source | Also the identity for source composition |
| `SourceSet((s₁, s₂, ...))` | Ordered source composition, also built with `+` | Use at most one solar/blackbody beam; see the limitations below |

`DiffuseBoundary` and `LidarPulse` are reserved extension types, not complete
RT capabilities. Source AD traits describe an extension seam; they do not
promise that all downstream tangent kernels exist.

## Units and coordinates

`SolarBeam.F₀` and prescribed `SurfaceSIF.SIF₀` are spectral irradiances.
With inputs in mW per m² per cm⁻¹, output Stokes values are spectral radiances
in mW per m² per sr per cm⁻¹. They are not dimensionless reflectances.
`ThermalEmission.B_layer` and retrievable SIF amplitudes are already radiances.
The source preparer converts retrievable SIF to hemispheric irradiance with π.

The wavelength coordinate in the `SIF755` form does not convert a spectral
radiance density from per cm⁻¹ to per nm. That density conversion is a separate
operation, and must be applied consistently to all sources and outputs.

A nonzero retrievable SIF amplitude **or slope** requires an explicitly supplied,
non-unit `SolarBeam.F₀` spectrum. The guard checks that a spectrum is supplied
and is not identically one; the caller must supply the physically calibrated
Fraunhofer structure. It does not validate calibration or spectral structure.

## Composition limits

`+` flattens source tuples, but the current solar extraction uses only the
first beam. Do not compose two solar/blackbody beams; for the same geometry,
sum their irradiance matrices into one `SolarBeam(F₀=F₁+F₂)`.

Prescribed SIF on non-Lambertian surfaces currently reaches a no-op surface
method. Treat that combination as unsupported, even if the solve completes.
Strict-interior observers reject surface emission explicitly. Use the
[full source guide](../../extending/sources.md) for the supported tangent and
surface extension hooks.

The main `rt_run` path ignores legacy `RS_type.SIF₀`; configure `SurfaceSIF`
through `sources=`. Historical diagnostic entry points may retain the old field.

To change solar geometry, set `params.sza` before model construction or use
`remake_geometry`. The `SolarBeam.sza` field has no numerical effect today.
