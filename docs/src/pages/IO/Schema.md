# Configuration schema (IO)

This page documents the expected configuration structure accepted by `read_parameters`, `parameters_from_file`, `parameters_from_dict`, and `parameters_from_source`. YAML, TOML, and in-memory `Dict` inputs use the same schema.

Top-level keys:
- radiative_transfer (required)
- geometry (required)
- atmospheric_profile (required)
- absorption (optional)
- scattering (optional)

## radiative_transfer
- spec_bands: Vector{String}
  - Each string is a Julia range expression like "(1e7/777):0.015:(1e7/757)".
    Units accepted; values converted to cm⁻¹. Use explicit units when not `cm⁻¹`.
  - **Note**: This field uses `eval()` internally to support flexible range syntax and Unitful conversions.
    All other fields use safe parsing without eval.
- surface: Vector{String}
  - Each is a surface constructor call-like string, parsed safely (no eval):
    - LambertianSurfaceScalar(0.15)
    - LambertianSurfaceSpectrum([0.12, 0.13])
    - LambertianSurfaceLegendre([0.2, 0.05, 0.01])
    - rpvSurfaceScalar(ρ₀, ρ_c, k, Θ)
    - RossLiSurfaceScalar(fvol, fgeo, fiso)
    - CoxMunkSurface(wind_speed=U)
  - `CoxMunkSurface` supports full-polarization ocean BRDFs. Optional
    keywords include `n_water`, `whitecap_albedo` (default `0.22`),
    `include_whitecaps` (default `true`), and `shadowing` (default `true`).
  - Vegetation canopies are configured through a top-level `canopy` section;
    the parser wraps each band surface as the soil BRDF inside a
    `CanopySurface`.
- quadrature_type: String in {RadauQuad(), GaussLegQuad()}
- polarization_type: String in {Stokes_I(), Stokes_IQ(), Stokes_IQU(), Stokes_IQUV()}
- max_m: Int
- Δ_angle: Real (degrees)
- l_trunc: Int
- depol: Real
- float_type: String in {Float32, Float64}
- architecture: String in {default_architecture, CPU(), GPU(), Architectures.CPU(), Architectures.GPU()}

Accuracy and cost are controlled mainly through `quadrature_type`, `max_m`,
`Δ_angle`, `l_trunc`, and `depol`. Increase quadrature order and Fourier
moments for strongly anisotropic or highly polarized scenes; increase
`l_trunc` when aerosol phase functions have sharp forward peaks. Use consistent
units for spectral bands and optical inputs.

## geometry
- sza: Real (deg)
- vza: Vector{Real} (deg)
- vaz: Vector{Real} (deg)
- obs_alt: Real (Pa)
  - **Note**: Currently not used. Only Top Of Atmosphere (TOA/satellite) observations are supported at the moment.

## atmospheric_profile
- T: Vector{Real} (K), mid-level, TOA→BOA
- p: Vector{Real} (hPa), boundaries, TOA→BOA (length(T)+1)
- q: Vector{Real} (optional)
- profile_reduction: Int, or -1 for no reduction

## absorption (optional)
- molecules: Vector{Vector{String}} (per band)
- vmr: Dict{String,Real|Vector}
- broadening: String in {Voigt(), Lorentz(), Doppler()}
- CEF: String in {HumlicekWeidemann32SDErrorFunction()}
- wing_cutoff: Real (cm⁻¹)
- LUTfiles: Vector{Vector{String}} (optional, per molecule per band). Paths may begin with `${ENV:NAME}` to resolve large local LUT directories outside the repository.

## scattering (optional)
- aerosols: Vector of Dicts. Two aerosol forms are supported:
  - Mie/lognormal aerosol:
    - Required optical keys: `τ_ref`, `μ`, `σ`, `nᵣ`, `nᵢ`
    - Required vertical-profile keys: either `p₀`, `σp` or `z₀`, `σ₀`
  - Analytic phase-function aerosol:
    - Required keys: `τ_ref`, `phase_function`
    - Optional keys: `ssa` or `ϖ` for single-scattering albedo
      (default `1`)
    - Required vertical-profile keys: either `p₀`, `σp` or `z₀`, `σ₀`
    - Optional Mie/lognormal keys: `μ`, `σ`, `nᵣ`, `nᵢ`; if one is
      provided, all four must be provided. They are retained as metadata, but
      the analytic `phase_function` supplies the scattering optics.
- r_max: Real (µm)
- nquad_radius: Int
- λ_ref: Real (µm)
- decomp_type: String in {NAI2(), PCW()}
- n_ref: Complex (optional). If not provided, computed from first aerosol.

Analytic `phase_function` may be written as a constructor-like string:

```yaml
phase_function: HenyeyGreensteinPhaseFunction(g=0.7)
```

or as a mapping:

```yaml
phase_function:
  type: SyntheticPolarizedHenyeyGreensteinPhaseFunction
  g: 0.45
  polarization_fraction: 0.6
```

Supported analytic phase-function types are
`HenyeyGreensteinPhaseFunction` and
`SyntheticPolarizedHenyeyGreensteinPhaseFunction`. They are converted to
`GreekCoefs` and use the same CoreRT layer-optics path as Mie-derived
aerosols, so they are available to both full MOM and StandaloneSS. Linearized
Mie-parameter Jacobians are currently defined only for Mie/lognormal aerosols;
analytic phase-function aerosols throw in `LinMode` until their parameter
layout is defined.

## Minimal YAML example

```yaml
radiative_transfer:
  spec_bands:
    - (1e7/777):0.015:(1e7/757)
  surface:
    - LambertianSurfaceScalar(0.15)
  quadrature_type: GaussLegQuad()
  polarization_type: Stokes_I()
  max_m: 3
  Δ_angle: 2.0
  l_trunc: 20
  depol: 0.0
  float_type: Float64
  architecture: default_architecture

geometry:
  sza: 60
  vza: [60, 45, 30, 15, 0, 15, 30, 45, 60]
  vaz: [180, 180, 180, 180, 0, 0, 0, 0, 0]
  obs_alt: 1000.0

atmospheric_profile:
  T: [260, 262, 264]
  p: [0.1, 100.0, 800.0, 1005.0]
  profile_reduction: -1
```
