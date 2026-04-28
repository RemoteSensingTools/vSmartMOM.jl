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
- quadrature_type: String in {RadauQuad(), GaussQuadHemisphere(), GaussQuadFullSphere()}
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
- aerosols: Vector of Dicts with keys: τ_ref, μ, σ, nᵣ, nᵢ, p₀, σp
- r_max: Real (µm)
- nquad_radius: Int
- λ_ref: Real (µm)
- decomp_type: String in {NAI2(), PCW()}
- n_ref: Complex (optional). If not provided, computed from first aerosol.

## Minimal YAML example

```yaml
radiative_transfer:
  spec_bands:
    - (1e7/777):0.015:(1e7/757)
  surface:
    - LambertianSurfaceScalar(0.15)
  quadrature_type: GaussQuadFullSphere()
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
