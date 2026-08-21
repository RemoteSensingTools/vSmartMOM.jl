# RRS–XCO2 work summary

Last updated: 2026-08-19  
Working branch: `suniti_multi_sensor`

## Objective

Build and validate a three-band OCO-2 forward/linearized simulation framework
covering the O2 A band and the weak and strong CO2 bands, including aerosols,
surface spectral reflectance, solar-induced fluorescence (SIF), rotational
Raman scattering (RRS), and a 64-scene truth map at fine spectral resolution.

## Three-band basis configuration

- Spectral calculations use a wavenumber spacing of `0.1 cm^-1`.
- Bands are O2 A, weak CO2, and strong CO2.
- The atmospheric representation was developed first with 12 layers and then
  refined to 16 layers for a better UTLS aerosol representation.
- Geometry for the truth-map calculations is nadir viewing with SZA = 30°.
- Lambertian surface spectra are represented independently in each band by the
  first three Legendre polynomials in wavelength/wavenumber coordinates.
- Forward and analytic-linearized workflows, timing scripts, output tables,
  and plotting scripts are under `RRS_XCO2/scripts/`.
- Basis plots and serialized results are under `RRS_XCO2/output/`.

The OCO-2 slit-function convolution and resampling stage remains downstream
work. The high-resolution simulations and band wavelength grids have been
saved in a form intended for that step.

## State-vector and Jacobian organization

The retrieval-state ordering is:

1. Surface parameters
2. Aerosol parameters
3. SIF parameters
4. XCO2 / gas parameters

Gas derivatives are shared physical parameters across bands rather than
separate band-specific state elements. For example,

```text
dI/dCO2(z) = [dI_A/dCO2(z), dI_weak/dCO2(z), dI_strong/dCO2(z)].
```

CO2 is absent from the A-band absorber configuration, so its A-band Jacobian
entries are naturally zero while the weak and strong-band entries contribute
to the same CO2 state element. SIF is active only in the O2 A band; its weak
and strong CO2-band derivatives are zero.

Every Jacobian row is plotted separately with parameter-aware filenames and
titles. Three vertical panels correspond to the three bands.

## SIF treatment

- The SIF spectral shape was taken from the source used by the historical
  `creategrid*.jl` workflow.
- The state variables are `SIF760`, referenced at 760 nm, and the local slope
  with respect to wavenumber.
- A total wavelength-integrated SIF of `0.5 mW nm^-1 sr^-1 m^-2` is converted
  consistently into wavenumber space before determining `SIF760` and its
  slope.
- Analytic SIF Jacobians were tested against finite differences.
- SIF is included only in the O2 A band for this experiment.

Relevant scripts include:

- `scripts/plot_creategrid_sif_shape.py`
- `scripts/validate_sif_jacobian.jl`
- `scripts/run_forward.jl`
- `scripts/run_linearized.jl`

## Surface-albedo library

ECOSTRESS spectra under `~/data/ECOSTRESSlib/` were used to create four
Lambertian surface cases:

- Urban: 20% trees, 50% roofing, 30% concrete paving
- Rural: 50% grassland, 30% roofing, 20% trees
- Desert: 100% Sahara dust surface
- Forest: 100% green vegetation

For each surface, the full spectrum is stored in NetCDF and the three
band-specific Legendre coefficients are stored in tabular form. Provenance,
source spectra, mixture recipes, and fitting details are documented in
`surface_albedos/PROVENANCE.md`.

Key products:

- `surface_albedos/*_surface_albedo.nc`
- `surface_albedos/lambertian_legendre_inputs.dat`
- `surface_albedos/surface_legendre_fits.png`

## Aerosol scenario

Three spherical, lognormally distributed aerosol species are used:

| Species | Median radius | Geometric sigma | Refractive index | AOD at 760 nm |
|---|---:|---:|---:|---:|
| Tropospheric sulfate | 0.10 micron | 1.60 | `1.43 + 1e-8i` | 0.1935471100 |
| Organic carbon | 0.20 micron | 1.80 | `1.53 + 1e-2i` | 0.0807084200 |
| UTLS sulfate | 0.15 micron | 1.50 | `1.45 + 1e-8i` | 0.0057444777 |

The total AOD is 0.28 at 760 nm. The component ratio is inherited from the
requested 8:1.8:0.2 ratio and normalized at 760 nm, not 550 nm.

Vertical profiles are normalized lognormal functions in geometric altitude,
integrated exactly between layer boundaries using CDF differences. Adopted
modes are approximately 1.2 km, 1.8 km, and 12 km. Tropospheric widths were
narrowed so their continuous tails are negligible by 10 km. The UTLS mode is
an intentionally enhanced sensitivity case, not a background Junge layer.

The literature reasoning, composition, size distributions, vertical-profile
parameters, and plots are documented in `truth_map_aerosols/README.md`.

## Truth-map simulations

The 64 scenes span:

- Uniform CO2 VMR: 380, 400, 420, or 440 ppm
- SIF: absent or total 0.5 mW nm^-1 sr^-1 m^-2
- Aerosol: absent or the three-component 0.28-AOD case
- Surface: urban, rural, desert, or forest

This gives `4 × 2 × 2 × 4 = 64` states. Exact state vectors are stored in
`truth_map/true_states.dat`; the more readable component definitions are in
`truth_map/scene_components.dat`.

Each `truth_map/hiressim_NNN.nc` contains the requested high-resolution
components:

- O2 A band: Cabannes, Rayleigh, and RRS-related output components
- Weak and strong CO2 bands: elastic/Rayleigh calculation only

The common spectral grids are stored in `truth_map/sim_wavelength.nc`.
Corrected 16-layer simulations are the active top-level dataset; the earlier
12-layer set was retained separately for comparison.

## Vertical-grid correction

The original layer-center interpretation was checked against the actual layer
boundaries and thicknesses. Aerosol optical depth is now integrated over
boundaries rather than assigned from a profile sampled only at layer centers.
Both 12- and 16-layer grids were validated and plotted side by side with:

- Exact horizontal layer boundaries
- Continuous aerosol profiles
- Per-layer integrated aerosol optical depths

Relevant products include:

- `truth_map/aerosol_vertical_distribution_12_vs_16_layers.png`
- `truth_map_aerosols/continuous_vertical_profiles.png`
- `truth_map_aerosols/layer_integrated_profiles_12_vs_16.png`

## External-solar `Z0`, `R0`, and `T0` design

The direct solar direction is outside the diffuse angular operator by default.
Instead of inserting SZA as a zero-weight quadrature node, the code evaluates
rectangular phase columns directly:

```text
Z0: NquadN × Npol × nSpec
R0/T0: NquadN × Npol × nSpec
Rdot0/Tdot0: NquadN × Npol × nSpec × nParams
```

Raman direct-source operators use the analogous representation with an
additional Raman dimension. This design is implemented for elastic forward,
elastic linearized, and forward rotational-Raman SFI paths. The legacy
embedded-solar representation remains available through
`external_solar=false` for unsupported paths and diagnostic callbacks.

Validation against the older square-operator implementation found:

- Rayleigh direct `Z0` columns: machine-precision agreement
- Direct tangent `Zdot0` columns: machine-precision agreement
- Elastic forward SFI: 70/70 regression checks passed
- Rotational Raman: 9/9 checks passed
- Spectrally varying aerosol `Z0`: maximum absolute difference exactly 0
- Aerosol forward SFI: maximum absolute and relative differences exactly 0
- Aerosol linearized TOA maximum relative difference: `1.96e-16`
- Aerosol Jacobian maximum relative difference: `3.39e-11`

For five weighted streams, one distinct VZA, and IQU polarization, the diffuse
operator was reduced from `21 × 21` to `18 × 18`. The solar direction is no
longer counted as a diffuse stream.

## Spectral optical-property optimization

The forward and linearized optical-property construction was reorganized so
that quantities independent of Fourier order `m` are cached outside the
`m` loop. Rayleigh phase blocks depend on `m` but not wavelength and therefore
are not recomputed for every spectral point. Aerosol phase matrices remain
inside the `m` loop because their angular structure is strongly wavelength
dependent.

Aerosol `Z`, optical depth, single-scattering albedo, and their derivatives are
computed only at band endpoint/reference nodes and interpolated across the
band:

- Two nodes: linear interpolation
- Endpoints plus an interior reference wavelength: natural cubic spline
- Phase interpolation occurs after BGE truncation

The same applicability was checked for the inelastic path.

## Delta-BGE polarization normalization

A normalization asymmetry was identified in `truncate_phase.jl` and its
linearized counterpart: beta-derived families were normalized by the retained
fraction `c0 = 1-ft`, while fitted gamma and epsilon were not. This diluted
polarization for coarse particles while leaving intensity largely unchanged.
The implementation and regression coverage were updated so all relevant Greek
coefficient families use the consistent retained-fraction normalization.

## Single-scattering interaction fix

The separate `interaction_ss!` issue discovered during diagnostics was fixed
and tested. It is not on the multiple-scattering path used for the truth map,
but resolving it prevents the known defect from remaining latent.

## Fourier-order and stream investigation

The original untruncated Mie allocation at 757.001655 nm used a conservative
maximum size parameter and produced `m_max = 220`, requiring 111 streams under
`stream_l_cap = 2*nstreams - 1`.

This overweights extremely weak tails of the particle-size distributions. A
new optional YAML parameter was implemented:

```yaml
radiative_transfer:
  greek_beta_cutoff: 1.0e-5
```

The evaluation is explicitly beta-only:

1. Maximum size parameter establishes the allocated Mie-series ceiling.
2. After integration over the full size distribution, retain through the last
   degree satisfying `abs(beta_l) >= greek_beta_cutoff` at any wavelength in
   the band.
3. Evaluate each aerosol separately and take their maximum.
4. Combine with Rayleigh, surface, and source traits.
5. Apply explicit user and stream caps.

Other Greek coefficients are not used because they may legitimately be zero
or negative and are not the scalar phase-function support benchmark.

At the shortest O2 A-band wavelength, the retained degrees are:

- Tropospheric sulfate: 24
- Organic carbon: 99
- UTLS sulfate: 22
- Combined aerosol effective support: `m_max = 99`

Thus an uncapped calculation needs 50 rather than 111 streams. The cutoff
shortens the Fourier loop; normal production configuration still treats
`nstreams` as an independent angular-resolution choice.

Subsequent Float64 diagnostics showed that radiances are already converged by
approximately `l_trunc=16` for this aerosol mixture. Production aerosol scenes
therefore use explicit `δBGE(16, 0°)` with 9 weighted streams
(`stream_l_cap=17`). Aerosol-free scenes retain 5 streams and no truncation.
`Nquad` is derived from weighted streams plus requested viewing nodes and is
not set directly.

The parser, JSON schema, default YAML, code docstrings, tests, user docs,
developer notes, `AGENTS.md`, and `CLAUDE.md` document this rule. Focused
component and schema tests passed.

## Untruncated moment-resolved datasets

At 757.001655 nm, isolated contributions were saved for every Fourier order
`m = 0:99`, every VZA from 0° through 70° in 5° increments, and relative
azimuths 0° and 180°. Negative signed VZA denotes relative azimuth 180°.
The files also contain cumulative sums over `m`.

| Dataset | Status | Total recorded compute time |
|---|---|---:|
| `truth_map_aerosols/untruncated_stokes_by_m_nstreams50_mmax99.jld2` | Complete, 15/15 VZAs | 2.899 h |
| `truth_map_aerosols/untruncated_stokes_by_m_nstreams111_mmax99.jld2` | Complete, 15/15 VZAs | 9.419 h |

Both datasets above are Float32 diagnostic artifacts and must not be used as
untruncated references. A replacement
`untruncated_stokes_by_m_float64_nstreams50_mmax99.jld2` was started on
2026-08-19 with the same VZA, azimuth, and per-m layout.

ASCII companions with the isolated and cumulative I/Q/U values are stored
beside both JLD2 files. These diagnostic runs use the embedded-solar layout
because the per-m callback is intentionally unavailable from the external-
solar TOA entry point; validated radiances are equivalent.

The 50-stream comparison against delta-BGE truncation levels 8, 12, 16, 24,
and 32 is plotted in:

- `truth_map_aerosols/truncation_vs_untruncated_50stream_IQ.png`

The left column contains I and Q. The right column contains percentage
difference for I and absolute difference for Q. U was omitted because it is
zero in the principal-plane geometries used here.

## Principal diagnostic finding

> **⚠ SUPERSEDED (2026-08-19).** The untruncated references described in this
> and the preceding two sections are numerically defective, not references.
> The discrepancy is Float32 precision amplified by stream count; the Legendre
> truncation and the phase function play no role. See
> [`truth_map_aerosols/HANDOFF_untruncated_artifact.md`](truth_map_aerosols/HANDOFF_untruncated_artifact.md).
> The text below records the earlier interpretation.


The earlier discontinuous untruncated angular result was localized largely to
the `m=0` contribution rather than the sum of nonzero Fourier moments. Phase
functions themselves were smooth, and individual beta coefficients decayed
smoothly. This motivated the detailed per-m, per-stream diagnostics and the
new beta-tail support criterion. The new 50- and 111-stream moment-resolved
datasets now permit a clean quadrature-resolution comparison at identical
`m=0:99` support.

## Important files

- `config/oco_grass_3aerosol.yaml`: principal three-band configuration
- `scripts/generate_truth_states.jl`: builds the 64-state table
- `scripts/generate_truth_map.jl`: generates high-resolution truth spectra
- `scripts/diagnose_untruncated_stokes_by_m.jl`: resumable per-m/VZA benchmark
- `scripts/plot_aerosol_truncation_vs_untruncated_50stream.py`: latest I/Q plot
- `surface_albedos/PROVENANCE.md`: surface-data provenance and recipes
- `truth_map/README.md`: truth-map schema and simulation contents
- `truth_map_aerosols/README.md`: aerosol provenance and parameter rationale

## Remaining work

- ~~Compare the completed 50- and 111-stream moment-resolved datasets directly,
  especially `m=0`.~~ Done 2026-08-19: the disagreement is entirely in `m=0`
  and is a Float32 defect, not a quadrature-resolution question. See
  [`truth_map_aerosols/HANDOFF_untruncated_artifact.md`](truth_map_aerosols/HANDOFF_untruncated_artifact.md).
- Fix the three numerics bugs identified there (`_convert_numerics` casting
  precision-dependent defaults; `rt_kernel_ss.jl` bypassing `get_dtau_ndoubl`;
  `gfct` inconsistency in the directional `get_dtau_ndoubl`), and add a
  Float32/Float64 stream-sweep regression test.
- Complete and plot the Float64 50-stream, `m=0:99` untruncated replacement.
- Complete OCO-2 slit-function convolution and resampling once the weak and
  strong CO2-band calibration/slit resources are available.
- Re-run the final production truth map if the stream/truncation comparison
  changes the selected production angular resolution.
