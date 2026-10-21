# Changelog

## 2.2.0 — 2026-10-21

Scientific authorship of the linearization and Raman work, including the
formulations, analytic derivatives, and retrieval methodology, is credited to
**Suniti Sanghavi**. Christian Frankenberg assists with performance optimization,
backend engineering, integration, packaging, and release tooling. See the
release manual's contribution statement
and the published references in `CITATION.bib`.

### Consolidated release

This release combines the surface/multisensor integration, Suniti's updated
retrieval workflows, and the validated Jacobian improvements. It preserves the
published `v2.1.0` tag. The migration baseline for registry users is **1.1.0**;
2.0/2.1 were not registered in General.

- Named endpoint/interior radiance results retain tuple-style iteration.
- Retrieval-selected, local-optical, and equivalent-source analytic Jacobians
  support efficient elastic inversions, including surface SIF.
- Historical Sanghavi inversion logs versus optimized replays show about
  **12–20×** lower retrieval wall time in the tested cases. These are workload
  measurements, not a universal Jacobian/forward cost ratio. See
  `docs/dev_notes/jacobian_batched/suniti_inversions.md` for timing and precision
  conditions, and the retained convergence/convolution accuracy evidence.
- Cox–Munk glint now uses the incident spectrum, both atmospheric transits,
  and consistent Fourier/BRDF normalization. Forward and analytic Jacobian
  runs apply the same correction, including wind and optical-depth tangents.
  **Ocean radiances change**; revalidate prior ocean baselines. The direct
  glint correction supports unpolarized incident light.
- Cox–Munk slope statistics corrected: `0.003 + 0.00512U` is now used as the
  *total* mean-square slope (per-axis variance half of it; previously the
  glint was √2 wider — twice the slope variance — and half as bright at the
  specular point), and the
  Smith shadowing function uses its correct `√π` normalization (shadowing at
  grazing angles was underestimated or zero). Both change ocean radiances.
- Polarized Cox–Munk reflection corrected: Stokes vectors are transformed
  explicitly between each beam's meridian plane and the facet's plane of
  reflection; all callers pass (incident, reflected) as documented; the
  surface Fourier moments follow the atmospheric phase-matrix storage
  convention (the I,Q ← U,V block had the opposite sign, which affected
  polarized skylight reflected at the surface); and the circular (V) terms of
  the Fresnel Mueller matrix follow vSmartMOM's convention. The full Mueller
  matrix now agrees with VLIDORT 2.8.3's Cox–Munk kernels. **Polarized ocean
  runs change in every Stokes component, including I when there is a
  scattering atmosphere** (up to ~3e-3 relative in a Rayleigh test case);
  `Stokes_I` runs and vacuum intensity are unchanged by this item.
- Cox–Munk water refractive index: `n_water = nothing` now uses the
  Segelstein (1981) pure-water index (25 °C) at every spectral point, with n
  and k linear in wavelength (table in `src/CoreRT/Surfaces/data/`, from the
  refractiveindex.info database, CC0). It used to take n = 1.328 at 550 nm
  from a built-in table that departs from Segelstein and Hale & Querry above
  about 900 nm. A vector `n_water` is used element-wise and must have one
  entry per spectral point (it used to use its first element); a scalar is
  unchanged bit for bit. The per-point reflection operator is computed
  exactly at a few knot indices and interpolated in n (within 1e-7 of the
  reflectance in Float64 for water at 17 bands from 0.25 to 3.9 µm, and
  within 1e-6 in Float32 at 0.55, 1.6, 2.06 and 2.8 µm), and the TOA exact
  glint uses each point's own index. Single-wavelength helpers
  (`coxmunk_brdf_mueller`, `reflectance`, …) now throw when the surface's
  `n_water` is `nothing` or a vector with differing entries (they fell back
  to 550 nm or to the first element), and `water_refractive_index` throws
  outside 34 nm–10 m instead of clamping. **Ocean radiances change**: TOA I
  in the glint view (U = 5 m/s, SZA 30°) changes by +3.6 % at 0.55 µm,
  −0.5 % at 0.76 µm, −9.3 % at 1.6 µm and −16.9 % at 2.06 µm. A standalone
  `CoxMunkSSSurface(n_water = nothing)` built without a model keeps its
  documented 1.34.
- **Ross–Li results change against v2.1.0**: the kernels cap both zenith
  angles at 80° (`e2c1ad5b`; v2.1.0 had no cap), which changes I by up to 0.8
  % at TOA and 1.1 % at BOA in the VLIDORT Case D comparison.
  `RossLiSurfaceScalar(fvol, fgeo, fiso; zenith_cap_deg = 80)` makes the cap a
  parameter: `zenith_cap_deg = nothing` (`nothing` or `null` in YAML/TOML
  surface strings) gives the uncapped MODIS/VLIDORT kernels, which match
  uncapped VLIDORT 2.8.3 (Case D) to ≤ 7e-5 relative in I. The default and the
  three-argument form are unchanged against the v2.2 pre-releases. Ross–Li
  surface strings that combine positional and keyword arguments now keep the
  positional values (they were dropped).
- Documentation: the exact Cox–Munk glint correction is applied at TOA only;
  interior observer levels keep the Fourier-truncated glint. For a
  canopy, `bhr_dw` is diffuse only, and the forward BOA radiance at the solar
  ordinate has no direct beam (the canopy sets j₀⁺ = 0; the linearized solver
  adds it).
- Mie radius quadrature fixed per run: the nodes lie on `[r_min, r_max]`, with
  a new optional `scattering.r_min` (µm; default `Scattering.DEFAULT_R_MIN` =
  1e-3 µm), instead of starting at each size distribution's 1e-8 quantile.
  That quantile moved the nodes with μ_logr/σ_logr while the analytic size
  tangents assume fixed nodes. The μ_logr/σ_logr Jacobians now equal finite
  differences of the forward model at any `nquad_radius`. Before, max |analytic
  − FD| / max |FD| of the TOA Jacobian was up to 18 % for a 0.2 µm,
  σ_g = 1.6 mode on JacobianTestFast settings and 26 % for the organic mode of
  an OCO-like three-mode scene at 300 nodes (O₂ A band). Forward aerosol
  optics change by the radius-quadrature error at the configured node count:
  TOA I by ≤ 1e-5 relative in that OCO-like scene and by 3.4e-4 in a two-mode
  test scene at 100 nodes. NAI2 (forward, linearized, GPU) and PCW use the
  same nodes. At a fixed `nquad_radius` the nodes span a wider range of ln r
  (1.2–1.6× for the tested modes), so check convergence, in particular for
  coarse modes. `MieModel` has an `r_min` field, `make_mie_model` an `r_min`
  keyword and `ScatteringParameters` a trailing `r_min`; the old constructor
  forms use the default.
- The standalone line-by-line model (`make_hitran_model`, `HitranModel`) is
  deprecated (removal in 3.0) with a one-time notice giving the
  AtmosphericAbsorption calls and their convention differences. The LUT
  builder `make_interpolation_model`, still the route for YAML `LUTfiles`, now
  uses its `vmr` argument (previously ignored) and rejects wavelength grids
  early with a uniform-wavenumber recipe.
- Combined-band solves (`rt_run(RS_type, model, iBand)` with several bands,
  and the single-scatter and multisensor equivalents) apply one lower
  boundary. They now require every selected band to hold the same surface
  (an identical plain-value surface counts); a shared spectral or spline
  albedo must cover all selected bands' spectral points, and Legendre albedos
  and canopies are rejected across bands. Previously the first band's surface
  was used silently.
- Pure-gas molecular-Rayleigh API in `Scattering` (also exported from
  `vSmartMOM`): `molecular_rayleigh_properties` and
  `molecular_rayleigh_cross_section_ratio` dispatch on
  `AbstractMolecularRayleighMethod`. `UsingRefInd()` is the refractive-index
  plus King-factor method for He, Ar, N₂, O₂ and CO₂ (`RefractiveIndexRayleigh`
  is kept as an alias, and the forms without a method argument use it);
  `UsingPolTen()` uses the polarizability tensors of N₂, O₂ and H₂.
  `molecular_rayleigh_provenance` covers `UsingRefInd()` for Ar and CO₂. The
  refractive-index part is ported from Suniti Sanghavi's Habitable Worlds
  Observatory branch. The Rayleigh optical depth is still Bodhaine (1999),
  rescaled by the King factor of the depolarization in use.
- N₂ and O₂ refractive-index fits extend through the NIR: Peck & Khanna
  (1966) for N₂ beyond 436 nm and Egan (2026) for O₂ beyond 467 nm, joined to
  the Bates (1984) fits where the two agree, so σ is continuous. σ(N₂) rises
  by up to 0.21 % and σ(O₂) falls by up to 0.89 % at 2 µm; air changes by
  less than 0.03 %. The O₂ King factor stays 1.09 (He et al. 2021). Sources,
  joins and numbers: `docs/dev_notes/molecular_rayleigh_dispersion.md`.
- **Results change for every `depol: -1` model** built by
  `model_from_parameters(params)`, `model_from_parameters(LinMode(), params)`
  or `update_model!`. The VS `_plus` builder follows the same selection; its
  `depol: -1` handling is new (see the VS `_plus` builder item below). The
  N₂/O₂ polarizability anisotropy now varies with wavelength
  (`DispersiveAnisotropy()`, implied by the refractive-index fits and King
  factors). It is the default of `UsingPolTen()` and of the automatic
  `depol: -1` Rayleigh depolarization, which sets the Rayleigh Greek
  coefficients and τ_Rayleigh. For air the depolarization is 0.0277 at
  400 nm and 0.0254 at 1.6 µm (v2.2 pre-releases 0.0263 and 0.0280; v2.1.0
  about 0.0281 at both), and τ_Rayleigh changes by +0.25 % at 400 nm and −0.44 %
  at 1.6 µm against the pre-releases (−0.06 % and −0.46 % against v2.1.0).
  `radiative_transfer.rayleigh_anisotropy: ConstantAnisotropy()` (params
  field `rayleigh_anisotropy`) keeps the constant anisotropy, the v2.2
  pre-release depolarization bit for bit, and
  `UsingPolTen(ConstantAnisotropy())` the pre-release `UsingPolTen()` values;
  neither restores v2.1.0 (see the next item), and the new `UsingRefInd()`
  fits have no switch. The
  Raman and Cabannes quantities (ϖ_Cabannes, the Cabannes depolarization,
  RRS and VS cross sections) keep the constant Raman-module tensor of the
  published tables. Paired noRS/Raman runs (RRS or VS; Ring or filling-in
  spectra) should set `ConstantAnisotropy()` until the Raman side is decided.
- Molecular constants of `depol: -1` and of the Cabannes fraction (Suniti
  Sanghavi, `817d93cb`): the N₂/O₂ constants are evaluated at the
  dry-column-weighted mean temperature of the profile instead of 300 K, and
  the Buldakov et al. dispersion pole of the mean polarizability uses the
  angular frequency 2πcν̃ (v2.1.0 used cν̃). The 2π dominates. With the
  constant anisotropy, the air depolarization becomes 0.0263 at 400 nm,
  0.0276 at 760 nm and 0.0280 at 1.6 µm (v2.1.0: about 0.0281 from 400 nm
  to 2 µm), τ_Rayleigh changes by −0.31 % at 400 nm and −0.09 % at 760 nm,
  and the RRS share 1 − ϖ_Cabannes falls by about 6 % at 400 nm and 2 % at
  760 nm.
- **Raman (RRS and VS `_plus`) results change in absorbing bands** (Suniti
  Sanghavi, `817d93cb`): the Rayleigh fraction that weights the inelastic
  source (`fscattRayl`) is the Rayleigh share of the total layer extinction,
  τ_Rayl/(τ_Rayl + τ_aer + τ_abs), where τ_aer is the δ-M-scaled aerosol
  extinction ((1 − fᵗϖ) times the unscaled one). The RRS and VS `_plus` elemental kernels
  attenuate with the total dτ and do not multiply by ϖ, so the fraction must
  be taken over the total extinction. v2.1.0 left gas absorption out of the
  denominator: at incident wavelengths inside absorption lines each layer's
  Raman source was too large by (τ_Rayl + τ_aer + τ_abs)/(τ_Rayl + τ_aer).
  noRS runs, and Raman runs without gas absorption, are unchanged.
- **Raman runs with `depol: -1` change slightly**: the Cabannes depolarization
  inverts Z_Rayl = ϖ_Cab·Z_Cab + (1 − ϖ_Cab)·Z_RRS (γ_RRS = 3/4) exactly
  (Suniti Sanghavi, `817d93cb`). v2.1.0's `compute_γ_mol_Cabannes!` used the
  denominator (1 + 2γ_R)(1 − ϖ) + 5 instead of (3 + 2ϖ) − 4γ_R(1 − ϖ), which
  agrees only as ϖ_Cab → 1. With the current molecular constants (2π pole) at
  300 K, the Cabannes γ rises by 1.5 % (N₂), 4.2 % (O₂) and 2.5 % (air) at 395
  nm, and by 1.6 %, 4.5 % and 2.7 % at 760 nm. The air Cabannes depolarization
  2γ/(1+γ) at 760 nm becomes 0.00737 (v2.1.0 formula, same constants:
  0.00718), and the factor (1−γ)/(1+2γ) of its Greek coefficients changes by
  −2.8e-4. It enters only the Cabannes Greek coefficients of RRS and VS runs;
  noRS runs do not use it.
- Raman: each per-J Q-branch line (ΔJ = 0) weights the anisotropy by the
  Placzek–Teller coefficient b_JJ instead of b_JJ², so Q + S + O lines sum to
  the Rayleigh (or vibrational-band) total for every J and the summed Cabannes
  lines agree with Sanghavi (2026, Table 1) to within 1.5 %. The vibrational Raman
  (VS) Q-branch cross section rises by 13.3 % (O₂) and 5.9 % (N₂), which
  changes VS / VS_plus runs, their Cabannes fractions and the stellar H₂ VS
  path; RRS, Rayleigh and all noRS/RRS results are unchanged. Derivation and
  references:
  `docs/dev_notes/raman_q_branch_placzek_teller.md`.
- Raman (RRS): `getRamanSSProp!` throws an `ArgumentError` for a band that
  holds no N₂/O₂ rotational line (the smallest shifts are ±11.9 cm⁻¹ for N₂
  and ±14.3 cm⁻¹ for O₂). Such a band used to get `n_Raman = 0` and lose the
  Raman share 1 − ϖ_Cabannes (3.2 % at 395 nm) without notice. RRS weights
  that sum to zero (an RRS object never filled) are rejected instead of
  giving NaN.
- Raman (RRS): each rotational line is shared between its two neighbouring
  grid points by linear interpolation, with shifts measured from the grid
  point nearest the reference wavelength and the source offsets explicitly
  negated; a reference wavelength outside the band grid is an
  `ArgumentError`. Results no longer depend on the parity of the band grid.
  Elastic results are unchanged. In the Phase 1b band (762–765 nm, 103
  points, flat F₀) the inelastic radiance changes by ≤ 2.1e-6 (relative)
  where every source lies inside the band, and by up to 12 % (≤ 0.17 % of I)
  at the six points where a line's source crosses the band edge. With a
  structured F₀ (60 %-deep features every 3.7 cm⁻¹) points whose sources all lie
  inside the band change by up to 3.9 % of the inelastic part on that
  odd-length grid (0.20 % of I) and up to 12 % on an even-length one (0.38 %
  of I). The VS, VS `_plus` and stellar H₂ paths share the line-splitting
  routine, so their sub-grid placement changes in the same way.
- Removed `getRamanSSProp!(::RRS, depol, λ, ν)`, which swapped Stokes and
  anti-Stokes (use `getRamanSSProp!(rs, λ, ν)`), and
  `getRamanSSProp!(::Union{VS_0to1, VS_1to0}, …)`, which could not run.
  Nothing in vSmartMOM, uni_vSmartMOM or ExoOptics called either.
- `VS_0to1_plus` depletes its source band (band 1, the incident wavelength)
  by f = (σ_VRS + σ_RVRS)/σ_Rayleigh (8.7e-4 at 395 nm) and no longer the
  target bands (it had ϖ_Cabannes = [1, 1 − f, 1 − f]), as `VS_1to0_plus`
  already did. Band-1 I changes by −1.3e-3 to −1.5e-3 (relative) at 395 nm
  over four viewing directions.
- The VS Q-branch phase matrix uses the population-weighted ⟨b_JJ⟩ of its
  lines instead of the whole-band depolarization: ρ goes 0.157 → 0.046 (N₂)
  and 0.290 → 0.100 (O₂) at 300 K. The band-summed inelastic I changes by
  −4.1 % to +4.1 % over four viewing directions (`VS_0to1_plus`, 395 nm,
  SZA 45°).
- VS `_plus` model builder: the Rayleigh and Cabannes depolarizations are
  resolved as in the standard builder. With `depol: -1` it used to pass −1
  itself, which gave τ_Rayleigh 0.22× the correct value and β₂ = 2 in
  `greek_rayleigh` (negative radiances from noRS on a VS model); band-1 I at
  nadir goes 7.59e-3 → 3.34e-2 (395 nm, pure Rayleigh, SZA 45°, constant
  anisotropy). With an
  explicit `depol` ≥ 0 the Cabannes Greek coefficients now use that value
  too (before, they came from the molecular constants for any `depol`), so
  those VS runs change as well. `params` is no longer modified, and line-by-line
  absorption and Float32 builds work. The VS `_plus` forward run works again
  (its kernel and interaction had fallen behind workspace changes; no physics
  change).
- Documentation: the Greek matrix B_l is shown as printed by Sanghavi &
  Stephens (2015, Eq. 22) in Siewert's form, which the vSmartMOM paper uses by
  reference (−ε at [3,4]), with the exact mapping to the code's stored γ and ε;
  no numerical change.
- Documentation: the surface-pressure Jacobian page states that, under
  `FixedPressureGrid()`, the temperature dependence of the automatic
  (`depol: -1`) Rayleigh depolarization (Greek coefficients and τ_Rayleigh)
  has no tangent (about 1e-5 relative, neglected by decision).
- Tests: `test/vlidort_baseline/extract_references.jl` now emits the
  committed `siewert2000_IIA_greek.jl` byte for byte, including the γ import
  flip it previously left out.
- HITRAN cache: `artifact(molecule)` under a downloaded edition now uses the
  cached file only if it is a verified full dataset (0–150000 cm⁻¹, all
  isotopologues, SHA-256 matching its metadata); otherwise it downloads the
  full dataset and, if that fails (e.g. offline), raises an error explaining
  why instead of silently using a narrower or unverified file. Downloads are
  validated before an atomic replacement that keeps the previous data on
  failure. `fetch_hitran`/`fetch_hitran_by_ids` reuse a cache only if its
  molecule, edition, spectral range, isotopologues, size and SHA-256 match
  (missing metadata forces a re-download). `hitran_is_cached(molecule[,
  edition])` now means "`artifact()` would use a verified full dataset without
  downloading".
- The positional linearized form `rt_run(RS_type, model, lin_model, NAer,
  NGas, NSurf, iBand)` accepts a one-element band vector and rejects several
  bands with a clear message (it solves one band per call).
- Polarized Jacobians of layers that are not doubled (scattering optical
  depth below `max(dτ_min_floor, dτ_max_threshold·μ_min)`) are corrected: the
  phase-matrix term of the reverse operators r⁺⁻ and t⁻⁻ had the wrong sign in
  the off-diagonal (I,Q)↔(U,V) blocks, so aerosol, size/profile and
  surface-pressure columns of `Stokes_IQU`/`Stokes_IQUV` runs with such layers
  were wrong (7.5% in the affected layers' columns in a 34-layer test); their
  gas columns are unchanged. The same kernel also treated Q as a U/V row for
  `Stokes_IQ`, which gave wrong I and Q Jacobians, gas columns included, in
  both doubled and undoubled layers. `Stokes_I` Jacobians and all forward
  radiances are unchanged.
- `model_from_parameters(LinMode(), params)` no longer adds water-vapour line
  absorption when the configuration has no `absorption` block; as documented
  and as the forward constructor does, omitting the block means no gas
  absorption (`q` still sets the dry-air column). Previously the linearized
  model described a different atmosphere from the forward model.
- `DefaultParameters.yaml` gives the aerosol pressure peak and width in hPa
  (`p₀: 900`, `σp: 50`), the unit of the profile. They were given in Pa,
  which placed the aerosol far below the surface: an exponential tail toward
  the ground in Float64 and NaN optical depths in Float32. Aerosol results of
  `default_parameters()` and of GEOS-Chem-generated configurations
  (`geoschem_to_dict`, `GeosChemSource`), which start from these defaults,
  change.
- **Aerosol results change for altitude-form profiles** (`z₀`, `σ₀`; Suniti
  Sanghavi, `b5f6e425`): the log-normal height distribution is integrated
  over the geometric layer interfaces (CDF differences, normalized over the
  model column) in the forward and linearized constructors. v2.1.0's
  forward constructor evaluated the density at the layer pressures in hPa as
  if they were heights in km (`pdf(LogNormal(log z₀, σ₀), p_full)·Δz`),
  which put the aerosol in the topmost model layers: for
  `ParamsEMIT_MODTRANcomp.yaml` (z₀ = 2 km, σ₀ = 0.4, model top 20 hPa) the
  forward AOD-weighted mean height moves from 24.7 km, with 99.997 % of the
  AOD in the top layer (20–40 hPa), to 2.2 km.
- **Linearized aerosol profiles match the forward model** (Suniti Sanghavi,
  `ee21b857`): v2.1.0's linearized constructor rebuilt every profile as a
  pressure Gaussian from the distribution's mean and standard deviation and
  weighted it by Δp, while the forward constructor weighted the layer-centre
  density by Δz. Linearized radiances and profile Jacobians of pressure-form
  profiles (`p₀`, `σp`) therefore change and now equal the forward model's.
  For altitude-form profiles the log-normal's moments were used as `p₀` and
  `σp` in hPa (in the Float32 configuration above every weight underflowed
  to zero and the profile was NaN); the aerosol's two profile Jacobian
  columns are now ∂/∂z₀ and ∂/∂σ₀.
- **Results change for runs with `profile_reduction: n`** (Suniti Sanghavi,
  `4c0758e7`): T, q and the VMR profiles are interpolated in log pressure on
  the actual input layer pressures (flat extrapolation), and the H₂O ratio is
  recomputed from the interpolated q. v2.1.0 interpolated them as if the input
  layer centres were evenly spaced in pressure between the shallowest and the
  deepest one. Without interior observer heights the reduced pressure levels
  are the same as in v2.1.0; the reduced layers' geometric thicknesses now
  come from the full-resolution profile (interpolated in log p) instead of a
  hydrostatic recomputation with the reduced T. For the
  `DefaultParameters.yaml` profile reduced to 5, 12 or 16 layers, layer
  temperatures differ by up to 53–55 K (n = 12: 274.0 K in v2.1.0 against
  218.6 K now at 209.5 hPa). A target of at least the input layer count now
  leaves the profile unchanged instead of failing an assertion.
- `remake_geometry` keeps the resolved observer heights, the TOA/BOA output
  selection, and the external-solar representation when it changes the solar
  or viewing angles. Previously it dropped interior observers and the TOA
  output metadata, and returned external-solar models to the embedded
  representation.
- `external_solar` defaults to `false` (embedded μ₀: the solar direction is
  an extra quadrature node) in the forward and linearized
  `model_from_parameters`. Code written for cores that defaulted to `true`,
  such as the HWO/ExoOptics core, must pass `external_solar = true` when
  building the model before calling `rt_run_toa` on the external-solar path.
- **Radiances change for runs with a non-unit or spectral solar `F₀`** (Suniti
  Sanghavi, `817d93cb`): the direct solar beam at the surface is the solve's
  `F₀` attenuated by exp(−τ/μ₀) for every Lambertian flavour and the analytic
  BRDFs (RPV, Ross–Li, Cox–Munk), in `rt_run`, `rt_run_ss` and the Raman
  solvers. v2.1.0 applied `F₀` to atmospheric scattering but used a unit,
  I-only beam for the surface's direct-beam source terms, so the directly
  reflected surface term and the BOA direct beam were computed for F₀ = 1:
  with a spectrally flat F₀ = c they were a factor c too small relative to the
  rest of the radiance, and with a solar spectrum they had none of its
  Fraunhofer structure. (`LambertianSurfaceLegendre`/`Spline` used the unit
  beam for their reflected term but had no BOA direct beam, and
  `LambertianSurfaceSpectrum` failed; see the Lambertian scaffold item below.)
  Runs with the default unit beam are unchanged by this item. The v2.1.0
  linearized Lambertian and Cox–Munk layers already used `F₀`, so forward and
  linearized radiances now agree for non-unit F₀.
- `CanopySurface` leaf scattering now uses the solar Stokes vector `F₀` of
  the solve, as its soil now does (previous item). It used a unit beam, so
  canopy radiances were not proportional to a spectral solar irradiance or a
  `SolarBeam(F₀=…)` scaling (e.g. R(2F₀)/R(F₀) was 1.34–1.45 instead of 2,
  measured with the soil already using F₀), and a zero solar beam still gave
  leaf-scattered radiance. Runs with the default unit beam are unchanged.
- **BOA downwelling radiance and the downward BHR change for
  `LambertianSurfaceLegendre` and `LambertianSurfaceSpline`** (shared
  Lambertian scaffold, `e6af91f1`; `t⁺⁺ = I` also Suniti Sanghavi,
  `817d93cb`): for m > 0 the surface layer passes the downwelling field
  (`t⁺⁺ = I`; v2.1.0 used 0, so BOA output kept only the m = 0 azimuthal
  mean), and `j₀⁺` carries the attenuated direct beam (v2.1.0 used 0, so BOA
  output and the downward BHR omitted it), as for `LambertianSurfaceScalar`.
  TOA radiances are unchanged. A single-coefficient Legendre albedo no longer
  fails.
- StandaloneSS `truncated_ss_path2` no longer halves the m = 0 Cox–Munk
  Fourier term (the azimuthal mean), matching the CoreRT correction to the
  glint residual; the truncated series now converges to `run_exact_ss`, and
  Cox–Munk `apply_back_correction!` results change accordingly. Its azimuth
  rule grows with the requested number of moments (at least 100 nodes,
  2·max_m above 50 moments).
- `truncate_phase(δBGE, …)` (and its linearized form) no longer fails for
  phase matrices with an identically zero f₁₂ or f₃₄, such as the analytic
  Henyey–Greenstein and synthetic polarized Henyey–Greenstein phase
  functions; those elements truncate to zero. The linearized form rejects a
  nonzero tangent of such an element, for which the fit has no derivative.
- **Polarized results with δBGE-truncated aerosols change** (`85e46b94`):
  the truncated γ and ϵ Greek coefficients are divided by c₀ = 1 − fᵗ, as β,
  α, δ and ζ already were (Sanghavi & Stephens 2015, Eq. 8), in the forward
  and linearized truncation. v2.1.0 omitted the division, so f₁₂/f₁₁ and
  f₃₄/f₁₁ of the truncated aerosol phase matrix were too small by the factor
  1 − fᵗ (29 % at fᵗ = 0.29). `Stokes_I` runs and untruncated aerosols are
  unchanged.
- `update_aerosol_microphysics!` (BatchContext) uses the constructor's
  anchor wavelengths, reference-wavelength node rule and interpolation
  endpoints for multi-wavelength bands. In Float32 it previously inserted an
  extra knot one ulp from the band edge, so even an identical update changed
  k, fᵗ and the radiance (~1e-4) instead of reproducing a fresh build.
- The parameter parser warns (once per key and session) about keys it does
  not read at the top level and in the `radiative_transfer` and `geometry`
  blocks, the sections the schema closes; a misspelled optional key such as
  `nstream` previously left the default in effect silently. Such keys remain
  non-fatal.
- Configuration and API changes:
  - New `radiative_transfer.greek_beta_cutoff` (Suniti Sanghavi,
    `b5f6e425`). Omitting it keeps the full Mie series, but configurations
    with aerosols that omit it now log a migration warning.
  - Setting both legacy `max_m` and `m_max` raises an error (`90d8788c`;
    v2.1.0 accepted both).
  - `compute_aerosol_optical_properties(model; autodiff=true)` is retired
    and raises an error (`a332d6d8`); use
    `compute_aerosol_optical_properties(LinMode(), …)`.
  - Raman solvers warn once when an explicit `depol ≥ 0` is set, because it
    is then used for Rayleigh, Cabannes and τ_Rayleigh alike (Suniti
    Sanghavi, `817d93cb`).
  - `build_sif_source(SIF₀::AbstractMatrix, ν_model, ν_sif, jSIF)` fills a
    matrix for `SurfaceSIF` (Suniti Sanghavi, `817d93cb`).
  - Parsed HITRAN line lists are cached per process (`ee0a8722`, `505dffcf`):
    after replacing a file at the same path by means other than
    `fetch_hitran`/`artifact`, call `clear_spectroscopy_cache!()`.
  - New forward `ThermalEmission` per-layer Planck source (Suniti Sanghavi,
    `3212b3f4`).
  - New `rt_run_streams` / `streams_callback` per-Fourier-moment
    quadrature-stream export (Suniti Sanghavi, `51307d14`, `d6051b7b`).
- **Breaking for positional Jacobian indexing** (Suniti Sanghavi,
  `5ba80c51`, `f1efc5ff`): the column order is `[p_surf, 7 per aerosol,
  gas columns per layer (species-major, TOA to BOA), surface, optional SIF,
  canopy]`; v2.1.0 had `[7 per aerosol, one column per gas, surface,
  canopy]`. The positional `NGas` of the linearized `rt_run`
  (`size(lin_model.τ̇_abs[1], 1)`) now counts one column per species and
  layer. Code that hard-codes offsets (e.g. `1:7` for the first aerosol) must
  use `ParameterLayout` (`psurf_index`, `aerosol_range`, `gas_profile_range`,
  `surface_range`) or `result.layout`.
- **Breaking for some linearized callers:** gas-abundance Jacobian columns
  are the native partial ∂τ/∂VMR = σ·N_dry
  with cross sections, CIA and MT_CKD held fixed (except the ABSCO
  H₂O-broadener term below). Solves and Jacobian plans
  that request a column whose forward absorption also responds to that
  abundance now raise an `ArgumentError` naming the omitted terms (for H₂O:
  line-by-line self broadening, MT_CKD and CIA; for a variable gas: CIA pairs
  with it as a collider). Pass `abundance_jacobians=:frozen_cross_section` to
  accept the partial (values are unchanged); fixed humidity
  (`compute_h2o_jacobians=false`), `fixed_molecules` and fixed-broadener
  tables need no change. Plans may not select q-driven H₂O columns when
  q ≡ 0. In one HITRAN2016 check the native H₂O column differed from the
  complete response by 4.4% (q = 0.0125) and 8.1% (q = 0.025).
- The q-driven H₂O Jacobian columns include the exact response of every
  ABSCO table with a multi-node H₂O-broadener axis (O₂, CO₂ and H₂O itself)
  under both `abundance_jacobians` settings, so these tables are accepted
  under `:complete`. In one check on the ABSCO v5.2 tables (RRS_XCO2
  16-layer profile), the term's share of the q-driven H₂O column, at its
  largest over the layers, reached 99 % (O₂ A), 15 % (weak CO₂) and 5 %
  (strong CO₂).
- `CoreRT.loadAbsco` reads the full 4-D ABSCO table again (NCDatasets 0.14
  flattens a variable read with `[:]`) and checks its dimensions; the `.jld2` tables built
  from it by `make_interpolation_model` hold only the H₂O VMR = 0 slice.
- A nonzero `radiative_transfer.Δ_angle` (the retired δBGE exclusion cone)
  is still accepted but is now stored as 0 with a one-time warning, so the
  parsed parameters and their display match what the fit uses.
- Documentation: `default_solar_spectrum_at_earth` states its units
  (photons s⁻¹ m⁻² µm⁻¹, not `SolarBeam` F₀ units), and the `BlackbodySource`
  `factor = π` is described as the exitance πB of a Lambertian emitter (use
  `scale = (R/d)²` for a distant disk).
- Layers that do not scatter (τϖ ≤ 2eps, common in Float32 and wherever
  Rayleigh is zeroed) no longer re-add the previous layer's sources, and the
  surface is coupled to the whole atmosphere above it (previously the lowest
  layer's scattering trait was reused, dropping surface–atmosphere coupling),
  in the forward, split, multisensor and linearized paths. In the forward
  path such layers now also emit their thermal source (1−ϖ)B(1−e^{−τ/μ}).
  **Radiances change for columns with non-scattering layers**: BOA was
  overestimated by ≈ 51 % at one view in a three-layer test, and Float32
  thermal-IR TOA was ≈ 3× too large (it now matches Float64).
- `RadauQuad` at SZA = 0 uses correctly mapped, halved weights (they summed
  to 2 and gave negative radiances); near-overhead suns no longer throw:
  coincident nodes are merged and μ₀ is kept.
- Cox–Munk BRDFs return 0 instead of NaN at extreme grazing in Float32,
  where the slope PDF and the weight denominator underflow together.
- Edge-case numerical fixes:
  - RPV and Ross–Li surfaces no longer produce intermittent NaNs (`311198c2`,
    `faf31d2d`). The Fourier accumulator is zero-filled (v2.1.0 multiplied
    uninitialized memory by 0), and `acos`/`sqrt` arguments are clamped.
    Finite results are unchanged.
  - Linearized elemental layers keep their ϖ and phase partials where ϖ or
    the projected Z·F₀ entry is zero (Suniti Sanghavi, `35a3495d`); v2.1.0
    set them to zero there.
  - Interpolation-model cross sections requested on a wavelength grid
    (`wavelength_flag=true`) are returned in the caller's order (Suniti
    Sanghavi, `092ff2c1`); v2.1.0 returned them reversed. The RT path uses
    wavenumber grids and is unaffected.
  - `mol_name(mol, iso)` and `HumlicekErrorFunction` no longer raise a
    `MethodError` (`3d86651d`).
- Tests: the Jacobian/forward helpers are one `VSmartMOMTestSupport` module
  (`run_linearized_rt`, `linearized_R`, `run_forward_rt`, `forward_R`,
  fail-closed `rel_errors`); NaN-tolerant accuracy checks now require finite
  values; tiled CUDA Jacobian tests run once and global switches are restored
  in `finally`; the local data runner fails on unexpected errors and pins
  documented legacy failures; `test/README.md` lists every test file.
- New surface-pressure Jacobian convention `FixedAltitudeGrid()` (opt-in,
  `model_from_parameters(LinMode(), params; vertical_grid=...)`): every
  interface stays at its geometric altitude, so every layer's pressure and
  column scale with p_surf; `surface_pressure_interfaces(grid, p_half,
  p_surf)` rebuilds the matching interface pressures. `FixedPressureGrid()` (only the
  bottom boundary moves) remains the default and is unchanged. The Round
  7_altgrd pilot campaign and its report are under
  `sandbox/workflows/RRS_XCO2/inversion/`.
- RRS_XCO2 workflow: six scripts that assumed data next to the code (the
  Round 4/5 prior builders and their tests, `validate_bottom_layer_retrievals`,
  and the `RetrievalCases` `FULL_COLUMN_TRUTH_ROOT` fallback) now resolve it
  below `RRS_XCO2_DATA_ROOT`, and the OCORaman workflow project and scripts find the
  repository from `sandbox/workflows/OCORaman`. The workflow keeps the
  converted dry ABSCO `.jld2` tables as its default; the original 4-D tables
  enter through the v2.2 rerun campaign's own pinned environment.
- TOML configurations accept arrays that mix integers and floats (e.g.
  `vza = [60, 45.5]`), as YAML already did; `${ENV:NAME}` markers in
  `cia_files` and `mtckd_file` now resolve like those in `LUTfiles`; surplus
  `surface` entries beyond the number of bands are reported once (they were
  silently ignored), and fewer surfaces than bands raise a clear error at
  model construction instead of a `BoundsError`; and the linearized gas
  tangents use the model float type unless the CO2 VMR is a ForwardDiff Dual
  (a Float64 or integer `CO2` VMR left Float64, abstract or failing tangent
  arrays in Float32 models).
- A reused `RS_type` no longer keeps the solar F₀ that an earlier call
  wrote into it (a call with `sources=SolarBeam(F₀=2F)` made later default
  calls return 2×). Precedence is the `sources=` keyword, then a user-preset
  `RS_type.F₀`, then `model.sources`. A `sources=` keyword applies to its
  call only: `rt_run`, `rt_run_ss`, `rt_run_toa(rs::RRS, …)` and the
  linearized `rt_run` restore `RS_type.F₀` and its driver-write record
  afterwards, also when the call throws, so a preset F₀ survives. The call's
  own results are unchanged. Split solves with TMS take the TMS beam from the
  captured `surface_F₀` instead of reading `RS_type.F₀` back (same values; no
  result changes).
- HDR, BHR (up and down) and `rt_run_ss`'s `hem_T` no longer double-count
  the direct beam when μ₀ is a weighted quadrature node (Gauss with the SZA
  on a node, every Radau run): e.g. vacuum Lambertian HDR was 1 + w_μ₀ too
  large at SZA 60° with Gauss. `hem_R`/`hem_T` are documented as diffuse
  fluxes in F₀ units, not normalized reflectances.
- An SZA or VZA whose cosine lies within 4·eps of a quadrature ordinate
  (e.g. `sza = acosd(node)`) no longer adds a numerically coincident second
  ordinate (in the worst sampled case TOA jumped by ≈ 9 % and HDR by ≈ 14 %);
  within
  that window the output now equals the on-node result, and other grids are
  unchanged. Just outside it the elemental factors are now evaluated in
  limit form (next item).
- The elemental t⁺⁺/J₀⁺ factor μⱼ/(μᵢ − μⱼ)(e^{−τ/μᵢ} − e^{−τ/μⱼ}) and its
  Raman and canopy analogues are evaluated as a divided difference of e^{−x}
  (`expdiff_ratio`, exact at μᵢ = μⱼ) in every elemental kernel (forward,
  linearized, canopy and Raman; CPU and GPU). In the tested cases away from
  coincident ordinates (pure Rayleigh IQUV, Radau, `rt_run_ss`, the
  JacobianTestFast linearized and forward runs, a canopy) results change by
  about 5e-15 relative in Float64 and up to 6e-6 in Float32; the Phase 1b RRS
  band (Float32) changes by up to 1.8e-6 elementwise, because the former
  Raman 1e-8/1e-6 coincidence windows are gone. The VS and VS `_plus`
  kernels were checked against BigFloat at kernel level only (no RT-level
  comparison). Near a node results are now
  accurate: in the quickstart Rayleigh scene (Stokes_IQU, Gauss nstreams 3,
  μ₀ = 0.5 + δ) the TOA radiance differed from the on-node result by 7.7e-5
  (max-normalized) at δ = 5·eps in Float64 and by 1.5e-4 at δ ≈ 1e-6 in
  Float32; with the limit form these are 2.1e-15 and 1.7e-6, the physical
  change (≈ 1.7·δ). With one aerosol in that scene and δ = 1e-12, the
  aerosol τ_ref column of the TOA Jacobian differs from a central difference
  by 2.5e-9 relative (was 1.2e-2), and a Radau nadir view at SZA 1e-5° from
  2.2e-6 to 8.8e-12 against SZA 0. Canopies with horizontal leaves (G ∝ μ) no
  longer hit 0/0. `rt_close_tol` and `rt_loose_tol` are kept but no longer
  used by live code.
- TMS (`ss_correction`) is evaluated per band in combined-band solves (it
  raised a `BoundsError`) and is rejected with `rt_run_streams`/stream
  callbacks, whose reconstruction it broke. Stream exports do not include
  the Cox–Munk glint residual (documented).
- Interior-height (multisensor) runs apply the Cox–Munk glint correction to
  TOA, matching the column path (Stokes I differed by up to 2.2 %, Q by
  more). The glint
  correction also works when optical depths are promoted beyond the model
  precision (Float32 model, Float64 aerosol inputs).
- **New rejections** with an `ArgumentError` instead of a silently wrong
  result: incident F₀ with nonzero U or V (azimuthal reconstruction assumes
  I/Q incidence); active `ThermalEmission` with `rt_run_ss` or the Raman
  solvers (it was ignored; the empty `ThermalEmission()` placeholder is still
  accepted); a `ThermalEmission.B_layer` whose row count differs from the
  model layer count (the one-row zero placeholder is still accepted).
  Linearized runs on surfaces without a linearized layer (Spectrum, Spline,
  RPV, Ross–Li) now fail up front with a clear `ArgumentError` instead of a
  deep `MethodError` (new trait `supports_linearized_surface`).
- **H₂O continuum, CIA and humid-layer columns change** (Suniti Sanghavi,
  `4774afd4`): the H₂O/dry-air molar ratio r is converted to the moist-air
  mole fractions x_H₂O = r/(1+r) and x_dry = 1/(1+r) before the hydrostatic
  column is split (v2.1.0 used r and 1 − r). At q = 0.0125 (0.025) the H₂O
  column falls by 2.0 % (4.0 %) and the dry column rises by 0.03 % (0.11 %).
  MT_CKD uses the hydrostatic H₂O column with the partial pressures p·r/(1+r)
  and p/(1+r), instead of n·Δz with p·r and p(1 − r). Before the separate
  T_ref/T change below, the self continuum falls by 4.0 % (7.8 %) and the
  foreign continuum by 2.0 % (3.9 %). Dry–dry CIA uses the dry number density
  and the hydrostatic dry column instead of n_air²·Δz: −4.0 % (−7.8 %). The
  continuum and CIA percentages are for thin layers: v2.1.0's midpoint column
  n·Δz exceeds the hydrostatic one by the factor p̄·ln(p₂/p₁)/Δp (1.0009 for
  900–1000 hPa, 1.040 for 500–1000 hPa, where the self continuum then falls by
  7.6 % at q = 0.0125). Line-by-line H₂O receives the layer's moist mole
  fraction, so its self-broadening is included (v2.1.0 passed no H₂O
  abundance).
- CIA and MT_CKD input (`4774afd4`): CIA files can select HITRAN reference
  families (`reference_codes`) and a `negative_policy`. Overlapping reference
  families without `reference_codes`, and negative cross sections on the
  model grid, now raise errors. `load_mtckd` requires the file's `Title` and
  `Version_description` attributes.
- MT_CKD H₂O continuum: the self and foreign terms now carry the
  number-density ratio `(p/p_ref)(T_ref/T)` of AER's `mt_ckd_h2o` (4.3);
  only the pressure ratio was applied. **Continuum optical depths increase
  by T_ref/T − 1 per layer** (≈ +3 % at 288 K, ≈ +36 % at a 217 K
  tropopause; ≈ +6.7 % for a 33-layer US-standard column at 1.6 and 2.0 µm).
  This factor is 1 at 296 K; against v2.1.0 the continuum also changes
  through the moist-air columns and partial pressures of the previous item.
- **Breaking:** the RT solvers reject `Stokes_IQ` (forward, linearized,
  `remake_geometry` and Raman model builds raise an `ArgumentError`; use
  `Stokes_IQU`). Neglecting U in multiple scattering can cause intensity
  errors of several percent: in one 770 nm scene with aerosol τ = 0.3,
  `Stokes_IQ` differed from `Stokes_IQU` by up to 4.2 % in I and 7.6 % of I
  in Q, largest in the principal plane (scalar `Stokes_I`: up to 5.7 %). Its phase-matrix
  moments also omitted the `T ζ T′` term from the Q←Q element. The parser,
  the Scattering-level type and StandaloneSS `ExactSSConfig` single-scattering
  runs still accept it (an RTModel passed to `run_exact_ss` cannot carry it).
- HITRAN, partition-sum, HAPI, CIA and MT_CKD references added to
  `CITATION.bib`, the References page and the README, with guidance on which
  HITRAN edition a calculation used.
- **Forward Mie aerosol optics change** (`e46281d2`; anchor grids, the λ_ref
  node and per-node truncation by Suniti Sanghavi, `b524a0d8`, `b5f6e425`,
  `7ef84106`): each
  aerosol's Mie calculation uses its own refractive index in every band.
  v2.1.0 wrote `n_ref` into the `Aerosol` object it then reused for the band
  loop, so every band's Mie used `n_ref`, by default the first aerosol's
  index. Scenes whose aerosols have different indices, or an explicit
  `n_ref` that differs from an aerosol's index, change at first order; a
  single aerosol with the default `n_ref` does not change through this. In
  bands with more than one spectral point the optics now vary across the
  band: Mie runs at the band (or `aerosol_anchor_bands`) endpoints, and at
  λ_ref when it lies inside the band; each node is truncated separately; and
  k_ext, k_sca (ϖ = k_sca/k_ext), fᵗ and the Greek coefficients are
  interpolated in wavenumber (v2.1.0: one Mie calculation at the band
  centre). **Breaking for code that indexes `model.τ_aer[iBand]`**: it is
  `[iAer, nSpec, iZ]` (v2.1.0: `[iAer, iZ]`).
- **Aerosol nᵣ, nᵢ and size Jacobian columns change** (Suniti Sanghavi,
  `7ef84106`): the reference extinction k_ref is computed with the configured
  `n_ref` and held fixed with respect to nᵣ and nᵢ (its size derivatives are
  kept), so the nᵣ and nᵢ columns no longer contain ∂k_ref/∂n. v2.1.0 computed
  each aerosol's k_ref with its own index and differentiated it, which held
  that aerosol's AOD at λ_ref fixed, and its linearized AOD normalization
  differed from the forward one whenever `n_ref` was not the aerosol's index.
  Forward and linearized constructors now use the same k_ref. The
  size-distribution tangents are divided by the quadrature sum S of the size
  distribution over the radius nodes: dw = (du − w·Σdu)/S. v2.1.0 divided by
  the already normalized sum, so its size columns were scaled by S (S ≈ 1
  unless `r_max` cuts the distribution; the radius nodes themselves changed
  later, see the fixed radius grid item).
- The spectral angular cache is sized to avoid spectral-quadratic growth.
- The experimental aerosol ingestion adapter now rejects calls that previously
  returned placeholder optics. Production Mie calculations remain in `Scattering`.
- Numerical conversion and CUDA inverse dispatch ambiguities are resolved;
  Aqua's ambiguity check is enabled.
- RT batched multiplication is package-owned, preserving NNlib's behavior in
  other packages. GPU tests fail on compilation/scalar-indexing errors when
  hardware is functional.
- Taplo 0.10.0 is pinned and checksum-verified; CI checks maintained TOML
  formatting, scene schemas, and negative configuration cases.
- Julia 1.13 is supported (`julia = "1.10, 1.11, 1.12, 1.13"`; `SHA` compat
  admits the 1.x standard library shipped with 1.13) and tested in CI.
- A HITRAN download rejected by validation no longer leaves its temporary file
  in the cache directory on Windows: the validator held the file open after
  the error, so it could not be deleted. The HITRAN metadata reader and
  `geoschem_to_dict` likewise close their files when reading fails, and
  quoted non-ASCII metadata values no longer raise a `StringIndexError`.
- Fourier convergence options are represented in the schema. Unknown numerics
  keys and convergence thresholds supplied with `all` now raise errors.
- Documentation distinguishes supported elastic Jacobians, forward Raman,
  forward thermal sources, and data-dependent experimental retrieval campaigns.

### Earlier integrated changes

- **RRS GPU performance (~90× on the O2A Raman LUT; ~120× whole-`rt_run`).**
  The rotational-Raman adding–doubling GPU path was rewritten: fused
  KernelAbstractions kernels replace the per-Δn `batched_mul` loops for the
  doubling n-loop, the `ScatteringInterface_11` matrix outputs, and the SFI
  source vectors (`src/CoreRT/CoreKernel/fused_raman_kernels.jl`); `interaction!`
  now defaults to a GPU-only (non-staged) workspace with a VRAM-aware fallback;
  and inelastic postprocessing reduces over Raman shifts once instead of per
  view-angle. Net on one L40S (Float32): a two-band O2A Raman LUT scene drops
  from ~1.08 h to ~43 s; peak GPU memory ~7× lower (substantially larger
  `nSpec` runnable). Results are bit-equivalent up to Float32 reduction reorder
  (elastic R/T bit-identical, inelastic ~1 ULP in Float64) — **GPU-fused configs
  need RRS goldens regenerated**. Also includes band-local `expandBandScalars`
  and spectral-F0 flat-Z correctness fixes, and a GPU regression suite locking
  each fused kernel against the `batched_mul` reference. CPU and large-`NquadN`
  configs fall back to the unchanged path.

- **AtmosphericAbsorption.jl engine swap.** `model_from_parameters` now
  delegates LBL gas absorption to the external
  [AtmosphericAbsorption.jl](https://github.com/RemoteSensingTools/AtmosphericAbsorption.jl)
  package (TIPS-2021 partition sums).  The internal `Absorption` module remains
  the standalone σ(ν, T, p) API (TIPS-2017).  High-temperature results
  (> ~500 K) may differ slightly from earlier releases. The parsed
  `broadening_function`/`CEF` fields hold AtmosphericAbsorption types
  (`Voigt()`, `Lorentz()`, `Doppler()`, `HumlicekWeideman32()`; `cfa632b6`);
  parameters built in code must use these instead of the
  `vSmartMOM.Absorption` types.

- **GPU Mie dispatch (10.7–12.9×).** `make_mie_model(...; architecture=GPU())`
  now runs NAI-2 Mie on GPU automatically. Metal supports the native
  Float32 policy; PCW retains its CPU path. With `architecture: GPU()` the
  forward `model_from_parameters` computes the in-band NAI2 Mie on the GPU
  (`995ea3e0`, `13bfb968`; the reference extinction stays on the CPU), with
  the precision policy following the float type: `NativeFloat64` for
  Float64, `DSEmulated` for Float32. The linearized constructor keeps CPU
  Mie, so forward and linearized aerosol optics from one GPU configuration
  can differ at GPU/CPU Mie precision.

- **`BatchContext` batch API.** `BatchContext` / `update_model!` /
  `update_aerosol_loading!` / `update_aerosol_microphysics!` enable efficient
  multi-scene loops with shared Mie/Fourier/HITRAN caches.

- **Six correctness fixes.**
  - Wing-cutoff window: wavenumber-space (was wavelength-space) — standalone
    Absorption API.
  - Canopy soil-albedo Jacobian axis fixed.
  - VS Raman call-chain + atomic-mass factor corrected (vibrational Raman
    results change; RRS unaffected).
  - Nonzero δBGE `Δ_angle` is retired: it warns and is forced to zero so
    forward and linearized truncation share one angular domain.
  - `LambertianSurfaceSpectrum` surface layer now works (missing
    `create_surface_layer!` method added, `e6af91f1`; see the Lambertian
    scaffold item under Consolidated release).
  - Mie Dₙ recurrence: `ForwardDiff.Dual` refractive indices run the
    downward recursion in native arithmetic (v2.1.0 cast it to
    `Complex{Float64}`, a `MethodError` for Duals); plain complex indices
    keep the Float64-stabilized recursion. The `autodiff=true` keyword of
    `compute_aerosol_optical_properties` is retired (see the configuration
    and API item).

- **Performance work.** Flat-Z zero-copy path; m-invariant Fourier-loop cache;
  GPU sync stripping (~19.5% faster forward RT on O₂-A test case).

- **~2,300-line dead-code purge** (retired `rt_run_canopy`, `RRS_plus`,
  legacy lin entry points).

- **Compat widenings.** LogExpFunctions `0.3`; Parameters.jl `0.13`.

- **Folded PRs #222–#225.**

---

## v2.1.0 — Fourier-stream resolution + source-term refactor

> **vSmartMOM** = **V**ector **S**imulated **M**easurements of the **A**tmosphere using **R**adiative **T**ransfer based on the **M**atrix **O**perator **M**ethod.

**Schema migration release.** Treat as breaking for downstream code
that hard-codes `params.max_m` / `params.l_trunc` field reads, raw
access to `SolarBeam`'s `F₀` on `RS_type`, or `Δ_angle = 2.0` as the
implicit default.  The legacy YAML parser branch (`max_m`/`l_trunc`)
is preserved for tactical reuse, but every in-tree YAML has been
migrated to `nstreams` + `truncation: auto`; users following examples
will land on the v0.7 schema.

### What's new since v2.0.0 — at-a-glance

If you last looked at `vSmartMOM.jl` on `main`, this is what changed.
Most items have their own subsection further down, or a dedicated
docs page in `docs/src/pages/`.

| Theme | Summary |
|-------|---------|
| **Source-term abstraction (v0.6)** | First-class `AbstractSource` types — `SolarBeam`, `BlackbodySource`, `SurfaceSIF`, `NoSource`, composed via `+` into `SourceSet`. Replaces the implicit `RS_type.F₀` ownership. Forward/linearized RT now flow through `prepare_source` → `contribute!` / `source_tangent!` dispatchers. AD seam at `prepare_source` (user-parameter space ⇄ kernel space). See `docs/src/pages/extending/sources.md`. |
| **Fourier / stream resolution schema (v0.7)** | `nstreams` is the primary user input; `stream_l_cap = 2·nstreams - 1` regardless of quadrature. `truncation: auto` (conservative in v2.1, see below) and `truncation: "δBGE(L, Δ)"` are the recommended idioms. Per-component `component_m_max(c, ctx)` traits — Cox-Munk / RPV / RossLi / canopy now run to their full `user_l_cap` instead of being silently half-truncated. |
| **VLIDORT 2.8.3 baseline suite** | New `test/vlidort_baseline/` with Siewert 2000 Problem IIA (Stokes-I), and `solar_tester` scalar (Task 1) + vector (Task 1, IQU) reference data committed in tree. No PyVLIDORT / Fortran runtime needed for regression. |
| **Stokes IQ polarization** | New `Stokes_IQ()` polarization type — 2-channel I/Q, between scalar I and IQU, useful for instruments with linear polarization but no Q/U separation. Forward, single-scatter, surface BRDFs (incl. Cox-Munk SS wind Jacobian), and Greek single-scatter all covered. |
| **Standalone exact single-scatter** | `StandaloneSS` driver (`run_exact_ss`, `ExactSSConfig`, `SSGeometry`) for high-precision SS-only computations independent of the adding-doubling solver. RTModel SS adapter included. |
| **Aerosols module** | Multi-moment + TOMAS-15 aerosol size-distribution support; analytic phase-function aerosols (Henyey-Greenstein); `compute_aerosol_optical_properties` exported; vector phase + surface BRDF chain rules wired through. |
| **Canopy surface (PROSAIL/4SAIL)** | `CanopySurface` BRDF integrated into the standard `rt_run` pipeline (legacy `rt_run_canopy` removed). Canopy block in YAML schema; per-band soil + leaf reflectance/transmittance. |
| **HITRAN edition selection** | `set_hitran_edition!`, `get_hitran_edition`, `available_hitran_editions`, `fetch_hitran`, `fetch_hitran_by_ids` — fetches from hitran.org with SHA-256 provenance. Defaults to the legacy 2016 artifacts; switch with one call. |
| **Numerics tuning surface** | `RTNumericalParameters` carries `dτ_max_threshold`, `dτ_min_floor`, `blas_threads` knobs that previously lived as hard-coded `0.001` / `1024·eps(FT)`. YAML key `radiative_transfer.numerics` exposes them. |
| **Documentation rewrite** | DocumenterVitepress migration, narrative-thread restructure (Concepts arc), per-block schema docs under `docs/src/pages/IO/Schema/`, JSON-Schema autocomplete for TOML/YAML, auto-generated benchmarks page (Siewert IIA et al). README + landing page now expand the **vSmartMOM** acronym (Vector Simulated Measurements of the Atmosphere using Radiative Transfer based on the Matrix Operator Method). |
| **`RTModel` hierarchy (v2.0 finalised)** | `model_from_parameters` returns a hierarchical `RTModel{ARCH, FT}` with physics-based sub-structs (`SolverConfig`, `Atmosphere`, `Optics`, `Surfaces`, `Sources`). `Base.getproperty` shim provides backward-compat field access. |
| **CI hardening** | Test matrix now Julia 1.10 / 1.11 / 1.12; `actions/checkout@v4`; `julia-actions/cache@v2` step; Aqua quality gate (stale deps, ambiguities, piracy, undocumented exports) included in `runtests.jl`. |

### Breaking changes vs v2.0.0

These will require user action if you have downstream code:

- `params.max_m` and `params.l_trunc` are still parsed but no longer
  the canonical knobs.  In-tree YAMLs all use `nstreams` now.
- `SolverConfig.max_m` / `get_max_m(model)` have been removed in
  favour of `m_max_bands` (per-band order). If you read these,
  replace with `model.solver.m_max_bands[iBand]`.
- `RS_type.F₀` is no longer the public interface for solar irradiance
  — pass `sources = SolarBeam(F₀ = …)` to `model_from_parameters` or
  `rt_run` instead.  Existing code that wrote to `RS_type.F₀`
  directly still works in the legacy code path but the field is
  scheduled for removal in a follow-up.
- `Δ_angle` is now an optional YAML field, default `0.0` (was a
  required field with implicit `2.0` examples).  Set explicitly when
  reproducing VLIDORT `DO_DELTAM_SCALING = 2°` benchmarks.
- `rt_run_canopy` (canopy-only entry point) was removed — call
  `rt_run` with `CanopySurface` as the BRDF.
- `parameters_from_yaml` is now YAML-only; use `parameters_from_file`
  / `read_parameters` for TOML or registry-dispatched inputs.
- `RRS_plus` has been retired. Rotational Raman remains available via
  `RRS`; concatenated-grid inelastic workflows should use the
  vibrational plus modes.

### Numerical changes

- The doubling-step `dτ_max` formula now filters zero-weight quadrature
  streams (user-supplied VZA / SZA nodes are appended for output
  postprocessing but should not constrain `dτ_initial`) and applies an
  absolute `dτ_min_floor` of `1024·eps(FT)`. This affects results at
  grazing geometries (e.g. `vza ≈ 89.99°`) at the level of single-bit
  rounding — large enough to flag in floating-point comparisons against
  pre-v2.1 reference outputs, but well within published-benchmark
  precision.  Threaded into `rt_kernel!`, `rt_kernel_ss!`,
  `rt_kernel_multisensor!`, and `rt_kernel_lin!` so the user-set
  `numerics.dτ_max_threshold` / `dτ_min_floor` reach all four paths.

### Highlights

First-class `AbstractSource` type vocabulary (`SolarBeam`,
`BlackbodySource`, `SurfaceSIF`); user-facing `nstreams` /
`truncation: auto` Fourier-resolution schema; per-model BLAS thread
cap; VLIDORT 2.8.3 baseline validation suite (Siewert 2000 Problem IIA,
solar_tester scalar/vector); `Δ_angle` is now optional and defaults to
`0.0` (legacy `2.0` matches VLIDORT `DO_DELTAM_SCALING` and can still
be set explicitly).

### New user-facing schema (Phase D)

- **`nstreams` is the primary resolution knob.** Public contract:
  `stream_l_cap = 2·nstreams - 1` regardless of `quadrature_type`.
  Minimum 3 (Rayleigh contributes through m=2). Default 8 when
  legacy `max_m`/`l_trunc` are absent — gives `stream_l_cap = 15`,
  which fits typical aerosol δ-fit truncations and gives ~4-decimal
  Rayleigh convergence at standard atmospheric depth.
- **`truncation: auto`** — VLIDORT-`DO_DELTAM_SCALING`-style.
  `NoTruncation()` for Rayleigh-only scenes; `δBGE(stream_l_cap,
  Δ_angle)` once aerosols are configured. The per-band Mie loop
  applies δBGE only when `length(greek.β) > stream_l_cap`, so a
  small-particle Mie series shorter than the cap automatically
  falls back to `NoTruncation()` — no crash, no user action
  needed. Logs the chosen branch via `@info`.
- **`quadrature_type` defaults to `GaussLegQuad()`** when omitted
  (Sanghavi: 5–50× more accurate per stream than `RadauQuad` on
  Rayleigh-only). `RadauQuad` remains supported but is documented
  as expert/legacy.
- **JSON Schema updated** —
  `schemas/vsmartmom-parameters.schema.json` covers the new fields.
  The repo's `.taplo.toml` already wires it to TOML configs; YAML
  editors honour it via `# yaml-language-server: $schema=...`. See
  `docs/src/pages/IO/Schema.md` for the per-block reference and
  editor-setup recipes.

### Internal architecture (Phases A–C)

- **`Nstreams` field on `QuadPoints`** (Phase A) — count of nonzero
  weights, distinct from the augmented `Nquad` total. Surfaces in
  `Base.show` and `conventions.md` §6.
- **`m_max_bands` order semantics on `SolverConfig`** (Phase B) —
  replaces the scalar `max_m` and the count-form `max_m_bands`.
  Forward and lin RT paths now flow through `_derive_m_max_bands`,
  fixing a latent precedence bug at even `l_max` in the lin path
  (`Int(ceil(l_max+1)/2)` → `Int(ceil((l_max+1)/2))`).
- **`component_m_max(c, ctx)` traits** (Phase C) — per-component
  Fourier support dispatch (Rayleigh → 2, Lambertian → 0,
  CoxMunk/RPV/RossLi/Canopy → `user_l_cap`, SolarBeam → 0,
  SurfaceSIF → 0). Resolves Phase B's regression for Cox-Munk
  forward; the Phase C P2 follow-up clamps trait-derived `m_max`
  against the actual quadrature `Nstreams` cap and includes active
  sources in the aggregation.

### Bit-equality

- All 64 in-tree YAMLs run byte-equal through Phase D — they carry
  explicit `max_m` / `l_trunc` so the parser hits the legacy branch.
- Phase B fixed the lin-only precedence bug at even `l_max`; lin
  Jacobians shift by one Fourier moment for affected configs
  (intentional silent-bug fix).
- Phase C+P2 unified forward and lin Fourier-loop bounds. Cox-Munk
  forward goes from the half-truncated count-aggregator output
  back to its full `user_l_cap` resolution.

### See also

- Plan: `docs/dev_notes/fourier_stream_resolution_plan.md`
- Schema index: `docs/src/pages/IO/Schema.md`

## v2.0.0

### Breaking Changes

- **`RTModel` replaces `vSmartMOM_Model`**: `model_from_parameters()` now returns an `RTModel{ARCH, FT}` with physics-based sub-structs (`SolverConfig`, `Atmosphere`, `Optics`, surfaces). The old flat `vSmartMOM_Model` struct has been removed. A `Base.getproperty` shim provides backward-compatible field access (e.g. `model.τ_abs`, `model.profile`).

- **`RTModelLin` replaces `vSmartMOM_Lin`**: The linearized model container returned by `model_from_parameters(LinMode(), params)` is now named `RTModelLin`.

- **`rt_run_canopy` removed**: Canopy RT is now handled through the standard `rt_run()` pipeline using `CanopySurface` as the surface type. The standalone `rt_run_canopy()` function has been removed.

- **`parameters_from_yaml` is YAML-only**: Use `parameters_from_file` or `read_parameters` for TOML or registry-dispatched inputs. Passing a non-`.yaml`/`.yml` path to `parameters_from_yaml` now raises `ArgumentError`.

### New Features

- **Hierarchical model architecture**: `RTModel` organizes state into `solver`, `geometry`, `quad_points`, `atmosphere`, `optics`, and `surfaces` sub-structs, clearly separating fixed configuration from differentiable state.

- **`ParameterLayout`**: New struct for Jacobian index arithmetic — provides `aerosol_range()`, `gas_range()`, `surface_range()`, `n_total()` accessors.

- **Task-based documentation**: The manual now starts from runnable tasks (quick start, scene configuration, Jacobians, GPU, and extension guides) and keeps module pages as reference material.

- **Module-grouped Library**: Public docstrings are organized by module in a dedicated Library section, with a generated function index for lookup.

- **Repository citation file**: `CITATION.bib` now ships at the repository root with the JOSS software paper and core method references.

- **Cox-Munk ocean surface**: Full polarized BRDF with Fresnel reflection, wind-speed-dependent roughness, whitecap correction, and Smith (1967) shadowing. Linearized variant included.

- **Accessor functions**: `architecture(model)`, `CoreRT.polarization_type(model)`, `CoreRT.float_type(model)`, `CoreRT.n_aerosols(model)`, `get_surface(model, iBand)`, `get_spec_bands(model)`.

### Improvements

- Julia compat is declared for Julia 1.10, 1.11, and 1.12. CUDA remains a weak dependency and now allows CUDA.jl 6.
- `read_parameters` is the preferred public loader for file paths, dictionaries, and `IOSource` values. `parameters_from_file`, `parameters_from_dict`, and `parameters_from_source` remain available as explicit aliases.
- Parser validation now throws `ArgumentError` instead of relying on `@assert`, so configuration checks are not elided by optimized Julia runs.
- The docs include a shipped CPU quick-start scene and a small docs smoke runner for the examples most likely to drift.
- The docs build runs with `checkdocs = :exports` and no `warnonly` suppressions.
- Canopy surface integration now targets CanopyOptics.jl 0.2, including the
  analytic Stokes-aware canopy Z-matrix API used by `CanopySurface`.
- Replaced `eval(Meta.parse(...))` in YAML spec_bands parsing with a safe literal parser.
- Fixed `Vector{Any}` fields in Raman types to use proper parametric types.
- Cleaned up dead code: removed unused prototype files, plot scripts, and old test utilities.
- Cleaned up dated and stale TODO comments across inelastic scattering and surface modules.

### Known Limitations

- Linearized RT (Jacobians) does not yet support Raman scattering (`RS_type = noRS()` only).
- SIF helper exports exist, but the end-to-end SIF workflow remains experimental until fixture/data policy is settled.
- EMIT/MODTRAN comparison scripts are parsed in CI but not executed by default because they require external LUTs and machine-local scene inputs.
- Aerosol scene input is documented, but the higher-level aerosol API is still being stabilized.
- Wigner 3j symbol validation test remains disabled by default (~60s runtime).
