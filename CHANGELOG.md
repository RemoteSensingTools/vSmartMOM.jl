# Changelog

## 2.2.0 — release candidate (unreleased)

Scientific authorship of the linearization and Raman work, including the
formulations, analytic derivatives, and retrieval methodology, is credited to
**Suniti Sanghavi**. Christian Frankenberg assists with performance optimization,
backend engineering, integration, packaging, and release tooling. See the
release manual's contribution statement
and the published references in `CITATION.bib`.

### Consolidated release

This candidate combines the surface/multisensor integration, Suniti's updated
retrieval workflows, and the validated Jacobian improvements. It preserves the
published `v2.1.0` tag. The migration baseline for registry users is **1.1.0**;
2.0/2.1 were not registered in General. Registration of 2.2.0 will require
review of General's sequential-version rule.

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
  `UsingPolTen()` uses the Raman-module polarizability tensors for N₂, O₂ and
  H₂. `molecular_rayleigh_provenance` covers `UsingRefInd()` for Ar and CO₂. The refractive-index part is ported from Suniti Sanghavi's
  Habitable Worlds Observatory branch. Model defaults are unchanged: Bodhaine
  (1999) optical depth, with depolarization from the polarizability constants
  when `depol: -1`.
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
- `CanopySurface` leaf scattering now uses the solar Stokes vector `F₀` of
  the solve, as the soil below it already did. It used a unit beam, so canopy
  radiances were not proportional to a spectral solar irradiance or a
  `SolarBeam(F₀=…)` scaling (e.g. R(2F₀)/R(F₀) was 1.34–1.45 instead of 2),
  and a zero solar beam still gave leaf-scattered radiance. Runs with the
  default unit beam are unchanged.
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
- **Breaking for some linearized callers:** gas-abundance Jacobian columns
  are the native partial ∂τ/∂VMR = σ·N_dry
  with cross sections, CIA and MT_CKD held fixed. Solves and Jacobian plans
  that request a column whose forward absorption also responds to that
  abundance now raise an `ArgumentError` naming the omitted terms (for H₂O:
  self broadening, multi-node ABSCO H₂O-broadener axes, MT_CKD and CIA; for
  a variable gas: CIA pairs with it as a collider). Pass
  `abundance_jacobians=:frozen_cross_section` to accept the partial (values
  are unchanged); fixed humidity (`compute_h2o_jacobians=false`),
  `fixed_molecules` and fixed-broadener tables need no change. Plans may not
  select q-driven H₂O columns when q ≡ 0. In one HITRAN2016 check the native
  H₂O column differed from the complete response by 4.4% (q = 0.0125) and
  8.1% (q = 0.025).
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
  `RS_type.F₀`, then `model.sources`.
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
  unchanged. Just outside it the elemental 1/(μᵢ − μ₀) cancellation remains
  (≈ 4e-4 relative at 5·eps in Float64).
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
- HITRAN, partition-sum, HAPI, CIA and MT_CKD references added to
  `CITATION.bib`, the References page and the README, with guidance on which
  HITRAN edition a calculation used.
- Shared aerosol reference normalization and spectral angular-cache sizing
  fixes restore constructor parity and avoid spectral-quadratic cache growth.
- The experimental aerosol ingestion adapter now rejects calls that previously
  returned placeholder optics. Production Mie calculations remain in `Scattering`.
- Numerical conversion and CUDA inverse dispatch ambiguities are resolved;
  Aqua's ambiguity check is enabled.
- RT batched multiplication is package-owned, preserving NNlib's behavior in
  other packages. GPU tests fail on compilation/scalar-indexing errors when
  hardware is functional.
- Taplo 0.10.0 is pinned and checksum-verified; CI checks maintained TOML
  formatting, scene schemas, and negative configuration cases.
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
  (> ~500 K) may differ slightly from earlier releases.

- **GPU Mie dispatch (10.7–12.9×).** `make_mie_model(...; architecture=GPU())`
  now runs NAI-2 Mie on GPU automatically. Metal supports the native
  Float32 policy; PCW retains its CPU path.

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
    `create_surface_layer!` method added, including Lambertian scaffold).
  - Mie Dₙ recurrence: Dual-number-safe initialization restores correct
    `ForwardDiff` derivatives through Mie.

- **Performance work.** Flat-Z zero-copy path; m-invariant Fourier-loop cache;
  GPU sync stripping (~19.5% faster forward RT on O₂-A test case).

- **~2,300-line dead-code purge** (retired `rt_run_canopy`, `RRS_plus`,
  legacy lin entry points).

- **Compat widenings.** LogExpFunctions `0.3`; Parameters.jl `0.13`.

- **Folded PRs #222–#225.**

Release publication is to be coordinated with S. Sanghavi after the candidate
passes the release gates; no tag or registration is implied by this entry.

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
