# vSmartMOM v2.2 consolidation — decision log and open items

Maintained during the v2.2 consolidation (branch `release/v2.2-consolidation`,
started 2026-10-17 from `7acab570`, the Round 5 integration). Each entry says
*what* was decided, *why*, and *what to look at* if the decision needs to be
revisited. Decisions are Suniti Sanghavi's unless stated otherwise.

The full release plan is Christian Frankenberg's handoff
`RELEASE_READINESS_AND_CLAUDE_HANDOFF_2026-10-16.md`; this file records only
decisions taken on top of it.

## 1. Decisions

### 1.1 Molecular Rayleigh: two methods, selected by dispatch

**Implemented (merge a405d87c).** The HWO refractive-index
Rayleigh code (`src/Scattering/molecular_rayleigh.jl` in the HWO worktree)
is ported with two method types:

| Type | Physics | Species available now |
|---|---|---|
| `UsingRefInd()` | refractive index + (effective) King factor | He, Ar, N₂, O₂, CO₂ |
| `UsingPolTen()` | polarizability tensor (mean polarizability ᾱ, anisotropy γ) | N₂, O₂, H₂ (constants from the Raman module) |

`RefractiveIndexRayleigh` stays as an alias of `UsingRefInd` so existing
ExoOptics code keeps working. The two-argument call
`molecular_rayleigh_properties(gas, λ)` remains the refractive-index method
(ExoOptics calls it for He, which has no polarizability tensor yet).

Note on what the RT model uses today: the default *elastic* Rayleigh optical
depth is the Bodhaine et al. (1999) Eq. 30 air formula scaled by
p_surf/1013.25 and rescaled by the ratio of King factors that takes
Bodhaine's built-in depolarization (0.0279) to the configured one
(`getRayleighLayerOptProp`).
With the default `depol: -1`, that depolarization — which also sets the
Rayleigh Greek coefficients — comes from the N₂/O₂ polarizability constants
of the Raman module (`model_from_parameters.jl`, `getRamanAtmoConstants`).
The Cabannes/RRS/VRS paths use the polarizability tensors directly.

### 1.2 Common names between vSmartMOM and ExoOptics

Methods live in their most general form in vSmartMOM; ExoOptics derives from
them. Where ExoOptics duplicates a vSmartMOM method under another name, the
ExoOptics name changes to the vSmartMOM one (or delegates to it). Where no
exact parallel exists, a vSmartMOM method is added via multiple dispatch.
An audit of ExoOptics is in progress (see §3).

### 1.3 Solar direction: embedded μ₀ is the default everywhere

`external_solar` controls the *direction* of the solar beam in the
quadrature, not the solar spectrum:

- `external_solar=false` (**default**, embedded): μ₀ is an extra quadrature
  node — zero-weight with Gauss quadrature, a weighted node with Radau.
  Supports every output (TOA, BOA, interior sensors, full R/T, Radau,
  atmosphere/surface split).
- `external_solar=true` (external): μ₀ stays outside the stream set and the
  exact direct-beam source is evaluated at the requested μ₀ through
  rectangular source columns; TOA-only, Gauss-only, through `rt_run_toa`.
  It saves one operator node per solve (unless μ₀ already coincides with a
  node). Supported only for forward `noRS` and `RRS` through `rt_run_toa`,
  linearized `noRS` through `rt_run(model, lin_model, …)`, Lambertian
  surfaces, no `rt_run_ss`, no interior sensors. Retained for quick-and-dirty use only.

**Not to be confused with stream interpolation.** Reusing one solve for many
solar geometries by interpolating between quadrature nodes is a *separate*
algorithm, `rt_run_streams` (per-Fourier-moment results at all streams). It
runs with the embedded representation (`external_solar=false`). ExoOptics'
production stellar path (`PerPixelRT._per_pixel_from_streams`) uses it; its
documented interpolation error (≈9 % DOLP at α = 90° with 5 streams) belongs
to that path, not to `external_solar`.

The solar *spectrum* is a separate setting: `SolarBeam(F₀=...)` (default
F₀ = 1) or the Raman `RS_type.F₀`. Raman runs need the real solar spectrum;
ExoOptics uses F₀ = 1 so that incident starlight stays generic. Both rules
are independent of `external_solar`.

Embedded is the default for the release *and* for ExoOptics. The HWO work
should also have used embedded; the ExoOptics production paths that force
the external mode are being inventoried (§3) before switching.
The Horak pure-gas diagnostic keeps `external_solar=true` explicitly (one
keyword added to its `model_from_parameters` call).

### 1.4 Workflow layout

All workflow source lives in `sandbox/workflows/<NAME>/` with its own
documentation. Generated data stays outside git and is located through
environment variables (`RRS_XCO2_DATA_ROOT`, `VSMARTMOM_SIF_DATA_DIR`).
Scripts that read workflow configs/data are updated to these locations.

### 1.5 `Δ_angle` retired permanently

The δBGE forward-peak exclusion angle is redundant: the fit is effectively
bounded by the smallest angle of the quadrature that represents the phase
function, so a separate exclusion cone is at best redundant and at worst
error-prone (the linearized path never honoured it, so non-zero values made
Jacobians inconsistent with the forward model). Its numerical effect is
already retired in this branch (the δBGE constructor forces it to 0 with a
one-time warning, and since `83676497` the parser also stores a nonzero
`radiative_transfer.Δ_angle` as 0 with a one-time warning). Compatibility
metadata (the YAML key, the `δBGE(l, Δ)` signature, the stored field) stays
accepted and ignored, so existing scripts keep running. The user-facing
documentation and docs examples no longer set it (2.2.0 documentation pass);
scripts under `sandbox/` still carry `Δ_angle: 2.0` and only trigger the
warning.

### 1.6 `Aerosols.compute_optical_properties` placeholder

The old TOMAS-15→optics placeholder returned made-up numbers (approximate
Mie efficiencies, fixed g = 0.7, a rough number density). The release raises
an explicit "not implemented" error; production Mie optics live in
`Scattering`. Remnant tests/examples that relied on the placeholder are
removed or rewritten.

### 1.7 Other

- Release worktree: `/home/sanghavi/code/github/vSmartMOM-v2.2-release`,
  local commits only; push/PR/tag/registration need explicit approval.
- HITRAN: citation only for now (no change to data hosting or the HITRAN2016
  default).
- TMS (Nakajima–Tanaka) GPU efficiency: deferred to Christian.
- Compute: free A100 (GPU 1 on curry) may be used.
- Review gate before every commit: focused tests → independent Claude review
  → Codex review (`gpt-6.1-sol`, xhigh) → fixes → re-review.

## 2. Explanations kept for later decisions

### 2.1 Why the Round 6/7 runners are pinned to commit `7acab570`

RRS_XCO2 Round 7 is complete (80 scenes × 11 members = 880 retrievals,
finished 2026-09-28 on curry and wurst, no failures; deployment tarball in
`bottom_layer_XCO2_retrievals/round7_fixed_sif_imperfect_correction/deployment/`).

The runners (`Round7RetrievalCampaign.jl`, `round6_environment.sh`) refuse to
start unless the loaded vSmartMOM checkout is exactly `7acab570`, its tree is
clean, and its file-tree hash matches. This is a reproducibility lock: the
campaign identity hash written into every output file then proves which code
produced it, and nobody can accidentally extend Round 7 with different code.

**When this matters to you:** only if Round 7 is re-run or extended. Keep the
clean checkout `/home/sanghavi/code/github/uni_vSmartMOM_round5_jacobian`
(at `7acab570`) as the archive environment. New campaigns get their own pin
(e.g. Round 7_altgrd pins its own commit; later rounds would pin v2.2.0).

### 2.2 Surface-pressure Jacobian: what the code actually does

Both the older core lineage (`suniti_multi_sensor`, `f14371a8` — dated
2026-10-09, but branched from older core code) and the Round 5+
core (`7acab570`) use the same vertical convention: on the *final* model
grid, all interface pressures except the bottom one are held fixed, together
with T, humidity and VMR on that grid, and only `p_half[end] = p_surf`
moves. (The analytic tangent does not include how profile reduction or
observer-interface insertion would move other interfaces if the model were
rebuilt from perturbed *input* parameters; a finite-difference check must
perturb the final grid, not the input file.) Upper layers are not entirely
inert: the Rayleigh optical depth is redistributed by column fraction, and
aerosol vertical allocations are renormalized over all layers (for a profile
defined in altitude all interfaces above BOA move), so upper layers carry
small nonzero tangents. Everything above the bottom layer keeps its pressure levels and
thickness; only the bottom layer gains or loses air (its column, and its
thickness Δz). Neither version holds levels fixed in altitude.

What changed in commit `36e7d1e0` (2026-09-08, "harden v2.2 pressure, source
and batch contracts"): the bottom layer's gas *cross-section* σ(p, T) now
also responds to its shifted mid-layer pressure (dσ/dp via the LUT/line
model, `absorption_pressure.jl`). Without it, the Jacobian is not the
derivative of the model's own forward calculation. The CIA and MT_CKD
density factors were already present in the older core.

Timeline: Round 4 outputs are dated 2026-09-07 and the Round 3-style
`retrievals_acos_*` outputs 2026-09-05, i.e. before this term existed;
Rounds 5–7 ran on `7acab570`, which includes it. (To be confirmed against
each round's provenance — task in progress.)

**Fixed-altitude grid (new):** for a fixed-height flux tower or a
fixed-altitude airborne observer the layering must not move around the
observer when p_surf changes. A `FixedAltitudeGrid()` convention (levels
fixed in altitude; a p_surf change propagates hydrostatically to every
level) is implemented as an opt-in alongside the unchanged
`FixedPressureGrid()` default (`model_from_parameters(LinMode(), params;
vertical_grid=FixedAltitudeGrid())`, with `surface_pressure_interfaces` to
rebuild the interface pressures; merged in `ca6dc61c`), and Round 7 has been
repeated as a pilot ("Round 7_altgrd": one noiseless and one noisy member
per scene). Whether Rounds 3–4 (and possibly all rounds) need re-running is
decided after comparing the pilot with Round 7 (§3.1).

## 3. Open items / TODO

- [ ] **For Suniti (molecular Rayleigh, from the port review):**
  (1) O₂ King-factor constant 1.09 (He et al. 2021, as ported) vs 1.096
  (Bates 1984/Bodhaine 1999; +0.54 % in σ); (2) the RI fits are extrapolated
  beyond 468 nm (N₂) and 546 nm (O₂); (3) the Raman-module anisotropy γ has no
  wavelength or temperature dependence, which also affects the default
  `depol: -1`; (4) UsingPolTen default temperature 288.15 K; (5) H₂ has no
  independent refractive-index check here.
- [ ] **Polarizability tensors for more species.** Populate ᾱ and γ (and
  their wavelength dependence) for as many species as possible — at least
  Ar, CO₂, He, H₂O, CH₄ — for use in Rayleigh (`UsingPolTen`) and Raman
  computations, and add intercomparisons between `UsingPolTen` and
  `UsingRefInd` for every species that has both.
- [ ] ExoOptics/vSmartMOM name unification (audit in progress).
- [ ] Switch ExoOptics production paths to the embedded solar direction;
  quantify embedded vs external differences on the A–D case first.
- [ ] Fix ExoOptics provenance: `ExoEarthRunner.jl:946` records and requires
  a clean `../uni_vSmartMOM`, which is not the vSmartMOM actually loaded.
  Use the loaded package's own path.
- [ ] Move remaining workflows to `sandbox/workflows/` and update the scripts
  that reference old paths (e.g. `generate_truth_map.jl:42`).
- [x] Fixed-altitude grid + Round 7_altgrd pilot: completed (merged in
  `ca6dc61c`); see §3.1 for the results and the decisions it needs
  (re-running Rounds 3–4 / all rounds).
- [ ] TMS GPU efficiency — Christian.
- [x] **For Suniti (convention, documented 2026-10-17):** `conventions.md` §1
  said Q > 0 means polarization parallel to the meridian plane. The code's
  convention (calibrated against its own Rayleigh phase matrix and VLIDORT)
  is Q = |E_φ|² − |E_θ|², U = 2Re(E_φE_θ*), V = 2Im(E_φE_θ*), i.e. Q > 0
  *perpendicular* to the meridian plane. The page now states this; please
  confirm the wording.
- [ ] **For Suniti (surface definition):** vSmartMOM's Ross–Li caps both zenith
  angles at 80° inside the kernels; VLIDORT does not (its Li-sparse kernel
  turns strongly negative at grazing quadrature nodes). The cap changes TOA
  intensity by 0.35–0.84 % in the VLIDORT comparison (BOA by up to 1.1 %).
  Keep the cap (physical, non-negative BRDF) or follow the MODIS definition?
- [ ] ExoOptics audit item 16 (`_read_hitran2024`): delegate to
  AtmosphericAbsorption's `load_lines(HitranPort(path))`, not to vSmartMOM's
  legacy `read_hitran` (the standalone line engine is deprecated in 2.2).
- [ ] `CanopySurface` caches its optics without keying them on the spectral
  grid: one canopy object reused for bands with different grids (or after a
  combined-band solve) reuses a stale cache. Combined-band canopy solves are
  rejected for now (B8); key the cache on the grid in a follow-up. Related:
  the linearized canopy wrapper calls the forward canopy layer without
  `spec_bands_wn`/`m_max`, so a cache it builds uses the default `m_max = 2`
  (later moments are clamped to the last stored one) and no spectral leaf
  optics; fix both together (key = quadrature, spectral grid, `m_max`,
  Stokes count).
- [ ] Generic analytic surfaces (`analytic_surface.jl`, RPV and Ross–Li) pass
  the reflected (row) direction as `μᵢ`; harmless for these reciprocal scalar
  BRDFs, but fix the order and docstring in a separate commit.
- [ ] Extremely grazing Float32 Cox–Munk inputs can overflow the
  1/(4μᵢμᵣcos⁴β) weight (pre-existing; reviewer finding) — add a tail guard.
  Partly done in `2707057a`: where the slope PDF underflows to 0 the glint
  term (and its wind derivative) is now exactly 0 instead of NaN; an overflow
  with a nonzero PDF is not guarded.
- [ ] Ocean/glint products made before the Cox–Munk corrections (B1 slope
  variance, B2 shadowing, B3 polarization) are stale, e.g. the outputs of
  `sandbox/workflows/RRS_XCO2/scripts/generate_truth_map_land_glint.jl`;
  regenerate before use.

### 3.1 Added by the 2026-10-17 bug audit (results-changing; need Suniti)

The audit findings that are code defects with a clear reference were fixed
in reviewed commits (see `CHANGELOG.md`). The items below change scientific
results or definitions and wait for a decision.

- [ ] **Solver-entry conventions (fix/rt-entry-guards).** (1) Solar F₀
  precedence is kwarg > user-preset `RS_type.F₀` > `model.sources` (as the
  code always did; older comments said `model.sources` beats a preset); a
  `sources=` call still overwrites a preset permanently. (2) `rt_run_ss`
  `hem_T` is now diffuse-only (the alternative is to always add
  μ₀F₀e^(−τ/μ₀)). (3) Interior upwelling levels keep only the
  Fourier-truncated glint (forward and linearized consistent). (4) Canopy
  BHR-down omits the direct beam (canopy sets j₀⁺ = 0). (5) The
  1/(μᵢ − μ₀) cancellation outside the 4·eps merge window remains (≈ 5e-19/δ
  in Float64, ≈ 2.4e-10/δ in Float32); a full fix needs a limit form in the
  elemental kernels. Radau's own near-coincident sub-interval nodes are
  still merged only when exactly equal.
- [ ] **Sign of ε in the Greek matrix B_l (wording only).** The code stores
  `B[3,4] = +ε`, `B[4,3] = −ε` (`construct_B_matrix`), and
  `concepts/03b_scattering.md` now shows that matrix (2.2.0 docs pass; it
  previously showed −ε at [3,4]) with a note that Siewert 2000 / de Haan et
  al. 1987 print −ε at [3,4]. `conventions.md` §2 calls `+ε` at [3,4] "the
  standard Hovenier convention". Please confirm the matrix against
  Sanghavi (2014) Mie-Fourier Eq. 16 and whether that wording is right.

- [x] **MT_CKD density factor (approved by Suniti 2026-10-18, merged).**
  AER's `mt_ckd_h2o` scales both continuum terms with the number-density
  ratio `(p/p_ref)(T_ref/T)`; `compute_τ_h2o_continuum!` applied only
  `p/p_ref`, so the continuum was low by `T/T_ref` (−26 % at 220 K, 0 at
  296 K). The fix from branch `proposal/mtckd-density-ratio` is merged into
  the release. MT_CKD is used by the
  `ParamsEMIT_MODTRANcomp_newLUT*` local configs and the OCORaman LUT script
  (when `MTCKD_FILE` is set); ExoOptics computes its MT_CKD continuum
  separately. uni_vSmartMOM's own copy has the same omission (not changed
  here).
- [x] **`Stokes_IQ` (decided by Suniti 2026-10-18: disabled in the RT
  solvers).** The I/Q truncation of Π and B drops the U-mediated `T ζ T′`
  term from the Q←Q element for m ≥ 1, so IQ Z moments are not the I/Q
  block of the IQU moments (m = 0 and Rayleigh-like cases without ζ are
  unaffected). A comparison (770 nm, aerosol τ = 0.3, μ = 0.2, SZA 40°,
  VZA 35°, 12 streams, albedo 0.05) showed the larger problem is the IQ
  approximation itself: with the exact I/Q block (identical to `Stokes_IQU`
  with its U coupling removed), I still differs from `Stokes_IQU` by up to
  4.2 % and Q by up to 6.4 % of I, largest in the principal plane (scalar
  `Stokes_I`: up to 5.7 %). A principal-plane restriction would therefore
  not help. Every RT model build now rejects `Stokes_IQ`
  (`_require_rt_polarization` in `rt_set_streams`); the parser, the
  Scattering-level type and StandaloneSS `ExactSSConfig` runs keep it.
- [ ] **Raman.** R1: a narrow RRS band with no rotational line inside gets
  `n_Raman = 0` and loses the Raman energy silently (warn/error?). R2: the VS
  Q-branch Placzek–Teller factor uses `(α′/(b γ′))²` where `α′²/(b γ′²)`
  appears intended (Q-branch low by 5.6 % N₂, 11.8 % O₂). R3: VS uses the
  RRS-free Cabannes Greek coefficients for the elastic part while keeping RRS
  in the elastic SSA (probable inconsistency). R4–R8 minor (band-depolarized
  Q-branch phase, VS_0to1_plus applies no source depletion (~8e-4 of the
  energy goes to VRS with ϖ_Cabannes = 1), unused/dead overloads,
  RRS line offsets from an off-grid mean ν).
- [ ] **Water refractive index for Cox–Munk.** The built-in table departs from
  Hale & Querry/Segelstein above ~900 nm (n 1.255 vs 1.317 at 1.6 µm), and
  `n_water = nothing` always uses the 550 nm value instead of a spectral
  lookup; a vector `n_water` uses its first element. Glint in the SWIR is
  ~6 % high. Choose the reference data and whether to look up per band or
  per spectral point.
- [ ] **NAI2 size-parameter Jacobians.** `r_min = quantile(dist, 1e-8)` moves
  the radius nodes with μ_logr/σ_logr while the analytic tangent assumes
  fixed nodes; size Jacobians are accurate only when `nquad_radius` is
  converged. Radiance Jacobians were off by 33–115 % with 100 nodes
  (`r_max = 50`, the JacobianTestFast settings); the optical-property
  derivatives agree to ≈1e-7 with 1600 nodes (`r_max = 5`). Fix the radius
  grid, differentiate it, or document a minimum `nquad_radius`.
- [ ] **Surface-pressure Jacobian with `depol: -1`.** The automatic
  depolarization depends on the column-mean temperature, which moves with
  `p_surf`; that term (~1e-5 relative) has no tangent.
- [ ] **VS `_plus` builder** uses the raw `params.depol` (with `depol: -1`
  τ_Rayleigh ≈ 0.22× the correct value) and mutates `params.spec_bands` and
  `params.brdf` (probable, from reading the code).
- [ ] **Round 7_altgrd pilot** (report in the RRS_XCO2 round7_altgrd output
  directory): no-aerosol scenes match Round 7 (mean ΔXCO2 +0.001 ppm), aerosol
  scenes are 0.49 ppm lower on average (up to 1.52 ppm); five desert solves
  hit the iteration limit. Decide: (1) the aerosol rule (altitude-defined
  aerosols fixed in altitude), (2) the comparability criterion and which
  rounds/members to re-run, (3) the iteration limit, (4) reconcile the note
  that the release "adds the bottom-layer dσ/dp term" (among Rounds 4–7, only
  Round 4 lacked it; Round 3 SIF-on used `cb90b83b`, also without it, and
  earlier unpinned runs have uncertain provenance).
