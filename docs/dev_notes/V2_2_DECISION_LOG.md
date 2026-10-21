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
| `UsingPolTen()` | polarizability tensor (mean polarizability ᾱ, anisotropy γ) | N₂, O₂ (dispersive by default, decision 13c), H₂ (Raman-module constants) |

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
Rayleigh Greek coefficients — comes from the N₂/O₂ polarizability tensors
selected by `params.rayleigh_anisotropy` (default `DispersiveAnisotropy()`,
decision 13c; `_auto_molecular_depolarizations` in `model_from_parameters.jl`).
ϖ_Cabannes, the Cabannes depolarization and the RRS/VRS paths use the
Raman-module constants (`getRamanAtmoConstants`).

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

### 1.8 Decisions of 2026-10-19 (open items of §3 and §3.1)

Taken one by one by Suniti; implementation status is tracked in §3/§3.1.

| # | Item | Decision |
|---|---|---|
| 1 | ABSCO 4-D H₂O-broadener axis in the q Jacobian | Implement the exact ∂σ/∂VMR_H₂O term (linear in the tabulated VMR). |
| 2 | Archived RRS_XCO2 products (dry ABSCO LUTs) | The 4-D tables apply to new campaigns; archives stay, labelled "ABSCO v5.2, dry broadener". |
| 3a | Aerosols under `FixedAltitudeGrid` | Altitude-defined aerosols stay fixed in altitude (as in the pilot). |
| 3b | Re-running earlier rounds | Re-run all rounds on v2.2, with 4-D ABSCO and truth maps regenerated with the 4-D tables. |
| 3c | Desert solves at the 7-iteration limit | Raise the limit for the reruns and repeat the five solves. |
| 3d | "Release adds dσ/dp" wording | Correct it: dσ/dp is in the core since `36e7d1e0`; it is confirmed absent in Round 4 and Round 3 SIF-on, and very likely absent in the unpinned earlier runs (done; the claim was only in the working notes). |
| 4 | R1: RRS band without rotational lines | Raise an `ArgumentError` (done). |
| 5 | R2: VS Q-branch Placzek–Teller factor | After the reference check: Placzek–Teller weight b_JJ to first power (applied; `raman_q_branch_placzek_teller.md`). |
| 6 | R3: VS elastic phase matrix vs SSA | Make them consistent (Cabannes in both). Not implemented: no VS mode redistributes rotational Raman, so a Cabannes-only SSA would lose its energy; options in §3. |
| 7 | R4–R8 | Fix all as cleanup (done). |
| 8 | VS `_plus` builder | Use the resolved depolarization; do not mutate `params` (done; follows `rayleigh_anisotropy` like the standard builder). |
| 9 | Cox–Munk water refractive index | Segelstein (1981), interpolated per spectral point; a vector `n_water` is used element-wise (done). |
| 10 | Ross–Li 80° zenith cap | Make it a parameter, default 80°; off reproduces MODIS/VLIDORT (done: `zenith_cap_deg`). |
| 11 | NAI2 size-parameter Jacobians | Fixed radius grid (nodes independent of the distribution parameters) (done: `scattering.r_min`, default 1e-3 µm). |
| 12 | p_surf tangent of `depol: -1` | Document as neglected (~1e-5 relative) (done, `jacobians.md`). |
| 13a | O₂ King factor | Keep 1.09 (He et al. 2021) (done; cited, with the 1.096 alternative). |
| 13b | Refractive-index fits beyond 468/546 nm | Switch to fits valid through the NIR where the current ones are extrapolated (done: Peck & Khanna 1966 for N₂, Egan 2026 for O₂, continuous joins). |
| 13c | Anisotropy γ | Add the wavelength dependence (N₂, O₂) (done: `DispersiveAnisotropy()`, default for `depol: -1` and `UsingPolTen()`; Raman side unchanged, see §3). |
| 13d | `UsingPolTen` default T | Keep 288.15 K (done; reasons documented). |
| 13e | H₂ refractive-index check | Add a test (done: Peck & Huang 1977 and Raj et al. 2018). |
| 14a | Solar F₀ precedence | Keep kwarg > preset > `model.sources`; a `sources=` call no longer overwrites the preset (done). |
| 14b | `rt_run_ss` `hem_T` | Keep diffuse-only (documented). |
| 14c | Interior upwelling glint | Keep the Fourier-truncated glint; document that the exact glint is TOA-only (documented). |
| 14d | Canopy BHR-down | Keep (no direct beam); document (documented). |
| 14e | 1/(μᵢ − μ₀) cancellation | Implement a limit form in the elemental kernels (done: `expdiff_ratio`). |
| 15 | Sign of ε in B_l | Docs follow the papers' printed form (Sanghavi & Stephens 2015, Eq. 22; Siewert's form, used by the vSmartMOM paper): −ε at [3,4], with the mapping to the code's stored γ, ε; no numerical change. |
| 16 | Q-sign wording in `conventions.md` §1 | Confirmed. |

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

Provenance (`ROUND7_ALTGRD_PILOT_REPORT.md`). Confirmed from the outputs:
Rounds 5–7 ran on `7acab570`, which includes the term. Round 4 (`cb90b83b`,
`630c9190`) and Round 3 SIF-on (`cb90b83b`) ran from `suniti_multi_sensor`
commits that never contained it, whatever their output dates. Inferred:
Rounds 1–2 and Round 3 no-SIF record no commit (the reflog suggests Round 3
no-SIF ran at `cb90b83b`); they finished before `36e7d1e0` existed, so the
term is very unlikely there, and only uncommitted local code could change
that.

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

- [x] **Molecular Rayleigh port review** (decisions 13a–13e): done; formulas,
  sources, joins and before/after numbers in
  `docs/dev_notes/molecular_rayleigh_dispersion.md`. The golden test in
  `test/test_vertical_grid_psurf.jl` and the Phase 1b RRS reference test
  select `ConstantAnisotropy()`, the formulation their stored values were
  made with; the stored files are unchanged.
- [ ] **For Suniti (follow-ups of decision 13):**
  (1) Raman side: with the dispersive Rayleigh default, an elastic-only run
  and a Raman run no longer share one elastic + RRS decomposition (implied
  depolarization differs by up to ≈10 % in the NIR). Options: keep the Raman
  constants and recommend `ConstantAnisotropy()` for paired noRS/RRS runs;
  switch the Raman cross sections too (σ_RRS −6.6 % at 760 nm, +6.3 % at
  395 nm); or make the Rayleigh default constant again.
  (2) The dispersive O₂ anisotropy follows the 1.09 King factor; Hohm (1994)
  supports 1.096 (γ 3 % higher). (3) Bodhaine Eq. 30 τ_Rayleigh is +1.0 %
  high at 1.6 µm and +3.6 % at 2.06 µm against Peck–Reeder/Ciddor
  refractivity, and the King rescaling uses a constant ρ₀ = 0.0279 (τ 0.3 %
  high at 400 nm); option: an RI-based air σ for τ_Rayleigh. (4) The
  Raman-module H₂ anisotropy is 8.5 % below the ab initio value (σ_RRS(H₂)
  ≈16 % low). (5) N₂ join kept at 436.3 nm (continuous) rather than 468 nm
  (2.85×10⁻⁴ jump in σ).
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
- [x] **Ross–Li zenith cap (decision 10).** `RossLiSurfaceScalar(...;
  zenith_cap_deg = 80)`, `nothing` gives the MODIS/VLIDORT kernels. VLIDORT
  Case D also checks the uncapped references (≤ 7e-5 in I); the cap changes
  I by 0.34–0.81 % at TOA and 0.39–1.07 % at BOA in that comparison.
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
- [x] **ABSCO H₂O-broadener tangent (decision 1, §1.8).** The q-driven H₂O
  columns include `N_dry·VMR·∂σ/∂b·db/dr` for every ABSCO table with a
  multi-node broadener axis (`_accumulate_broadener_tangent!`) under both
  `abundance_jacobians` settings, so `:complete` accepts these models and
  `:frozen_cross_section` changes no values. Checked against central differences on the real v5.2 tables
  (≤ 4e-8 relative). Note: the forward call passes `r = N_H2O/N_dry` to the
  O₂/CO₂ tables but the moist mole fraction `r/(1+r)` to the H₂O table; which
  one the ABSCO `Broadener_01_VMR` axis means is still to be confirmed.
- [ ] **4-D ABSCO in RRS_XCO2 (decisions 2 and 3b; scope decided by Suniti
  2026-10-19).** The workflow default stays the dry `.jld2` tables, so
  existing and archived products never mix with 4-D ones. The 4-D tables enter
  through the v2.2 rerun campaign: a fresh data root with its own pinned
  environment that fixes and hashes the five table paths, as
  `round7_altgrd_environment.sh` does. A generic per-data-root marker and
  provenance system was drafted and reviewed three times, but left out because
  each review found new mixing paths (overrides, resumes, concurrency,
  launchers); the draft is saved in
  `vsmartmom-release-evidence/2026-10-19/absco_round3_full_WIP.patch`.

### 3.1 Added by the 2026-10-17 bug audit (results-changing; need Suniti)

The audit findings that are code defects with a clear reference were fixed
in reviewed commits (see `CHANGELOG.md`). The items below change scientific
results or definitions and wait for a decision.

- [x] **Solver-entry conventions (decisions 14a–14e).** (1) Solar F₀
  precedence stays `sources=` > user-preset `RS_type.F₀` > `model.sources`;
  a `sources=` keyword now applies to its call only (14a, done:
  `_with_call_scoped_F₀` around `rt_run`, `rt_run_ss`, `rt_run_toa(rs::RRS)`
  and the linearized `rt_run`). (2) `rt_run_ss` `hem_T` stays diffuse-only
  (14b, documented). (3) Interior upwelling levels keep only the
  Fourier-truncated glint (14c, documented: the exact glint is TOA-only).
  (4) Canopy BHR-down omits the direct beam (14d, documented). (5) The
  elemental factors are evaluated in limit form, `expdiff_ratio` /
  `elem_trans_factor[_dτ]`, in every live forward, linearized, canopy and
  Raman kernel (14e, done; `StandaloneSS._τ_integral` is not an elemental
  kernel, see the follow-ups). The 4·eps merge window is
  kept for exact on-node results but is no longer needed for accuracy, and
  Radau's exact-equality node merge is sufficient (nadir view at SZA 1e-5°:
  2.2e-6 → 8.8e-12).
- [ ] **For Suniti (solver follow-ups of decision 14):**
  (1) `StandaloneSS._τ_integral` (`src/StandaloneSS/kernels.jl:62-71`) has
  the same cancellation with a 1e-10 trapezoid fallback (a Float32 case
  returns 0 instead of 0.0082): use Δτ·`expdiff_ratio` in a follow-up
  (recommended) or leave it?
  (2) Removing the VS_plus window also replaced an assignment (`=`) by an
  accumulation (`+=`) in that branch; the arrays are zeroed first and each
  entry is written once per launch, so results do not change (note only).
  (3) The VS (non-plus) and canopy kernels set t⁺⁺[i,j] = 0 for i ≠ j at
  equal μ (Stokes cross terms), while the elastic, RRS and VS_plus kernels
  use the limit; pre-existing, left untouched. Separate review
  (recommended) or leave?
  (4) Metal (no Float64) is not verified with the new series helpers: one
  Metal smoke test before release (recommended)?
- [x] **Sign of ε in the Greek matrix B_l (decision 15).** `conventions.md`
  §2 and `concepts/03b_scattering.md` show B_l and Π_l^m as printed by
  Sanghavi & Stephens (2015, Eqs. 21–22), in Siewert's form, which the
  vSmartMOM paper uses by reference: −ε at [3,4]. Sanghavi (2014) Eq. 16 was
  not checked (no copy available). Relative to the Siewert/VLIDORT convention
  the code's Stokes vector is (I, −Q, −U, −V). The code therefore stores γ and
  ε with the opposite sign and writes +ε at [3,4], which is numerically
  Siewert's −ε. Coefficients from the code's own Mie code need no flips.
  A sign-flip experiment on Siewert's Problem IIA confirms this: the papers'
  B_l with raw Siewert coefficients reproduces VLIDORT 2.8.3 without flips,
  and negating ε alone reverses V only. (Evidence:
  `vsmartmom-release-evidence/2026-10-19/greek_sign_flip_experiment.md`.)
  `test/vlidort_baseline/extract_references.jl` emitted `SIEWERT_γ` without
  the import flip that the committed `siewert2000_IIA_greek.jl` carries (the
  file had been corrected by hand). The generator now emits that file byte
  for byte, so a re-extraction keeps the flip.

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
- [x] **Raman** (decisions 4, 5, 7; decision 6 open below, decision 8 has its
  own item). R2 done 2026-10-19 (Q-branch weight b_JJ,
  `raman_q_branch_placzek_teller.md`). R1 (`ArgumentError` for a band without
  rotational lines), R4 (VS Q-branch phase from ⟨b_JJ⟩), R5 (`VS_0to1_plus`
  source depletion), R6/R8 (dead overloads removed) and R7 (linear
  interpolation from the reference grid point) done 2026-10-20. The Phase 1b
  RRS reference is corrected in-test for the line placement (data
  unchanged).
- [ ] **For Suniti (Raman follow-ups):**
  (1) R3 is not implemented: no VS mode carries rotational Raman
  inelastically (VS_plus weights are VRS/RVRS only; `VS_0to1`/`VS_1to0` cannot
  run), so a Cabannes-only SSA would lose 1 − ϖ_Cab = 3.2 % at 395 nm.
  Options, with effects vs noRS measured for `VS_0to1_plus` at 395 nm (pure
  Rayleigh, one layer, SZA 45°, flat F₀; TOA at nadir, VZA 45° Δφ 0° and
  180°, VZA 60° Δφ 90°; constant anisotropy, i.e. before `f8211b00`, which
  moves the noRS phase by about 0.2 % in β₂; scripts and logs in
  `vsmartmom-release-evidence/2026-10-20/raman_vs_r3/`):
  (a) Rayleigh phase with the Rayleigh-minus-VRS SSA (closure ≤ 4.9e-4
  relative in TOA, ≤ 1.4e-4 in BOA; recommended), (b) Cabannes-only SSA as
  decided (I −3.7 % to −6.4 % relative), (c) a combined RRS+VS mode (new
  development), (d) keep (DOLP 3.6–5.3 % relative too high over the four
  directions; 0.2776 vs 0.2681 at nadir).
  (2) RRS bands that hold only some lines: energy is conserved by scaling the
  captured lines (×13.4 in Phase 1b, 7.4 % of σ_RRS), which makes the nearest
  ghosts 13× too strong. Keep and document (done in `08_inelastic.md`), warn
  below ~95 % capture (recommended), error, or extend the grid internally?
  (3) Stellar H₂ paths (`sol_RRS`, `sol_VS_plus`; not run by CoreRT today)
  have the R6/R7/R4/R8 analogues (Stokes/anti-Stokes swap, off-grid ν₀ with
  one array reused for both line sets, whole-band ρ, dead method): apply the
  same fixes in a follow-up?
  (4) Accept the in-test Phase 1b correction, or regenerate the reference?
  (5) R7 changes RRS_XCO2 inelastic radiances in structured spectra: check one
  truth-map member on v2.2 before mixing new and archived products.
- [x] **Water refractive index for Cox–Munk (decision 9).** Segelstein (1981)
  pure-water table, n and k linear in wavelength at every spectral point; a
  vector `n_water` element-wise; the per-point operator from a few knot
  indices interpolated in n (≤ 1e-7 in Float64 for water at 17 bands from
  0.25 to 3.9 µm, ≤ 1e-6 in Float32 at 0.55, 1.6, 2.06 and 2.8 µm); exact TOA
  glint per point.
  Glint TOA I (U = 5 m/s, SZA 30°): +3.6 % at 0.55 µm, −9.3 % at 1.6 µm,
  −16.9 % at 2.06 µm (the audit's −6 % estimate at 1.6 µm used Hale &
  Querry's n = 1.317; Segelstein gives 1.3096).
- [ ] **For Suniti (surface follow-ups of decisions 9, 10, 14d):**
  (1) Segelstein is pure water at 25 °C, with no salinity: keep it
  (recommended) or add a seawater correction later? (For reference, Hale &
  Querry's pure-water n is 1.317 at 1.6 µm, Segelstein's 1.3096.)
  (2) Knot interpolation in n as the reading of "per spectral point"
  (≤ 1e-7 in Float64, ≤ 1e-6 in Float32 at the tested bands; exact
  evaluation at every point
  costs about nSpec/knots times the surface work): keep (recommended)?
  (3) A standalone `CoxMunkSSSurface(n_water = nothing)` built without a model
  keeps its documented 1.34 (the model adapter passes the per-point index):
  keep documented (recommended) or throw?
  (4) Canopy: the forward BOA radiance at the solar ordinate has no direct
  beam (j₀⁺ = 0) while the linearized solver adds it; documented per 14d.
  Add the direct beam to the forward solver in a follow-up (recommended)?
  (5) Single-wavelength Cox–Munk helpers now throw for `n_water = nothing` or
  a vector with differing entries: keep (recommended) or return to a
  documented 550 nm fallback?
- [ ] Bug (pre-existing, found during decision 9): `rt_run(model; i_band = 4)`
  fails on a model with more than three bands, because the default `noRS`
  `ϖ_Cabannes` has three entries (`Inelastic/types.jl:191`, indexed per band
  in `compEffectiveLayerProperties.jl`).
- [x] **NAI2 size-parameter Jacobians (decision 11).** Fixed `[r_min, r_max]`
  radius grid, `scattering.r_min` default 1e-3 µm, also for PCW. Size
  tangents equal central differences to ~1e-10 (optics); TOA Jacobians to
  3e-8 (audit mode on JacobianTestFast settings) and ≤ 6e-6 (OCO-like scene),
  metric max |analytic − FD| / max |FD|. Before: 7–18 % with the audit's
  0.2 µm, σ_g = 1.6 mode and 26 % for the OCO organic mode at 300 nodes.
  Accuracy against the converged integral still depends on `nquad_radius` and
  needs a convergence check (the TMS ordering test moved to 400 nodes for its
  1 µm mode). The PCW reference test and the 7acab570 golden test pin
  `r_min` to the old quantile.
- [ ] **For Suniti (radius-grid follow-ups of decision 11):**
  (1) Default `r_min`: 1e-3 µm (recommended; truncates nothing physical) or
  5e-3 µm (node density closer to the old grids)?
  (2) A legacy quantile option for exact reproduction of archived multi-mode
  runs (RRS_XCO2 moves ≤ 1e-5 in TOA): no (recommended)?
  (3) Coarse modes at about 100 nodes: raise `nquad_radius` in the
  JacobianTestFast-style configs, or document only (done, recommended)?
  (4) Downstream: `uni_vSmartMOM_v2-2/.../generate_truth_map.jl:219` rebuilds
  `ScatteringParameters` with 6 arguments and would drop a configured
  `r_min` (none is set today); pass `FT(sp.r_min)` when that workflow is next
  edited.
- [x] **Surface-pressure Jacobian with `depol: -1` (decision 12).** The
  automatic depolarization (Rayleigh Greek coefficients and, through the King
  factor, τ_Rayleigh) depends on the column-mean temperature, which moves
  with `p_surf` under `FixedPressureGrid()` (not under `FixedAltitudeGrid()`);
  that term (~1e-5 relative) has no tangent. Documented as neglected in
  `jacobians.md`.
- [x] **VS `_plus` builder** (decision 8): resolves the depolarization as the
  standard builder does (including `rayleigh_anisotropy`), works on a copy of
  `params`; line-by-line absorption and Float32 builds work.
- [ ] **Round 7_altgrd pilot** (report in the RRS_XCO2 round7_altgrd output
  directory): no-aerosol scenes match Round 7 (mean ΔXCO2 +0.001 ppm), aerosol
  scenes are 0.49 ppm lower on average (up to 1.52 ppm); five desert solves
  hit the iteration limit. Decide: (1) the aerosol rule (altitude-defined
  aerosols fixed in altitude), (2) the comparability criterion and which
  rounds/members to re-run, (3) the iteration limit, (4) reconcile the note
  that the release "adds the bottom-layer dσ/dp term" (among Rounds 4–7, only
  Round 4 lacked it; Round 3 SIF-on used `cb90b83b`, also without it, and
  earlier unpinned runs have uncertain provenance). Decided 2026-10-19
  (decisions 3a–3d, §1.8); (4) is done. (1)–(3) are carried out by the v2.2
  rerun campaign.
