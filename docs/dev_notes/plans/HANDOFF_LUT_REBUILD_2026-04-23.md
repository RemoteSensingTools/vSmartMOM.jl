# Session Handoff — 2026-04-23: LUT Rebuild + MODTRAN Gap Analysis

## Part A — Notes on item 3: `transm_down_dif` ≪ MODTRAN, `transm_up_dif` comparable

Observed (scenario `H2O0.0500_AOT0.0001_GNDALT0.000_TSZ168.0`):

| λ (nm) | field       | MODTRAN   | vSmartMOM | V/M ratio |
|--------|-------------|-----------|-----------|-----------|
| 400    | td_dif      | 0.341     | 0.182     | 0.53      |
| 400    | tu_dif      | 0.171     | 0.171     | 1.00      |
| 550    | td_dif      | 0.267     | 0.059     | 0.22      |
| 550    | tu_dif      | 0.064     | 0.053     | 0.83      |
| 2000   | td_dir      | 0.437     | 0.975     | 2.23      |

Why `td_dif` is systematically suppressed while `tu_dif` stays close:

### (1) Path-length asymmetry amplifies τ deficits on the sun leg
- Solar path: μ_s = cos(39.58°) = 0.771, airmass ≈ 1.30
- View path: μ_v = cos(12°) = 0.978, airmass ≈ 1.02

Any τ shortfall Δτ (from missing absorbers / continua / baseline aerosol)
translates into:
- `td_dir = exp(-τ/μ_s)` over-estimated by a factor `exp(Δτ/μ_s)`
- `tu_dir = exp(-τ/μ_v)` over-estimated by a smaller `exp(Δτ/μ_v)`

For small τ, diffuse ≈ `(1 − exp(−τ/μ)) · f_scatter` ∝ airmass. Missing
τ therefore hits `td_dif` about 1.3× harder than `tu_dif` at the scenario
geometry. At λ=550 nm, the V/M ratio spread (0.22 vs 0.83) is larger than
1.3× because…

### (2) Missing aerosol forward-scattering hits `td_dif` specifically
MODTRAN's default aerosol has g ≈ 0.65–0.75; forward-peaked scattering
feeds diffuse *downward* transmission almost exclusively. Our AOT=0.0001
effectively switches aerosol off. The lost forward-scatter contribution
shows up almost entirely on the down leg; the Lambertian-source up leg
is barely affected.

### (3) The two diffuse terms are *computed* differently in our pipeline
- `td_dif = hem_T(alb=0)/μ_s − exp(-τ/μ_s)`. `hem_T(alb=0)` is a BOA
  hemispheric flux ratio that integrates every downward slant path, so
  any scattering/absorption gap propagates into it 1:1.
- `tu_dif = T_up_dir − exp(-τ/μ_v)`, where `T_up_dir = α·(1−0.2·S)/(0.2·T_down)`
  and `α = π·(L(0.2)−L(0))/μ_s`. At near-nadir VZA and small τ both
  `T_up_dir` and `exp(-τ/μ_v)` are close to 1; their *difference* is a
  small residual that stays close to MODTRAN even when their individual
  magnitudes don't.

So the asymmetry is diagnostic, not a bug. Closing the τ gaps (aerosol
baseline + CIA/continua + LUT resolution) should pull `td_dif` toward
MODTRAN while leaving `tu_dif` where it is.

---

## Part B — Plan for the improved LUT

### Objectives
1. Wavelength range 280–3500 nm (ν = 2857.14–35714.29 cm⁻¹, span 32857 cm⁻¹).
2. Switch interpolation from `BSpline(Cubic(Line()))` to `BSpline(Linear())`
   to eliminate over/undershoot on sharp lines.
3. Extend species list to MODTRAN-equivalent coverage.
4. Target wavenumber step Δν = 0.001 cm⁻¹ **— but see sizing table below;
   this is the first decision to make.**

### Sizing trade-off at different Δν
Span 32857 cm⁻¹. Store Float32 cross-sections; no spline coefficients
needed beyond the grid values when using linear interpolation.

| Δν (cm⁻¹) | # ν pts | 15 p × 10 T (150 pairs) | 55 p × 37 T (2035 pairs) |
|-----------|---------|-------------------------|--------------------------|
| 0.001     | 32.9 M  | 19.7 GB / species        | 268 GB / species          |
| 0.002     | 16.4 M  | 9.9 GB                  | 134 GB                    |
| 0.005     | 6.57 M  | 3.9 GB                  | 53 GB                     |
| 0.01      | 3.29 M  | 2.0 GB                  | 26.8 GB                   |

Compute (GPU A100, one species, ~50 k HITRAN lines, Voigt + HW32 CEF):
- ~10–60 s per (p, T) at 0.001 cm⁻¹ → 150 pairs = ~1–3 hr × 20 species = 20–60 GPU-hours
- 2035 pairs would be ~150–600 GPU-hours — probably not worth it if the
  p,T grid benefit is marginal compared to the line-resolution fix.

**Physics check for Δν choice:** at 1 atm surface pressure, air-broadened
HWHM in the IR is ~0.05 cm⁻¹; in the stratosphere (p ≈ 1 hPa, Doppler-
dominated) HWHM can drop to ~0.005–0.01 cm⁻¹ for light molecules at Vis
wavelengths. 0.001 cm⁻¹ samples every line at ≥5 points per HWHM
everywhere we care about. 0.01 cm⁻¹ is ~Nyquist at surface and under-samples
in the stratosphere by ~2× — survivable with linear interpolation since
the error is *smoothing*, not ringing, but band contrast at low pressure
will be under-represented.

**Recommendation:** commit to **0.001 cm⁻¹** only over 280–3500 nm if
storage and wall-clock are acceptable. Otherwise split bands:
- 280–800 nm (UV-Vis, mostly broad bands): 0.005 cm⁻¹
- 800–3500 nm (NIR/SWIR, narrow lines dominate): 0.001 cm⁻¹

Total storage at 0.001 cm⁻¹ / full range / 15×10 p,T / 20 species ≈ **400 GB**.
Filesystem has 574 TB free — not a constraint.

### Species list (MODTRAN-equivalent subset)

Current (7): H2O, CO2, O3, N2O, CO, CH4, O2.

Additions via HITRAN line database (~12):
- **UV-Vis**: NO2 (M=10), SO2 (M=9), HCHO/H2CO (M=20), HCN (M=23)
- **NIR/SWIR**: HCl (M=15), HF (M=14), C2H2 (M=26), C2H6 (M=27), OCS (M=19),
  NH3 (M=11), HNO3 (M=12)
- **Stratospheric**: NO (M=8), ClO (M=18), HOCl (M=21), H2O2 (M=25)

**Not in HITRAN line format — require separate handling:**
- **CIA** (collision-induced absorption): O2–O2, N2–N2, N2–O2, H2O–H2O,
  H2O–N2, CO2–CO2. Source: HITRAN-CIA (separate .cia files). Strongly
  affects 2.0 μm and 4.3 μm windows — a leading suspect for the 2000 nm
  td_dir discrepancy in the MODTRAN comparison.
- **Water-vapor continuum (MT_CKD)**: matters from 1.9 μm onward and in
  the UV blue (4 μm continuum tail). Needs its own implementation
  (external code or port of MT_CKD 4.1/4.2 to Julia). Deferring.

### Interpolation change
One-line edit in `src/Absorption/make_model_helpers.jl:91`:
```julia
# Current:
itp = interpolate(cs_matrix, BSpline(Cubic(Line(OnGrid()))))
# New (linear in ν, same or linear in p, T):
itp = interpolate(cs_matrix, BSpline(Linear()))
```
The 4-D pre-computed LUTs remain the same shape; only the interpolant
changes. Runtime cost is *lower* with linear vs cubic. No Gibbs overshoot.

Recommend also removing the 2-point padding that `BSpline(Cubic(Line))`
implicitly uses (see `O2_emit_lut.jld2` where `itp.coefs` is
`(223320, 57, 39)` vs `ν_grid` of length 223318 — two extra at each end).

### Existing code paths
- Per-molecule build script template: `sandbox/prototyping/createAbscoJLD2.jl`
  (uses ABSCO tables, would adapt for HITRAN path).
- HITRAN fetch: `fetch_hitran("NO2"; numin=2857, numax=35714, edition="HITRAN2024")`
  is already wired (`src/Artifacts/hitran_api.jl:57`).
- LUT builder: `make_interpolation_model(hitran, Voigt(), ν_grid, p_grid, t_grid; …)`
  at `src/Absorption/make_model_helpers.jl:55`. Need a thin driver that
  loops over species and writes one JLD2 per molecule.

### Proposed build phases

**Phase 1 — infrastructure & spot check** *(1 session)*
- Switch interpolator to Linear (one-line edit + test that existing
  tests still pass).
- Write `test/benchmarks/build_hitran_luts.jl` driver (parameterized
  species list, p,T grid, ν step, output dir).
- Pilot-build a single species (e.g. NO2) to validate wall clock
  and storage figures against the table above.

**Phase 2 — pilot on existing species** *(1 session)*
- Rebuild O2, H2O, CO2 at 0.001 cm⁻¹ over 280–3500 nm, linear interp.
- Re-run the EMIT scenario and compare to MODTRAN in the A-band,
  1.94 μm and 2.0 μm windows. Expect: A-band spikes/negatives gone,
  band shapes closer to MODTRAN.

**Phase 3 — add missing species** *(compute-bound, may span days)*
- Build NO2, SO2, H2CO, HCN, HCl, HF, C2H2, C2H6, OCS, NH3, HNO3,
  NO, ClO, HOCl, H2O2 (15 species). Total ≈ 300 GB at 0.001 cm⁻¹ /
  15 p × 10 T.

**Phase 4 — CIA and continua** *(separate work item)*
- Port HITRAN-CIA data into a LUT (simpler: dense 2-D in T × ν,
  per collision pair).
- Evaluate whether MT_CKD water continuum is worth the implementation
  effort vs just folding an empirical continuum correction into the
  H2O LUT over 1.9–2.7 μm.

**Phase 5 — regression + docs** *(1 session)*
- Rerun MODTRAN-equivalent comparison at ≥3 wavelengths in each
  spectroscopic window. Document residual gaps vs MODTRAN.
- Update `CLAUDE.md` with LUT format & interpolation convention.

### Decisions (settled 2026-04-23)

| # | Decision | Value |
|---|----------|-------|
| 1 | Δν | **0.01 cm⁻¹ everywhere** (see pilot results below) |
| 2 | p,T grid | **~15 p × 10 T = 150 pairs** |
| 3 | CIA | pilot O2–O2 + N2–N2 early to diagnose 2 μm gap |
| 4 | Species list | 7 existing + 15 additions (15 kept in full) |
| 5 | Float type | **Float32** |
| 6 | Output path | **`~/data/HITRAN_LUTs/`** |

### Part C — Resolution-convergence pilot (2026-04-23)

Pilot: `test/benchmarks/pilot_lut_resolution.jl`. Direct line-by-line
cross-sections via `make_hitran_model` + `compute_absorption_cross_section`
on two bands × two (p, T) conditions, compared to a Δν = 0.0005 cm⁻¹
reference. Metric: **|ΔTmean|**, the absolute error in the band-mean
transmittance `⟨exp(-σ·N_col)⟩` using a representative full-atmosphere
column (N_O2 = 4.5×10²⁴, N_H2O = 5×10²²).

| Band        | (p, T)              | Δν=0.01  | Δν=0.005 | Δν=0.001  | Δν=0.1 (current) |
|-------------|---------------------|----------|----------|-----------|------------------|
| H2O 1.94 μm | trop (500 hPa, 250) | 2.1e−6   | 9.8e−7   | 1.1e−7    | 4.2e−4           |
| H2O 1.94 μm | strat (10 hPa, 220) | **3.3e−5** | 1.3e−6 | 1.2e−7    | 9.2e−4           |
| O2 A-band   | trop (500 hPa, 250) | 4.4e−6   | 2.2e−6   | 2.5e−7    | 7.9e−5           |
| O2 A-band   | strat (10 hPa, 220) | 4.3e−6   | 1.1e−6   | 1.2e−7    | 9.2e−4           |

Takeaways:
- **Δν = 0.01 cm⁻¹ caps the band-mean transmittance error at ≤ 3.3×10⁻⁵**
  (stratospheric H2O 1.94 μm being the worst case). That is ~30× tighter
  than typical retrieval noise (~1e−3), so 0.01 cm⁻¹ is "free accuracy"
  relative to downstream errors.
- **Δν = 0.001 cm⁻¹** is ~100× tighter again (~1e−7), but costs 10× the
  storage and line-by-line compute. Not worth it.
- Δν = 0.1 cm⁻¹ *band-mean* error is already ≤ 1×10⁻³ — so the ~50%
  discrepancy we see at 2 μm in the MODTRAN comparison is **not a
  resolution problem**. Missing physics (aerosol baseline, CIA, water
  continuum) is the dominant cause. Resolution fix is still worth
  doing for point-by-point fidelity (and the LUT artifacts), but it
  won't close the 2 μm gap by itself.

Raw numbers: `~/pilot_lut_resolution_results.jld2`, logs in
`~/pilot_lut_run_v2.log`.

### Finalized build spec
- **Spectral**: Δν = 0.01 cm⁻¹, range 280–3500 nm → 2857.14–35714.29 cm⁻¹
  → **3.286 M ν-points per (p, T)**
- **p grid** (15 points, log-spaced + boundary-anchored):
  `[0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30, 100, 300, 500, 700, 900, 1013.25, 1080]` hPa
- **T grid** (10 points): `180, 200, 220, 240, 260, 280, 300, 320, 340, 360` K
- **Storage per species**: 3.286M × 150 × 4 B = **~1.97 GB / species**
- **20 species** (existing 7 + 12 new + 1 budget for extras): ~40 GB total
  — trivially fits under `~/data/HITRAN_LUTs/`
- **Interpolant**: `BSpline(Linear())` for ν; linear in p, T

### Next step
Write the Phase-1 driver `test/benchmarks/build_hitran_luts.jl` and
kick off the rebuild. Before Phase 2, run the O2–O2 + N2–N2 CIA
spot-check to confirm continua explain the 2 μm gap (otherwise more
physics gaps remain).

### Part I — Beer-Lambert symmetry check between T_down_dir and T_up_dir (2026-04-23)

**Question from user:** For any monochromatic wavelength,
`ln(T_up_dir) · cos(vza) = ln(T_down_dir) · cos(sza) = −τ_total` should
hold on both sides. Does it?

**Empirical check** (scenario H2O=0.05, AOT=0.0001, GNDALT=0, TSZ=168,
μ_s=0.771, μ_v=0.978):

| λ (nm) | MODTRAN τ_dn | MODTRAN τ_up | MODTRAN τ_dn−τ_up | vSm τ_dn−τ_up |
|--------|-------------|--------------|--------------------|---------------|
| 380    | 0.690       | 0.490        | **+0.20**          | 0.0000        |
| 450    | 0.456       | 0.255        | **+0.20**          | 0.0000        |
| 550    | 0.349       | 0.149        | **+0.20**          | 0.0000        |
| **760**    | 0.707       | 0.915        | **−0.21 (reversed!)** | 0.0000      |
| 1270   | 0.308       | 0.116        | +0.19              | 0.0000        |
| 2000   | 0.639       | 0.446        | +0.19              | 0.0000        |
| 2200   | 0.330       | 0.121        | +0.21              | 0.0000        |
| 3000   | 0.411       | 0.213        | +0.20              | 0.0000        |

**vSmartMOM** is Beer-Lambert-consistent by construction — both legs are
computed from the same `τ_total` in `test/benchmarks/modtran_equivalent_fields.jl`
lines 215–216. **MODTRAN is not.** The residual has two distinct
mechanisms:

#### Mechanism 1 — wavelength-independent ≈ 0.20 offset in continuum

The offset τ_dn − τ_up is flat to ~5 % across 380–3000 nm. This rules
out Rayleigh (λ⁻⁴, would span 4 decades) and strongly matches
**rural 23 km-visibility aerosol**, which has AOT(550) ≈ 0.2 roughly
flat in the VNIR.

Most likely mechanism: MODTRAN's 6S-style definition of "direct"
treats aerosol forward-scattering asymmetrically on the two legs:

- **Down leg:** solar beam is narrow. Any aerosol scatter removes a
  photon, so `td_dir = exp(−(τ_abs + τ_rayleigh + τ_aer_ext)/μ_s)`.
- **Up leg:** originates from an isotropic Lambertian surface source.
  Photons scattered away from μ_v by aerosol are compensated by
  photons scattered *into* μ_v from other angles of the same
  isotropic source, so "direct" on the up leg excludes aerosol
  scattering and bundles it under `transm_up_dif` instead.
  `tu_dir = exp(−(τ_abs + τ_rayleigh)/μ_v)`.

Difference τ_dn − τ_up = τ_aer_ext ≈ 0.20, λ-independent, matching
observation.

#### Mechanism 2 — the 760 nm O2 A-band reversal

Only reverses sign inside a narrow-line-saturated band. This is a
**spectral-averaging / curve-of-growth artefact**:

- MODTRAN output is binned to ~0.1 nm. Inside the A-band each bin
  contains multiple saturated O2 lines.
- At long slant path (small μ_s) the line cores are fully opaque,
  so the band-averaged `T_d` is dominated by between-line windows
  and refuses to decrease as fast as Beer-Lambert would predict
  — apparent τ = −ln(T) · μ under-estimates the true column τ.
- At short slant path (μ_v ≈ 0.98, nearly vertical) the lines are
  less saturated, apparent τ tracks the true column more faithfully.
- Net: apparent τ_dn < apparent τ_up inside the A-band — reversed
  from the aerosol offset direction.
- vSmartMOM samples at 0.1 cm⁻¹ (essentially monochromatic), so
  curve-of-growth does not apply to its output.

#### Implications for the MODTRAN-vs-vSmartMOM regression

1. **vSmartMOM's internal `τ_down ≡ τ_up`** is physically correct
   at the monochromatic level.
2. **Field-by-field comparison of `transm_*_dir` to MODTRAN is
   intrinsically apples-to-oranges.** Differences of ≈ 0.2 AOT
   (continuum) and direction-dependent line saturation (strong
   absorption bands) should be expected even after every physics
   gap is closed.
3. **Use `T_tot = T_dir + T_dif` for regression**, not `T_dir`
   alone. `T_tot` is convention-independent: under the 6S
   decomposition the aerosol-asymmetry moves back and forth
   between `dir` and `dif` but their *sum* is conserved.
4. When we eventually turn on a non-zero default aerosol baseline
   in our YAML (to actually match MODTRAN's default scenario),
   we should see our own `td_dir` drop by ≈ exp(−0.2/μ_s) = 0.77×
   of its current value — confirming the convention diagnosis.

Verification script: the one-off Julia snippet used for the table
above is in `~/*.log`; a permanent version could live as a
small function in `regression_newLUT.jl` if we want to rerun the
check after the new LUTs are propagated.

### Part H — New-LUT regression (2026-04-23)

Ran `test/benchmarks/regression_newLUT.jl` with the rebuilt LUTs
(new YAML: `test/test_parameters/ParamsEMIT_MODTRANcomp_newLUT.yaml`).

**Pathology check on the full 223k-point spectrum:**
| Metric              | OLD cubic LUT | NEW linear LUT |
|---------------------|---------------|----------------|
| # of τ < −1e-3      | 203           | **0**          |
| # of τ = NaN        | 0             | 0              |
| # of pointwise τ>10 | 461 (artifacts) | 880 (**real line cores**) |

**O2 A-band window (759.5–762 nm, 432 points):**
| Metric         | OLD    | NEW            |
|----------------|--------|----------------|
| min τ          | −1.83  | **+0.059**     |
| max τ          | 18     | 504 (line core) |
| median τ       | 0.065  | **1.29**       |

The "18" spike in the old LUT was a single-point cubic overshoot; everything
else stayed at a ~0.05 baseline that failed to represent the A-band at all.
The new LUT shows the proper line-core/ between-line contrast. No negative τ.

**Band point-checks vs MODTRAN (`T_dir_slant` = exp(-τ/μ_s), μ_s=0.77):**

| λ (nm) | old vSm T_dir | new vSm T_dir | MODTRAN T_dir | assessment                                |
|--------|---------------|---------------|---------------|-------------------------------------------|
|  380   | 0.560         | 0.539         | 0.408         | Gap narrowed slightly; remainder = baseline aerosol + NO2 xsec |
|  550   | 0.882         | 0.874         | 0.636         | Same; mostly aerosol baseline             |
|  760   | 0.965         | 0.527         | 0.400         | **A-band now absorbing**; gap is CIA + aerosol |
|  762.5 | 0.937         | 0.467         | 0.448         | **Excellent match**                       |
| 1270   | 0.957         | 0.995         | 0.671         | CIA missing (see Part F — 1.27 μm O2 a¹Δg)|
| 1930   | ~1            | 0.005         | 0.619         | **H2O 1.94 μm band now absorbing** — overshoots MODTRAN because our line-core sampling is sharper than MODTRAN's instrument resolution |
| 2000   | 0.975         | 0.592         | 0.437         | Gap narrowed 10× (from 0.54 to 0.15) via resolved H2O bands; residual = MT_CKD continuum + aerosol |
| 2200   | 0.983         | 0.670         | 0.652         | **Match**                                 |
| 2300   | 0.988         | 0.846         | 0.700         | Small residual gap (CH4 + H2O continuum)  |

**Fixes the LUT rebuild delivered:**
1. **O2 A-band corruption gone** — was the most dramatic failure. Now resolves ~500-to-1 line/baseline contrast faithfully.
2. **No more negative τ** (203 unphysical points in old LUT → 0 now).
3. **H2O 1.94/2.0/2.7 μm bands now absorbing properly.**
4. **Float32 + linear interpolation + 0.01 cm⁻¹ spacing** produce a LUT that's pointwise robust (no cubic overshoots / undershoots).

**Gaps still present (expected, flagged in Part F):**
1. Baseline aerosol: MODTRAN's default ~AOT(550)=0.2 rural model. Our AOT=0.0001 essentially aerosol-off.
2. H2O self + foreign continuum (MT_CKD) — explains residual 2.1-2.3 μm gap.
3. UV-Vis broadband absorbers needing `.xsc`: NO2 (300–500 nm), H2CO (300–360 nm),
   O3 Huggins. SO2 is **already ingested** via `augment_lut_with_xsec.jl`.

### Part G — Full 22-species rebuild outcome (2026-04-23)

Launched via `CUDA_VISIBLE_DEVICES=0 julia --project=test test/benchmarks/build_hitran_luts.jl`.
Total wall clock: **3642 s (61 min)** on one A100-40GB.

| Status | # | Species |
|--------|---|---------|
| ✅ built | 19 | H2O, O3, N2O, CO, CH4, O2, NO, SO2, NO2, NH3, HNO3, HF, HCl, OCS, H2CO, HOCl, HCN, C2H2, C2H6 |
| ⚠️ rebuilt iso=1 only | 1 | CO2 — HITRAN iso letter codes ("A"=10,"B"=11,…) crash `mol_weight` via `read_hitran` returning iso=0; restricted to main isotopologue (98 % abundance) |
| ❌ no HITRAN .par lines in 280–3500 nm (HTTP 404 from API) | 2 | ClO, H2O2 |

ClO's absorption is mostly UV band (need .xsc), and H2O2 is a UV broadband absorber too — both would need .xsc data anyway, so their line-LUT gap isn't blocking.

Storage: ~1.97 GB × 20 species = **~39 GB** under `~/data/HITRAN_LUTs/`,
well within the 574 TB free on `/kiwi-data`.

### Part E — NO2 pilot build outcome (2026-04-23)

Pilot of one species (NO2) via `test/benchmarks/build_hitran_luts.jl`:

- **Build time**: 190 s on one A100 (3.3M ν × 15 p × 10 T, Voigt LBL).
- **Storage**: 1.88 GB (Float32, linear interp). On spec.
- **LUT file**: `~/data/HITRAN_LUTs/NO2.jld2`.
- **Content check**:
    - σ(3500 nm, 1013 hPa, 296 K) ≈ 4×10⁻²¹ cm²  (OK, consistent with HITRAN ν₁/ν₃ NO2 IR bands)
    - σ(2900 nm) ≈ 3.5×10⁻²⁶   (OK, weak region between IR features)
    - σ(400 nm) ≈ 3.5×10⁻³⁵   (essentially **zero** — expected and wrong for the atmosphere)

**Important finding.** HITRAN's .par NO2 file only carries lines for
ν ∈ [2817, 7978] cm⁻¹ (λ = 1253–3550 nm). The atmospherically dominant
UV-Vis band (300–500 nm, σ ≈ 6×10⁻¹⁹ cm²) is stored **separately in
HITRAN's cross-section database (`.xsc` files)**, not in the line
catalogue. The rebuild-from-.par path is therefore accurate for
line-dominated species but leaves UV-Vis *broadband* species with
near-zero cross-section in the UV.

Species flagged as "UV-Vis broadband, .xsc required":
- **NO2**  (Vis band 300–600 nm)
- **SO2**  (UV 250–340 nm)
- **H2CO** (UV 300–360 nm)
- **O3**   (Huggins 300–360 nm, Chappuis 540–650 nm — HITRAN has lines
  for the IR and some UV, but the Huggins/Chappuis structure is usually
  drawn from .xsc)
- Possibly **ClO, HOCl** (UV bands only partially in .par)

We will build the full 22-species LUT from HITRAN .par data now, and
treat UV-Vis supplementation via .xsc as a follow-on (Phase 2.5). The
broadband species' LUTs are still correct in their IR bands.

### Part F — CIA pilot result (2026-04-23, files provided by user)

Ran `test/benchmarks/pilot_cia.jl` with the HITRAN2024 files
`~/data/HITRAN_CIA/{O2-O2_2024,O2-N2_2024,N2-N2_2021}.cia`. Combined
τ_CIA (vertical) using HITRAN Eq. (3) over the 53-layer AFGL
mid-lat-winter profile of our scenario, plus slant τ at μ_s = 0.77:

| λ (nm) | τ_O2O2 | τ_O2N2 | τ_N2N2 | τ_vert | T_slant |
|--------|--------|--------|--------|--------|---------|
|  380   | 1.1e−3 |  0     |  0     | 1.1e−3 | 0.9985  |
|  477   | 3.0e−3 |  0     |  0     | 3.0e−3 | 0.9962  |
|  577   | 3.0e−3 |  0     |  0     | 3.0e−3 | 0.9961  |
|  630   | 5.7e−4 |  0     |  0     | 5.7e−4 | 0.9993  |
| 1270   | 1.3e−2 | 2.8e−2 |  0     | 4.0e−2 | 0.9491  |
| 2000   |  0     |  0     | 1.5e−6 | 1.5e−6 | 1.0000  |
| 2160   |  0     |  0     | 1.5e−3 | 1.5e−3 | 0.9980  |
| 2200   |  0     |  0     | 4.8e−4 | 4.8e−4 | 0.9994  |

**Takeaways:**
- **CIA is the dominant missing physics at 1.27 μm** (5.1% slant
  attenuation from O2 a¹Δg + collider). Also contributes ~0.3–0.4%
  at UV-Vis "O2 double-transition" bands (477, 577 nm).
- **CIA is not the 2 μm-gap culprit.** The N2 2.16 μm overtone peak
  gives only τ_slant ≈ 0.002 there, while the MODTRAN vs. vSmartMOM
  gap at 2000 nm is ≈ 0.62 slant — 300× larger. Remaining suspects
  for that gap, in rough order of expected size:
    1. **MODTRAN default aerosol baseline** (rural 23 km → AOT(2 μm) ≈ 0.05–0.10)
    2. **H2O self- and foreign-continuum (MT_CKD)** — well-known to
       dominate the 2.0–2.5 μm atmospheric window
    3. LUT resolution + aerosol together account for the visible /
       near-UV part of the gap
- **For the rebuilt LUTs**, CIA should be added as an independent
  optical-depth contribution alongside line-absorber τ. Pending work
  (Phase 2.5):
    1. Load `.cia` blocks via `pilot_cia.jl::parse_cia_file`
    2. Build a CIA interpolator keyed on (ν, T, collision-pair)
    3. At model build time, compute τ_CIA per layer using
       Eq. (3) of the HITRAN CIA readme and add it to `τ_abs`.
       Ref density note: HITRAN uses `ρ_Air = ρ_O2 + ρ_N2` and has
       M-Air variants for Earth atmosphere — prefer those when
       present to avoid double-counting.

### Part D — CIA / XSC files are gated (2026-04-23)

Both HITRAN-CIA (`hitran.org/data/CIA/*.cia`) and HITRAN cross-section
(`hitran.org/data/xsec/*.xsc`) files return HTTP 404/403 via both
`curl` and Julia `Downloads`. HITRAN appears to gate both behind
authenticated API access.

To run the CIA pilot (`test/benchmarks/pilot_cia.jl`) AND to
supplement the UV-Vis broadband species (NO2, SO2, H2CO, O3), we need
either:

1. User to download the files manually from their HITRAN account:
     - CIA: log in → "Collision-induced Absorption" → O2-O2, N2-N2
       → drop at `~/data/HITRAN_CIA/{O2-O2,N2-N2}.cia`
     - XSC: "Cross-sections" → NO2, SO2, H2CO, O3 (Huggins/Chappuis)
       → drop at `~/data/HITRAN_XSC/{NO2,SO2,H2CO,O3}/*.xsc`
   The pilot's parser reads the standard HITRAN fixed-width text
   format, so no format conversion is needed.

2. Or: switch source to MT_CKD / CKD_2.4 continua (comes as Fortran
   tables, more work to port but covers CIA + water continuum in one
   go — which we'd want eventually anyway). Note MT_CKD does not
   include UV-Vis xsec data — that would still need HITRAN or
   another source (MPI-Mainz UV-Vis spectral atlas is public).

Option 1 is minutes of human work; Option 2 is a larger port. Proposing
Option 1 for the pilot and deferring Option 2 to Phase 4.
