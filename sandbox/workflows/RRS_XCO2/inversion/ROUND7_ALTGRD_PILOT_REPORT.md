# Round 7_altgrd pilot: fixed-altitude grid vs Round 7

Prepared 2026-10-18 for Suniti Sanghavi. This report describes and measures
the pilot; it does **not** decide whether a full re-run is needed.

## 1. Bottom line

* **No-aerosol scenes (80 solves): comparable.** XCO2 differs by mean
  +0.0007 ppm (SD 0.037, max |Δ| 0.12 ppm). That is ≤ 0.55 of the posterior
  XCO2 σ in every solve, and ≤ 0.61 of the Round 7 noise scatter. Δp_surf is
  ≤ 0.19 hPa. Convergence and fit-quality outcomes are identical. χ² is
  unchanged (|Δχ²_O2A| ≤ 0.005) except for the desert solves: there, 16 of 20
  took one more iteration (4 instead of 3), and their O₂A reduced χ² fell by
  up to 0.13 (e.g. 041/11: 0.335 → 0.209). Round 7 may have stopped slightly
  early there (it took one iteration fewer); the noisy Round 7 desert O₂A χ²
  values span 1.15–1.40.
* **Aerosol scenes (80 solves): not comparable by the criterion in §6.**
  XCO2 is systematically lower under the fixed-altitude grid: mean −0.49 ppm,
  SD 0.42 ppm. Only 39 % of these solves agree within 1σ.
  * Largest ΔXCO2: −1.52 ppm (078/01), i.e. 4.7 posterior σ and 3.5 × that
    scene's Round 7 noise scatter. The largest σ-normalized difference is
    6.0σ (059/01, −1.01 ppm). The largest relative to noise scatter is
    8.6 × (052/11, −0.89 ppm).
  * The retrieved aerosol partition also moves (mean Δ ln AOD_organic
    +0.37), and Δp_surf reaches 3.2 hPa.
  * By surface the mean difference is −0.81 ppm (desert), −0.71 (forest),
    −0.24 (rural) and −0.21 (urban).
  * Five of the 20 desert aerosol solves did not converge within the
    7-iteration limit (still fit_ok = 1). Two divergent trials raised the
    damping γ to 10⁴, and XCO2 was still drifting when they stopped. Their
    ΔXCO2 (−0.51 to −0.60 ppm) is therefore a snapshot, not a converged
    difference. Excluding them leaves the aerosol mean at −0.49 ppm.
* **Against truth** (noiseless member 11) the no-aerosol bias is unchanged
  (+0.225 ppm in both rounds). The aerosol bias changes sign and grows:
  Round 7 +0.09 ± 0.18 ppm, fixed-altitude grid −0.27 ± 0.29 ppm. Pooled
  over classes these partly cancel; the pooled numbers in §4 should not be
  read as an improvement.

## 2. What was run

| Item | Round 7 (reference) | Round 7_altgrd (pilot) |
|---|---|---|
| Grid convention for p_surf | `FixedPressureGrid` (interfaces above the surface fixed in pressure) | `FixedAltitudeGrid` (all interfaces fixed in altitude; p_half ∝ p_surf) |
| Code | `7acab57000dae259207a6760faae156cfa1734f6` | `4075fa16a0adfa5bf524b32d44b075854ab4e2d5` (branch `experiment/round7-altgrd` = 7acab570 + this change; local only) |
| State model | `round7_fixed_sif_imperfect_correction` | `round7_altgrd_fixed_sif_imperfect_correction` |
| Observations, noise, Se, seeds | sealed Round 7 release (manifest `a281bedd…`) | the same files, read-only |
| Priors, truth table, analyzer, components, SIF template, solar, RT config, Julia Manifest | Round 7 release `round7_fixed_sif_imperfect_correction_v1` (SHA256SUMS `7014e804…`) | the same files, verified against that manifest |
| Spectroscopy | ABSCO v5.2 + rebuilt HITRAN A-band H₂O LUT | same paths; contents hashed into the identity and sealed in `deployment/round7_altgrd_release.json` |
| Settings | Julia 1.12.5, Float32, GPU, 9 streams, 16 layers, 28 active coordinates, fixed SIF, tight UTLS | same |
| Members | 01–11 (880 solves) | 01 (noisy) and 11 (noiseless) for all 80 scenes: 160 solves |
| Hardware | curry GPU1 (all 440 no-aerosol solves and 165 SIF-on aerosol solves) and wurst GPU1 (the remaining 275 aerosol solves) | curry physical GPU1 only (UUID `GPU-c22ccdf7…`) |
| Campaign identity | off `98b4378d…`, on `4df29f6c…` | off `1edb9348…`, on `4bdecee2…` |

The only intended difference is the grid convention. Both conventions give a
bit-identical forward model at the truth surface pressure of 1000 hPa (see
the validation note `ROUND7_ALTGRD.md`), so the Round 7 observations are
valid for both.

Two genuine but secondary differences come with the convention. They are
part of the change, not separate choices:

* under `FixedAltitudeGrid`, the XCO2 dry-air weights do not depend on the
  retrieved p_surf (all layer columns scale together);
* the altitude-form aerosol layers stay fixed in altitude, so their
  p_surf tangent is zero.

One operational difference is not physics. Compiled Julia caches for this
branch were written to a private first depot
(`.julia_depot_altgrd_curry` in this output tree). The Round 7 curry depot
and `~/.julia` were read, never written. Package versions are fixed by the
identical Manifest.

### Execution log

* Smoke test (states 001 and 036, member 11): both converged and passed the
  fit gate. These two outputs are the official pilot outputs for 001/11 and
  036/11; the full pilot validated and skipped them.
* First launch: refused, because another session's process held GPU1
  (correct behaviour, no work done). Relaunched once GPU1 had been idle for
  60 s.
* The SIF-off block completed all 79 remaining solves (02:24–03:03 UTC).
* The SIF-on block failed twice at its first GPU kernel compile with a
  `GPUCompiler InvalidIRError`. This happened in a trivial `Float32`
  broadcast inside the aerosol optical-property cache, before any result was
  written. Both times, "module … missing from the cache" warnings appeared
  and packages were recompiled. At the same time another session was
  rewriting the shared `~/.julia/compiled` caches (writes from about
  02:58 UTC, a mass rewrite of about 600 files from 04:45 UTC).
  * The evidence is circumstantial. The SIF-off block finished cleanly while
    shared-cache writes were already happening. The second failure (about
    04:34 UTC) came before the mass rewrite. The successful diagnostic solve
    ran during it.
  * The same code then worked in three separate processes: a GPU probe, the
    Jacobian driver on scene 036, and a single-state runner solve of scene
    006. That diagnostic solve produced the official 006/11 output, which the
    final launch validated and skipped.
  * A third launch completed every remaining solve without failure.
  * Suspected cause: the shared-depot cache churn. The root cause was not
    isolated. A code/compiler interaction is not excluded, and the failure
    was not reproduced on demand.
  * No output was overwritten. The runner validated every existing file
    before skipping it.
  * The private first depot did not isolate the runs from `~/.julia/compiled`,
    because Julia still reads that depot. A full re-run should use a depot
    whose compiled caches nobody else rewrites.
* Final: 160 / 160 outputs written, 0 failed solves. Logs are in `logs/` of
  this tree.

Timing per forward+Jacobian evaluation (non-first trials, mean): Round 7
4.24 s, altgrd 4.86 s. Model build: 0.36 s vs 0.47 s. This is not a
controlled comparison: most Round 7 aerosol solves ran on wurst, and the
solve paths differ. Of the 0.62 s difference, 0.11 s is in model build
(where the extra per-layer pressure derivatives are computed); nothing was
profiled, so the split is not attributed further.

## 3. p_surf Jacobian columns at the truth state

These are instrument-space columns per hPa, from the production evaluator
(GPU, Float32) at the truth states of scene 011 (urban, aerosol, SIF off) and
scene 036 (rural, aerosol, SIF on), with p_surf = 1000 hPa. The forward
spectra are bit-identical between conventions; all other 29 columns are also
bit-identical.


| state | band | ‖K₁‖ pressure | ‖K₁‖ altitude | ratio | cos angle | ‖ΔK₁‖/‖K₁ pressure‖ |
|---|---|---:|---:|---:|---:|---:|
| s011_truth | O2A | 3.8086e-01 | 3.9242e-01 | 1.0304 | 0.991521 | 0.1356 |
| s011_truth | weak CO2 | 1.9504e-02 | 1.9741e-02 | 1.0122 | 0.966145 | 0.2621 |
| s011_truth | strong CO2 | 3.5842e-02 | 3.2146e-02 | 0.8969 | 0.967834 | 0.2614 |
| s011_truth | all | 3.8304e-01 | 3.9423e-01 | 1.0292 | 0.991203 | 0.1377 |
| s036_truth | O2A | 6.0850e-01 | 6.2664e-01 | 1.0298 | 0.991972 | 0.1320 |
| s036_truth | weak CO2 | 2.0495e-02 | 2.0753e-02 | 1.0126 | 0.966295 | 0.2616 |
| s036_truth | strong CO2 | 2.2104e-02 | 1.9786e-02 | 0.8951 | 0.968041 | 0.2612 |
| s036_truth | all | 6.0925e-01 | 6.2729e-01 | 1.0296 | 0.991905 | 0.1325 |

The altitude-grid column is about 3 % larger in the O₂ A band. It is about
10 % smaller in the strong CO₂ band. Its spectral shape differs by
13–26 % (cos angle 0.97 in the CO₂ bands). The reason is physical: every
layer's cross sections now respond to p_surf, instead of only the bottom
layer's. Both columns are confirmed by Float64 finite differences of the full
instrument-space forward model: at the best step, to 6e-9 (altitude) and
2e-7 (pressure). In Float32, which is what production uses, FD agreement
at a 2 hPa step is 1.2 % (altitude) and 2.3 % (pressure) raw, and 0.5–1.0 %
noise-weighted. At 0.5 hPa it is worse (up to 6.9 % raw), consistent with
Float32 cancellation in the FD. Details are in `ROUND7_ALTGRD.md` §5.

## 4. Summary statistics (altgrd − Round 7)

The tables below are written by
`retrieval_setup/compare_round7_altgrd_pilot.jl`. A few numbers quoted in §1
and §5 were computed by hand from the CSV with short Python checks: the
26/80 χ² count, the seven-largest list, the −0.490 ppm mean without the
non-converged solves, and the desert iteration/χ² counts.

σ_XCO2 is the Round 7 posterior XCO2 uncertainty of the same solve:
√(hᵀŜh), with h the gradient of each round's own XCO2 diagnostic. SD₁₋₁₀ is
the standard deviation of Round 7 XCO2 over noisy members 01–10 of the same
scene. XCO2 recomputed from the stored final states matches the stored XCO2
exactly in both rounds.

Coverage: pairs compared=160/160; solves whose Round 7 scene has all ten noisy members=160/160

| subset | n | ΔXCO2 mean | SD | max\|Δ\| | Δpsurf mean | SD | max\|Δ\| | max \|ΔXCO2\|/σ_XCO2 | mean \|ΔXCO2\|/σ_XCO2 | max \|ΔXCO2\|/SD₁₋₁₀ | mean \|ΔXCO2\|/SD₁₋₁₀ |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| member 11 (noiseless) | 80 | -0.1828 | 0.3135 | 0.9042 | -0.195 | 0.740 | 1.349 | 5.321 | 0.950 | 8.608 | 1.154 |
| member 01 (noisy) | 80 | -0.3104 | 0.4374 | 1.5235 | +0.112 | 0.924 | 3.157 | 5.978 | 1.387 | 8.571 | 1.652 |
| all | 160 | -0.2466 | 0.3847 | 1.5235 | -0.041 | 0.849 | 3.157 | 5.978 | 1.168 | 8.608 | 1.403 |
| no aerosol | 80 | +0.0007 | 0.0368 | 0.1185 | +0.009 | 0.043 | 0.187 | 0.553 | 0.147 | 0.608 | 0.177 |
| aerosol | 80 | -0.4939 | 0.4155 | 1.5235 | -0.091 | 1.201 | 3.157 | 5.978 | 2.190 | 8.608 | 2.629 |
| SIF off | 80 | -0.2487 | 0.3899 | 1.4815 | -0.041 | 0.862 | 2.529 | 5.975 | 1.176 | 8.608 | 1.497 |
| SIF on | 80 | -0.2445 | 0.3819 | 1.5235 | -0.042 | 0.841 | 3.157 | 5.978 | 1.161 | 6.711 | 1.309 |

Agreement by class (fractions of solves with |ΔXCO2| ≤ 0.5σ and ≤ 1σ, σ = Round 7 posterior):

| class | n | mean ΔXCO2 | max \|ΔXCO2\| | ≤0.5σ | ≤1σ | max \|Δ\|/σ | outcomes agree |
|---|---:|---:|---:|---:|---:|---:|---:|
| no aerosol all | 80 | +0.001 | 0.118 | 98.8 % | 100.0 % | 0.55 | 80/80 |
| no aerosol off | 40 | -0.002 | 0.118 | 97.5 % | 100.0 % | 0.55 | 40/40 |
| no aerosol on | 40 | +0.004 | 0.082 | 100.0 % | 100.0 % | 0.38 | 40/40 |
| no aerosol urban | 20 | -0.059 | 0.118 | 95.0 % | 100.0 % | 0.55 | 20/20 |
| no aerosol rural | 20 | +0.020 | 0.040 | 100.0 % | 100.0 % | 0.17 | 20/20 |
| no aerosol desert | 20 | +0.020 | 0.033 | 100.0 % | 100.0 % | 0.22 | 20/20 |
| no aerosol forest | 20 | +0.022 | 0.037 | 100.0 % | 100.0 % | 0.15 | 20/20 |
| aerosol all | 80 | -0.494 | 1.523 | 27.5 % | 38.8 % | 5.98 | 75/80 |
| aerosol off | 40 | -0.495 | 1.481 | 27.5 % | 40.0 % | 5.97 | 37/40 |
| aerosol on | 40 | -0.493 | 1.523 | 27.5 % | 37.5 % | 5.98 | 38/40 |
| aerosol urban | 20 | -0.215 | 1.481 | 45.0 % | 75.0 % | 5.97 | 20/20 |
| aerosol rural | 20 | -0.237 | 0.905 | 60.0 % | 75.0 % | 2.96 | 20/20 |
| aerosol desert | 20 | -0.810 | 1.014 | 0.0 % | 0.0 % | 5.98 | 15/20 |
| aerosol forest | 20 | -0.714 | 1.523 | 5.0 % | 5.0 % | 4.68 | 20/20 |

Truth XCO2 bias by class (retrieved − truth, ppm, mean ± SD):

| class | member | Round 7 | altgrd |
|---|---|---:|---:|
| no aerosol | 11 | +0.225 ± 0.190 | +0.225 ± 0.215 |
| no aerosol | 01 | +0.157 ± 0.258 | +0.158 ± 0.277 |
| aerosol | 11 | +0.092 ± 0.175 | -0.274 ± 0.292 |
| aerosol | 01 | +0.212 ± 0.313 | -0.410 ± 0.351 |

Solves that did not converge under altgrd:

- scene 054 member 01 (desert, aod760_0p28, SIF off): outcome=3, divergences=2, iterations 7 (R7 4), ΔXCO2 -0.562 ppm at stop, fit_ok=1
- scene 055 member 11 (desert, aod760_0p28, SIF off): outcome=3, divergences=2, iterations 7 (R7 4), ΔXCO2 -0.541 ppm at stop, fit_ok=1
- scene 055 member 01 (desert, aod760_0p28, SIF off): outcome=3, divergences=2, iterations 7 (R7 4), ΔXCO2 -0.510 ppm at stop, fit_ok=1
- scene 058 member 01 (desert, aod760_0p28, SIF angular_integral760_0p5): outcome=3, divergences=2, iterations 7 (R7 4), ΔXCO2 -0.597 ppm at stop, fit_ok=1
- scene 060 member 11 (desert, aod760_0p28, SIF angular_integral760_0p5): outcome=3, divergences=2, iterations 7 (R7 4), ΔXCO2 -0.537 ppm at stop, fit_ok=1

Round 7 noise-driven scatter (SD over members 01–10, per scene): XCO2 median 0.1874 ppm (range 0.0715–0.5059); psurf median 0.134 hPa

Truth errors (retrieved − truth XCO2, ppm), pooled over aerosol and no-aerosol scenes — see the by-class table above:

| subset | R7 mean | R7 SD | altgrd mean | altgrd SD |
|---|---:|---:|---:|---:|
| member 11 (noiseless) | +0.1585 | 0.1934 | -0.0243 | 0.3573 |
| member 01 (noisy) | +0.1846 | 0.2863 | -0.1258 | 0.4250 |

Truth errors in psurf (retrieved − 1000 hPa):

| subset | R7 mean | R7 SD | altgrd mean | altgrd SD |
|---|---:|---:|---:|---:|
| member 11 (noiseless) | +1.168 | 1.355 | +0.973 | 1.295 |
| member 01 (noisy) | +1.087 | 1.307 | +1.200 | 1.526 |
| all | +1.128 | 1.328 | +1.087 | 1.415 |
| no aerosol | -0.123 | 0.353 | -0.114 | 0.384 |
| aerosol | +2.378 | 0.506 | +2.287 | 0.981 |
| SIF off | +1.150 | 1.344 | +1.109 | 1.467 |
| SIF on | +1.106 | 1.320 | +1.064 | 1.370 |

Other state differences (altgrd − R7), aerosol scenes, all members:

| parameter | mean Δ | SD Δ | max \|Δ\| |
|---|---:|---:|---:|
| co2_vmr_layer16 | -5.885e-06 | 3.582e-06 | 1.446e-05 |
| ln_sulfate_aod760 | -0.08009 | 0.1616 | 0.6527 |
| ln_organic_carbon_aod760 | +0.3689 | 0.1234 | 0.7185 |
| ln_utls_sulfate_aod760 | -0.09479 | 0.1281 | 0.255 |
| ln_sulfate_z0 | -0.01363 | 0.08477 | 0.1868 |
| ln_organic_carbon_z0 | +0.0009949 | 0.01153 | 0.03408 |
| ln_utls_sulfate_z0 | -0.003722 | 0.004286 | 0.009012 |

Convergence/fit/iterations:

- member 11 (noiseless): converged R7=80 altgrd=78; fit_ok R7=80 altgrd=80; both converged+fit_ok=78/80; iterations equal=40, altgrd more=30, fewer=10; mean iterations R7=4.12 altgrd=4.30
- member 01 (noisy): converged R7=80 altgrd=77; fit_ok R7=80 altgrd=80; both converged+fit_ok=77/80; iterations equal=49, altgrd more=17, fewer=14; mean iterations R7=4.30 altgrd=4.21

Chi² per band (altgrd − R7):

| subset | band | mean Δχ² | max \|Δχ²\| | mean χ² R7 | mean χ² altgrd |
|---|---|---:|---:|---:|---:|
| member 11 (noiseless) | O2A | -0.0405 | 0.2535 | 0.1944 | 0.1540 |
| member 11 (noiseless) | weak CO2 | -0.0009 | 0.0039 | 0.0058 | 0.0049 |
| member 11 (noiseless) | strong CO2 | -0.0125 | 0.0545 | 0.0289 | 0.0164 |
| member 01 (noisy) | O2A | -0.0191 | 0.2538 | 1.1824 | 1.1632 |
| member 01 (noisy) | weak CO2 | -0.0009 | 0.0211 | 0.9946 | 0.9937 |
| member 01 (noisy) | strong CO2 | -0.0052 | 0.0548 | 1.0136 | 1.0084 |

Posterior σ (altgrd/R7 ratio):

- member 11 (noiseless): σ_XCO2 R7 mean 0.2378 ppm, altgrd mean 0.2543 ppm, ratio 1.0046–1.1960; σ_psurf R7 mean 0.234 hPa, altgrd mean 0.284 hPa, ratio 0.9918–1.3641
- member 01 (noisy): σ_XCO2 R7 mean 0.2377 ppm, altgrd mean 0.2542 ppm, ratio 1.0041–1.2061; σ_psurf R7 mean 0.233 hPa, altgrd mean 0.285 hPa, ratio 0.9718–1.3673


## 5. Interpretation

* The two conventions share the forward model only at 1000 hPa. The first
  step differs only through the p_surf Jacobian column. Every later trial
  pressure gives different interface pressures, spectra and Jacobian columns,
  because the forward pressure response changes together with its
  (consistent) tangent. In scenes without aerosol this changes the retrieved
  p_surf and XCO2 by small fractions of their uncertainty.
* In the aerosol scenes, p_surf, the sulfate/organic AOD split, aerosol
  heights and bottom-layer CO₂ are strongly correlated. Round 7 already
  retrieves p_surf about +2.4 hPa high there, and the fixed-altitude grid
  does not remove that bias (+2.3 hPa).
  * The state differences point at that degenerate direction: ln
    AOD_organic +0.37 and layer-16 CO₂ −5.9 ppm on average. I did not do a
    correlation analysis, so this is a reading of the state differences, not
    a demonstrated mechanism.
  * The fits are of similar quality. Mean χ² is slightly lower under the
    fixed-altitude grid, but 26 of the 80 aerosol solves have a higher total
    χ². χ² from the two forward models is not strictly comparable away from
    1000 hPa.
  * The posterior σ_XCO2 of the aerosol solves grows by 9–21 % (no aerosol:
    0.4–3 %).
* **Convergence path also matters.** The seven largest |ΔXCO2| (078/01,
  012/01, 072/01, 077/01, 074/01, 059/01, 080/01) are all noisy member 01.
  Six have a higher total χ² under the fixed-altitude grid, and five stopped
  after fewer iterations than Round 7 (4 vs 5–6). Part of the largest
  differences may therefore reflect where each solve stopped, not only a
  different optimum.
* The five non-converged desert solves (§1) fit either reading. Their
  aerosol solutions may be sensitive to the p_surf response, or this may be
  a plain optimizer-damping limit: two divergences, γ up to 10⁴, then too
  few iterations to recover.

## 6. Are the two rounds comparable?

Proposed criterion (for Suniti to accept or change). The rounds are called
comparable if all of the following hold:

1. convergence and fit-quality outcomes agree in ≥ 95 % of solves, with no
   new fit failures;
2. |ΔXCO2| ≤ 0.5 σ_XCO2 in ≥ 95 % of solves and ≤ 1 σ_XCO2 in all solves;
3. the mean ΔXCO2 of every scene class (aerosol/no aerosol, SIF off/on,
   surface) is within ±0.1 ppm, about half the median Round 7 noise scatter
   (0.19 ppm).

The per-class numbers are in the "Agreement by class" table in §4.

| Class | (1) outcomes agree | (2) ≤0.5σ / ≤1σ | (3) class means of ΔXCO2 | Verdict under the criterion |
|---|---|---|---|---|
| No aerosol (80) | 80/80 | 98.8 % / 100 % (max 0.55σ) | −0.06 to +0.02 ppm (SIF off/on, all surfaces) | comparable |
| Aerosol (80) | 75/80 (5 non-converged, desert) | 27.5 % / 38.8 % (max 6.0σ) | −0.81 to −0.21 ppm | not comparable |

In plain language: for clear-sky scenes the grid convention does not matter
at the precision of this study. For aerosol scenes it changes the retrieved
XCO2 by about half a ppm on average, and by up to 1.5 ppm. That is larger
than both the retrieval's own uncertainty and the noise scatter between
Round 7 members, and it reverses the sign of the aerosol-scene XCO2 bias
against truth. Whether this calls for a full re-run of the aerosol half, or
of the whole campaign, is your decision. The pilot covers only members 01
and 11, so member-to-member noise statistics under the fixed-altitude grid
(the SD over members 01–10) were not measured.

## 7. Which code earlier rounds ran on (output provenance)

The round-to-directory mapping comes from `RRS_XCO2/visualization/README.md`
(retrieval-type table) and `CHEATSHEET.md`. "Round 1/2" as names for indices
1/2 are my reading; the files do not use those words. Sources: NetCDF global
attributes (every file), `logs/*provenance*`, `.control/` identity files, and
`git log` of `uni_vSmartMOM_round5_jacobian`.

| Round | Outputs | Code recorded in outputs | created_utc range |
|---|---|---|---|
| 1 (full column) | `RRS_XCO2/inversion/{corrected,uncorrected}` (704) | none | 2026-08-31 15:54 – 09-04 05:09 |
| 2 (bottom layer) | `bottom_layer_XCO2_retrievals/retrievals/` (1078) | none | 09-03 22:33 – 09-05 11:08 |
| 3 no-SIF | `retrievals_acos_mapped_tapered_vertical_correlation_nosif/` | none (checkout HEAD was cb90b83b throughout, by reflog — inferred) | 09-05 01:54 – 09-06 21:53 |
| 3 SIF-on | `~/RRS_XCO2_private/results/bottom_layer_sif_acos_mapped_tapered_vertical_correlation_v1/` | `source_sha=cb90b83b…` in production provenance logs | 09-05 01:39 – 09-07 08:56 |
| 4 no-SIF | `round4_known_sif759/retrievals_nosif/` | `round4_code_checkpoint_sha=cb90b83be9101209d2e23292d7085dc432e02d3e` (all 880) | 09-07 02:20 – 09-09 07:53 |
| 4 SIF-on | `~/RRS_XCO2_private/results/bottom_layer_round4_known_sif759_sif_on_…/retrievals` | `round4_code_checkpoint_sha=630c9190e82632660636d87ce99afd8935084bd3` (all 880) | 09-07 16:31 – 09-09 21:47 |
| 5, 6, 7 | `round5_fixed_sif`, `round6_fixed_sif`, `round7_fixed_sif_imperfect_correction` | `7acab57000dae259207a6760faae156cfa1734f6` | R5 09-09/10; R6 09-10/11; R7 09-26 – 09-28 |

**The bottom-layer dσ/dp term.**

* `36e7d1e0` (cfranken, 2026-09-08 14:29 −0700) adds
  `src/CoreRT/tools/absorption_pressure.jl`, with
  `_absorption_pressure_derivative` and `_accumulate_absorption_pressure!`.
  `git log --all -S` finds no earlier occurrence. It is an ancestor of
  7acab570, so **Rounds 5–7 include the term**.
* **Round 4 did not include it**, confirmed by the outputs:
  * Both recorded Round 4 commits (cb90b83b of 09-04, and 630c9190 of 09-07,
    whose `src/` is identical) lie on `suniti_multi_sensor`. They are neither
    ancestors nor descendants of 36e7d1e0; the merge base is 55a0a73f
    (2026-08-27).
  * Neither contains `absorption_pressure.jl`. The comment in cb90b83b's
    `lin_model_from_parameters.jl` states that the bottom-layer pressure
    dependence is "through molecular column only".
* The statement "Round 4 outputs are dated 2026-09-07" is only partly
  accurate. Round 4 outputs span 2026-09-07 02:20 to 09-09 21:47 UTC, so
  about half were written after 36e7d1e0 was committed. The dates do not
  matter, though: Round 4 ran from a different branch that never contained
  the term.
* Rounds 1–2 and Round 3 no-SIF record no code SHA. Rounds 1–2 ran from the
  moving working checkout `uni_vSmartMOM` and are unpinned. They finished by
  2026-09-05, before 36e7d1e0 existed (09-08), and the term first appears in
  git history in that commit. That makes the term's presence very unlikely,
  but it cannot exclude uncommitted local code. Round 3 SIF-on records
  cb90b83b, which does not contain the term.

## 8. Files

* This report, a per-solve CSV (`round7_altgrd_vs_round7_pilot.csv`), and the
  outputs `retrievals_{nosif,sif}/corrected/retrieval_stateNNN_perturbation{01,11}.nc`
  are in this tree. Logs are in `logs/`, and the sealed pin is in
  `deployment/round7_altgrd_release.json`.
* The comparison script is
  `sandbox/workflows/RRS_XCO2/inversion/retrieval_setup/compare_round7_altgrd_pilot.jl`
  on the branch. The Jacobian drivers are `vertical_grid_jacobian_states.jl`
  and `vertical_grid_psurf_fd.jl` in the same directory.
* Round 7 was read only. Nothing in its tree was modified.

## Appendix: every scene and member

χ² is the final reduced χ² per band. σXCO2 is each round's posterior XCO2 σ.
The Round 7 noise scatter and posterior p_surf σ per solve are in the CSV
(`round7_altgrd_vs_round7_pilot.csv`), which has the same rows.

| scene | mem | surface | aerosol | SIF | XCO2 R7 | XCO2 altgrd | ΔXCO2 | psurf R7 | psurf altgrd | Δpsurf | σXCO2 R7 / altgrd | χ² R7 (O2A, wCO2, sCO2) | χ² altgrd (O2A, wCO2, sCO2) | conv R7/alt | fit R7/alt | iter R7/alt |
|---|---|---|---|---|---:|---:|---:|---:|---:|---:|---|---|---|---|---|---|
| 001 | 11 | urban | none | off | 397.5221 | 397.4634 | -0.0587 | 1000.439 | 1000.505 | +0.067 | 0.212 / 0.216 | 0.212, 0.003, 0.009 | 0.212, 0.002, 0.009 | 1/1 | 1/1 | 4/4 |
| 001 | 01 | urban | none | off | 397.3772 | 397.3247 | -0.0525 | 1000.630 | 1000.652 | +0.023 | 0.212 / 0.216 | 1.196, 0.959, 0.969 | 1.196, 0.957, 0.972 | 1/1 | 1/1 | 4/4 |
| 002 | 11 | urban | none | off | 398.7701 | 398.7104 | -0.0597 | 1000.440 | 1000.508 | +0.069 | 0.213 / 0.216 | 0.212, 0.003, 0.009 | 0.212, 0.003, 0.009 | 1/1 | 1/1 | 4/4 |
| 002 | 01 | urban | none | off | 398.5560 | 398.5000 | -0.0560 | 1000.397 | 1000.463 | +0.066 | 0.213 / 0.216 | 1.244, 0.975, 0.977 | 1.244, 0.975, 0.977 | 1/1 | 1/1 | 4/4 |
| 003 | 11 | urban | none | off | 400.0183 | 399.9575 | -0.0607 | 1000.441 | 1000.511 | +0.070 | 0.213 / 0.217 | 0.212, 0.003, 0.009 | 0.212, 0.003, 0.009 | 1/1 | 1/1 | 4/4 |
| 003 | 01 | urban | none | off | 400.3210 | 400.2418 | -0.0792 | 1000.471 | 1000.578 | +0.107 | 0.214 / 0.217 | 1.202, 0.956, 1.012 | 1.201, 0.955, 1.012 | 1/1 | 1/1 | 4/4 |
| 004 | 11 | urban | none | off | 401.2663 | 401.2046 | -0.0617 | 1000.442 | 1000.515 | +0.072 | 0.214 / 0.217 | 0.212, 0.003, 0.009 | 0.212, 0.003, 0.009 | 1/1 | 1/1 | 4/4 |
| 004 | 01 | urban | none | off | 401.3037 | 401.1853 | -0.1185 | 1000.548 | 1000.735 | +0.187 | 0.214 / 0.218 | 1.185, 1.047, 0.985 | 1.184, 1.046, 0.982 | 1/1 | 1/1 | 4/4 |
| 005 | 11 | urban | none | off | 402.5144 | 402.4517 | -0.0627 | 1000.444 | 1000.518 | +0.074 | 0.214 / 0.218 | 0.212, 0.003, 0.009 | 0.212, 0.003, 0.009 | 1/1 | 1/1 | 4/4 |
| 005 | 01 | urban | none | off | 402.3940 | 402.3627 | -0.0313 | 1000.508 | 1000.501 | -0.007 | 0.214 / 0.218 | 1.253, 0.988, 0.991 | 1.251, 0.987, 0.996 | 1/1 | 1/1 | 4/4 |
| 006 | 11 | urban | none | on | 397.5244 | 397.4659 | -0.0585 | 1000.434 | 1000.500 | +0.067 | 0.212 / 0.216 | 0.212, 0.003, 0.009 | 0.212, 0.002, 0.009 | 1/1 | 1/1 | 4/4 |
| 006 | 01 | urban | none | on | 397.3939 | 397.3568 | -0.0372 | 1000.396 | 1000.418 | +0.022 | 0.212 / 0.216 | 1.228, 0.916, 0.945 | 1.227, 0.916, 0.947 | 1/1 | 1/1 | 4/4 |
| 007 | 11 | urban | none | on | 398.7723 | 398.7130 | -0.0593 | 1000.435 | 1000.504 | +0.068 | 0.213 / 0.216 | 0.212, 0.003, 0.009 | 0.212, 0.003, 0.009 | 1/1 | 1/1 | 4/4 |
| 007 | 01 | urban | none | on | 398.8047 | 398.7613 | -0.0434 | 1000.540 | 1000.555 | +0.015 | 0.213 / 0.217 | 1.198, 0.997, 1.000 | 1.198, 0.997, 1.003 | 1/1 | 1/1 | 4/4 |
| 008 | 11 | urban | none | on | 400.0203 | 399.9600 | -0.0603 | 1000.437 | 1000.507 | +0.070 | 0.213 / 0.217 | 0.212, 0.003, 0.009 | 0.212, 0.003, 0.009 | 1/1 | 1/1 | 4/4 |
| 008 | 01 | urban | none | on | 399.5624 | 399.5103 | -0.0520 | 1000.424 | 1000.479 | +0.055 | 0.213 / 0.217 | 1.250, 0.988, 1.014 | 1.251, 0.989, 1.013 | 1/1 | 1/1 | 4/4 |
| 009 | 11 | urban | none | on | 401.2685 | 401.2072 | -0.0613 | 1000.438 | 1000.510 | +0.072 | 0.214 / 0.217 | 0.211, 0.003, 0.009 | 0.212, 0.003, 0.009 | 1/1 | 1/1 | 4/4 |
| 009 | 01 | urban | none | on | 401.1093 | 401.0834 | -0.0259 | 1000.536 | 1000.513 | -0.023 | 0.213 / 0.217 | 1.184, 1.005, 0.973 | 1.183, 1.004, 0.979 | 1/1 | 1/1 | 4/4 |
| 010 | 11 | urban | none | on | 402.5166 | 402.4543 | -0.0623 | 1000.439 | 1000.513 | +0.074 | 0.214 / 0.218 | 0.211, 0.003, 0.009 | 0.212, 0.003, 0.009 | 1/1 | 1/1 | 4/4 |
| 010 | 01 | urban | none | on | 402.3841 | 402.3021 | -0.0821 | 1000.383 | 1000.510 | +0.127 | 0.215 / 0.218 | 1.178, 0.988, 0.983 | 1.177, 0.987, 0.983 | 1/1 | 1/1 | 4/4 |
| 011 | 11 | urban | aer | off | 397.3805 | 397.4327 | +0.0522 | 1002.913 | 1001.565 | -1.349 | 0.248 / 0.276 | 0.200, 0.008, 0.059 | 0.051, 0.005, 0.008 | 1/1 | 1/1 | 6/4 |
| 011 | 01 | urban | aer | off | 397.4486 | 397.5692 | +0.1206 | 1002.862 | 1001.441 | -1.421 | 0.248 / 0.276 | 1.243, 1.005, 1.100 | 1.099, 1.003, 1.046 | 1/1 | 1/1 | 6/4 |
| 012 | 11 | urban | aer | off | 398.6796 | 398.6757 | -0.0039 | 1002.918 | 1001.580 | -1.338 | 0.249 / 0.277 | 0.200, 0.007, 0.060 | 0.051, 0.005, 0.008 | 1/1 | 1/1 | 6/4 |
| 012 | 01 | urban | aer | off | 398.8198 | 397.3383 | -1.4815 | 1002.735 | 1005.264 | +2.529 | 0.248 / 0.270 | 1.226, 0.905, 1.067 | 1.183, 0.926, 1.110 | 1/1 | 1/1 | 6/4 |
| 013 | 11 | urban | aer | off | 399.9789 | 399.9194 | -0.0596 | 1002.922 | 1001.595 | -1.328 | 0.250 / 0.278 | 0.200, 0.007, 0.060 | 0.051, 0.005, 0.007 | 1/1 | 1/1 | 6/4 |
| 013 | 01 | urban | aer | off | 399.8727 | 399.6510 | -0.2217 | 1002.948 | 1001.655 | -1.293 | 0.250 / 0.279 | 1.176, 0.977, 1.014 | 1.035, 0.971, 0.966 | 1/1 | 1/1 | 6/4 |
| 014 | 11 | urban | aer | off | 401.2778 | 401.1622 | -0.1156 | 1002.927 | 1001.608 | -1.318 | 0.251 / 0.279 | 0.200, 0.006, 0.061 | 0.052, 0.005, 0.007 | 1/1 | 1/1 | 6/4 |
| 014 | 01 | urban | aer | off | 401.4649 | 401.1841 | -0.2809 | 1002.996 | 1002.398 | -0.598 | 0.251 / 0.278 | 1.188, 1.009, 1.042 | 1.049, 1.008, 0.996 | 1/1 | 1/1 | 6/3 |
| 015 | 11 | urban | aer | off | 402.5767 | 402.4054 | -0.1713 | 1002.932 | 1001.620 | -1.312 | 0.252 / 0.280 | 0.200, 0.005, 0.062 | 0.051, 0.005, 0.007 | 1/1 | 1/1 | 6/4 |
| 015 | 01 | urban | aer | off | 402.8098 | 402.6577 | -0.1520 | 1002.831 | 1001.572 | -1.259 | 0.252 / 0.280 | 1.092, 1.001, 1.078 | 0.958, 0.998, 1.025 | 1/1 | 1/1 | 6/4 |
| 016 | 11 | urban | aer | on | 397.3811 | 397.4321 | +0.0509 | 1002.912 | 1001.570 | -1.342 | 0.248 / 0.276 | 0.200, 0.008, 0.059 | 0.051, 0.005, 0.008 | 1/1 | 1/1 | 6/4 |
| 016 | 01 | urban | aer | on | 397.5436 | 397.2813 | -0.2624 | 1002.849 | 1002.649 | -0.200 | 0.248 / 0.274 | 1.190, 0.967, 1.033 | 1.097, 0.966, 1.006 | 1/1 | 1/1 | 6/4 |
| 017 | 11 | urban | aer | on | 398.6798 | 398.6753 | -0.0045 | 1002.916 | 1001.583 | -1.333 | 0.249 / 0.277 | 0.200, 0.007, 0.060 | 0.051, 0.005, 0.008 | 1/1 | 1/1 | 6/4 |
| 017 | 01 | urban | aer | on | 398.6187 | 398.0889 | -0.5298 | 1002.862 | 1002.120 | -0.741 | 0.248 / 0.277 | 1.165, 0.996, 1.036 | 1.072, 0.990, 1.012 | 1/1 | 1/1 | 7/3 |
| 018 | 11 | urban | aer | on | 399.9788 | 399.9189 | -0.0600 | 1002.920 | 1001.597 | -1.323 | 0.250 / 0.278 | 0.200, 0.007, 0.060 | 0.051, 0.005, 0.007 | 1/1 | 1/1 | 6/4 |
| 018 | 01 | urban | aer | on | 399.8085 | 399.6675 | -0.1411 | 1003.063 | 1001.887 | -1.176 | 0.250 / 0.278 | 1.210, 1.013, 1.041 | 1.058, 1.012, 0.990 | 1/1 | 1/1 | 6/4 |
| 019 | 11 | urban | aer | on | 401.2776 | 401.1617 | -0.1159 | 1002.925 | 1001.612 | -1.313 | 0.251 / 0.279 | 0.200, 0.006, 0.061 | 0.051, 0.005, 0.007 | 1/1 | 1/1 | 6/4 |
| 019 | 01 | urban | aer | on | 401.3772 | 400.8357 | -0.5415 | 1002.664 | 1002.570 | -0.094 | 0.250 / 0.277 | 1.193, 1.004, 1.074 | 1.098, 1.003, 1.048 | 1/1 | 1/1 | 7/3 |
| 020 | 11 | urban | aer | on | 402.5769 | 402.4048 | -0.1721 | 1002.930 | 1001.624 | -1.306 | 0.252 / 0.280 | 0.200, 0.005, 0.062 | 0.051, 0.005, 0.007 | 1/1 | 1/1 | 6/4 |
| 020 | 01 | urban | aer | on | 402.7175 | 402.5172 | -0.2003 | 1002.789 | 1001.675 | -1.114 | 0.252 / 0.280 | 1.152, 0.966, 1.034 | 1.020, 0.962, 0.985 | 1/1 | 1/1 | 6/4 |
| 021 | 11 | rural | none | off | 397.8750 | 397.8957 | +0.0207 | 999.695 | 999.683 | -0.013 | 0.238 / 0.239 | 0.203, 0.003, 0.021 | 0.205, 0.003, 0.019 | 1/1 | 1/1 | 4/4 |
| 021 | 01 | rural | none | off | 397.9242 | 397.9435 | +0.0193 | 999.713 | 999.701 | -0.011 | 0.238 / 0.240 | 1.240, 0.992, 1.006 | 1.242, 0.992, 1.004 | 1/1 | 1/1 | 4/4 |
| 022 | 11 | rural | none | off | 399.1227 | 399.1430 | +0.0204 | 999.696 | 999.684 | -0.012 | 0.238 / 0.240 | 0.202, 0.003, 0.021 | 0.205, 0.003, 0.019 | 1/1 | 1/1 | 4/4 |
| 022 | 01 | rural | none | off | 399.1229 | 399.1390 | +0.0161 | 999.493 | 999.516 | +0.023 | 0.238 / 0.240 | 1.203, 1.012, 0.982 | 1.207, 1.012, 0.980 | 1/1 | 1/1 | 4/4 |
| 023 | 11 | rural | none | off | 400.3703 | 400.3903 | +0.0200 | 999.697 | 999.686 | -0.011 | 0.239 / 0.240 | 0.202, 0.003, 0.021 | 0.205, 0.003, 0.019 | 1/1 | 1/1 | 4/4 |
| 023 | 01 | rural | none | off | 400.4352 | 400.4482 | +0.0130 | 999.586 | 999.604 | +0.018 | 0.239 / 0.240 | 1.285, 1.048, 1.026 | 1.289, 1.048, 1.023 | 1/1 | 1/1 | 4/4 |
| 024 | 11 | rural | none | off | 401.6181 | 401.6377 | +0.0196 | 999.698 | 999.687 | -0.011 | 0.239 / 0.241 | 0.202, 0.003, 0.021 | 0.204, 0.004, 0.019 | 1/1 | 1/1 | 4/4 |
| 024 | 01 | rural | none | off | 400.9658 | 400.9780 | +0.0122 | 999.868 | 999.854 | -0.015 | 0.239 / 0.240 | 1.212, 0.979, 1.063 | 1.214, 0.979, 1.062 | 1/1 | 1/1 | 4/4 |
| 025 | 11 | rural | none | off | 402.8659 | 402.8850 | +0.0192 | 999.698 | 999.688 | -0.010 | 0.240 / 0.241 | 0.202, 0.004, 0.022 | 0.204, 0.004, 0.020 | 1/1 | 1/1 | 4/4 |
| 025 | 01 | rural | none | off | 402.7003 | 402.7017 | +0.0014 | 999.664 | 999.699 | +0.035 | 0.240 / 0.241 | 1.202, 0.954, 0.993 | 1.206, 0.954, 0.992 | 1/1 | 1/1 | 4/4 |
| 026 | 11 | rural | none | on | 397.8765 | 397.8974 | +0.0209 | 999.692 | 999.679 | -0.013 | 0.238 / 0.239 | 0.202, 0.003, 0.021 | 0.205, 0.003, 0.019 | 1/1 | 1/1 | 4/4 |
| 026 | 01 | rural | none | on | 397.7112 | 397.7369 | +0.0256 | 999.586 | 999.575 | -0.011 | 0.238 / 0.239 | 1.159, 1.008, 1.039 | 1.164, 1.008, 1.034 | 1/1 | 1/1 | 4/4 |
| 027 | 11 | rural | none | on | 399.1242 | 399.1447 | +0.0205 | 999.693 | 999.681 | -0.012 | 0.238 / 0.240 | 0.202, 0.003, 0.021 | 0.205, 0.003, 0.019 | 1/1 | 1/1 | 4/4 |
| 027 | 01 | rural | none | on | 399.2255 | 399.2655 | +0.0400 | 999.716 | 999.656 | -0.060 | 0.238 / 0.240 | 1.146, 0.999, 1.089 | 1.147, 0.999, 1.087 | 1/1 | 1/1 | 4/4 |
| 028 | 11 | rural | none | on | 400.3720 | 400.3920 | +0.0200 | 999.693 | 999.682 | -0.011 | 0.239 / 0.240 | 0.202, 0.003, 0.021 | 0.205, 0.003, 0.019 | 1/1 | 1/1 | 4/4 |
| 028 | 01 | rural | none | on | 399.7143 | 399.7417 | +0.0274 | 999.754 | 999.718 | -0.036 | 0.239 / 0.240 | 1.140, 1.049, 1.015 | 1.142, 1.049, 1.013 | 1/1 | 1/1 | 4/4 |
| 029 | 11 | rural | none | on | 401.6197 | 401.6392 | +0.0196 | 999.694 | 999.684 | -0.010 | 0.239 / 0.241 | 0.202, 0.003, 0.021 | 0.204, 0.004, 0.019 | 1/1 | 1/1 | 4/4 |
| 029 | 01 | rural | none | on | 401.6544 | 401.6730 | +0.0186 | 999.752 | 999.737 | -0.015 | 0.239 / 0.241 | 1.192, 1.020, 1.006 | 1.194, 1.020, 1.005 | 1/1 | 1/1 | 4/4 |
| 030 | 11 | rural | none | on | 402.8675 | 402.8867 | +0.0192 | 999.695 | 999.686 | -0.009 | 0.240 / 0.241 | 0.202, 0.004, 0.022 | 0.204, 0.004, 0.020 | 1/1 | 1/1 | 4/4 |
| 030 | 01 | rural | none | on | 402.8173 | 402.8504 | +0.0330 | 999.771 | 999.721 | -0.050 | 0.240 / 0.241 | 1.221, 1.012, 1.013 | 1.222, 1.013, 1.010 | 1/1 | 1/1 | 4/4 |
| 031 | 11 | rural | aer | off | 397.5312 | 397.5691 | +0.0379 | 1002.403 | 1001.163 | -1.240 | 0.302 / 0.334 | 0.138, 0.003, 0.036 | 0.039, 0.002, 0.005 | 1/1 | 1/1 | 4/5 |
| 031 | 01 | rural | aer | off | 398.0042 | 397.7762 | -0.2280 | 1001.663 | 1001.147 | -0.515 | 0.300 / 0.333 | 1.093, 0.976, 0.971 | 1.048, 0.972, 0.963 | 1/1 | 1/1 | 5/5 |
| 032 | 11 | rural | aer | off | 398.8075 | 398.8101 | +0.0026 | 1002.410 | 1001.178 | -1.232 | 0.303 / 0.335 | 0.138, 0.003, 0.036 | 0.039, 0.002, 0.005 | 1/1 | 1/1 | 4/5 |
| 032 | 01 | rural | aer | off | 399.0813 | 398.2679 | -0.8133 | 1002.255 | 1003.318 | +1.062 | 0.304 / 0.334 | 1.102, 1.043, 1.002 | 1.077, 1.040, 0.993 | 1/1 | 1/1 | 4/4 |
| 033 | 11 | rural | aer | off | 400.0856 | 400.0540 | -0.0316 | 1002.416 | 1001.184 | -1.232 | 0.304 / 0.337 | 0.138, 0.003, 0.036 | 0.039, 0.002, 0.005 | 1/1 | 1/1 | 4/5 |
| 033 | 01 | rural | aer | off | 399.7942 | 398.9925 | -0.8016 | 1002.344 | 1003.207 | +0.863 | 0.303 / 0.335 | 1.089, 0.971, 1.006 | 1.058, 0.972, 1.003 | 1/1 | 1/1 | 4/4 |
| 034 | 11 | rural | aer | off | 401.3624 | 401.2965 | -0.0659 | 1002.423 | 1001.195 | -1.228 | 0.305 / 0.338 | 0.137, 0.003, 0.036 | 0.039, 0.002, 0.005 | 1/1 | 1/1 | 4/5 |
| 034 | 01 | rural | aer | off | 401.5584 | 400.6533 | -0.9051 | 1002.565 | 1003.434 | +0.869 | 0.306 / 0.337 | 1.229, 0.957, 1.017 | 1.194, 0.957, 1.011 | 1/1 | 1/1 | 4/4 |
| 035 | 11 | rural | aer | off | 402.6392 | 402.5365 | -0.1027 | 1002.430 | 1001.213 | -1.217 | 0.306 / 0.339 | 0.137, 0.003, 0.037 | 0.039, 0.002, 0.005 | 1/1 | 1/1 | 4/5 |
| 035 | 01 | rural | aer | off | 402.4383 | 402.2349 | -0.2034 | 1002.567 | 1001.318 | -1.250 | 0.307 / 0.340 | 1.159, 1.023, 1.011 | 1.054, 1.021, 0.977 | 1/1 | 1/1 | 4/5 |
| 036 | 11 | rural | aer | on | 397.5309 | 397.5689 | +0.0379 | 1002.402 | 1001.167 | -1.235 | 0.302 / 0.334 | 0.138, 0.003, 0.036 | 0.039, 0.002, 0.005 | 1/1 | 1/1 | 4/5 |
| 036 | 01 | rural | aer | on | 397.5339 | 397.4314 | -0.1026 | 1002.227 | 1001.484 | -0.743 | 0.302 / 0.334 | 1.149, 1.022, 1.031 | 1.044, 1.023, 1.002 | 1/1 | 1/1 | 4/5 |
| 037 | 11 | rural | aer | on | 398.8089 | 398.8099 | +0.0011 | 1002.408 | 1001.181 | -1.227 | 0.303 / 0.335 | 0.138, 0.003, 0.036 | 0.039, 0.002, 0.005 | 1/1 | 1/1 | 4/5 |
| 037 | 01 | rural | aer | on | 398.9257 | 398.8905 | -0.0352 | 1002.531 | 1001.078 | -1.453 | 0.304 / 0.337 | 1.130, 1.025, 1.010 | 1.039, 1.023, 0.976 | 1/1 | 1/1 | 4/5 |
| 038 | 11 | rural | aer | on | 400.0854 | 400.0554 | -0.0300 | 1002.414 | 1001.182 | -1.232 | 0.304 / 0.337 | 0.137, 0.003, 0.036 | 0.039, 0.002, 0.005 | 1/1 | 1/1 | 4/5 |
| 038 | 01 | rural | aer | on | 399.6835 | 399.5116 | -0.1718 | 1002.469 | 1001.498 | -0.972 | 0.304 / 0.336 | 1.155, 1.019, 0.985 | 1.057, 1.016, 0.960 | 1/1 | 1/1 | 4/5 |
| 039 | 11 | rural | aer | on | 401.3631 | 401.2972 | -0.0658 | 1002.422 | 1001.198 | -1.224 | 0.305 / 0.338 | 0.137, 0.003, 0.036 | 0.039, 0.002, 0.005 | 1/1 | 1/1 | 4/5 |
| 039 | 01 | rural | aer | on | 400.9701 | 400.6497 | -0.3204 | 1002.783 | 1001.857 | -0.926 | 0.305 / 0.339 | 1.073, 1.000, 1.036 | 0.974, 0.999, 1.001 | 1/1 | 1/1 | 4/5 |
| 040 | 11 | rural | aer | on | 402.6391 | 402.5354 | -0.1038 | 1002.430 | 1001.217 | -1.212 | 0.306 / 0.339 | 0.137, 0.003, 0.037 | 0.039, 0.002, 0.005 | 1/1 | 1/1 | 4/5 |
| 040 | 01 | rural | aer | on | 402.8780 | 402.0330 | -0.8450 | 1002.210 | 1003.214 | +1.004 | 0.306 / 0.338 | 1.121, 0.956, 1.046 | 1.097, 0.956, 1.042 | 1/1 | 1/1 | 4/4 |
| 041 | 11 | desert | none | off | 397.5645 | 397.5844 | +0.0199 | 999.771 | 999.746 | -0.025 | 0.146 / 0.151 | 0.335, 0.009, 0.018 | 0.209, 0.006, 0.010 | 1/1 | 1/1 | 3/4 |
| 041 | 01 | desert | none | off | 397.6617 | 397.6836 | +0.0220 | 999.824 | 999.800 | -0.024 | 0.146 / 0.150 | 1.340, 1.033, 1.037 | 1.214, 1.030, 1.031 | 1/1 | 1/1 | 3/4 |
| 042 | 11 | desert | none | off | 398.8112 | 398.8331 | +0.0219 | 999.773 | 999.753 | -0.019 | 0.147 / 0.151 | 0.335, 0.010, 0.018 | 0.250, 0.008, 0.012 | 1/1 | 1/1 | 3/4 |
| 042 | 01 | desert | none | off | 398.7922 | 398.7892 | -0.0030 | 999.934 | 999.950 | +0.016 | 0.147 / 0.151 | 1.232, 1.032, 0.954 | 1.234, 1.032, 0.953 | 1/1 | 1/1 | 4/4 |
| 043 | 11 | desert | none | off | 400.0580 | 400.0791 | +0.0211 | 999.774 | 999.758 | -0.017 | 0.147 / 0.151 | 0.334, 0.010, 0.018 | 0.250, 0.008, 0.013 | 1/1 | 1/1 | 3/4 |
| 043 | 01 | desert | none | off | 400.2918 | 400.3021 | +0.0103 | 999.825 | 999.832 | +0.008 | 0.147 / 0.151 | 1.381, 0.953, 0.953 | 1.298, 0.951, 0.950 | 1/1 | 1/1 | 3/4 |
| 044 | 11 | desert | none | off | 401.3048 | 401.3253 | +0.0204 | 999.776 | 999.762 | -0.015 | 0.147 / 0.151 | 0.334, 0.010, 0.018 | 0.250, 0.008, 0.013 | 1/1 | 1/1 | 3/4 |
| 044 | 01 | desert | none | off | 401.3144 | 401.3473 | +0.0329 | 999.693 | 999.671 | -0.022 | 0.147 / 0.151 | 1.346, 1.018, 0.953 | 1.262, 1.017, 0.947 | 1/1 | 1/1 | 3/4 |
| 045 | 11 | desert | none | off | 402.5517 | 402.5714 | +0.0196 | 999.778 | 999.766 | -0.012 | 0.148 / 0.152 | 0.334, 0.011, 0.019 | 0.250, 0.009, 0.013 | 1/1 | 1/1 | 3/4 |
| 045 | 01 | desert | none | off | 402.5335 | 402.5621 | +0.0287 | 999.787 | 999.758 | -0.029 | 0.147 / 0.152 | 1.252, 0.994, 1.062 | 1.254, 0.996, 1.060 | 1/1 | 1/1 | 4/4 |
| 046 | 11 | desert | none | on | 397.5656 | 397.5860 | +0.0204 | 999.769 | 999.743 | -0.026 | 0.146 / 0.151 | 0.334, 0.009, 0.018 | 0.208, 0.006, 0.011 | 1/1 | 1/1 | 3/4 |
| 046 | 01 | desert | none | on | 397.5429 | 397.5701 | +0.0272 | 999.744 | 999.727 | -0.017 | 0.146 / 0.151 | 1.148, 1.010, 1.019 | 1.151, 1.010, 1.017 | 1/1 | 1/1 | 4/4 |
| 047 | 11 | desert | none | on | 398.8124 | 398.8346 | +0.0222 | 999.770 | 999.750 | -0.020 | 0.147 / 0.151 | 0.334, 0.010, 0.018 | 0.250, 0.008, 0.012 | 1/1 | 1/1 | 3/4 |
| 047 | 01 | desert | none | on | 398.9264 | 398.9484 | +0.0220 | 999.762 | 999.752 | -0.010 | 0.147 / 0.151 | 1.236, 1.025, 0.969 | 1.239, 1.026, 0.967 | 1/1 | 1/1 | 4/4 |
| 048 | 11 | desert | none | on | 400.0592 | 400.0806 | +0.0215 | 999.772 | 999.755 | -0.017 | 0.147 / 0.151 | 0.334, 0.010, 0.018 | 0.250, 0.008, 0.013 | 1/1 | 1/1 | 3/4 |
| 048 | 01 | desert | none | on | 400.2919 | 400.3125 | +0.0206 | 999.760 | 999.743 | -0.017 | 0.147 / 0.151 | 1.361, 0.981, 0.982 | 1.274, 0.979, 0.977 | 1/1 | 1/1 | 3/4 |
| 049 | 11 | desert | none | on | 401.3060 | 401.3267 | +0.0206 | 999.773 | 999.759 | -0.015 | 0.147 / 0.151 | 0.334, 0.010, 0.018 | 0.250, 0.008, 0.013 | 1/1 | 1/1 | 3/4 |
| 049 | 01 | desert | none | on | 401.3035 | 401.3367 | +0.0332 | 999.689 | 999.663 | -0.026 | 0.147 / 0.151 | 1.396, 0.974, 0.982 | 1.309, 0.973, 0.978 | 1/1 | 1/1 | 3/4 |
| 050 | 11 | desert | none | on | 402.5529 | 402.5728 | +0.0199 | 999.775 | 999.763 | -0.012 | 0.148 / 0.152 | 0.334, 0.011, 0.019 | 0.250, 0.009, 0.013 | 1/1 | 1/1 | 3/4 |
| 050 | 01 | desert | none | on | 402.5545 | 402.5504 | -0.0041 | 999.862 | 999.891 | +0.029 | 0.148 / 0.152 | 1.304, 1.065, 1.006 | 1.223, 1.062, 1.002 | 1/1 | 1/1 | 3/4 |
| 051 | 11 | desert | aer | off | 397.8428 | 396.9659 | -0.8769 | 1001.955 | 1002.801 | +0.847 | 0.168 / 0.196 | 0.131, 0.015, 0.019 | 0.133, 0.011, 0.034 | 1/1 | 1/1 | 4/4 |
| 051 | 01 | desert | aer | off | 397.8332 | 396.8785 | -0.9547 | 1002.137 | 1003.084 | +0.947 | 0.168 / 0.196 | 1.140, 1.075, 1.052 | 1.139, 1.071, 1.066 | 1/1 | 1/1 | 4/4 |
| 052 | 11 | desert | aer | off | 399.1028 | 398.2169 | -0.8859 | 1001.949 | 1002.820 | +0.871 | 0.169 / 0.196 | 0.131, 0.015, 0.019 | 0.133, 0.011, 0.032 | 1/1 | 1/1 | 4/4 |
| 052 | 01 | desert | aer | off | 399.1767 | 398.3335 | -0.8432 | 1001.991 | 1002.697 | +0.707 | 0.168 / 0.196 | 1.101, 0.978, 1.014 | 1.105, 0.975, 1.031 | 1/1 | 1/1 | 4/4 |
| 053 | 11 | desert | aer | off | 400.3629 | 399.4675 | -0.8954 | 1001.946 | 1002.839 | +0.894 | 0.169 / 0.196 | 0.131, 0.015, 0.019 | 0.133, 0.011, 0.031 | 1/1 | 1/1 | 4/4 |
| 053 | 01 | desert | aer | off | 400.6156 | 399.7529 | -0.8627 | 1001.961 | 1002.709 | +0.748 | 0.169 / 0.197 | 1.160, 0.928, 0.973 | 1.170, 0.923, 0.987 | 1/1 | 1/1 | 4/4 |
| 054 | 11 | desert | aer | off | 401.6234 | 400.7192 | -0.9042 | 1001.941 | 1002.857 | +0.916 | 0.170 / 0.196 | 0.131, 0.015, 0.019 | 0.132, 0.011, 0.030 | 1/1 | 1/1 | 4/4 |
| 054 | 01 | desert | aer | off | 401.3696 | 400.8080 | -0.5615 | 1001.989 | 1000.783 | -1.206 | 0.170 / 0.204 | 1.133, 0.973, 1.023 | 1.387, 0.971, 1.029 | 1/0 | 1/1 | 4/7 |
| 055 | 11 | desert | aer | off | 402.8842 | 402.3432 | -0.5410 | 1001.936 | 1000.690 | -1.246 | 0.171 / 0.204 | 0.132, 0.015, 0.019 | 0.385, 0.018, 0.025 | 1/0 | 1/1 | 4/7 |
| 055 | 01 | desert | aer | off | 402.6923 | 402.1825 | -0.5099 | 1001.756 | 1000.571 | -1.185 | 0.171 / 0.205 | 1.148, 1.091, 0.980 | 1.383, 1.078, 1.003 | 1/0 | 1/1 | 4/7 |
| 056 | 11 | desert | aer | on | 397.8429 | 396.9649 | -0.8780 | 1001.953 | 1002.802 | +0.850 | 0.168 / 0.196 | 0.131, 0.015, 0.019 | 0.133, 0.011, 0.034 | 1/1 | 1/1 | 4/4 |
| 056 | 01 | desert | aer | on | 397.5577 | 396.6761 | -0.8816 | 1002.056 | 1002.816 | +0.760 | 0.168 / 0.196 | 1.109, 0.989, 1.059 | 1.108, 0.986, 1.077 | 1/1 | 1/1 | 4/4 |
| 057 | 11 | desert | aer | on | 399.1033 | 398.2165 | -0.8868 | 1001.947 | 1002.819 | +0.872 | 0.169 / 0.196 | 0.131, 0.015, 0.019 | 0.132, 0.011, 0.032 | 1/1 | 1/1 | 4/4 |
| 057 | 01 | desert | aer | on | 398.9839 | 398.1504 | -0.8335 | 1001.854 | 1002.644 | +0.790 | 0.169 / 0.197 | 1.113, 0.973, 0.984 | 1.112, 0.967, 0.997 | 1/1 | 1/1 | 4/4 |
| 058 | 11 | desert | aer | on | 400.3633 | 399.4677 | -0.8956 | 1001.944 | 1002.840 | +0.896 | 0.169 / 0.196 | 0.131, 0.015, 0.019 | 0.132, 0.011, 0.031 | 1/1 | 1/1 | 4/4 |
| 058 | 01 | desert | aer | on | 400.4610 | 399.8640 | -0.5970 | 1002.030 | 1000.859 | -1.171 | 0.169 / 0.204 | 1.096, 1.006, 1.009 | 1.348, 1.004, 1.010 | 1/0 | 1/1 | 4/7 |
| 059 | 11 | desert | aer | on | 401.6228 | 400.7190 | -0.9038 | 1001.939 | 1002.856 | +0.917 | 0.170 / 0.196 | 0.131, 0.015, 0.019 | 0.132, 0.011, 0.030 | 1/1 | 1/1 | 4/4 |
| 059 | 01 | desert | aer | on | 401.7271 | 400.7135 | -1.0137 | 1002.212 | 1003.282 | +1.070 | 0.170 / 0.194 | 1.178, 1.007, 1.018 | 1.180, 1.000, 1.028 | 1/1 | 1/1 | 4/4 |
| 060 | 11 | desert | aer | on | 402.8830 | 402.3457 | -0.5373 | 1001.935 | 1000.684 | -1.251 | 0.171 / 0.204 | 0.131, 0.015, 0.019 | 0.385, 0.018, 0.025 | 1/0 | 1/1 | 4/7 |
| 060 | 01 | desert | aer | on | 402.9929 | 402.0650 | -0.9279 | 1002.077 | 1002.960 | +0.883 | 0.170 / 0.195 | 1.058, 0.982, 0.998 | 1.062, 0.975, 1.014 | 1/1 | 1/1 | 4/4 |
| 061 | 11 | forest | none | off | 397.9528 | 397.9744 | +0.0216 | 999.562 | 999.559 | -0.003 | 0.251 / 0.252 | 0.201, 0.003, 0.025 | 0.204, 0.003, 0.021 | 1/1 | 1/1 | 4/4 |
| 061 | 01 | forest | none | off | 397.6278 | 397.6629 | +0.0351 | 999.511 | 999.482 | -0.028 | 0.250 / 0.252 | 1.229, 1.035, 1.107 | 1.232, 1.035, 1.104 | 1/1 | 1/1 | 4/4 |
| 062 | 11 | forest | none | off | 399.2003 | 399.2218 | +0.0215 | 999.563 | 999.560 | -0.002 | 0.251 / 0.252 | 0.201, 0.003, 0.025 | 0.204, 0.003, 0.022 | 1/1 | 1/1 | 4/4 |
| 062 | 01 | forest | none | off | 399.2221 | 399.2515 | +0.0295 | 999.659 | 999.627 | -0.032 | 0.251 / 0.252 | 1.200, 1.002, 1.005 | 1.202, 1.002, 1.002 | 1/1 | 1/1 | 4/4 |
| 063 | 11 | forest | none | off | 400.4482 | 400.4692 | +0.0210 | 999.563 | 999.562 | -0.001 | 0.252 / 0.253 | 0.201, 0.003, 0.025 | 0.204, 0.004, 0.022 | 1/1 | 1/1 | 4/4 |
| 063 | 01 | forest | none | off | 400.3480 | 400.3500 | +0.0020 | 999.503 | 999.552 | +0.049 | 0.252 / 0.253 | 1.216, 0.934, 1.073 | 1.221, 0.935, 1.071 | 1/1 | 1/1 | 4/4 |
| 064 | 11 | forest | none | off | 401.6959 | 401.7168 | +0.0208 | 999.564 | 999.563 | -0.001 | 0.252 / 0.253 | 0.201, 0.003, 0.025 | 0.204, 0.004, 0.022 | 1/1 | 1/1 | 4/4 |
| 064 | 01 | forest | none | off | 401.3029 | 401.3093 | +0.0064 | 999.745 | 999.757 | +0.012 | 0.252 / 0.253 | 1.218, 0.994, 1.013 | 1.221, 0.994, 1.011 | 1/1 | 1/1 | 4/4 |
| 065 | 11 | forest | none | off | 402.9438 | 402.9642 | +0.0204 | 999.564 | 999.564 | -0.000 | 0.253 / 0.254 | 0.200, 0.004, 0.025 | 0.203, 0.004, 0.022 | 1/1 | 1/1 | 4/4 |
| 065 | 01 | forest | none | off | 403.1419 | 403.1576 | +0.0157 | 999.551 | 999.563 | +0.012 | 0.253 / 0.254 | 1.173, 1.000, 1.028 | 1.178, 1.000, 1.026 | 1/1 | 1/1 | 4/4 |
| 066 | 11 | forest | none | on | 397.9542 | 397.9760 | +0.0218 | 999.559 | 999.556 | -0.003 | 0.251 / 0.252 | 0.201, 0.003, 0.025 | 0.204, 0.003, 0.021 | 1/1 | 1/1 | 4/4 |
| 066 | 01 | forest | none | on | 397.9125 | 397.9483 | +0.0357 | 999.553 | 999.519 | -0.034 | 0.250 / 0.252 | 1.213, 1.039, 0.982 | 1.216, 1.039, 0.977 | 1/1 | 1/1 | 4/4 |
| 067 | 11 | forest | none | on | 399.2019 | 399.2233 | +0.0215 | 999.559 | 999.557 | -0.002 | 0.251 / 0.252 | 0.201, 0.003, 0.025 | 0.204, 0.003, 0.022 | 1/1 | 1/1 | 4/4 |
| 067 | 01 | forest | none | on | 398.8244 | 398.8501 | +0.0257 | 999.812 | 999.771 | -0.040 | 0.251 / 0.252 | 1.177, 1.014, 1.037 | 1.179, 1.014, 1.035 | 1/1 | 1/1 | 4/4 |
| 068 | 11 | forest | none | on | 400.4497 | 400.4708 | +0.0212 | 999.560 | 999.559 | -0.001 | 0.252 / 0.253 | 0.200, 0.003, 0.025 | 0.204, 0.004, 0.022 | 1/1 | 1/1 | 4/4 |
| 068 | 01 | forest | none | on | 400.3997 | 400.4366 | +0.0369 | 999.527 | 999.493 | -0.034 | 0.252 / 0.253 | 1.134, 0.979, 0.974 | 1.136, 0.979, 0.970 | 1/1 | 1/1 | 4/4 |
| 069 | 11 | forest | none | on | 401.6974 | 401.7183 | +0.0209 | 999.561 | 999.560 | -0.001 | 0.252 / 0.253 | 0.200, 0.003, 0.025 | 0.203, 0.004, 0.022 | 1/1 | 1/1 | 4/4 |
| 069 | 01 | forest | none | on | 401.7520 | 401.7791 | +0.0271 | 999.442 | 999.441 | -0.001 | 0.252 / 0.254 | 1.199, 0.918, 1.035 | 1.203, 0.918, 1.031 | 1/1 | 1/1 | 4/4 |
| 070 | 11 | forest | none | on | 402.9453 | 402.9657 | +0.0205 | 999.561 | 999.561 | -0.000 | 0.253 / 0.254 | 0.200, 0.004, 0.025 | 0.203, 0.004, 0.022 | 1/1 | 1/1 | 4/4 |
| 070 | 01 | forest | none | on | 403.0712 | 403.0802 | +0.0090 | 999.576 | 999.601 | +0.025 | 0.253 / 0.254 | 1.238, 0.961, 0.972 | 1.242, 0.961, 0.970 | 1/1 | 1/1 | 4/4 |
| 071 | 11 | forest | aer | off | 397.3837 | 396.9654 | -0.4183 | 1002.579 | 1002.874 | +0.294 | 0.326 / 0.359 | 0.139, 0.002, 0.042 | 0.087, 0.002, 0.022 | 1/1 | 1/1 | 4/5 |
| 071 | 01 | forest | aer | off | 397.1625 | 397.2954 | +0.1330 | 1003.379 | 1002.147 | -1.232 | 0.326 / 0.360 | 1.092, 0.919, 1.036 | 1.008, 0.918, 1.000 | 1/1 | 1/1 | 4/5 |
| 072 | 11 | forest | aer | off | 398.6627 | 398.2113 | -0.4515 | 1002.589 | 1002.876 | +0.288 | 0.328 / 0.360 | 0.139, 0.002, 0.042 | 0.087, 0.002, 0.022 | 1/1 | 1/1 | 4/5 |
| 072 | 01 | forest | aer | off | 399.2775 | 397.9650 | -1.3125 | 1001.134 | 1003.655 | +2.522 | 0.328 / 0.359 | 1.138, 1.021, 0.962 | 1.183, 1.022, 0.993 | 1/1 | 1/1 | 5/4 |
| 073 | 11 | forest | aer | off | 399.9423 | 399.4573 | -0.4850 | 1002.597 | 1002.881 | +0.283 | 0.329 / 0.361 | 0.138, 0.002, 0.043 | 0.087, 0.002, 0.022 | 1/1 | 1/1 | 4/5 |
| 073 | 01 | forest | aer | off | 400.0693 | 399.4731 | -0.5962 | 1002.972 | 1003.495 | +0.522 | 0.329 / 0.362 | 1.127, 1.056, 0.997 | 1.065, 1.057, 0.977 | 1/1 | 1/1 | 4/5 |
| 074 | 11 | forest | aer | off | 401.2218 | 400.7019 | -0.5199 | 1002.605 | 1002.890 | +0.285 | 0.330 / 0.363 | 0.139, 0.002, 0.043 | 0.087, 0.002, 0.022 | 1/1 | 1/1 | 4/5 |
| 074 | 01 | forest | aer | off | 402.0258 | 400.8364 | -1.1893 | 1001.612 | 1003.766 | +2.155 | 0.326 / 0.361 | 1.091, 0.962, 0.958 | 1.143, 0.965, 0.976 | 1/1 | 1/1 | 5/4 |
| 075 | 11 | forest | aer | off | 402.5005 | 401.6172 | -0.8832 | 1002.614 | 1003.955 | +1.340 | 0.331 / 0.363 | 0.138, 0.002, 0.043 | 0.132, 0.004, 0.039 | 1/1 | 1/1 | 4/4 |
| 075 | 01 | forest | aer | off | 402.4683 | 401.6499 | -0.8184 | 1002.593 | 1003.807 | +1.214 | 0.331 / 0.363 | 1.138, 1.037, 1.061 | 1.135, 1.038, 1.057 | 1/1 | 1/1 | 4/4 |
| 076 | 11 | forest | aer | on | 397.3839 | 396.9691 | -0.4147 | 1002.579 | 1002.866 | +0.288 | 0.326 / 0.359 | 0.139, 0.002, 0.042 | 0.087, 0.002, 0.022 | 1/1 | 1/1 | 4/5 |
| 076 | 01 | forest | aer | on | 397.6977 | 397.3278 | -0.3698 | 1002.595 | 1003.025 | +0.430 | 0.327 / 0.358 | 1.147, 0.898, 1.048 | 1.103, 0.898, 1.023 | 1/1 | 1/1 | 4/5 |
| 077 | 11 | forest | aer | on | 398.6632 | 398.2115 | -0.4517 | 1002.587 | 1002.876 | +0.289 | 0.328 / 0.360 | 0.139, 0.002, 0.042 | 0.087, 0.002, 0.022 | 1/1 | 1/1 | 4/5 |
| 077 | 01 | forest | aer | on | 399.1671 | 397.8552 | -1.3119 | 1001.232 | 1003.708 | +2.477 | 0.329 / 0.360 | 1.133, 0.994, 1.031 | 1.173, 0.997, 1.058 | 1/1 | 1/1 | 5/4 |
| 078 | 11 | forest | aer | on | 399.9420 | 399.4586 | -0.4833 | 1002.597 | 1002.878 | +0.281 | 0.329 / 0.361 | 0.138, 0.002, 0.043 | 0.087, 0.002, 0.022 | 1/1 | 1/1 | 4/5 |
| 078 | 01 | forest | aer | on | 400.9328 | 399.4093 | -1.5235 | 1000.324 | 1003.482 | +3.157 | 0.326 / 0.359 | 1.031, 0.992, 1.000 | 1.125, 0.994, 1.047 | 1/1 | 1/1 | 6/4 |
| 079 | 11 | forest | aer | on | 401.2217 | 400.7017 | -0.5199 | 1002.605 | 1002.888 | +0.283 | 0.330 / 0.363 | 0.138, 0.002, 0.043 | 0.087, 0.002, 0.022 | 1/1 | 1/1 | 4/5 |
| 079 | 01 | forest | aer | on | 402.1650 | 401.3750 | -0.7901 | 1001.761 | 1002.884 | +1.122 | 0.327 / 0.362 | 1.082, 0.980, 0.973 | 1.080, 0.982, 0.975 | 1/1 | 1/1 | 5/5 |
| 080 | 11 | forest | aer | on | 402.5011 | 401.6168 | -0.8844 | 1002.613 | 1003.955 | +1.342 | 0.331 / 0.363 | 0.138, 0.002, 0.043 | 0.132, 0.004, 0.039 | 1/1 | 1/1 | 4/4 |
| 080 | 01 | forest | aer | on | 402.7511 | 401.7607 | -0.9904 | 1001.604 | 1003.040 | +1.436 | 0.330 / 0.364 | 1.126, 1.018, 1.056 | 1.111, 1.014, 1.064 | 1/1 | 1/1 | 5/5 |
