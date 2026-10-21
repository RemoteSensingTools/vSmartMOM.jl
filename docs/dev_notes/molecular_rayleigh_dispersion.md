# Molecular Rayleigh: N₂/O₂ dispersion, tensor anisotropy, H₂ check (v2.2)

Implements decisions 13a–13e of 2026-10-19 (`V2_2_DECISION_LOG.md` §1.8).
Code: `src/Scattering/molecular_rayleigh.jl` (fits, tensor types),
`src/Inelastic/src/inelastic_cross_section.jl` (`compute_effective_coefficents!`
dispatch), `src/Inelastic/inelastic_helper.jl` (`getRamanAtmoConstants`,
`rayleigh_air_depolarization`), `src/CoreRT/tools/model_from_parameters.jl`
(`_auto_molecular_depolarizations`, used by the forward, linearized and batch
builders). Tests: `test/test_molecular_rayleigh.jl`,
`test/test_rayleigh_anisotropy.jl`.

Units: ν in cm⁻¹, s = 1/λ in µm⁻¹, refractivity r = (n − 1)·10⁸. All
wavelengths are vacuum wavelengths. Cross section of a fit given at number
density N:

σ = 24π³ν⁴ [(n² − 1)/(n² + 2)]² F_K / N²,  ρ = 6(F_K − 1)/(7F_K + 3).

## 13a. O₂ King factor: 1.09 (He et al. 2021), kept

F_K(O₂) = 1.09 + 1.385×10⁻¹¹ν² + 1.448×10⁻²⁰ν⁴ (He et al. 2021, ACP 21, 14927,
Table 1, verified from the article's Table 1 XLSX). Bates (1984) and Bodhaine et
al. (1999) use the constant 1.096 with the same dispersion terms: +0.543 % in
σ_O₂ at 400 nm, +0.550 % at 1.6 µm. Cited in the `UsingRefInd` docstring and
the King-factor source comment.

Independent evidence, for the record (not acted on): the two-oscillator
anisotropy dispersion of Hohm (1994, Chem. Phys. 179, 533) for O₂, as
transcribed in PICASO (`picaso/rayleigh.py`; primary not accessible here),
gives F_K = 1.1028 at 550 nm and 1.0949 at 1.6 µm, i.e. γ within ±0.9 % of the
Bates-1.096 anisotropy over 400–2060 nm and 2.3–3.7 % above the He-1.09 one.
The same transcription reproduces the H₂ ab initio King factor (below) to
7×10⁻⁵, which supports it. He et al. themselves found measured O₂ cross
sections 2.8 ± 1.2 % above the n-based values.

## 13b. Refractive-index fits valid through the NIR

| Gas | Fit | Formula (r = (n−1)·10⁸) | State, N (cm⁻³) | Range | Used |
|---|---|---|---|---|---|
| N₂ | Bates (1984) | 5677.465 + 318.81874×10¹²/(1.44×10¹⁰ − ν²) | 288.15 K, 1013.25 hPa; 2.546899×10¹⁹ | 254–468 nm (fit) | λ ≤ 436.309 nm |
| N₂ | Peck & Khanna (1966) | 6855.200 + 324.3157×10¹²/(1.44×10¹⁰ − ν²) = 10⁸[6.8552×10⁻⁵ + 3.243157×10⁻²/(144 − s²)] | 0 °C, 760 Torr; 2.68678×10¹⁹ | 467.9–2058.7 nm (measured) | λ > 436.309 nm |
| O₂ | Bates (1984) | 20564.8 + 2.480899×10¹³/(4.09×10⁹ − ν²) | 273.15 K, 1013.25 hPa; 2.68678×10¹⁹ | 288–546 nm (fit) | λ ≤ 466.678 nm |
| O₂ | Egan (2026), Eq. (8) | 5627886.8/(281.0 − s²) + 197258.6/(45.4 − s²) | 20 °C, 100 kPa; 1×10⁵ Pa/(k_B·293.15 K) = 2.470739×10¹⁹ | constrained 0.254–1.542 µm | λ > 466.678 nm |

Sources and verification:

- Bates (1984) formulas: He et al. (2021) Table 1 (XLSX), which also gives the
  ν ranges and marks the O₂ fit as 273.15 K / 2.68678×10¹⁹.
- Peck & Khanna (1966), JOSA 56, 1059 (paywalled here). Formula checked
  against the refractiveindex.info records `N2/nk/Peck-0C.yml` (0 °C, the form
  used) and `Peck-15C.yml` (15 °C: 6.497378×10⁻⁵ + 3.0738649×10⁻²/(144 − s²)),
  database commit `37773dce8128ccabcd9d5108d9958d6041244661`. The 15 °C record
  is the 0 °C one × 0.947802 (real-gas density ratio), 1.4×10⁻⁴ below the
  ideal-gas ratio 273.15/288.15; the 15 °C form of Sneep & Ubachs (2005,
  Sec. 4.2) and Bates, 6498.2 + 307.43305×10¹²/(1.44×10¹⁰ − ν²), is the 0 °C
  formula scaled by the ideal-gas ratio (agreement 4×10⁻⁶; seen in the PICASO
  and TauREx-2 transcriptions). We use the 0 °C form with the ideal-gas
  Loschmidt density, which is the same convention.
- Egan (2026), Meas. Sci. Technol. 37, 255003, CC BY 4.0
  (https://doi.org/10.1088/1361-6501/ae6e94); Eq. (8) and Secs. 4.2–4.5 read
  from the published PDF. Data behind Eq. (8): Ladenburg & Wolfsohn (1932)
  UV–visible (10⁻³ uncertainty), Hohm visible points, the author's two-colour
  refractivity at 0.6329908 and 1.542383 µm (relative uncertainty 7.2×10⁻⁶),
  and a 0.436/0.254 µm ratio (Keeling et al.). Egan states Eq. (8) "should be
  accurate throughout the optical frequency range"; the 0.76 µm A band perturbs
  n − 1 by ~10⁻⁵ locally. Test: our ideal-gas per-molecule polarizability from
  Eq. (8) equals Egan's measured molar refractivity A_R(293.15 K) at both
  wavelengths times (1 + b_Rρ)/Z with B(293.15 K) = −17.0 cm³/mol (his
  Table 2), to 5×10⁻⁶.
- O₂ alternatives not used: Zhang et al. (2008, Appl. Opt. 47, 3143) — Egan
  finds their O₂ value at 0.8 µm 4.5×10⁻⁴ (n−1) low, > 5× the mutual expanded
  uncertainty; the refractiveindex.info `O2/nk/Zhang.yml` record
  (1.181494×10⁻⁴ + 9.708931×10⁻³/(75.4 − s²), 20 °C) is Křen's (2011) inferred
  single-oscillator formula, which agrees with Eq. (8) to 6×10⁻⁵ at 633 nm and
  1.5×10⁻⁴ at 1.54 µm (Egan quotes 8×10⁻⁵ and 1.7×10⁻⁴), i.e. our Egan
  implementation reproduces Egan's own comparison. Smith et al. (1976, Phys.
  Rev. A 13, 1422) cover only 168–288 nm.
- Number densities are ideal-gas values at each fit's reference state, as in
  Bates (1984), Bodhaine et al. (1999), Sneep & Ubachs (2005) and He et al.
  (2021). Compressibility makes the per-molecule polarizability high by
  1 − Z ≈ 4.6×10⁻⁴ (N₂, 0 °C), 9.8×10⁻⁴ (O₂, 0 °C) and 7.0×10⁻⁴ (O₂, 20 °C,
  100 kPa); σ by twice that. Not corrected (convention of all compilations).

Join. The previous fits ended at 468 nm (N₂) and 546 nm (O₂). Joined at those
limits they would jump by:

| | n − 1 (old/new − 1) | σ (old/new − 1) |
|---|---|---|
| N₂ at 21 360 cm⁻¹ (468.2 nm) | −1.42×10⁻⁴ | −2.85×10⁻⁴ |
| O₂ at 18 315 cm⁻¹ (546.0 nm) | +1.04×10⁻³ | +2.08×10⁻³ |

The O₂ jump is not negligible (it exceeds the Egan and Peck–Khanna precision by
two orders of magnitude and would sit inside the 466–546 nm region). Both
gases are therefore switched where the two fits give the same per-molecule
Lorentz–Lorenz volume [(n²−1)/(n²+2)]/N (bisection): ν× = 22 919.551 cm⁻¹
(436.309 nm) for N₂ and 21 428.071 cm⁻¹ (466.678 nm) for O₂. n − 1 and σ are
continuous (residual mismatch at the rounded joins 6×10⁻¹² and 8×10⁻¹¹); only
the slope has a kink. The O₂ join lies inside the validity of both fits. The
N₂ join uses Peck & Khanna 31 nm below their shortest line (467.9 nm), where it
differs from the Bates fit by < 1.4×10⁻⁴ in n − 1; the alternative (joining
inside a common range) does not exist for N₂. The Bates fits are kept unchanged
below the joins, so all values at λ ≤ 436 nm (N₂) and ≤ 466 nm (O₂) are
bit-identical to the previous fits (the 360 and 400 nm regression anchors
still pass).

`refractive_index` is still reported at the species reference density
(N₂ 2.546899×10¹⁹, O₂ 2.68678×10¹⁹); a fit at another density is converted
through L ∝ N.

## 13c. Wavelength-dependent anisotropy γ(λ) and mean polarizability ᾱ(λ)

Two valid formulations, selected by dispatch (`AbstractPolarizabilityDispersion`):

- `ConstantAnisotropy()` — the Raman-module tensor: Buldakov
  ᾱ(ν,T) = ᾱ₀₀(1 + α_bT + α_cT²)/(1 − (2πcν/ω₀)²) and γ = γ̄₀₀ (N₂ 0.71,
  O₂ 1.080, H₂ 0.288 ×10⁻²⁴ cm³). Default of the Raman and Cabannes cross
  sections (`getRamanAtmoConstants(ν, T)`, the published Sanghavi 2022 and
  Sanghavi 2026 Table 1 values; `test_ray_cab_raman.jl`,
  `test_raman_placzek_teller.jl` unchanged).
- `DispersiveAnisotropy()` — default of `UsingPolTen()` and of the automatic
  `depol: -1` Rayleigh depolarization:
  ᾱ(ν,T) = (3/4π)·L(ν)/N · f(T)/f(T_fit), γ(ν) = (3/4π)·L(ν)/N·√(9(F_K(ν) − 1)/2),
  with L, N, T_fit of the 13b fit used at ν, F_K of `UsingRefInd`, and the
  Buldakov factor f(T) = 1 + α_bT + α_cT² (Egan measures the O₂ temperature
  coefficient 2.11×10⁻⁶ K⁻¹ at 633 nm; Buldakov's f gives 2.6×10⁻⁶ K⁻¹ near
  288 K). γ is the anisotropy implied by the measured depolarization ratios
  behind the King-factor fits (Bates 1984; as in Bodhaine et al. 1999 and
  Tomasi et al. 2005, Appl. Opt. 44, 3320), and ᾱ is Lorentz–Lorenz-consistent
  with the RI fits by construction. At T = T_fit the cross section and ρ equal
  `UsingRefInd` exactly; for 200–300 K they differ by < 4×10⁻⁴ (largest for
  O₂ at 200 K: −3.3×10⁻⁴ in σ, +3.3×10⁻⁴ in ρ). Validity: as
  the fits (N₂ 254–2059 nm, O₂ 288–1542 nm), with the King factors fitted to
  visible depolarization data and extrapolated elsewhere. At T ≠ T_fit the
  different reference temperatures of the two fits leave steps ≤ 1×10⁻⁴ in σ at
  the joins. H₂: no `UsingRefInd` fit; keeps the Raman-module tensor. The
  vibrational-Raman derivatives α′, γ′ stay the Raman-module constants.

Cross-checks of γ(λ): the constant γ̄₀₀ equals the King-factor anisotropy at
488 nm (Ar⁺ line) to within 0.03 % (N₂) and 0.32 % (O₂), the King-factor value
being the lower one; relative to the constant,
γ_disp is +2.1 % (N₂) / +4.2 % (O₂) at 400 nm and −3.9 % / −7.2 % at 2060 nm.
The dispersive ρ falls with λ, the constant-γ ρ rises.

| γ (10⁻²⁴ cm³) | 400 | 488 | 550 | 760 | 1600 | 2060 nm |
|---|---|---|---|---|---|---|
| N₂ dispersive | 0.7247 | 0.7098 | 0.7035 | 0.6925 | 0.6833 | 0.6823 |
| O₂ dispersive (He 1.09) | 1.1252 | 1.0765 | 1.0573 | 1.0270 | 1.0049 | 1.0026 |
| O₂ with Bates 1.096 | 1.1571 | 1.1089 | 1.0898 | 1.0597 | 1.0377 | 1.0353 |
| O₂ Hohm (1994) via PICASO | 1.1668 | 1.1161 | 1.0949 | 1.0588 | 1.0292 | 1.0259 |

ᾱ consistency: over 400–2060 nm the Buldakov ᾱ of the constant formulation
differs from the Lorentz–Lorenz ᾱ of the joined fits by −3.7×10⁻⁴ to
+7.8×10⁻⁴ (N₂) and +1×10⁻⁴ to −2.5×10⁻³ (O₂, ideal-gas basis; −0.4×10⁻³ to
−1.7×10⁻³ against Egan's real-gas polarizability).

Selection in the RT model: `radiative_transfer.rayleigh_anisotropy`
(`"DispersiveAnisotropy()"` default, `"ConstantAnisotropy()"`), field
`vSmartMOM_Parameters.rayleigh_anisotropy` (the positional constructor
without the new field keeps working and selects the default). With `depol: -1` the
Rayleigh depolarization (greek_rayleigh, King rescaling of τ_rayl) uses it;
the Cabannes depolarization (greek_cabannes) and ϖ_Cabannes always use the
Raman-module constants. `ConstantAnisotropy()` reproduces the v2.2
pre-release values bit for bit. It does not reproduce v2.1.0, which evaluated
the constants at 300 K and left 2π out of the Buldakov dispersion pole
(changed in `817d93cb`; air ρ about 0.0281 from 400 nm to 2 µm).

## 13d. `UsingPolTen` default temperature 288.15 K, kept

It is the temperature of the default number density 2.546899×10¹⁹ cm⁻³ and of
the standard state of the He, Ar, N₂ (short-wavelength) and CO₂ RI fits and of
the Bates/Bodhaine cross sections, so `UsingPolTen` is compared with them at
their own state. Per-molecule results depend on T only through f(T):
200→300 K changes σ by 3.5×10⁻⁴ (N₂), 3.6×10⁻⁴ (O₂), 1.9×10⁻³ (H₂). The RT
model never uses this default (it takes the column-weighted mean T).

## 13e. H₂ check

Refractivity: Peck & Huang (1977), JOSA 67, 1550 (authors E. R. Peck and
Shiezen Huang per the Crossref/Optica record; refractiveindex.info spells the
second author "S. Hung"), 0 °C/760 Torr, 168–1695 nm:
n − 1 = 0.0148956/(180.7 − s²) + 0.0049037/(92 − s²) (refractiveindex.info
`H2/nk/Peck.yml`; same coefficients in an independent transcription,
github.com/julie-dujardin/space-map). King factor: Raj, Hamaguchi & Witek
(2018), J. Chem. Phys. 148, 104308 — ab initio ⟨v=0,J|ᾱ|v=0,J⟩ and ⟨γ⟩ obtained
here by integrating their published r- and λ-dependent polarizabilities over
the v=0 wavefunctions (data repository github.com/ankit7540/
H2-PolarizabilityMatrixElements; own integration code, not theirs). J = 1
(most populated near 288 K) is used; J = 0 gives F_K lower by 1.4×10⁻⁴.

| λ (nm) | ᾱ_PH/ᾱ_Raj(J=1) − 1 | F_K Raj | F_K Raman module | ρ module/Raj − 1 | σ_PolTen/σ(PH+F_Raj) − 1 |
|---|---|---|---|---|---|
| 400 | −3.2×10⁻⁴ | 1.03328 | 1.02575 | −22 % | −0.37 % |
| 550 | −3.2×10⁻⁴ | 1.03235 | 1.02700 | −16 % | −0.19 % |
| 760 | −3.2×10⁻⁴ | 1.03187 | 1.02768 | −13 % | −0.11 % |
| 1064 | −3.2×10⁻⁴ | 1.03162 | 1.02805 | −11 % | −0.07 % |

The Peck & Huang and ab initio mean polarizabilities agree to 3×10⁻⁴; the
Raman-module H₂ ᾱ is 0.15 % above Peck & Huang; the Raman-module H₂ anisotropy
0.288×10⁻²⁴ cm³ (Asawaroengchai & Rosenblatt 1980) is 4.7 % below the static
and 8.5 % below the 550 nm ab initio value, so its F_K is 0.35–0.73 % low and ρ 11–22 %
low. The test requires σ within 0.5 %, ᾱ_PH vs ab initio within 10⁻³.

## Effect on results (before = v2.2 branch base f3ae5bdb, after = this change)

### `UsingRefInd` (n − 1 at the species reference density; ρ unchanged)

| | λ (nm) | n − 1 before | n − 1 after | σ before (cm²) | σ after (cm²) | Δσ | ρ |
|---|---|---|---|---|---|---|---|
| N₂ | 400 | 2.8822201e-4 | same | 1.7138759e-26 | same | 0 | 0.021058 |
| N₂ | 550 | 2.8337867e-4 | 2.8349228e-4 | 4.6308196e-27 | 4.6345334e-27 | +0.080 % | 0.020525 |
| N₂ | 760 | 2.8087085e-4 | 2.8107405e-4 | 1.2471657e-27 | 1.2489709e-27 | +0.145 % | 0.020240 |
| N₂ | 1600 | 2.7877878e-4 | 2.7905672e-4 | 6.2521092e-29 | 6.2645813e-29 | +0.199 % | 0.019997 |
| N₂ | 2060 | 2.7853946e-4 | 2.7882595e-4 | 2.2712819e-29 | 2.2759563e-29 | +0.206 % | 0.019969 |
| O₂ | 400 | 2.7724682e-4 | same | 1.5190048e-26 | same | 0 | 0.058328 |
| O₂ | 550 | 2.7163951e-4 | 2.7134565e-4 | 4.0493398e-27 | 4.0405836e-27 | −0.216 % | 0.054058 |
| O₂ | 760 | 2.6898682e-4 | 2.6826579e-4 | 1.0857658e-27 | 1.0799530e-27 | −0.535 % | 0.052301 |
| O₂ | 1600 | 2.6689059e-4 | 2.6575560e-4 | 5.4301630e-29 | 5.3840782e-29 | −0.849 % | 0.051099 |
| O₂ | 2060 | 2.6665719e-4 | 2.6547168e-4 | 1.9722998e-29 | 1.9548026e-29 | −0.887 % | 0.050977 |

Dry air from the `UsingRefInd` gases (Bodhaine composition: N₂ 78.084 %, O₂
20.946 %, Ar 0.934 %, CO₂ 360 ppm; ρ from the isotropic/anisotropic sums):

| λ | n − 1 after | vs Peck–Reeder before → after | vs Ciddor before → after | σ after | vs Bodhaine Eq. 29 before → after | ρ before → after |
|---|---|---|---|---|---|---|
| 400 | 2.8280116e-4 | +1.86e-4 → same | +1.39e-4 → same | 1.6717605e-26 | −0.125 % → same | 0.027927 |
| 550 | 2.7789269e-4 | +1.42e-4 → +2.51e-4 | +0.89e-4 → +1.98e-4 | 4.5066343e-27 | −0.109 % → −0.085 % | 0.026634 → 0.026619 |
| 760 | 2.7537034e-4 | +2.00e-4 → +2.56e-4 | +1.48e-4 → +2.04e-4 | 1.2126280e-27 | −0.084 % → −0.068 % | 0.026058 → 0.026025 |
| 1600 | 2.7327788e-4 | +2.87e-4 → +2.56e-4 | +2.39e-4 → +2.09e-4 | 6.0754360e-29 | −1.10 % → −1.10 % | 0.025629 → 0.025580 |
| 2060 | 2.7303916e-4 | +2.98e-4 → +2.56e-4 | +2.52e-4 → +2.09e-4 | 2.2069744e-29 | −3.50 % → −3.50 % | 0.025583 → 0.025532 |

With the new fits the air refractivity departs from Peck & Reeder (1972) and
Ciddor (1996) by a nearly constant 2.51–2.56×10⁻⁴ and 1.98–2.09×10⁻⁴ from 550
to 2060 nm (before: 1.4–3.0×10⁻⁴ and 0.9–2.5×10⁻⁴, growing with λ), i.e. the
dispersion shape now matches both independent air formulas (the constant
offset is the CO₂ content and compressibility convention). Bodhaine's Eq. 29 is a fit valid to ~1 µm: it
reproduces Peck–Reeder + Bodhaine's King factor to ≤ 7×10⁻⁵ over 250–760 nm,
but is +0.05 % (1.0 µm), +0.29 % (1.27 µm), +1.05 % (1.6 µm), +3.6 % (2.06 µm)
and +5.9 % (2.3 µm) high beyond it. The RT model's τ_rayl uses this fit (Eq. 30);
see open question 3.

### `UsingPolTen()` default (288.15 K): constant → dispersive

| | λ | σ constant | σ dispersive | ρ constant | ρ dispersive |
|---|---|---|---|---|---|
| N₂ | 400 | 1.7102660e-26 | 1.7138759e-26 (+0.21 %) | 0.020247 | 0.021058 |
| N₂ | 550 | 4.6409084e-27 | 4.6347750e-27 (−0.13 %) | 0.020881 | 0.020524 |
| N₂ | 760 | 1.2525770e-27 | 1.2490360e-27 (−0.28 %) | 0.021224 | 0.020239 |
| N₂ | 1600 | 6.2903021e-29 | 6.2649083e-29 (−0.40 %) | 0.021517 | 0.019996 |
| N₂ | 2060 | 2.2856123e-29 | 2.2760751e-29 (−0.42 %) | 0.021551 | 0.019967 |
| O₂ | 400 | 1.5081644e-26 | 1.5191082e-26 (+0.73 %) | 0.054014 | 0.058324 |
| O₂ | 550 | 4.0477317e-27 | 4.0404848e-27 (−0.18 %) | 0.056368 | 0.054059 |
| O₂ | 760 | 1.0863358e-27 | 1.0799265e-27 (−0.59 %) | 0.057643 | 0.052302 |
| O₂ | 1600 | 5.4300210e-29 | 5.3839459e-29 (−0.85 %) | 0.058737 | 0.051100 |
| O₂ | 2060 | 1.9719740e-29 | 1.9547546e-29 (−0.87 %) | 0.058864 | 0.050978 |

H₂ unchanged (σ 3.7446817e-27 at 400 nm).

### `depol: -1` in the RT model (air VMR 0.8/0.2, column-mean T = 257.1 K of `PureRayleighParameters.yaml`)

| λ | ρ_Rayleigh before | after | τ_R(1013.25 hPa) before | after | Δτ |
|---|---|---|---|---|---|
| 400 | 0.026267 | 0.027721 | 0.3592161 | 0.3601038 | +0.247 % |
| 550 | 0.027146 | 0.026447 | 0.0969411 | 0.0968261 | −0.119 % |
| 760 | 0.027621 | 0.025867 | 0.0261011 | 0.0260235 | −0.297 % |
| 1600 | 0.028029 | 0.025430 | 0.0013223 | 0.0013165 | −0.44 % |
| 2060 | 0.028076 | 0.025383 | 0.0004923 | 0.0004901 | −0.45 % |

ρ_Cabannes (0.006974–0.007641) and ϖ_Cabannes (0.96800–0.96611) are
unchanged. τ_R is Bodhaine Eq. 30 rescaled by F(ρ)/F(0.0279); for reference
Bodhaine Eq. 30 alone gives 0.3602134, 0.0970652, 0.0261134, 0.0013220,
0.0004922. The τ change follows the King factor F(ρ) only; the Rayleigh phase
matrix changes through ρ: ρ changes by +5.5 %, −2.6 %, −6.4 %, −9.3 % and −9.6 %
at the five wavelengths, the l = 2 Greek coefficients (∝ (1 − ρ)/(1 + ρ/2)) by
−0.22 %, +0.11 %, +0.27 %, +0.40 % and +0.41 %.

### Raman cross sections if they also switched (report only; not changed)

`getRamanAtmoConstants(ν, 250 K, DispersiveAnisotropy())` vs the default
(constant), air VMR 0.8/0.2:

| λ | σ_RRS constant | σ_RRS dispersive | Δσ_RRS | ϖ_Cab constant | ϖ_Cab dispersive | Δϖ_Cab | ρ_Cab const → disp |
|---|---|---|---|---|---|---|---|
| 395 | 5.6175556e-28 | 5.9701767e-28 | +6.28 % | 0.96806196 | 0.96617027 | −1.89×10⁻³ | 0.006965 → 0.007395 |
| 550 | 1.4931945e-28 | 1.4532528e-28 | −2.67 % | 0.96697679 | 0.96781489 | +8.4×10⁻⁴ | 0.007234 → 0.007042 |
| 760 | 4.0909401e-29 | 3.8217505e-29 | −6.58 % | 0.96644464 | 0.96854647 | +2.10×10⁻³ | 0.007387 → 0.006906 |
| 1600 | 2.0737282e-30 | 1.8749607e-30 | −9.58 % | 0.96610119 | 0.96920173 | +3.10×10⁻³ | 0.007588 → 0.006867 |

σ_RRS ∝ γ², so it moves by about twice the γ change; the Rayleigh total moves
by +0.33 % (395 nm) to −0.48 % (1600 nm). In a Raman run the RRS optical depth
is τ_R(1 − ϖ_Cab): −6 % at 760 nm and +6 % at 395 nm, which directly scales
Fraunhofer/O₂-A filling-in.

### Elastic-only vs Raman runs under the default

With `depol: -1` and the dispersive default, a pure-elastic (`noRS`) run uses
`greek_rayleigh` from the dispersive ρ, while a Raman run (`RRS`, `VS_*`) uses
`greek_cabannes` and ϖ_Cabannes from the constant Raman-module tensor plus RRS
lines from the same constant tensor. The elastic-only phase matrix is then no
longer ϖ_Cab·Z_Cab + (1 − ϖ_Cab)·Z_RRS of the Raman run: the implied ρ
differs by +5.5 % (400 nm), −6 % (760 nm) and −9 % (1.6 µm), about 0.2–0.4 %
in the l = 2 Greek coefficients. Both runs use the same τ_rayl. Workflows that
difference a `noRS` and an RRS run (Ring/filling-in spectra, OCORaman,
RRS_XCO2 elastic references) should set `rayleigh_anisotropy:
ConstantAnisotropy()` until the Raman cross sections are decided (open
question 1).

## Tests that move

- `test_vertical_grid_psurf.jl`, "default convention is FixedPressureGrid,
  bit for bit": the golden file (7acab570) was made with the constant
  anisotropy; with the dispersive default τ_rayl at 1613 nm moves by −0.44 %
  (computed for the fixture profile; it would fail the rtol 1e-12 check). Handling: the test now sets
  `params.rayleigh_anisotropy = ConstantAnisotropy()`; golden data untouched.
- `test_forward_raman_phase1b.jl`: its RRS reference
  (`phase1b_RRS_sanghavi_q0.jld2`) predates the dispersive anisotropy; with
  the default, the `depol: -1` τ_rayl at 763 nm moves by −0.3 % (inside its
  rtol 0.02, but a mixed comparison). Handling: the test sets
  `params.rayleigh_anisotropy = ConstantAnisotropy()`; reference untouched.
- `test_molecular_rayleigh.jl`: the checks that are specific to the
  Raman-module tensor (equality with `getRamanAtmoConstants`, closed-form
  anchors, 400 nm pinned σ/ρ, intercomparison) now name
  `UsingPolTen(ConstantAnisotropy())`; values unchanged. The intercomparison
  comment numbers are updated for the new RI fits (thresholds unchanged).

## Downstream

- ExoOptics `exo_Earth/src/AToD/OneDimensionalDiagnostic.jl` calls
  `molecular_rayleigh_properties(RefractiveIndexRayleigh(), species, λ)` for
  N₂, O₂, CO₂, Ar: N₂ and O₂ change beyond 436/466 nm (table above; ρ
  unchanged). Its windows 740–780, 1550–1700, 1950–2100, 2200–2400 nm see
  N₂ +0.14…+0.21 % and O₂ −0.53…−0.89 %. Its `_rayleigh_extrapolation_masks`
  flags only CO₂/Ar. The N₂ and O₂ dispersion fits are now measured-range fits
  up to 2059 nm (N₂) and 1542 nm (O₂); the 1550–2400 nm windows still
  extrapolate O₂ (and N₂ beyond 2059 nm), smoothly, and the King factors of
  both gases are extrapolations beyond the visible. ExoOptics may want to flag
  these.
- ExoOptics production (`production_base.yaml`, `depol: -1`) gets the
  dispersive ρ and the τ_rayl rescaling. For its five-layer profile
  (column-mean 246.2 K): ρ +5.5 % / τ_rayl +0.25 % at 400 nm, −6.4 % / −0.30 %
  at 760 nm, −9.3 % / −0.44 % at 1625 nm, −9.6 % / −0.46 % at 2025 nm,
  −9.7 % / −0.46 % at 2300 nm. `ExoEarthRunner._tau400_diagnostic`
  recomputes ρ and τ at 400 nm with the two-argument (constant)
  `getRamanAtmoConstants`, so its recorded `rayleigh_depolarization_400nm`
  and `tau_rayleigh_400nm` no longer equal what the RT used (ρ ≈ +5.5 %,
  τ ≈ +0.25 % at 400 nm). Suggested ExoOptics change: call
  `vSmartMOM.InelasticScattering.rayleigh_air_depolarization(params.rayleigh_anisotropy, 400, T)`
  or read the model's band depol. No syntax change is needed.
- uni_vSmartMOM has its own copy of the Raman module; nothing it calls is
  removed or renamed.

## Open questions (for the scientific owner)

1. Raman/Cabannes cross sections: keep the constant tensor (published tables,
   current default) or switch to `DispersiveAnisotropy()` (σ_RRS −6.6 % at
   760 nm, +6.3 % at 395 nm, −9.6 % at 1.6 µm; ϖ_Cab +2.1×10⁻³ at 760 nm).
   The dispersive one is physically consistent with the Rayleigh default and
   removes the elastic/Raman mismatch above; it should first be checked
   against an independent RRS reference (measured N₂/O₂ rotational Raman cross
   sections or Ring spectra).
2. O₂ anisotropy: the dispersive γ follows decision 13a (King factor 1.09).
   The independent Hohm (1994) dispersion supports the 1.096 level (γ 3 %
   higher, ρ_O₂ 6 % higher).
3. Bodhaine Eq. 30 beyond 1 µm: τ_rayl is +1.0 % (1.6 µm), +3.6 % (2.06 µm)
   and +5.9 % (2.3 µm) above Peck–Reeder/Ciddor refractivity with Bodhaine's
   King factor; the King rescaling also uses a constant ρ₀ = 0.0279 instead of
   Bodhaine's F_air(λ), which leaves τ_rayl 0.3 % high at 400 nm and 0.1 % low
   at 1.6 µm. Options: an
   RI-based air σ (e.g. from the `UsingRefInd` gases) for τ_rayl, or Bodhaine's
   F_air(λ) as the rescaling reference.
4. H₂ Raman-module anisotropy (0.288×10⁻²⁴ cm³) is 8.5 % below the ab initio
   550 nm value (ρ 11–22 % low; σ_RRS(H₂) ∝ γ² would be ~16 % low). Relevant
   for H₂-rich (exoplanet) Raman runs.
5. N₂ join at 436.3 nm (continuous) vs the conventional 468 nm switch of
   Bates/Bodhaine (jump 2.85×10⁻⁴ in σ).
6. Compressibility: all fits use ideal-gas densities (σ high by 1–2×10⁻³);
   Egan gives the virial coefficients needed to correct O₂.

## Validation record

Tests run on this branch (CPU, Julia 1.12.5): `test_molecular_rayleigh.jl`
(1111 pass), `test_rayleigh_anisotropy.jl` (82), `test_parameters_parser.jl`,
`test_vertical_grid_psurf.jl` (175, golden comparison with
`ConstantAnisotropy()`), `test_ray_cab_raman.jl` (38), `test_raman_placzek_teller.jl`,
`test_update_model.jl`, `test_update_model_transaction.jl`, `test_float32.jl`,
`test_forward_raman_phase1b.jl` with `PHASE1B_CPU=1` (12, sanghavi reference,
rtol 0.02), `test_io_exports.jl`, `test_canopy.jl`, `test_jacobians_unit.jl`,
`test_lin_band_selection.jl`, `test_julia_interfaces.jl`,
`test_external_solar_raman.jl`, `test_forward_ss.jl`,
`test_inelastic_vs_smoke.jl`, `test_rt_entry_guards.jl`,
`test_type_stability.jl`: all pass (one pre-existing `broken` in
`test_update_model.jl`). Reviews: independent Claude reviewer and Codex
(gpt-6.1-sol, xhigh) both APPROVE-WITH-FIXES; fixes applied; Codex re-review
APPROVE.
