# Per-J Q-branch lines and the VS Q-branch phase matrix: Placzek–Teller weight b_JJ (decisions 5 and 7; R2, R4)

Decided by Suniti Sanghavi on 2026-10-19 after the reference check below.
Code: `compute_σ_Rayl_VibRaman_coeff_hires!` in
`src/Inelastic/src/inelastic_cross_section.jl`; test:
`test/test_raman_placzek_teller.jl`. The VS Q-branch phase matrix (R4) is in
its own section below.

## What changed

Each per-J Q-branch (ΔJ = 0) line uses the Chandrasekhar depolarization
factor of a line whose anisotropic part carries the Placzek–Teller weight
b_JJ:

| | Before | After |
|---|---|---|
| γ_C(J) | `3/(4+45(α/(b_JJ·γ))^2)` = 3b_JJ²γ²/(45α² + 4b_JJ²γ²) | `_q_branch_depolarization(α, γ, b_JJ)` = `3/(4 + 45(α/γ)^2/b_JJ)` = 3b_JJγ²/(45α² + 4b_JJγ²) |

with b_JJ = J(J+1)/((2J−1)(2J+3)) (Sanghavi 2022, Eq. 31). The change applies
to the elastic Cabannes line, with (ᾱ, γ̄), and to the vibrational Q branch,
with (α′, γ′). The unchanged prefactor α²(1+2γ_C)/(3−4γ_C) equals α²/3 + (2/27)X
for anisotropic weight X, so each Q line is now

    σ_Q(J) ∝ g_J (2J+1) N_J [α²/3 + (2/27) b_JJ γ²].

The old form entered without a citation in commit `66f5b3c1` (2022-02-10). It
squared b_JJ, apparently a parenthesization slip.

## Why b_JJ, to the first power

1. The O/S lines of RRS and RVRS carry b_{J±2,J}·γ² to first power
   (Sanghavi 2022, Eqs. 17–22 and 33–34, from Long 1977). The coefficients
   b_{J±2,J} (Eq. 19) and b_JJ (Eq. 31) are the Placzek–Teller coefficients
   (Placzek & Teller 1933), i.e. the squared Clebsch–Gordan coefficients
   ⟨J0;20|J′0⟩². By completeness, b_{J−2,J} + b_{J,J} + b_{J+2,J} = 1 for
   every J. The isotropic part α² (rank 0) contributes only to ΔJ = 0.
2. Rayleigh scattering is the sum of the Cabannes line and RRS (Sanghavi &
   Frankenberg 2023, Eq. 3; Sanghavi 2026, Eqs. 1–3). With (1), the Cabannes Q
   line is left with exactly α²/3 + (2/27)·b_JJ·γ².
3. The same holds for the vibrational band. Its total α′²/3 + (2/27)γ′² (the
   summed shape of Sanghavi 2022, Eq. 28) is shared between the VS Q line and
   the RVRS O/S lines only if the Q line carries b_JJ·γ′².
4. Independent check: Young (1981), Eq. 3, ρ_n^t = 2ρ_v^t/(1+ρ_v^t) =
   8ρ_v^c/(1+8ρ_v^c). This holds exactly when the Cabannes line carries ¼ of the
   anisotropy, which is the large-J limit of b_JJ.
5. Comparison with Sanghavi (2026), Table 1. The code's energy-conservation
   route (`compute_γ_mol_Cabannes!`) reproduces the table at 300 K (O₂ 0.00746 /
   0.00805 and N₂ 0.00265 / 0.00280 at 395 / 760 nm; `test_ray_cab_raman.jl`).
   The J-summed per-J Cabannes lines (⟨b_JJ⟩ ≈ 0.26) give 0.00739 / 0.00793 (O₂)
   and 0.00263 / 0.00277 (N₂) at 300 K, 0.8–1.5 % lower. The difference is the
   (ν̃ − Δν̃)⁴ weighting of the RRS lines, which the per-J sum leaves out.
   b_JJ² would give about 0.002 (O₂) and 0.0007 (N₂), and the J-independent
   form of Sanghavi (2022), Eq. 29, gives γ_C,Ray itself (0.0277 / 0.0102 at
   395 nm).

## Where this departs from Sanghavi (2022)

- **Eq. 29** uses a J-independent γ_C,vib with the full γ′². That gives the
  VS Q branch all of γ′², while Eqs. 33–34 give the O/S branches another ≈ ¾,
  so the anisotropy is counted twice.
- **Eq. 32** multiplies the whole line, including α′², by b_JJ. That scales
  the isotropic part by ≈ ¼, which Placzek–Teller does not do.
- **The code's per-J line strengths** now follow Placzek–Teller (points 1–3
  above), and so does the VS Q-branch *phase matrix* (R4, next section). It
  used the J-independent depolarization of Eq. 29 with the full γ′²
  (γ_C = 0.085 for N₂ and 0.170 for O₂; Sanghavi 2022, Table 4).

## VS Q-branch phase matrix: ⟨b_JJ⟩ (decision 7, R4)

Decided by Suniti Sanghavi on 2026-10-19 (decision 7, "fix all as cleanup").
Code: `get_greek_raman_VS` in `src/Inelastic/inelastic_helper.jl` and
`q_branch_mean_b` in `src/Inelastic/src/inelastic_cross_section.jl`; test:
`test/test_raman_placzek_teller.jl`, testset "VS Q-branch phase matrix uses
the line-weighted ⟨b_JJ⟩".

The VS runs treat the vibrational Q branch (ΔJ = 0) of each molecule as one
set of lines with one phase matrix; the O and S branches are separate RVRS
lines with the RRS phase matrix (γ_C = 3/4). Each per-J Q line is a
Rayleigh-type line with isotropic part α′²/3 and anisotropic part
(2/27) b_JJ γ′² (section "What changed"). The phase matrix of such a line is
linear in the two parts, so the sum of the lines is again a Rayleigh-type
line, with anisotropic part (2/27) ⟨b_JJ⟩ γ′²:

    ⟨b_JJ⟩ = Σ_J w_J b_JJ / Σ_J w_J,   w_J = g_J (2J+1) exp(−hc E(v,J)/(k T)),
    γ_C,Q = 3⟨b_JJ⟩γ′² / (45α′² + 4⟨b_JJ⟩γ′²) = _q_branch_depolarization(α′, γ′, ⟨b_JJ⟩),
    ρ_Q = 2γ_C,Q / (1 + γ_C,Q),

summed over the levels J = 0…Jmax of the per-J lines, with v = 0 for
`VS_0to1(_plus)` and v = 1 for `VS_1to0(_plus)`. T is the temperature at which
the per-J lines were computed: `compute_σ_Rayl_VibRaman_coeff_hires!` now
stores it in `mol.effCoeff.T` (the field existed but was never set), so the
phase matrix and the line strengths describe the same population. The ν⁴
factors of the individual lines (the Q-branch lines spread over a few cm⁻¹)
are left out of the weights; they change ⟨b_JJ⟩ by much less than 10⁻⁴.

| 395 nm, 300 K | ⟨b_JJ⟩ | γ_C before → after | ρ before → after | β₂ before → after |
|---|---|---|---|---|
| N₂ | 0.2552 | 0.0850 → 0.0237 | 0.1566 → 0.0463 | 0.391 → 0.466 |
| O₂ | 0.2595 | 0.1695 → 0.0528 | 0.2899 → 0.1004 | 0.310 → 0.428 |

⟨b_JJ⟩ hardly depends on T (N₂ 0.2563 at 220 K, O₂ 0.2624) and v = 1 differs
from v = 0 by less than 5 × 10⁻⁴. `γ_C_VibRaman` and `rho_depol_VibRaman` keep
the whole-band value of Eq. 29; it still enters `compute_σ_VibRaman_coeff!`
(whole-band cross section) and the stellar H₂ `get_greek_raman_VS`, which is
not changed here.

**Effect** (`VS_0to1_plus` at 395 nm, pure Rayleigh, one layer with
τ_Rayleigh = 0.378, SZA 45°, flat F₀; measured with the corrected builder of
decision 8). The elastic radiances do not change, and neither does the
vibrational Raman energy (the line strengths are the same). Its angular and
polarization distribution does: the band-summed inelastic I changes by
+4.1 % at nadir, −4.1 % at VZA 45° (Δφ = 0°), +3.1 % at VZA 45° (Δφ = 180°)
and −1.9 % at VZA 60° (Δφ = 90°); its Q/I goes from 0.208 to 0.243 at nadir
and from 0.658 to 0.850 at VZA 45° (Δφ = 0°).

## Effect of the b_JJ line weights

- **σ_VRS rises.** σ_VRS (vibrational Q branch, `inelastic_helper.jl`) rises by
  13.3 % (O₂) and 5.9 % (N₂) relative to the old code, at 250–300 K and
  395–770 nm. Put the other way, the old values were 11.7 % and 5.6 % low.
- **Other consumers of these lines change too:**
  - `compute_optical_RS!` / `compute_optical_VRS!` for VS and VS_plus runs;
  - the VS-specific Cabannes fractions (`compute_ϖ_Cabannes_VS` and the
    VS_plus band ϖ_Cabannes, e.g. 0.999176 → 0.999128 at 395 nm);
  - the stellar H₂ `sol_VS_*` path, by +5.4 to +5.6 % from about 250 K upward
    (+4.2 % at 100 K).
- **Stored but unused:** the per-J elastic line `σ_Rayl_coeff_hires`.
- **Unchanged:** RRS, RVRS, Rayleigh, the RRS-based Cabannes quantities
  (ϖ_Cab, γ_C,Cab) and every noRS/RRS result.

## Tests

`test/test_raman_placzek_teller.jl` checks the following.

- **Per-J identity.** Q + S + O = 128π⁵ g_J(2J+1) N_J (α²/3 + 2γ²/27) holds to
  1e-10 for every J ≤ Jmax−2. It is checked for the elastic lines and the
  0→1 and 1→0 vibrational lines, for O₂ and N₂, at 250 and 296 K.
  - **Populations.** Each branch is normalized by its own population factor:
    v = 0–2 for the elastic Q line and v = 0 for the RRS lines. That is how the
    code defines them; the small difference between the two, about 6e-4 for
    O₂, is left as it is.
- **Float32.** `_q_branch_depolarization` stays finite in Float32.
- **J-summed Cabannes line.** Its depolarization agrees with the
  energy-conservation route `compute_γ_mol_Cabannes!` to within 2 % at 395 and
  760 nm, and its anisotropy fraction lies in 0.24–0.28.
- **Discrimination.** With the old b_JJ², 176 of the original 210 checks fail.
- **VS Q-branch phase matrix (R4).** `q_branch_mean_b` equals the anisotropy
  fraction X/γ′² of the summed per-J Q lines to 1e-10 (v = 0 and 1, N₂ and O₂,
  220 and 300 K) and lies in 0.24–0.28; the depolarization carried by the
  Greek coefficients of `get_greek_raman_VS` equals
  3⟨b_JJ⟩γ′²/(45α′² + 4⟨b_JJ⟩γ′²) to 1e-10 and is below a third of the
  whole-band value.

## References

- Long, D.A. (1977). *Raman Spectroscopy.* McGraw-Hill.
- Placzek, G. & Teller, E. (1933). Die Rotationsstruktur der Ramanbanden
  mehratomiger Moleküle. *Z. Phys.* 81, 209–258.
- Sanghavi, S. (2022). Raman scattering in the Earth's atmosphere, Part I:
  Optical properties. *JQSRT* 291, 108328.
- Sanghavi, S. & Frankenberg, C. (2023). Raman scattering in the Earth's
  atmosphere, Part II: Radiative transfer modeling for remote sensing
  applications. *JQSRT* 311, 108791.
- Sanghavi, S. (2026). Isolating Cabannes and rotational Raman scattering to
  determine forward modeling error due to Rayleigh scattering in the Earth's
  atmosphere (manuscript, JQSRT).
- Young, A.T. (1981). On the Rayleigh-scattering optical depth of the
  atmosphere. *J. Appl. Meteor.* 20, 328–330.
