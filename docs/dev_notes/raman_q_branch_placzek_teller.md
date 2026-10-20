# Per-J Q-branch line strengths: Placzek–Teller weight b_JJ (decision 5, R2)

Decided by Suniti Sanghavi on 2026-10-19 after the reference check below.
Code: `compute_σ_Rayl_VibRaman_coeff_hires!` in
`src/Inelastic/src/inelastic_cross_section.jl`; test:
`test/test_raman_placzek_teller.jl`.

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
  above). The VS Q-branch *phase matrix* (`get_greek_raman_VS`, through
  `rho_depol_VibRaman`) still uses the J-independent depolarization of Eq. 29
  with the full γ′² (γ_C = 0.085 for N₂ and 0.170 for O₂; Sanghavi 2022,
  Table 4). The value consistent with these line strengths, using ⟨b_JJ⟩, is
  γ_C ≈ 0.024 / 0.053. That is the
  open Raman item R4 ("band-depolarized Q-branch phase"), not part of this change.

## Effect

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

`test/test_raman_placzek_teller.jl` checks two things.

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
