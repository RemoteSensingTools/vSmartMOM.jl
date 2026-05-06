# Case C (vector solar_tester) — open Q sign + U principal-plane conventions

**Status**: Surfaced 2026-05-06 by `test/vlidort_baseline/cases/case_C_solar_tester_vector.jl`. Branch `SS-exact`. Open.

## Symptom

Vector RT regression test against VLIDORT 2.8.3 `solar_tester_IQU0.all`
(NSTOKES=3, Problem III aerosol + Rayleigh + Lambertian). At sza=35°,
raz=0°, vza=10° (geom 1, level 1, dir up, task 1):

| Stokes | modeled | truth | rel-err |
| ------ | ------- | ----- | ------- |
| I      | 0.0695  | 0.0696 | 1.4e-3 |
| Q      | +0.0134 | −0.0128 | ~2.05 (sign flipped) |
| U      | 0.0     | −0.0037 | 1.0 (modeled exactly 0) |

I matches well (~7e-3 max across VZAs). Q has the **wrong sign**. U is
**exactly zero** in vSmartMOM at vaz=0°.

## Diagnosis

### Q sign — Rayleigh γ convention difference

Compared B[1,2] = γ at L=2 between codes:

- **vSmartMOM** (`get_greek_rayleigh` in `src/Scattering/mie_helper_functions.jl:454`):
  γ[3] = `0.5 · sqrt(6) · dpl_p` ≈ **+1.205** (with dpl_p = 0.984)
- **VLIDORT** (`V2p8p3_solar_tester.f90:336`):
  PROBLEM_RAY(5,2) = `−sqrt(6) · β_2` ≈ **−1.205** (with β_2 = 0.4920)

These have **opposite signs** for the same Rayleigh L=2 coefficient.

In Case A (Siewert IIA, τ_rayl ≡ 0 — pure aerosol), Q matches at ~1e-6
because Rayleigh never enters. In Case C, the Rayleigh γ contribution
flows into the Z matrix with the opposite sign vs VLIDORT, so the I↔Q
coupling has the wrong sign for any Rayleigh-influenced layer — hence
the global Q sign flip.

Aerosol γ mapping (via PROBLEMIII_b1 → vSmartMOM γ, **no sign flip**) is
consistent with Case A's Siewert PROBLEM_IIA(3,:) → vSmartMOM γ
convention. Tried mirroring VLIDORT's SMASK = −1 on b1 (γ = −b1 in
vSmartMOM): I rel-err got worse (vza=40 went from 4e-4 to 7e-3) and Q
sign was unchanged. So keep the Siewert-style aerosol mapping; the Q
sign flip is purely from Rayleigh γ.

### U at vaz = 0° / 180° — postprocessing weight by sin(m·φ)

In `src/CoreRT/tools/postprocessing_vza.jl:33`:
```julia
w = if n == 1
    weight * cos_m_phi
else
    weight * Diagonal([cos_m_phi, cos_m_phi, sin_m_phi, sin_m_phi][1:n])
end
```

For Stokes_IQU (n=3): U is weighted by `sin(m·φ)`. At φ = 0° or 180°,
`sin(m·φ) = 0` for all m → U is identically 0 in the principal plane.

VLIDORT's `solar_tester_IQU0.all` reports a non-zero U at raz=0° (e.g.
−4e-3 at vza=10°), which means VLIDORT uses a different sign/rotation
convention that produces U ≠ 0 in this geometry. Likely a different
Stokes-vector axis orientation.

## Knock-on effect: BOA-dn I

The Rayleigh γ sign mismatch also leaks into the I↔Q multi-scatter
coupling for downwelling. Observed at sza=35°, raz=0° (geom 1, level
5/BOA, dir=2):

| Direction | I rel-err | comment |
| --------- | --------- | ------- |
| TOA-up I  | ~2e-3 max | OK |
| BOA-dn I  | ~19% max  | dominated by multi-scatter; γ-driven |

So Case C currently gates TOA-up I tightly (5e-3) and BOA-dn I loosely
(3e-1). Reconciling Rayleigh γ should bring BOA-dn I back to ~1e-3.

## How Case C currently handles it

The committed Case C tests gate **only on I**, with TOA-up tight and
BOA-dn loose (see above). Q and U are computed and printed for
visibility but **not asserted**. The docstring at the top of
`case_C_solar_tester_vector.jl` lists the open conventions.

## Resolution paths

1. **Q (Rayleigh sign)**:
   - Either flip vSmartMOM's `get_greek_rayleigh` γ in source (changes
     wider behavior — would need to verify Q sign convention with the
     full vSmartMOM test suite)
   - Or override Rayleigh greek_coefs in the test setup, computing γ
     explicitly with VLIDORT's sign and assigning to
     `model.optics.rayleigh.greek_rayleigh.γ`. Test-only fix.
2. **U (principal-plane)**: investigate VLIDORT's U convention at
   raz=0°. May involve a 180° rotation/sign convention. If VLIDORT's
   non-zero U at raz=0° is due to a phase convention rather than
   physics, document and skip in test. Otherwise add a follow-up
   raz=90° run that exercises U non-trivially.

## Why this isn't a Case A issue

Case A uses Siewert PROBLEM_IIA aerosol with **τ_rayl explicitly set to
zero** in `case_A_siewert2000.jl:34`. So no Rayleigh contribution → no
Q sign flip. Also Case A's Siewert truth tables only report U at az=90°
(Table 6, ~1e-9 — essentially zero), so the U principal-plane issue
never appears.
