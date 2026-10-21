# Case C vector solar_tester conventions

**Status**: Resolved as of 2026-05-06 for the checked
`GaussLegQuad/Float64/CPU` axis. The large Case C residuals were from
comparing vSmartMOM at `raz = 0°` against VLIDORT vector truth at `raz = 10°`,
plus a transient `l_trunc = 16` mismatch. With the corrected setup, I/Q/U all
match VLIDORT at about `1e-3` relative error or better.

## Current setup

VLIDORT source facts from `V2p8p3_solar_tester.f90`:

- `NSTREAMS = 8`, so Task 1 uses `NMOMENTS = 2*NSTREAMS - 1 = 15`.
- The vector saved results use `USER_RELAZMS = [10, 90, 170]`, unlike scalar
  solar_tester which uses `[0, 90, 180]`.
- Task 1 has `DO_DELTAM_SCALING = .false.`, so vSmartMOM uses
  `NoTruncation()`.
- The vector aerosol moments come from `ProblemIII.Moms` as file columns
  `a1, b1, a2, a3, b2, a4`.
- The VLIDORT driver applies `SMASK = (/ 1, -1, -1, 1, 1, -1, 1, 1 /)` when
  forming `GREEKMAT`; for the I/Q coupling entries this means final VLIDORT
  `GREEKMAT(1,2)` uses `-b1`.
- VLIDORT Rayleigh uses `PROBLEM_RAY(5,2) = -sqrt(6)*beta2`, opposite in sign
  to vSmartMOM's native `get_greek_rayleigh` γ.

The vSmartMOM Case C comparison currently keeps a single internal convention:

- `solar_tester_vector.yaml`: `GaussLegQuad()`, `l_trunc = 15`,
  `max_m = 16`, `NoTruncation()`, `Stokes_IQU()`, one spectral point, and
  the first baseline comparison at `raz = 10°`.
- Aerosol import: `PROBLEMIII_b1 -> γ` directly. Do **not** apply VLIDORT's
  `SMASK` here; those file values are pre-driver values.
- Truth comparison: flip VLIDORT truth Q/U signs to compare against
  vSmartMOM's internal convention.
- Exact principal-plane U would be zero by symmetry, but the shipped vector
  saved-results geometry is not exact principal plane for the first/last RAZ
  entries (`10°` and `170°`).

Current `GaussLegQuad/Float64/CPU` residuals at
`sza = 35°`, `raz = 10°`, Task 1:

```text
TOA-up I: max = 0.000487
BOA-dn I: max = 0.000608
TOA-up Q: max = 0.000350
BOA-dn Q: max = 0.00114
TOA-up U: max = 0.000356
BOA-dn U: max = 0.00121
```

## Regression caught

Claude's attempted fix changed Case C to `γ = -PROBLEMIII_b1` and tightened
the BOA/Q gates. That mixed post-SMASK aerosol γ with native vSmartMOM
Rayleigh γ. The same sweep also compared vSmartMOM at `raz = 0°` against
VLIDORT vector truth at `raz = 10°`, so the printed residuals below include a
geometry mismatch:

```text
Case C [GaussLegQuad/Float64/CPU], vSmartMOM raz=0 vs VLIDORT raz=10:
  TOA-up I max ≈ 0.00715
  BOA-dn I max ≈ 0.189
  TOA-up Q max ≈ 0.137
```

The same edit also left `solar_tester_vector.yaml` at `l_trunc = 16`, which is
an off-by-one mismatch with VLIDORT's `NMOMENTS = 15`.

## Remaining question

No Case C source-path discrepancy is currently indicated by the corrected
comparison. A black-surface mixed Rayleigh+aerosol vector case would still be a
useful future diagnostic to isolate atmospheric coupling from Lambertian
feedback, but it is no longer needed to explain the old residuals.
