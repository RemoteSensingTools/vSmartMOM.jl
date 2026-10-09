# Round-5 fixed-SIF and tight-UTLS prior

Round 5 is designed as an apples-to-apples comparison with round 4. It uses
the same truth spectra, noise realizations, forward model, analytical
Jacobians, surface priors, tropospheric-aerosol priors, and CO2 prior. Its
scientific changes are deliberately limited to the SIF and stratospheric
aerosol coordinates described below.

## Exact changes from round 4

Both SIF quantities are known and fixed:

```text
SIF off: Lnu(759 nm) = 0; mSIF = 0
SIF on:  Lnu(759 nm) = 0.004818031987713776
         mSIF        = 1.2291230681458325e-5
```

The SIF-on values are in the native vSmartMOM units
`mW m-2 sr-1 (cm-1)-1` and `mW m-2 sr-1 (cm-1)-2`, respectively. The core
forward model remains referenced at 760 nm and receives

```text
SIF760 = Lnu759 + mSIF * (nu760 - nu759).
```

Neither SIF coefficient is an active retrieval coordinate. Both SIF-on and
SIF-off retrievals therefore have the same 28-coordinate numerical state:

```text
full-state indices: [1; 6:32]
core-state indices: 1:28
```

The UTLS/stratospheric sulfate AOD and profile-height log-coordinate standard
deviations are each reduced by exactly a factor of ten:

| Coordinate | Round 4 sigma | Round 5 sigma | Variance ratio |
|---|---:|---:|---:|
| `ln(AOD760_UTLS)` | 0.75 | 0.075 | 0.01 |
| `ln(z0_UTLS/km)` | 0.10 | 0.01 | 0.01 |

The sulfate and organic-carbon tropospheric values stay at `0.75` for
`ln(AOD760)` and `0.10` for `ln(z0/km)`. The dedicated prior builder applies
the tightening as a covariance congruence transform, so any future
cross-covariances involving those coordinates retain their correlations.

## CO2 is intentionally unchanged

Round 5 does **not** alter CO2 flexibility. The prior builder copies the full
16-element CO2 mean and full `16 x 16` covariance block from the approved
round-3/round-4 tapered prior and tests them for exact equality. This includes
the four fixed upper layers and every off-diagonal term among the active lower
layers. Consequently, differences between rounds 4 and 5 cannot be attributed
to a changed CO2 prior.

A proposed CO2 precision-matrix modification remains deferred to a separate,
unnamed experiment. The approved round 6 instead fixes SIF while restoring
the original UTLS prior; see [ROUND6_FIXED_SIF_PRIOR.md](ROUND6_FIXED_SIF_PRIOR.md).
Neither round changes the CO2 covariance.

## Implementation and validation

- `build_round5_fixed_sif_apriori.jl` generates separate SIF-on and SIF-off
  NetCDF priors plus human-readable `.dat` audits.
- `Round5FixedSIF.jl` injects the two fixed SIF values into the established
  30-coordinate `OCO_RRS_synth` core and removes their Jacobian columns.
- `Round5RetrievalCampaign.jl` validates the state mapping, prior identity,
  truth ownership, and output provenance.
- `run_round5_fixed_sif_retrievals.jl` reads existing truth/noise products and
  writes into a visibly isolated `round5/fixed_sif` output namespace.

The validation suite checks the 28-to-30 state mapping, analytical-Jacobian
column reduction against two-sided finite differences, exact preservation of
the complete CO2 block, factor-of-ten UTLS standard deviations, positive
definiteness of every active covariance, and round-5 provenance barriers.
