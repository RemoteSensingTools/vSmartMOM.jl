# Log-AOD prior sensitivity test

This directory records the isolated test of
`sigma(ln(AOD760)) = 0.75` for all three aerosol species that preceded its
adoption as the production prior on 2026-09-03.

The test keeps the production a-priori state means and every other covariance
entry unchanged. Only the three independent log-AOD variances change:

```text
parameter                         production variance   test variance
ln_sulfate_aod760                 4.0                   0.5625
ln_organic_carbon_aod760          4.0                   0.5625
ln_utls_sulfate_aod760            4.0                   0.5625
```

The authoritative production prior remains
`../../retrieval_setup/apriori_states.nc`, SHA-256
`6497919cd35ab8d05f79e96d827e05e567576931c9eae1b8cf534151c1f0ea2c`.
The isolated test prior is `retrieval_setup/apriori_states.nc`, SHA-256
`147d42466735bb33e390b437754c94cba0362bb6288b61cacac04b4c4b0f2612`.

The retrieval comparison uses the noiseless member (perturbation 11) of two
no-SIF truth states:

- state 001: urban, no aerosol, uniform 400 ppm background with a 360 ppm
  bottom layer;
- state 013: urban, aerosol AOD760 = 0.28, uniform 400 ppm CO2.

Both corrected and uncorrected measurements are retrieved with Float32,
nine streams, `noRS`, and the `OCO_RRS_synth` analytic Jacobian. Outputs are
written only below `retrievals/`; they never replace production retrievals or
publish a production smoke-release gate. The Wurst GPU0 launcher is
`../../../inversion/run_ln_aod_sigma_0p75_sensitivity.sh`, and its combined
log is `run_wurst_gpu0.log`.

## Status

The test was launched on Wurst GPU0 on 2026-09-03. Both state-001 retrievals
converged in three accepted trials with no divergences. The state-013 corrected
retrieval converged with the tighter prior, whereas its sigma-2 baseline had
stopped after the maximum number of divergent LM trials. The state-013
uncorrected sensitivity run was deliberately interrupted when the production
restart was approved, so this directory contains three of the four planned
products and is not a complete comparison suite.

The canonical production prior and restart record are now
`../../retrieval_setup/apriori_states.nc` and
`../../RETRIEVAL_RESTART_SIGMA_LNAOD_0P75.md`, respectively.
